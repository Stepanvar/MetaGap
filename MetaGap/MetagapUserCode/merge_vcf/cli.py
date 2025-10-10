"""Command-line entrypoint for the VCF merging workflow."""

from __future__ import annotations

import argparse
import datetime
import glob
import os
import shutil

from .logging_utils import (
    MergeConflictError,
    MergeVCFError,
    ValidationError,
    handle_critical_error,
    log_message,
)
from .metadata import append_metadata_to_merged_vcf, load_metadata_lines
from .merging import merge_vcfs
from .validation import (
    find_first_vcf_with_header,
    normalize_vcf_version,
    validate_all_vcfs,
    validate_merged_vcf,
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Consolidate multiple VCF files into a single merged VCF file."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the VCF files to merge.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        "--output",
        dest="output_dir",
        help="Directory for the merged VCF. Defaults to the input directory when omitted.",
    )
    parser.add_argument(
        "-t",
        "--metadata-file",
        required=True,
        dest="metadata_file",
        help=(
            "Path to a metadata template file containing header lines to apply to the merged VCF."
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose console logging in addition to the log file.",
    )
    parser.add_argument(
        "--qual-threshold",
        type=float,
        default=30.0,
        help=(
            "Minimum QUAL value required to keep a variant after merging. "
            "Set to a negative value to disable the filter."
        ),
    )
    parser.add_argument(
        "--an-threshold",
        type=float,
        default=50.0,
        help=(
            "Minimum INFO/AN value required to keep a variant after tag recalculation. "
            "Set to a negative value to disable the filter."
        ),
    )
    return parser.parse_args()


def summarize_produced_vcfs(output_dir: str, fallback_vcf: str):
    """Return a directory, representative filename, and count for output VCFs."""

    discovered = []
    output_dir = os.path.abspath(output_dir) if output_dir else None
    if output_dir and os.path.isdir(output_dir):
        for pattern in ("SAMPLE_*.vcf.gz", "SAMPLE_*.vcf"):
            discovered.extend(
                sorted(
                    glob.glob(os.path.join(output_dir, pattern)),
                    key=lambda path: os.path.basename(path),
                )
            )

    if discovered:
        discovered.sort(
            key=lambda path: (
                0 if path.endswith(".vcf.gz") else 1,
                os.path.basename(path),
            )
        )
        representative_path = discovered[0]
        directory = os.path.dirname(os.path.abspath(representative_path))
        filename = os.path.basename(representative_path)
        count = len(discovered)
        return directory, filename, count

    fallback_abs = os.path.abspath(fallback_vcf)
    return os.path.dirname(fallback_abs), os.path.basename(fallback_abs), 1


def main():
    args = parse_arguments()
    verbose = args.verbose

    try:
        metadata_header_lines = load_metadata_lines(args.metadata_file, verbose)

        log_message(
            "Script Execution Log - "
            + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            verbose,
        )

        input_dir = os.path.abspath(args.input_dir)
        if not os.path.isdir(input_dir):
            handle_critical_error(
                f"Input directory does not exist: {input_dir}",
                exc_cls=ValidationError,
            )
        log_message("Input directory: " + input_dir, verbose)

        output_dir = os.path.abspath(args.output_dir) if args.output_dir else input_dir
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        log_message("Output directory: " + output_dir, verbose)

        (
            detected_file,
            detected_fileformat,
            detected_reference,
        ) = find_first_vcf_with_header(input_dir, verbose)
        ref_genome = detected_reference
        vcf_version = normalize_vcf_version(detected_fileformat)

        if not ref_genome:
            handle_critical_error(
                "Reference genome build must be auto-detectable from input files.",
                exc_cls=ValidationError,
            )
        if not vcf_version:
            handle_critical_error(
                "VCF version must be auto-detectable from input files.",
                exc_cls=ValidationError,
            )

        if detected_file:
            log_message(
                f"Auto-detected metadata from {detected_file} -> reference={detected_reference or 'unknown'}, "
                f"version={detected_fileformat or 'unknown'}",
                verbose,
            )

        log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}", verbose)

        valid_files, sample_order = validate_all_vcfs(
            input_dir,
            ref_genome,
            vcf_version,
            verbose,
        )

        merged_vcf = merge_vcfs(
            valid_files,
            output_dir,
            verbose,
            sample_order=sample_order,
        )

        final_vcf = append_metadata_to_merged_vcf(
            merged_vcf,
            header_metadata_lines=metadata_header_lines,
            verbose=verbose,
        )

        validate_merged_vcf(final_vcf, verbose)
            merged_vcf = merge_vcfs(
            valid_files,
            output_dir,
            verbose,
            sample_order=sample_order,
            qual_threshold=args.qual_threshold if args.qual_threshold >= 0 else None,
            an_threshold=args.an_threshold if args.an_threshold >= 0 else None,
        )

        summary_dir, sample_filename, produced_count = summarize_produced_vcfs(
            output_dir, final_vcf
        )
        summary_path = os.path.join(summary_dir, sample_filename)
        log_message(
            f"Discovered {produced_count} per-sample shard file(s) in {summary_dir} (e.g., {sample_filename}).",
            verbose,
        )

        final_target = os.path.join(output_dir, "cohort_final.vcf.gz")
        if final_vcf != final_target:
            os.makedirs(os.path.dirname(final_target), exist_ok=True)
            if os.path.exists(final_target):
                os.remove(final_target)
            shutil.move(final_vcf, final_target)
            for ext in (".tbi", ".csi"):
                source_index = final_vcf + ext
                target_index = final_target + ext
                if os.path.exists(source_index):
                    if os.path.exists(target_index):
                        os.remove(target_index)
                    shutil.move(source_index, target_index)
            final_vcf = final_target
        else:
            final_target = final_vcf

        log_message(
            f"Script execution completed successfully. Final cohort VCF: {final_target}",
            verbose,
        )
        print(f"Wrote: {final_target} x 1.")
    except (ValidationError, MergeConflictError, MergeVCFError) as exc:
        print(f"ERROR: {exc}")
        raise SystemExit(1)


__all__ = [
    "parse_arguments",
    "summarize_produced_vcfs",
    "main",
]
