"""Command-line interface for the standalone ``merge_vcf`` helper."""

from __future__ import annotations

import argparse
import datetime
import os
from types import SimpleNamespace
from typing import Any, List, Optional, Sequence, Tuple

from MetaGap.MetagapUserCode.merge_vcf.cli import summarize_produced_vcfs
from MetaGap.MetagapUserCode.merge_vcf.logging_utils import handle_critical_error, log_message
from MetaGap.MetagapUserCode.merge_vcf.metadata import (
    _parse_sample_metadata_line,
    append_metadata_to_merged_vcf,
    parse_metadata_arguments,
)
from MetaGap.MetagapUserCode.merge_vcf.merging import merge_vcfs
from MetaGap.MetagapUserCode.merge_vcf.validation import (
    find_first_vcf_with_header,
    normalize_vcf_version,
    validate_all_vcfs,
    validate_merged_vcf,
)


def _build_parser() -> argparse.ArgumentParser:
    """Return the argument parser for the standalone CLI."""

    parser = argparse.ArgumentParser(
        description=(
            "Merge a directory of VCF shards into an anonymised cohort file while "
            "injecting metadata from a template file."
        )
    )
    parser.add_argument(
        "input_dir",
        help="Directory containing the VCF files to merge.",
    )
    parser.add_argument(
        "output_dir",
        help="Destination directory for merged and anonymised VCF artefacts.",
    )
    parser.add_argument(
        "metadata_template",
        help="Path to a metadata template file containing ## header definitions.",
    )
    parser.add_argument(
        "--ref",
        dest="ref",
        help=(
            "Expected reference genome build. When omitted the workflow attempts to "
            "infer it from the first VCF header."
        ),
    )
    parser.add_argument(
        "--vcf-version",
        dest="vcf_version",
        help=(
            "Expected VCF version (for example 4.2). When omitted the workflow attempts "
            "to infer it from the first VCF header."
        ),
    )
    parser.add_argument(
        "--allow-gvcf",
        action="store_true",
        help="Permit input shards that include gVCF specific annotations.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Echo progress messages to stdout in addition to the log file.",
    )
    return parser


def _load_template_metadata(
    template_path: str,
) -> Tuple[List[str], List[str], Optional[str]]:
    """Return sample metadata entries and header lines parsed from *template_path*."""

    if not template_path:
        handle_critical_error("A metadata template path must be supplied.")

    if not os.path.exists(template_path):
        handle_critical_error(f"Metadata template does not exist: {template_path}")

    if not os.path.isfile(template_path):
        handle_critical_error(f"Metadata template is not a file: {template_path}")

    sample_entries: List[str] = []
    header_lines: List[str] = []
    serialized_sample_line: Optional[str] = None

    try:
        with open(template_path, "r", encoding="utf-8") as handle:
            for raw_line in handle:
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#CHROM"):
                    continue
                if stripped.startswith("##SAMPLE="):
                    if serialized_sample_line is not None:
                        handle_critical_error(
                            "Metadata template contains multiple ##SAMPLE definitions."
                        )
                    serialized_sample_line = stripped
                    mapping = _parse_sample_metadata_line(stripped)
                    sample_entries = [
                        f"{key}={value}"
                        for key, value in mapping.items()
                        if value is not None
                    ]
                    continue
                if stripped.startswith("#") and not stripped.startswith("##"):
                    continue
                header_lines.append(stripped)
    except OSError as exc:  # pragma: no cover - defensive error handling
        handle_critical_error(f"Failed to read metadata template {template_path}: {exc}")

    return sample_entries, header_lines, serialized_sample_line


def _prepare_metadata_from_template(
    template_path: str,
    verbose: bool,
) -> Tuple[Optional[Any], Sequence[Any], Optional[Sequence[str]], Sequence[str], Optional[str]]:
    """Load metadata from *template_path* and normalise it for the workflow."""

    sample_entries, header_lines, serialized_sample_line = _load_template_metadata(template_path)

    namespace = SimpleNamespace(
        sample_metadata_entries=sample_entries,
        header_metadata_lines=header_lines,
    )

    (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries_mapping,
        sanitized_header_lines,
        serialized_sample_line_from_args,
    ) = parse_metadata_arguments(namespace, verbose)

    serialized_line = serialized_sample_line_from_args or serialized_sample_line

    if sample_entries:
        log_message(
            "Loaded sample metadata definition from template: "
            + (serialized_line or "<invalid sample metadata>"),
            verbose,
        )

    if header_lines:
        log_message(
            "Loaded additional header metadata lines from template: "
            + ", ".join(header_lines),
            verbose,
        )

    return (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries_mapping,
        sanitized_header_lines,
        serialized_line,
    )


def run(
    input_dir: str,
    output_dir: str,
    metadata_template: str,
    *,
    ref: Optional[str] = None,
    vcf_version: Optional[str] = None,
    allow_gvcf: bool = False,
    verbose: bool = False,
) -> None:
    """Execute the merge workflow with parameters supplied by the CLI."""

    (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries,
        header_metadata_lines,
        serialized_sample_line,
    ) = _prepare_metadata_from_template(metadata_template, verbose)

    log_message(
        "Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        verbose,
    )

    input_dir = os.path.abspath(input_dir)
    if not os.path.isdir(input_dir):
        handle_critical_error(f"Input directory does not exist: {input_dir}")
    log_message("Input directory: " + input_dir, verbose)

    output_dir = os.path.abspath(output_dir)
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    log_message("Output directory: " + output_dir, verbose)

    ref_genome = ref.strip() if ref else None
    vcf_version_normalized = normalize_vcf_version(vcf_version) if vcf_version else None

    detection_needed = not ref_genome or not vcf_version_normalized
    detected_file = detected_fileformat = detected_reference = None
    if detection_needed:
        detected_file, detected_fileformat, detected_reference = find_first_vcf_with_header(
            input_dir, verbose
        )
        if not ref_genome and detected_reference:
            ref_genome = detected_reference
        if not vcf_version_normalized and detected_fileformat:
            vcf_version_normalized = normalize_vcf_version(detected_fileformat)

    if not ref_genome:
        handle_critical_error(
            "Reference genome build must be provided via --ref or detectable from input headers."
        )
    if not vcf_version_normalized:
        handle_critical_error(
            "VCF version must be provided via --vcf-version or detectable from input headers."
        )

    if detection_needed and detected_file:
        log_message(
            f"Auto-detected metadata from {detected_file} -> reference="
            f"{detected_reference or 'unknown'}, version={detected_fileformat or 'unknown'}",
            verbose,
        )

    log_message(
        f"Reference genome: {ref_genome}, VCF version: {vcf_version_normalized}",
        verbose,
    )

    valid_files, sample_order = validate_all_vcfs(
        input_dir,
        ref_genome,
        vcf_version_normalized,
        verbose=verbose,
        allow_gvcf=allow_gvcf,
    )

    merged_vcf = merge_vcfs(
        valid_files,
        output_dir,
        verbose,
        sample_order=sample_order,
        sample_header_line=sample_header_line,
        simple_header_lines=simple_header_lines,
    )

    final_vcf = append_metadata_to_merged_vcf(
        merged_vcf,
        sample_metadata_entries=sample_metadata_entries,
        header_metadata_lines=header_metadata_lines,
        serialized_sample_line=serialized_sample_line,
        verbose=verbose,
    )

    validate_merged_vcf(final_vcf, verbose)

    summary_dir, sample_filename, produced_count = summarize_produced_vcfs(
        output_dir, final_vcf
    )
    summary_path = os.path.join(summary_dir, sample_filename)
    print(f"Wrote: {summary_path} x {produced_count}")
    log_message(
        "Script execution completed successfully. "
        f"Wrote {produced_count} VCF file(s) (e.g., {summary_path}).",
        verbose,
    )


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Parse arguments and execute the merge workflow."""

    parser = _build_parser()
    args = parser.parse_args(argv)
    run(
        args.input_dir,
        args.output_dir,
        args.metadata_template,
        ref=args.ref,
        vcf_version=args.vcf_version,
        allow_gvcf=args.allow_gvcf,
        verbose=args.verbose,
    )


__all__ = ["main", "run"]


if __name__ == "__main__":  # pragma: no cover - manual execution helper
    main()
