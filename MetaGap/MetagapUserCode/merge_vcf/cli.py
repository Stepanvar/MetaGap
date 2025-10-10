"""Command-line entrypoint for the VCF merging workflow."""
from __future__ import annotations

import argparse
import datetime
import logging
import os
import sys
from pathlib import Path

from . import merging, metadata as metadata_module, validation
from .filtering import DEFAULT_ALLOWED_FILTER_VALUES
from .logging_utils import (
    LOG_FILE,
    MergeConflictError,
    MergeVCFError,
    ValidationError,
    configure_logging,
    handle_critical_error,
    log_message,
)


# Re-export frequently patched helpers for unit tests.
validate_all_vcfs = validation.validate_all_vcfs
merge_vcfs = merging.merge_vcfs
append_metadata_to_merged_vcf = metadata_module.append_metadata_to_merged_vcf
validate_merged_vcf = validation.validate_merged_vcf


def _validate_metadata_entry(arg: str) -> str:
    """Validate KEY=VALUE and return sanitized string."""
    if arg is None:
        raise argparse.ArgumentTypeError("Metadata entry cannot be empty")
    s = arg.strip()
    if s == "":
        return s
    if "=" not in s:
        raise argparse.ArgumentTypeError("Use KEY=VALUE (e.g., ID=VALUE)")
    key, value = s.split("=", 1)
    key, value = key.strip(), value.strip()
    if key == "" or (key.lower() == "id" and value == ""):
        raise argparse.ArgumentTypeError("Metadata must include a non-empty ID")
    if key.startswith("-"):
        raise argparse.ArgumentTypeError("Metadata key cannot start with '-'")
    return f"{key}={value}"


def parse_arguments():
    """Parse CLI args for the VCF merging tool."""
    parser = argparse.ArgumentParser(
        prog="merge_vcf",
        description="Consolidate multiple VCF files into a single merged VCF.",
    )
    parser.add_argument("input_vcfs", nargs="*", help="Input VCF files.")
    parser.add_argument("--inputs", nargs="+", dest="input_list", help="Alternate way to pass input VCFs.")
    parser.add_argument(
        "-o",
        "--out-dir",
        "--output-dir",
        "--output",
        dest="output_dir",
        help="Output directory for the merged VCF. Defaults to current directory.",
    )
    parser.add_argument("--qual-threshold", type=float, default=30.0, help="Minimum QUAL to keep a variant. <0 disables.")
    parser.add_argument("--an-threshold", type=float, default=50.0, help="Minimum INFO/AN to keep a variant. <0 disables.")
    parser.add_argument(
        "--sample-header",
        action="append",
        dest="sample_header_entries",
        type=_validate_metadata_entry,
        default=[],
        help="SAMPLE header KEY=VALUE entry (repeatable). Must include ID=...",
    )
    parser.add_argument(
        "--simple-header",
        action="append",
        dest="simple_header_lines",
        type=_validate_metadata_entry,
        default=[],
        help="Additional header KEY=VALUE (repeatable). '##' is auto-prefixed.",
    )
    parser.add_argument(
        "--metadata-file",
        dest="metadata_file",
        help="Optional file with extra header lines. One entry per line. Lines without '##' will be prefixed.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose console logging.")
    args = parser.parse_args()

    # Collect inputs
    input_files: list[str] = []
    if args.input_list:
        input_files.extend(args.input_list)
    if args.input_vcfs:
        input_files.extend(args.input_vcfs)
    if not input_files:
        parser.error("No input VCF files specified.")

    # Clean metadata
    args.sample_header_entries = [s for s in args.sample_header_entries if s and s.strip()]
    args.simple_header_lines = [
        s if s.startswith("##") else f"##{s}" for s in args.simple_header_lines if s and s.strip()
    ]

    # Require SAMPLE ID if any SAMPLE metadata present
    has_id = any(
        (kv := ent.split("=", 1)) and kv[0].strip().lower() == "id" and kv[1].strip()
        for ent in args.sample_header_entries
    )
    if args.sample_header_entries and not has_id:
        parser.error("Provide a non-empty SAMPLE ID (e.g., --sample-header ID=SampleName).")

    # Normalize paths
    args.output_dir = str(Path(args.output_dir)) if args.output_dir else "."
    args.input_files = [str(Path(p)) for p in input_files]
    return args


def summarize_produced_vcfs(output_dir: str, fallback_vcf: str):
    """Return (dir, representative_filename, count) for SAMPLE_*.vcf[.gz] in output_dir, else fallback file."""
    directory = Path(output_dir).resolve()
    if directory.is_dir():
        vcf_candidates = sorted(directory.glob("SAMPLE_*.vcf.gz")) + sorted(directory.glob("SAMPLE_*.vcf"))
        if vcf_candidates:
            representative_path = vcf_candidates[0].resolve()
            return str(representative_path.parent), representative_path.name, len(vcf_candidates)
    fallback_path = Path(fallback_vcf).resolve()
    return str(fallback_path.parent), fallback_path.name, 1


def main():
    args = parse_arguments()
    verbose = getattr(args, "verbose", False)
    allow_gvcf = getattr(args, "allow_gvcf", False)
    input_dir = None
    output_dir = None
    extra_file_lines: list[str] = []
    ref_genome = getattr(args, "ref", None)
    vcf_version = getattr(args, "vcf_version", None)

    # Phase 1: environment + autodetect metadata
    try:
        # Prefer an explicit input directory but fall back to legacy positional files.
        input_dir = getattr(args, "input_dir", None)
        if input_dir:
            input_dir = os.path.abspath(input_dir)
        else:
            input_files = getattr(args, "input_files", [])
            base_dirs = [os.path.dirname(os.path.abspath(p)) or "." for p in input_files]
            input_dir = os.path.commonpath(base_dirs) if base_dirs else "."

        if not os.path.isdir(input_dir):
            raise ValidationError(f"Input directory does not exist: {input_dir}")

        output_dir_arg = getattr(args, "output_dir", None)
        output_dir = os.path.abspath(output_dir_arg) if output_dir_arg else input_dir
        os.makedirs(output_dir, exist_ok=True)

        log_path = os.path.join(output_dir, LOG_FILE)
        configure_logging(
            log_level=logging.DEBUG if verbose else logging.INFO,
            log_file=log_path,
            enable_file_logging=True,
            enable_console=True,
        )
        if not verbose:
            for handler in logging.getLogger("vcf_merger").handlers:
                if isinstance(handler, logging.StreamHandler):
                    handler.setLevel(logging.WARNING)

        metadata_file = getattr(args, "metadata_file", None)
        extra_file_lines = (
            metadata_module.load_metadata_lines(metadata_file, verbose)
            if metadata_file
            else []
        )

        log_message("Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        log_message(f"Input directory: {input_dir}")
        log_message(f"Output directory: {output_dir}")

        detected_file = None
        detected_fileformat = None
        detected_reference = None

        if not ref_genome or not vcf_version:
            (
                detected_file,
                detected_fileformat,
                detected_reference,
            ) = validation.find_first_vcf_with_header(input_dir, verbose)

            if not ref_genome:
                ref_genome = detected_reference
            if not vcf_version:
                vcf_version = validation.normalize_vcf_version(detected_fileformat)

        if not ref_genome:
            raise ValidationError("Reference genome build must be auto-detectable from input files.")
        if not vcf_version:
            raise ValidationError("VCF version must be auto-detectable from input files.")

        if detected_file:
            log_message(
                f"Auto-detected metadata from {detected_file} -> "
                f"reference={detected_reference or 'unknown'}, version={detected_fileformat or 'unknown'}"
            )
        log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}")

    except ValidationError as e:
        handle_critical_error(str(e), exc_cls=ValidationError)
    except FileNotFoundError as e:
        handle_critical_error(f"File not found: {e}", exc_cls=ValidationError)
    except PermissionError as e:
        handle_critical_error(f"Permission error: {e}", exc_cls=MergeVCFError)
    except OSError as e:
        handle_critical_error(f"Filesystem error: {e}", exc_cls=MergeVCFError)
    except Exception as e:
        handle_critical_error(f"Unexpected error: {e}", exc_cls=MergeVCFError)

    # Phase 2: validate, merge, annotate, finalize
    try:
        out_dir = Path(output_dir or (input_dir or "."))
        out_dir.mkdir(parents=True, exist_ok=True)
        log_message(f"Output directory: {out_dir}", verbose)

        try:
            valid_files, sample_order = validate_all_vcfs(
                input_dir,
                ref_genome,
                vcf_version,
                verbose=verbose,
                allow_gvcf=allow_gvcf,
            )
        except ValidationError:
            if not allow_gvcf:
                print("Use --allow-gvcf to include gVCF inputs.")
            raise

        qual_raw = getattr(args, "qual_threshold", None)
        an_raw = getattr(args, "an_threshold", None)
        qual_threshold = qual_raw if isinstance(qual_raw, (int, float)) and qual_raw >= 0 else None
        an_threshold = an_raw if isinstance(an_raw, (int, float)) and an_raw >= 0 else None
        allowed_filter_values = DEFAULT_ALLOWED_FILTER_VALUES

        merged_vcf_path = merge_vcfs(
            valid_files,
            str(out_dir),
            verbose=verbose,
            sample_order=sample_order,
            qual_threshold=qual_threshold,
            an_threshold=an_threshold,
            allowed_filter_values=allowed_filter_values,
        )

        # Combine file-provided header lines with CLI-provided metadata entries
        simple_header_lines = getattr(args, "simple_header_lines", None)
        header_metadata_lines = getattr(args, "header_metadata_lines", None)
        meta_entries = getattr(args, "meta", None)
        header_metadata_lines_combined: list[str] = []
        for source in (
            header_metadata_lines,
            simple_header_lines,
            meta_entries,
            extra_file_lines,
        ):
            if source:
                header_metadata_lines_combined.extend(source)
        header_metadata_lines = header_metadata_lines_combined

        sample_header_entries = getattr(args, "sample_header_entries", None)
        sample_metadata_entries = getattr(args, "sample_metadata_entries", None)
        sample_entries = sample_metadata_entries or sample_header_entries

        final_vcf_path = append_metadata_to_merged_vcf(
            merged_vcf_path,
            sample_metadata_entries=sample_entries,
            sample_header_entries=sample_header_entries,
            header_metadata_lines=header_metadata_lines,
            qual_threshold=qual_threshold,
            an_threshold=an_threshold,
            allowed_filter_values=allowed_filter_values,
            verbose=verbose,
        )

        validate_merged_vcf(final_vcf_path, verbose=verbose)

        final_target = out_dir / "cohort_final.vcf.gz"
        final_vcf_path = Path(final_vcf_path)

        if final_vcf_path.resolve() != final_target.resolve():
            if final_target.exists():
                final_target.unlink()

            if final_vcf_path.suffix != ".gz":
                try:
                    validation.pysam.tabix_compress(
                        str(final_vcf_path), str(final_vcf_path) + ".gz", force=True
                    )
                except Exception as e:
                    raise MergeVCFError(f"Failed to compress final VCF: {e}")
                try:
                    final_vcf_path.unlink()
                except Exception:
                    pass
                final_vcf_path = Path(str(final_vcf_path) + ".gz")

            try:
                final_vcf_path.replace(final_target)
            except Exception as e:
                raise MergeVCFError(f"Failed to rename final VCF: {e}")

            try:
                validation.pysam.tabix_index(str(final_target), preset="vcf", force=True)
            except Exception as e:
                raise MergeVCFError(f"Failed to index final VCF: {e}")

            final_vcf_path = final_target

        log_message(
            f"Script execution completed successfully. Final cohort VCF: {final_vcf_path}",
            verbose,
        )
        print(f"Wrote: {final_vcf_path} x 1.")
    except (ValidationError, MergeConflictError, MergeVCFError) as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)
