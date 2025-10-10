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


# Re-export frequently patched helpers for the test suite.
validate_all_vcfs = validation.validate_all_vcfs
validate_merged_vcf = validation.validate_merged_vcf
merge_vcfs = merging.merge_vcfs
append_metadata_to_merged_vcf = metadata_module.append_metadata_to_merged_vcf

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

    def _metadata_list(values: list[str]) -> list[str]:
        """Return sanitized metadata entries, discarding blanks."""

        sanitized: list[str] = []
        for entry in values or []:
            if entry is None:
                continue
            cleaned = entry.strip()
            if not cleaned:
                continue
            sanitized.append(cleaned)
        return sanitized

    parser = argparse.ArgumentParser(
        prog="merge_vcf",
        description="Consolidate multiple VCF files into a single merged VCF.",
    )
    parser.add_argument(
        "--input-dir",
        dest="input_dir",
        required=True,
        help="Directory containing the VCF shards to merge.",
    )
    parser.add_argument(
        "-o",
        "--output",
        "--output-dir",
        dest="output_dir",
        help="Directory where merged VCF artifacts will be written.",
    )
    parser.add_argument("--ref", dest="ref", help="Reference genome build override.")
    parser.add_argument(
        "--vcf-version",
        dest="vcf_version",
        help="VCF specification version override (e.g., 4.2).",
    )
    parser.add_argument(
        "--allow-gvcf",
        dest="allow_gvcf",
        action="store_true",
        help="Permit gVCF inputs that contain <NON_REF> alleles.",
    )
    parser.add_argument(
        "--meta",
        dest="meta",
        action="append",
        default=[],
        type=_validate_metadata_entry,
        help="Generic KEY=VALUE metadata entries to add to the header.",
    )
    parser.add_argument(
        "--sample-metadata",
        dest="sample_metadata_entries",
        action="append",
        default=[],
        type=_validate_metadata_entry,
        help="SAMPLE metadata KEY=VALUE entries (repeatable).",
    )
    parser.add_argument(
        "--header-metadata",
        dest="header_metadata_lines",
        action="append",
        default=[],
        type=_validate_metadata_entry,
        help="Header KEY=VALUE entries (repeatable).",
    )
    parser.add_argument(
        "--metadata-template",
        dest="metadata_template_path",
        help="Optional metadata template file to seed header entries.",
    )
    parser.add_argument(
        "--metadata-file",
        dest="metadata_file",
        help="Optional file with additional header metadata lines.",
    )
    parser.add_argument(
        "--qual-threshold",
        type=float,
        default=30.0,
        help="Minimum QUAL to keep a variant. <0 disables.",
    )
    parser.add_argument(
        "--an-threshold",
        type=float,
        default=50.0,
        help="Minimum INFO/AN to keep a variant. <0 disables.",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose console logging.")

    args = parser.parse_args()

    try:
        input_dir_path = Path(args.input_dir).expanduser().resolve()
    except OSError as exc:  # pragma: no cover - defensive
        parser.error(f"Unable to normalize input directory: {exc}")
        raise

    if not input_dir_path.exists() or not input_dir_path.is_dir():
        parser.error(f"Input directory does not exist: {args.input_dir}")

    args.input_dir = str(input_dir_path)

    output_dir_value = args.output_dir
    if output_dir_value:
        try:
            output_dir_path = Path(output_dir_value).expanduser().resolve()
        except OSError as exc:  # pragma: no cover - defensive
            parser.error(f"Unable to normalize output directory: {exc}")
            raise
    else:
        output_dir_path = input_dir_path
    args.output_dir = str(output_dir_path)

    if args.metadata_template_path:
        try:
            template_path = Path(args.metadata_template_path).expanduser().resolve()
        except OSError as exc:  # pragma: no cover - defensive
            parser.error(f"Unable to normalize metadata template path: {exc}")
            raise
        args.metadata_template_path = str(template_path)

    if args.metadata_file:
        try:
            metadata_file_path = Path(args.metadata_file).expanduser().resolve()
        except OSError as exc:  # pragma: no cover - defensive
            parser.error(f"Unable to normalize metadata file path: {exc}")
            raise
        args.metadata_file = str(metadata_file_path)

    args.meta = _metadata_list(args.meta)
    args.sample_metadata_entries = _metadata_list(args.sample_metadata_entries)
    args.header_metadata_lines = _metadata_list(args.header_metadata_lines)

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
    verbose = args.verbose
    ref_genome = None
    vcf_version = None

    # Phase 1: environment + autodetect metadata
    try:
        input_dir = getattr(args, "input_dir", None)
        if not input_dir:
            raise ValidationError("An input directory must be provided.")
        if not os.path.isdir(input_dir):
            raise ValidationError(f"Input directory does not exist: {input_dir}")

        output_dir = getattr(args, "output_dir", None) or input_dir
        os.makedirs(output_dir, exist_ok=True)

        log_path = os.path.join(output_dir, LOG_FILE)
        configure_logging(
            log_level=logging.DEBUG if verbose else logging.INFO,
            log_file=log_path,
            enable_file_logging=True,
            enable_console=verbose,
        )

        metadata_file = getattr(args, "metadata_file", None)
        extra_file_lines = (
            metadata_module.load_metadata_lines(metadata_file, verbose) if metadata_file else []
        )

        log_message("Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        log_message(f"Input directory: {input_dir}")
        log_message(f"Output directory: {output_dir}")

        detected_file, detected_fileformat, detected_reference = validation.find_first_vcf_with_header(
            input_dir, verbose
        )
        ref_override = getattr(args, "ref", None)
        vcf_override = getattr(args, "vcf_version", None)
        ref_genome = ref_override or detected_reference
        vcf_version = vcf_override or validation.normalize_vcf_version(detected_fileformat)

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
        out_dir = Path(getattr(args, "output_dir", input_dir))
        out_dir.mkdir(parents=True, exist_ok=True)
        log_message(f"Output directory: {out_dir}", verbose)

        valid_files, sample_order = validate_all_vcfs(
            input_dir,
            ref_genome,
            vcf_version,
            verbose=verbose,
            allow_gvcf=getattr(args, "allow_gvcf", False),
        )

        qual_threshold = getattr(args, "qual_threshold", 30.0)
        qual_threshold = qual_threshold if qual_threshold >= 0 else None
        an_threshold = getattr(args, "an_threshold", 50.0)
        an_threshold = an_threshold if an_threshold >= 0 else None
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

        # Combine file-provided header lines with CLI-provided simple headers
        header_metadata_lines = (
            (getattr(args, "header_metadata_lines", []) or []) + (extra_file_lines or [])
        )

        final_vcf_path = append_metadata_to_merged_vcf(
            merged_vcf_path,
            sample_metadata_entries=getattr(args, "sample_metadata_entries", []),
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
                    validation.pysam.tabix_compress(str(final_vcf_path), str(final_vcf_path) + ".gz", force=True)
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

        log_message(f"Script execution completed successfully. Final cohort VCF: {final_vcf_path}", verbose)
        print(f"Wrote: {final_vcf_path} x 1.")
    except (ValidationError, MergeConflictError, MergeVCFError) as exc:
        print(f"ERROR: {exc}")
        sys.exit(1)
