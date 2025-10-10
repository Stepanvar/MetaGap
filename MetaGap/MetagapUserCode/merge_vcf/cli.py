"""Command-line entrypoint for the VCF merging workflow.

The CLI can be executed via ``python -m merge_vcf.cli`` or
``python merge_vcf/cli.py``.
"""
from __future__ import annotations

import argparse
import datetime
import logging
import os
import sys
from pathlib import Path

if __package__ in {None, ""} and __name__ == "__main__":
    # Ensure relative imports succeed when executed as a script.
    parent_dir = Path(__file__).resolve().parent.parent
    if str(parent_dir) not in sys.path:
        sys.path.insert(0, str(parent_dir))
    __package__ = "merge_vcf"

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

# Re-export helpers for tests/patching.
merge_vcfs = merging.merge_vcfs
append_metadata_to_merged_vcf = metadata_module.append_metadata_to_merged_vcf
validate_all_vcfs = validation.validate_all_vcfs
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

    def _metadata_list(values: list[str] | None) -> list[str]:
        """Return sanitized metadata entries, discarding blanks."""
        sanitized: list[str] = []
        for entry in values or []:
            cleaned = entry.strip() if entry is not None else ""
            if cleaned:
                sanitized.append(cleaned)
        return sanitized

    parser = argparse.ArgumentParser(
        prog="merge_vcf",
        description="Consolidate multiple VCF files into a single merged VCF.",
    )
    parser.add_argument("input_vcfs", nargs="*", help="Input VCF files.")

    input_group = parser.add_mutually_exclusive_group()
    input_group.add_argument(
        "-i",
        "--input-dir",
        dest="input_dir",
        help="Directory containing input VCF files. All *.vcf and *.vcf.gz files will be merged.",
    )
    input_group.add_argument(
        "--inputs",
        nargs="+",
        dest="input_list",
        help="Alternate way to pass input VCFs.",
    )
    parser.add_argument(
        "-o",
        "--output",
        "--output-dir",
        dest="output_dir",
        help="Directory where merged VCF artifacts will be written.",
    )
    parser.add_argument(
        "-m",
        "--metadata-template",
        dest="metadata_template_path",
        help="Optional file with extra header lines. One entry per line. Lines without '##' will be prefixed.",
    )
    parser.add_argument(
        "--allow-gvcf",
        action="store_true",
        help="Allow gVCF inputs during validation.",
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

    # Collect inputs
    explicit_inputs: list[str] = []
    if args.input_list:
        explicit_inputs.extend(args.input_list)
    if args.input_vcfs:
        explicit_inputs.extend(args.input_vcfs)

    if args.input_dir and explicit_inputs:
        parser.error("Provide either --input-dir or explicit VCF paths, not both.")
    if not args.input_dir and not explicit_inputs:
        parser.error("No input VCF files specified. Provide paths or --input-dir.")

    # Resolve input files
    input_files: list[str] = []
    if args.input_dir:
        try:
            input_dir_path = Path(args.input_dir).expanduser().resolve()
        except OSError as exc:  # defensive
            parser.error(f"Unable to normalize input directory: {exc}")
            raise
        if not input_dir_path.is_dir():
            parser.error(f"Input directory does not exist: {args.input_dir}")
        globbed_files = list(input_dir_path.glob("*.vcf")) + list(input_dir_path.glob("*.vcf.gz"))
        input_files = [str(p) for p in sorted(globbed_files)]
        args.input_dir = str(input_dir_path)
    else:
        input_files = [str(Path(p).expanduser().resolve()) for p in explicit_inputs]
        # Derive a common directory for logging if possible
        base_dirs = [os.path.dirname(p) or "." for p in input_files]
        if base_dirs:
            try:
                args.input_dir = os.path.commonpath(base_dirs)
            except ValueError:
                args.input_dir = os.path.dirname(input_files[0]) if input_files else "."

    if not input_files:
        parser.error("No input VCF files specified.")

    # Output directory
    output_dir_value = args.output_dir or args.input_dir
    try:
        output_dir_path = Path(output_dir_value).expanduser().resolve()
    except OSError as exc:  # defensive
        parser.error(f"Unable to normalize output directory: {exc}")
        raise
    args.output_dir = str(output_dir_path)

    # Metadata template path
    if args.metadata_template_path:
        try:
            template_path = Path(args.metadata_template_path).expanduser().resolve()
        except OSError as exc:  # defensive
            parser.error(f"Unable to normalize metadata template path: {exc}")
            raise
        args.metadata_template_path = str(template_path)

    # Sanitize lists
    args.header_metadata_lines = _metadata_list(getattr(args, "header_metadata_lines", None))
    args.sample_metadata_entries = _metadata_list(getattr(args, "sample_metadata_entries", None))

    # Finalize input list
    args.input_files = [str(Path(p)) for p in input_files]
    return args


def _normalize_argument_namespace(args):
    """Ensure the parsed args namespace exposes a consistent attribute surface."""

    # Sample metadata: convert list of KEY=VALUE to dict
    sample_mapping: dict[str, str] = {}
    for entry in getattr(args, "sample_metadata_entries", []) or []:
        if "=" not in entry:
            continue
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key:
            sample_mapping[key] = value
    args.sample_metadata_entries = sample_mapping

    # Normalize header metadata lines, ensure '##' prefix, drop empties and fileformat
    normalized_headers: list[str] = []
    for line in getattr(args, "header_metadata_lines", []) or []:
        s = line.strip()
        if not s:
            continue
        if not s.startswith("##"):
            s = f"##{s}"
        if s.lower().startswith("##fileformat"):
            continue
        normalized_headers.append(s)
    # de-duplicate preserving order
    seen = set()
    args.header_metadata_lines = [x for x in normalized_headers if not (x in seen or seen.add(x))]

    # Unify metadata file arg name
    metadata_file = getattr(args, "metadata_file", None)
    if metadata_file is None:
        metadata_file = getattr(args, "metadata_template_path", None)
    args.metadata_file = metadata_file

    # Ensure defaults
    if not hasattr(args, "allow_gvcf"):
        args.allow_gvcf = False
    if not hasattr(args, "verbose"):
        args.verbose = False
    if not hasattr(args, "qual_threshold"):
        args.qual_threshold = 30.0
    if not hasattr(args, "an_threshold"):
        args.an_threshold = 50.0

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
    args = _normalize_argument_namespace(parse_arguments())
    verbose = args.verbose
    ref_genome = None
    vcf_version = None

    # Phase 1: environment + autodetect metadata
    try:
        input_files = list(getattr(args, "input_files", []) or [])
        explicit_input_dir = getattr(args, "input_dir", None)

        if input_files:
            base_dirs = [os.path.dirname(os.path.abspath(p)) or "." for p in input_files]
            input_dir = os.path.commonpath(base_dirs) if base_dirs else "."
        elif explicit_input_dir:
            input_dir = os.path.abspath(explicit_input_dir)
            if not os.path.isdir(input_dir):
                raise ValidationError(f"Input directory does not exist: {input_dir}")
            input_files = [
                str(path)
                for path in sorted(Path(input_dir).glob("*.vcf"))
                + sorted(Path(input_dir).glob("*.vcf.gz"))
            ]
        else:
            raise ValidationError("No input VCF files or input directory specified.")

        if not os.path.isdir(input_dir):
            raise ValidationError(f"Input directory does not exist: {input_dir}")

        # Output directory
        output_dir_arg = getattr(args, "output_dir", None)
        output_dir = os.path.abspath(output_dir_arg) if output_dir_arg else input_dir
        if not output_dir:
            raise ValueError("output_dir is not set and input_dir is empty")
        os.makedirs(output_dir, exist_ok=True)

        # Logging
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

        # Load extra header lines from template file
        metadata_file = getattr(args, "metadata_file", None)
        extra_file_lines = (
            metadata_module.load_metadata_lines(metadata_file, verbose) if metadata_file else []
        )

        log_message("Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        log_message(f"Input directory: {input_dir}")
        log_message(f"Output directory: {output_dir}")

        # Detect reference and VCF version from first readable VCF
        detected_file = None
        detected_fileformat = None
        detected_reference = None

        (
            detected_file,
            detected_fileformat,
            detected_reference,
        ) = validation.find_first_vcf_with_header(input_dir, verbose)

        ref_genome = detected_reference
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

        # Validate inputs
        args_input_dir = getattr(args, "input_dir", None) or input_dir
        setattr(args, "input_dir", args_input_dir)

        try:
            valid_files, sample_order = validate_all_vcfs(
                args_input_dir,
                ref_genome,
                vcf_version,
                verbose=verbose,
                allow_gvcf=getattr(args, "allow_gvcf", False),
            )
        except ValidationError as exc:
            msg = str(exc)
            if (
                not getattr(args, "allow_gvcf", False)
                and "No valid VCF files remain after validation" in msg
            ):
                msg = msg.rstrip(".") + ". Use --allow-gvcf to include gVCF inputs."
            raise ValidationError(msg) from exc

        def _coerce_threshold(val):
            if isinstance(val, (int, float)) and val >= 0:
                return float(val)
            return None

        qual_threshold = _coerce_threshold(getattr(args, "qual_threshold", None))
        an_threshold = _coerce_threshold(getattr(args, "an_threshold", None))
        allowed_filter_values = DEFAULT_ALLOWED_FILTER_VALUES

        # Merge
        merged_vcf_path = merge_vcfs(
            valid_files,
            str(out_dir),
            verbose=verbose,
            sample_order=sample_order,
            qual_threshold=qual_threshold,
            an_threshold=an_threshold,
            allowed_filter_values=allowed_filter_values,
        )

        # Build header metadata (combine CLI + file), sanitize, de-dup
        def _norm(s: str) -> str: return s.strip()
        def _ensure_hashes(s: str) -> str: return s if s.startswith("##") else f"##{s}"

        raw_header_sources = [
            getattr(args, "header_metadata_lines", None),
            extra_file_lines,
        ]
        header_metadata_lines: list[str] = []
        for src in raw_header_sources:
            if not src:
                continue
            for s in src:
                s = _norm(s)
                if not s or s.lower().startswith("##fileformat"):
                    continue
                header_metadata_lines.append(_ensure_hashes(s))
        seen = set()
        header_metadata_lines = [x for x in header_metadata_lines if not (x in seen or seen.add(x))]

        # Sample metadata: build SAMPLE=<...> if provided as dict
        sample_metadata_entries = getattr(args, "sample_metadata_entries", None)
        serialized_sample_line = None
        if isinstance(sample_metadata_entries, dict) and sample_metadata_entries:
            try:
                serialized_sample_line = metadata_module.build_sample_metadata_line(
                    sample_metadata_entries
                )
            except ValueError as exc:
                raise ValidationError(str(exc)) from exc

        # Append metadata
        final_vcf_path = append_metadata_to_merged_vcf(
            merged_vcf_path,
            sample_metadata_entries=sample_metadata_entries,
            header_metadata_lines=header_metadata_lines,
            serialized_sample_line=serialized_sample_line,
            qual_threshold=qual_threshold,
            an_threshold=an_threshold,
            allowed_filter_values=allowed_filter_values,
            verbose=verbose,
        )

        # Validate and place final artifact at standard path
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


if __name__ == "__main__":
    if __package__ in {None, ""}:
        parent_dir = Path(__file__).resolve().parent.parent
        if str(parent_dir) not in sys.path:
            sys.path.insert(0, str(parent_dir))
        __package__ = "merge_vcf"
    main()
