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

# Re-export frequently patched helpers for easier test monkeypatching.
merge_vcfs = merging.merge_vcfs
append_metadata_to_merged_vcf = metadata_module.append_metadata_to_merged_vcf
validate_all_vcfs = validation.validate_all_vcfs
validate_merged_vcf = validation.validate_merged_vcf
from .logging_utils import (
    LOG_FILE,
    MergeConflictError,
    MergeVCFError,
    ValidationError,
    configure_logging,
    handle_critical_error,
    log_message,
)


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
        "--out-dir",
        "--output-dir",
        "--output",
        dest="output_dir",
        help="Output directory for the merged VCF. Defaults to current directory.",
    )
    parser.add_argument("--ref", dest="ref", help="Reference genome build to use.")
    parser.add_argument("--vcf-version", dest="vcf_version", help="VCF version for the merged output.")
    parser.add_argument(
        "--allow-gvcf",
        dest="allow_gvcf",
        action="store_true",
        help="Allow gVCF inputs when validating files.",
    )
    parser.add_argument("--qual-threshold", type=float, default=30.0, help="Minimum QUAL to keep a variant. <0 disables.")
    parser.add_argument("--an-threshold", type=float, default=50.0, help="Minimum INFO/AN to keep a variant. <0 disables.")
    parser.add_argument(
        "--sample-metadata",
        action="append",
        dest="sample_metadata_entries",
        type=_validate_metadata_entry,
        default=[],
        help="SAMPLE header KEY=VALUE entry (repeatable). Must include ID=...",
    )
    parser.add_argument(
        "--header-metadata",
        action="append",
        dest="header_metadata_lines",
        type=_validate_metadata_entry,
        default=[],
        help="Additional header KEY=VALUE (repeatable). '##' is auto-prefixed.",
    )
    parser.add_argument(
        "-m",
        "--metadata-template",
        dest="metadata_template_path",
        help="Optional file with extra header lines. One entry per line. Lines without '##' will be prefixed.",
    )
    parser.add_argument(
        "--meta",
        action="append",
        dest="meta_entries",
        type=_validate_metadata_entry,
        default=[],
        help="Generic metadata KEY=VALUE entry. ID keys route to sample metadata; others to header metadata.",
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

    input_files: list[str] = []
    if args.input_dir:
        input_dir_path = Path(args.input_dir)
        if not input_dir_path.is_dir():
            parser.error(f"Input directory does not exist: {args.input_dir}")
        globbed_files = list(input_dir_path.glob("*.vcf")) + list(input_dir_path.glob("*.vcf.gz"))
        input_files = [str(p) for p in sorted(globbed_files)]
    else:
        input_files = explicit_inputs

    if not input_files:
        parser.error("No input VCF files specified.")

# Clean metadata (conflict-resolved)
	def _norm(s: str) -> str:
		return s.strip()

	def _ensure_hashes(s: str) -> str:
		return s if s.startswith("##") else f"##{s}"

	def _route_meta(entry: str) -> tuple[str, str]:
		# SAMPLE=... goes to sample metadata; everything else goes to header metadata
		body = entry.lstrip("#").strip()
		if body.upper().startswith("SAMPLE="):
			return "sample", _ensure_hashes(body)
		return "header", _ensure_hashes(body)

	def _dedupe(seq: list[str]) -> list[str]:
		seen = set()
		out = []
		for x in seq:
			if x not in seen:
				seen.add(x)
				out.append(x)
		return out

	sample_src = (getattr(args, "sample_header_entries", []) or []) + \
				 (getattr(args, "sample_metadata_entries", []) or [])
	header_src = (getattr(args, "simple_header_lines", []) or []) + \
				 (getattr(args, "header_metadata_lines", []) or [])
	meta_src   = (getattr(args, "meta_entries", []) or [])

	sample_list: list[str] = []
	header_list: list[str] = []

	for s in sample_src:
		s = _norm(s)
		if s:
			sample_list.append(_ensure_hashes(s))

	for h in header_src:
		h = _norm(h)
		if h:
			header_list.append(_ensure_hashes(h))

	for m in meta_src:
		m = _norm(m)
		if not m:
			continue
		which, v = _route_meta(m)
		(sample_list if which == "sample" else header_list).append(v)

	# drop fileformat lines; dedupe; set canonical fields
	args.sample_metadata_entries = _dedupe([x for x in sample_list if not x.startswith("##fileformat")])
	args.header_metadata_lines   = _dedupe([x for x in header_list if not x.startswith("##fileformat")])

	# optional: clear legacy fields to avoid later confusion
	args.sample_header_entries = []
	args.simple_header_lines = []
	args.meta_entries = []

    # Require SAMPLE ID if any SAMPLE metadata present
    has_id = any(
        (kv := ent.split("=", 1)) and kv[0].strip().lower() == "id" and kv[1].strip()
        for ent in args.sample_metadata_entries
    )
    if args.sample_metadata_entries and not has_id:
        parser.error("Provide a non-empty SAMPLE ID (e.g., --sample-metadata ID=SampleName).")

    # Normalize paths
    args.output_dir = str(Path(args.output_dir)) if args.output_dir else "."
    args.input_dir = str(Path(args.input_dir)) if args.input_dir else None
    args.metadata_template_path = (
        str(Path(args.metadata_template_path)) if args.metadata_template_path else None
    )
    args.input_files = [str(Path(p)) for p in input_files]
    args.meta_entries = meta_entries
    return args


def _normalize_argument_namespace(args):
    """Ensure the parsed args namespace exposes a consistent attribute surface."""

    # Normalize sample metadata entries, accepting legacy attribute names.
    sample_entries: list[str] = []
    for attr in ("sample_header_entries", "sample_metadata_entries", "meta"):
        values = getattr(args, attr, None)
        if values:
            sample_entries.extend(list(values))
    args.sample_header_entries = [s for s in sample_entries if s]

    sample_mapping: dict[str, str] = {}
    for entry in args.sample_header_entries:
        if "=" not in entry:
            continue
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key:
            sample_mapping[key] = value
    args.sample_metadata_entries = sample_mapping

    # Normalize simple header metadata lines.
    simple_lines: list[str] = list(getattr(args, "simple_header_lines", []) or [])
    fallback_lines = getattr(args, "header_metadata_lines", None)
    if fallback_lines:
        for line in fallback_lines:
            if not line:
                continue
            simple_lines.append(line if line.startswith("##") else f"##{line}")
    args.simple_header_lines = [s for s in simple_lines if s]

    # Prefer explicit metadata file/template if provided via alternate argument name.
    metadata_file = getattr(args, "metadata_file", None)
    if metadata_file is None:
        metadata_file = getattr(args, "metadata_template_path", None)
    args.metadata_file = metadata_file

    # Ensure commonly accessed attributes exist with sensible defaults.
    if not hasattr(args, "allow_gvcf"):
        args.allow_gvcf = False
    if not hasattr(args, "verbose"):
        args.verbose = False
    if not hasattr(args, "output_dir"):
        args.output_dir = getattr(args, "out_dir", ".")
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

    # Phase 1: environment + autodetect metadata
    try:
        # Infer a base input directory from provided files or explicit directory.
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

        args.input_files = [str(Path(p)) for p in input_files]

        output_dir = os.path.abspath(args.output_dir) if args.output_dir else input_dir
        os.makedirs(output_dir, exist_ok=True)

        log_path = os.path.join(output_dir, LOG_FILE)
        configure_logging(
            log_level=logging.DEBUG if verbose else logging.INFO,
            log_file=log_path,
            enable_file_logging=True,
            enable_console=verbose,
        )

        extra_file_lines = (
            metadata_module.load_metadata_lines(args.metadata_template_path, verbose)
            if args.metadata_template_path
            else []
        )

        log_message("Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        log_message(f"Input directory: {input_dir}")
        log_message(f"Output directory: {output_dir}")

        detected_file, detected_fileformat, detected_reference = validation.find_first_vcf_with_header(input_dir, verbose)
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
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        log_message(f"Output directory: {out_dir}", verbose)

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
            message = str(exc)
            if (
                not getattr(args, "allow_gvcf", False)
                and "No valid VCF files remain after validation" in message
            ):
                message = message.rstrip(".") + ". Use --allow-gvcf to include gVCF inputs."
            raise ValidationError(message) from exc

        qual_threshold = args.qual_threshold if args.qual_threshold >= 0 else None
        an_threshold = args.an_threshold if args.an_threshold >= 0 else None
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
        cli_header_metadata = getattr(args, "header_metadata_lines", None)
        header_metadata_lines = (cli_header_metadata or []) + (extra_file_lines or [])

        sample_metadata_entries = getattr(args, "sample_metadata_entries", {}) or {}
        serialized_sample_line = None
        if sample_metadata_entries:
            try:
                serialized_sample_line = metadata_module.build_sample_metadata_line(
                    sample_metadata_entries
                )
            except ValueError as exc:
                raise ValidationError(str(exc)) from exc

        final_vcf_path = append_metadata_to_merged_vcf(
            merged_vcf_path,
            sample_metadata_entries=sample_metadata_entries,
            header_metadata_lines=header_metadata_lines,
            serialized_sample_line=serialized_sample_line,
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


if __name__ == "__main__":
    if __package__ in {None, ""}:
        parent_dir = Path(__file__).resolve().parent.parent
        if str(parent_dir) not in sys.path:
            sys.path.insert(0, str(parent_dir))
        __package__ = "merge_vcf"
    main()
