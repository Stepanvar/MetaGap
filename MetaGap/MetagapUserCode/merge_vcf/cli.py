"""Command-line entrypoint for the VCF merging workflow."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .logging_utils import MergeConflictError, MergeVCFError, ValidationError, log_message
from . import metadata as metadata_module
from . import validation
from . import merging

def _validate_metadata_entry(arg: str) -> str:
    """Validate that a metadata argument is in KEY=VALUE format and return a sanitized string."""
    if arg is None:
        raise argparse.ArgumentTypeError("Metadata entry cannot be empty")
    s = arg.strip()
    if s == "":
        # Allow empty strings (they will be filtered out later)
        return s
    if "=" not in s:
        raise argparse.ArgumentTypeError("Metadata entries must be in KEY=VALUE format (e.g., ID=VALUE)")
    key, value = s.split("=", 1)
    key = key.strip()
    value = value.strip()
    if key == "" or (key.lower() == "id" and value == ""):
        raise argparse.ArgumentTypeError("Metadata entries must have a non-empty ID")
    if key.startswith("-"):
        raise argparse.ArgumentTypeError("Metadata key cannot start with '-'")
    return f"{key}={value}"

def parse_arguments():
    """Parse and return command-line arguments for the VCF merging CLI."""
    parser = argparse.ArgumentParser(
        prog="merge_vcf",
        description="Consolidate multiple VCF files into a single merged VCF file."
    )
    # Input VCFs: allow multiple via --inputs or as positional arguments
    parser.add_argument(
        "input_vcfs",
        nargs="*",
        help="Input VCF files to merge (provide one or more file paths)."
    )
    parser.add_argument(
        "--inputs",
        nargs="+",
        dest="input_list",
        help="List of input VCF files to merge (alternative to providing them positionally)."
    )
    parser.add_argument(
        "-o", "--out-dir", "--output-dir", "--output",
        dest="output_dir",
        help="Directory for the merged output VCF. Defaults to current directory if omitted."
    )
    parser.add_argument(
        "--qual-threshold",
        type=float,
        default=30.0,
        help="Minimum QUAL value required to keep a variant (set negative to disable QUAL filtering)."
    )
    parser.add_argument(
        "--an-threshold",
        type=float,
        default=50.0,
        help="Minimum INFO/AN (allele count) required to keep a variant (set negative to disable AN filtering)."
    )
    parser.add_argument(
        "--sample-header",
        action="append",
        dest="sample_header_entries",
        type=_validate_metadata_entry,
        default=[],
        help="Sample metadata KEY=VALUE entry to include in the output VCF's SAMPLE header line (repeatable)."
    )
    parser.add_argument(
        "--simple-header",
        action="append",
        dest="simple_header_lines",
        type=_validate_metadata_entry,
        default=[],
        help="Additional header KEY=VALUE line to include in output VCF (repeatable, '##' will be prefixed automatically)."
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable verbose logging to console."
    )
    args = parser.parse_args()
    # Determine input file list (combine positional and --inputs)
    input_files = []
    if args.input_list:
        input_files.extend(args.input_list)
    if args.input_vcfs:
        input_files.extend(args.input_vcfs)
    if not input_files:
        parser.error("No input VCF files specified. Provide file paths or use --inputs.")
    # Remove any empty metadata entries that might have been added
    args.sample_header_entries = [s for s in args.sample_header_entries if s and s.strip()]
    args.simple_header_lines = [s if s.startswith("##") else f"##{s}" for s in args.simple_header_lines if s and s.strip()]
    # Ensure at least one sample ID is provided if any sample metadata given
    has_id = any(entry.split("=", 1)[0].strip().lower() == "id" and entry.split("=", 1)[1].strip() for entry in args.sample_header_entries)
    if not has_id and args.sample_header_entries:
        parser.error("A non-empty ID must be provided in sample header entries (e.g., --sample-header ID=SampleName).")
    # Normalize output directory
    args.output_dir = str(Path(args.output_dir)) if args.output_dir else "."
    # Finalize input file paths
    args.input_files = [str(Path(p)) for p in input_files]
    return args

def summarize_produced_vcfs(output_dir: str, fallback_vcf: str):
    """Return a tuple (directory, representative_filename, count) for output VCF shards in output_dir.
    If no per-sample VCF shards are found, returns the fallback_vcf's directory, filename, and count=1."""
    directory = Path(output_dir).resolve()
    if directory.is_dir():
        # Look for any per-sample VCF files in the output directory
        vcf_candidates = sorted(directory.glob("SAMPLE_*.vcf.gz")) + sorted(directory.glob("SAMPLE_*.vcf"))
        if vcf_candidates:
            # Pick the first (alphabetically) as representative
            representative_path = vcf_candidates[0].resolve()
            count = len(vcf_candidates)
            return str(representative_path.parent), representative_path.name, count
    # If no per-sample files found, fall back to the final VCF itself
    fallback_path = Path(fallback_vcf).resolve()
    return str(fallback_path.parent), fallback_path.name, 1

def main():
    """Main entry point for the VCF merging CLI. Parses arguments, runs the merging pipeline, and handles output."""
    args = parse_arguments()
    verbose = args.verbose
    try:
        out_dir = Path(args.output_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        log_message(f"Output directory: {out_dir}", verbose)
        # Validate and prepare input VCF files
        valid_files, sample_order = validation.validate_all_vcfs(args.input_files, verbose=verbose)
        # Perform merge and apply QUAL/AN filters
        merged_vcf_path = merging.merge_vcfs(
            valid_files,
            str(out_dir),
            verbose=verbose,
            sample_order=sample_order,
            qual_threshold=(args.qual_threshold if args.qual_threshold >= 0 else None),
            an_threshold=(args.an_threshold if args.an_threshold >= 0 else None)
        )
        # Append sample and simple header metadata (if provided) to the merged VCF
        final_vcf_path = metadata_module.append_metadata_to_merged_vcf(
            merged_vcf_path,
            sample_header_entries=args.sample_header_entries,
            simple_header_lines=args.simple_header_lines,
            verbose=verbose
        )
        # Validate the final merged VCF file structure and content
        validation.validate_merged_vcf(final_vcf_path, verbose=verbose)
        # Determine final output path and compress/index if needed
        final_target = out_dir / "cohort_final.vcf.gz"
        final_vcf_path = Path(final_vcf_path)
        if final_vcf_path.resolve() != final_target.resolve():
            if final_target.exists():
                final_target.unlink()
            # Compress final output if not already compressed
            if final_vcf_path.suffix != ".gz":
                try:
                    validation.pysam.tabix_compress(str(final_vcf_path), str(final_vcf_path) + ".gz", force=True)
                except Exception as e:
                    raise MergeVCFError(f"Failed to compress final VCF: {e}")
                try:
                    final_vcf_path.unlink()  # remove uncompressed file
                except Exception:
                    pass
                final_vcf_path = Path(str(final_vcf_path) + ".gz")
            # Rename to standard final file name
            try:
                final_vcf_path.replace(final_target)
            except Exception as e:
                raise MergeVCFError(f"Failed to rename final VCF: {e}")
            # Index the final VCF (create .tbi index for the bgzipped file)
            try:
                validation.pysam.tabix_index(str(final_target), preset="vcf", force=True)
            except Exception as e:
                raise MergeVCFError(f"Failed to index final VCF: {e}")
            final_vcf_path = final_target
        log_message(f"Script execution completed successfully. Final cohort VCF: {final_vcf_path}", verbose)
        print(f"Wrote: {final_vcf_path} x 1.")
    except (ValidationError, MergeConflictError, MergeVCFError) as exc:
        # Print error message and exit with code 1 on known exceptions
        print(f"ERROR: {exc}")
        sys.exit(1)
