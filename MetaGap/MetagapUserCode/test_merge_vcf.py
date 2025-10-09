#!/usr/bin/env python3
"""Public CLI surface for the MetaGap VCF merging workflow."""

import datetime
import os
import subprocess
import sys

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from merge_vcf import VCFPY_AVAILABLE, vcfpy
from merge_vcf import metadata as _metadata
from merge_vcf import cli as _cli
from merge_vcf.logging_utils import (
    LOG_FILE,
    handle_critical_error,
    handle_non_critical_error,
    log_message,
    logger,
)
from merge_vcf.metadata import (
    _format_sample_metadata_value,
    _parse_sample_metadata_line,
    _parse_simple_metadata_line,
    _validate_anonymized_vcf_header as _metadata_validate_anonymized_vcf_header,
    append_metadata_to_merged_vcf,
    apply_metadata_to_header,
    build_sample_metadata_line,
    parse_metadata_arguments,
    recalculate_cohort_info_tags,
)
from merge_vcf.merging import (
    _create_missing_call_factory,
    _pad_record_samples,
    _remove_format_and_sample_definitions,
    merge_vcfs,
    preprocess_vcf,
    union_headers,
)
from merge_vcf.validation import (
    find_first_vcf_with_header,
    is_gvcf_header,
    normalize_vcf_version,
    validate_all_vcfs,
    validate_merged_vcf,
    validate_vcf,
)


def parse_arguments():
    """Proxy to the CLI parser that remains patchable by tests."""

    return _cli.parse_arguments()


def summarize_produced_vcfs(output_dir, fallback_vcf):
    return _cli.summarize_produced_vcfs(output_dir, fallback_vcf)


def _validate_anonymized_vcf_header(
    final_vcf_path, ensure_for_uncompressed: bool = False
):
    original = _metadata.handle_critical_error
    _metadata.handle_critical_error = handle_critical_error
    try:
        return _metadata_validate_anonymized_vcf_header(
            final_vcf_path, ensure_for_uncompressed=ensure_for_uncompressed
        )
    finally:
        _metadata.handle_critical_error = original


def main():
    args = parse_arguments()
    verbose = args.verbose
    (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries,
        sanitized_header_lines,
        serialized_sample_line,
    ) = parse_metadata_arguments(args, verbose)

    log_message(
        "Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        verbose,
    )

    input_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(input_dir):
        handle_critical_error(f"Input directory does not exist: {input_dir}")
    log_message("Input directory: " + input_dir, verbose)

    output_dir = os.path.abspath(args.output_dir) if args.output_dir else input_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    log_message("Output directory: " + output_dir, verbose)

    ref_genome = args.ref.strip() if args.ref else None
    vcf_version = normalize_vcf_version(args.vcf_version) if args.vcf_version else None

    detection_needed = not ref_genome or not vcf_version
    detected_file = detected_fileformat = detected_reference = None
    if detection_needed:
        detected_file, detected_fileformat, detected_reference = find_first_vcf_with_header(
            input_dir, verbose
        )
        if not ref_genome and detected_reference:
            ref_genome = detected_reference
        if not vcf_version and detected_fileformat:
            vcf_version = normalize_vcf_version(detected_fileformat)

    if not ref_genome:
        handle_critical_error(
            "Reference genome build must be provided via --ref or auto-detectable from input files."
        )
    if not vcf_version:
        handle_critical_error(
            "VCF version must be provided via --vcf-version or auto-detectable from input files."
        )

    if detection_needed and detected_file:
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
        allow_gvcf=getattr(args, "allow_gvcf", False),
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
        header_metadata_lines=sanitized_header_lines,
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

__all__ = [
    "VCFPY_AVAILABLE",
    "vcfpy",
    "LOG_FILE",
    "logger",
    "log_message",
    "handle_critical_error",
    "handle_non_critical_error",
    "subprocess",
    "is_gvcf_header",
    "parse_arguments",
    "_format_sample_metadata_value",
    "build_sample_metadata_line",
    "_parse_sample_metadata_line",
    "parse_metadata_arguments",
    "apply_metadata_to_header",
    "_validate_anonymized_vcf_header",
    "append_metadata_to_merged_vcf",
    "recalculate_cohort_info_tags",
    "normalize_vcf_version",
    "find_first_vcf_with_header",
    "validate_vcf",
    "validate_all_vcfs",
    "merge_vcfs",
    "union_headers",
    "preprocess_vcf",
    "_create_missing_call_factory",
    "_remove_format_and_sample_definitions",
    "_pad_record_samples",
    "_parse_simple_metadata_line",
    "validate_merged_vcf",
    "summarize_produced_vcfs",
    "main",
]


if __name__ == "__main__":
    main()
