#!/usr/bin/env python3
"""Backward compatible wrapper for the :mod:`merge_vcf` package."""

from __future__ import annotations

import sys
import glob
import argparse
import csv
import datetime
import re
import copy
import tempfile
import shutil
import subprocess
import gzip
from collections import OrderedDict

MODULE_DIR = os.path.dirname(os.path.abspath(__file__))
if MODULE_DIR not in sys.path:
    sys.path.insert(0, MODULE_DIR)

from merging import (
    MergeCallbacks,
    append_metadata_to_merged_vcf,
    configure_callbacks as configure_merging_callbacks,
    merge_vcfs,
    union_headers,
    _create_missing_call_factory,
    _pad_record_samples,
    _remove_format_and_sample_definitions,
)

try:
    import vcfpy  # type: ignore
    VCFPY_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised when dependency is missing
    VCFPY_AVAILABLE = False

    class _MissingVcfpyModule:
        """Placeholder that raises a helpful error when vcfpy is unavailable."""

        __slots__ = ()

        def __getattr__(self, name):  # pragma: no cover - defensive, attribute driven
            raise ModuleNotFoundError(
                "Error: vcfpy package is required. Please install it with 'pip install vcfpy'."
            )

    vcfpy = _MissingVcfpyModule()  # type: ignore

# Configure logging
LOG_FILE = "script_execution.log"
logger = logging.getLogger("vcf_merger")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
fh = logging.FileHandler(LOG_FILE)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)


def log_message(message, verbose=False):
    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message):
    log_message("CRITICAL ERROR: " + message)
    print("A critical error occurred. Check {} for details.".format(LOG_FILE))
    sys.exit(1)


def handle_non_critical_error(message):
    log_message("WARNING: " + message)
    print("Warning: " + message)


def is_gvcf_header(header_lines):
    """Return True if the provided header metadata indicates a gVCF file."""

    for raw_line in header_lines:
        line_text = None
        if hasattr(raw_line, "to_line"):
            try:
                line_text = raw_line.to_line()
            except Exception:
                line_text = None
        if line_text is None and isinstance(raw_line, str):
            line_text = raw_line

        if line_text:
            stripped = line_text.strip()
            if stripped.startswith("##"):
                normalized = stripped.upper().replace(" ", "")
                if normalized.startswith("##GVCF"):
                    return True
                if normalized.startswith("##ALT=") and (
                    "ID=<NON_REF>" in normalized or "ID=NON_REF" in normalized
                ):
                    return True

        key = getattr(raw_line, "key", None)
        if isinstance(key, str):
            upper_key = key.upper()
            if upper_key == "GVCF":
                return True
            if upper_key == "ALT":
                candidates = [
                    getattr(raw_line, "id", None),
                    getattr(raw_line, "value", None),
                ]
                mapping = getattr(raw_line, "mapping", None)
                if isinstance(mapping, dict):
                    candidates.extend(mapping.values())
                for candidate in candidates:
                    if not isinstance(candidate, str):
                        continue
                    normalized_candidate = candidate.upper().replace(" ", "")
                    if "<NON_REF>" in normalized_candidate or normalized_candidate == "NON_REF":
                        return True
    return False


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
        "--ref",
        help="Expected reference genome build. When omitted the script attempts to auto-detect it.",
    )
    parser.add_argument(
        "--vcf-version",
        help="Expected VCF version (e.g., 4.2). When omitted the script attempts to auto-detect it.",
    )
    parser.add_argument(
        "--meta",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Sample metadata in KEY=VALUE form. Repeat for multiple keys. "
            "An ID entry is required when metadata is provided."
        ),
    )
    parser.add_argument(
        "--allow-gvcf",
        action="store_true",
        help="Allow input files that contain gVCF annotations such as <NON_REF> ALT alleles.",
    )
    parser.add_argument(
        "--sample-metadata",
        dest="sample_metadata_entries",
        action="append",
        metavar="KEY=VALUE",
        default=[],
        help="Add a key/value pair to the ##SAMPLE metadata line. Provide at least ID=... to emit the line.",
    )
    parser.add_argument(
        "--header-metadata",
        dest="header_metadata_lines",
        action="append",
        metavar="LINE",
        default=[],
        help="Add an arbitrary metadata header line (##key=value). The '##' prefix is optional.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose console logging in addition to the log file.",
    )
    return parser.parse_args()


def _validate_anonymized_vcf_header(final_vcf_path, ensure_for_uncompressed=False):
    """Use ``bcftools`` to confirm that ``final_vcf_path`` has a readable header."""

    if not final_vcf_path:
        return

    requires_validation = ensure_for_uncompressed or final_vcf_path.endswith(".vcf.gz")
    if not requires_validation:
        return

    try:
        subprocess.run(
            ["bcftools", "view", "-h", final_vcf_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (subprocess.CalledProcessError, OSError) as exc:
        handle_critical_error(
            f"Failed to read header from anonymized VCF using bcftools: {exc}"
        )




def normalize_vcf_version(version_value):
    """Extract the numeric portion of the VCF version if present."""
    if not version_value:
        return version_value
    match = re.search(r"(\d+\.\d+)", version_value)
    if match:
        return match.group(1)
    return version_value.strip()


def find_first_vcf_with_header(input_dir, verbose=False):
    """Locate the first readable VCF file and return its header metadata."""
    candidates = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    for file_path in candidates:
        try:
            reader = vcfpy.Reader.from_path(file_path)
        except Exception as e:
            handle_non_critical_error(
                f"Unable to open {file_path} for auto-detection: {str(e)}. Trying next file.",
            )
            continue

        header = reader.header
        fileformat = None
        reference = None
        for line in header.lines:
            if line.key == "fileformat" and fileformat is None:
                fileformat = line.value
            elif line.key == "reference" and reference is None:
                reference = line.value
            if fileformat and reference:
                break

        if fileformat and reference:
            log_message(
                f"Auto-detection using {file_path}: fileformat={fileformat}, reference={reference}",
                verbose,
            )
            return file_path, fileformat, reference

        handle_non_critical_error(
            f"Required metadata missing in {file_path} during auto-detection. Trying next file.",
        )

    return None, None, None


def validate_vcf(file_path, ref_genome, vcf_version, verbose=False, allow_gvcf=False):
    log_message(f"Validating file: {file_path}", verbose)
    if not os.path.isfile(file_path):
        handle_non_critical_error(f"File {file_path} does not exist. Skipping.")
        return False

    preprocessed_file = preprocess_vcf(file_path)
    try:
        reader = vcfpy.Reader.from_path(preprocessed_file)
    except Exception as e:
        handle_non_critical_error(f"Could not open {file_path}: {str(e)}. Skipping.")
        return False
    if preprocessed_file != file_path:
        os.remove(preprocessed_file)


    header = reader.header

    if is_gvcf_header(header.lines) and not allow_gvcf:
        handle_non_critical_error(
            f"{file_path} appears to be a gVCF based on header metadata. Use --allow-gvcf to include gVCFs. Skipping."
        )
        return False

    fileformat = None
    for line in header.lines:
        if line.key == "fileformat":
            fileformat = line.value
            break

    expected_version = normalize_vcf_version(vcf_version)
    actual_version = normalize_vcf_version(fileformat)

    if fileformat is None or expected_version != actual_version:
        handle_non_critical_error(
            f"VCF version mismatch in {file_path}. Expected: {vcf_version}, Found: {fileformat}. Skipping.",
        )
        return False

    ref = None
    for line in header.lines:
        if line.key == "reference":
            ref = line.value
            break
    if ref is None or ref != ref_genome:
        handle_non_critical_error(f"Reference genome mismatch in {file_path}. Expected: {ref_genome}, Found: {ref}. Skipping.")
        return False


    try:
        from itertools import islice

        if not allow_gvcf:
            for record in islice(reader, 5):
                if any(
                    (
                        (alt_value := getattr(alt, "value", "")) in {"<NON_REF>", "NON_REF"}
                        or str(alt) in {"<NON_REF>", "SymbolicAllele('NON_REF')"}
                    )
                    for alt in getattr(record, "ALT", [])
                ):
                    handle_non_critical_error(
                        f"{file_path} appears to be a gVCF (found <NON_REF> ALT allele). Use --allow-gvcf to include gVCFs. Skipping."
                    )
                    return False

        _ = next(reader)
    except StopIteration:
        pass
    except Exception as e:
        handle_non_critical_error(f"Structural integrity check failed for {file_path}: {str(e)}. Skipping.")
        return False

    log_message(f"Validation passed for {file_path}", verbose)
    return True


def _extract_header_definitions(header, key):
    """Return an OrderedDict mapping header line IDs to their metadata."""

    try:
        header_lines = list(header.get_lines(key))
    except Exception:  # pragma: no cover - defensive; vcfpy may raise AttributeError
        header_lines = []

    definitions = OrderedDict()
    for line in header_lines:
        line_id = getattr(line, "id", None)
        if line_id is None and hasattr(line, "mapping"):
            line_id = line.mapping.get("ID")
        if line_id is None:
            continue
        metadata = copy.deepcopy(getattr(line, "mapping", {}))
        definitions[line_id] = metadata
    return definitions


def validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose=False, allow_gvcf=False):
    valid_vcfs = []
    log_message(f"Validating all VCF files in {input_dir}", verbose)
    reference_samples = []
    reference_contigs = None
    reference_info_defs = None
    reference_format_defs = None
    for file_path in glob.glob(os.path.join(input_dir, "*.vcf")):
        try:
            with open(file_path, "r", encoding="utf-8") as raw_vcf:
                header_lines = []
                for raw_line in raw_vcf:
                    if not raw_line.startswith("#"):
                        break
                    if raw_line.startswith("##"):
                        header_lines.append(raw_line.rstrip("\n"))
        except OSError as exc:
            handle_non_critical_error(f"Could not read header from {file_path}: {exc}. Skipping.")
            continue

        if is_gvcf_header(header_lines) and not allow_gvcf:
            handle_non_critical_error(
                f"{file_path} appears to be a gVCF based on header metadata. Use --allow-gvcf to include gVCFs. Skipping."
            )
            continue

        if validate_vcf(file_path, ref_genome, vcf_version, verbose=verbose, allow_gvcf=allow_gvcf):
            preprocessed_file = preprocess_vcf(file_path)
            reader = None
            try:
                reader = vcfpy.Reader.from_path(preprocessed_file)
            except Exception as exc:
                handle_critical_error(
                    f"Failed to reopen {file_path} after validation: {exc}."
                )
            try:
                header = reader.header

                current_samples = list(header.samples.names)
                current_contigs = _extract_header_definitions(header, "contig")
                current_info_defs = _extract_header_definitions(header, "INFO")
                current_format_defs = _extract_header_definitions(header, "FORMAT")

                if not reference_samples:
                    reference_samples = list(current_samples)
                    reference_contigs = current_contigs
                    reference_info_defs = current_info_defs
                    reference_format_defs = current_format_defs
                else:
                    new_samples = [
                        sample for sample in current_samples if sample not in reference_samples
                    ]
                    if new_samples:
                        reference_samples.extend(new_samples)
                        handle_non_critical_error(
                            "Sample columns differ between VCF shards. The merged VCF will "
                            "include the union of all sample names. Newly observed samples "
                            f"{new_samples} were appended to the merged header."
                        )

                    missing_samples = [
                        sample for sample in reference_samples if sample not in current_samples
                    ]
                    if missing_samples:
                        handle_non_critical_error(
                            f"File {file_path} is missing sample columns {missing_samples}. "
                            "These entries will be filled with missing ('.') values during merge."
                        )
                    if current_contigs != reference_contigs:
                        handle_critical_error(
                            "Contig definitions differ between VCF shards. MetaGap assumes vertical "
                            "concatenation of shards, so headers must match across shards. "
                            f"Expected contigs: {list(reference_contigs.keys())}; found in {file_path}: "
                            f"{list(current_contigs.keys())}."
                        )
                    if current_info_defs != reference_info_defs:
                        handle_critical_error(
                            "INFO field definitions differ between VCF shards. MetaGap assumes "
                            "vertical concatenation of shards, so headers must match across shards. "
                            f"Expected INFO IDs: {list(reference_info_defs.keys())}; found in {file_path}: "
                            f"{list(current_info_defs.keys())}."
                        )
                    if current_format_defs != reference_format_defs:
                        handle_critical_error(
                            "FORMAT field definitions differ between VCF shards. MetaGap assumes "
                            "vertical concatenation of shards, so headers must match across shards. "
                            f"Expected FORMAT IDs: {list(reference_format_defs.keys())}; found in {file_path}: "
                            f"{list(current_format_defs.keys())}."
                        )

                contig_order = {name: idx for idx, name in enumerate(current_contigs.keys())}
                last_contig_index = None
                last_position = None
                last_contig_name = None
                for record in reader:
                    chrom = record.CHROM
                    pos = record.POS
                    if contig_order:
                        if chrom not in contig_order:
                            handle_critical_error(
                                f"Encountered unknown contig '{chrom}' in {file_path}. MetaGap assumes vertical "
                                "concatenation of shards with consistent headers."
                            )
                        contig_index = contig_order[chrom]
                        if last_contig_index is None:
                            last_contig_index = contig_index
                            last_position = pos
                            last_contig_name = chrom
                            continue

                        if contig_index < last_contig_index or (
                            contig_index == last_contig_index and pos < last_position
                        ):
                            handle_critical_error(
                                f"VCF shard {file_path} is not coordinate-sorted: encountered {chrom}:{pos} after "
                                f"{last_contig_name}:{last_position}. MetaGap assumes vertical concatenation of shards, "
                                "so input shards must already be sorted."
                            )

                        last_contig_index = contig_index
                    else:
                        if last_contig_name is None:
                            last_contig_name = chrom
                            last_position = pos
                            continue

                        if chrom == last_contig_name:
                            if pos < last_position:
                                handle_critical_error(
                                    f"VCF shard {file_path} is not coordinate-sorted: encountered {chrom}:{pos} after "
                                    f"{last_contig_name}:{last_position}. MetaGap assumes vertical concatenation of shards, "
                                    "so input shards must already be sorted."
                                )
                        else:
                            last_contig_name = chrom
                            last_position = pos
                            continue

                    last_contig_name = chrom
                    last_position = pos
                valid_vcfs.append(file_path)
            finally:
                if reader is not None and hasattr(reader, "close"):
                    try:
                        reader.close()
                    except Exception:
                        pass
                if (
                    preprocessed_file != file_path
                    and preprocessed_file is not None
                    and os.path.exists(preprocessed_file)
                ):
                    os.remove(preprocessed_file)
        else:
            log_message(f"File {file_path} failed validation and is skipped.", verbose)
    if not valid_vcfs:
        handle_critical_error("No valid VCF files remain after validation. Aborting.")
    log_message("Validation completed. Valid VCF files: " + ", ".join(valid_vcfs), verbose)
    return valid_vcfs, reference_samples


import os, shutil, subprocess, datetime

def _parse_simple_metadata_line(line):
    stripped = line.strip()
    if not stripped.startswith("##") or "=" not in stripped:
        return None
    key, value = stripped[2:].split("=", 1)
    key = key.strip()
    value = value.strip()
    if not key:
        return None
    return key, value


def preprocess_vcf(file_path):
    """
    Check if the VCF file uses spaces instead of tabs for the column header and data lines.
    If so, create a temporary file where the column header (#CHROM) and all subsequent lines 
    are converted to be tab-delimited. Lines starting with "##" (metadata) are left unchanged.
    Returns the path to the file to be used (original or temporary).
    """
    import re
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    modified = False
    new_lines = []
    header_found = False  # indicates when the column header has been encountered
    for line in lines:
        if line.startswith("##"):
            # Do not change metadata lines (they may contain spaces that are part of the value)
            new_lines.append(line)
        elif line.startswith("#"):
            # This is the column header line (e.g. "#CHROM ...")
            new_line = re.sub(r'\s+', '\t', line.rstrip()) + "\n"
            new_lines.append(new_line)
            header_found = True
            if new_line != line:
                modified = True
        else:
            # Data lines: once header_found is True, convert spaces to tabs.
            if header_found:
                new_line = re.sub(r'\s+', '\t', line.rstrip()) + "\n"
                new_lines.append(new_line)
                if new_line != line:
                    modified = True
            else:
                new_lines.append(line)
    
    if modified:
        temp_file = file_path + ".tmp"
        with open(temp_file, 'w') as f:
            f.writelines(new_lines)
        return temp_file
    else:
        return file_path



def summarize_produced_vcfs(output_dir, fallback_vcf):
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
        # Prefer gzipped outputs when both compressed and uncompressed variants exist.
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
    (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries,
        sanitized_header_lines,
        serialized_sample_line,
    ) = parse_metadata_arguments(
        args,
        verbose,
        log_message=log_message,
        handle_critical_error=handle_critical_error,
    )

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
        detected_file, detected_fileformat, detected_reference = find_first_vcf_with_header(input_dir, verbose)
        if not ref_genome and detected_reference:
            ref_genome = detected_reference
        if not vcf_version and detected_fileformat:
            vcf_version = normalize_vcf_version(detected_fileformat)

    if not ref_genome:
        handle_critical_error("Reference genome build must be provided via --ref or auto-detectable from input files.")
    if not vcf_version:
        handle_critical_error("VCF version must be provided via --vcf-version or auto-detectable from input files.")

    if detection_needed and detected_file:
        log_message(
            f"Auto-detected metadata from {detected_file} -> reference={detected_reference or 'unknown'}, "
            f"version={detected_fileformat or 'unknown'}",
            verbose,
        )

    log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}", verbose)

    valid_files, sample_order = validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose)

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

# Ensure the repository root is on ``sys.path`` so the ``merge_vcf`` package
# can be imported when this file is executed directly.
_THIS_FILE = Path(__file__).resolve()
_REPO_ROOT = _THIS_FILE.parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from merge_vcf.main import *  # noqa: F401,F403 - re-export legacy symbols
from merge_vcf import main as _main


configure_merging_callbacks(
    MergeCallbacks(
        preprocess_vcf=preprocess_vcf,
        log_message=log_message,
        handle_critical_error=handle_critical_error,
        build_sample_metadata_line=build_sample_metadata_line,
        parse_sample_metadata_line=_parse_sample_metadata_line,
        apply_metadata_to_header=apply_metadata_to_header,
        vcfpy=vcfpy,
    )
)


if __name__ == "__main__":
    main()
