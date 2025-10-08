#!/usr/bin/env python3
"""
This script consolidates multiple VCF files into one merged VCF file.
It replicates the bash script functionality:
  - Parses command-line options for the input directory, output directory, reference genome, VCF version, and optional metadata.
  - Validates individual VCF files for header fileformat and reference genome.
  - Merges valid VCFs using vcfpy.
  - Appends metadata (if provided via CLI) to the merged VCF header.
  - Performs a final validation of the merged VCF.
  - Logs execution details to a log file.

Requirements:
pip install vcfpy
"""

import os
import sys
import glob
import argparse
import csv
import logging
import datetime
import re
import copy
from collections import OrderedDict

try:
    import vcfpy
except ImportError:
    sys.exit("Error: vcfpy package is required. Please install it with 'pip install vcfpy'.")

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
        "--output-dir",
        required=False,
        help="Directory to write the merged VCF file (defaults to the input directory).",
    )
    parser.add_argument(
        "--ref",
        required=True,
        help="Reference genome build expected in the VCF headers.",
    )
    parser.add_argument(
        "--vcf-version",
        required=True,
        help="VCF version expected in the input files.",
    )
    parser.add_argument(
        "--meta",
        required=False,
        help="Optional CSV list of key=value pairs describing SAMPLE metadata.",
    )

    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = args.input_dir
    return args

def _format_sample_metadata_value(value: str) -> str:
    """Return a VCF-safe representation of the provided metadata value."""

    value = str(value)
    needs_quotes = any(char in value for char in [" ", ",", "\t", "\"", "<", ">", "="])
    if not needs_quotes:
        return value

    escaped = value.replace("\\", "\\\\").replace("\"", "\\\"")
    return f'"{escaped}"'


def build_sample_metadata_line(entries: "OrderedDict[str, str]") -> str:
    """Serialize an ordered mapping into a single ##SAMPLE metadata line."""

    id_value = entries.get("ID", "").strip()
    if not id_value:
        raise ValueError("Sample metadata must include a non-empty ID value.")

    parts = []
    for key, raw_value in entries.items():
        if raw_value is None:
            continue
        value = str(raw_value).strip()
        if not value and key != "ID":
            continue
        parts.append(f"{key}={_format_sample_metadata_value(value)}")

    serialized = ",".join(parts)
    return f"##SAMPLE=<{serialized}>"


def parse_metadata_string(raw_metadata):
    """Parse a comma-separated list of key=value pairs into an ordered mapping."""

    if not raw_metadata:
        return None

    try:
        entries = next(csv.reader([raw_metadata]))
    except StopIteration:
        return None

    mapping = OrderedDict()
    for entry in entries:
        if not entry:
            continue
        if "=" not in entry:
            raise ValueError(f"Metadata entry '{entry}' is missing an '=' separator.")
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError("Metadata keys cannot be empty.")
        mapping[key] = value

    if not mapping:
        return None

    if "ID" not in mapping or not mapping["ID"].strip():
        raise ValueError("Metadata must include a non-empty ID entry.")

    return mapping


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


def validate_vcf(file_path, ref_genome, vcf_version, verbose=False):
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

    for line in header.lines:
        if isinstance(getattr(line, "key", None), str) and line.key.upper().startswith("GVCF"):
            handle_non_critical_error(
                f"{file_path} appears to be a gVCF (found header line {line.key}). Skipping."
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

        for record in islice(reader, 5):
            if any(
                (
                    (alt_value := getattr(alt, "value", "")) in {"<NON_REF>", "NON_REF"}
                    or str(alt) in {"<NON_REF>", "SymbolicAllele('NON_REF')"}
                )
                for alt in getattr(record, "ALT", [])
            ):
                handle_non_critical_error(
                    f"{file_path} appears to be a gVCF (found <NON_REF> ALT allele). Skipping."
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


def validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose=False):
    valid_vcfs = []
    log_message(f"Validating all VCF files in {input_dir}", verbose)
    reference_samples = None
    reference_contigs = None
    reference_info_defs = None
    reference_format_defs = None
    for file_path in glob.glob(os.path.join(input_dir, "*.vcf")):
        if validate_vcf(file_path, ref_genome, vcf_version, verbose):
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
                current_contigs = copy.deepcopy(header.contigs)
                current_info_defs = copy.deepcopy(header.info_defs)
                current_format_defs = copy.deepcopy(header.format_defs)

                if reference_samples is None:
                    reference_samples = current_samples
                    reference_contigs = current_contigs
                    reference_info_defs = current_info_defs
                    reference_format_defs = current_format_defs
                else:
                    if current_samples != reference_samples:
                        handle_critical_error(
                            "Sample names differ between VCF shards. MetaGap assumes vertical "
                            "concatenation of shards, so headers and sample ordering must match. "
                            f"Expected samples: {reference_samples}; found in {file_path}: {current_samples}."
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

                contig_order = {name: idx for idx, name in enumerate(header.contigs.keys())}
                last_contig_index = None
                last_position = None
                last_contig_name = None
                for record in reader:
                    chrom = record.CHROM
                    pos = record.POS
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

                    if contig_index == last_contig_index:
                        last_position = pos
                        last_contig_name = chrom
                    else:
                        last_contig_index = contig_index
                        last_position = pos
                        last_contig_name = chrom
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
    return valid_vcfs


def union_headers(valid_files):
    """Return a unified header containing all INFO/FORMAT/FILTER/contig lines."""

    if not valid_files:
        handle_critical_error("No valid VCF files remain to merge. Aborting.")

    temporary_paths = []

    def _cleanup_temp_paths():
        for temp_path in temporary_paths:
            try:
                if temp_path and os.path.exists(temp_path):
                    os.remove(temp_path)
            except OSError:
                pass

    def _get_index_mapping(header_obj, key):
        indices = getattr(header_obj, "_indices", {})
        mapping = indices.get(key, {})
        if isinstance(mapping, dict):
            return mapping
        return {}

    first_path = preprocess_vcf(valid_files[0])
    if first_path != valid_files[0]:
        temporary_paths.append(first_path)

    try:
        reader = vcfpy.Reader.from_path(first_path)
    except Exception as exc:
        _cleanup_temp_paths()
        handle_critical_error(f"Failed to open the first valid VCF for header union: {exc}")

    try:
        unified_header = reader.header.copy()
    finally:
        reader.close()

    contig_mapping = _get_index_mapping(unified_header, "contig")
    unified_header.contigs = dict(contig_mapping)

    existing_info_ids = set(unified_header.info_ids())
    existing_format_ids = set(unified_header.format_ids())
    existing_filter_ids = set(unified_header.filter_ids())
    existing_contig_ids = set(contig_mapping.keys())

    try:
        for file_path in valid_files[1:]:
            preprocessed_path = preprocess_vcf(file_path)
            if preprocessed_path != file_path:
                temporary_paths.append(preprocessed_path)
            try:
                reader = vcfpy.Reader.from_path(preprocessed_path)
            except Exception as exc:
                _cleanup_temp_paths()
                handle_critical_error(
                    f"Failed to open {file_path} while unifying headers: {exc}"
                )

            try:
                header = reader.header
                info_lines = _get_index_mapping(header, "INFO")
                for info_id, info_line in info_lines.items():
                    if info_id not in existing_info_ids:
                        unified_header.add_line(info_line.copy())
                        existing_info_ids.add(info_id)

                format_lines = _get_index_mapping(header, "FORMAT")
                for format_id, format_line in format_lines.items():
                    if format_id not in existing_format_ids:
                        unified_header.add_line(format_line.copy())
                        existing_format_ids.add(format_id)

                filter_lines = _get_index_mapping(header, "FILTER")
                for filter_id, filter_line in filter_lines.items():
                    if filter_id not in existing_filter_ids:
                        unified_header.add_line(filter_line.copy())
                        existing_filter_ids.add(filter_id)

                contig_lines = _get_index_mapping(header, "contig")
                for contig_id, contig_line in contig_lines.items():
                    if contig_id not in existing_contig_ids:
                        unified_header.add_line(contig_line.copy())
                        existing_contig_ids.add(contig_id)
            finally:
                reader.close()

        unified_header.contigs = _get_index_mapping(unified_header, "contig")
        return unified_header
    finally:
        _cleanup_temp_paths()


def merge_vcfs(valid_files, output_dir, verbose=False):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    merged_filename = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")
    log_message("Merging VCF files...", verbose)
    header = union_headers(valid_files)

    writer = vcfpy.Writer.from_path(merged_filename, header)
    for file_path in valid_files:
        log_message(f"Processing file: {file_path}", verbose)
        preprocessed_file = preprocess_vcf(file_path)
        try:
            reader = vcfpy.Reader.from_path(preprocessed_file)
            for record in reader:
                writer.write_record(record)
        except Exception as e:
            handle_non_critical_error(f"Error processing {file_path}: {str(e)}. Skipping remaining records.")
        if preprocessed_file != file_path:
            os.remove(preprocessed_file)

    writer.close()
    log_message(f"Merged VCF file created successfully: {merged_filename}", verbose)
    return merged_filename

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


def _parse_sample_metadata_line(line):
    stripped = line.strip()
    if not stripped.startswith("##SAMPLE="):
        return None

    payload = stripped[len("##SAMPLE=") :].strip()
    if not payload.startswith("<") or not payload.endswith(">"):
        raise ValueError("Malformed ##SAMPLE metadata detected; expected angle bracket encapsulation.")

    inner = payload[1:-1]
    entries = []
    current = []
    in_quotes = False
    escape = False

    for char in inner:
        if escape:
            current.append(char)
            escape = False
            continue

        if char == "\\":
            if in_quotes:
                escape = True
                continue
            current.append(char)
            continue

        if char == '"':
            in_quotes = not in_quotes
            continue

        if char == "," and not in_quotes:
            entry = "".join(current).strip()
            if entry:
                entries.append(entry)
            current = []
            continue

        current.append(char)

    if current:
        entry = "".join(current).strip()
        if entry:
            entries.append(entry)

    mapping = OrderedDict()
    for entry in entries:
        if "=" not in entry:
            continue
        key, raw_value = entry.split("=", 1)
        key = key.strip()
        value = raw_value.strip()
        mapping[key] = value

    if "ID" not in mapping or not mapping["ID"]:
        raise ValueError("Sample metadata must include an ID entry.")

    return mapping


def append_metadata_to_merged_vcf(merged_vcf, metadata, verbose=False):
    if not metadata:
        log_message("No metadata provided. Skipping metadata append step.", verbose)
        return

    log_message("Appending metadata to merged VCF header.", verbose)

    try:
        reader = vcfpy.Reader.from_path(merged_vcf)
    except Exception as exc:
        handle_non_critical_error(f"Could not open {merged_vcf}: {exc}. Skipping metadata append.")
        return

    header = reader.header

    try:
        metadata_line = build_sample_metadata_line(metadata)
        metadata_mapping = _parse_sample_metadata_line(metadata_line)
    except ValueError as exc:
        reader.close()
        handle_non_critical_error(str(exc))
        return

    header.lines = [line for line in header.lines if getattr(line, "key", None) != "SAMPLE"]

    sample_line = vcfpy.SampleHeaderLine.from_mapping(metadata_mapping)
    header.add_line(sample_line)

    temp_path = merged_vcf + ".tmp"
    writer = vcfpy.Writer.from_path(temp_path, header)

    for record in reader:
        writer.write_record(record)

    writer.close()
    reader.close()

    os.replace(temp_path, merged_vcf)

    log_message("Metadata appended successfully.", verbose)

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



def _info_field_requires_value(number):
    """Return True when the INFO definition expects an accompanying value."""

    if number is None:
        return False
    if isinstance(number, str):
        stripped = number.strip()
        if stripped == "":
            return False
        if stripped == "0":
            return False
        # symbolic numbers (., A, G, R) all expect a value when present
        if stripped in {".", "A", "G", "R"}:
            return True
        number = stripped
    try:
        return int(number) > 0
    except (TypeError, ValueError):
        # Unknown token, err on the side of requiring a value
        return True


def _has_null_value(value):
    """Return True if *value* or any nested value is ``None``."""

    if value is None:
        return True
    if isinstance(value, (list, tuple)):
        return any(_has_null_value(v) for v in value)
    return False


def validate_merged_vcf(merged_vcf, verbose=False):
    log_message(f"Starting validation of merged VCF: {merged_vcf}", verbose)
    if not os.path.isfile(merged_vcf):
        handle_critical_error(f"Merged VCF file {merged_vcf} does not exist.")
        
    try:
        reader = vcfpy.Reader.from_path(merged_vcf)
    except Exception as e:
        handle_critical_error(f"Could not open {merged_vcf}: {str(e)}.")


    header = reader.header
    required_meta = ["fileformat", "reference"]
    for meta in required_meta:
        found = False
        for line in header.lines:
            if line.key == meta:
                found = True
                break
        if not found:
            handle_critical_error(f"Missing required meta-information: ##{meta} in {merged_vcf} header.")


    defined_info_ids = OrderedDict()
    info_ids_iterable = []
    if hasattr(header, "info_ids"):
        candidate = header.info_ids
        if callable(candidate):
            try:
                info_ids_iterable = list(candidate())
            except Exception:
                info_ids_iterable = []
        else:
            info_ids_iterable = list(candidate)

    for info_id in info_ids_iterable:
        try:
            info_def = header.get_info_field_info(info_id)
        except Exception:
            info_def = None
        if isinstance(info_def, vcfpy.header.InfoHeaderLine):
            defined_info_ids[info_id] = info_def

    if not defined_info_ids:
        for line in header.lines:
            if isinstance(line, vcfpy.header.InfoHeaderLine):
                defined_info_ids[line.id] = line

    required_info_ids = {
        info_id
        for info_id, info_def in defined_info_ids.items()
        if _info_field_requires_value(getattr(info_def, "number", None))
    }
    required_info_ids = {"AC", "AN", "AF"}.intersection(defined_info_ids)

    encountered_exception = False
    try:
        for record in reader:
            info_map = getattr(record, "INFO", None)
            record_label = f"{getattr(record, 'CHROM', '?')}:{getattr(record, 'POS', '?')}"

            if info_map is None:
                handle_non_critical_error(
                    f"Record {record_label} in {merged_vcf} is missing INFO data. Continuing without INFO validation for this record."
                )
                continue

            if not isinstance(info_map, dict):
                handle_non_critical_error(
                    f"Record {record_label} in {merged_vcf} has an unexpected INFO type: {type(info_map).__name__}. Continuing without INFO validation for this record."
                )
                continue

            record_info_keys = set(info_map.keys())

            undefined_keys = sorted(record_info_keys.difference(defined_info_ids))
            if undefined_keys:
                logger.warning(
                    "Record %s in %s has INFO fields not present in header definitions: %s.",
                    record_label,
                    merged_vcf,
                    ", ".join(undefined_keys),
                )

            missing_keys = sorted(required_info_ids.difference(record_info_keys))
            if missing_keys:
                handle_non_critical_error(
                    f"Record {record_label} in {merged_vcf} is missing required INFO fields: {', '.join(missing_keys)}."
                )

            null_keys = [key for key, value in info_map.items() if _has_null_value(value)]
            if null_keys:
                handle_non_critical_error(
                    f"Record {record_label} in {merged_vcf} has INFO fields with null values: {', '.join(sorted(null_keys))}."
                )

    except SystemExit:
        raise
    except Exception as exc:
        encountered_exception = True
        logger.warning("Error while parsing records in %s: %s", merged_vcf, exc)

    if not encountered_exception:
        log_message(f"Validation completed successfully for merged VCF: {merged_vcf}", verbose)
        print(f"Validation completed successfully for merged VCF: {merged_vcf}")


def main():
    args = parse_arguments()
    verbose = False

    log_message(
        "Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        verbose,
    )

    input_dir = args.input_dir
    if not os.path.isdir(input_dir):
        handle_critical_error(f"Input directory {input_dir} does not exist or is not a directory.")
    log_message("Input directory: " + input_dir, verbose)

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)
    log_message("Output directory: " + output_dir, verbose)

    ref_genome = args.ref
    vcf_version = normalize_vcf_version(args.vcf_version)
    log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}", verbose)

    try:
        metadata_entries = parse_metadata_string(args.meta)
    except ValueError as exc:
        handle_critical_error(str(exc))

    valid_files = validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose)
    merged_vcf = merge_vcfs(valid_files, output_dir, verbose)
    append_metadata_to_merged_vcf(merged_vcf, metadata_entries, verbose)
    validate_merged_vcf(merged_vcf, verbose)

    print("----------------------------------------")
    print("Script Execution Summary:")
    print("Merged VCF File: " + merged_vcf)
    print("Log File: " + LOG_FILE)
    print("For detailed logs, refer to " + LOG_FILE)
    print("----------------------------------------")
    log_message("Script execution completed successfully.", verbose)


if __name__ == "__main__":
    main()
