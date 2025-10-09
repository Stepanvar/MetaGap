"""Utility routines for validating VCF inputs and merged outputs.

This module centralises the validation helpers that were previously embedded in
the command-line merge script so that they can be re-used independently (e.g.,
from tests).
"""

from __future__ import annotations

import copy
import glob
import gzip
import logging
import os
import re
import sys
from collections import OrderedDict
from itertools import islice
from typing import Iterable, List, Optional, Sequence, Tuple

try:  # pragma: no cover - exercised when dependency is missing
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


LOG_FILE = "script_execution.log"
logger = logging.getLogger("vcf_merger")


def _configure_logger() -> None:
    if logger.handlers:
        return

    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(LOG_FILE)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    ch = logging.StreamHandler()
    ch.setFormatter(formatter)
    logger.addHandler(ch)


_configure_logger()


def log_message(message: str, verbose: bool = False) -> None:
    """Log *message* and optionally echo it to stdout when ``verbose`` is True."""

    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message: str) -> None:
    """Log a critical error and terminate the process."""

    log_message("CRITICAL ERROR: " + message)
    print("A critical error occurred. Check {} for details.".format(LOG_FILE))
    sys.exit(1)


def handle_non_critical_error(message: str) -> None:
    """Log a warning without terminating the process."""

    log_message("WARNING: " + message)
    print("Warning: " + message)


def is_gvcf_header(header_lines: Sequence[object]) -> bool:
    """Return True if the provided header metadata indicates a gVCF file."""

    for raw_line in header_lines:
        line_text: Optional[str] = None
        if hasattr(raw_line, "to_line"):
            try:
                line_text = raw_line.to_line()
            except Exception:  # pragma: no cover - defensive fallback
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


def normalize_vcf_version(version_value: Optional[str]) -> Optional[str]:
    """Extract the numeric portion of the VCF version if present."""

    if not version_value:
        return version_value
    match = re.search(r"(\d+\.\d+)", version_value)
    if match:
        return match.group(1)
    return version_value.strip()


def preprocess_vcf(file_path: str) -> str:
    """
    Check if the VCF file uses spaces instead of tabs for the column header and data lines.
    If so, create a temporary file where the column header (#CHROM) and all subsequent lines
    are converted to be tab-delimited. Lines starting with ``##`` (metadata) are left unchanged.
    Returns the path to the file to be used (original or temporary).
    """

    with open(file_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    modified = False
    new_lines: List[str] = []
    header_found = False  # indicates when the column header has been encountered
    for line in lines:
        if line.startswith("##"):
            # Do not change metadata lines (they may contain spaces that are part of the value)
            new_lines.append(line)
        elif line.startswith("#"):
            # This is the column header line (e.g. "#CHROM ...")
            new_line = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
            new_lines.append(new_line)
            header_found = True
            if new_line != line:
                modified = True
        else:
            # Data lines: once header_found is True, convert spaces to tabs.
            if header_found:
                new_line = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
                new_lines.append(new_line)
                if new_line != line:
                    modified = True
            else:
                new_lines.append(line)

    if modified:
        temp_file = file_path + ".tmp"
        with open(temp_file, "w", encoding="utf-8") as handle:
            handle.writelines(new_lines)
        return temp_file
    return file_path


def find_first_vcf_with_header(
    input_dir: str, verbose: bool = False
) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Locate the first readable VCF file and return its header metadata."""

    candidates = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    for file_path in candidates:
        try:
            reader = vcfpy.Reader.from_path(file_path)
        except Exception as exc:
            handle_non_critical_error(
                f"Unable to open {file_path} for auto-detection: {str(exc)}. Trying next file."
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
            f"Required metadata missing in {file_path} during auto-detection. Trying next file."
        )

    return None, None, None


def _extract_header_definitions(header, key: str) -> "OrderedDict[str, dict]":
    """Return an OrderedDict mapping header line IDs to their metadata."""

    try:
        header_lines = list(header.get_lines(key))
    except Exception:  # pragma: no cover - defensive; vcfpy may raise AttributeError
        header_lines = []

    definitions: "OrderedDict[str, dict]" = OrderedDict()
    for line in header_lines:
        line_id = getattr(line, "id", None)
        if line_id is None and hasattr(line, "mapping"):
            line_id = line.mapping.get("ID")
        if line_id is None:
            continue
        metadata = copy.deepcopy(getattr(line, "mapping", {}))
        definitions[line_id] = metadata
    return definitions


def validate_vcf(
    file_path: str,
    ref_genome: str,
    vcf_version: str,
    verbose: bool = False,
    allow_gvcf: bool = False,
) -> bool:
    """Validate that *file_path* matches the provided reference and VCF version."""

    log_message(f"Validating file: {file_path}", verbose)
    if not os.path.isfile(file_path):
        handle_non_critical_error(f"File {file_path} does not exist. Skipping.")
        return False

    preprocessed_file = preprocess_vcf(file_path)
    try:
        reader = vcfpy.Reader.from_path(preprocessed_file)
    except Exception as exc:
        handle_non_critical_error(f"Could not open {file_path}: {str(exc)}. Skipping.")
        return False
    if preprocessed_file != file_path:
        os.remove(preprocessed_file)

    header = reader.header
    fileformat = None
    for line in header.lines:
        if line.key == "fileformat":
            fileformat = line.value
            break

    expected_version = normalize_vcf_version(vcf_version)
    actual_version = normalize_vcf_version(fileformat)

    if fileformat is None or expected_version != actual_version:
        handle_non_critical_error(
            f"VCF version mismatch in {file_path}. Expected: {vcf_version}, Found: {fileformat}. Skipping."
        )
        return False

    ref = None
    for line in header.lines:
        if line.key == "reference":
            ref = line.value
            break
    if ref is None or ref != ref_genome:
        handle_non_critical_error(
            f"Reference genome mismatch in {file_path}. Expected: {ref_genome}, Found: {ref}. Skipping."
        )
        return False

    try:
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
    except Exception as exc:
        handle_non_critical_error(
            f"Structural integrity check failed for {file_path}: {str(exc)}. Skipping."
        )
        return False

    log_message(f"Validation passed for {file_path}", verbose)
    return True


def validate_all_vcfs(
    input_dir: str,
    ref_genome: str,
    vcf_version: str,
    verbose: bool = False,
    allow_gvcf: bool = False,
) -> Tuple[List[str], List[str]]:
    """Validate all VCF files inside *input_dir* and return the successful ones."""

    valid_vcfs: List[str] = []
    log_message(f"Validating all VCF files in {input_dir}", verbose)
    reference_samples: List[str] = []
    reference_contigs = None
    reference_info_defs = None
    reference_format_defs = None
    for file_path in glob.glob(os.path.join(input_dir, "*.vcf")):
        try:
            with open(file_path, "r", encoding="utf-8") as raw_vcf:
                header_lines: List[str] = []
                for raw_line in raw_vcf:
                    if not raw_line.startswith("#"):
                        break
                    if raw_line.startswith("##"):
                        header_lines.append(raw_line.rstrip("\n"))
        except OSError as exc:
            handle_non_critical_error(
                f"Could not read header from {file_path}: {exc}. Skipping."
            )
            continue

        if is_gvcf_header(header_lines) and not allow_gvcf:
            handle_non_critical_error(
                f"{file_path} appears to be a gVCF based on header metadata. Use --allow-gvcf to include gVCFs. Skipping."
            )
            continue

        if validate_vcf(
            file_path,
            ref_genome,
            vcf_version,
            verbose=verbose,
            allow_gvcf=allow_gvcf,
        ):
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


def _info_field_requires_value(number) -> bool:
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


def _has_null_value(value) -> bool:
    """Return True if *value* or any nested value is considered missing."""

    if value is None:
        return True
    if isinstance(value, str):
        stripped = value.strip()
        return stripped == "" or stripped == "."
    if isinstance(value, (list, tuple)):
        return any(_has_null_value(v) for v in value)
    return False


def _parse_header_info_definition(line: str) -> Tuple[Optional[str], dict]:
    """Return ``(id, metadata)`` parsed from a ``##INFO`` header line."""

    lower = line.lower()
    if not lower.startswith("##info"):
        return None, {}
    start = line.find("<")
    end = line.rfind(">")
    if start == -1 or end == -1 or end <= start + 1:
        return None, {}

    payload = line[start + 1 : end]
    entries = {}
    for part in payload.split(","):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        entries[key.strip().upper()] = value.strip()

    info_id = entries.get("ID")
    if not info_id:
        return None, {}
    return info_id, entries


def _read_vcf_without_vcfpy(file_path: str) -> Tuple[List[str], List[str]]:
    """Yield header metadata and record lines for a VCF file without vcfpy."""

    opener = gzip.open if str(file_path).endswith(".gz") else open
    column_header_seen = False
    try:
        with opener(file_path, "rt", encoding="utf-8") as handle:
            header_lines: List[str] = []
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if not line:
                    continue
                if line.startswith("##"):
                    header_lines.append(line)
                    continue
                if line.startswith("#CHROM"):
                    column_header_seen = True
                    break

            records: List[str] = []
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if line.strip():
                    records.append(line)
    except OSError as exc:
        handle_critical_error(f"Could not open {file_path}: {exc}.")

    if not column_header_seen:
        handle_critical_error(
            f"Could not open {file_path}: Missing #CHROM header line."
        )

    return header_lines, records


def _validate_merged_vcf_without_vcfpy(merged_vcf: str, verbose: bool = False) -> None:
    """Fallback validator that performs minimal checks without vcfpy."""

    header_lines, records = _read_vcf_without_vcfpy(merged_vcf)

    normalized_headers = [line.lower() for line in header_lines]
    if not any(line.startswith("##fileformat=") for line in normalized_headers):
        handle_critical_error(
            f"Missing required meta-information: ##fileformat in {merged_vcf} header."
        )
    if not any(line.startswith("##reference=") for line in normalized_headers):
        handle_critical_error(
            f"Missing required meta-information: ##reference in {merged_vcf} header."
        )

    info_definitions = {}
    for line in header_lines:
        info_id, entries = _parse_header_info_definition(line)
        if not info_id:
            continue
        info_definitions[info_id] = entries

    required_info_ids = {
        info_id
        for info_id, entries in info_definitions.items()
        if _info_field_requires_value(entries.get("NUMBER"))
    }

    for record_line in records:
        if record_line.startswith("#"):
            continue
        fields = record_line.split("\t")
        if len(fields) < 8:
            continue
        chrom, pos = fields[0], fields[1]
        record_label = f"{chrom}:{pos}"
        info_field = fields[7].strip()
        if info_field in {"", "."}:
            handle_non_critical_error(
                f"Record {record_label} in {merged_vcf} is missing INFO data. Continuing without INFO validation for this record."
            )
            continue

        info_map = {}
        for entry in info_field.split(";"):
            if not entry:
                continue
            if "=" in entry:
                key, value = entry.split("=", 1)
                info_map[key] = value
            else:
                info_map[entry] = ""

        undefined_keys = sorted(
            key for key in info_map.keys() if key not in info_definitions
        )
        if undefined_keys:
            logger.warning(
                "Record %s in %s has INFO fields not present in header definitions: %s.",
                record_label,
                merged_vcf,
                ", ".join(undefined_keys),
            )

        missing_keys = sorted(required_info_ids.difference(info_map.keys()))
        if missing_keys:
            handle_non_critical_error(
                f"Record {record_label} in {merged_vcf} is missing required INFO fields: {', '.join(missing_keys)}."
            )

        null_keys = [key for key, value in info_map.items() if _has_null_value(value)]
        if null_keys:
            handle_non_critical_error(
                f"Record {record_label} in {merged_vcf} has INFO fields with null values: {', '.join(sorted(null_keys))}."
            )

    log_message(f"Validation completed successfully for merged VCF: {merged_vcf}", verbose)
    print(f"Validation completed successfully for merged VCF: {merged_vcf}")


def validate_merged_vcf(merged_vcf: str, verbose: bool = False) -> None:
    """Validate the merged VCF using vcfpy when available."""

    log_message(f"Starting validation of merged VCF: {merged_vcf}", verbose)
    if not os.path.isfile(merged_vcf):
        handle_critical_error(f"Merged VCF file {merged_vcf} does not exist.")

    if not VCFPY_AVAILABLE:
        _validate_merged_vcf_without_vcfpy(merged_vcf, verbose=verbose)
        return

    reader = None
    gz_stream = None
    try:
        reader = vcfpy.Reader.from_path(merged_vcf)
    except Exception as exc:
        if merged_vcf.endswith(".gz"):
            try:
                gz_stream = gzip.open(merged_vcf, "rt")
                reader = vcfpy.Reader.from_stream(gz_stream, merged_vcf)
            except Exception as gzip_exc:
                if gz_stream is not None:
                    gz_stream.close()
                handle_critical_error(
                    f"Could not open {merged_vcf}: {str(gzip_exc)}."
                )
        else:
            handle_critical_error(f"Could not open {merged_vcf}: {str(exc)}.")

    if reader is None:
        handle_critical_error(f"Could not open {merged_vcf}: Unknown error.")

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
    info_ids_iterable: Iterable = []
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
    except Exception as exc:  # pragma: no cover - defensive logging
        encountered_exception = True
        logger.warning("Error while parsing records in %s: %s", merged_vcf, exc)

    if not encountered_exception:
        log_message(f"Validation completed successfully for merged VCF: {merged_vcf}", verbose)
        print(f"Validation completed successfully for merged VCF: {merged_vcf}")

    try:
        reader.close()
    except Exception:
        pass
    if gz_stream is not None:
        try:
            gz_stream.close()
        except Exception:
            pass


__all__ = [
    "LOG_FILE",
    "VCFPY_AVAILABLE",
    "handle_critical_error",
    "handle_non_critical_error",
    "is_gvcf_header",
    "log_message",
    "logger",
    "normalize_vcf_version",
    "preprocess_vcf",
    "find_first_vcf_with_header",
    "validate_vcf",
    "validate_all_vcfs",
    "validate_merged_vcf",
    "vcfpy",
]

