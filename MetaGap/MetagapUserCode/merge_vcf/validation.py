"""Validation routines for VCF shards and merged outputs."""

from __future__ import annotations

import copy
import glob
import gzip
import os
from collections import OrderedDict
from typing import Iterable, List, Optional

from . import VCFPY_AVAILABLE, vcfpy
from .logging_utils import (
    ValidationError,
    handle_critical_error,
    handle_non_critical_error,
    log_message,
    logger,
)
from .merging import preprocess_vcf


def is_gvcf_header(header_lines: Iterable) -> bool:
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


def normalize_vcf_version(version_value: Optional[str]):
    """Extract the numeric portion of the VCF version if present."""

    if not version_value:
        return version_value
    import re

    match = re.search(r"(\d+\.\d+)", version_value)
    if match:
        return match.group(1)
    return version_value.strip()


def find_first_vcf_with_header(input_dir: str, verbose: bool = False):
    """Locate the first readable VCF file and return its header metadata."""

    candidates = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    for file_path in candidates:
        try:
            reader = vcfpy.Reader.from_path(file_path)
        except Exception as exc:
            handle_non_critical_error(
                f"Unable to open {file_path} for auto-detection: {exc}. Trying next file.",
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


def validate_vcf(
    file_path: str,
    ref_genome: str,
    vcf_version: str,
    verbose: bool = False,
    allow_gvcf: bool = False,
) -> bool:
    log_message(f"Validating file: {file_path}", verbose)
    if not os.path.isfile(file_path):
        handle_non_critical_error(f"File {file_path} does not exist. Skipping.")
        return False

    preprocessed_file = preprocess_vcf(file_path)
    try:
        reader = vcfpy.Reader.from_path(preprocessed_file)
    except Exception as exc:
        handle_non_critical_error(f"Could not open {file_path}: {exc}. Skipping.")
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
        from itertools import islice

        if not allow_gvcf:
            for record in islice(reader, 5):
                if any(
                    (
                        (alt_value := getattr(alt, "value", ""))
                        in {"<NON_REF>", "NON_REF"}
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
            f"Structural integrity check failed for {file_path}: {exc}. Skipping."
        )
        return False

    log_message(f"Validation passed for {file_path}", verbose)
    return True


def _extract_header_definitions(header, key: str) -> "OrderedDict[str, dict]":
    """Return an OrderedDict mapping header line IDs to their metadata."""

    try:
        header_lines = list(header.get_lines(key))
    except Exception:
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


def validate_all_vcfs(
    input_dir: str,
    ref_genome: str,
    vcf_version: str,
    verbose: bool = False,
    allow_gvcf: bool = False,
):
    valid_vcfs: List[str] = []
    log_message(f"Validating all VCF files in {input_dir}", verbose)
    all_candidates = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    log_message(
        f"Discovered {len(all_candidates)} VCF shard(s) prior to validation.",
        verbose,
    )
    reference_samples: List[str] = []
    reference_contigs = None
    reference_info_defs = None
    reference_format_defs = None
    for file_path in all_candidates:
        try:
            with open(file_path, "r", encoding="utf-8") as raw_vcf:
                header_lines = []
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
                    f"Failed to reopen {file_path} after validation: {exc}.",
                    exc_cls=ValidationError,
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
                            f"{list(current_contigs.keys())}.",
                            exc_cls=ValidationError,
                        )
                    if current_info_defs != reference_info_defs:
                        handle_critical_error(
                            "INFO field definitions differ between VCF shards. MetaGap assumes "
                            "vertical concatenation of shards, so headers must match across shards. "
                            f"Expected INFO IDs: {list(reference_info_defs.keys())}; found in {file_path}: "
                            f"{list(current_info_defs.keys())}.",
                            exc_cls=ValidationError,
                        )
                    if current_format_defs != reference_format_defs:
                        handle_critical_error(
                            "FORMAT field definitions differ between VCF shards. MetaGap assumes "
                            "vertical concatenation of shards, so headers must match across shards. "
                            f"Expected FORMAT IDs: {list(reference_format_defs.keys())}; found in {file_path}: "
                            f"{list(current_format_defs.keys())}.",
                            exc_cls=ValidationError,
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
                                "concatenation of shards with consistent headers.",
                                exc_cls=ValidationError,
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
                                "so input shards must already be sorted.",
                                exc_cls=ValidationError,
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
                                    "so input shards must already be sorted.",
                                    exc_cls=ValidationError,
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
        handle_critical_error(
            "No valid VCF files remain after validation. Aborting.",
            exc_cls=ValidationError,
        )
    log_message("Validation completed. Valid VCF files: " + ", ".join(valid_vcfs), verbose)
    log_message(
        f"Validation summary: {len(valid_vcfs)} of {len(all_candidates)} shard(s) passed.",
        verbose,
    )
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
        if stripped in {".", "A", "G", "R"}:
            return True
        number = stripped
    try:
        return int(number) > 0
    except (TypeError, ValueError):
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


def _parse_header_info_definition(line: str):
    """Return (id, metadata) parsed from a ##INFO header line."""

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


def _read_vcf_without_vcfpy(file_path: str):
    """Yield header metadata and record lines for a VCF file without vcfpy."""

    opener = gzip.open if str(file_path).endswith(".gz") else open
    try:
        with opener(file_path, "rt", encoding="utf-8") as handle:
            header_lines = []
            column_header_seen = False
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

            records = []
            for raw_line in handle:
                line = raw_line.rstrip("\n")
                if line.strip():
                    records.append(line)
    except OSError as exc:
        handle_critical_error(
            f"Could not open {file_path}: {exc}.",
            exc_cls=ValidationError,
        )

    if not column_header_seen:
        handle_critical_error(
            f"Could not open {file_path}: Missing #CHROM header line.",
            exc_cls=ValidationError,
        )

    return header_lines, records


def _validate_merged_vcf_without_vcfpy(merged_vcf: str, verbose: bool = False):
    """Fallback validator that performs minimal checks without vcfpy."""

    header_lines, records = _read_vcf_without_vcfpy(merged_vcf)

    normalized_headers = [line.lower() for line in header_lines]
    if not any(line.startswith("##fileformat=") for line in normalized_headers):
        handle_critical_error(
            f"Missing required meta-information: ##fileformat in {merged_vcf} header.",
            exc_cls=ValidationError,
        )
    if not any(line.startswith("##reference=") for line in normalized_headers):
        handle_critical_error(
            f"Missing required meta-information: ##reference in {merged_vcf} header.",
            exc_cls=ValidationError,
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


def validate_merged_vcf(merged_vcf: str, verbose: bool = False):
    log_message(f"Starting validation of merged VCF: {merged_vcf}", verbose)
    if not os.path.isfile(merged_vcf):
        handle_critical_error(
            f"Merged VCF file {merged_vcf} does not exist.",
            exc_cls=ValidationError,
        )

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
                    f"Could not open {merged_vcf}: {gzip_exc}.",
                    exc_cls=ValidationError,
                )
        else:
            handle_critical_error(
                f"Could not open {merged_vcf}: {exc}.",
                exc_cls=ValidationError,
            )

    if reader is None:
        handle_critical_error(
            f"Could not open {merged_vcf}: Unknown error.",
            exc_cls=ValidationError,
        )

    header = reader.header
    required_meta = ["fileformat", "reference"]
    for meta in required_meta:
        found = False
        for line in header.lines:
            if line.key == meta:
                found = True
                break
        if not found:
            handle_critical_error(
                f"Missing required meta-information: ##{meta} in {merged_vcf} header.",
                exc_cls=ValidationError,
            )

    defined_info_ids = OrderedDict()
    try:
        info_ids_iterable = list(header.info_ids)
    except Exception:
        info_ids_iterable = []
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
    "is_gvcf_header",
    "normalize_vcf_version",
    "find_first_vcf_with_header",
    "validate_vcf",
    "validate_all_vcfs",
    "validate_merged_vcf",
]
