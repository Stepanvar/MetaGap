"""Validation routines for VCF shards and merged outputs."""

from __future__ import annotations

import copy
import gzip
import os
from collections import OrderedDict
from dataclasses import dataclass
from typing import Iterable, List, Optional

try:  # pragma: no cover - dependency availability checked dynamically
    import pysam  # type: ignore

    PYSAM_AVAILABLE = True
except ImportError:  # pragma: no cover - exercised when dependency is missing
    pysam = None  # type: ignore
    PYSAM_AVAILABLE = False

from . import VCFPY_AVAILABLE, vcfpy
from .logging_utils import (
    handle_critical_error,
    handle_non_critical_error,
    log_message,
    logger,
)
from .merging import preprocess_vcf


@dataclass(frozen=True)
class PreparedVCFInput:
    """Metadata describing a discovered VCF input file."""

    source_path: str
    """The path that was originally discovered (may be ``.vcf`` or ``.vcf.gz``)."""

    compressed_path: str
    """Path to the BGZF-compressed ``.vcf.gz`` representation."""

    index_path: str
    """Path to the corresponding ``.tbi`` index file."""


def _iter_candidate_vcfs(input_dir: str) -> List[str]:
    """Return sorted candidate VCF paths discovered under *input_dir*."""

    candidates: List[str] = []
    for root, _, files in os.walk(input_dir):
        for name in files:
            lower = name.lower()
            if lower.endswith(".vcf") or lower.endswith(".vcf.gz"):
                candidates.append(os.path.join(root, name))
    candidates.sort()
    return candidates


def _contains_non_ref_alt(record) -> bool:
    """Return True if *record* includes a symbolic NON_REF alternate allele."""

    for alt in getattr(record, "ALT", []) or []:
        alt_value = getattr(alt, "value", None)
        if alt_value in {"<NON_REF>", "NON_REF"}:
            return True
        alt_str = str(alt)
        if alt_str in {"<NON_REF>", "SymbolicAllele('NON_REF')"}:
            return True
    return False


def discover_and_prepare_inputs(
    input_dir: str, verbose: bool = False, allow_gvcf: bool = False
) -> List[PreparedVCFInput]:
    """Discover VCF inputs, compress/index them, and screen for gVCF markers."""

    if not os.path.isdir(input_dir):
        handle_critical_error(f"Input directory does not exist: {input_dir}")

    if not PYSAM_AVAILABLE:
        handle_critical_error(
            "pysam dependency is required to prepare VCF inputs. Please install pysam."
        )

    if not VCFPY_AVAILABLE:
        handle_critical_error(
            "vcfpy dependency is required to inspect VCF headers. Please install vcfpy."
        )

    prepared: List[PreparedVCFInput] = []
    discovered_paths = _iter_candidate_vcfs(input_dir)
    log_message(
        f"Discovered {len(discovered_paths)} potential VCF input file(s) in {input_dir}",
        verbose,
    )

    for path in discovered_paths:
        source_path = path
        compressed_path = path
        created_compressed = False

        if not path.lower().endswith(".vcf.gz"):
            compressed_path = path + ".gz"
            try:
                pysam.tabix_compress(path, compressed_path, force=True)
                created_compressed = True
                log_message(f"Compressed {path} -> {compressed_path}", verbose)
            except Exception as exc:  # pragma: no cover - defensive logging path
                handle_non_critical_error(
                    f"Failed to BGZF-compress {path}: {exc}. Skipping."
                )
                continue

        index_path = compressed_path + ".tbi"
        try:
            pysam.tabix_index(compressed_path, preset="vcf", force=True)
            log_message(f"Indexed {compressed_path} -> {index_path}", verbose)
        except Exception as exc:  # pragma: no cover - defensive logging path
            handle_non_critical_error(
                f"Failed to create tabix index for {compressed_path}: {exc}. Skipping."
            )
            if created_compressed:
                try:
                    if os.path.exists(compressed_path):
                        os.remove(compressed_path)
                except OSError:
                    pass
            continue

        reader = None
        skip_file = False
        try:
            reader = vcfpy.Reader.from_path(compressed_path)
            if not allow_gvcf:
                if is_gvcf_header(reader.header.lines):
                    handle_non_critical_error(
                        f"{source_path} appears to be a gVCF based on header metadata. "
                        "Use --allow-gvcf to include gVCFs. Skipping."
                    )
                    skip_file = True
                else:
                    try:
                        from itertools import islice

                        for record in islice(reader, 5):
                            if _contains_non_ref_alt(record):
                                handle_non_critical_error(
                                    f"{source_path} appears to be a gVCF (found <NON_REF> ALT allele). "
                                    "Use --allow-gvcf to include gVCFs. Skipping."
                                )
                                skip_file = True
                                break
                    except Exception as exc:
                        handle_non_critical_error(
                            f"Failed to inspect records in {source_path}: {exc}. Skipping."
                        )
                        skip_file = True
        except Exception as exc:
            handle_non_critical_error(
                f"Could not open {compressed_path} for validation: {exc}. Skipping."
            )
            skip_file = True
        finally:
            if reader is not None:
                try:
                    reader.close()
                except Exception:
                    pass

        if skip_file:
            if created_compressed:
                for cleanup_path in (compressed_path, index_path):
                    try:
                        if cleanup_path and os.path.exists(cleanup_path):
                            os.remove(cleanup_path)
                    except OSError:
                        pass
            continue

        prepared.append(
            PreparedVCFInput(
                source_path=source_path,
                compressed_path=compressed_path,
                index_path=index_path,
            )
        )

    log_message(
        f"Prepared {len(prepared)} VCF input file(s) for downstream validation.",
        verbose,
    )
    return prepared


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

    candidates = _iter_candidate_vcfs(input_dir)
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
    reference_samples: List[str] = []
    reference_contigs: "OrderedDict[str, dict]" | None = None
    reference_info_defs: "OrderedDict[str, dict]" | None = None
    reference_filter_defs: "OrderedDict[str, dict]" | None = None
    reference_format_defs: "OrderedDict[str, dict]" | None = None
    prepared_inputs = discover_and_prepare_inputs(
        input_dir, verbose=verbose, allow_gvcf=allow_gvcf
    )

    for prepared in prepared_inputs:
        file_path = prepared.source_path
        opener = gzip.open if file_path.endswith(".gz") else open
        mode = "rt" if opener is gzip.open else "r"
        try:
            with opener(file_path, mode, encoding="utf-8") as raw_vcf:
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
                    f"Failed to reopen {file_path} after validation: {exc}."
                )
            try:
                header = reader.header

                current_samples = list(header.samples.names)
                current_contigs = _extract_header_definitions(header, "contig")
                current_filter_defs = _extract_header_definitions(header, "FILTER")
                current_info_defs = _extract_header_definitions(header, "INFO")
                current_format_defs = _extract_header_definitions(header, "FORMAT")

                if not reference_samples:
                    reference_samples = list(current_samples)
                    reference_contigs = current_contigs
                    reference_filter_defs = current_filter_defs
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
                    if reference_contigs is None:
                        reference_contigs = OrderedDict()
                    for contig_id, mapping in current_contigs.items():
                        existing = reference_contigs.get(contig_id)
                        if existing is None:
                            reference_contigs[contig_id] = mapping
                            continue
                        if existing != mapping:
                            handle_critical_error(
                                "Contig header definitions conflict across shards. "
                                f"For contig '{contig_id}', previous definition was {existing!r} "
                                f"but {file_path} defines {mapping!r}."
                            )

                    if reference_filter_defs is None:
                        reference_filter_defs = OrderedDict()
                    for filter_id, mapping in current_filter_defs.items():
                        existing = reference_filter_defs.get(filter_id)
                        if existing is None:
                            reference_filter_defs[filter_id] = mapping
                            continue
                        if existing.get("Description") != mapping.get("Description"):
                            handle_critical_error(
                                "FILTER header definitions conflict across shards. "
                                f"For filter '{filter_id}', previous Description was "
                                f"{existing.get('Description')!r} but {file_path} defines "
                                f"{mapping.get('Description')!r}."
                            )

                    if reference_info_defs is None:
                        reference_info_defs = OrderedDict()
                    for info_id, mapping in current_info_defs.items():
                        existing = reference_info_defs.get(info_id)
                        if existing is None:
                            reference_info_defs[info_id] = mapping
                            continue
                        for key in ("Number", "Type", "Description"):
                            if existing.get(key) != mapping.get(key):
                                handle_critical_error(
                                    "INFO header definitions conflict across shards. "
                                    f"For INFO '{info_id}', field '{key}' differed: "
                                    f"{existing.get(key)!r} vs {mapping.get(key)!r} in {file_path}."
                                )

                    if reference_format_defs is None:
                        reference_format_defs = OrderedDict()
                    for format_id, mapping in current_format_defs.items():
                        existing = reference_format_defs.get(format_id)
                        if existing is None:
                            reference_format_defs[format_id] = mapping
                            continue
                        for key in ("Number", "Type", "Description"):
                            if existing.get(key) != mapping.get(key):
                                handle_critical_error(
                                    "FORMAT header definitions conflict across shards. "
                                    f"For FORMAT '{format_id}', field '{key}' differed: "
                                    f"{existing.get(key)!r} vs {mapping.get(key)!r} in {file_path}."
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
        handle_critical_error(f"Could not open {file_path}: {exc}.")

    if not column_header_seen:
        handle_critical_error(
            f"Could not open {file_path}: Missing #CHROM header line."
        )

    return header_lines, records


def _validate_merged_vcf_without_vcfpy(merged_vcf: str, verbose: bool = False):
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


def validate_merged_vcf(merged_vcf: str, verbose: bool = False):
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
                    f"Could not open {merged_vcf}: {gzip_exc}."
                )
        else:
            handle_critical_error(f"Could not open {merged_vcf}: {exc}.")

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
            handle_critical_error(
                f"Missing required meta-information: ##{meta} in {merged_vcf} header."
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
    "PreparedVCFInput",
    "discover_and_prepare_inputs",
    "is_gvcf_header",
    "normalize_vcf_version",
    "find_first_vcf_with_header",
    "validate_vcf",
    "validate_all_vcfs",
    "validate_merged_vcf",
]
