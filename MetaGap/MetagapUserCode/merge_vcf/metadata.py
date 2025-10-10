"""Metadata helpers for augmenting merged VCF files."""

from __future__ import annotations

import logging
import copy
import gzip
import re
from contextlib import ExitStack, closing
from pathlib import Path
from collections import OrderedDict
from typing import List, Optional, Tuple, Union

import logging
import subprocess

from . import PYSAM_AVAILABLE, VCFPY_AVAILABLE, pysam, vcfpy
from .logging_utils import (
    MergeConflictError,
    ValidationError,
    handle_critical_error,
    log_message,
)

try:  # pragma: no cover - optional dependency
    import pysam
except ImportError:  # pragma: no cover - exercised in environments without pysam
    pysam = None


STANDARD_INFO_DEFINITIONS = OrderedDict(
    [
        (
            "AC",
            OrderedDict(
                [
                    ("Number", "A"),
                    ("Type", "Integer"),
                    (
                        "Description",
                        "Alternate allele count in genotypes, for each ALT allele",
                    ),
                ]
            ),
        ),
        (
            "AN",
            OrderedDict(
                [
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Total number of alleles in called genotypes"),
                ]
            ),
        ),
        (
            "AF",
            OrderedDict(
                [
                    ("Number", "A"),
                    ("Type", "Float"),
                    ("Description", "Alternate allele frequency"),
                ]
            ),
        ),
    ]
)


INFO_HEADER_PATTERN = re.compile(r"^##INFO=<ID=([^,>]+)")
FILTER_HEADER_PATTERN = re.compile(r"^##FILTER=<ID=([^,>]+)")
CONTIG_HEADER_PATTERN = re.compile(r"^##contig=<ID=([^,>]+)", re.IGNORECASE)


def _format_info_definition(info_id: str, definition_mapping: OrderedDict) -> str:
    parts = [f"ID={info_id}"]
    for key, value in definition_mapping.items():
        if value is None:
            continue
        if key == "Description":
            escaped_value = str(value).replace('"', '\\"')
            parts.append(f'Description="{escaped_value}"')
            continue
        parts.append(f"{key}={value}")
    return "##INFO=<" + ",".join(parts) + ">"


def ensure_standard_info_definitions(header, verbose: bool = False):
    if not VCFPY_AVAILABLE:
        return header

    info_module = getattr(vcfpy, "header", None)
    info_cls = getattr(info_module, "InfoHeaderLine", None) if info_module else None
    if info_cls is None:
        return header

    existing_ids = set()
    for line in getattr(header, "lines", []):
        if isinstance(line, info_cls):
            existing_ids.add(getattr(line, "id", None))

    added_ids = []
    for info_id, definition in STANDARD_INFO_DEFINITIONS.items():
        if info_id in existing_ids:
            continue

        mapping = OrderedDict([("ID", info_id)])
        mapping.update(definition)

        try:
            info_line = info_cls.from_mapping(mapping)
        except Exception as exc:  # pragma: no cover - defensive logging
            log_message(
                f"Unable to create INFO header definition for {info_id}: {exc}",
                verbose,
                level=logging.WARNING,
            )
            continue

        try:
            header.add_line(info_line)
        except Exception as exc:  # pragma: no cover - defensive logging
            log_message(
                f"Failed to append INFO header definition for {info_id}: {exc}",
                verbose,
                level=logging.WARNING,
            )
            continue

        existing_ids.add(info_id)
        added_ids.append(info_id)

    if added_ids:
        log_message(
            "Inserted INFO header definitions for missing fields: "
            + ", ".join(added_ids),
            verbose,
            level=logging.INFO,
        )

    return header


def ensure_standard_info_header_lines(
    final_header_lines, existing_header_lines, verbose: bool = False
):
    current_ids = set()
    for line in final_header_lines:
        if not isinstance(line, str):
            continue
        match = INFO_HEADER_PATTERN.match(line.strip())
        if match:
            current_ids.add(match.group(1))

    added_ids = []
    for info_id, definition in STANDARD_INFO_DEFINITIONS.items():
        if info_id in current_ids:
            continue

        formatted = _format_info_definition(info_id, definition)
        if formatted in existing_header_lines:
            current_ids.add(info_id)
            continue

        final_header_lines.append(formatted)
        existing_header_lines.add(formatted)
        current_ids.add(info_id)
        added_ids.append(info_id)

    if added_ids:
        log_message(
            "Inserted INFO header definitions for missing fields: "
            + ", ".join(added_ids),
            verbose,
            level=logging.INFO,
        )

    return final_header_lines


def _format_sample_metadata_value(value: str) -> str:
    """Return a VCF-safe representation of the provided metadata value."""

    value = str(value)
    needs_quotes = any(char in value for char in [" ", ",", "\t", "\"", "<", ">", "="])
    if not needs_quotes:
        return value

    escaped = value.replace("\\", "\\\\").replace("\"", "\\\"")
    return f'"{escaped}"'


def build_sample_metadata_line(entries: "OrderedDict[str, str]") -> str:
    """Serialize an ordered mapping into a single ``##SAMPLE`` metadata line."""

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


def _parse_sample_metadata_line(serialized: str) -> "OrderedDict[str, str]":
    """Return an ordered mapping extracted from a serialized ``##SAMPLE`` line."""

    if not isinstance(serialized, str):
        raise TypeError("Serialized sample metadata must be provided as a string.")

    text = serialized.strip()
    prefix = "##SAMPLE=<"
    suffix = ">"
    if not text.startswith(prefix) or not text.endswith(suffix):
        raise ValueError("Serialized sample metadata must be in '##SAMPLE=<...>' format.")

    body = text[len(prefix) : -len(suffix)]
    entries = OrderedDict()

    token: List[str] = []
    stack: List[str] = []
    in_quotes = False
    escape = False

    def flush_token() -> None:
        raw = "".join(token).strip()
        token.clear()
        if not raw:
            return
        if "=" not in raw:
            raise ValueError(f"Invalid SAMPLE metadata entry: '{raw}'")
        key, value = raw.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            raise ValueError("SAMPLE metadata keys cannot be empty.")
        if len(value) >= 2 and value[0] == value[-1] == '"':
            inner = value[1:-1]
            unescaped: List[str] = []
            i = 0
            while i < len(inner):
                ch = inner[i]
                if ch == "\\" and i + 1 < len(inner):
                    next_ch = inner[i + 1]
                    if next_ch in {'\\', '"'}:
                        unescaped.append(next_ch)
                        i += 2
                        continue
                unescaped.append(ch)
                i += 1
            value_to_store = "".join(unescaped)
        else:
            value_to_store = value
        entries[key] = value_to_store

    for ch in body:
        if escape:
            token.append(ch)
            escape = False
            continue

        if ch == "\\":
            token.append(ch)
            escape = True
            continue

        if in_quotes:
            if ch == '"':
                in_quotes = False
            token.append(ch)
            continue

        if ch == '"':
            in_quotes = True
            token.append(ch)
            continue

        if ch in "{[":
            stack.append(ch)
            token.append(ch)
            continue

        if ch in "}]":
            if stack:
                opener = stack[-1]
                if (opener == "{" and ch == "}") or (opener == "[" and ch == "]"):
                    stack.pop()
            token.append(ch)
            continue

        if ch == "," and not stack:
            flush_token()
            continue

        token.append(ch)

    flush_token()

    if "ID" in entries and isinstance(entries["ID"], str):
        entries["ID"] = entries["ID"].strip()
    return entries


def _parse_simple_metadata_line(line: str) -> Optional[Tuple[str, str]]:
    stripped = line.strip()
    if not stripped.startswith("##") or "=" not in stripped:
        return None
    key, value = stripped[2:].split("=", 1)
    key = key.strip()
    value = value.strip()
    if not key:
        return None
    return key, value


def _load_metadata_template(
    template_path: Optional[str], verbose: bool = False
) -> Tuple[Optional["OrderedDict[str, str]"], List, List[str], Optional[str]]:
    """Return metadata derived from a user-supplied template header file."""

    if not template_path:
        return None, [], [], None

    normalized_path = Path(template_path).expanduser().resolve()
    if not normalized_path.exists():
        handle_critical_error(
            f"Metadata template header file does not exist: {template_path}",
            exc_cls=ValidationError,
        )
    if not normalized_path.is_file():
        handle_critical_error(
            f"Metadata template header path is not a file: {template_path}",
            exc_cls=ValidationError,
        )

    template_sample_mapping: Optional["OrderedDict[str, str]"] = None
    template_serialized_sample: Optional[str] = None
    sanitized_lines: List[str] = []
    simple_header_lines = []
    existing_simple = set()
    existing_sanitized = set()

    try:
        with normalized_path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#CHROM"):
                    continue
                if not stripped.startswith("##"):
                    log_message(
                        (
                            "Skipping non-metadata line in template "
                            f"{normalized_path} line {line_number}: {stripped}"
                        ),
                        verbose,
                        level=logging.DEBUG,
                    )
                    continue
                if stripped.startswith("##SAMPLE="):
                    try:
                        parsed_sample = _parse_sample_metadata_line(stripped)
                    except ValueError as exc:
                        handle_critical_error(
                            "Invalid SAMPLE metadata line in template "
                            f"{normalized_path} line {line_number}: {exc}",
                            exc_cls=ValidationError,
                        )
                    if template_sample_mapping is not None:
                        log_message(
                            (
                                "Multiple ##SAMPLE lines found in metadata "
                                f"template {normalized_path}; using the last occurrence"
                            ),
                            verbose,
                            level=logging.WARNING,
                        )
                    template_sample_mapping = OrderedDict(parsed_sample.items())
                    template_serialized_sample = stripped
                    continue

                parsed = _parse_simple_metadata_line(stripped)
                if not parsed:
                    handle_critical_error(
                        "Metadata template lines must be in '##key=value' format. "
                        f"Problematic entry at {normalized_path} line {line_number}: {stripped}",
                        exc_cls=ValidationError,
                    )
                key, value = parsed
                sanitized = f"##{key}={value}"
                if sanitized not in existing_sanitized:
                    sanitized_lines.append(sanitized)
                    existing_sanitized.add(sanitized)

                try:
                    simple_line = vcfpy.SimpleHeaderLine(key, value)
                except Exception as exc:  # pragma: no cover - defensive logging
                    log_message(
                        f"Failed to parse metadata template line '{sanitized}': {exc}",
                        verbose,
                        level=logging.WARNING,
                    )
                    continue

                identifier = (simple_line.key, getattr(simple_line, "value", None))
                if identifier in existing_simple:
                    continue
                simple_header_lines.append(simple_line)
                existing_simple.add(identifier)
    except OSError as exc:
        handle_critical_error(
            f"Unable to read metadata template header file {template_path}: {exc}",
            exc_cls=ValidationError,
        )

    log_message(
        f"Loaded metadata template header from {normalized_path}",
        verbose,
        level=logging.INFO,
    )
    if sanitized_lines:
        log_message(
            "Metadata template header lines: " + ", ".join(sanitized_lines),
            verbose,
            level=logging.DEBUG,
        )
    if template_serialized_sample:
        log_message(
            f"Metadata template sample line: {template_serialized_sample}",
            verbose,
            level=logging.DEBUG,
        )

    return (
        template_sample_mapping,
        simple_header_lines,
        sanitized_lines,
        template_serialized_sample,
    )


def load_metadata_lines(metadata_file: str, verbose: bool = False) -> List[str]:
    """Return sanitized ``##key=value`` metadata lines from ``metadata_file``."""

    if not metadata_file:
        handle_critical_error(
            "A metadata file path must be provided.",
            exc_cls=ValidationError,
        )

    normalized_path = Path(metadata_file).expanduser().resolve()
    if not normalized_path.exists():
        handle_critical_error(
            f"Metadata file does not exist: {metadata_file}",
            exc_cls=ValidationError,
        )
    if not normalized_path.is_file():
        handle_critical_error(
            f"Metadata file path is not a file: {metadata_file}",
            exc_cls=ValidationError,
        )

    sanitized_lines: List[str] = []
    seen_lines = set()

    try:
        with normalized_path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                stripped = raw_line.strip()
                if not stripped:
                    continue
                if not stripped.startswith("##"):
                    stripped = "##" + stripped
                parsed = _parse_simple_metadata_line(stripped)
                if not parsed:
                    handle_critical_error(
                        "Metadata file lines must be in '##key=value' format. "
                        f"Problematic entry at {normalized_path} line {line_number}: {raw_line.strip()}",
                        exc_cls=ValidationError,
                    )
                key, value = parsed
                sanitized = f"##{key}={value}"
                if sanitized in seen_lines:
                    continue
                sanitized_lines.append(sanitized)
                seen_lines.add(sanitized)
    except OSError as exc:
        handle_critical_error(
            f"Unable to read metadata file {metadata_file}: {exc}",
            exc_cls=ValidationError,
        )

    log_message(
        f"Loaded {len(sanitized_lines)} metadata header line(s) from {normalized_path}",
        verbose,
        level=logging.INFO,
    )

    return sanitized_lines


def apply_metadata_to_header(
    header,
    sample_header_line=None,
    simple_header_lines=None,
    verbose: bool = False,
):
    """Return ``header`` with CLI metadata applied."""

    if sample_header_line is None and not simple_header_lines:
        return header

    simple_header_lines = simple_header_lines or []

    if simple_header_lines:
        existing_simple = {
            (line.key, getattr(line, "value", None))
            for line in header.lines
            if isinstance(line, vcfpy.SimpleHeaderLine)
        }
        for simple_line in simple_header_lines:
            identifier = (simple_line.key, simple_line.value)
            if identifier in existing_simple:
                continue
            header.add_line(simple_line)
            existing_simple.add(identifier)

    if sample_header_line is not None:
        header.lines = [
            line for line in header.lines if getattr(line, "key", None) != "SAMPLE"
        ]
        header.add_line(sample_header_line)

    header = ensure_standard_info_definitions(header, verbose=verbose)

    log_message("Applied CLI metadata to merged header.", verbose, level=logging.INFO)
    return header


def _validate_anonymized_vcf_header(
    final_vcf_path: Union[str, Path], ensure_for_uncompressed: bool = False
):
    """Ensure the anonymized VCF header can be parsed by ``vcfpy``."""

    if not final_vcf_path:
        return

    resolved_path = Path(final_vcf_path)
    suffixes = resolved_path.suffixes
    requires_validation = ensure_for_uncompressed or suffixes[-2:] == [
        ".vcf",
        ".gz",
    ]
    if not requires_validation:
        return

    opener = gzip.open if resolved_path.suffix == ".gz" else open
    try:
        with opener(resolved_path, "rt", encoding="utf-8") as stream, closing(
            vcfpy.Reader.from_stream(stream)
        ) as reader:
            _ = reader.header
    except Exception as exc:
        handle_critical_error(
            f"Failed to read header from anonymized VCF using vcfpy: {exc}"
        )


def append_metadata_to_merged_vcf(
    merged_vcf: str,
    sample_metadata_entries=None,
    header_metadata_lines=None,
    serialized_sample_line=None,
    verbose: bool = False,
):
    """Finalize the merged VCF by applying in-Python processing and metadata."""

    if not PYSAM_AVAILABLE or not VCFPY_AVAILABLE:  # pragma: no cover - defensive
        handle_critical_error(
            "vcfpy and pysam must be available to finalize the merged VCF output."
        )

    header_metadata_lines = header_metadata_lines or []
    merged_path = Path(merged_vcf)
    joint_temp = merged_path.with_name(f"{merged_path.name}.joint.temp.vcf")
    filtered_temp = merged_path.with_name(f"{merged_path.name}.filtered.temp.vcf")

    expects_gzip = merged_path.suffix == ".gz"
    final_plain_vcf_path: Optional[Path] = None
    final_vcf_path: Optional[Path] = None

    name = merged_path.name
    parent = merged_path.parent
    if name.endswith(".vcf.gz"):
        base_name = name[: -len(".vcf.gz")]
        final_plain_vcf_path = parent / f"{base_name}.anonymized.vcf"
        final_vcf_path = final_plain_vcf_path.with_name(
            f"{final_plain_vcf_path.name}.gz"
        )
    elif name.endswith(".vcf"):
        base_name = name[: -len(".vcf")]
        final_plain_vcf_path = parent / f"{base_name}.anonymized.vcf"
        final_vcf_path = final_plain_vcf_path
    elif name.endswith(".gz"):
        base_name = name[: -len(".gz")]
        final_plain_vcf_path = parent / f"{base_name}.anonymized"
        final_vcf_path = final_plain_vcf_path.with_name(
            f"{final_plain_vcf_path.name}.gz"
        )
    else:
        extension = merged_path.suffix or ".vcf"
        base_name = name[: -len(extension)] if merged_path.suffix else name
        final_plain_vcf_path = parent / f"{base_name}.anonymized{extension}"
        final_vcf_path = final_plain_vcf_path

    def _open_vcf(path: Path):
        if path.suffix == ".gz":
            return gzip.open(path, "rt", encoding="utf-8")
        return path.open("r", encoding="utf-8")

    def _read_header_lines(path: Path) -> List[str]:
        header_lines: List[str] = []
        with _open_vcf(path) as handle:
            for raw in handle:
                if not raw.startswith("#"):
                    break
                header_lines.append(raw.rstrip("\n"))
        return header_lines

    def _format_scalar(value) -> str:
        if value is None:
            return "."
        if isinstance(value, float):
            return f"{value:g}"
        return str(value)

    def _serialize_info(info: OrderedDict) -> str:
        entries = []
        for key, value in info.items():
            if value is None:
                continue
            if isinstance(value, list):
                if not value:
                    serialized = "."
                else:
                    serialized = ",".join(_format_scalar(item) for item in value)
            else:
                serialized = _format_scalar(value)
            entries.append(f"{key}={serialized}")
        return ";".join(entries) if entries else "."

    def _normalize_filter(field) -> str:
        if field in {None, [], ["PASS"], "PASS"}:
            return "PASS"
        if isinstance(field, list):
            return ";".join(str(entry) for entry in field)
        return str(field)

    def _normalize_quality(qual_value):
        if qual_value in {None, ".", ""}:
            return None
        try:
            return float(qual_value)
        except (TypeError, ValueError):
            return None

    def _extract_called_alleles(gt_value) -> List[int]:
        if gt_value is None:
            return []
        if isinstance(gt_value, (list, tuple)):
            tokens = []
            for item in gt_value:
                if item is None:
                    continue
                tokens.extend(str(item).replace("|", "/").split("/"))
        else:
            tokens = str(gt_value).replace("|", "/").split("/")
        allele_indices: List[int] = []
        for token in tokens:
            cleaned = token.strip()
            if not cleaned or cleaned == ".":
                continue
            try:
                allele_indices.append(int(cleaned))
            except ValueError:
                continue
        return allele_indices

    def _cleanup_temp_files():
        for temp_path in (joint_temp, filtered_temp):
            try:
                temp_path.unlink()
            except FileNotFoundError:
                continue
            except Exception:
                pass

    def _recalculate_info_fields(record) -> int:
        alt_count = len(getattr(record, "ALT", []) or [])
        ac = [0] * alt_count
        an = 0
        for call in getattr(record, "calls", []):
            data = getattr(call, "data", {})
            alleles = _extract_called_alleles(data.get("GT"))
            for allele in alleles:
                if allele is None:
                    continue
                if allele >= 0:
                    an += 1
                if 1 <= allele <= alt_count:
                    ac[allele - 1] += 1
        info = getattr(record, "INFO", OrderedDict())
        info["AC"] = ac if alt_count else []
        info["AN"] = an
        info["AF"] = [count / an for count in ac] if an else []
        record.INFO = info
        return an

    def _record_passes_filters(record, allele_number: int) -> bool:
        qual_value = _normalize_quality(getattr(record, "QUAL", None))
        if qual_value is None or qual_value <= 30:
            return False
        if allele_number <= 50:
            return False
        filters = getattr(record, "FILTER", None)
        filter_text = _normalize_filter(filters)
        return filter_text == "PASS"

    def _serialize_record(record) -> str:
        alt_field = (
            ",".join(str(alt) for alt in getattr(record, "ALT", []) or [])
            if getattr(record, "ALT", [])
            else "."
        )
        qual_value = _normalize_quality(getattr(record, "QUAL", None))
        qual_field = _format_scalar(qual_value) if qual_value is not None else "."
        info_field = _serialize_info(getattr(record, "INFO", OrderedDict()))
        filter_field = _normalize_filter(getattr(record, "FILTER", None))
        record_id = getattr(record, "ID", None)
        if isinstance(record_id, list):
            record_id_field = ";".join(str(entry) for entry in record_id)
        else:
            record_id_field = str(record_id or ".")
        return "\t".join(
            [
                str(getattr(record, "CHROM", ".")),
                str(getattr(record, "POS", ".")),
                record_id_field,
                str(getattr(record, "REF", ".")),
                alt_field,
                qual_field,
                filter_field,
                info_field,
            ]
        )

    log_message(
        "Recalculating AC, AN, AF tags across the merged cohort...",
        verbose,
        level=logging.DEBUG,
    )

    header_lines = _read_header_lines(merged_path)
    fileformat_line = None
    remaining_header_lines: List[str] = []
    for line in header_lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#CHROM"):
            continue
        if stripped.startswith("##fileformat") and fileformat_line is None:
            fileformat_line = stripped
            continue
        if stripped.startswith("##SAMPLE="):
            continue
        remaining_header_lines.append(stripped)

    final_header_lines: List[str] = []
    if fileformat_line is not None:
        final_header_lines.append(fileformat_line)
    final_header_lines.extend(remaining_header_lines)
    existing_header_lines = set(final_header_lines)

    ensure_standard_info_header_lines(
        final_header_lines, existing_header_lines, verbose=verbose
    )

    if sample_metadata_entries:
        try:
            serialized_line = (
                serialized_sample_line
                if serialized_sample_line is not None
                else build_sample_metadata_line(sample_metadata_entries)
            )
            if serialized_line not in existing_header_lines:
                final_header_lines.append(serialized_line)
                existing_header_lines.add(serialized_line)
        except ValueError as exc:
            handle_critical_error(str(exc))

    for metadata_line in header_metadata_lines:
        normalized = metadata_line.strip()
        if not normalized:
            continue
        if not normalized.startswith("##"):
            normalized = "##" + normalized
        if normalized in existing_header_lines:
            continue
        final_header_lines.append(normalized)
        existing_header_lines.add(normalized)

    ensure_standard_info_header_lines(
        final_header_lines, existing_header_lines, verbose=verbose
    )

    final_header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    body_lines: List[str] = []
    with ExitStack() as stack:
        reader = None
        try:
            stream = stack.enter_context(_open_vcf(merged_path))
            reader = stack.enter_context(closing(vcfpy.Reader.from_stream(stream)))
        except NotImplementedError:
            reader = stack.enter_context(
                closing(vcfpy.Reader.from_path(str(merged_path)))
            )

        for record in reader:
            allele_number = _recalculate_info_fields(record)
            if not _record_passes_filters(record, allele_number):
                continue
            body_lines.append(_serialize_record(record))

    if final_plain_vcf_path is None:
        extension = merged_path.suffix or ".vcf"
        base_name = merged_path.stem if merged_path.suffix else merged_path.name
        final_plain_vcf_path = merged_path.with_name(
            f"{base_name}.anonymized{extension}"
        )

    try:
        with final_plain_vcf_path.open("w", encoding="utf-8") as output_handle:
            for line in final_header_lines:
                output_handle.write(line + "\n")
            for line in body_lines:
                output_handle.write(line + "\n")
        final_vcf_path = final_plain_vcf_path
    except Exception as exc:
        if final_plain_vcf_path and final_plain_vcf_path.exists():
            try:
                final_plain_vcf_path.unlink()
            except Exception:
                pass
        handle_critical_error(
            f"Failed to assemble final anonymized VCF contents: {exc}",
            exc_cls=MergeConflictError,
        )

    if expects_gzip and final_plain_vcf_path and final_plain_vcf_path.exists():
        if pysam is None:
            _cleanup_temp_files()
            handle_critical_error(
                "pysam is required to BGZF-compress and index the anonymized VCF."
            )

        try:
            final_vcf_path = final_plain_vcf_path.with_name(
                f"{final_plain_vcf_path.name}.gz"
            )
            pysam.bgzip(str(final_plain_vcf_path), str(final_vcf_path), force=True)
            pysam.tabix_index(str(final_vcf_path), preset="vcf", force=True)
        except (OSError, ValueError) as exc:
            _cleanup_temp_files()
            if final_plain_vcf_path.exists():
                try:
                    final_plain_vcf_path.unlink()
                except Exception:
                    pass
            if final_vcf_path and final_vcf_path.exists():
                try:
                    final_vcf_path.unlink()
                except Exception:
                    pass
            handle_critical_error(
                f"Failed to compress or index anonymized VCF: {exc}",
                exc_cls=MergeConflictError,
            )
    elif not expects_gzip:
        final_vcf_path = final_plain_vcf_path

    final_vcf_str = str(final_vcf_path) if final_vcf_path is not None else ""
    _validate_anonymized_vcf_header(final_vcf_str, ensure_for_uncompressed=True)

    log_message(
        f"Anonymized merged VCF written to: {final_vcf_str}",
        verbose,
        level=logging.INFO,
    )
    return final_vcf_str


def recalculate_cohort_info_tags(vcf_path: str, verbose: bool = False) -> None:
    """Recompute AC, AN, and AF INFO tags in-place for ``vcf_path``."""

    log_message(
        "Recalculating AC, AN, AF tags across the merged cohort...",
        verbose,
        level=logging.DEBUG,
    )

    if VCFPY_AVAILABLE:
        header = None
        records: List = []
        try:
            with closing(vcfpy.Reader.from_path(vcf_path)) as reader:
                header = reader.header

                for record in reader:
                    info = dict(getattr(record, "INFO", {}) or {})
                    alt_alleles = list(getattr(record, "ALT", []) or [])
                    alt_counts = [0] * len(alt_alleles)
                    allele_number = 0

                    for call in getattr(record, "calls", []) or []:
                        data = getattr(call, "data", {}) or {}
                        genotype = data.get("GT")
                        if not genotype:
                            continue
                        tokens = genotype.replace("|", "/").split("/")
                        for token in tokens:
                            token = token.strip()
                            if not token or token == ".":
                                continue
                            try:
                                allele_index = int(token)
                            except ValueError:
                                continue
                            if allele_index < 0:
                                continue
                            allele_number += 1
                            if allele_index == 0:
                                continue
                            normalized_index = allele_index - 1
                            if 0 <= normalized_index < len(alt_counts):
                                alt_counts[normalized_index] += 1

                    info["AN"] = allele_number
                    if alt_counts:
                        info["AC"] = alt_counts
                        if allele_number > 0:
                            info["AF"] = [
                                count / allele_number for count in alt_counts
                            ]
                        else:
                            info["AF"] = [0.0 for _ in alt_counts]
                    else:
                        info["AC"] = []
                        info["AF"] = []

                    record.INFO = info
                    records.append(record)
        except Exception as exc:
            handle_critical_error(
                f"Failed to open {vcf_path} for AC/AN/AF recalculation: {exc}",
                exc_cls=MergeConflictError,
            )

        try:
            with closing(vcfpy.Writer.from_path(vcf_path, header)) as writer:
                for record in records:
                    writer.write_record(record)
        except Exception as exc:
            handle_critical_error(
                f"Failed to reopen {vcf_path} for writing: {exc}",
                exc_cls=MergeConflictError,
            )
        return

    _recalculate_cohort_info_tags_without_vcfpy(vcf_path)


def _recalculate_cohort_info_tags_without_vcfpy(vcf_path: str) -> None:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    try:
        with opener(vcf_path, "rt", encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle]
    except OSError as exc:
        handle_critical_error(
            f"Failed to open {vcf_path}: {exc}",
            exc_cls=MergeConflictError,
        )

    header_lines = []
    records = []
    column_header = None
    for line in lines:
        if line.startswith("##"):
            header_lines.append(line)
            continue
        if line.startswith("#CHROM"):
            column_header = line
            header_lines.append(line)
            continue
        if not line:
            continue
        records.append(line)

    if column_header is None:
        handle_critical_error(
            f"Failed to parse {vcf_path}: Missing #CHROM header line.",
            exc_cls=MergeConflictError,
        )

    updated_records = []
    for record_line in records:
        fields = record_line.split("\t")
        if len(fields) < 8:
            updated_records.append(record_line)
            continue

        alt_alleles = [allele for allele in fields[4].split(",") if allele]
        alt_counts = [0] * len(alt_alleles)
        allele_number = 0

        if len(fields) >= 10:
            format_keys = fields[8].split(":") if fields[8] else []
            try:
                gt_index = format_keys.index("GT")
            except ValueError:
                gt_index = None

            if gt_index is not None:
                for sample_entry in fields[9:]:
                    parts = sample_entry.split(":")
                    if gt_index >= len(parts):
                        continue
                    genotype = parts[gt_index]
                    if not genotype:
                        continue
                    tokens = genotype.replace("|", "/").split("/")
                    for token in tokens:
                        token = token.strip()
                        if not token or token == ".":
                            continue
                        try:
                            allele_index = int(token)
                        except ValueError:
                            continue
                        if allele_index < 0:
                            continue
                        allele_number += 1
                        if allele_index == 0:
                            continue
                        normalized_index = allele_index - 1
                        if 0 <= normalized_index < len(alt_counts):
                            alt_counts[normalized_index] += 1

        info_field = fields[7]
        info_map = {}
        order = []
        for entry in info_field.split(";"):
            if not entry:
                continue
            if "=" in entry:
                key, value = entry.split("=", 1)
            else:
                key, value = entry, ""
            if key not in info_map:
                order.append(key)
            info_map[key] = value

        info_map["AN"] = str(allele_number)
        if alt_counts:
            info_map["AC"] = ",".join(str(count) for count in alt_counts)
            if allele_number > 0:
                info_map["AF"] = ",".join(
                    f"{count / allele_number:.6g}" for count in alt_counts
                )
            else:
                info_map["AF"] = ",".join("0" for _ in alt_counts)
        else:
            info_map["AC"] = "."
            info_map["AF"] = "."

        for key in ("AC", "AN", "AF"):
            if key not in order:
                order.append(key)

        new_entries = []
        for key in order:
            value = info_map.get(key, "")
            if value == "":
                new_entries.append(f"{key}=")
            else:
                new_entries.append(f"{key}={value}")
        fields[7] = ";".join(new_entries)
        updated_records.append("\t".join(fields))

    try:
        with opener(vcf_path, "wt", encoding="utf-8") as handle:
            for line in header_lines:
                handle.write(line + "\n")
            for line in updated_records:
                handle.write(line + "\n")
    except OSError as exc:
        handle_critical_error(
            f"Failed to rewrite {vcf_path}: {exc}",
            exc_cls=MergeConflictError,
        )


__all__ = [
    "_format_sample_metadata_value",
    "build_sample_metadata_line",
    "_parse_sample_metadata_line",
    "_parse_simple_metadata_line",
    "load_metadata_lines",
    "apply_metadata_to_header",
    "_validate_anonymized_vcf_header",
    "append_metadata_to_merged_vcf",
    "recalculate_cohort_info_tags",
]
