"""Metadata helpers for augmenting merged VCF files.

This module orchestrates the construction of VCF headers by combining CLI
metadata with auto-generated INFO definitions and helper utilities for reading
and writing VCF text streams. The helpers centralize tricky logic such as
escaping header descriptions and de-duplicating metadata contributed by
multiple sources so that downstream workflows can reuse the same guarantees.
"""

from __future__ import annotations

import collections
from collections.abc import MutableMapping
import gzip
import re
from contextlib import ExitStack, closing
from pathlib import Path
import io
import os
import re
import shutil
import tempfile
from collections import OrderedDict
from typing import Iterable, List, Optional, Tuple, Mapping, Union

import logging
import subprocess

from . import PYSAM_AVAILABLE, VCFPY_AVAILABLE, pysam, vcfpy
from .filtering import DEFAULT_ALLOWED_FILTER_VALUES, prepare_allowed_filter_set
from .logging_utils import (
    MergeConflictError,
    ValidationError,
    handle_critical_error,
    log_message,
)

DEFAULT_QUAL_THRESHOLD = 30.0
"""Minimum QUAL value required for a record to pass default filtering."""

DEFAULT_AN_THRESHOLD = 50.0
"""Minimum allele number required for a record to pass default filtering."""

try:  # pragma: no cover - optional dependency
    import pysam
except ImportError:  # pragma: no cover - exercised in environments without pysam
    pysam = None


DEFAULT_QUAL_THRESHOLD: Optional[float] = 30.0
DEFAULT_AN_THRESHOLD: Optional[float] = 50.0
DEFAULT_ALLOWED_FILTER_VALUES: Tuple[str, ...] = ("PASS",)


STANDARD_INFO_DEFINITIONS: dict[str, dict[str, str]] = {
    "AC": {
        "Number": "A",
        "Type": "Integer",
        "Description": "Alternate allele count in genotypes, for each ALT allele",
    },
    "AN": {
        "Number": "1",
        "Type": "Integer",
        "Description": "Total number of alleles in called genotypes",
    },
    "AF": {
        "Number": "A",
        "Type": "Float",
        "Description": "Alternate allele frequency",
    },
}


INFO_HEADER_PATTERN = re.compile(r"^##INFO=<ID=([^,>]+)")
FILTER_HEADER_PATTERN = re.compile(r"^##FILTER=<ID=([^,>]+)")
CONTIG_HEADER_PATTERN = re.compile(r"^##contig=<ID=([^,>]+)", re.IGNORECASE)


def _metadata_record_passes_filters(
    record,
    allele_number,
    qual_threshold: Optional[float],
    an_threshold: Optional[float],
    allowed_filter_values: Optional[Sequence[str]],
) -> bool:
    """Mirror merging._record_passes_filters using recalculated allele counts."""

    if qual_threshold is not None:
        q = getattr(record, "QUAL", None)
        try:
            if q is None or float(q) <= qual_threshold:
                return False
        except (TypeError, ValueError):
            return False

    if an_threshold is not None:
        an_val = allele_number
        if isinstance(an_val, list):
            an_val = an_val[0] if an_val else None
        try:
            if an_val is None or float(an_val) <= an_threshold:
                return False
        except (TypeError, ValueError):
            return False

    filters_field = getattr(record, "FILTER", None)
    filters = filters_field or []
    if not isinstance(filters, (list, tuple)):
        filters = [filters]
    vals = [f for f in filters if f not in {None, "", "."}]
    if not vals:
        return True

    allowed_set = prepare_allowed_filter_set(allowed_filter_values)
    return bool(allowed_set) and all(value in allowed_set for value in vals)


def _format_info_definition(info_id: str, definition_mapping: Mapping[str, object]) -> str:
    parts = [f"ID={info_id}"]
    for key, value in definition_mapping.items():
        if value is None:
            continue
        if key == "Description":
            # Escape embedded quotes so the serialized INFO line survives a
            # round-trip through vcfpy and pysam without corrupting the
            # comma-delimited payload defined by the VCF specification.
            escaped_value = str(value).replace('"', '\\"')
            parts.append(f'Description="{escaped_value}"')
            continue
        parts.append(f"{key}={value}")
    return "##INFO=<" + ",".join(parts) + ">"


def ensure_standard_info_definitions(header, verbose: bool = False):
    """Ensure AC/AN/AF INFO definitions exist on a ``vcfpy`` header object.

    The helper inspects the mutable ``header`` supplied by ``vcfpy`` and
    inserts definitions for the canonical AC, AN, and AF INFO tags when they are
    absent. All operations are performed defensively so that environments
    lacking ``vcfpy`` support can continue without modification.
    """

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

        mapping: dict[str, object] = {"ID": info_id}
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
    """Inject AC/AN/AF INFO lines into serialized header text when needed.

    ``final_header_lines`` contains the mutable list of header strings that will
    be written to disk, while ``existing_header_lines`` tracks the set of lines
    already present. The helper mirrors :func:`ensure_standard_info_definitions`
    but operates on raw strings so that code paths without ``vcfpy`` access can
    benefit from the same guarantees.
    """

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
            # Header lines supplied by CLI metadata or pre-existing files must
            # not be duplicated; mark the ID so later iterations skip it.
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


def _open_vcf(path: str):
    return (
        gzip.open(path, "rt", encoding="utf-8")
        if path.endswith(".gz")
        else open(path, "r", encoding="utf-8")
    )


def _read_header_lines(path: str) -> List[str]:
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
        return "%g" % value
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


def _serialize_record(record) -> str:
    alt_values = getattr(record, "ALT", []) or []
    alt_field = ",".join(str(alt) for alt in alt_values) if alt_values else "."
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


def _cleanup_temp_files(
    joint_temp: str, filtered_temp: str, verbose: bool = False
) -> None:
    for path in (joint_temp, filtered_temp):
        if not path:
            continue
        if not os.path.exists(path):
            continue
        try:
            os.remove(path)
        except OSError as exc:  # pragma: no cover - defensive logging
            log_message(
                f"Failed to remove temporary file {path}: {exc}",
                verbose,
                level=logging.WARNING,
            )


def _format_sample_metadata_value(value: str) -> str:
    """Return a VCF-safe representation of the provided metadata value."""

    value = str(value)
    needs_quotes = any(char in value for char in [" ", ",", "\t", "\"", "<", ">", "="])
    if not needs_quotes:
        return value

    escaped = value.replace("\\", "\\\\").replace("\"", "\\\"")
    return f'"{escaped}"'


def build_sample_metadata_line(entries: Mapping[str, str]) -> str:
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


def _parse_sample_metadata_line(serialized: str) -> dict[str, str]:
    """Return an ordered mapping extracted from a serialized ``##SAMPLE`` line."""

    if not isinstance(serialized, str):
        raise TypeError("Serialized sample metadata must be provided as a string.")

    text = serialized.strip()
    prefix = "##SAMPLE=<"
    suffix = ">"
    if not text.startswith(prefix) or not text.endswith(suffix):
        raise ValueError("Serialized sample metadata must be in '##SAMPLE=<...>' format.")

    body = text[len(prefix) : -len(suffix)]
    entries: dict[str, str] = {}

    token: list[str] = []
    stack: list[str] = []
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
            unescaped: list[str] = []
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


def _parse_simple_metadata_line(line: str) -> tuple[str, str] | None:
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
) -> Tuple[Optional[dict[str, str]], List, List[str], Optional[str]]:
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

    template_sample_mapping: Optional[dict[str, str]] = None
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
                        raise ValidationError(
                            "Invalid SAMPLE metadata line in template "
                            f"{normalized_path} line {line_number}: {exc}"
                        ) from exc
                    if template_sample_mapping is not None:
                        log_message(
                            (
                                "Multiple ##SAMPLE lines found in metadata "
                                f"template {normalized_path}; using the last occurrence"
                            ),
                            verbose,
                            level=logging.WARNING,
                        )
                    template_sample_mapping = dict(parsed_sample.items())
                    template_serialized_sample = stripped
                    continue

                parsed = _parse_simple_metadata_line(stripped)
                if not parsed:
                    raise ValidationError(
                        "Metadata template lines must be in '##key=value' format. "
                        f"Problematic entry at {normalized_path} line {line_number}: {stripped}"
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
        raise ValidationError(
            f"Unable to read metadata template header file {template_path}: {exc}"
        ) from exc

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


def _split_metadata_segments(raw_line: str) -> list[str]:
    """Split concatenated metadata segments into individual ``##`` records."""

    if raw_line is None:
        return []

    text = raw_line.strip()
    if not text:
        return []

    # Ensure the string begins with ``##`` before scanning for additional segments.
    if not text.startswith("##"):
        text = "##" + text.lstrip("#")

    segments: list[str] = []
    start = 0
    length = len(text)

    while start < length:
        if not text.startswith("##", start):
            # Residual text without a ``##`` prefix is coerced into a valid segment.
            remainder = text[start:].strip()
            if remainder:
                if not remainder.startswith("##"):
                    remainder = "##" + remainder.lstrip("#")
                segments.append(remainder)
            break

        next_start = text.find("##", start + 2)
        segment = text[start:] if next_start == -1 else text[start:next_start]
        segment = segment.strip()
        if segment:
            if not segment.startswith("##"):
                segment = "##" + segment.lstrip("#")
            segments.append(segment)

        if next_start == -1:
            break
        start = next_start

    return segments


def load_metadata_lines(metadata_file: str, verbose: bool = False) -> list[str]:
    """Return sanitized ``##key=value`` metadata lines from ``metadata_file``."""

    if not metadata_file:
        raise ValidationError("A metadata file path must be provided.")

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

    sanitized_lines: list[str] = []
    seen_lines = set()

    try:
        with normalized_path.open("r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, 1):
                segments = _split_metadata_segments(raw_line)
                for segment in segments:
                    parsed = _parse_simple_metadata_line(segment)
                    if not parsed:
                        raise ValidationError(
                            "Metadata file lines must be in '##key=value' format. "
                            f"Problematic entry at {normalized_path} line {line_number}: {raw_line.strip()}"
                        )
                    key, value = parsed
                    sanitized = f"##{key}={value}"
                    if sanitized in seen_lines:
                        continue
                    sanitized_lines.append(sanitized)
                    seen_lines.add(sanitized)
    except OSError as exc:
        raise ValidationError(
            f"Unable to read metadata file {metadata_file}: {exc}"
        ) from exc

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
        with opener(final_vcf_path, "rt", encoding="utf-8") as stream:
            try:
                reader = vcfpy.Reader.from_stream(stream)
            except NotImplementedError:
                reader = vcfpy.Reader.from_path(final_vcf_path)
            # Force header consumption
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
    sample_header_entries=None,
    simple_header_lines=None,
    *,
    qual_threshold: Optional[float] = DEFAULT_QUAL_THRESHOLD,
    an_threshold: Optional[float] = DEFAULT_AN_THRESHOLD,
    allowed_filter_values: Optional[Sequence[str]] = DEFAULT_ALLOWED_FILTER_VALUES,
    verbose: bool = False,
):
    """Finalize the merged VCF by applying in-Python processing and metadata."""

    if not PYSAM_AVAILABLE or not VCFPY_AVAILABLE:  # pragma: no cover - defensive
        handle_critical_error(
            "vcfpy and pysam must be available to finalize the merged VCF output."
        )

    if sample_metadata_entries is None and sample_header_entries is not None:
        sample_metadata_entries = sample_header_entries
    if header_metadata_lines is None and simple_header_lines is not None:
        header_metadata_lines = simple_header_lines

    header_metadata_lines = header_metadata_lines or []
    effective_allowed_filters: Sequence[str] = (
        allowed_filter_values if allowed_filter_values is not None else DEFAULT_ALLOWED_FILTER_VALUES
    )
    joint_temp = f"{merged_vcf}.joint.temp.vcf"
    filtered_temp = f"{merged_vcf}.filtered.temp.vcf"

    expects_gzip = merged_vcf.endswith(".gz")
    final_plain_vcf: Optional[str]
    final_vcf: Optional[str]

    if merged_vcf.endswith(".vcf.gz"):
        base_path = merged_vcf[: -len(".vcf.gz")]
        final_plain_vcf = f"{base_path}.anonymized.vcf"
        final_vcf = f"{final_plain_vcf}.gz"
    elif merged_vcf.endswith(".vcf"):
        base_path = merged_vcf[: -len(".vcf")]
        final_plain_vcf = f"{base_path}.anonymized.vcf"
        final_vcf = final_plain_vcf
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

    def _ensure_standard_info_lines(lines: list[str]) -> None:
        need = {
            "AC": '##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count in genotypes, for each ALT allele">',
            "AN": '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
            "AF": '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, for each ALT allele">',
        }
        present = set()
        for line in lines:
            if line.startswith("##INFO=<ID="):
                try:
                    kid = line.split("##INFO=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
                    present.add(kid)
                except Exception:
                    continue
        for key, value in need.items():
            if key not in present:
                lines.append(value)

    header_lines = _read_header_lines(Path(merged_vcf))
    fileformat_line = None
    remaining_header_lines: list[str] = []
    for line in header_lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#CHROM"):
            continue
        if stripped.startswith("##fileformat") and fileformat_line is None:
            fileformat_line = stripped
            continue
        if stripped.startswith("##SAMPLE="):
            continue  # drop per-sample meta; we'll add cohort SAMPLE
        remaining_header_lines.append(stripped)

    final_header_lines: list[str] = [
        fileformat_line or "##fileformat=VCFv4.2",
        *remaining_header_lines,
    ]
    existing = set(final_header_lines)

    _ensure_standard_info_lines(final_header_lines)

    if sample_metadata_entries or serialized_sample_line:
        try:
            smp = serialized_sample_line or build_sample_metadata_line(sample_metadata_entries)
            if not smp.startswith("##"):
                smp = "##" + smp
            if smp not in existing:
                final_header_lines.append(smp)
                existing.add(smp)
        except ValueError as exc:
            handle_critical_error(str(exc))

    normalized_metadata: list[str] = []
    for meta in header_metadata_lines:
        for segment in _split_metadata_segments(meta):
            normalized_metadata.append(segment)

    deduplicated_metadata: list[str] = []
    seen_normalized: set[str] = set()
    for nm in normalized_metadata:
        cleaned = (nm or "").strip()
        if not cleaned:
            continue
        if not cleaned.startswith("##"):
            cleaned = "##" + cleaned.lstrip("#")
        if cleaned in existing or cleaned in seen_normalized:
            continue
        deduplicated_metadata.append(cleaned)
        seen_normalized.add(cleaned)

    final_header_lines.extend(deduplicated_metadata)
    existing.update(seen_normalized)

    _ensure_standard_info_lines(final_header_lines)
    final_header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # ----- Serializer for 8-column anonymized records -----
    def _fmt_filter(filt) -> str:
        if not filt:
            return "."
        if isinstance(filt, (list, tuple)):
            if not filt:
                return "."
            if filt == ["PASS"] or filt == ("PASS",):
                return "PASS"
            return ",".join(str(x) for x in filt)
        return str(filt)

    def _fmt_info(info: dict) -> str:
        if not info:
            return "."
        parts = []
        for key, value in info.items():
            if value is None:
                parts.append(key)
            elif isinstance(value, (list, tuple)):
                parts.append(
                    f"{key}=" + ",".join("" if entry is None else str(entry) for entry in value)
                )
            else:
                parts.append(f"{key}={value}")
        return ";".join(parts) if parts else "."

    def _serialize_record_8col(rec) -> str:
        chrom = str(getattr(rec, "CHROM", "."))
        pos = str(getattr(rec, "POS", "."))
        rid = getattr(rec, "ID", None)
        rid_field = (
            ";".join(rid)
            if isinstance(rid, (list, tuple))
            else (rid if rid else ".")
        )
        ref = str(getattr(rec, "REF", "."))
        alts = getattr(rec, "ALT", []) or []
        alt_field = ",".join(str(a) for a in alts) if alts else "."
        qual = getattr(rec, "QUAL", None)
        qual_field = "." if qual is None else str(qual)
        filt_field = _fmt_filter(getattr(rec, "FILTER", None))
        info_field = _fmt_info(getattr(rec, "INFO", {}) or {})
        return "\t".join(
            [chrom, pos, rid_field, ref, alt_field, qual_field, filt_field, info_field]
        )

    def _format_scalar(value) -> str:
        if value is None:
            return "."
        if isinstance(value, float):
            return f"{value:g}"
        return str(value)

    def _serialize_info(info: Mapping[str, object]) -> str:
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
        for path in [joint_temp, filtered_temp]:
            if not path:
                continue
            try:
                if os.path.exists(path):
                    os.remove(path)
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
        info = getattr(record, "INFO", {})
        info["AC"] = ac if alt_count else []
        info["AN"] = an
        info["AF"] = af if an else []
        record.INFO = info
        return an

    def _record_passes_filters(record, allele_number: int) -> bool:
        """Return whether a record satisfies the default QUAL and AN thresholds."""

        qual_value = _normalize_quality(getattr(record, "QUAL", None))
        if qual_threshold is not None:
            if qual_value is None or qual_value <= qual_threshold:
                return False
        if an_threshold is not None:
            if allele_number is None or allele_number <= an_threshold:
                return False

        filters = getattr(record, "FILTER", None)

        if allowed_filter_values is None:
            return True

        if filters in {None, [], (), "", "."}:
            return True

        if isinstance(filters, (list, tuple)):
            filter_tokens = [token for token in filters if token not in {None, "", "."}]
        else:
            filter_tokens = [filters] if filters not in {None, "", "."} else []

        if not filter_tokens:
            return True

        allowed = {
            str(value)
            for value in allowed_filter_values
            if value not in {None, "", "."}
        }

        if not allowed:
            return False

        return all(str(token) in allowed for token in filter_tokens)

    def _serialize_record(record) -> str:
        alt_field = (
            ",".join(str(alt) for alt in getattr(record, "ALT", []) or [])
            if getattr(record, "ALT", [])
            else "."
        )
        qual_value = _normalize_quality(getattr(record, "QUAL", None))
        qual_field = _format_scalar(qual_value) if qual_value is not None else "."
        info_field = _serialize_info(getattr(record, "INFO", {}))
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

    def append_metadata_to_merged_vcf(
    merged_vcf: str,
    sample_metadata_entries=None,
    header_metadata_lines=None,
    serialized_sample_line=None,
    qual_threshold: Optional[float] = DEFAULT_QUAL_THRESHOLD,
    an_threshold: Optional[float] = DEFAULT_AN_THRESHOLD,
    allowed_filter_values: Optional[Sequence[str]] = DEFAULT_ALLOWED_FILTER_VALUES,
    verbose: bool = False,
    ):
        """Finalize the merged VCF by applying metadata, AC/AN/AF recomputation, filtering, anonymization, and optional BGZF+Tabix."""
        if not PYSAM_AVAILABLE or not VCFPY_AVAILABLE:  # pragma: no cover
            handle_critical_error("vcfpy and pysam must be available to finalize the merged VCF output.")

        header_metadata_lines = header_metadata_lines or []

        # Decide final output names from input suffix.
        expects_gzip = merged_vcf.endswith(".gz")
        if merged_vcf.endswith(".vcf.gz"):
            base = merged_vcf[: -len(".vcf.gz")]
            final_plain_vcf = f"{base}.anonymized.vcf"
            final_vcf = f"{final_plain_vcf}.gz"
        elif merged_vcf.endswith(".vcf"):
            base = merged_vcf[: -len(".vcf")]
            final_plain_vcf = f"{base}.anonymized.vcf"
            final_vcf = final_plain_vcf
        elif merged_vcf.endswith(".gz"):
            base = merged_vcf[: -len(".gz")]
            final_plain_vcf = f"{base}.anonymized.vcf"
            final_vcf = f"{final_plain_vcf}.gz"
        else:
            final_plain_vcf = f"{merged_vcf}.anonymized.vcf"
            final_vcf = final_plain_vcf

    # ----- Build final header (line-based) -----
    def _read_header_lines(path: str) -> list[str]:
        op = gzip.open if path.endswith(".gz") else open
        out = []
        with op(path, "rt", encoding="utf-8") as fh:
            for ln in fh:
                if ln.startswith("#"):
                    out.append(ln.rstrip("\n"))
                    if ln.startswith("#CHROM"):
                        break
                else:
                    break
        return out

    def _ensure_standard_info_lines(lines: list[str]) -> None:
        need = {
            "AC": '##INFO=<ID=AC,Number=A,Type=Integer,Description="Alternate allele count in genotypes, for each ALT allele">',
            "AN": '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">',
            "AF": '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency, for each ALT allele">',
        }
        present = set()
        for l in lines:
            if l.startswith("##INFO=<ID="):
                try:
                    kid = l.split("##INFO=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
                    present.add(kid)
                except Exception:
                    pass
        for k, v in need.items():
            if k not in present:
                lines.append(v)

    header_lines = _read_header_lines(merged_vcf)
    fileformat_line = None
    remaining_header_lines: list[str] = []
    for line in header_lines:
        s = line.strip()
        if not s or s.startswith("#CHROM"):
            continue
        if s.startswith("##fileformat") and fileformat_line is None:
            fileformat_line = s
            continue
        if s.startswith("##SAMPLE="):
            continue  # drop per-sample meta; we'll add cohort SAMPLE
        remaining_header_lines.append(s)

    final_header_lines: list[str] = [fileformat_line or "##fileformat=VCFv4.2", *remaining_header_lines]
    existing = set(final_header_lines)

    _ensure_standard_info_lines(final_header_lines)

    if sample_metadata_entries or serialized_sample_line:
        try:
            smp = serialized_sample_line or build_sample_metadata_line(sample_metadata_entries)
            if not smp.startswith("##"):
                smp = "##" + smp
            if smp not in existing:
                final_header_lines.append(smp)
                existing.add(smp)
        except ValueError as exc:
            handle_critical_error(str(exc))

    for meta in header_metadata_lines:
        nm = (meta or "").strip()
        if not nm:
            continue
        if not nm.startswith("##"):
            nm = "##" + nm
        if nm in existing:
            continue
        final_header_lines.append(nm)
        existing.add(nm)

    _ensure_standard_info_lines(final_header_lines)
    final_header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    # ----- Serializer for 8-column anonymized records -----
    def _fmt_filter(filt) -> str:
        if not filt:
            return "."
        if isinstance(filt, (list, tuple)):
            if not filt:
                return "."
            if filt == ["PASS"] or filt == ("PASS",):
                return "PASS"
            return ",".join(str(x) for x in filt)
        return str(filt)

    def _fmt_info(info: dict) -> str:
        if not info:
            return "."
        parts = []
        for k, v in info.items():
            if v is None:
                parts.append(k)
            elif isinstance(v, (list, tuple)):
                parts.append(f"{k}=" + ",".join("" if vv is None else str(vv) for vv in v))
            else:
                parts.append(f"{k}={v}")
        return ";".join(parts) if parts else "."

    def _serialize_record_8col(rec) -> str:
        chrom = str(getattr(rec, "CHROM", "."))
        pos = str(getattr(rec, "POS", "."))
        rid = getattr(rec, "ID", None)
        rid_field = ";".join(rid) if isinstance(rid, (list, tuple)) else (rid if rid else ".")
        ref = str(getattr(rec, "REF", "."))
        alts = getattr(rec, "ALT", []) or []
        alt_field = ",".join(str(a) for a in alts) if alts else "."
        qual = getattr(rec, "QUAL", None)
        qual_field = "." if qual is None else str(qual)
        filt_field = _fmt_filter(getattr(rec, "FILTER", None))
        info_field = _fmt_info(getattr(rec, "INFO", {}) or {})
        return "\t".join([chrom, pos, rid_field, ref, alt_field, qual_field, filt_field, info_field])

    # ----- Stream, recompute AC/AN/AF, filter, anonymize, write -----
    reader = None
    temp_out: Optional[str] = None
    try:
        import tempfile
        out_dir = os.path.dirname(os.path.abspath(final_plain_vcf)) or None
        with tempfile.NamedTemporaryFile("w", encoding="utf-8", delete=False, dir=out_dir, suffix=".tmp") as outfh:
            temp_out = outfh.name
            # header
            for ln in final_header_lines:
                outfh.write(ln + "\n")

            # body
            try:
                reader = vcfpy.Reader.from_path(merged_vcf)
            except Exception as exc:
                handle_critical_error(f"Failed to open merged VCF for finalization: {exc}", exc_cls=MergeConflictError)

            try:
                for record in reader:
                    _recompute_ac_an_af(record)
                    if not _record_passes_filters(
                        record,
                        qual_threshold=qual_threshold,
                        an_threshold=an_threshold,
                        allowed_filter_values=allowed_filter_values,
                    ):
                        continue
                    # anonymize: drop FORMAT and calls
                    record.FORMAT = []
                    record.calls = []
                    outfh.write(_serialize_record_8col(record) + "\n")
            finally:
                try:
                    reader.close()
                except Exception:
                    pass

        shutil.move(temp_out, final_plain_vcf)
        temp_out = None
    except Exception as exc:
        if temp_out and os.path.exists(temp_out):
            try:
                os.remove(temp_out)
            except Exception:
                pass
        if final_plain_vcf and os.path.exists(final_plain_vcf):
            try:
                os.remove(final_plain_vcf)
            except Exception:
                pass
        handle_critical_error(f"Failed to assemble final anonymized VCF contents: {exc}", exc_cls=MergeConflictError)

    # ----- Optional compression + index -----
    if expects_gzip and final_plain_vcf and os.path.exists(final_plain_vcf):
        try:
            pysam.tabix_compress(final_plain_vcf, final_vcf, force=True)
            pysam.tabix_index(final_vcf, preset="vcf", force=True)
        except Exception as exc:
            for p in (final_plain_vcf, final_vcf):
                if p and os.path.exists(p):
                    try:
                        os.remove(p)
                    except Exception:
                        pass
            handle_critical_error(f"Failed to compress or index anonymized VCF: {exc}", exc_cls=MergeConflictError)
    else:
        final_vcf = final_plain_vcf

    log_message(f"Anonymized merged VCF written to: {final_vcf}", verbose)
    return final_vcf

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
    writer_opener = opener

    try:
        input_handle = opener(vcf_path, "rt", encoding="utf-8")
    except OSError as exc:
        handle_critical_error(f"Failed to open {vcf_path}: {exc}", exc_cls=MergeConflictError)
        return

    import tempfile

    temp_path: Optional[str] = None
    try:
        temp_dir = os.path.dirname(os.path.abspath(vcf_path)) or None
        try:
            tmp = tempfile.NamedTemporaryFile(delete=False, dir=temp_dir, suffix=".tmp")
            temp_path = tmp.name
            tmp.close()
        except OSError as exc:
            try:
                input_handle.close()
            except Exception:
                pass
            handle_critical_error(
                f"Failed to create temporary file for {vcf_path}: {exc}",
                exc_cls=MergeConflictError,
            )
            return

        column_header_found = False
        try:
            with input_handle, writer_opener(temp_path, "wt", encoding="utf-8") as out:
                for raw_line in input_handle:
                    line = raw_line.rstrip("\n")
                    if not line:
                        continue
                    if line.startswith("##"):
                        out.write(line + "\n")
                        continue
                    if line.startswith("#CHROM"):
                        column_header_found = True
                        out.write(line + "\n")
                        continue

                    fields = line.split("\t")
                    if len(fields) < 8:
                        out.write(line + "\n")
                        continue

                    # ALT allele list
                    alt_alleles = [a for a in fields[4].split(",") if a]
                    alt_counts = [0] * len(alt_alleles)
                    allele_number = 0

                    # Extract GT from FORMAT and per-sample fields
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
                                for token in genotype.replace("|", "/").split("/"):
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
                                    idx = allele_index - 1
                                    if 0 <= idx < len(alt_counts):
                                        alt_counts[idx] += 1

                    # Parse existing INFO preserving order
                    info_raw = fields[7] if len(fields) >= 8 else ""
                    info_map: dict = {}
                    order: List[str] = []
                    if info_raw and info_raw != ".":
                        for entry in info_raw.split(";"):
                            if not entry:
                                continue
                            if "=" in entry:
                                k, v = entry.split("=", 1)
                            else:
                                k, v = entry, ""
                            if k not in info_map:
                                order.append(k)
                            info_map[k] = v

                    # Update AN/AC/AF
                    info_map["AN"] = str(allele_number)
                    if alt_counts:
                        info_map["AC"] = ",".join(str(c) for c in alt_counts)
                        if allele_number > 0:
                            info_map["AF"] = ",".join(f"{c / allele_number:.6g}" for c in alt_counts)
                        else:
                            info_map["AF"] = ",".join("0" for _ in alt_counts)
                    else:
                        info_map["AC"] = "."
                        info_map["AF"] = "."

                    for k in ("AC", "AN", "AF"):
                        if k not in order:
                            order.append(k)

                    new_entries: List[str] = []
                    for k in order:
                        v = info_map.get(k, "")
                        new_entries.append(f"{k}={v}" if v != "" else f"{k}=")

                    fields[7] = ";".join(new_entries) if new_entries else "."
                    out.write("\t".join(fields) + "\n")
        except OSError as exc:
            if temp_path and os.path.exists(temp_path):
                try:
                    os.remove(temp_path)
                except Exception:
                    pass
            handle_critical_error(f"Failed to rewrite {vcf_path}: {exc}", exc_cls=MergeConflictError)
            return

        if not column_header_found:
            if temp_path and os.path.exists(temp_path):
                try:
                    os.remove(temp_path)
                except Exception:
                    pass
            handle_critical_error(
                f"Failed to parse {vcf_path}: Missing #CHROM header line.",
                exc_cls=MergeConflictError,
            )
            return

        try:
            shutil.move(temp_path, vcf_path)
        except OSError as exc:
            if temp_path and os.path.exists(temp_path):
                try:
                    os.remove(temp_path)
                except Exception:
                    pass
            handle_critical_error(f"Failed to rewrite {vcf_path}: {exc}", exc_cls=MergeConflictError)
            return
        finally:
            temp_path = None
    finally:
        if temp_path and os.path.exists(temp_path):
            try:
                os.remove(temp_path)
            except Exception:
                pass


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
