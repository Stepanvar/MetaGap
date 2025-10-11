"""File helpers used by the VCF importer."""

from __future__ import annotations

import gzip
import logging
import re
from contextlib import contextmanager
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Optional, TextIO

from django.core.exceptions import ValidationError

from .vcf_database import VCFDatabaseWriter
from .vcf_metadata import (
    METADATA_SECTION_MAP,
    VCFMetadataParser,
    normalize_metadata_key,
    normalize_metadata_value,
)

logger = logging.getLogger(__name__)


@contextmanager
def open_vcf_text(file_path: str) -> Iterator[TextIO]:
    """Open a VCF file as text, transparently handling gzip compression."""

    path = Path(file_path)
    suffixes = {suffix.lower() for suffix in path.suffixes}
    if {".gz", ".bgz"} & suffixes:
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8") as handle:
            yield handle


def extract_metadata_text_fallback(
    file_path: str, warnings: Optional[list[str]] = None
) -> Dict[str, Any]:
    """Extract metadata from a VCF file via text parsing."""

    metadata: Dict[str, Any] = {}
    parser = VCFMetadataParser(warnings)
    try:
        with open_vcf_text(file_path) as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#CHROM"):
                    break
                if not stripped.startswith("##"):
                    continue

                for entry_key, entry_items in _iter_metadata_entries(stripped):
                    parser.ingest_metadata_items(metadata, entry_key, entry_items)

                    normalized_entry_key = (entry_key or "").upper()
                    skip_raw_capture = (
                        normalized_entry_key in METADATA_SECTION_MAP
                        or normalized_entry_key == "SEQUENCING_PLATFORM"
                    )

                    if skip_raw_capture:
                        continue

                    for raw_key, value in entry_items.items():
                        if raw_key.lower() == "value":
                            continue
                        normalized_raw = normalize_metadata_key(raw_key)
                        if normalized_raw and normalized_raw not in metadata:
                            metadata[normalized_raw] = value
        if "name" not in metadata and "sample_group_name" in metadata:
            metadata["name"] = metadata["sample_group_name"]
    except UnicodeDecodeError as exc:  # pragma: no cover - defensive fallback
        raise ValidationError(
            "The uploaded VCF file could not be decoded using UTF-8. "
            "Please ensure the file is UTF-8 encoded or convert it before retrying."
        ) from exc
    return metadata


def parse_vcf_text_fallback(
    file_path: str,
    sample_group: Any,
    writer: VCFDatabaseWriter,
    warnings: Optional[list[str]] = None,
) -> None:
    """Parse a VCF via text processing when pysam is unavailable."""

    header_sample: Optional[str] = None
    header_seen = False
    warnings_list = warnings if warnings is not None else None

    try:
        with open_vcf_text(file_path) as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue

                # Header parsing
                if stripped.startswith("#CHROM"):
                    header_seen = True
                    columns = stripped.lstrip("#").split("\t")
                    if len(columns) > 9:
                        header_sample = columns[9]
                    continue
                if stripped.startswith("#"):
                    continue

                # Data line
                fields = stripped.split("\t")
                if len(fields) < 8:
                    msg = (
                        "Encountered a truncated data line at line "
                        f"{line_number}: '{stripped}'."
                    )
                    if warnings_list is not None:
                        warnings_list.append(msg)
                    logger.warning("%s", msg)
                    continue

                chrom, pos, variant_id, ref, alt, qual, filter_value, info_field = fields[:8]

                # INFO map
                info_mapping: Dict[str, Any] = {}
                for entry in info_field.split(";"):
                    if not entry:
                        continue
                    if "=" in entry:
                        key, value = entry.split("=", 1)
                        info_mapping[key] = value
                    else:
                        info_mapping[entry] = True

                info_instance = writer.create_info_instance(info_mapping)

                # FORMAT / sample
                format_instance: Optional[Any] = None
                sample_identifier: Optional[str] = None
                if len(fields) > 9:
                    format_keys = fields[8].split(":") if len(fields) > 8 else []
                    format_values = fields[9].split(":")
                    sample_identifier = header_sample or "Sample"
                    samples = (
                        {sample_identifier: dict(zip(format_keys, format_values))}
                        if format_keys
                        else {}
                    )
                    if samples:
                        format_instance, sample_identifier = writer.create_format_instance(samples)

                # POS
                try:
                    pos_value = int(pos)
                except ValueError as exc:
                    raise ValidationError(
                        f"Invalid position value '{pos}' encountered on line "
                        f"{line_number}: '{stripped}'."
                    ) from exc

                # QUAL
                if qual in {".", ""}:
                    qual_value: Optional[float] = None
                else:
                    try:
                        qual_value = float(qual)
                    except ValueError as exc:
                        raise ValidationError(
                            f"Invalid quality value '{qual}' encountered on line "
                            f"{line_number}: '{stripped}'."
                        ) from exc

                # Persist
                writer.create_allele_frequency(
                    sample_group,
                    chrom=chrom,
                    pos=pos_value,
                    variant_id=None if variant_id == "." else variant_id,
                    ref=ref,
                    alt=writer.serialize_alt((alt,)),
                    qual=qual_value,
                    filter_value=writer.serialize_filter(filter_value),
                    info=info_instance,
                    format_instance=format_instance,
                    format_sample=sample_identifier,
                )
    except UnicodeDecodeError as exc:  # pragma: no cover - defensive fallback
        raise ValidationError(
            "The uploaded VCF file could not be decoded using UTF-8. "
            "Please ensure the file is UTF-8 encoded or convert it before retrying."
        ) from exc

    if not header_seen:
        raise ValidationError(
            "The uploaded VCF file is missing a '#CHROM' header line and cannot be parsed."
        )

def split_sample_attributes(content: str) -> Iterable[str]:
    items: list[str] = []
    current: list[str] = []
    quote_char: Optional[str] = None
    escape = False
    bracket_stack: list[str] = []
    opening = {"{": "}", "[": "]", "(": ")"}
    closing = {value: key for key, value in opening.items()}

    for char in content:
        if quote_char:
            current.append(char)
            if escape:
                escape = False
                continue
            if char == "\\":
                escape = True
                continue
            if char == quote_char:
                quote_char = None
            continue

        if char in {'"', "'"}:
            quote_char = char
            current.append(char)
            continue

        if char in opening:
            bracket_stack.append(char)
            current.append(char)
            continue

        if char in closing:
            if bracket_stack and bracket_stack[-1] == closing[char]:
                bracket_stack.pop()
            current.append(char)
            continue

        if char == "," and not bracket_stack:
            item = "".join(current).strip()
            if item:
                items.append(item)
            current = []
            continue

        current.append(char)

    tail = "".join(current).strip()
    if tail:
        items.append(tail)

    return items


def _iter_metadata_entries(line: str) -> Iterable[Tuple[str, Dict[str, Any]]]:
    if not line.startswith("##"):
        return []

    normalized = "##" + re.sub(r"(?<=.)##", "\n##", line[2:])
    entries: list[Tuple[str, Dict[str, Any]]] = []
    for candidate in normalized.splitlines():
        if not candidate.startswith("##"):
            continue
        key, items = _parse_metadata_entry(candidate[2:])
        if key:
            entries.append((key, items))
    return entries


def _parse_metadata_entry(entry: str) -> Tuple[str, Dict[str, Any]]:
    key, _, remainder = entry.partition("=")
    key = key.strip()
    remainder = remainder.strip()

    if remainder.startswith("<") and remainder.endswith(">"):
        items = _parse_angle_bracket_items(remainder[1:-1])
    elif remainder:
        items = {"value": normalize_metadata_value(remainder)}
    else:
        items = {}
    return key, items


def _parse_angle_bracket_items(content: str) -> Dict[str, Any]:
    parsed: Dict[str, Any] = {}
    for item in split_sample_attributes(content):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        parsed[key] = normalize_metadata_value(value)
    return parsed


__all__ = [
    "extract_metadata_text_fallback",
    "parse_vcf_text_fallback",
    "split_sample_attributes",
    "open_vcf_text",
]
