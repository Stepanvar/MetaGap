"""File helpers used by the VCF importer."""

from __future__ import annotations

from typing import Any, Dict, Iterable, Optional

from .vcf_database import VCFDatabaseWriter
from .vcf_metadata import normalize_metadata_value


def extract_metadata_text_fallback(file_path: str) -> Dict[str, Any]:
    """Extract metadata from a VCF file via text parsing."""

    metadata: Dict[str, Any] = {}
    with open(file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("##SAMPLE"):
                start = stripped.find("<")
                end = stripped.rfind(">")
                if start == -1 or end == -1 or end <= start:
                    continue
                content = stripped[start + 1 : end]
                for item in split_sample_attributes(content):
                    if "=" not in item:
                        continue
                    key, value = item.split("=", 1)
                    metadata[key.lower()] = normalize_metadata_value(value)
                metadata.setdefault("name", metadata.get("id"))
            if stripped.startswith("#CHROM"):
                break
    return metadata


def parse_vcf_text_fallback(
    file_path: str,
    sample_group: Any,
    writer: VCFDatabaseWriter,
) -> None:
    """Parse a VCF via text processing when pysam is unavailable."""

    header_sample: Optional[str] = None

    with open(file_path, "r", encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#CHROM"):
                columns = stripped.lstrip("#").split("\t")
                if len(columns) > 9:
                    header_sample = columns[9]
                continue
            if stripped.startswith("#"):
                continue

            fields = stripped.split("\t")
            if len(fields) < 8:
                continue

            chrom, pos, variant_id, ref, alt, qual, filter_value, info_field = fields[:8]
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

            writer.create_allele_frequency(
                sample_group,
                chrom=chrom,
                pos=int(pos),
                variant_id=None if variant_id == "." else variant_id,
                ref=ref,
                alt=writer.serialize_alt((alt,)),
                qual=None if qual in {".", ""} else float(qual),
                filter_value=writer.serialize_filter(filter_value),
                info=info_instance,
                format_instance=format_instance,
                format_sample=sample_identifier,
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


__all__ = [
    "extract_metadata_text_fallback",
    "parse_vcf_text_fallback",
    "split_sample_attributes",
]
