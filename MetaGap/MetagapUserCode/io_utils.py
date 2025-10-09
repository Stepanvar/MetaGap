"""Utility helpers shared across VCF processing modules."""

from __future__ import annotations

from typing import Optional, Tuple

import re


def parse_simple_metadata_line(line: str) -> Optional[Tuple[str, str]]:
    """Parse a VCF metadata line of the form ``##key=value``.

    Parameters
    ----------
    line:
        The raw metadata text line to parse.

    Returns
    -------
    Optional[Tuple[str, str]]
        ``None`` when the line is not a simple ``##key=value`` metadata entry.
        Otherwise a ``(key, value)`` tuple with leading/trailing whitespace stripped.
    """

    stripped = line.strip()
    if not stripped.startswith("##") or "=" not in stripped:
        return None

    key, value = stripped[2:].split("=", 1)
    key = key.strip()
    value = value.strip()
    if not key:
        return None
    return key, value


def preprocess_vcf(file_path: str) -> str:
    """Normalise spacing issues in a VCF file and return the usable path.

    Some VCF producers emit files that use spaces instead of tabs for the
    ``#CHROM`` header row and subsequent data rows.  ``vcfpy`` expects
    tab-delimited columns, so we rewrite these rows when needed.  Metadata
    lines (``##``) are preserved verbatim.

    Parameters
    ----------
    file_path:
        The path to the original VCF file on disk.

    Returns
    -------
    str
        The path to a VCF file safe for consumption by ``vcfpy``.  If no
        modifications were required the original path is returned; otherwise a
        temporary sibling file with a ``.tmp`` suffix is created and returned.
    """

    with open(file_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    modified = False
    new_lines = []
    header_found = False

    for line in lines:
        if line.startswith("##"):
            new_lines.append(line)
            continue

        if line.startswith("#"):
            new_line = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
            new_lines.append(new_line)
            header_found = True
            if new_line != line:
                modified = True
            continue

        if header_found:
            new_line = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
            new_lines.append(new_line)
            if new_line != line:
                modified = True
        else:
            new_lines.append(line)

    if not modified:
        return file_path

    temp_file = f"{file_path}.tmp"
    with open(temp_file, "w", encoding="utf-8") as handle:
        handle.writelines(new_lines)
    return temp_file


__all__ = ["parse_simple_metadata_line", "preprocess_vcf"]
