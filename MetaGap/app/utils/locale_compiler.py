"""Utilities for compiling gettext catalogs without external dependencies."""
from __future__ import annotations

import array
import ast
import codecs
import logging
import struct
from email.parser import HeaderParser
from pathlib import Path
from typing import Dict, Optional

_logger = logging.getLogger(__name__)

_ID_SECTION = 1
_STR_SECTION = 2
_CTX_SECTION = 3


class CatalogCompilationError(RuntimeError):
    """Raised when a ``.po`` file cannot be compiled into a ``.mo`` catalog."""


def _add_entry(
    messages: Dict[bytes, bytes],
    context: Optional[bytes],
    msgid: bytes,
    msgstr: bytes,
    fuzzy: bool,
) -> None:
    if fuzzy:
        return
    key = msgid if context is None else b"%b\x04%b" % (context, msgid)
    messages[key] = msgstr


def _generate_mo(messages: Dict[bytes, bytes]) -> bytes:
    keys = sorted(messages.keys())
    offsets = []
    ids = b""
    strings = b""

    for msgid in keys:
        msgstr = messages[msgid]
        offsets.append((len(ids), len(msgid), len(strings), len(msgstr)))
        ids += msgid + b"\0"
        strings += msgstr + b"\0"

    header_size = 7 * 4
    index_size = 16 * len(keys)
    keystart = header_size + index_size
    valuestart = keystart + len(ids)

    key_offsets = []
    value_offsets = []
    for id_offset, id_length, str_offset, str_length in offsets:
        key_offsets.extend([id_length, id_offset + keystart])
        value_offsets.extend([str_length, str_offset + valuestart])

    all_offsets = key_offsets + value_offsets

    header = struct.pack(
        "Iiiiiii",
        0x950412DE,
        0,
        len(keys),
        header_size,
        header_size + len(keys) * 8,
        0,
        0,
    )

    return header + array.array("i", all_offsets).tobytes() + ids + strings


def _parse_po(po_path: Path) -> Dict[bytes, bytes]:
    try:
        lines = po_path.read_bytes().splitlines(True)
    except OSError as exc:
        raise CatalogCompilationError(str(exc)) from exc

    if lines and lines[0].startswith(codecs.BOM_UTF8):
        raise CatalogCompilationError(
            "Po file starts with a UTF-8 BOM; save without BOM and retry."
        )

    encoding = "latin-1"
    messages: Dict[bytes, bytes] = {}

    section = None
    context: Optional[bytes] = None
    msgid = b""
    msgstr = b""
    fuzzy = False
    is_plural = False

    header_parser = HeaderParser()

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.decode(encoding)

        if line.startswith("#") and section == _STR_SECTION:
            _add_entry(messages, context, msgid, msgstr, fuzzy)
            section = None
            context = None
            fuzzy = False

        if line.startswith("#,") and "fuzzy" in line:
            fuzzy = True

        if line.startswith("#"):
            continue

        if line.startswith("msgctxt"):
            if section == _STR_SECTION:
                _add_entry(messages, context, msgid, msgstr, fuzzy)
            section = _CTX_SECTION
            context = b""
            msgid = b""
            msgstr = b""
            fuzzy = False
            is_plural = False
            content = line[7:]
        elif line.startswith("msgid") and not line.startswith("msgid_plural"):
            if section == _STR_SECTION:
                if not msgid:
                    parsed_header = header_parser.parsestr(msgstr.decode(encoding))
                    charset = parsed_header.get_content_charset()
                    if charset:
                        encoding = charset
                _add_entry(messages, context, msgid, msgstr, fuzzy)
                context = None
            section = _ID_SECTION
            msgid = b""
            msgstr = b""
            fuzzy = False
            is_plural = False
            content = line[5:]
        elif line.startswith("msgid_plural"):
            if section != _ID_SECTION:
                raise CatalogCompilationError(
                    f"msgid_plural not preceded by msgid on line {line_number}"
                )
            section = _ID_SECTION
            content = line[12:]
            msgid += b"\0"
            is_plural = True
        elif line.startswith("msgstr"):
            section = _STR_SECTION
            if line.startswith("msgstr["):
                if not is_plural:
                    raise CatalogCompilationError(
                        f"plural without msgid_plural on line {line_number}"
                    )
                content = line.split("]", 1)[1]
                if msgstr:
                    msgstr += b"\0"
            else:
                if is_plural:
                    raise CatalogCompilationError(
                        f"indexed msgstr required for plural on line {line_number}"
                    )
                content = line[6:]
        else:
            content = line

        content = content.strip()
        if not content:
            continue

        try:
            parsed = ast.literal_eval(content)
        except SyntaxError as exc:
            raise CatalogCompilationError(
                f"Syntax error on line {line_number}: {content}"
            ) from exc

        if section == _CTX_SECTION:
            context = (context or b"") + parsed.encode(encoding)
        elif section == _ID_SECTION:
            msgid += parsed.encode(encoding)
        elif section == _STR_SECTION:
            msgstr += parsed.encode(encoding)
        else:
            raise CatalogCompilationError(
                f"Syntax error before line {line_number}: {content}"
            )

    if section == _STR_SECTION:
        _add_entry(messages, context, msgid, msgstr, fuzzy)

    return messages


def compile_po_to_mo(po_path: Path, mo_path: Path) -> None:
    messages = _parse_po(po_path)
    output = _generate_mo(messages)
    mo_path.parent.mkdir(parents=True, exist_ok=True)
    mo_path.write_bytes(output)


def ensure_compiled_catalogs(locale_root: Path) -> None:
    for po_path in locale_root.rglob("*.po"):
        mo_path = po_path.with_suffix(".mo")
        try:
            if mo_path.exists() and mo_path.stat().st_mtime >= po_path.stat().st_mtime:
                continue
            compile_po_to_mo(po_path, mo_path)
        except CatalogCompilationError as exc:
            _logger.warning("Failed to compile %s: %s", po_path, exc)
        except OSError as exc:
            _logger.warning("Unable to write compiled catalog %s: %s", mo_path, exc)

