"""Utilities for handling VCF header metadata configuration."""
from __future__ import annotations

from collections import OrderedDict
from typing import Callable, Optional, Sequence, Tuple

try:
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


_LogMessageFn = Callable[[str, bool], None]
_HandleErrorFn = Callable[[str], None]


def _safe_log(log_fn: Optional[_LogMessageFn], message: str, verbose: bool) -> None:
    if log_fn is None:
        return
    log_fn(message, verbose)


def _handle_error(
    handler: Optional[_HandleErrorFn], message: str
) -> None:  # pragma: no cover - control flow depends on caller
    if handler is not None:
        handler(message)
    raise ValueError(message)


def _format_sample_metadata_value(value: str) -> str:
    """Return a VCF-safe representation of the provided metadata value."""

    value = str(value)
    needs_quotes = any(char in value for char in [" ", ",", "\t", '"', "<", ">", "="])
    if not needs_quotes:
        return value

    escaped = value.replace("\\", "\\\\").replace('"', '\\"')
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
    entries: "OrderedDict[str, str]" = OrderedDict()

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


def parse_metadata_arguments(
    args,
    verbose: bool = False,
    *,
    log_message: Optional[_LogMessageFn] = None,
    handle_critical_error: Optional[_HandleErrorFn] = None,
):
    """Return header metadata derived from CLI arguments."""

    sample_entries = getattr(args, "sample_metadata_entries", None) or []
    additional_lines = getattr(args, "header_metadata_lines", None) or []

    sample_mapping: "OrderedDict[str, str]" = OrderedDict()
    for raw_entry in sample_entries:
        entry = raw_entry.strip()
        if not entry:
            continue
        if "=" not in entry:
            _handle_error(
                handle_critical_error,
                f"Invalid sample metadata entry '{raw_entry}'. Expected KEY=VALUE format.",
            )
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            _handle_error(handle_critical_error, "Sample metadata keys cannot be empty.")
        sample_mapping[key] = value

    sample_header_line = None
    serialized_sample_line = None
    if sample_mapping:
        try:
            sample_line = build_sample_metadata_line(sample_mapping)
        except ValueError as exc:
            _handle_error(handle_critical_error, str(exc))
        _safe_log(log_message, f"Using CLI sample metadata: {sample_line}", verbose)
        serialized_sample_line = sample_line
        sample_header_line = vcfpy.SampleHeaderLine.from_mapping(sample_mapping)

    simple_header_lines = []
    sanitized_header_lines = []
    for raw_line in additional_lines:
        normalized = raw_line.strip()
        if not normalized:
            continue
        if not normalized.startswith("##"):
            normalized = "##" + normalized
        parsed = _parse_simple_metadata_line(normalized)
        if not parsed:
            _handle_error(
                handle_critical_error,
                f"Additional metadata '{raw_line}' must be in '##key=value' format.",
            )
        key, value = parsed
        simple_header_lines.append(vcfpy.SimpleHeaderLine(key, value))
        sanitized_header_lines.append(f"##{key}={value}")

    if simple_header_lines:
        summary = ", ".join(f"##{line.key}={line.value}" for line in simple_header_lines)
        _safe_log(log_message, f"Using CLI header metadata lines: {summary}", verbose)

    sample_mapping_copy = OrderedDict(sample_mapping) if sample_mapping else None

    return (
        sample_header_line,
        simple_header_lines,
        sample_mapping_copy,
        sanitized_header_lines,
        serialized_sample_line,
    )


def apply_metadata_to_header(
    header,
    sample_header_line=None,
    simple_header_lines: Optional[Sequence[vcfpy.HeaderLine]] = None,
    verbose: bool = False,
    *,
    log_message: Optional[_LogMessageFn] = None,
):
    """Return ``header`` with CLI metadata applied."""

    if sample_header_line is None and not simple_header_lines:
        return header

    simple_header_lines = list(simple_header_lines or [])

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

    _safe_log(log_message, "Applied CLI metadata to merged header.", verbose)
    return header


__all__ = [
    "_format_sample_metadata_value",
    "build_sample_metadata_line",
    "_parse_sample_metadata_line",
    "_parse_simple_metadata_line",
    "parse_metadata_arguments",
    "apply_metadata_to_header",
]
