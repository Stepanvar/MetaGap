"""Metadata helpers for augmenting merged VCF files."""

from __future__ import annotations

import copy
import gzip
import os
import re
import subprocess
from collections import OrderedDict
from typing import List, Optional, Tuple

from . import VCFPY_AVAILABLE, vcfpy
from .logging_utils import handle_critical_error, log_message

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
                f"WARNING: Unable to create INFO header definition for {info_id}: {exc}",
                verbose,
            )
            continue

        try:
            header.add_line(info_line)
        except Exception as exc:  # pragma: no cover - defensive logging
            log_message(
                f"WARNING: Failed to append INFO header definition for {info_id}: {exc}",
                verbose,
            )
            continue

        existing_ids.add(info_id)
        added_ids.append(info_id)

    if added_ids:
        log_message(
            "Inserted INFO header definitions for missing fields: "
            + ", ".join(added_ids),
            verbose,
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

    normalized_path = os.path.abspath(template_path)
    if not os.path.exists(normalized_path):
        handle_critical_error(
            f"Metadata template header file does not exist: {template_path}"
        )
    if not os.path.isfile(normalized_path):
        handle_critical_error(
            f"Metadata template header path is not a file: {template_path}"
        )

    template_sample_mapping: Optional["OrderedDict[str, str]"] = None
    template_serialized_sample: Optional[str] = None
    sanitized_lines: List[str] = []
    simple_header_lines = []
    existing_simple = set()
    existing_sanitized = set()

    try:
        with open(normalized_path, "r", encoding="utf-8") as handle:
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
                    )
                    continue
                if stripped.startswith("##SAMPLE="):
                    try:
                        parsed_sample = _parse_sample_metadata_line(stripped)
                    except ValueError as exc:
                        handle_critical_error(
                            "Invalid SAMPLE metadata line in template "
                            f"{normalized_path} line {line_number}: {exc}"
                        )
                    if template_sample_mapping is not None:
                        log_message(
                            (
                                "WARNING: Multiple ##SAMPLE lines found in metadata "
                                f"template {normalized_path}; using the last occurrence"
                            ),
                            verbose,
                        )
                    template_sample_mapping = OrderedDict(parsed_sample.items())
                    template_serialized_sample = stripped
                    continue

                parsed = _parse_simple_metadata_line(stripped)
                if not parsed:
                    handle_critical_error(
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
                        f"WARNING: Failed to parse metadata template line '{sanitized}': {exc}",
                        verbose,
                    )
                    continue

                identifier = (simple_line.key, getattr(simple_line, "value", None))
                if identifier in existing_simple:
                    continue
                simple_header_lines.append(simple_line)
                existing_simple.add(identifier)
    except OSError as exc:
        handle_critical_error(
            f"Unable to read metadata template header file {template_path}: {exc}"
        )

    log_message(
        f"Loaded metadata template header from {normalized_path}",
        verbose,
    )
    if sanitized_lines:
        log_message(
            "Metadata template header lines: " + ", ".join(sanitized_lines),
            verbose,
        )
    if template_serialized_sample:
        log_message(
            f"Metadata template sample line: {template_serialized_sample}",
            verbose,
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
        handle_critical_error("A metadata file path must be provided.")

    normalized_path = os.path.abspath(metadata_file)
    if not os.path.exists(normalized_path):
        handle_critical_error(f"Metadata file does not exist: {metadata_file}")
    if not os.path.isfile(normalized_path):
        handle_critical_error(f"Metadata file path is not a file: {metadata_file}")

    sanitized_lines: List[str] = []
    seen_lines = set()

    try:
        with open(normalized_path, "r", encoding="utf-8") as handle:
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
                        f"Problematic entry at {normalized_path} line {line_number}: {raw_line.strip()}"
                    )
                key, value = parsed
                sanitized = f"##{key}={value}"
                if sanitized in seen_lines:
                    continue
                sanitized_lines.append(sanitized)
                seen_lines.add(sanitized)
    except OSError as exc:
        handle_critical_error(f"Unable to read metadata file {metadata_file}: {exc}")

    log_message(
        f"Loaded {len(sanitized_lines)} metadata header line(s) from {normalized_path}",
        verbose,
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

    log_message("Applied CLI metadata to merged header.", verbose)
    return header


def _validate_anonymized_vcf_header(final_vcf_path: str, ensure_for_uncompressed: bool = False):
    """Use ``bcftools`` to confirm that ``final_vcf_path`` has a readable header."""

    if not final_vcf_path:
        return

    requires_validation = ensure_for_uncompressed or final_vcf_path.endswith(".vcf.gz")
    if not requires_validation:
        return

    try:
        subprocess.run(
            ["bcftools", "view", "-h", final_vcf_path],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except (subprocess.CalledProcessError, OSError) as exc:
        handle_critical_error(
            f"Failed to read header from anonymized VCF using bcftools: {exc}"
        )


def append_metadata_to_merged_vcf(
    merged_vcf: str,
    sample_metadata_entries=None,
    header_metadata_lines=None,
    serialized_sample_line=None,
    verbose: bool = False,
):
    """Finalize the merged VCF by applying bcftools processing and metadata."""

    header_metadata_lines = header_metadata_lines or []

    joint_temp = f"{merged_vcf}.joint.temp.vcf"
    filtered_temp = f"{merged_vcf}.filtered.temp.vcf"

    expects_gzip = merged_vcf.endswith(".gz")
    final_plain_vcf = None
    final_vcf = None

    if merged_vcf.endswith(".vcf.gz"):
        base_path = merged_vcf[: -len(".vcf.gz")]
        final_plain_vcf = f"{base_path}.anonymized.vcf"
        final_vcf = f"{final_plain_vcf}.gz"
    elif merged_vcf.endswith(".vcf"):
        base_path = merged_vcf[: -len(".vcf")]
        final_plain_vcf = f"{base_path}.anonymized.vcf"
        final_vcf = final_plain_vcf
    else:
        base_path, ext = os.path.splitext(merged_vcf)
        if ext == ".gz":
            base_path = base_path or merged_vcf[:-3]
            final_plain_vcf = f"{base_path}.anonymized"
            final_vcf = f"{final_plain_vcf}.gz"
        else:
            final_plain_vcf = f"{base_path}.anonymized{ext or '.vcf'}"
            final_vcf = final_plain_vcf

    def _cleanup_temp_files():
        for path in [joint_temp, filtered_temp]:
            try:
                if path and os.path.exists(path):
                    os.remove(path)
            except Exception:
                pass

    log_message(
        "Recalculating AC, AN, AF tags across the merged cohort...",
        verbose,
    )
    try:
        subprocess.run(
            [
                "bcftools",
                "+fill-tags",
                merged_vcf,
                "-O",
                "v",
                "-o",
                joint_temp,
                "--",
                "-t",
                "AC,AN,AF",
            ],
            check=True,
        )
    except (subprocess.CalledProcessError, OSError) as exc:
        _cleanup_temp_files()
        handle_critical_error(
            f"Failed to recalculate INFO tags with bcftools +fill-tags: {exc}"
        )

    log_message("Applying quality filters to remove low-confidence variants...", verbose)
    try:
        subprocess.run(
            [
                "bcftools",
                "view",
                joint_temp,
                "-i",
                'QUAL>30 && INFO/AN>50 && FILTER="PASS"',
                "-O",
                "v",
                "-o",
                filtered_temp,
            ],
            check=True,
        )
    except (subprocess.CalledProcessError, OSError) as exc:
        _cleanup_temp_files()
        handle_critical_error(
            f"Failed to apply bcftools filtering to merged VCF: {exc}"
        )

    log_message("Removing individual sample genotype columns (anonymizing data)...", verbose)
    try:
        reader = vcfpy.Reader.from_path(filtered_temp)
    except Exception as exc:
        _cleanup_temp_files()
        handle_critical_error(
            f"Failed to open filtered VCF for anonymization: {exc}"
        )

    try:
        anonymized_header = reader.header.copy()

        # Remove FORMAT definitions and sample columns from the header
        filtered_lines = []
        for line in getattr(anonymized_header, "lines", []):
            if isinstance(line, vcfpy.header.FormatHeaderLine):
                continue
            if getattr(line, "key", None) == "FORMAT":
                continue
            filtered_lines.append(line)
        anonymized_header.lines = filtered_lines

        if hasattr(anonymized_header, "formats"):
            try:
                anonymized_header.formats.clear()
            except AttributeError:
                anonymized_header.formats = OrderedDict()

        if hasattr(anonymized_header, "samples") and hasattr(
            anonymized_header.samples, "names"
        ):
            anonymized_header.samples.names = []

        # Attach sample metadata if supplied
        if sample_metadata_entries or serialized_sample_line is not None:
            try:
                mapping = (
                    OrderedDict(sample_metadata_entries)
                    if sample_metadata_entries
                    else _parse_sample_metadata_line(serialized_sample_line)
                )
                sample_line = vcfpy.header.SampleHeaderLine.from_mapping(mapping)
            except Exception as exc:
                _cleanup_temp_files()
                reader.close()
                handle_critical_error(str(exc))

            anonymized_header.lines = [
                line
                for line in anonymized_header.lines
                if getattr(line, "key", None) != "SAMPLE"
            ]
            anonymized_header.add_line(sample_line)

        # Add simple metadata header lines supplied via CLI/template
        existing_simple = {
            (line.key, getattr(line, "value", None))
            for line in getattr(anonymized_header, "lines", [])
            if isinstance(line, vcfpy.SimpleHeaderLine)
        }

        for metadata_line in header_metadata_lines:
            normalized = metadata_line.strip()
            if not normalized:
                continue
            if not normalized.startswith("##"):
                normalized = "##" + normalized
            parsed = _parse_simple_metadata_line(normalized)
            if not parsed:
                continue
            key, value = parsed
            identifier = (key, value)
            if identifier in existing_simple:
                continue
            try:
                anonymized_header.add_line(vcfpy.SimpleHeaderLine(key, value))
            except Exception as exc:  # pragma: no cover - defensive logging
                log_message(
                    f"WARNING: Failed to add metadata header line '{normalized}': {exc}",
                    verbose,
                )
                continue
            existing_simple.add(identifier)

        anonymized_header = ensure_standard_info_definitions(
            anonymized_header, verbose=verbose
        )

        if not final_plain_vcf:
            base, ext = os.path.splitext(merged_vcf)
            if not ext:
                ext = ".vcf"
            final_plain_vcf = f"{base}.anonymized{ext}"

        try:
            writer = vcfpy.Writer.from_path(final_plain_vcf, anonymized_header)
        except Exception as exc:
            reader.close()
            _cleanup_temp_files()
            if final_plain_vcf and os.path.exists(final_plain_vcf):
                try:
                    os.remove(final_plain_vcf)
                except Exception:
                    pass
            handle_critical_error(
                f"Failed to open anonymized VCF writer: {exc}"
            )

        with writer:
            for record in reader:
                anonymized_record = vcfpy.Record(
                    CHROM=record.CHROM,
                    POS=record.POS,
                    ID=record.ID,
                    REF=record.REF,
                    ALT=copy.deepcopy(record.ALT),
                    QUAL=record.QUAL,
                    FILTER=copy.deepcopy(record.FILTER),
                    INFO=copy.deepcopy(record.INFO),
                    FORMAT=None,
                    calls=[],
                )
                writer.write_record(anonymized_record)

        final_vcf = final_plain_vcf
    except Exception as exc:
        _cleanup_temp_files()
        if final_plain_vcf and os.path.exists(final_plain_vcf):
            try:
                os.remove(final_plain_vcf)
            except Exception:
                pass
        handle_critical_error(
            f"Failed to assemble final anonymized VCF contents: {exc}"
        )
    finally:
        try:
            reader.close()
        except Exception:
            pass

    if expects_gzip and final_plain_vcf and os.path.exists(final_plain_vcf):
        if pysam is None:
            _cleanup_temp_files()
            handle_critical_error(
                "pysam is required to BGZF-compress and index the anonymized VCF."
            )

        try:
            final_vcf = f"{final_plain_vcf}.gz"
            pysam.bgzip(final_plain_vcf, final_vcf, force=True)
            pysam.tabix_index(final_vcf, preset="vcf", force=True)
        except (OSError, ValueError) as exc:
            _cleanup_temp_files()
            if os.path.exists(final_plain_vcf):
                try:
                    os.remove(final_plain_vcf)
                except Exception:
                    pass
            if os.path.exists(final_vcf):
                try:
                    os.remove(final_vcf)
                except Exception:
                    pass
            handle_critical_error(
                f"Failed to compress or index anonymized VCF: {exc}"
            )
        finally:
            if os.path.exists(final_plain_vcf):
                try:
                    os.remove(final_plain_vcf)
                except Exception:
                    pass
    elif not expects_gzip:
        final_vcf = final_plain_vcf

    _cleanup_temp_files()

    _validate_anonymized_vcf_header(final_vcf, ensure_for_uncompressed=True)

    log_message(f"Anonymized merged VCF written to: {final_vcf}", verbose)
    return final_vcf


def recalculate_cohort_info_tags(vcf_path: str, verbose: bool = False) -> None:
    """Recompute AC, AN, and AF INFO tags in-place for ``vcf_path``."""

    log_message(
        "Recalculating AC, AN, AF tags across the merged cohort...",
        verbose,
    )

    if VCFPY_AVAILABLE:
        try:
            reader = vcfpy.Reader.from_path(vcf_path)
        except Exception as exc:
            handle_critical_error(
                f"Failed to open {vcf_path} for AC/AN/AF recalculation: {exc}"
            )

        records = []
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
                    info["AF"] = [count / allele_number for count in alt_counts]
                else:
                    info["AF"] = [0.0 for _ in alt_counts]
            else:
                info["AC"] = []
                info["AF"] = []

            record.INFO = info
            records.append(record)

        try:
            reader.close()
        except Exception:
            pass

        try:
            writer = vcfpy.Writer.from_path(vcf_path, header)
        except Exception as exc:
            handle_critical_error(f"Failed to reopen {vcf_path} for writing: {exc}")

        try:
            for record in records:
                writer.write_record(record)
        finally:
            writer.close()
        return

    _recalculate_cohort_info_tags_without_vcfpy(vcf_path)


def _recalculate_cohort_info_tags_without_vcfpy(vcf_path: str) -> None:
    opener = gzip.open if vcf_path.endswith(".gz") else open
    try:
        with opener(vcf_path, "rt", encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle]
    except OSError as exc:
        handle_critical_error(f"Failed to open {vcf_path}: {exc}")

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
            f"Failed to parse {vcf_path}: Missing #CHROM header line."
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
        handle_critical_error(f"Failed to rewrite {vcf_path}: {exc}")


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
