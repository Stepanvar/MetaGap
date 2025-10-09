#!/usr/bin/env python3
"""
This script consolidates multiple VCF files into one merged VCF file.
It replicates the bash script functionality:
  - Parses command-line options for the input directory, output directory, reference genome, VCF version, and optional metadata.
  - Validates individual VCF files for header fileformat and reference genome.
  - Merges valid VCFs using vcfpy.
  - Appends metadata (if provided via CLI) to the merged VCF header.
  - Performs a final validation of the merged VCF.
  - Logs execution details to a log file.

At completion the script prints a compact summary in the form
``Wrote: <DIR>/<SAMPLE_XXXX>.vcf[.gz] x N`` where the directory, a
representative sample filename, and the number of produced cohort VCF
shards are derived from the output directory contents.

Requirements:
pip install vcfpy
"""

import os
import sys
import glob
import argparse
import csv
import datetime
import re
import copy
import tempfile
import shutil
import subprocess
import gzip
from collections import OrderedDict

from validation import (
    LOG_FILE,
    VCFPY_AVAILABLE,
    handle_critical_error,
    handle_non_critical_error,
    is_gvcf_header,
    log_message,
    normalize_vcf_version,
    preprocess_vcf,
    find_first_vcf_with_header,
    validate_vcf,
    validate_all_vcfs,
    validate_merged_vcf,
    vcfpy,
)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Consolidate multiple VCF files into a single merged VCF file."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the VCF files to merge.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        "--output",
        dest="output_dir",
        help="Directory for the merged VCF. Defaults to the input directory when omitted.",
    )
    parser.add_argument(
        "--ref",
        help="Expected reference genome build. When omitted the script attempts to auto-detect it.",
    )
    parser.add_argument(
        "--vcf-version",
        help="Expected VCF version (e.g., 4.2). When omitted the script attempts to auto-detect it.",
    )
    parser.add_argument(
        "--meta",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Sample metadata in KEY=VALUE form. Repeat for multiple keys. "
            "An ID entry is required when metadata is provided."
        ),
    )
    parser.add_argument(
        "--allow-gvcf",
        action="store_true",
        help="Allow input files that contain gVCF annotations such as <NON_REF> ALT alleles.",
    )
    parser.add_argument(
        "--sample-metadata",
        dest="sample_metadata_entries",
        action="append",
        metavar="KEY=VALUE",
        default=[],
        help="Add a key/value pair to the ##SAMPLE metadata line. Provide at least ID=... to emit the line.",
    )
    parser.add_argument(
        "--header-metadata",
        dest="header_metadata_lines",
        action="append",
        metavar="LINE",
        default=[],
        help="Add an arbitrary metadata header line (##key=value). The '##' prefix is optional.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose console logging in addition to the log file.",
    )
    return parser.parse_args()

def _format_sample_metadata_value(value: str) -> str:
    """Return a VCF-safe representation of the provided metadata value."""

    value = str(value)
    needs_quotes = any(char in value for char in [" ", ",", "\t", "\"", "<", ">", "="])
    if not needs_quotes:
        return value

    escaped = value.replace("\\", "\\\\").replace("\"", "\\\"")
    return f'"{escaped}"'


def build_sample_metadata_line(entries: "OrderedDict[str, str]") -> str:
    """Serialize an ordered mapping into a single ##SAMPLE metadata line."""

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

    token = []
    stack = []
    in_quotes = False
    escape = False

    def flush_token():
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
            unescaped = []
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


def parse_metadata_arguments(args, verbose=False):
    """Return header metadata derived from CLI arguments."""

    sample_entries = getattr(args, "sample_metadata_entries", None) or []
    additional_lines = getattr(args, "header_metadata_lines", None) or []

    sample_mapping = OrderedDict()
    for raw_entry in sample_entries:
        entry = raw_entry.strip()
        if not entry:
            continue
        if "=" not in entry:
            handle_critical_error(
                f"Invalid sample metadata entry '{raw_entry}'. Expected KEY=VALUE format."
            )
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if not key:
            handle_critical_error("Sample metadata keys cannot be empty.")
        sample_mapping[key] = value

    sample_header_line = None
    serialized_sample_line = None
    if sample_mapping:
        try:
            sample_line = build_sample_metadata_line(sample_mapping)
        except ValueError as exc:
            handle_critical_error(str(exc))
        log_message(f"Using CLI sample metadata: {sample_line}", verbose)
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
            handle_critical_error(
                f"Additional metadata '{raw_line}' must be in '##key=value' format."
            )
        key, value = parsed
        simple_header_lines.append(vcfpy.SimpleHeaderLine(key, value))
        sanitized_header_lines.append(f"##{key}={value}")

    if simple_header_lines:
        log_message(
            "Using CLI header metadata lines: "
            + ", ".join(f"##{line.key}={line.value}" for line in simple_header_lines),
            verbose,
        )

    sample_mapping_copy = OrderedDict(sample_mapping) if sample_mapping else None

    return (
        sample_header_line,
        simple_header_lines,
        sample_mapping_copy,
        sanitized_header_lines,
        serialized_sample_line,
    )


def apply_metadata_to_header(
    header, sample_header_line=None, simple_header_lines=None, verbose=False
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

    log_message("Applied CLI metadata to merged header.", verbose)
    return header


def _validate_anonymized_vcf_header(final_vcf_path, ensure_for_uncompressed=False):
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
    merged_vcf,
    sample_metadata_entries=None,
    header_metadata_lines=None,
    serialized_sample_line=None,
    verbose=False,
):
    """Finalize the merged VCF by applying bcftools processing and metadata."""

    header_metadata_lines = header_metadata_lines or []

    joint_temp = f"{merged_vcf}.joint.temp.vcf"
    filtered_temp = f"{merged_vcf}.filtered.temp.vcf"
    header_temp = f"{merged_vcf}.header.temp"
    body_temp = f"{merged_vcf}.body.temp"

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
        for path in [joint_temp, filtered_temp, header_temp, body_temp]:
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
        with open(header_temp, "w", encoding="utf-8") as header_handle:
            subprocess.run(
                ["bcftools", "view", "-h", filtered_temp],
                check=True,
                stdout=header_handle,
            )

        with open(body_temp, "w", encoding="utf-8") as body_handle:
            view_proc = subprocess.Popen(
                ["bcftools", "view", "-H", filtered_temp],
                stdout=subprocess.PIPE,
            )
            try:
                subprocess.run(
                    ["cut", "-f1-8"],
                    check=True,
                    stdin=view_proc.stdout,
                    stdout=body_handle,
                )
            finally:
                if view_proc.stdout is not None:
                    view_proc.stdout.close()
                return_code = view_proc.wait()
            if return_code != 0:
                raise subprocess.CalledProcessError(
                    return_code, ["bcftools", "view", "-H", filtered_temp]
                )
    except (subprocess.CalledProcessError, OSError) as exc:
        _cleanup_temp_files()
        handle_critical_error(
            f"Failed to anonymize merged VCF columns using bcftools: {exc}"
        )

    log_message("Combining custom metadata header with VCF header...", verbose)
    try:
        with open(header_temp, "r", encoding="utf-8") as header_handle:
            header_lines = [line.rstrip("\n") for line in header_handle]

        fileformat_line = None
        remaining_header_lines = []
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

        final_header_lines = []
        if fileformat_line is not None:
            final_header_lines.append(fileformat_line)

        final_header_lines.extend(remaining_header_lines)
        existing_header_lines = set(final_header_lines)

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
                _cleanup_temp_files()
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

        final_header_lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

        body_lines = []
        with open(body_temp, "r", encoding="utf-8") as body_handle:
            for line in body_handle:
                body_lines.append(line.rstrip("\n"))

        is_gzipped_output = False
        if merged_vcf.endswith(".vcf.gz"):
            base = merged_vcf[: -len(".vcf.gz")]
            final_vcf = f"{base}.anonymized.vcf.gz"
            is_gzipped_output = True
        else:
            base, ext = os.path.splitext(merged_vcf)
            if not ext:
                ext = ".vcf"
            final_vcf = f"{base}.anonymized{ext}"

        writer = gzip.open if is_gzipped_output else open
        open_kwargs = {"mode": "wt", "encoding": "utf-8"} if is_gzipped_output else {"mode": "w", "encoding": "utf-8"}
        with writer(final_vcf, **open_kwargs) as output_handle:
            for line in final_header_lines:
                output_handle.write(line + "\n")
            for line in body_lines:
                output_handle.write(line + "\n")
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

    if expects_gzip:
        try:
            subprocess.run(["bgzip", "-f", final_plain_vcf], check=True)
            subprocess.run(["tabix", "-p", "vcf", "-f", final_vcf], check=True)
        except (subprocess.CalledProcessError, OSError) as exc:
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
                f"Failed to finalize anonymized VCF compression/indexing: {exc}"
            )
        finally:
            if os.path.exists(final_plain_vcf):
                try:
                    os.remove(final_plain_vcf)
                except Exception:
                    pass
    else:
        final_vcf = final_plain_vcf

    _cleanup_temp_files()

    _validate_anonymized_vcf_header(final_vcf, ensure_for_uncompressed=True)

    log_message(f"Anonymized merged VCF written to: {final_vcf}", verbose)
    return final_vcf


def union_headers(valid_files, sample_order=None):
    """Return a merged header with combined metadata from *valid_files*."""

    combined_header = None
    info_ids = set()
    filter_ids = set()
    contig_ids = set()
    computed_sample_order = []
    merged_sample_metadata = None

    for file_path in valid_files:
        preprocessed_file = preprocess_vcf(file_path)
        reader = None
        try:
            reader = vcfpy.Reader.from_path(preprocessed_file)
            header = reader.header

            if combined_header is None:
                combined_header = header.copy()
                for line in combined_header.lines:
                    if isinstance(line, vcfpy.header.InfoHeaderLine):
                        info_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.FilterHeaderLine):
                        filter_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.ContigHeaderLine):
                        contig_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.SampleHeaderLine):
                        mapping = OrderedDict(getattr(line, "mapping", {}))
                        if mapping.get("ID"):
                            merged_sample_metadata = OrderedDict(mapping)

                _remove_format_and_sample_definitions(combined_header)
                if hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
                    computed_sample_order = list(combined_header.samples.names)
                else:
                    computed_sample_order = []
                continue

            if hasattr(header, "samples") and hasattr(header.samples, "names"):
                for sample_name in header.samples.names:
                    if sample_name not in computed_sample_order:
                        computed_sample_order.append(sample_name)

            for line in header.lines:
                if isinstance(line, vcfpy.header.InfoHeaderLine):
                    if line.id not in info_ids:
                        combined_header.add_line(copy.deepcopy(line))
                        info_ids.add(line.id)
                elif isinstance(line, vcfpy.header.FilterHeaderLine):
                    if line.id not in filter_ids:
                        combined_header.add_line(copy.deepcopy(line))
                        filter_ids.add(line.id)
                elif isinstance(line, vcfpy.header.ContigHeaderLine):
                    if line.id not in contig_ids:
                        combined_header.add_line(copy.deepcopy(line))
                        contig_ids.add(line.id)
                elif isinstance(line, vcfpy.header.SampleHeaderLine):
                    mapping = OrderedDict(getattr(line, "mapping", {}))
                    sample_id = mapping.get("ID")
                    if not sample_id:
                        continue
                    if merged_sample_metadata is None:
                        merged_sample_metadata = OrderedDict(mapping)
                        continue
                    for key, value in mapping.items():
                        if key not in merged_sample_metadata or not merged_sample_metadata[key]:
                            merged_sample_metadata[key] = value
        except Exception as exc:
            handle_critical_error(f"Failed to read VCF header from {file_path}: {exc}")
        finally:
            if reader is not None:
                try:
                    reader.close()
                except Exception:
                    pass
            if preprocessed_file != file_path and os.path.exists(preprocessed_file):
                os.remove(preprocessed_file)

    if combined_header is None:
        handle_critical_error("Unable to construct a merged VCF header.")

    target_sample_order = sample_order if sample_order is not None else computed_sample_order
    if target_sample_order and hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
        combined_header.samples.names = list(target_sample_order)

    _remove_format_and_sample_definitions(combined_header)

    if merged_sample_metadata:
        serialized = build_sample_metadata_line(merged_sample_metadata)
        parsed_mapping = _parse_sample_metadata_line(serialized)
        sample_line = vcfpy.header.SampleHeaderLine.from_mapping(parsed_mapping)
        combined_header.lines = [
            line for line in combined_header.lines if getattr(line, "key", None) != "SAMPLE"
        ]
        combined_header.add_line(sample_line)

    return combined_header


def _create_missing_call_factory(format_keys, header):
    if not format_keys:
        return lambda: {}

    template = {}
    for key in format_keys:
        if key == "GT":
            template[key] = "./."
            continue

        field_info = None
        try:
            field_info = header.get_format_field_info(key)
        except Exception:
            field_info = None

        number = getattr(field_info, "number", None)
        if isinstance(number, str):
            stripped = number.strip()
            if not stripped:
                number = None
            elif stripped in {".", "A", "G", "R"}:
                template[key] = []
                continue
            else:
                try:
                    number = int(stripped)
                except ValueError:
                    number = None

        if isinstance(number, int):
            if number <= 1:
                default_value = None
            else:
                default_value = [None] * number
        elif number in {".", "A", "G", "R"}:
            default_value = []
        else:
            default_value = None

        template.setdefault(key, default_value)

    def factory():
        data = {}
        for fmt_key, default in template.items():
            if isinstance(default, list):
                data[fmt_key] = [None] * len(default)
            else:
                data[fmt_key] = default
        return data

    return factory


def _remove_format_and_sample_definitions(header):
    """Strip FORMAT definitions and sample columns from a VCF header."""

    if header is None:
        return

    if hasattr(header, "lines"):
        filtered_lines = []
        for line in header.lines:
            if isinstance(line, vcfpy.header.FormatHeaderLine):
                continue
            key = getattr(line, "key", None)
            if isinstance(key, str) and key.upper() == "FORMAT":
                continue
            filtered_lines.append(line)
        header.lines = filtered_lines

    if hasattr(header, "formats"):
        try:
            header.formats.clear()
        except AttributeError:
            header.formats = OrderedDict()

    if hasattr(header, "samples") and hasattr(header.samples, "names"):
        header.samples.names = []


def _pad_record_samples(record, header, sample_order):
    if not sample_order or not hasattr(record, "call_for_sample"):
        return

    format_keys = list(record.FORMAT or [])
    factory = _create_missing_call_factory(format_keys, header)
    updated_calls = []
    for name in sample_order:
        call = record.call_for_sample.get(name)
        if call is None:
            call = vcfpy.Call(name, factory())
        else:
            defaults = factory()
            for key in format_keys:
                if key not in call.data:
                    call.data[key] = defaults[key]
                    continue

                default_value = defaults[key]
                current_value = call.data[key]
                if isinstance(default_value, list) and not isinstance(current_value, list):
                    if current_value is None:
                        call.data[key] = default_value
                    elif isinstance(current_value, tuple):
                        call.data[key] = list(current_value)
                    else:
                        call.data[key] = [current_value]

            if "GT" in defaults:
                genotype_value = call.data.get("GT")
                if genotype_value in {None, "", "."}:
                    call.data["GT"] = defaults["GT"]
        updated_calls.append(call)
    record.update_calls(updated_calls)


import os, shutil, subprocess, datetime

def merge_vcfs(
    valid_files,
    output_dir,
    verbose=False,
    sample_order=None,
    sample_header_line=None,
    simple_header_lines=None,
):
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_vcf = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")   # plain VCF
    filled_vcf = base_vcf + ".filled.vcf"                                 # after +fill-tags
    gz_vcf   = base_vcf + ".gz"

    file_count = len(valid_files)
    log_message(f"Found {file_count} VCF files. Merging them into a combined VCF...", verbose)

    temp_files, preprocessed_files = [], []
    for file_path in valid_files:
        try:
            pre = preprocess_vcf(file_path)   # your helper prepares bgzipped+indexed or plain inputs
        except Exception as exc:
            handle_critical_error(f"Failed to preprocess {file_path}: {exc}")
        preprocessed_files.append(pre)
        if pre != file_path:
            temp_files.append(pre)

    log_message("Merging VCF files with bcftools...", verbose)
    try:
        # Write UNCOMPRESSED VCF so vcfpy can adjust the header easily
        result = subprocess.run(
            ["bcftools", "merge", "-m", "all", "-Ov", "-o", base_vcf, *preprocessed_files],
            capture_output=True, text=True,
        )
    finally:
        for tmp in temp_files:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass
    if result.returncode != 0:
        handle_critical_error((result.stderr or "bcftools merge failed").strip())

    # Open merged VCF, apply metadata/header edits, rewrite
    try:
        reader = vcfpy.Reader.from_path(base_vcf)
    except Exception as exc:
        handle_critical_error(f"Failed to open merged VCF for post-processing: {exc}")

    try:
        header = apply_metadata_to_header(
            reader.header,
            sample_header_line=sample_header_line,
            simple_header_lines=simple_header_lines,
            verbose=verbose,
        )
    except SystemExit:
        reader.close()
        raise
    except Exception as exc:
        reader.close()
        handle_critical_error(f"Failed to apply metadata to merged VCF: {exc}")

    tmp_out = base_vcf + ".tmp"
    try:
        writer = vcfpy.Writer.from_path(tmp_out, header)  # writes full header immediately
    except Exception as exc:
        reader.close()
        handle_critical_error(f"Failed to open temporary writer for merged VCF: {exc}")

    try:
        for record in reader:
            # optional: _pad_record_samples(record, header, sample_order)
            writer.write_record(record)
    finally:
        reader.close()
        writer.close()

    shutil.move(tmp_out, base_vcf)

    # Recalculate AC/AN/AF using bcftools +fill-tags
    log_message("Recomputing AC, AN, AF with bcftools +fill-tags...", verbose)
    res2 = subprocess.run(
        ["bcftools", "+fill-tags", base_vcf, "-Ov", "-o", filled_vcf, "--", "-t", "AC,AN,AF"],
        capture_output=True, text=True,
    )
    if res2.returncode != 0:
        handle_critical_error((res2.stderr or "bcftools +fill-tags failed").strip())
    shutil.move(filled_vcf, base_vcf)

    # Compress (BGZF) and index (tabix)
    log_message("Compressing and indexing the final VCF...", verbose)
    try:
        subprocess.run(["bgzip", "-f", base_vcf], check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", gz_vcf], check=True)
    except subprocess.CalledProcessError as exc:
        handle_critical_error(f"Failed to compress or index merged VCF ({base_vcf}): {exc}")

    log_message(f"Merged VCF file created and indexed successfully: {gz_vcf}", verbose)
    return gz_vcf

def _parse_simple_metadata_line(line):
    stripped = line.strip()
    if not stripped.startswith("##") or "=" not in stripped:
        return None
    key, value = stripped[2:].split("=", 1)
    key = key.strip()
    value = value.strip()
    if not key:
        return None
    return key, value


def preprocess_vcf(file_path):
    """
    Check if the VCF file uses spaces instead of tabs for the column header and data lines.
    If so, create a temporary file where the column header (#CHROM) and all subsequent lines 
    are converted to be tab-delimited. Lines starting with "##" (metadata) are left unchanged.
    Returns the path to the file to be used (original or temporary).
    """
    import re
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    modified = False
    new_lines = []
    header_found = False  # indicates when the column header has been encountered
    for line in lines:
        if line.startswith("##"):
            # Do not change metadata lines (they may contain spaces that are part of the value)
            new_lines.append(line)
        elif line.startswith("#"):
            # This is the column header line (e.g. "#CHROM ...")
            new_line = re.sub(r'\s+', '\t', line.rstrip()) + "\n"
            new_lines.append(new_line)
            header_found = True
            if new_line != line:
                modified = True
        else:
            # Data lines: once header_found is True, convert spaces to tabs.
            if header_found:
                new_line = re.sub(r'\s+', '\t', line.rstrip()) + "\n"
                new_lines.append(new_line)
                if new_line != line:
                    modified = True
            else:
                new_lines.append(line)
    
    if modified:
        temp_file = file_path + ".tmp"
        with open(temp_file, 'w') as f:
            f.writelines(new_lines)
        return temp_file
    else:
        return file_path



def summarize_produced_vcfs(output_dir, fallback_vcf):
    """Return a directory, representative filename, and count for output VCFs."""

    discovered = []
    output_dir = os.path.abspath(output_dir) if output_dir else None
    if output_dir and os.path.isdir(output_dir):
        for pattern in ("SAMPLE_*.vcf.gz", "SAMPLE_*.vcf"):
            discovered.extend(
                sorted(
                    glob.glob(os.path.join(output_dir, pattern)),
                    key=lambda path: os.path.basename(path),
                )
            )

    if discovered:
        # Prefer gzipped outputs when both compressed and uncompressed variants exist.
        discovered.sort(
            key=lambda path: (
                0 if path.endswith(".vcf.gz") else 1,
                os.path.basename(path),
            )
        )
        representative_path = discovered[0]
        directory = os.path.dirname(os.path.abspath(representative_path))
        filename = os.path.basename(representative_path)
        count = len(discovered)
        return directory, filename, count

    fallback_abs = os.path.abspath(fallback_vcf)
    return os.path.dirname(fallback_abs), os.path.basename(fallback_abs), 1


def main():
    args = parse_arguments()
    verbose = args.verbose
    (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries,
        sanitized_header_lines,
        serialized_sample_line,
    ) = parse_metadata_arguments(args, verbose)

    log_message(
        "Script Execution Log - " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        verbose,
    )

    input_dir = os.path.abspath(args.input_dir)
    if not os.path.isdir(input_dir):
        handle_critical_error(f"Input directory does not exist: {input_dir}")
    log_message("Input directory: " + input_dir, verbose)

    output_dir = os.path.abspath(args.output_dir) if args.output_dir else input_dir
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    log_message("Output directory: " + output_dir, verbose)

    ref_genome = args.ref.strip() if args.ref else None
    vcf_version = normalize_vcf_version(args.vcf_version) if args.vcf_version else None

    detection_needed = not ref_genome or not vcf_version
    detected_file = detected_fileformat = detected_reference = None
    if detection_needed:
        detected_file, detected_fileformat, detected_reference = find_first_vcf_with_header(input_dir, verbose)
        if not ref_genome and detected_reference:
            ref_genome = detected_reference
        if not vcf_version and detected_fileformat:
            vcf_version = normalize_vcf_version(detected_fileformat)

    if not ref_genome:
        handle_critical_error("Reference genome build must be provided via --ref or auto-detectable from input files.")
    if not vcf_version:
        handle_critical_error("VCF version must be provided via --vcf-version or auto-detectable from input files.")

    if detection_needed and detected_file:
        log_message(
            f"Auto-detected metadata from {detected_file} -> reference={detected_reference or 'unknown'}, "
            f"version={detected_fileformat or 'unknown'}",
            verbose,
        )

    log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}", verbose)

    valid_files, sample_order = validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose)

    merged_vcf = merge_vcfs(
        valid_files,
        output_dir,
        verbose,
        sample_order=sample_order,
        sample_header_line=sample_header_line,
        simple_header_lines=simple_header_lines,
    )

    final_vcf = append_metadata_to_merged_vcf(
        merged_vcf,
        sample_metadata_entries=sample_metadata_entries,
        header_metadata_lines=sanitized_header_lines,
        serialized_sample_line=serialized_sample_line,
        verbose=verbose,
    )

    validate_merged_vcf(final_vcf, verbose)

    summary_dir, sample_filename, produced_count = summarize_produced_vcfs(
        output_dir, final_vcf
    )
    summary_path = os.path.join(summary_dir, sample_filename)
    print(f"Wrote: {summary_path} x {produced_count}")
    log_message(
        "Script execution completed successfully. "
        f"Wrote {produced_count} VCF file(s) (e.g., {summary_path}).",
        verbose,
    )



if __name__ == "__main__":
    main()
