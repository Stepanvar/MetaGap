"""Routines for combining VCF shards into a merged cohort file."""

from __future__ import annotations

import copy
import datetime
import os
import re
import shutil
from collections import OrderedDict
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import pysam, vcfpy
from .logging_utils import handle_critical_error, log_message
from .metadata import (
    _parse_sample_metadata_line,
    apply_metadata_to_header,
    build_sample_metadata_line,
)


def preprocess_vcf(file_path: str) -> str:
    """Normalize whitespace delimiters in ``file_path`` if necessary."""

    with open(file_path, "r", encoding="utf-8") as handle:
        lines = handle.readlines()

    modified = False
    new_lines: List[str] = []
    header_found = False
    for line in lines:
        if line.startswith("##"):
            new_lines.append(line)
        elif line.startswith("#"):
            new_line = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
            new_lines.append(new_line)
            header_found = True
            if new_line != line:
                modified = True
        else:
            if header_found:
                new_line = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
                new_lines.append(new_line)
                if new_line != line:
                    modified = True
            else:
                new_lines.append(line)

    if modified:
        temp_file = file_path + ".tmp"
        with open(temp_file, "w", encoding="utf-8") as handle:
            handle.writelines(new_lines)
        return temp_file
    return file_path


def _create_missing_call_factory(format_keys: Sequence[str], header) -> Callable[[], dict]:
    if not format_keys:
        return lambda: {}

    template = {}
    for key in format_keys:
        if key == "GT":
            template[key] = "./."
            continue

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

    def factory() -> dict:
        data = {}
        for fmt_key, default in template.items():
            if isinstance(default, list):
                data[fmt_key] = [None] * len(default)
            else:
                data[fmt_key] = default
        return data

    return factory


def _remove_format_and_sample_definitions(header) -> None:
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


def _pad_record_samples(record, header, sample_order: Sequence[str]) -> None:
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
                    call.data[key] = defaults.get(key)
                else:
                    current_value = call.data[key]
                    if isinstance(current_value, list):
                        continue
                    if isinstance(current_value, tuple):
                        call.data[key] = list(current_value)
                    else:
                        call.data[key] = [current_value]

            if "GT" in defaults:
                genotype_value = call.data.get("GT")
                if genotype_value in {None, "", "."}:
                    call.data["GT"] = defaults["GT"]
        updated_calls.append(call)
    record.update_calls(updated_calls)


def merge_vcfs(
    valid_files: Sequence[str],
    output_dir: str,
    verbose: bool = False,
    sample_order: Optional[Sequence[str]] = None,
    sample_header_line=None,
    simple_header_lines=None,
) -> str:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_vcf = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")
    gz_vcf = base_vcf + ".gz"

    file_count = len(valid_files)
    log_message(
        f"Found {file_count} VCF files. Merging them into a combined VCF...",
        verbose,
    )

    temp_files: List[str] = []
    preprocessed_files: List[str] = []
    for file_path in valid_files:
        try:
            pre = preprocess_vcf(file_path)
        except Exception as exc:
            handle_critical_error(f"Failed to preprocess {file_path}: {exc}")
        preprocessed_files.append(pre)
        if pre != file_path:
            temp_files.append(pre)

    merged_header = None
    info_ids: set[str] = set()
    format_ids: set[str] = set()
    filter_ids: set[str] = set()
    contig_ids: set[str] = set()
    discovered_samples: List[str] = []

    record_store: "OrderedDict[Tuple[str, int, str, Tuple[str, ...]], Dict[str, object]]" = OrderedDict()

    def _record_key(record) -> Tuple[str, int, str, Tuple[str, ...]]:
        alts = tuple(str(alt) for alt in getattr(record, "ALT", []) or [])
        return (str(record.CHROM), int(record.POS), str(record.REF), alts)

    def _initialise_store(record) -> Dict[str, object]:
        return {
            "template": copy.deepcopy(record),
            "format_keys": list(getattr(record, "FORMAT", []) or []),
            "calls": {},
        }

    try:
        for file_path in preprocessed_files:
            try:
                reader = vcfpy.Reader.from_path(file_path)
            except Exception as exc:
                handle_critical_error(f"Failed to read VCF file {file_path}: {exc}")

            header = reader.header

            if merged_header is None:
                merged_header = copy.deepcopy(header)
                for line in getattr(merged_header, "lines", []):
                    if isinstance(line, vcfpy.header.InfoHeaderLine) and line.id:
                        info_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.FormatHeaderLine) and line.id:
                        format_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.FilterHeaderLine) and line.id:
                        filter_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.ContigHeaderLine) and line.id:
                        contig_ids.add(line.id)
                discovered_samples.extend(getattr(header.samples, "names", []))
            else:
                for line in getattr(header, "lines", []):
                    try:
                        cloned = copy.deepcopy(line)
                    except Exception:
                        cloned = line
                    if isinstance(line, vcfpy.header.InfoHeaderLine) and line.id:
                        if line.id not in info_ids:
                            merged_header.add_line(cloned)
                            info_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.FormatHeaderLine) and line.id:
                        if line.id not in format_ids:
                            merged_header.add_line(cloned)
                            format_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.FilterHeaderLine) and line.id:
                        if line.id not in filter_ids:
                            merged_header.add_line(cloned)
                            filter_ids.add(line.id)
                    elif isinstance(line, vcfpy.header.ContigHeaderLine) and line.id:
                        if line.id not in contig_ids:
                            merged_header.add_line(cloned)
                            contig_ids.add(line.id)

                for sample_name in getattr(header.samples, "names", []):
                    if sample_name not in discovered_samples:
                        discovered_samples.append(sample_name)

            for record in reader:
                key = _record_key(record)
                container = record_store.get(key)
                if container is None:
                    container = _initialise_store(record)
                    record_store[key] = container
                else:
                    format_keys = container["format_keys"]
                    for fmt_key in getattr(record, "FORMAT", []) or []:
                        if fmt_key not in format_keys:
                            format_keys.append(fmt_key)

                call_mapping: Dict[str, dict] = container["calls"]  # type: ignore[assignment]
                for call in getattr(record, "calls", []):
                    sample_name = getattr(call, "sample", getattr(call, "name", None))
                    if sample_name is None:
                        continue
                    call_mapping[sample_name] = copy.deepcopy(getattr(call, "data", {}))

            reader.close()
    finally:
        for tmp in temp_files:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass

    if merged_header is None:
        handle_critical_error("Failed to construct a merged VCF header from the provided files.")

    final_sample_order = list(sample_order) if sample_order is not None else discovered_samples
    if hasattr(merged_header, "samples") and hasattr(merged_header.samples, "names"):
        merged_header.samples.names = list(final_sample_order)

    try:
        header = apply_metadata_to_header(
            merged_header,
            sample_header_line=sample_header_line,
            simple_header_lines=simple_header_lines,
            verbose=verbose,
        )
    except SystemExit:
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to apply metadata to merged VCF header: {exc}")

    tmp_out = base_vcf + ".tmp"
    try:
        writer = vcfpy.Writer.from_path(tmp_out, header)
    except Exception as exc:
        handle_critical_error(f"Failed to create VCF writer for merged output: {exc}")

    try:
        for container in record_store.values():
            template = copy.deepcopy(container["template"])  # type: ignore[index]
            format_keys = list(container["format_keys"])  # type: ignore[index]
            template.FORMAT = format_keys

            factory = _create_missing_call_factory(format_keys, header)
            default_values = factory()
            updated_calls = []
            call_mapping = container["calls"]  # type: ignore[index]
            for sample_name in final_sample_order:
                sample_data = call_mapping.get(sample_name)
                data = {}
                for fmt_key in format_keys:
                    default_value = copy.deepcopy(default_values.get(fmt_key))
                    if sample_data is not None and fmt_key in sample_data:
                        value = sample_data[fmt_key]
                        if isinstance(value, tuple):
                            value = list(value)
                        elif isinstance(value, list):
                            value = list(value)
                        data[fmt_key] = value
                    else:
                        data[fmt_key] = default_value
                call = vcfpy.Call(sample_name, data)
                updated_calls.append(call)

            template.update_calls(updated_calls)
            writer.write_record(template)
    finally:
        writer.close()

    shutil.move(tmp_out, base_vcf)

    log_message("Compressing and indexing the final VCF with pysam...", verbose)
    try:
        pysam.tabix_compress(base_vcf, gz_vcf, force=True)
        pysam.tabix_index(gz_vcf, preset="vcf", force=True)
        if os.path.exists(base_vcf):
            os.remove(base_vcf)
    except Exception as exc:
        handle_critical_error(f"Failed to compress or index merged VCF ({base_vcf}): {exc}")

    log_message(f"Merged VCF file created and indexed successfully: {gz_vcf}", verbose)
    return gz_vcf


def union_headers(valid_files: Sequence[str], sample_order: Optional[Sequence[str]] = None):
    """Return a merged header with combined metadata from *valid_files*."""

    combined_header = None
    info_ids = set()
    filter_ids = set()
    contig_ids = set()
    computed_sample_order: List[str] = []
    merged_sample_metadata = None

    initial_sample_names: List[str] = []

    for file_path in valid_files:
        preprocessed_file = preprocess_vcf(file_path)
        reader = None
        try:
            reader = vcfpy.Reader.from_path(preprocessed_file)
            header = reader.header

            if combined_header is None:
                combined_header = header.copy()
                if hasattr(combined_header, "samples") and hasattr(
                    combined_header.samples, "names"
                ):
                    initial_sample_names = list(combined_header.samples.names)
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
                computed_sample_order = list(initial_sample_names)
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
    samples_attr_present = hasattr(combined_header, "samples") and hasattr(
        combined_header.samples, "names"
    )
    if target_sample_order and samples_attr_present:
        combined_header.samples.names = list(target_sample_order)

    _remove_format_and_sample_definitions(combined_header)

    if target_sample_order and samples_attr_present:
        combined_header.samples.names = list(target_sample_order)

    if merged_sample_metadata:
        serialized = build_sample_metadata_line(merged_sample_metadata)
        parsed_mapping = _parse_sample_metadata_line(serialized)
        sample_line = vcfpy.header.SampleHeaderLine.from_mapping(parsed_mapping)
        combined_header.lines = [
            line for line in combined_header.lines if getattr(line, "key", None) != "SAMPLE"
        ]
        combined_header.add_line(sample_line)

    return combined_header


__all__ = [
    "preprocess_vcf",
    "merge_vcfs",
    "union_headers",
    "_create_missing_call_factory",
    "_remove_format_and_sample_definitions",
    "_pad_record_samples",
]
