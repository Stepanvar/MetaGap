"""Routines for combining VCF shards into a merged cohort file."""

from __future__ import annotations

import copy
import datetime
import os
import re
import shutil
import subprocess
from collections import OrderedDict
from typing import Callable, List, Optional, Sequence

from . import vcfpy
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


def _record_passes_filters(
    record,
    qual_threshold: Optional[float],
    an_threshold: Optional[float],
    allowed_filter_values: Optional[Sequence[str]],
) -> bool:
    if qual_threshold is not None:
        qual_value = getattr(record, "QUAL", None)
        if qual_value is None:
            return False
        try:
            numeric_qual = float(qual_value)
        except (TypeError, ValueError):
            return False
        if numeric_qual <= qual_threshold:
            return False

    if an_threshold is not None:
        an_value = record.INFO.get("AN")
        if isinstance(an_value, list):
            an_value = an_value[0] if an_value else None
        try:
            numeric_an = float(an_value) if an_value is not None else None
        except (TypeError, ValueError):
            numeric_an = None
        if numeric_an is None or numeric_an <= an_threshold:
            return False

    filters = record.FILTER or []
    if not filters:
        return True

    normalized_filters = []
    for entry in filters:
        if entry in {None, "", "."}:
            continue
        normalized_filters.append(entry)

    if not normalized_filters:
        return True

    if allowed_filter_values is None:
        allowed_filter_values = ("PASS",)

    allowed_set = {value for value in allowed_filter_values if value not in {None, "", "."}}
    return bool(allowed_set) and all(value in allowed_set for value in normalized_filters)


def _filter_vcf_records(
    input_path: str,
    qual_threshold: Optional[float],
    an_threshold: Optional[float],
    allowed_filter_values: Optional[Sequence[str]],
    verbose: bool,
) -> None:
    try:
        reader = vcfpy.Reader.from_path(input_path)
    except Exception as exc:
        handle_critical_error(f"Failed to read VCF for filtering ({input_path}): {exc}")

    tmp_filtered = input_path + ".filtered"
    kept_records = 0
    total_records = 0

    try:
        writer = vcfpy.Writer.from_path(tmp_filtered, reader.header)
    except Exception as exc:
        reader.close()
        handle_critical_error(f"Failed to open filtered VCF writer ({input_path}): {exc}")

    try:
        for record in reader:
            total_records += 1
            if _record_passes_filters(
                record,
                qual_threshold=qual_threshold,
                an_threshold=an_threshold,
                allowed_filter_values=allowed_filter_values,
            ):
                writer.write_record(record)
                kept_records += 1
    finally:
        reader.close()
        writer.close()

    shutil.move(tmp_filtered, input_path)

    log_message(
        (
            "Applied in-memory variant quality filter: "
            f"kept {kept_records} of {total_records} records."
        ),
        verbose,
    )


def merge_vcfs(
    valid_files: Sequence[str],
    output_dir: str,
    verbose: bool = False,
    sample_order: Optional[Sequence[str]] = None,
    sample_header_line=None,
    simple_header_lines=None,
    qual_threshold: Optional[float] = 30.0,
    an_threshold: Optional[float] = 50.0,
    allowed_filter_values: Optional[Sequence[str]] = ("PASS",),
) -> str:
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_vcf = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")
    filled_vcf = base_vcf + ".filled.vcf"
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

    log_message("Merging VCF files with bcftools...", verbose)
    try:
        result = subprocess.run(
            ["bcftools", "merge", "-m", "all", "-Ov", "-o", base_vcf, *preprocessed_files],
            capture_output=True,
            text=True,
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
        writer = vcfpy.Writer.from_path(tmp_out, header)
    except Exception as exc:
        reader.close()
        handle_critical_error(f"Failed to open temporary writer for merged VCF: {exc}")

    try:
        for record in reader:
            writer.write_record(record)
    finally:
        reader.close()
        writer.close()

    shutil.move(tmp_out, base_vcf)

    log_message("Recomputing AC, AN, AF with bcftools +fill-tags...", verbose)
    res2 = subprocess.run(
        ["bcftools", "+fill-tags", base_vcf, "-Ov", "-o", filled_vcf, "--", "-t", "AC,AN,AF"],
        capture_output=True,
        text=True,
    )
    if res2.returncode != 0:
        handle_critical_error((res2.stderr or "bcftools +fill-tags failed").strip())
    shutil.move(filled_vcf, base_vcf)

    _filter_vcf_records(
        base_vcf,
        qual_threshold=qual_threshold,
        an_threshold=an_threshold,
        allowed_filter_values=allowed_filter_values,
        verbose=verbose,
    )

    log_message("Compressing and indexing the final VCF...", verbose)
    try:
        subprocess.run(["bgzip", "-f", base_vcf], check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", gz_vcf], check=True)
    except subprocess.CalledProcessError as exc:
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
    "_filter_vcf_records",
    "_record_passes_filters",
]
