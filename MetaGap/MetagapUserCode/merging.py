"""Utilities for merging VCF files using shared callbacks from the CLI script."""
from __future__ import annotations

import copy
import datetime
import os
import re
import shutil
import subprocess
from collections import OrderedDict
from dataclasses import dataclass
from typing import Callable, Mapping, MutableMapping, Optional


@dataclass
class MergeCallbacks:
    """Container for helper callables required by the merging helpers."""

    preprocess_vcf: Callable[[str], str]
    log_message: Callable[[str, bool], None]
    handle_critical_error: Callable[[str], None]
    build_sample_metadata_line: Callable[["OrderedDict[str, str]"], str]
    parse_sample_metadata_line: Callable[[str], Mapping[str, str]]
    apply_metadata_to_header: Callable[..., "object"]
    vcfpy: "object"


_callbacks: Optional[MergeCallbacks] = None


def configure_callbacks(callbacks: MergeCallbacks) -> None:
    """Register helper callbacks that are required for merging operations."""

    global _callbacks
    _callbacks = callbacks


def _require_callbacks() -> MergeCallbacks:
    if _callbacks is None:
        raise RuntimeError(
            "Merging callbacks were not configured. Call configure_callbacks() before using these helpers."
        )
    return _callbacks


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


def _format_info_definition(info_id, definition_mapping):
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


def ensure_standard_info_definitions(header, verbose=False):
    callbacks = _require_callbacks()
    vcfpy = callbacks.vcfpy
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
            callbacks.log_message(
                f"WARNING: Unable to create INFO header definition for {info_id}: {exc}",
                verbose,
            )
            continue

        try:
            header.add_line(info_line)
        except Exception as exc:  # pragma: no cover - defensive logging
            callbacks.log_message(
                f"WARNING: Failed to append INFO header definition for {info_id}: {exc}",
                verbose,
            )
            continue

        existing_ids.add(info_id)
        added_ids.append(info_id)

    if added_ids:
        callbacks.log_message(
            "Inserted INFO header definitions for missing fields: "
            + ", ".join(added_ids),
            verbose,
        )

    return header


def ensure_standard_info_header_lines(final_header_lines, existing_header_lines, verbose=False):
    callbacks = _require_callbacks()

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
        callbacks.log_message(
            "Inserted INFO header definitions for missing fields: "
            + ", ".join(added_ids),
            verbose,
        )

    return final_header_lines


def union_headers(valid_files, sample_order=None):
    """Return a merged header with combined metadata from *valid_files*."""

    callbacks = _require_callbacks()
    vcfpy = callbacks.vcfpy

    combined_header = None
    info_ids = set()
    filter_ids = set()
    contig_ids = set()
    computed_sample_order = []
    merged_sample_metadata: Optional[OrderedDict] = None

    for file_path in valid_files:
        preprocessed_file = callbacks.preprocess_vcf(file_path)
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
        except Exception as exc:  # pragma: no cover - propagated to CLI handler
            callbacks.handle_critical_error(f"Failed to read VCF header from {file_path}: {exc}")
        finally:
            if reader is not None:
                try:
                    reader.close()
                except Exception:  # pragma: no cover - defensive cleanup
                    pass
            if preprocessed_file != file_path and os.path.exists(preprocessed_file):
                os.remove(preprocessed_file)

    if combined_header is None:
        callbacks.handle_critical_error("Unable to construct a merged VCF header.")

    target_sample_order = sample_order if sample_order is not None else computed_sample_order
    if target_sample_order and hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
        combined_header.samples.names = list(target_sample_order)

    _remove_format_and_sample_definitions(combined_header)

    if merged_sample_metadata:
        serialized = callbacks.build_sample_metadata_line(merged_sample_metadata)
        parsed_mapping = callbacks.parse_sample_metadata_line(serialized)
        sample_line = vcfpy.header.SampleHeaderLine.from_mapping(parsed_mapping)
        combined_header.lines = [
            line for line in combined_header.lines if getattr(line, "key", None) != "SAMPLE"
        ]
        combined_header.add_line(sample_line)

    return combined_header


def _create_missing_call_factory(format_keys, header):
    callbacks = _require_callbacks()
    vcfpy = callbacks.vcfpy

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
        data: MutableMapping[str, object] = {}
        for fmt_key, default in template.items():
            if isinstance(default, list):
                data[fmt_key] = [None] * len(default)
            else:
                data[fmt_key] = default
        return data

    return factory


def _remove_format_and_sample_definitions(header):
    """Strip FORMAT definitions and sample columns from a VCF header."""

    callbacks = _require_callbacks()
    vcfpy = callbacks.vcfpy

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
    callbacks = _require_callbacks()
    vcfpy = callbacks.vcfpy

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


def merge_vcfs(
    valid_files,
    output_dir,
    verbose=False,
    sample_order=None,
    sample_header_line=None,
    simple_header_lines=None,
):
    callbacks = _require_callbacks()
    vcfpy = callbacks.vcfpy

    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_vcf = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")
    filled_vcf = base_vcf + ".filled.vcf"
    gz_vcf = base_vcf + ".gz"

    file_count = len(valid_files)
    callbacks.log_message(
        f"Found {file_count} VCF files. Merging them into a combined VCF...",
        verbose,
    )

    temp_files = []
    preprocessed_files = []
    for file_path in valid_files:
        try:
            pre = callbacks.preprocess_vcf(file_path)
        except Exception as exc:  # pragma: no cover - propagated to CLI
            callbacks.handle_critical_error(f"Failed to preprocess {file_path}: {exc}")
        preprocessed_files.append(pre)
        if pre != file_path:
            temp_files.append(pre)

    callbacks.log_message("Merging VCF files with bcftools...", verbose)
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
            except OSError:  # pragma: no cover - cleanup best effort
                pass
    if result.returncode != 0:
        callbacks.handle_critical_error((result.stderr or "bcftools merge failed").strip())

    try:
        reader = vcfpy.Reader.from_path(base_vcf)
    except Exception as exc:  # pragma: no cover - propagated to CLI
        callbacks.handle_critical_error(f"Failed to open merged VCF for post-processing: {exc}")

    try:
        header = callbacks.apply_metadata_to_header(
            reader.header,
            sample_header_line=sample_header_line,
            simple_header_lines=simple_header_lines,
            verbose=verbose,
        )
        header = ensure_standard_info_definitions(header, verbose=verbose)
    except SystemExit:
        reader.close()
        raise
    except Exception as exc:  # pragma: no cover - propagated to CLI
        reader.close()
        callbacks.handle_critical_error(f"Failed to apply metadata to merged VCF: {exc}")

    tmp_out = base_vcf + ".tmp"
    try:
        writer = vcfpy.Writer.from_path(tmp_out, header)
    except Exception as exc:  # pragma: no cover - propagated to CLI
        reader.close()
        callbacks.handle_critical_error(f"Failed to open temporary writer for merged VCF: {exc}")

    try:
        for record in reader:
            writer.write_record(record)
    finally:
        reader.close()
        writer.close()

    shutil.move(tmp_out, base_vcf)

    callbacks.log_message("Recomputing AC, AN, AF with bcftools +fill-tags...", verbose)
    res2 = subprocess.run(
        ["bcftools", "+fill-tags", base_vcf, "-Ov", "-o", filled_vcf, "--", "-t", "AC,AN,AF"],
        capture_output=True,
        text=True,
    )
    if res2.returncode != 0:
        callbacks.handle_critical_error((res2.stderr or "bcftools +fill-tags failed").strip())
    shutil.move(filled_vcf, base_vcf)

    callbacks.log_message("Compressing and indexing the final VCF...", verbose)
    try:
        subprocess.run(["bgzip", "-f", base_vcf], check=True)
        subprocess.run(["tabix", "-p", "vcf", "-f", gz_vcf], check=True)
    except subprocess.CalledProcessError as exc:  # pragma: no cover - propagated to CLI
        callbacks.handle_critical_error(
            f"Failed to compress or index merged VCF ({base_vcf}): {exc}"
        )

    callbacks.log_message(
        f"Merged VCF file created and indexed successfully: {gz_vcf}",
        verbose,
    )
    return gz_vcf


def append_metadata_to_merged_vcf(
    merged_vcf,
    sample_metadata_entries=None,
    header_metadata_lines=None,
    serialized_sample_line=None,
    verbose=False,
):
    """Finalize the merged VCF by applying bcftools processing and metadata."""

    callbacks = _require_callbacks()

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
            except Exception:  # pragma: no cover - cleanup best effort
                pass

    callbacks.log_message(
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
    except (subprocess.CalledProcessError, OSError) as exc:  # pragma: no cover - propagated to CLI
        _cleanup_temp_files()
        callbacks.handle_critical_error(
            f"Failed to recalculate INFO tags with bcftools +fill-tags: {exc}"
        )

    callbacks.log_message("Applying quality filters to remove low-confidence variants...", verbose)
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
    except (subprocess.CalledProcessError, OSError) as exc:  # pragma: no cover - propagated to CLI
        _cleanup_temp_files()
        callbacks.handle_critical_error(
            f"Failed to apply bcftools filtering to merged VCF: {exc}"
        )

    callbacks.log_message("Removing individual sample genotype columns (anonymizing data)...", verbose)
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
    except (subprocess.CalledProcessError, OSError) as exc:  # pragma: no cover - propagated to CLI
        _cleanup_temp_files()
        callbacks.handle_critical_error(
            f"Failed to anonymize merged VCF columns using bcftools: {exc}"
        )

    callbacks.log_message("Combining custom metadata header with VCF header...", verbose)
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

        for custom_line in header_metadata_lines:
            custom_line = custom_line.strip()
            if not custom_line:
                continue
            if not custom_line.startswith("##"):
                custom_line = f"##{custom_line}"
            final_header_lines.append(custom_line)

        if serialized_sample_line:
            if not serialized_sample_line.startswith("##"):
                serialized_sample_line = f"##{serialized_sample_line}"
            final_header_lines.append(serialized_sample_line)

        final_header_lines.extend(remaining_header_lines)
        existing_header_lines = set(final_header_lines)

        ensure_standard_info_header_lines(
            final_header_lines, existing_header_lines, verbose=verbose
        )

        with open(final_plain_vcf, "w", encoding="utf-8") as final_handle:
            for line in final_header_lines:
                final_handle.write(f"{line}\n")

            with open(body_temp, "r", encoding="utf-8") as body_handle:
                shutil.copyfileobj(body_handle, final_handle)
    except Exception as exc:  # pragma: no cover - propagated to CLI
        _cleanup_temp_files()
        callbacks.handle_critical_error(
            f"Failed to combine metadata headers for merged VCF: {exc}"
        )

    _cleanup_temp_files()

    callbacks.log_message("Compressing and indexing anonymized VCF...", verbose)
    try:
        if expects_gzip:
            subprocess.run(["bgzip", "-f", final_plain_vcf], check=True)
            subprocess.run(["tabix", "-p", "vcf", "-f", final_vcf], check=True)
        else:
            if final_vcf.endswith(".gz"):
                subprocess.run(["bgzip", "-f", final_plain_vcf], check=True)
            else:
                final_vcf = final_plain_vcf
    except (subprocess.CalledProcessError, OSError) as exc:  # pragma: no cover - propagated
        callbacks.handle_critical_error(
            f"Failed to compress or index anonymized VCF ({final_plain_vcf}): {exc}"
        )

    callbacks.log_message(
        f"Final anonymized VCF created successfully: {final_vcf}",
        verbose,
    )

    return final_vcf
