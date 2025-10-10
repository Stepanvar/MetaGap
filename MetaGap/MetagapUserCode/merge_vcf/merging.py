"""Routines for combining VCF shards into a merged cohort file."""

from __future__ import annotations

import copy
import datetime
import heapq
import itertools
import gzip
import os
import re
import shutil
import subprocess
from collections import OrderedDict
from typing import Callable, List, Optional, Sequence, Tuple, Iterable

from . import vcfpy
from .logging_utils import handle_critical_error, log_message
from .metadata import (
    _parse_sample_metadata_line,
    apply_metadata_to_header,
    build_sample_metadata_line,
)


def _alt_value(alt) -> str:
    """Return a comparable representation for an ALT allele object."""

    if alt is None:
        return ""
    return getattr(alt, "value", str(alt))


def _normalized_alt_key(record) -> Tuple[str, ...]:
    """Return a sorted tuple of ALT allele values for merge-key comparisons."""

    alts = record.ALT or []
    values = {_alt_value(alt) for alt in alts if alt is not None}
    if not values:
        return ()
    return tuple(sorted(values))


def _record_sort_key(record, contig_ranks: dict) -> Tuple[int, int, str, Tuple[str, ...]]:
    """Return a tuple suitable for ordering VCF records during merge."""

    chrom = getattr(record, "CHROM", "")
    rank = contig_ranks.get(chrom, len(contig_ranks))
    pos = getattr(record, "POS", 0) or 0
    ref = getattr(record, "REF", "") or ""
    return (rank, pos, ref, _normalized_alt_key(record))


def _remap_genotype(value: Optional[str], allele_map: dict) -> Optional[str]:
    """Return ``value`` with allele indices remapped according to ``allele_map``."""

    if value is None:
        return None

    text = str(value)
    if not text or text in {".", "./.", ".|."}:
        return text

    parts = re.split(r"([/|])", text)
    remapped: List[str] = []
    for token in parts:
        if token in {"/", "|"}:
            remapped.append(token)
            continue
        if not token:
            continue
        if token == ".":
            remapped.append(token)
            continue
        try:
            allele_index = int(token)
        except ValueError:
            remapped.append(token)
            continue
        if allele_index == 0:
            remapped.append("0")
            continue
        remapped.append(str(allele_map.get(allele_index, allele_index)))
    return "".join(remapped)


def _merge_colliding_records(
    grouped_records: Sequence[Tuple["vcfpy.Record", int]],
    header,
    sample_order: Sequence[str],
) -> "vcfpy.Record":
    """Combine records from distinct shards that reference the same locus."""

    base_record = copy.deepcopy(grouped_records[0][0])

    # Merge ID values while preserving the first occurrence order.
    merged_ids: List[str] = []
    for record, _ in grouped_records:
        raw_id = getattr(record, "ID", None)
        if not raw_id or raw_id == ".":
            continue
        for token in str(raw_id).split(";"):
            token = token.strip()
            if token and token not in merged_ids:
                merged_ids.append(token)
    base_record.ID = ";".join(merged_ids) if merged_ids else "."

    # Combine FILTER entries.
    merged_filters: List[str] = []
    for record, _ in grouped_records:
        filters = getattr(record, "FILTER", None) or []
        for entry in filters:
            if entry not in merged_filters:
                merged_filters.append(entry)
    if not merged_filters:
        merged_filters = []
    base_record.FILTER = merged_filters

    # Construct the unified ALT allele order.
    alt_objects: OrderedDict[str, object] = OrderedDict()
    for record, _ in grouped_records:
        for alt in getattr(record, "ALT", []) or []:
            value = _alt_value(alt)
            if value not in alt_objects:
                alt_objects[value] = copy.deepcopy(alt)
    base_record.ALT = list(alt_objects.values())

    allele_map = {value: index for index, value in enumerate(alt_objects.keys(), start=1)}

    # Consolidate FORMAT keys across records.
    format_keys: List[str] = []
    for record, _ in grouped_records:
        for key in record.FORMAT or []:
            if key not in format_keys:
                format_keys.append(key)
    base_record.FORMAT = format_keys

    # Merge per-sample genotype calls.
    sample_calls: dict = {}
    for record, _ in grouped_records:
        per_record_map = {}
        for idx, alt in enumerate(record.ALT or [], start=1):
            per_record_map[idx] = allele_map[_alt_value(alt)]
        for call in getattr(record, "calls", []) or []:
            data = copy.deepcopy(call.data) if call.data is not None else {}
            if "GT" in data:
                data["GT"] = _remap_genotype(data["GT"], per_record_map)
            sample_calls[call.sample] = vcfpy.Call(call.sample, data)

    ordered_calls: List["vcfpy.Call"] = []
    for name in sample_order:
        call = sample_calls.get(name)
        if call is None:
            call = vcfpy.Call(name, {})
        ordered_calls.append(call)

    base_record.update_calls(ordered_calls)
    _pad_record_samples(base_record, header, sample_order)

    return base_record


def preprocess_vcf(file_path: str, *, chunk_size: int = 1024) -> str:
    """Normalize whitespace delimiters in ``file_path`` if necessary."""

    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")

    opener = gzip.open if str(file_path).endswith(".gz") else open
    mode = "rt" if opener is gzip.open else "r"

    modified = False
    header_found = False

    with opener(file_path, mode, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header_found = True
                normalized = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
                if normalized != line:
                    modified = True
            elif header_found:
                normalized = re.sub(r"\s+", "\t", line.rstrip()) + "\n"
                if normalized != line:
                    modified = True

    if not modified:
        return file_path

    temp_file = f"{file_path}.tmp"

    with opener(file_path, mode, encoding="utf-8") as read_handle, open(
        temp_file, "w", encoding="utf-8"
    ) as write_handle:
        header_found = False
        buffer: List[str] = []

        def flush_buffer() -> None:
            if buffer:
                write_handle.writelines(buffer)
                buffer.clear()

        for line in read_handle:
            if line.startswith("##"):
                buffer.append(line)
            elif line.startswith("#"):
                header_found = True
                buffer.append(re.sub(r"\s+", "\t", line.rstrip()) + "\n")
            else:
                if header_found:
                    buffer.append(re.sub(r"\s+", "\t", line.rstrip()) + "\n")
                else:
                    buffer.append(line)

            if len(buffer) >= chunk_size:
                flush_buffer()

        flush_buffer()

    return temp_file


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


def _ensure_info_header_lines(header) -> None:
    """Ensure standard AC/AN/AF INFO definitions exist on *header*."""

    info_definitions = (
        (
            "AC",
            OrderedDict(
                (
                    ("ID", "AC"),
                    ("Number", "A"),
                    ("Type", "Integer"),
                    (
                        "Description",
                        "Alternate allele count in genotypes, for each ALT allele",
                    ),
                )
            ),
        ),
        (
            "AN",
            OrderedDict(
                (
                    ("ID", "AN"),
                    ("Number", "1"),
                    ("Type", "Integer"),
                    ("Description", "Total number of alleles in called genotypes"),
                )
            ),
        ),
        (
            "AF",
            OrderedDict(
                (
                    ("ID", "AF"),
                    ("Number", "A"),
                    ("Type", "Float"),
                    (
                        "Description",
                        "Allele frequency, for each ALT allele",
                    ),
                )
            ),
        ),
    )

    for key, mapping in info_definitions:
        try:
            existing = header.get_info_field_info(key)
        except Exception:
            existing = None
        if existing is None:
            header.add_info_line(mapping)


def _iter_called_genotype_alleles(value) -> Iterable[int]:
    """Yield allele indices from a genotype value, skipping missing alleles."""

    if value is None:
        return

    if isinstance(value, (list, tuple)):
        for item in value:
            yield from _iter_called_genotype_alleles(item)
        return

    if isinstance(value, int):
        yield value
        return

    if isinstance(value, str):
        stripped = value.strip()
        if not stripped or stripped == ".":
            return

        if "/" in stripped or "|" in stripped:
            parts = re.split(r"[\\/|]", stripped)
        else:
            parts = [stripped]

        for part in parts:
            token = part.strip()
            if not token or token == ".":
                continue
            try:
                yield int(token)
            except ValueError:
                continue
        return

    try:
        yield int(value)
    except (TypeError, ValueError):
        return


def _recompute_ac_an_af(record) -> None:
    """Recompute AC/AN/AF INFO fields for a record based on genotype calls."""

    alt_allele_count = len(getattr(record, "ALT", []) or [])
    ac = [0] * alt_allele_count
    an = 0

    calls = getattr(record, "calls", None) or []
    for call in calls:
        data = getattr(call, "data", {})
        genotype_value = data.get("GT")
        for allele_index in _iter_called_genotype_alleles(genotype_value):
            if isinstance(allele_index, bool):
                continue
            an += 1
            if allele_index is None or allele_index <= 0:
                continue
            if 1 <= allele_index <= alt_allele_count:
                ac[allele_index - 1] += 1

    record.INFO["AC"] = ac
    record.INFO["AN"] = an
    if an > 0:
        record.INFO["AF"] = [count / an for count in ac]
    else:
        record.INFO["AF"] = [0.0] * alt_allele_count


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

    if not preprocessed_files:
        handle_critical_error("No VCF files available for merging.")
    try:
        unified_header_template = union_headers(valid_files, sample_order=sample_order)
    except SystemExit:
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to construct unified header: {exc}")

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

    log_message("Constructing unified header for streaming merge...", verbose)
    try:
        combined_header = union_headers(preprocessed_files, sample_order=sample_order)
    except SystemExit:
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to construct unified VCF header: {exc}")

    # Ensure combined header includes all samples observed across inputs.
    sample_names: List[str] = []
    if hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
        sample_names = list(combined_header.samples.names or [])

    readers = []
    reader_iters = []
    format_ids = set()
    for path in preprocessed_files:
        try:
            reader = vcfpy.Reader.from_path(path)
        except Exception as exc:
            handle_critical_error(f"Failed to open VCF for merging ({path}): {exc}")
        readers.append(reader)
        reader_iters.append(iter(reader))

        file_samples = []
        if hasattr(reader.header, "samples") and hasattr(reader.header.samples, "names"):
            file_samples = list(reader.header.samples.names or [])
        for name in file_samples:
            if name not in sample_names:
                sample_names.append(name)

        for line in getattr(reader.header, "lines", []):
            if isinstance(line, vcfpy.header.FormatHeaderLine) and line.id not in format_ids:
                combined_header.add_line(copy.deepcopy(line))
                format_ids.add(line.id)

        header_formats = getattr(reader.header, "formats", None)
        if header_formats is not None:
            for fmt_id, fmt_line in header_formats.items():
                if fmt_id in format_ids:
                    continue
                combined_formats = getattr(combined_header, "formats", None)
                if combined_formats is not None:
                    try:
                        combined_formats[fmt_id] = copy.deepcopy(fmt_line)
                    except Exception:
                        pass
                format_ids.add(fmt_id)

    if hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
        combined_header.samples.names = list(sample_names)

    header = apply_metadata_to_header(
        combined_header,
        sample_header_line=sample_header_line,
        simple_header_lines=simple_header_lines,
        verbose=verbose,
    )

    if hasattr(header, "samples") and hasattr(header.samples, "names"):
        sample_names = list(header.samples.names or sample_names)

    contig_names: List[str] = []
    if hasattr(header, "contigs") and header.contigs:
        contig_names = list(header.contigs.keys())
    if not contig_names:
        for line in getattr(header, "lines", []):
            if isinstance(line, vcfpy.header.ContigHeaderLine):
                contig_names.append(line.id)
    contig_ranks = {name: idx for idx, name in enumerate(contig_names)}

    tmp_out = base_vcf + ".tmp"
    _ensure_info_header_lines(header)

    try:
        writer = vcfpy.Writer.from_path(base_vcf, header)
    except Exception as exc:
        for reader in readers:
            try:
                reader.close()
            except Exception:
                pass
        handle_critical_error(f"Failed to open writer for merged VCF: {exc}")

    def _advance(idx: int):
        iterator = reader_iters[idx]
        try:
            return next(iterator)
        except StopIteration:
            return None

    log_message("Performing heap-based k-way merge of VCF shards...", verbose)
    heap: List[Tuple[Tuple[int, int, str, Tuple[str, ...]], int, int, "vcfpy.Record"]] = []
    counter = itertools.count()
    for idx in range(len(reader_iters)):
        record = _advance(idx)
        if record is None:
            continue
        key = _record_sort_key(record, contig_ranks)
        heapq.heappush(heap, (key, next(counter), idx, record))

    target_sample_order = []
    if hasattr(header, "samples") and hasattr(header.samples, "names"):
        target_sample_order = list(header.samples.names)

    try:
        while heap:
            key, _, src_idx, record = heapq.heappop(heap)
            colliding: List[Tuple["vcfpy.Record", int]] = [(record, src_idx)]

            while heap and heap[0][0] == key:
                _, _, other_idx, other_record = heapq.heappop(heap)
                colliding.append((other_record, other_idx))

            pending: List[Tuple["vcfpy.Record", int]] = []
            idx = 0
            while idx < len(colliding):
                current_record, reader_index = colliding[idx]
                next_record = _advance(reader_index)
                while next_record is not None and _record_sort_key(next_record, contig_ranks) == key:
                    colliding.append((next_record, reader_index))
                    next_record = _advance(reader_index)
                if next_record is not None:
                    pending.append((next_record, reader_index))
                idx += 1

            merged_record = _merge_colliding_records(colliding, header, sample_names)
            _recompute_ac_an_af(merged_record)
            writer.write_record(merged_record)

            for next_record, reader_index in pending:
                heapq.heappush(
                    heap,
                    (
                        _record_sort_key(next_record, contig_ranks),
                        next(counter),
                        reader_index,
                        next_record,
                    ),
                )
    finally:
        writer.close()
        for reader in readers:
            try:
                reader.close()
            except Exception:
                pass

    for tmp in temp_files:
        try:
            if os.path.exists(tmp):
                os.remove(tmp)
        except OSError:
            pass

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
    info_lines: "OrderedDict[str, vcfpy.header.InfoHeaderLine]" = OrderedDict()
    filter_lines: "OrderedDict[str, vcfpy.header.FilterHeaderLine]" = OrderedDict()
    contig_lines: "OrderedDict[str, vcfpy.header.ContigHeaderLine]" = OrderedDict()
    format_lines: "OrderedDict[str, vcfpy.header.FormatHeaderLine]" = OrderedDict()
    computed_sample_order: List[str] = []
    merged_sample_metadata = None

    def _line_mapping(line) -> "OrderedDict[str, str]":
        mapping = getattr(line, "mapping", None)
        if isinstance(mapping, dict):
            return OrderedDict(mapping)
        return OrderedDict()

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
                    computed_sample_order = list(combined_header.samples.names)
                for line in combined_header.lines:
                    if isinstance(line, vcfpy.header.InfoHeaderLine):
                        info_lines[line.id] = copy.deepcopy(line)
                    elif isinstance(line, vcfpy.header.FilterHeaderLine):
                        filter_lines[line.id] = copy.deepcopy(line)
                    elif isinstance(line, vcfpy.header.ContigHeaderLine):
                        contig_lines[line.id] = copy.deepcopy(line)
                    elif isinstance(line, vcfpy.header.FormatHeaderLine):
                        format_lines[line.id] = copy.deepcopy(line)
                    elif isinstance(line, vcfpy.header.SampleHeaderLine):
                        mapping = _line_mapping(line)
                        if mapping.get("ID"):
                            merged_sample_metadata = OrderedDict(mapping)
                continue

            if hasattr(header, "samples") and hasattr(header.samples, "names"):
                for sample_name in header.samples.names:
                    if sample_name not in computed_sample_order:
                        computed_sample_order.append(sample_name)

            for line in header.lines:
                if isinstance(line, vcfpy.header.InfoHeaderLine):
                    existing = info_lines.get(line.id)
                    if existing is None:
                        info_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    existing_mapping = _line_mapping(existing)
                    new_mapping = _line_mapping(line)
                    for key in ("Number", "Type", "Description"):
                        if existing_mapping.get(key) != new_mapping.get(key):
                            handle_critical_error(
                                "INFO header definitions conflict across shards. "
                                f"Field '{key}' for INFO '{line.id}' differs between shards: "
                                f"{existing_mapping.get(key)!r} vs {new_mapping.get(key)!r} "
                                f"(encountered in {file_path})."
                            )
                elif isinstance(line, vcfpy.header.FilterHeaderLine):
                    existing = filter_lines.get(line.id)
                    if existing is None:
                        filter_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    existing_mapping = _line_mapping(existing)
                    new_mapping = _line_mapping(line)
                    if existing_mapping.get("Description") != new_mapping.get("Description"):
                        handle_critical_error(
                            "FILTER header definitions conflict across shards. "
                            f"Description for FILTER '{line.id}' differs between shards: "
                            f"{existing_mapping.get('Description')!r} vs {new_mapping.get('Description')!r} "
                            f"(encountered in {file_path})."
                        )
                elif isinstance(line, vcfpy.header.ContigHeaderLine):
                    existing = contig_lines.get(line.id)
                    if existing is None:
                        contig_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    if _line_mapping(existing) != _line_mapping(line):
                        handle_critical_error(
                            "Contig header definitions conflict across shards. "
                            f"Contig '{line.id}' differs between shards (conflict in {file_path})."
                        )
                elif isinstance(line, vcfpy.header.FormatHeaderLine):
                    existing = format_lines.get(line.id)
                    if existing is None:
                        format_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    existing_mapping = _line_mapping(existing)
                    new_mapping = _line_mapping(line)
                    for key in ("Number", "Type", "Description"):
                        if existing_mapping.get(key) != new_mapping.get(key):
                            handle_critical_error(
                                "FORMAT header definitions conflict across shards. "
                                f"Field '{key}' for FORMAT '{line.id}' differs between shards: "
                                f"{existing_mapping.get(key)!r} vs {new_mapping.get(key)!r} "
                                f"(encountered in {file_path})."
                            )
                elif isinstance(line, vcfpy.header.SampleHeaderLine):
                    mapping = _line_mapping(line)
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
