from __future__ import annotations

import copy
import datetime
import heapq
import itertools
import gzip
import os
import re
import shutil
from collections import OrderedDict
from typing import Callable, List, Optional, Sequence, Tuple, Iterable

import pysam  # BGZF + Tabix
import concurrent.futures
from . import vcfpy
from .logging_utils import (
    MergeConflictError,
    handle_critical_error,
    log_message,
)
from .metadata import (
    _parse_sample_metadata_line,
    apply_metadata_to_header,
    build_sample_metadata_line,
)

def _alt_value(alt) -> str:
    if alt is None:
        return ""
    return getattr(alt, "value", str(alt))

def _normalized_alt_key(record) -> Tuple[str, ...]:
    alts = record.ALT or []
    values = {_alt_value(alt) for alt in alts if alt is not None}
    return tuple(sorted(values)) if values else ()

def _record_sort_key(record, contig_ranks: dict) -> Tuple[int, int, str, Tuple[str, ...]]:
    chrom = getattr(record, "CHROM", "")
    rank = contig_ranks.get(chrom, len(contig_ranks))
    pos = getattr(record, "POS", 0) or 0
    ref = getattr(record, "REF", "") or ""
    return (rank, pos, ref, _normalized_alt_key(record))

def _sort_record_alts(record) -> None:
    """Sort the ALT alleles of *record* lexicographically in-place."""
    if record is None:
        return
    alts = getattr(record, "ALT", None) or []
    if len(alts) <= 1:
        return
    sorted_alts = sorted(alts, key=_alt_value)
    if list(alts) != sorted_alts:
        record.ALT = list(sorted_alts)

def _apply_contig_order(header, contig_order: Sequence[str]) -> None:
    """Ensure *header* preserves *contig_order* for contig definitions."""
    if header is None or not contig_order:
        return
    contig_lines = [
        line
        for line in getattr(header, "lines", [])
        if isinstance(line, vcfpy.header.ContigHeaderLine)
    ]
    if not contig_lines:
        return
    contig_lookup: "OrderedDict[str, object]" = OrderedDict()
    for line in contig_lines:
        identifier = getattr(line, "id", None)
        if identifier is None:
            continue
        if identifier not in contig_lookup:
            contig_lookup[identifier] = line
    ordered_ids: List[str] = []
    for name in contig_order:
        if name in contig_lookup:
            ordered_ids.append(name)
    for name in contig_lookup:
        if name not in ordered_ids:
            ordered_ids.append(name)
    ordered_lines = [contig_lookup[name] for name in ordered_ids]
    iterator = iter(ordered_lines)
    new_lines: List[object] = []
    for line in getattr(header, "lines", []):
        if isinstance(line, vcfpy.header.ContigHeaderLine):
            new_lines.append(next(iterator))
        else:
            new_lines.append(line)
    header.lines = new_lines
    contigs_attr = getattr(header, "contigs", None)
    if contigs_attr:
        reordered = OrderedDict()
        for name in ordered_ids:
            if name in contigs_attr:
                reordered[name] = contigs_attr[name]
        for name, value in contigs_attr.items():
            if name not in reordered:
                reordered[name] = value
        try:
            header.contigs = reordered
        except Exception:
            pass
    setattr(header, "_metagap_contig_order", list(ordered_ids))

def _remap_genotype(value: Optional[str], allele_map: dict) -> Optional[str]:
    if value is None:
        return None
    text = str(value)
    if not text or text in {".", "./.", ".|."}:
        return text
    parts = re.split(r"([/|])", text)
    out: List[str] = []
    for tok in parts:
        if tok in {"/", "|"}:
            out.append(tok)
            continue
        if tok == "" or tok == ".":
            out.append(tok or "")
            continue
        try:
            a = int(tok)
        except ValueError:
            out.append(tok)
            continue
        if a == 0:
            out.append("0")
        else:
            out.append(str(allele_map.get(a, a)))
    return "".join(out)

def _merge_colliding_records(
    grouped_records: Sequence[Tuple["vcfpy.Record", int]],
    header,
    sample_order: Sequence[str],
) -> "vcfpy.Record":
    base = copy.deepcopy(grouped_records[0][0])
    # IDs: merge unique IDs from all records
    merged_ids: List[str] = []
    for r, _ in grouped_records:
        rid = getattr(r, "ID", None)
        if not rid or rid == ".":
            continue
        for tok in str(rid).split(";"):
            tok = tok.strip()
            if tok and tok not in merged_ids:
                merged_ids.append(tok)
    base.ID = ";".join(merged_ids) if merged_ids else "."
    # FILTER: merge unique FILTERs from all records
    merged_filters: List[str] = []
    for r, _ in grouped_records:
        for f in (r.FILTER or []):
            if f not in merged_filters:
                merged_filters.append(f)
    base.FILTER = merged_filters
    # Unified ALT order across records
    alt_objects: OrderedDict[str, object] = OrderedDict()
    for r, _ in grouped_records:
        for alt in (r.ALT or []):
            val = _alt_value(alt)
            if val not in alt_objects:
                alt_objects[val] = copy.deepcopy(alt)
    base.ALT = list(alt_objects.values())
    allele_map_global = {val: i for i, val in enumerate(alt_objects.keys(), start=1)}
    # Union of all FORMAT keys
    fmt_keys: List[str] = []
    for r, _ in grouped_records:
        for k in (r.FORMAT or []):
            if k not in fmt_keys:
                fmt_keys.append(k)
    base.FORMAT = fmt_keys
    # Merge calls for all samples, with genotype allele remapping
    calls_by_sample: dict = {}
    for r, _ in grouped_records:
        per_map = {}
        for idx, alt in enumerate(r.ALT or [], start=1):
            per_map[idx] = allele_map_global[_alt_value(alt)]
        for call in (getattr(r, "calls", []) or []):
            data = copy.deepcopy(call.data) if call.data is not None else {}
            if "GT" in data:
                data["GT"] = _remap_genotype(data["GT"], per_map)
            calls_by_sample[call.sample] = vcfpy.Call(call.sample, data)
    # Order calls by final sample list, adding empty calls for missing samples
    ordered_calls: List[vcfpy.Call] = []
    for name in sample_order:
        call = calls_by_sample.get(name)
        if call is None:
            call = vcfpy.Call(name, {})
        ordered_calls.append(call)
    base.update_calls(ordered_calls)
    _pad_record_samples(base, header, sample_order)
    return base

def preprocess_vcf(file_path: str, *, chunk_size: int = 1024) -> str:
    """Normalize whitespace delimiters in `file_path` if necessary.

    Tabs are mandatory between VCF columns; meta-info lines (##...) must be preserved as-is.
    Returns original path if no changes were needed; otherwise writes a normalized .tmp file and returns its path.
    """
    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")
    opener = gzip.open if str(file_path).endswith(".gz") else open
    mode = "rt" if opener is gzip.open else "r"
    # First pass: detect if normalization is needed (do not buffer whole file)
    modified = False
    header_found = False
    with opener(file_path, mode, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header_found = True
                if re.sub(r"\s+", "\t", line.rstrip()) + "\n" != line:
                    modified = True
            elif header_found:
                if re.sub(r"\s+", "\t", line.rstrip()) + "\n" != line:
                    modified = True
            # keep scanning if header not reached yet
    if not modified:
        return file_path
    # Second pass: rewrite file with normalized tabs in a buffer
    temp_file = f"{file_path}.tmp"
    with opener(file_path, mode, encoding="utf-8") as read_handle, open(temp_file, "w", encoding="utf-8") as write_handle:
        header_found = False
        buffer: List[str] = []
        def flush_buffer():
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
    template: dict = {}
    for key in format_keys:
        if key == "GT":
            template[key] = "./."
            continue
        try:
            fi = header.get_format_field_info(key)
        except Exception:
            fi = None
        number = getattr(fi, "number", None)
        if isinstance(number, str):
            s = number.strip()
            if not s:
                number = None
            elif s in {".", "A", "G", "R"}:
                template[key] = []
                continue
            else:
                try:
                    number = int(s)
                except ValueError:
                    number = None
        if isinstance(number, int):
            default_value = None if number <= 1 else [None] * number
        elif number in {".", "A", "G", "R"}:
            default_value = []
        else:
            default_value = None
        template.setdefault(key, default_value)
    def factory() -> dict:
        data: dict = {}
        for k, d in template.items():
            data[k] = ([None] * len(d)) if isinstance(d, list) else d
        return data
    return factory

def _ensure_info_header_lines(header) -> None:
    # Ensure standard AC/AN/AF INFO definitions present
    info_defs = (
        ("AC", OrderedDict((("ID", "AC"), ("Number", "A"), ("Type", "Integer"),
            ("Description", "Alternate allele count in genotypes, for each ALT allele")))),
        ("AN", OrderedDict((("ID", "AN"), ("Number", "1"), ("Type", "Integer"),
            ("Description", "Total number of alleles in called genotypes")))),
        ("AF", OrderedDict((("ID", "AF"), ("Number", "A"), ("Type", "Float"),
            ("Description", "Allele frequency, for each ALT allele")))),
    )
    for key, mapping in info_defs:
        try:
            existing = header.get_info_field_info(key)
        except Exception:
            existing = None
        if existing is None:
            header.add_info_line(mapping)

def _iter_called_genotype_alleles(value) -> Iterable[int]:
    if value is None:
        return
    if isinstance(value, (list, tuple)):
        for item in value:
            yield from _iter_called_genotype_alleles(item)
        return
    if isinstance(value, int):
        yield value; return
    if isinstance(value, str):
        s = value.strip()
        if not s or s == ".":
            return
        parts = re.split(r"[\\/|]", s) if ("/" in s or "|" in s) else [s]
        for part in parts:
            t = part.strip()
            if not t or t == ".":
                continue
            try:
                yield int(t)
            except ValueError:
                continue
        return
    try:
        yield int(value)
    except (TypeError, ValueError):
        return

def _recompute_ac_an_af(record) -> None:
    alt_n = len(record.ALT or [])
    ac = [0] * alt_n
    an = 0
    for call in (record.calls or []):
        g = (call.data or {}).get("GT")
        for ai in _iter_called_genotype_alleles(g):
            if isinstance(ai, bool):
                continue
            an += 1
            if ai is None or ai <= 0:
                continue
            if 1 <= ai <= alt_n:
                ac[ai - 1] += 1
    record.INFO["AC"] = ac
    record.INFO["AN"] = an
    record.INFO["AF"] = ([c / an for c in ac] if an > 0 else [0.0] * alt_n)

def _remove_format_and_sample_definitions(header) -> None:
    if header is None:
        return
    if hasattr(header, "lines"):
        filtered = []
        for line in header.lines:
            if isinstance(line, vcfpy.header.FormatHeaderLine):
                continue
            key = getattr(line, "key", None)
            if isinstance(key, str) and key.upper() == "FORMAT":
                continue
            filtered.append(line)
        header.lines = filtered
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
    fmt_keys = list(record.FORMAT or [])
    factory = _create_missing_call_factory(fmt_keys, header)
    new_calls = []
    for name in sample_order:
        call = record.call_for_sample.get(name)
        if call is None:
            call = vcfpy.Call(name, factory())
        else:
            defaults = factory()
            for k in fmt_keys:
                if k not in call.data:
                    call.data[k] = defaults.get(k)
                else:
                    if isinstance(call.data[k], tuple):
                        call.data[k] = list(call.data[k])
        if "GT" in call.data and call.data["GT"] in {None, "", "."}:
            call.data["GT"] = defaults.get("GT", "./.")
        new_calls.append(call)
    record.update_calls(new_calls)

def _record_passes_filters(
    record,
    qual_threshold: Optional[float],
    an_threshold: Optional[float],
    allowed_filter_values: Optional[Sequence[str]],
) -> bool:
    if qual_threshold is not None:
        q = getattr(record, "QUAL", None)
        try:
            if q is None or float(q) <= qual_threshold:
                return False
        except (TypeError, ValueError):
            return False
    if an_threshold is not None:
        an_val = record.INFO.get("AN")
        if isinstance(an_val, list):
            an_val = an_val[0] if an_val else None
        try:
            if an_val is None or float(an_val) <= an_threshold:
                return False
        except (TypeError, ValueError):
            return False
    filters = record.FILTER or []
    if not filters:
        return True
    vals = [f for f in filters if f not in {None, "", "."}]
    if not vals:
        return True
    allowed = tuple(allowed_filter_values or ("PASS",))
    allowed_set = {v for v in allowed if v not in {None, "", "."}}
    return bool(allowed_set) and all(v in allowed_set for v in vals)

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
    tmp = input_path + ".filtered"
    kept = 0
    total = 0
    try:
        writer = vcfpy.Writer.from_path(tmp, reader.header)
    except Exception as exc:
        reader.close()
        handle_critical_error(f"Failed to open filtered VCF writer ({input_path}): {exc}")
    try:
        for rec in reader:
            total += 1
            if _record_passes_filters(rec, qual_threshold, an_threshold, allowed_filter_values):
                writer.write_record(rec)
                kept += 1
    finally:
        reader.close()
        writer.close()
    shutil.move(tmp, input_path)
    log_message(f"Applied variant filter: kept {kept} of {total}.", verbose)

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
    log_message(f"Discovered {file_count} validated VCF shard(s) for merging.", verbose)
    temp_files: List[str] = []
    preprocessed: List[str] = []
    # Preprocess all input files concurrently for performance
    try:
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(preprocess_vcf, valid_files))
    except Exception as exc:
        handle_critical_error(f"Failed to preprocess input VCF files: {exc}", exc_cls=MergeConflictError)
    for orig_path, pre_path in zip(valid_files, results):
        preprocessed.append(pre_path)
        if pre_path != orig_path:
            temp_files.append(pre_path)
    # Merge headers from all files to create a combined header
    try:
        merged_header = union_headers(preprocessed, sample_order=sample_order)
    except MergeConflictError:
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to merge VCF headers: {exc}", exc_cls=MergeConflictError)
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
    # Determine contig ordering and finalize header
    contigs: List[str] = []
    if getattr(header, "contigs", None):
        contigs = list(header.contigs.keys())
    if not contigs:
        for ln in getattr(header, "lines", []):
            if isinstance(ln, vcfpy.header.ContigHeaderLine):
                contigs.append(ln.id)
    contig_ranks = {name: i for i, name in enumerate(contigs)}
    _ensure_info_header_lines(header)
    # Ensure header.formats includes all FORMAT definitions
    if hasattr(header, "formats") and header.formats is not None:
        for line in getattr(header, "lines", []):
            if isinstance(line, vcfpy.header.FormatHeaderLine):
                header.formats[line.id] = line
    # Prepare final sample list
    sample_names: List[str] = list(getattr(header.samples, "names", []))
    # Open readers for each VCF shard and perform k-way merge
    readers = []
    iters = []
    for path in preprocessed:
        try:
            r = vcfpy.Reader.from_path(path)
        except Exception as exc:
            handle_critical_error(f"Failed to open VCF for merging ({path}): {exc}")
        readers.append(r)
        iters.append(iter(r))
        # Grow union of samples if any new ones discovered
        file_samps = list(getattr(r.header.samples, "names", []) or [])
        for n in file_samps:
            if n not in sample_names:
                sample_names.append(n)
    if hasattr(header, "samples") and hasattr(header.samples, "names"):
        header.samples.names = list(sample_names)
    # Write merged records
    try:
        writer = vcfpy.Writer.from_path(base_vcf, header)
    except Exception as exc:
        for r in readers:
            try:
                r.close()
            except Exception:
                pass
        handle_critical_error(f"Failed to open writer for merged VCF: {exc}", exc_cls=MergeConflictError)
    def _advance(i: int):
        try:
            rec = next(iters[i])
        except StopIteration:
            return None
        _sort_record_alts(rec)
        return rec
    log_message("Performing heap-based k-way merge of VCF shards...", verbose)
    heap: List[Tuple[Tuple[int, int, str, Tuple[str, ...]], int, int, vcfpy.Record]] = []
    counter = itertools.count()
    for i in range(len(iters)):
        rec = _advance(i)
        if rec is None:
            continue
        heapq.heappush(heap, (_record_sort_key(rec, contig_ranks), next(counter), i, rec))
    try:
        while heap:
            key, _, src_idx, rec = heapq.heappop(heap)
            colliding_records: List[Tuple[vcfpy.Record, int]] = [(rec, src_idx)]
            while heap and heap[0][0] == key:
                _, _, j, rec2 = heapq.heappop(heap)
                colliding_records.append((rec2, j))
            pending: List[Tuple[vcfpy.Record, int]] = []
            for cur_idx, ridx in enumerate(colliding_records):
                cur_rec, src_index = colliding_records[cur_idx]
                nxt = _advance(src_index)
                while nxt is not None and _record_sort_key(nxt, contig_ranks) == key:
                    colliding_records.append((nxt, src_index))
                    nxt = _advance(src_index)
                if nxt is not None:
                    pending.append((nxt, src_index))
            merged_rec = _merge_colliding_records(colliding_records, header, sample_names)
            _recompute_ac_an_af(merged_rec)
            writer.write_record(merged_rec)
            for nxt, ridx in pending:
                heapq.heappush(heap, (_record_sort_key(nxt, contig_ranks), next(counter), ridx, nxt))
    finally:
        writer.close()
        for r in readers:
            try:
                r.close()
            except Exception:
                pass
        for t in temp_files:
            try:
                if os.path.exists(t):
                    os.remove(t)
            except OSError:
                pass
    # Apply in-memory variant filtering thresholds
    _filter_vcf_records(
        base_vcf,
        qual_threshold=qual_threshold,
        an_threshold=an_threshold,
        allowed_filter_values=allowed_filter_values,
        verbose=verbose,
    )
    # Anonymize: drop all FORMAT fields and sample columns for privacy
    try:
        rdr = vcfpy.Reader.from_path(base_vcf)
    except Exception as exc:
        handle_critical_error(f"Failed to reopen VCF for anonymization: {exc}")
    _remove_format_and_sample_definitions(rdr.header)
    anon_tmp = base_vcf + ".anon"
    try:
        w = vcfpy.Writer.from_path(anon_tmp, rdr.header)
    except Exception as exc:
        rdr.close()
        handle_critical_error(f"Failed to open anonymized writer: {exc}")
    try:
        for rec in rdr:
            rec.FORMAT = []
            rec.calls = []
            w.write_record(rec)
    finally:
        rdr.close()
        w.close()
    shutil.move(anon_tmp, base_vcf)
    # Compress and index the final VCF using pysam
    log_message("Compressing and indexing the final VCF...", verbose)
    try:
        pysam.tabix_compress(base_vcf, gz_vcf, force=True)
        pysam.tabix_index(gz_vcf, preset="vcf", force=True)
        os.remove(base_vcf)
    except Exception as exc:
        handle_critical_error(f"Failed to compress/index final VCF: {exc}", exc_cls=MergeConflictError)
    log_message(f"Merged VCF created and indexed successfully: {gz_vcf}", verbose)
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
        return OrderedDict(mapping) if isinstance(mapping, dict) else OrderedDict()
    for file_path in valid_files:
        pre = preprocess_vcf(file_path)
        reader = None
        try:
            reader = vcfpy.Reader.from_path(pre)
            header = reader.header
            if combined_header is None:
                combined_header = header.copy()
                if hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
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
                        m = _line_mapping(line)
                        if m.get("ID"):
                            merged_sample_metadata = OrderedDict(m)
                continue
            if hasattr(header, "samples") and hasattr(header.samples, "names"):
                for s in header.samples.names:
                    if s not in computed_sample_order:
                        computed_sample_order.append(s)
            for line in header.lines:
                if isinstance(line, vcfpy.header.InfoHeaderLine):
                    exist = info_lines.get(line.id)
                    if exist is None:
                        info_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    em, nm = _line_mapping(exist), _line_mapping(line)
                    for k in ("Number", "Type", "Description"):
                        if em.get(k) != nm.get(k):
                            handle_critical_error(
                                "INFO header definitions conflict across shards. "
                                f"Field '{k}' for INFO '{line.id}' differs: {em.get(k)!r} vs {nm.get(k)!r} "
                                f"(in {file_path}).",
                                exc_cls=MergeConflictError,
                            )
                elif isinstance(line, vcfpy.header.FilterHeaderLine):
                    exist = filter_lines.get(line.id)
                    if exist is None:
                        filter_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    em, nm = _line_mapping(exist), _line_mapping(line)
                    if em.get("Description") != nm.get("Description"):
                        handle_critical_error(
                            "FILTER header definitions conflict across shards. "
                            f"Description for FILTER '{line.id}' differs: {em.get('Description')!r} vs {nm.get('Description')!r} "
                            f"(in {file_path}).",
                            exc_cls=MergeConflictError,
                        )
                elif isinstance(line, vcfpy.header.ContigHeaderLine):
                    exist = contig_lines.get(line.id)
                    if exist is None:
                        contig_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    if _line_mapping(exist) != _line_mapping(line):
                        handle_critical_error(
                            "Contig header definitions conflict across shards. "
                            f"Contig '{line.id}' differs between shards (conflict in {file_path}).",
                            exc_cls=MergeConflictError,
                        )
                elif isinstance(line, vcfpy.header.FormatHeaderLine):
                    exist = format_lines.get(line.id)
                    if exist is None:
                        format_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                        continue
                    em, nm = _line_mapping(exist), _line_mapping(line)
                    for k in ("Number", "Type", "Description"):
                        if em.get(k) != nm.get(k):
                            handle_critical_error(
                                "FORMAT header definitions conflict across shards. "
                                f"Field '{k}' for FORMAT '{line.id}' differs: {em.get(k)!r} vs {nm.get(k)!r} "
                                f"(in {file_path}).",
                                exc_cls=MergeConflictError,
                            )
                elif isinstance(line, vcfpy.header.SampleHeaderLine):
                    m = _line_mapping(line)
                    sid = m.get("ID")
                    if not sid:
                        continue
                    if merged_sample_metadata is None:
                        merged_sample_metadata = OrderedDict(m)
                        continue
                    for k, v in m.items():
                        if k not in merged_sample_metadata or not merged_sample_metadata[k]:
                            merged_sample_metadata[k] = v
        except Exception as exc:
            handle_critical_error(f"Failed to read VCF header from {file_path}: {exc}", exc_cls=MergeConflictError)
        finally:
            if reader is not None:
                try:
                    reader.close()
                except Exception:
                    pass
            if pre != file_path and os.path.exists(pre):
                os.remove(pre)
    if combined_header is None:
        handle_critical_error("Unable to construct a merged VCF header.", exc_cls=MergeConflictError)
    target_samples = sample_order if sample_order is not None else computed_sample_order
    if target_samples and hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
        combined_header.samples.names = list(target_samples)
    if merged_sample_metadata:
        serialized = build_sample_metadata_line(merged_sample_metadata)
        parsed = _parse_sample_metadata_line(serialized)
        sample_line = vcfpy.header.SampleHeaderLine.from_mapping(parsed)
        combined_header.lines = [ln for ln in combined_header.lines if getattr(ln, "key", None) != "SAMPLE"]
        combined_header.add_line(sample_line)
    contig_order = list(contig_lines.keys())
    if contig_order:
        _apply_contig_order(combined_header, contig_order)
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
