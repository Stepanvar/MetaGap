from __future__ import annotations

import logging
import copy
import datetime
import heapq
import itertools
import gzip
import os
import re
import shutil
from typing import Callable, Iterable, List, Optional, Sequence, Tuple

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
    """Sorting key for records: by contig rank, position, REF, and ALT alleles."""
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
    # Reorder contig header lines to match the given contig_order
    contig_lines = [
        line for line in getattr(header, "lines", [])
        if isinstance(line, vcfpy.header.ContigHeaderLine)
    ]
    if not contig_lines:
        return
    contig_lookup: dict[str, object] = {}
    for line in contig_lines:
        identifier = getattr(line, "id", None)
        if identifier is None:
            continue
        if identifier not in contig_lookup:
            contig_lookup[identifier] = line
    # Build ordered list of contig IDs: first those in contig_order, then the rest
    ordered_ids: List[str] = []
    for name in contig_order:
        if name in contig_lookup:
            ordered_ids.append(name)
    for name in contig_lookup:
        if name not in ordered_ids:
            ordered_ids.append(name)
    # Replace header lines with ordered contigs
    ordered_lines = [contig_lookup[name] for name in ordered_ids]
    iterator = iter(ordered_lines)
    new_lines: List[object] = []
    for line in getattr(header, "lines", []):
        if isinstance(line, vcfpy.header.ContigHeaderLine):
            new_lines.append(next(iterator))
        else:
            new_lines.append(line)
    header.lines = new_lines
    # Also reorder header.contigs attribute if present
    contigs_attr = getattr(header, "contigs", None)
    if contigs_attr:
        reordered: dict[str, object] = {}
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
    # Store the applied contig order (for reference/debugging)
    setattr(header, "_metagap_contig_order", list(ordered_ids))

def _remap_genotype(value: Optional[str], allele_map: dict) -> Optional[str]:
    """Remap genotype allele indices in a GT string according to allele_map."""
    if value is None:
        return None
    text = str(value)
    if not text or text in {".", "./.", ".|."}:
        return text
    # Split genotype string by separator ("/" or "|"), preserve separators in parts
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
    grouped_records: Sequence[Tuple[vcfpy.Record, int]],
    header,
    sample_order: Sequence[str],
) -> vcfpy.Record:
    """Merge multiple VCF records representing the same variant into one record."""
    # Start with a deep copy of the first record as the base
    base = copy.deepcopy(grouped_records[0][0])
    # Merge IDs: collect unique IDs from all records
    merged_ids: List[str] = []
    for record, _ in grouped_records:
        rid = getattr(record, "ID", None)
        if not rid or rid == ".":
            continue
        for tok in str(rid).split(";"):
            tok = tok.strip()
            if tok and tok not in merged_ids:
                merged_ids.append(tok)
    base.ID = ";".join(merged_ids) if merged_ids else "."
    # Merge FILTERs: combine unique filter flags from all records
    merged_filters: List[str] = []
    for record, _ in grouped_records:
        for f in (record.FILTER or []):
            if f not in merged_filters:
                merged_filters.append(f)
    base.FILTER = merged_filters
    # Unified ALT allele list across all records (ensure consistent indexing for genotypes)
    alt_objects: dict[str, object] = {}
    for record, _ in grouped_records:
        for alt in (record.ALT or []):
            val = _alt_value(alt)
            if val not in alt_objects:
                alt_objects[val] = copy.deepcopy(alt)
    base.ALT = list(alt_objects.values())
    # Mapping from original allele value to new allele index (for genotype remapping)
    allele_map_global = {val: i for i, val in enumerate(alt_objects.keys(), start=1)}
    # Union of all FORMAT keys across records
    fmt_keys: List[str] = []
    for record, _ in grouped_records:
        for key in (record.FORMAT or []):
            if key not in fmt_keys:
                fmt_keys.append(key)
    base.FORMAT = fmt_keys
    # Merge genotype calls for all samples, remapping allele indices in GT fields
    calls_by_sample: dict = {}
    for record, _ in grouped_records:
        # Map allele indices from this record's ALT set to global ALT indices
        per_map: dict = {}
        for idx, alt in enumerate(record.ALT or [], start=1):
            per_map[idx] = allele_map_global[_alt_value(alt)]
        for call in (record.calls or []):
            sample_name = call.sample  # vcfpy.Call.sample holds the sample name
            data = copy.deepcopy(call.data) if call.data is not None else {}
            if "GT" in data:
                data["GT"] = _remap_genotype(data["GT"], per_map)
            calls_by_sample[sample_name] = vcfpy.Call(sample_name, data)
    # Order merged calls by the final sample order, add empty calls for missing samples
    ordered_calls: List[vcfpy.Call] = []
    for name in sample_order:
        call = calls_by_sample.get(name)
        if call is None:
            call = vcfpy.Call(name, {})
        ordered_calls.append(call)
    base.update_calls(ordered_calls)
    # Pad any missing FORMAT fields with default values for all calls
    _pad_record_samples(base, header, sample_order)
    return base

def preprocess_vcf(file_path: str, *, chunk_size: int = 1024) -> str:
    """Normalize whitespace delimiters in `file_path` if necessary (tabs between columns). 
    Returns the original path if no changes were needed; otherwise writes a normalized `.tmp` file and returns its path."""
    if chunk_size <= 0:
        raise ValueError("chunk_size must be a positive integer")
    opener = gzip.open if str(file_path).endswith(".gz") else open
    mode = "rt" if opener is gzip.open else "r"
    # First pass: detect if normalization is needed without loading entire file
    modified = False
    header_found = False
    with opener(file_path, mode, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue  # meta-info lines can be left as-is
            if line.startswith("#"):
                header_found = True
                # Check if header line has correct tab separation
                if re.sub(r"\s+", "\t", line.rstrip()) + "\n" != line:
                    modified = True
            elif header_found:
                # For data lines after header, check delimiter normalization
                if re.sub(r"\s+", "\t", line.rstrip()) + "\n" != line:
                    modified = True
            # if header not reached, continue scanning without modification check
    if not modified:
        return file_path  # no changes needed
    # Second pass: rewrite file with normalized tabs
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
                buffer.append(line)  # preserve meta-info lines exactly
            elif line.startswith("#"):
                header_found = True
                buffer.append(re.sub(r"\s+", "\t", line.rstrip()) + "\n")
            else:
                if header_found:
                    # normalize delimiters in variant data lines
                    buffer.append(re.sub(r"\s+", "\t", line.rstrip()) + "\n")
                else:
                    # (In case file has weird content before header)
                    buffer.append(line)
            if len(buffer) >= chunk_size:
                flush_buffer()
        flush_buffer()
    return temp_file

def _create_missing_call_factory(format_keys: Sequence[str], header) -> Callable[[], dict]:
    """Create a factory for default call data given a set of FORMAT keys."""
    if not format_keys:
        return lambda: {}
    template: dict = {}
    for key in format_keys:
        if key == "GT":
            template[key] = "./."
            continue
        # Determine default value based on the FORMAT definition's number
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
    """Ensure that standard AC/AN/AF INFO header definitions are present in the VCF header."""
    info_defs = (
        (
            "AC",
            {
                "ID": "AC",
                "Number": "A",
                "Type": "Integer",
                "Description": "Alternate allele count in genotypes, for each ALT allele",
            },
        ),
        (
            "AN",
            {
                "ID": "AN",
                "Number": "1",
                "Type": "Integer",
                "Description": "Total number of alleles in called genotypes",
            },
        ),
        (
            "AF",
            {
                "ID": "AF",
                "Number": "A",
                "Type": "Float",
                "Description": "Allele frequency, for each ALT allele",
            },
        ),
    )
    for key, mapping in info_defs:
        try:
            existing = header.get_info_field_info(key)
        except Exception:
            existing = None
        if existing is None:
            header.add_info_line(mapping)

def _iter_called_genotype_alleles(value) -> Iterable[int]:
    """Yield all allele indices present in a genotype field value (handles complex structures)."""
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
    # Attempt to cast other types to int
    try:
        yield int(value)
    except (TypeError, ValueError):
        return

def _recompute_ac_an_af(record) -> None:
    """Recompute the AC, AN, AF INFO fields for a merged record based on its genotype calls."""
    alt_n = len(record.ALT or [])
    ac = [0] * alt_n
    an = 0
    for call in (record.calls or []):
        gt_value = (call.data or {}).get("GT")
        for allele_index in _iter_called_genotype_alleles(gt_value):
            if isinstance(allele_index, bool):
                continue  # skip boolean values (if any)
            an += 1
            if allele_index is None or allele_index <= 0:
                continue  # 0 or None = reference/missing allele, not counted in AC
            if 1 <= allele_index <= alt_n:
                ac[allele_index - 1] += 1
    record.INFO["AC"] = ac
    record.INFO["AN"] = an
    record.INFO["AF"] = ([c / an for c in ac] if an > 0 else [0.0] * alt_n)

def _remove_format_and_sample_definitions(header) -> None:
    """Strip all FORMAT definitions and sample names from the VCF header (for anonymization)."""
    if header is None:
        return
    # Remove FORMAT lines from header lines list
    if hasattr(header, "lines"):
        new_lines = []
        for line in header.lines:
            if isinstance(line, vcfpy.header.FormatHeaderLine):
                continue
            key = getattr(line, "key", None)
            if isinstance(key, str) and key.upper() == "FORMAT":
                continue
            new_lines.append(line)
        header.lines = new_lines
    # Clear the header.formats dictionary if present
    if hasattr(header, "formats"):
        try:
            header.formats.clear()
        except AttributeError:
            header.formats = {}
    # Remove all sample names from header
    if hasattr(header, "samples") and hasattr(header.samples, "names"):
        header.samples.names = []

def _pad_record_samples(record, header, sample_order: Sequence[str]) -> None:
    """Ensure that *record* has calls for all samples in sample_order and all FORMAT keys present."""
    if not sample_order or not hasattr(record, "call_for_sample"):
        return
    fmt_keys = list(record.FORMAT or [])
    factory = _create_missing_call_factory(fmt_keys, header)
    new_calls: List[vcfpy.Call] = []
    for name in sample_order:
        call = record.call_for_sample.get(name)
        if call is None:
            # If sample is missing, create a new call with default values
            call = vcfpy.Call(name, factory())
        else:
            # Fill missing format fields in existing calls with defaults
            defaults = factory()
            for key in fmt_keys:
                if key not in call.data:
                    call.data[key] = defaults.get(key)
                else:
                    # Convert tuple values to list for consistency
                    if isinstance(call.data[key], tuple):
                        call.data[key] = list(call.data[key])
        # Ensure GT is not None/empty for proper output (use default "./." if so)
        if "GT" in call.data and call.data["GT"] in {None, "", "."}:
            call.data["GT"] = (factory().get("GT", "./."))
        new_calls.append(call)
    record.update_calls(new_calls)

def _record_passes_filters(
    record,
    qual_threshold: Optional[float],
    an_threshold: Optional[float],
    allowed_filter_values: Optional[Sequence[str]],
) -> bool:
    """Check if a record meets the quality and AN thresholds and has allowed FILTER values."""
    # QUAL threshold
    if qual_threshold is not None:
        q = getattr(record, "QUAL", None)
        try:
            if q is None or float(q) <= qual_threshold:
                return False
        except (TypeError, ValueError):
            return False
    # AN (allele number) threshold
    if an_threshold is not None:
        an_val = record.INFO.get("AN")
        if isinstance(an_val, list):
            an_val = an_val[0] if an_val else None
        try:
            if an_val is None or float(an_val) <= an_threshold:
                return False
        except (TypeError, ValueError):
            return False
    # FILTER field check: only allow records with filter values in allowed_filter_values (default "PASS")
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
    """In-place filter of VCF records based on quality/AN thresholds and FILTER values."""
    import os

    tmp_path = input_path + ".filtered"
    kept = 0
    total = 0

    try:
        import vcfpy  # local import to avoid global dependency at import time

        # Read → write filtered → atomic replace
        try:
            with vcfpy.Reader.from_path(input_path) as reader:
                try:
                    with vcfpy.Writer.from_path(tmp_path, header=reader.header) as writer:
                        for rec in reader:
                            total += 1
                            if _record_passes_filters(rec, qual_threshold, an_threshold, allowed_filter_values):
                                writer.write_record(rec)
                                kept += 1
                except Exception as exc:
                    # Best-effort cleanup of partial tmp file
                    try:
                        if os.path.exists(tmp_path):
                            os.remove(tmp_path)
                    except OSError:
                        pass
                    handle_critical_error(f"Failed to write filtered VCF ({input_path}): {exc}", exc_cls=MergeVCFError)
        except Exception as exc:
            handle_critical_error(f"Failed to read VCF for filtering ({input_path}): {exc}", exc_cls=MergeVCFError)

        try:
            os.replace(tmp_path, input_path)  # atomic on same filesystem
        except Exception as exc:
            handle_critical_error(
                f"Failed to replace original with filtered VCF ({input_path}): {exc}", exc_cls=MergeVCFError
            )

        log_message(f"Applied variant filter: kept {kept} of {total}.", verbose, level=logging.INFO)

    except MergeVCFError:
        raise
    except Exception as exc:
        handle_critical_error(f"Unexpected filtering error for {input_path}: {exc}", exc_cls=MergeVCFError)


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
    """Merge multiple VCF files into one .vcf.gz and index it."""
    import os, copy, datetime, concurrent.futures
    import vcfpy, pysam

    os.makedirs(output_dir, exist_ok=True)
    timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    base_vcf = os.path.join(output_dir, f"merged_vcf_{timestamp}.vcf")
    gz_vcf = base_vcf + ".gz"

    file_count = len(valid_files)
    log_message(
        f"Discovered {file_count} validated VCF shard(s) for merging.",
        verbose,
        level=logging.INFO,
    )

    temp_files: List[str] = []
    preprocessed: List[str] = []

    # Preprocess inputs in parallel
    try:
        def _run_pre(p: str) -> str:
            return preprocess_vcf(
                p,
                qual_threshold=qual_threshold,
                an_threshold=an_threshold,
                allowed_filter_values=allowed_filter_values,
                verbose=verbose,
            )

        max_workers = min(8, max(1, len(valid_files)))
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(executor.map(_run_pre, valid_files))
    except Exception as exc:
        handle_critical_error(f"Failed to preprocess input VCF files: {exc}", exc_cls=MergeConflictError)

    for orig_path, pre_path in zip(valid_files, results):
        preprocessed.append(pre_path)
        if pre_path != orig_path:
            temp_files.append(pre_path)

    # Build merged header
    try:
        merged_header = union_headers(preprocessed, sample_order=sample_order)
    except MergeConflictError:
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to merge VCF headers: {exc}", exc_cls=MergeConflictError)

    # Apply custom metadata
    try:
        header = apply_metadata_to_header(
            merged_header,
            sample_header_line=sample_header_line,
            simple_header_lines=simple_header_lines,
            verbose=verbose,
        )
    except Exception as exc:
        handle_critical_error(f"Failed to apply metadata to header: {exc}", exc_cls=MergeConflictError)

    final_sample_order = list(sample_order) if sample_order else list(getattr(header.samples, "names", []))

    # Collect records keyed by variant
    def _alt_tuple(rec: vcfpy.Record) -> tuple[str, ...]:
        return tuple(getattr(a, "value", str(a)) for a in (rec.ALT or []))

    def _rec_key(rec: vcfpy.Record) -> tuple[str, int, str, tuple[str, ...]]:
        return (rec.CHROM, rec.POS, rec.REF, _alt_tuple(rec))

    record_store: dict[tuple, dict] = {}

    try:
        for path in preprocessed:
            with vcfpy.Reader.from_path(path) as reader:
                for rec in reader:
                    key = _rec_key(rec)
                    container = record_store.get(key)
                    if container is None:
                        container = {
                            "CHROM": rec.CHROM,
                            "POS": rec.POS,
                            "ID": copy.deepcopy(rec.ID),
                            "REF": rec.REF,
                            "ALT": copy.deepcopy(rec.ALT),
                            "QUAL": rec.QUAL,
                            "FILTER": copy.deepcopy(rec.FILTER),
                            "INFO": copy.deepcopy(rec.INFO),
                            "format_keys": list(rec.FORMAT or []),
                            "calls": {},
                        }
                        record_store[key] = container
                    else:
                        fmt = container["format_keys"]
                        for fk in rec.FORMAT or []:
                            if fk not in fmt:
                                fmt.append(fk)

                    calls_map: dict[str, dict] = container["calls"]
                    for call in rec.calls or []:
                        sname = getattr(call, "sample", None) or getattr(call, "name", None)
                        if sname:
                            calls_map[sname] = copy.deepcopy(call.data or {})
    except Exception as exc:
        handle_critical_error(f"Failed while aggregating records: {exc}", exc_cls=MergeConflictError)

    # Write merged VCF, bgzip, tabix
    try:
        def _order_format(fmt_keys: List[str]) -> List[str]:
            if "GT" in fmt_keys and fmt_keys[0] != "GT":
                fmt_keys = [k for k in fmt_keys if k != "GT"]
                fmt_keys.insert(0, "GT")
            seen, ordered = set(), []
            for k in fmt_keys:
                if k not in seen:
                    seen.add(k)
                    ordered.append(k)
            return ordered

        # sort by header contig order then position
        contig_order: dict[str, int] = {}
        try:
            from vcfpy import header as vcfhdr
            idx = 0
            for line in getattr(header, "lines", []):
                if isinstance(line, vcfhdr.ContigHeaderLine) and line.id:
                    contig_order[line.id] = idx
                    idx += 1
        except Exception:
            contig_order = {}

        def _sort_key(k: tuple[str, int, str, tuple[str, ...]]):
            chrom, pos, _, _ = k
            return (contig_order.get(chrom, 10**9), chrom, pos)

        with vcfpy.Writer.from_path(base_vcf, header=header) as writer:
            for key in sorted(record_store.keys(), key=_sort_key):
                cont = record_store[key]
                fmt_keys = _order_format(list(cont["format_keys"]))

                calls_out: List[vcfpy.Call] = []
                for s in final_sample_order:
                    data = cont["calls"].get(s, {})
                    filled: dict[str, object] = {}
                    for fk in fmt_keys:
                        filled[fk] = data.get(fk, "./." if fk == "GT" else ".")
                    calls_out.append(vcfpy.Call(sample=s, data=filled))

                out_rec = vcfpy.Record(
                    CHROM=cont["CHROM"],
                    POS=cont["POS"],
                    ID=cont["ID"],
                    REF=cont["REF"],
                    ALT=cont["ALT"],
                    QUAL=cont["QUAL"],
                    FILTER=cont["FILTER"],
                    INFO=cont["INFO"],
                    FORMAT=fmt_keys,
                    calls=calls_out,
                )
                writer.write_record(out_rec)

        pysam.tabix_compress(base_vcf, gz_vcf, force=True)
        pysam.tabix_index(gz_vcf, preset="vcf", force=True)
        try:
            os.remove(base_vcf)
        except OSError:
            pass
    except Exception as exc:
        handle_critical_error(f"Failed to write merged VCF: {exc}", exc_cls=MergeConflictError)
    finally:
        for tmp in temp_files:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass

    log_message(f"Merged VCF written and indexed: {gz_vcf}", verbose)
    return gz_vcf


    # Merge headers from all files to create a combined header
    try:
        merged_header = union_headers(preprocessed, sample_order=sample_order)
    except MergeConflictError:
        # Pass through known header merge conflicts as MergeConflictError
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to merge VCF headers: {exc}", exc_cls=MergeConflictError)

    # Apply additional metadata lines (if provided) to the merged header
    try:
        header = apply_metadata_to_header(
            merged_header,
            sample_header_line=sample_header_line,
            simple_header_lines=simple_header_lines,
            verbose=verbose,
        )
    except SystemExit:
        # allow SystemExit to propagate (if handle_critical_error called sys.exit)
        raise
    except Exception as exc:
        handle_critical_error(f"Failed to apply metadata to merged VCF header: {exc}")

    # Determine final contig order and ensure header contigs preserve this order
    contigs: List[str] = []
    if getattr(header, "contigs", None):
        contigs = list(header.contigs.keys())
    if not contigs:
        # If header.contigs is empty, gather contig names from any Contig header lines (if present)
        for ln in getattr(header, "lines", []):
            if isinstance(ln, vcfpy.header.ContigHeaderLine):
                contigs.append(ln.id)
    contig_ranks = {name: i for i, name in enumerate(contigs)}
    _ensure_info_header_lines(header)  # ensure AC/AN/AF are defined in INFO headers

    # Ensure header.formats dict is populated with all FORMAT definitions (vcfpy quirk)
    if hasattr(header, "formats") and header.formats is not None:
        for ln in getattr(header, "lines", []):
            if isinstance(ln, vcfpy.header.FormatHeaderLine):
                header.formats[ln.id] = ln

    # Initialize sample list from merged header; include all samples from all files
    sample_names: List[str] = list(getattr(header.samples, "names", []) or [])

    # Open each VCF file with vcfpy Reader and prepare for k-way merge
    readers: List[vcfpy.Reader] = []
    iters: List[Iterator] = []
    for path in preprocessed:
        try:
            reader = vcfpy.Reader.from_path(path)
        except Exception as exc:
            handle_critical_error(f"Failed to open VCF for merging ({path}): {exc}")
        readers.append(reader)
        iters.append(iter(reader))
        # Add any new sample names found in this file to the master sample list
        file_samples = list(getattr(reader.header.samples, "names", []) or [])
        for name in file_samples:
            if name not in sample_names:
                sample_names.append(name)
    # Update the merged header's sample names to include all discovered samples (or maintain given order)
    if hasattr(header, "samples") and hasattr(header.samples, "names"):
        header.samples.names = list(sample_names)

    # Create output VCF writer with the merged header
    try:
        writer = vcfpy.Writer.from_path(base_vcf, header)
    except Exception as exc:
        # Clean up readers before aborting
        for r in readers:
            try:
                r.close()
            except Exception:
                pass
        handle_critical_error(f"Failed to open writer for merged VCF: {exc}", exc_cls=MergeConflictError)

    def _advance(reader_index: int) -> Optional[vcfpy.Record]:
        """Fetch the next record from reader at index, or return None if exhausted."""
        try:
            rec = next(iters[reader_index])
        except StopIteration:
            return None
        _sort_record_alts(rec)
        return rec

    # Heap-based k-way merge: each heap entry is (sort_key, counter, source_index, record)
    log_message(
        "Performing heap-based k-way merge of VCF shards...",
        verbose,
        level=logging.DEBUG,
    )
    heap: List[Tuple[Tuple[int, int, str, Tuple[str, ...]], int, int, vcfpy.Record]] = []
    counter = itertools.count()  # unique counter to avoid comparison of records in heap
    # Initialize heap with the first record from each reader (if available)
    for i in range(len(iters)):
        rec = _advance(i)
        if rec is not None:
            heapq.heappush(heap, (_record_sort_key(rec, contig_ranks), next(counter), i, rec))

    # Merge loop
    try:
        while heap:
            key, _, src_idx, rec = heapq.heappop(heap)
            # Collect all records with the same sort key (colliding variants)
            colliding_records: List[Tuple[vcfpy.Record, int]] = [(rec, src_idx)]
            # Also gather any other entries from heap that have identical key
            while heap and heap[0][0] == key:
                _, _, j, rec2 = heapq.heappop(heap)
                colliding_records.append((rec2, j))
            # For each source that contributed a colliding record, advance it to find if more colliding records follow
            pending: List[Tuple[vcfpy.Record, int]] = []
            cur_index = 0
            while cur_index < len(colliding_records):
                _, source_idx = colliding_records[cur_index]
                nxt = _advance(source_idx)
                # If the next record from this source is still the same variant (same sort key), merge it as well
                while nxt is not None and _record_sort_key(nxt, contig_ranks) == key:
                    colliding_records.append((nxt, source_idx))
                    nxt = _advance(source_idx)
                if nxt is not None:
                    pending.append((nxt, source_idx))
                cur_index += 1
            # Merge all colliding records into one and write to output
            merged_rec = _merge_colliding_records(colliding_records, header, sample_names)
            _recompute_ac_an_af(merged_rec)
            writer.write_record(merged_rec)
            # Push any pending next records (that were not colliding) onto the heap
            for nxt_rec, ridx in pending:
                heapq.heappush(heap, (_record_sort_key(nxt_rec, contig_ranks), next(counter), ridx, nxt_rec))
    finally:
        # Close writer and all readers, remove any temporary files
        writer.close()
        for r in readers:
            try:
                r.close()
            except Exception:
                pass
        for tmp in temp_files:
            try:
                if os.path.exists(tmp):
                    os.remove(tmp)
            except OSError:
                pass

    # Apply in-place variant filtering on the merged VCF file (QUAL/AN thresholds, FILTER values)
    _filter_vcf_records(
        base_vcf,
        qual_threshold=qual_threshold,
        an_threshold=an_threshold,
        allowed_filter_values=allowed_filter_values,
        verbose=verbose,
    )

    # Anonymize the VCF: drop all sample genotype data and FORMAT definitions
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
        handle_critical_error(f"Failed to open anonymized VCF writer: {exc}")
    try:
        for rec in rdr:
            rec.FORMAT = []    # no format fields
            rec.calls = []     # no sample calls
            w.write_record(rec)
    finally:
        rdr.close()
        w.close()
    # Replace original merged VCF with the anonymized version
    shutil.move(anon_tmp, base_vcf)

    # BGZF compress + Tabix index using pysam
    log_message(
        "Compressing and indexing the final VCF...",
        verbose,
        level=logging.INFO,
    )
    try:
        pysam.tabix_compress(base_vcf, gz_vcf, force=True)
        pysam.tabix_index(gz_vcf, preset="vcf", force=True)
        os.remove(base_vcf)  # remove uncompressed VCF after successful compression
    except Exception as exc:
        handle_critical_error(
            f"Failed to compress/index final VCF: {exc}",
            exc_cls=MergeConflictError,
        )
    log_message(
        f"Merged VCF created and indexed successfully: {gz_vcf}",
        verbose,
        level=logging.INFO,
    )
    return gz_vcf

def union_headers(valid_files: Sequence[str], sample_order: Optional[Sequence[str]] = None):
    """Merge the headers of all VCF files in `valid_files` into a combined header."""
    combined_header = None
    # Dictionaries to keep track of seen header lines by ID
    info_lines: dict[str, vcfpy.header.InfoHeaderLine] = {}
    filter_lines: dict[str, vcfpy.header.FilterHeaderLine] = {}
    contig_lines: dict[str, vcfpy.header.ContigHeaderLine] = {}
    format_lines: dict[str, vcfpy.header.FormatHeaderLine] = {}
    computed_sample_order: List[str] = []
    merged_sample_metadata = None

    def _line_mapping(line) -> dict[str, str]:
        """Helper to get an ordered mapping of a header line's metadata (ID, Number, Type, etc)."""
        mapping = getattr(line, "mapping", None)
        return dict(mapping) if isinstance(mapping, dict) else {}

    for file_path in valid_files:
        pre = preprocess_vcf(file_path)  # ensure file is normalized (idempotent if already done)
        reader = None
        try:
            reader = vcfpy.Reader.from_path(pre)
            hdr = reader.header
            if combined_header is None:
                # Initialize combined_header as a copy of the first file's header
                combined_header = hdr.copy()
                # Initialize sample order and header line dictionaries from first header
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
                        # Keep sample metadata line for merging
                        mapping = _line_mapping(line)
                        if mapping.get("ID"):
                            merged_sample_metadata = dict(mapping)
                continue  # move to next file
            # For subsequent files, merge header lines:
            if hasattr(hdr, "samples") and hasattr(hdr.samples, "names"):
                for s in hdr.samples.names:
                    if s not in computed_sample_order:
                        computed_sample_order.append(s)
            for line in hdr.lines:
                if isinstance(line, vcfpy.header.InfoHeaderLine):
                    existing = info_lines.get(line.id)
                    if existing is None:
                        info_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                    else:
                        # Ensure no conflicting definitions for same INFO ID
                        exist_map, new_map = _line_mapping(existing), _line_mapping(line)
                        for field in ("Number", "Type", "Description"):
                            if exist_map.get(field) != new_map.get(field):
                                handle_critical_error(
                                    "INFO header definitions conflict across shards. "
                                    f"Field '{field}' for INFO '{line.id}' differs: "
                                    f"{exist_map.get(field)!r} vs {new_map.get(field)!r} (in {file_path}).",
                                    exc_cls=MergeConflictError,
                                )
                elif isinstance(line, vcfpy.header.FilterHeaderLine):
                    existing = filter_lines.get(line.id)
                    if existing is None:
                        filter_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                    else:
                        # Ensure no conflicting descriptions for same FILTER ID
                        exist_map, new_map = _line_mapping(existing), _line_mapping(line)
                        if exist_map.get("Description") != new_map.get("Description"):
                            handle_critical_error(
                                "FILTER header definitions conflict across shards. "
                                f"Description for FILTER '{line.id}' differs: "
                                f"{exist_map.get('Description')!r} vs {new_map.get('Description')!r} (in {file_path}).",
                                exc_cls=MergeConflictError,
                            )
                elif isinstance(line, vcfpy.header.ContigHeaderLine):
                    existing = contig_lines.get(line.id)
                    if existing is None:
                        contig_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                    else:
                        # Ensure no conflicting contig definitions for same ID
                        if _line_mapping(existing) != _line_mapping(line):
                            handle_critical_error(
                                "Contig header definitions conflict across shards. "
                                f"Contig '{line.id}' differs between shards (conflict in {file_path}).",
                                exc_cls=MergeConflictError,
                            )
                elif isinstance(line, vcfpy.header.FormatHeaderLine):
                    existing = format_lines.get(line.id)
                    if existing is None:
                        format_lines[line.id] = copy.deepcopy(line)
                        combined_header.add_line(copy.deepcopy(line))
                    else:
                        # Ensure no conflicting definitions for same FORMAT ID
                        exist_map, new_map = _line_mapping(existing), _line_mapping(line)
                        for field in ("Number", "Type", "Description"):
                            if exist_map.get(field) != new_map.get(field):
                                handle_critical_error(
                                    "FORMAT header definitions conflict across shards. "
                                    f"Field '{field}' for FORMAT '{line.id}' differs: "
                                    f"{exist_map.get(field)!r} vs {new_map.get(field)!r} (in {file_path}).",
                                    exc_cls=MergeConflictError,
                                )
                elif isinstance(line, vcfpy.header.SampleHeaderLine):
                    # Merge sample-level metadata (e.g., SAMPLE lines)
                    mapping = _line_mapping(line)
                    sid = mapping.get("ID")
                    if not sid:
                        continue
                    if merged_sample_metadata is None:
                        merged_sample_metadata = dict(mapping)
                    else:
                        for k, v in mapping.items():
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
            # Remove temp file if one was created for normalization
            if pre != file_path and os.path.exists(pre):
                os.remove(pre)
    if combined_header is None:
        handle_critical_error("Unable to construct a merged VCF header.", exc_cls=MergeConflictError)
    # Determine final sample ordering for the combined header
    target_samples = sample_order if sample_order is not None else computed_sample_order
    if target_samples and hasattr(combined_header, "samples") and hasattr(combined_header.samples, "names"):
        combined_header.samples.names = list(target_samples)
    # If sample metadata was merged, add a single combined SAMPLE header line
    if merged_sample_metadata:
        serialized = build_sample_metadata_line(merged_sample_metadata)
        parsed = _parse_sample_metadata_line(serialized)
        sample_line = vcfpy.header.SampleHeaderLine.from_mapping(parsed)
        # Remove any existing SAMPLE lines and replace with the merged one
        combined_header.lines = [ln for ln in combined_header.lines if getattr(ln, "key", None) != "SAMPLE"]
        combined_header.add_line(sample_line)
    # Ensure contig order is applied to combined header
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