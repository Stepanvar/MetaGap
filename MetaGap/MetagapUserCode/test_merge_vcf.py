#!/usr/bin/env python3
"""Backward compatible wrapper for the :mod:`merge_vcf` package."""

from __future__ import annotations

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

if __package__ in {None, ""}:
    current_dir = os.path.dirname(os.path.abspath(__file__))
    if current_dir not in sys.path:
        sys.path.append(current_dir)

try:  # pragma: no cover - import guard for script/module execution
    from .io_utils import parse_simple_metadata_line, preprocess_vcf
except ImportError:  # pragma: no cover - executed when running as a script
    from io_utils import parse_simple_metadata_line, preprocess_vcf  # type: ignore

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

# Configure logging
LOG_FILE = "script_execution.log"
logger = logging.getLogger("vcf_merger")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
fh = logging.FileHandler(LOG_FILE)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)


def log_message(message, verbose=False):
    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message):
    log_message("CRITICAL ERROR: " + message)
    print("A critical error occurred. Check {} for details.".format(LOG_FILE))
    sys.exit(1)


def handle_non_critical_error(message):
    log_message("WARNING: " + message)
    print("Warning: " + message)


def is_gvcf_header(header_lines):
    """Return True if the provided header metadata indicates a gVCF file."""

    for raw_line in header_lines:
        line_text = None
        if hasattr(raw_line, "to_line"):
            try:
                line_text = raw_line.to_line()
            except Exception:
                line_text = None
        if line_text is None and isinstance(raw_line, str):
            line_text = raw_line

        if line_text:
            stripped = line_text.strip()
            if stripped.startswith("##"):
                normalized = stripped.upper().replace(" ", "")
                if normalized.startswith("##GVCF"):
                    return True
                if normalized.startswith("##ALT=") and (
                    "ID=<NON_REF>" in normalized or "ID=NON_REF" in normalized
                ):
                    return True

        key = getattr(raw_line, "key", None)
        if isinstance(key, str):
            upper_key = key.upper()
            if upper_key == "GVCF":
                return True
            if upper_key == "ALT":
                candidates = [
                    getattr(raw_line, "id", None),
                    getattr(raw_line, "value", None),
                ]
                mapping = getattr(raw_line, "mapping", None)
                if isinstance(mapping, dict):
                    candidates.extend(mapping.values())
                for candidate in candidates:
                    if not isinstance(candidate, str):
                        continue
                    normalized_candidate = candidate.upper().replace(" ", "")
                    if "<NON_REF>" in normalized_candidate or normalized_candidate == "NON_REF":
                        return True
    return False


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
        parsed = parse_simple_metadata_line(normalized)
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




def normalize_vcf_version(version_value):
    """Extract the numeric portion of the VCF version if present."""
    if not version_value:
        return version_value
    match = re.search(r"(\d+\.\d+)", version_value)
    if match:
        return match.group(1)
    return version_value.strip()


def find_first_vcf_with_header(input_dir, verbose=False):
    """Locate the first readable VCF file and return its header metadata."""
    candidates = sorted(glob.glob(os.path.join(input_dir, "*.vcf")))
    for file_path in candidates:
        try:
            reader = vcfpy.Reader.from_path(file_path)
        except Exception as e:
            handle_non_critical_error(
                f"Unable to open {file_path} for auto-detection: {str(e)}. Trying next file.",
            )
            continue

        header = reader.header
        fileformat = None
        reference = None
        for line in header.lines:
            if line.key == "fileformat" and fileformat is None:
                fileformat = line.value
            elif line.key == "reference" and reference is None:
                reference = line.value
            if fileformat and reference:
                break

        if fileformat and reference:
            log_message(
                f"Auto-detection using {file_path}: fileformat={fileformat}, reference={reference}",
                verbose,
            )
            return file_path, fileformat, reference

        handle_non_critical_error(
            f"Required metadata missing in {file_path} during auto-detection. Trying next file.",
        )

    return None, None, None


def validate_vcf(file_path, ref_genome, vcf_version, verbose=False, allow_gvcf=False):
    log_message(f"Validating file: {file_path}", verbose)
    if not os.path.isfile(file_path):
        handle_non_critical_error(f"File {file_path} does not exist. Skipping.")
        return False

    preprocessed_file = preprocess_vcf(file_path)
    try:
        reader = vcfpy.Reader.from_path(preprocessed_file)
    except Exception as e:
        handle_non_critical_error(f"Could not open {file_path}: {str(e)}. Skipping.")
        return False
    if preprocessed_file != file_path:
        os.remove(preprocessed_file)


    header = reader.header

    if is_gvcf_header(header.lines) and not allow_gvcf:
        handle_non_critical_error(
            f"{file_path} appears to be a gVCF based on header metadata. Use --allow-gvcf to include gVCFs. Skipping."
        )
        return False

    fileformat = None
    for line in header.lines:
        if line.key == "fileformat":
            fileformat = line.value
            break

    expected_version = normalize_vcf_version(vcf_version)
    actual_version = normalize_vcf_version(fileformat)

    if fileformat is None or expected_version != actual_version:
        handle_non_critical_error(
            f"VCF version mismatch in {file_path}. Expected: {vcf_version}, Found: {fileformat}. Skipping.",
        )
        return False

    ref = None
    for line in header.lines:
        if line.key == "reference":
            ref = line.value
            break
    if ref is None or ref != ref_genome:
        handle_non_critical_error(f"Reference genome mismatch in {file_path}. Expected: {ref_genome}, Found: {ref}. Skipping.")
        return False


    try:
        from itertools import islice

        if not allow_gvcf:
            for record in islice(reader, 5):
                if any(
                    (
                        (alt_value := getattr(alt, "value", "")) in {"<NON_REF>", "NON_REF"}
                        or str(alt) in {"<NON_REF>", "SymbolicAllele('NON_REF')"}
                    )
                    for alt in getattr(record, "ALT", [])
                ):
                    handle_non_critical_error(
                        f"{file_path} appears to be a gVCF (found <NON_REF> ALT allele). Use --allow-gvcf to include gVCFs. Skipping."
                    )
                    return False

        _ = next(reader)
    except StopIteration:
        pass
    except Exception as e:
        handle_non_critical_error(f"Structural integrity check failed for {file_path}: {str(e)}. Skipping.")
        return False

    log_message(f"Validation passed for {file_path}", verbose)
    return True


def _extract_header_definitions(header, key):
    """Return an OrderedDict mapping header line IDs to their metadata."""

    try:
        header_lines = list(header.get_lines(key))
    except Exception:  # pragma: no cover - defensive; vcfpy may raise AttributeError
        header_lines = []

    definitions = OrderedDict()
    for line in header_lines:
        line_id = getattr(line, "id", None)
        if line_id is None and hasattr(line, "mapping"):
            line_id = line.mapping.get("ID")
        if line_id is None:
            continue
        metadata = copy.deepcopy(getattr(line, "mapping", {}))
        definitions[line_id] = metadata
    return definitions


def validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose=False, allow_gvcf=False):
    valid_vcfs = []
    log_message(f"Validating all VCF files in {input_dir}", verbose)
    reference_samples = []
    reference_contigs = None
    reference_info_defs = None
    reference_format_defs = None
    for file_path in glob.glob(os.path.join(input_dir, "*.vcf")):
        try:
            with open(file_path, "r", encoding="utf-8") as raw_vcf:
                header_lines = []
                for raw_line in raw_vcf:
                    if not raw_line.startswith("#"):
                        break
                    if raw_line.startswith("##"):
                        header_lines.append(raw_line.rstrip("\n"))
        except OSError as exc:
            handle_non_critical_error(f"Could not read header from {file_path}: {exc}. Skipping.")
            continue

        if is_gvcf_header(header_lines) and not allow_gvcf:
            handle_non_critical_error(
                f"{file_path} appears to be a gVCF based on header metadata. Use --allow-gvcf to include gVCFs. Skipping."
            )
            continue

        if validate_vcf(file_path, ref_genome, vcf_version, verbose=verbose, allow_gvcf=allow_gvcf):
            preprocessed_file = preprocess_vcf(file_path)
            reader = None
            try:
                reader = vcfpy.Reader.from_path(preprocessed_file)
            except Exception as exc:
                handle_critical_error(
                    f"Failed to reopen {file_path} after validation: {exc}."
                )
            try:
                header = reader.header

                current_samples = list(header.samples.names)
                current_contigs = _extract_header_definitions(header, "contig")
                current_info_defs = _extract_header_definitions(header, "INFO")
                current_format_defs = _extract_header_definitions(header, "FORMAT")

                if not reference_samples:
                    reference_samples = list(current_samples)
                    reference_contigs = current_contigs
                    reference_info_defs = current_info_defs
                    reference_format_defs = current_format_defs
                else:
                    new_samples = [
                        sample for sample in current_samples if sample not in reference_samples
                    ]
                    if new_samples:
                        reference_samples.extend(new_samples)
                        handle_non_critical_error(
                            "Sample columns differ between VCF shards. The merged VCF will "
                            "include the union of all sample names. Newly observed samples "
                            f"{new_samples} were appended to the merged header."
                        )

                    missing_samples = [
                        sample for sample in reference_samples if sample not in current_samples
                    ]
                    if missing_samples:
                        handle_non_critical_error(
                            f"File {file_path} is missing sample columns {missing_samples}. "
                            "These entries will be filled with missing ('.') values during merge."
                        )
                    if current_contigs != reference_contigs:
                        handle_critical_error(
                            "Contig definitions differ between VCF shards. MetaGap assumes vertical "
                            "concatenation of shards, so headers must match across shards. "
                            f"Expected contigs: {list(reference_contigs.keys())}; found in {file_path}: "
                            f"{list(current_contigs.keys())}."
                        )
                    if current_info_defs != reference_info_defs:
                        handle_critical_error(
                            "INFO field definitions differ between VCF shards. MetaGap assumes "
                            "vertical concatenation of shards, so headers must match across shards. "
                            f"Expected INFO IDs: {list(reference_info_defs.keys())}; found in {file_path}: "
                            f"{list(current_info_defs.keys())}."
                        )
                    if current_format_defs != reference_format_defs:
                        handle_critical_error(
                            "FORMAT field definitions differ between VCF shards. MetaGap assumes "
                            "vertical concatenation of shards, so headers must match across shards. "
                            f"Expected FORMAT IDs: {list(reference_format_defs.keys())}; found in {file_path}: "
                            f"{list(current_format_defs.keys())}."
                        )

                contig_order = {name: idx for idx, name in enumerate(current_contigs.keys())}
                last_contig_index = None
                last_position = None
                last_contig_name = None
                for record in reader:
                    chrom = record.CHROM
                    pos = record.POS
                    if contig_order:
                        if chrom not in contig_order:
                            handle_critical_error(
                                f"Encountered unknown contig '{chrom}' in {file_path}. MetaGap assumes vertical "
                                "concatenation of shards with consistent headers."
                            )
                        contig_index = contig_order[chrom]
                        if last_contig_index is None:
                            last_contig_index = contig_index
                            last_position = pos
                            last_contig_name = chrom
                            continue

                        if contig_index < last_contig_index or (
                            contig_index == last_contig_index and pos < last_position
                        ):
                            handle_critical_error(
                                f"VCF shard {file_path} is not coordinate-sorted: encountered {chrom}:{pos} after "
                                f"{last_contig_name}:{last_position}. MetaGap assumes vertical concatenation of shards, "
                                "so input shards must already be sorted."
                            )

                        last_contig_index = contig_index
                    else:
                        if last_contig_name is None:
                            last_contig_name = chrom
                            last_position = pos
                            continue

                        if chrom == last_contig_name:
                            if pos < last_position:
                                handle_critical_error(
                                    f"VCF shard {file_path} is not coordinate-sorted: encountered {chrom}:{pos} after "
                                    f"{last_contig_name}:{last_position}. MetaGap assumes vertical concatenation of shards, "
                                    "so input shards must already be sorted."
                                )
                        else:
                            last_contig_name = chrom
                            last_position = pos
                            continue

                    last_contig_name = chrom
                    last_position = pos
                valid_vcfs.append(file_path)
            finally:
                if reader is not None and hasattr(reader, "close"):
                    try:
                        reader.close()
                    except Exception:
                        pass
                if (
                    preprocessed_file != file_path
                    and preprocessed_file is not None
                    and os.path.exists(preprocessed_file)
                ):
                    os.remove(preprocessed_file)
        else:
            log_message(f"File {file_path} failed validation and is skipped.", verbose)
    if not valid_vcfs:
        handle_critical_error("No valid VCF files remain after validation. Aborting.")
    log_message("Validation completed. Valid VCF files: " + ", ".join(valid_vcfs), verbose)
    return valid_vcfs, reference_samples


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
    ) = parse_metadata_arguments(
        args,
        verbose,
        log_message=log_message,
        handle_critical_error=handle_critical_error,
    )

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

# Ensure the repository root is on ``sys.path`` so the ``merge_vcf`` package
# can be imported when this file is executed directly.
_THIS_FILE = Path(__file__).resolve()
_REPO_ROOT = _THIS_FILE.parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from merge_vcf.main import *  # noqa: F401,F403 - re-export legacy symbols
from merge_vcf import main as _main


configure_merging_callbacks(
    MergeCallbacks(
        preprocess_vcf=preprocess_vcf,
        log_message=log_message,
        handle_critical_error=handle_critical_error,
        build_sample_metadata_line=build_sample_metadata_line,
        parse_sample_metadata_line=_parse_sample_metadata_line,
        apply_metadata_to_header=apply_metadata_to_header,
        vcfpy=vcfpy,
    )
)


if __name__ == "__main__":
    main()
