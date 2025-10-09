#!/usr/bin/env python3
"""Public CLI surface for the MetaGap VCF merging workflow."""

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
import csv
import logging
import datetime
import os
import subprocess
import gzip
import importlib
import importlib.util
import types
from collections import OrderedDict

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

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from merge_vcf import VCFPY_AVAILABLE, vcfpy
from merge_vcf import metadata as _metadata
from merge_vcf import cli as _cli
from merge_vcf.logging_utils import (
    LOG_FILE,
    handle_critical_error,
    handle_non_critical_error,
    log_message,
    logger,
)
from merge_vcf.metadata import (
    _format_sample_metadata_value,
    _parse_sample_metadata_line,
    _parse_simple_metadata_line,
    _validate_anonymized_vcf_header as _metadata_validate_anonymized_vcf_header,
    append_metadata_to_merged_vcf,
    apply_metadata_to_header,
    build_sample_metadata_line,
    parse_metadata_arguments,
    recalculate_cohort_info_tags,
)
from merge_vcf.merging import (
    _create_missing_call_factory,
    _pad_record_samples,
    _remove_format_and_sample_definitions,
    merge_vcfs,
    preprocess_vcf,
    union_headers,
)
from merge_vcf.validation import (
    find_first_vcf_with_header,
    is_gvcf_header,
    normalize_vcf_version,
    validate_all_vcfs,
    validate_merged_vcf,
    validate_vcf,
)


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

    return _cli.parse_arguments()


def summarize_produced_vcfs(output_dir, fallback_vcf):
    return _cli.summarize_produced_vcfs(output_dir, fallback_vcf)


def _validate_anonymized_vcf_header(
    final_vcf_path, ensure_for_uncompressed: bool = False
):
    original = _metadata.handle_critical_error
    _metadata.handle_critical_error = handle_critical_error
    try:
        return _metadata_validate_anonymized_vcf_header(
            final_vcf_path, ensure_for_uncompressed=ensure_for_uncompressed
        )
    finally:
        _metadata.handle_critical_error = original


def run_workflow(args):
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
        detected_file, detected_fileformat, detected_reference = find_first_vcf_with_header(
            input_dir, verbose
        )
        if not ref_genome and detected_reference:
            ref_genome = detected_reference
        if not vcf_version and detected_fileformat:
            vcf_version = normalize_vcf_version(detected_fileformat)

    if not ref_genome:
        handle_critical_error(
            "Reference genome build must be provided via --ref or auto-detectable from input files."
        )
    if not vcf_version:
        handle_critical_error(
            "VCF version must be provided via --vcf-version or auto-detectable from input files."
        )

    if detection_needed and detected_file:
        log_message(
            f"Auto-detected metadata from {detected_file} -> reference={detected_reference or 'unknown'}, "
            f"version={detected_fileformat or 'unknown'}",
            verbose,
        )

    log_message(f"Reference genome: {ref_genome}, VCF version: {vcf_version}", verbose)

    valid_files, sample_order = validate_all_vcfs(
        input_dir,
        ref_genome,
        vcf_version,
        verbose,
        allow_gvcf=getattr(args, "allow_gvcf", False),
    )

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
    return summary_path


def _load_cli_module():
    current_module = sys.modules.get(__name__)
    if current_module is not None:
        for alias in ("MetagapUserCode.test_merge_vcf", "test_merge_vcf"):
            sys.modules.setdefault(alias, current_module)

    for name in ("MetagapUserCode.cli", "cli"):
        module = sys.modules.get(name)
        if module is not None:
            return module

    for name in ("MetagapUserCode.cli", "cli"):
        try:
            return importlib.import_module(name)
        except ModuleNotFoundError:
            continue

    cli_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "cli.py")
    spec = importlib.util.spec_from_file_location("cli", cli_path)
    if spec is None or spec.loader is None:
        raise ImportError("Unable to locate CLI module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules.setdefault("cli", module)
    return module


def parse_arguments():
    cli_module = _load_cli_module()
    return cli_module.parse_arguments()


def main():
    cli_module = _load_cli_module()
    args = parse_arguments()
    workflow_module = sys.modules.get(__name__)
    if workflow_module is None:
        workflow_module = types.SimpleNamespace(
            run_workflow=run_workflow,
            parse_arguments=parse_arguments,
        )
    return cli_module.main(args, workflow_module=workflow_module)


if __name__ == "__main__":
    main()
