import gzip
import importlib.util
import os
import sys
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest


MODULE_PATH = (
    Path(__file__).resolve().parents[2] / "MetagapUserCode" / "test_merge_vcf.py"
)


def load_user_module():
    spec = importlib.util.spec_from_file_location("user_test_merge_vcf", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@dataclass
class _HeaderDefinition:
    text: str

    def to_line(self):
        return self.text

    def copy(self):
        return _HeaderDefinition(self.text)


class _SimpleHeaderLine(_HeaderDefinition):
    def __init__(self, key: str, value: str):
        self.key = key
        self.value = value
        super().__init__(f"##{key}={value}")

    def copy(self):
        return _SimpleHeaderLine(self.key, self.value)


class _InfoHeaderLine(_HeaderDefinition):
    def __init__(self, mapping):
        self.key = "INFO"
        self.mapping = dict(mapping)
        self.id = self.mapping.get("ID")
        parts = ",".join(f"{k}={v}" for k, v in self.mapping.items())
        super().__init__(f"##INFO=<{parts}>")

    def copy(self):
        return _InfoHeaderLine(self.mapping)


class _FormatHeaderLine(_HeaderDefinition):
    def __init__(self, mapping):
        self.key = "FORMAT"
        self.mapping = dict(mapping)
        self.id = self.mapping.get("ID")
        self.number = self.mapping.get("Number") or self.mapping.get("NUMBER")
        parts = ",".join(f"{k}={v}" for k, v in self.mapping.items())
        super().__init__(f"##FORMAT=<{parts}>")

    def copy(self):
        return _FormatHeaderLine(self.mapping)


class _SampleHeaderLine(_HeaderDefinition):
    def __init__(self, mapping):
        self.key = "SAMPLE"
        self.mapping = OrderedDict(mapping)
        self.id = self.mapping.get("ID")
        parts = ",".join(f"{k}={v}" for k, v in self.mapping.items())
        super().__init__(f"##SAMPLE=<{parts}>")

    def copy(self):
        return _SampleHeaderLine(self.mapping)


class _FakeCall:
    def __init__(self, name, data):
        self.name = name
        self.data = {key: value for key, value in data.items()}

    def copy(self):
        return _FakeCall(self.name, self.data.copy())

    def to_string(self, format_keys):
        values = []
        for key in format_keys:
            raw = self.data.get(key)
            if isinstance(raw, list):
                formatted = ",".join("." if val is None else str(val) for val in raw)
            elif raw in {None, ""}:
                formatted = "."
            else:
                formatted = str(raw)
            values.append(formatted)
        return ":".join(values) if format_keys else "."


class _FakeRecord:
    def __init__(
        self,
        chrom,
        pos,
        record_id,
        ref,
        alts,
        qual,
        filt,
        info,
        format_keys,
        calls,
    ):
        self.CHROM = chrom
        self.POS = pos
        self.ID = record_id
        self.REF = ref
        self.ALT = list(alts)
        self.QUAL = qual
        self.FILTER = filt
        self.INFO = OrderedDict(info)
        self.FORMAT = list(format_keys)
        self.calls = list(calls)
        self.call_for_sample = {call.name: call for call in self.calls}

    def copy(self):
        return _FakeRecord(
            self.CHROM,
            self.POS,
            self.ID,
            self.REF,
            list(self.ALT),
            self.QUAL,
            self.FILTER,
            OrderedDict((k, v if not isinstance(v, list) else list(v)) for k, v in self.INFO.items()),
            list(self.FORMAT),
            [call.copy() for call in self.calls],
        )

    def update_calls(self, updated_calls):
        self.calls = list(updated_calls)
        self.call_for_sample = {call.name: call for call in self.calls}

    def to_line(self):
        info_parts = []
        for key, value in self.INFO.items():
            if isinstance(value, list):
                serialized = ",".join(str(entry) for entry in value)
            else:
                serialized = str(value)
            info_parts.append(f"{key}={serialized}")
        info_field = ";".join(info_parts) if info_parts else "."

        fields = [
            self.CHROM,
            str(self.POS),
            self.ID,
            self.REF,
            ",".join(self.ALT) if self.ALT else ".",
            self.QUAL,
            self.FILTER,
            info_field,
        ]

        if self.FORMAT:
            fields.append(":".join(self.FORMAT))
            fields.extend(call.to_string(self.FORMAT) for call in self.calls)

        return "\t".join(fields)


class _FakeHeader:
    def __init__(self, lines, samples, columns_line, format_fields):
        self.lines = [line.copy() for line in lines]
        self.samples = SimpleNamespace(names=list(samples))
        self._columns_line = columns_line
        self._format_fields = {key: value.copy() for key, value in format_fields.items()}

    def add_line(self, line):
        self.lines.append(line)
        if isinstance(line, _FormatHeaderLine) and line.id:
            self._format_fields[line.id] = line

    def copy(self):
        return _FakeHeader(self.lines, self.samples.names, self._columns_line, self._format_fields)

    def get_format_field_info(self, key):
        return self._format_fields.get(key)

    def to_columns_line(self):
        return self._columns_line

    def get_lines(self, key):
        key_upper = key.upper()
        return [line for line in self.lines if getattr(line, "key", "").upper() == key_upper]

    @property
    def info_ids(self):
        entries = OrderedDict()
        for line in self.lines:
            if isinstance(line, _InfoHeaderLine) and line.id:
                entries[line.id] = line.mapping
        return entries


_VCF_STORAGE = {}


def _parse_mapping_body(body: str):
    mapping = OrderedDict()
    token = []
    in_quotes = False
    escape = False

    def _flush():
        if not token:
            return
        raw = "".join(token).strip()
        token.clear()
        if not raw or "=" not in raw:
            return
        key, value = raw.split("=", 1)
        mapping[key.strip()] = value.strip()

    for char in body:
        if escape:
            token.append(char)
            escape = False
            continue
        if char == "\\":
            token.append(char)
            escape = True
            continue
        if char == '"':
            in_quotes = not in_quotes
            token.append(char)
            continue
        if char == "," and not in_quotes:
            _flush()
            continue
        token.append(char)
    _flush()
    return mapping


def _header_from_text(text: str):
    stripped = text.strip()
    if stripped.startswith("##INFO=<") and stripped.endswith(">"):
        mapping = _parse_mapping_body(stripped[len("##INFO=<") : -1])
        return _InfoHeaderLine(mapping)
    if stripped.startswith("##FORMAT=<") and stripped.endswith(">"):
        mapping = _parse_mapping_body(stripped[len("##FORMAT=<") : -1])
        return _FormatHeaderLine(mapping)
    if stripped.startswith("##SAMPLE=<") and stripped.endswith(">"):
        mapping = _parse_mapping_body(stripped[len("##SAMPLE=<") : -1])
        return _SampleHeaderLine(mapping)
    if stripped.startswith("##") and "=" in stripped[2:]:
        key, value = stripped[2:].split("=", 1)
        return _SimpleHeaderLine(key, value)
    raise ValueError(f"Unsupported header line: {text}")


def _ensure_vcf_loaded(path: str):
    if path in _VCF_STORAGE:
        return

    header_lines = []
    columns_line = None
    records = []
    format_fields = {}

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                header_line = _header_from_text(line)
                header_lines.append(header_line)
                if isinstance(header_line, _FormatHeaderLine) and header_line.id:
                    format_fields[header_line.id] = header_line
                continue
            if line.startswith("#CHROM"):
                columns_line = line
                continue

            fields = line.split("\t")
            info_entries = OrderedDict()
            for part in fields[7].split(";"):
                if "=" not in part:
                    continue
                key, value = part.split("=", 1)
                info_entries[key] = value
            format_keys = fields[8].split(":") if len(fields) > 8 else []
            calls = []
            for idx, sample_value in enumerate(fields[9:]):
                sample_name = columns_line.split("\t")[9 + idx]
                values = sample_value.split(":") if format_keys else []
                data = {fmt: val for fmt, val in zip(format_keys, values)}
                calls.append(_FakeCall(sample_name, data))
            record = _FakeRecord(
                chrom=fields[0],
                pos=fields[1],
                record_id=fields[2],
                ref=fields[3],
                alts=fields[4].split(",") if fields[4] != "." else [],
                qual=fields[5],
                filt=fields[6],
                info=info_entries,
                format_keys=format_keys,
                calls=calls,
            )
            records.append(record)

    if columns_line is None:
        raise ValueError(f"VCF file {path} is missing #CHROM header")

    _VCF_STORAGE[path] = {
        "header_lines": header_lines,
        "columns_line": columns_line,
        "records": records,
        "format_fields": format_fields,
        "samples": columns_line.split("\t")[9:],
    }


def _register_vcf(path, header_lines, columns_line, records, format_fields, *, write_to_disk=True):
    _VCF_STORAGE[path] = {
        "header_lines": [line.copy() for line in header_lines],
        "columns_line": columns_line,
        "records": [record.copy() for record in records],
        "format_fields": {key: value.copy() for key, value in format_fields.items()},
        "samples": columns_line.split("\t")[9:],
    }

    if not write_to_disk:
        return

    with open(path, "w", encoding="utf-8") as handle:
        for line in header_lines:
            handle.write(line.to_line() + "\n")
        handle.write(columns_line + "\n")
        for record in records:
            handle.write(record.to_line() + "\n")


class _FakeReader:
    def __init__(self, path):
        _ensure_vcf_loaded(path)
        data = _VCF_STORAGE[path]
        self.header = _FakeHeader(
            data["header_lines"],
            data["samples"],
            data["columns_line"],
            data["format_fields"],
        )
        self._records = [record.copy() for record in data["records"]]
        self._index = 0

    @classmethod
    def from_path(cls, path, *args, **kwargs):
        return cls(path)

    @classmethod
    def from_stream(cls, *args, **kwargs):  # pragma: no cover - defensive
        raise NotImplementedError("Stream-based reading is not implemented in the stub")

    def __iter__(self):
        self._index = 0
        return self

    def __next__(self):
        if self._index >= len(self._records):
            raise StopIteration
        record = self._records[self._index].copy()
        self._index += 1
        return record

    def close(self):
        return None


class _FakeWriter:
    def __init__(self, path, header):
        self._path = path
        self._header = header
        self._records = []
        self._handle = open(path, "w", encoding="utf-8")
        for line in header.lines:
            self._handle.write(line.to_line() + "\n")
        self._handle.write(header.to_columns_line() + "\n")

    @classmethod
    def from_path(cls, path, header, *args, **kwargs):
        return cls(path, header)

    def write_record(self, record):
        self._records.append(record.copy())
        self._handle.write(record.to_line() + "\n")

    def close(self):
        self._handle.close()
        _register_vcf(
            self._path,
            [line.copy() for line in self._header.lines],
            self._header.to_columns_line(),
            self._records,
            self._header._format_fields,
        )


class _FakeVcfpyModule:
    Reader = _FakeReader
    Writer = _FakeWriter
    Call = _FakeCall
    SimpleHeaderLine = _SimpleHeaderLine
    SampleHeaderLine = _SampleHeaderLine
    header = SimpleNamespace(
        InfoHeaderLine=_InfoHeaderLine,
        FilterHeaderLine=_HeaderDefinition,
        ContigHeaderLine=_HeaderDefinition,
        SampleHeaderLine=_SampleHeaderLine,
        FormatHeaderLine=_FormatHeaderLine,
    )


def _generate_demo_vcf(path: Path, sample_name: str, genotype: str) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
        '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">',
        '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}".format(
            sample=sample_name
        ),
        "1\t100\t.\tA\tG\t.\tPASS\tAC=.;AN=.;AF=.\tGT\t{gt}".format(gt=genotype),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _install_vcfpy_stub(monkeypatch, module):
    stub = _FakeVcfpyModule()
    merging_module = sys.modules[module.merge_vcfs.__module__]
    validation_module = sys.modules[module.validate_all_vcfs.__module__]

    monkeypatch.setattr(module, "vcfpy", stub, raising=False)
    monkeypatch.setattr(merging_module, "vcfpy", stub, raising=False)
    monkeypatch.setattr(merging_module, "VCFPY_AVAILABLE", True, raising=False)
    monkeypatch.setattr(validation_module, "vcfpy", stub, raising=False)
    monkeypatch.setattr(validation_module, "VCFPY_AVAILABLE", True, raising=False)


def _configure_fake_subprocess(monkeypatch, module):
    merging_module = sys.modules[module.merge_vcfs.__module__]
    real_shutil = merging_module.shutil
    original_move = real_shutil.move

    def fake_move(src, dst, *args, **kwargs):
        result = original_move(src, dst, *args, **kwargs)
        if src in _VCF_STORAGE:
            _VCF_STORAGE[dst] = _VCF_STORAGE.pop(src)
        return result

    monkeypatch.setattr(real_shutil, "move", fake_move)

    def handle_merge(cmd):
        output_path = cmd[cmd.index("-o") + 1]
        input_files = cmd[cmd.index("-o") + 2 :]

        header_lines = [
            _SimpleHeaderLine("fileformat", "VCFv4.2"),
            _SimpleHeaderLine("reference", "GRCh38"),
            _InfoHeaderLine(
                {
                    "ID": "AC",
                    "Number": "A",
                    "Type": "Integer",
                    "Description": "Allele count in genotypes",
                }
            ),
            _InfoHeaderLine(
                {
                    "ID": "AN",
                    "Number": "1",
                    "Type": "Integer",
                    "Description": "Total number of alleles",
                }
            ),
            _InfoHeaderLine(
                {
                    "ID": "AF",
                    "Number": "A",
                    "Type": "Float",
                    "Description": "Allele frequency",
                }
            ),
            _FormatHeaderLine(
                {
                    "ID": "GT",
                    "Number": "1",
                    "Type": "String",
                    "Description": "Genotype",
                }
            ),
        ]

        sample_names = []
        calls = []
        for file_path in input_files:
            name = Path(file_path).stem
            with open(file_path, "r", encoding="utf-8") as handle:
                for raw_line in handle:
                    line = raw_line.strip()
                    if not line or line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        parts = line.split("\t")
                        name = parts[-1]
                        continue
                    fields = line.split("\t")
                    format_keys = fields[8].split(":") if len(fields) > 8 else []
                    sample_values = fields[9].split(":") if len(fields) > 9 else []
                    data = {fmt: val for fmt, val in zip(format_keys, sample_values)}
                    calls.append(_FakeCall(name, data))
            sample_names.append(name)

        record = _FakeRecord(
            chrom="1",
            pos="100",
            record_id=".",
            ref="A",
            alts=["G"],
            qual=".",
            filt="PASS",
            info=OrderedDict([("AC", "."), ("AN", "."), ("AF", ".")]),
            format_keys=["GT"],
            calls=calls,
        )

        columns_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(
            sample_names
        )
        format_fields = {"GT": header_lines[-1]}
        _register_vcf(output_path, header_lines, columns_line, [record], format_fields)

    def handle_fill_tags(cmd):
        input_path = cmd[2]
        output_path = cmd[cmd.index("-o") + 1]
        data = _VCF_STORAGE[input_path]
        records = []
        for record in data["records"]:
            updated = record.copy()
            alt_counts = [0 for _ in updated.ALT]
            allele_total = 0
            for call in updated.calls:
                genotype = call.data.get("GT", "")
                tokens = genotype.replace("|", "/").split("/")
                for token in tokens:
                    if not token or token == ".":
                        continue
                    try:
                        allele_index = int(token)
                    except ValueError:
                        continue
                    if allele_index < 0:
                        continue
                    allele_total += 1
                    if allele_index == 0:
                        continue
                    normalized = allele_index - 1
                    if 0 <= normalized < len(alt_counts):
                        alt_counts[normalized] += 1

            updated.INFO["AN"] = allele_total
            updated.INFO["AC"] = alt_counts
            if allele_total:
                updated.INFO["AF"] = [count / allele_total for count in alt_counts]
            else:
                updated.INFO["AF"] = [0.0 for _ in alt_counts]
            records.append(updated)

        _register_vcf(
            output_path,
            [line.copy() for line in data["header_lines"]],
            data["columns_line"],
            records,
            data["format_fields"],
        )

    def handle_bgzip(cmd):
        target = cmd[-1]
        gz_path = target + ".gz"
        with open(target, "rb") as src, gzip.open(gz_path, "wb") as dest:
            dest.write(src.read())
        os.remove(target)
        if target in _VCF_STORAGE:
            _VCF_STORAGE[gz_path] = _VCF_STORAGE.pop(target)

    def handle_tabix(cmd):
        target = cmd[-1]
        index_path = target + ".tbi"
        with open(index_path, "wb") as handle:
            handle.write(b"")

    def fake_run(cmd, *args, **kwargs):
        if cmd[:2] == ["bcftools", "merge"]:
            handle_merge(cmd)
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:2] == ["bcftools", "+fill-tags"]:
            handle_fill_tags(cmd)
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:1] == ["bgzip"]:
            handle_bgzip(cmd)
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:1] == ["tabix"]:
            handle_tabix(cmd)
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        raise AssertionError(f"Unexpected command: {cmd}")

    monkeypatch.setattr(merging_module, "subprocess", module.subprocess, raising=False)
    monkeypatch.setattr(module.subprocess, "run", fake_run)


def test_large_cohort_merge_workflow(tmp_path, monkeypatch):
    module = load_user_module()
    _install_vcfpy_stub(monkeypatch, module)

    _VCF_STORAGE.clear()

    input_dir = tmp_path / "inputs"
    output_dir = tmp_path / "outputs"
    input_dir.mkdir()
    output_dir.mkdir()

    sample_count = 50
    expected_alt = 0
    sample_files = []

    for index in range(sample_count):
        sample_name = f"SAMPLE_{index + 1:04d}"
        genotype = "0/1" if index % 2 == 0 else "1/1"
        if genotype == "0/1":
            expected_alt += 1
        else:
            expected_alt += 2
        path = input_dir / f"{sample_name}.vcf"
        _generate_demo_vcf(path, sample_name, genotype)
        sample_files.append(str(path))

    _configure_fake_subprocess(monkeypatch, module)

    final_path_holder = {}

    def fake_append_metadata(merged_vcf, *args, **kwargs):
        final_path_holder["path"] = merged_vcf
        return merged_vcf

    monkeypatch.setattr(module, "append_metadata_to_merged_vcf", fake_append_metadata)

    original_parse_metadata = module.parse_metadata_arguments

    def patched_parse_metadata(args, verbose=False, **kwargs):
        return original_parse_metadata(args, verbose=verbose)

    monkeypatch.setattr(module, "parse_metadata_arguments", patched_parse_metadata)

    args = SimpleNamespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        ref="GRCh38",
        vcf_version="4.2",
        meta=[],
        sample_metadata_entries=[],
        header_metadata_lines=[],
        verbose=False,
        allow_gvcf=False,
    )

    monkeypatch.setattr(module, "parse_arguments", lambda: args)

    summary_path = module.run_workflow(args)

    assert summary_path
    final_vcf_path = Path(final_path_holder["path"])
    assert final_vcf_path.exists()

    with gzip.open(final_vcf_path, "rt", encoding="utf-8") as handle:
        contents = handle.read()

    assert "##fileformat=VCFv4.2" in contents
    assert "##reference=GRCh38" in contents
    assert "##INFO=<ID=AC" in contents
    assert "##INFO=<ID=AN" in contents
    assert "##INFO=<ID=AF" in contents

    raw_lines = [line.strip() for line in contents.splitlines() if line.strip()]
    data_line = next(line for line in raw_lines if not line.startswith("#"))
    info_field = data_line.split("\t")[7]
    info_entries = {}
    for part in info_field.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        info_entries[key] = value

    assert info_entries["AC"] not in {".", ""}
    assert info_entries["AN"] not in {".", ""}
    assert info_entries["AF"] not in {".", ""}

    observed_ac = [int(val) for val in info_entries["AC"].split(",")]
    observed_an = int(info_entries["AN"])
    observed_af = [float(val) for val in info_entries["AF"].split(",")]

    assert observed_an == sample_count * 2
    assert observed_ac == [expected_alt]
    assert pytest.approx(observed_af[0], rel=1e-6) == expected_alt / (sample_count * 2)

