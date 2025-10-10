"""Unit tests covering merging helpers and anonymization workflows."""

from __future__ import annotations

import copy
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest
from _pytest.monkeypatch import MonkeyPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from MetagapUserCode.merge_vcf import merging
from MetagapUserCode.merge_vcf import metadata as metadata_module
from MetagapUserCode.tests.test_append_metadata_cli_smoke import (
    SAMPLE_BODY_TRIMMED,
    _configure_fake_bcftools,
)
from MetagapUserCode.tests.test_large_merge_workflow import _FakeVcfpyModule, _VCF_STORAGE


class _CallData(dict):
    """Mapping that normalizes GT lookups to mimic vcfpy semantics."""

    def __setitem__(self, key, value):
        if key == "GT" and isinstance(value, list):
            value = value[0] if value else None
        super().__setitem__(key, value)

    def __getitem__(self, key):
        value = super().__getitem__(key)
        return value

    def get(self, key, default=None):
        if key not in self:
            return default
        return self[key]

    def copy(self):  # pragma: no cover - simple delegation
        return _CallData({key: copy.deepcopy(value) for key, value in self.items()})


class _StubCall:
    """Minimal stand-in for ``vcfpy.Call`` used in unit tests."""

    def __init__(self, sample: str, data: dict[str, object]):
        self.sample = sample
        self.name = sample
        self.data = _CallData(data)

    def copy(self) -> "_StubCall":
        return _StubCall(self.sample, self.data.copy())

    def __deepcopy__(self, memo) -> "_StubCall":  # pragma: no cover - delegation
        return self.copy()


class _StubRecord:
    """Simple record carrying FORMAT/call data for padding tests."""

    def __init__(self, sample: str, genotype):
        self.CHROM = "1"
        self.POS = 100
        self.ID = "."
        self.REF = "A"
        self.ALT = ["G"]
        self.QUAL = "."
        self.FILTER = []
        self.INFO = {}
        self.FORMAT = ["GT"]
        self.calls = [_StubCall(sample, {"GT": genotype})]
        self.call_for_sample = {sample: self.calls[0]}

    def copy(self) -> "_StubRecord":
        duplicate = _StubRecord.__new__(_StubRecord)
        duplicate.CHROM = self.CHROM
        duplicate.POS = self.POS
        duplicate.ID = self.ID
        duplicate.REF = self.REF
        duplicate.ALT = list(self.ALT)
        duplicate.QUAL = self.QUAL
        duplicate.FILTER = list(self.FILTER)
        duplicate.INFO = dict(self.INFO)
        duplicate.FORMAT = list(self.FORMAT)
        duplicate.calls = [call.copy() for call in self.calls]
        duplicate.call_for_sample = {call.sample: call for call in duplicate.calls}
        return duplicate

    def __deepcopy__(self, memo) -> "_StubRecord":  # pragma: no cover - delegation
        return self.copy()

    def update_calls(self, updated_calls):
        self.calls = [call for call in updated_calls]
        self.call_for_sample = {call.sample: call for call in self.calls}


def _header_with_formats(sample_names):
    header = SimpleNamespace()
    header.samples = SimpleNamespace(names=list(sample_names))
    header.lines = []
    header.formats = {}

    def _get_format_field_info(key):
        if key == "GT":
            return SimpleNamespace(number="1")
        return None

    header.get_format_field_info = _get_format_field_info
    header.copy = lambda: _header_with_formats(header.samples.names)
    return header


def test_merge_colliding_records_recomputes_allele_metrics(monkeypatch):
    header = _header_with_formats(["S1", "S2"])
    sample_order = ["S1", "S2"]

    monkeypatch.setattr(merging.vcfpy, "Call", _StubCall)

    record_a = _StubRecord("S1", "0/1")
    record_a.INFO.update({"AC": [99], "AN": 99, "AF": [0.99]})

    record_b = _StubRecord("S2", "1/1")
    record_b.INFO.update({"AC": [1], "AN": 2, "AF": [0.5]})

    merged = merging._merge_colliding_records(
        [(record_a, 0), (record_b, 1)], header, sample_order
    )

    merging._recompute_ac_an_af(merged)

    assert merged.INFO["AC"] == [3]
    assert merged.INFO["AN"] == 4
    assert merged.INFO["AF"] == [pytest.approx(0.75)]

    zero_alt = merged.copy()
    zero_alt.INFO.update({"AC": [7], "AN": 8, "AF": [1.0]})
    for call in zero_alt.calls:
        call.data["GT"] = "./."

    merging._recompute_ac_an_af(zero_alt)

    assert zero_alt.INFO["AC"] == [0]
    assert zero_alt.INFO["AN"] == 0
    assert zero_alt.INFO["AF"] == [0.0]


def test_pad_record_samples_adds_missing_calls(monkeypatch):
    header = _header_with_formats(["S1", "S2"])

    record = _StubRecord("S1", "0/1")
    record.FORMAT.append("DP")
    record.call_for_sample["S1"].data["DP"] = 5

    def _format_info(key):
        if key == "DP":
            return SimpleNamespace(number="1")
        if key == "GT":
            return SimpleNamespace(number="1")
        return None

    header.get_format_field_info = _format_info

    merging._pad_record_samples(record, header, ["S1", "S2"])

    assert set(record.call_for_sample) == {"S1", "S2"}
    assert record.call_for_sample["S2"].data["GT"] == "./."
    assert record.call_for_sample["S2"].data["DP"] is None

    # Existing numeric scalars are normalized to single-element lists
    assert record.call_for_sample["S1"].data["DP"] == [5]


def test_pad_record_samples_normalizes_blank_genotypes():
    header = _header_with_formats(["S1"])
    record = _StubRecord("S1", ".")

    merging._pad_record_samples(record, header, ["S1"])

    assert record.call_for_sample["S1"].data["GT"] == "./."


def test_merge_vcfs_pads_missing_samples(monkeypatch, tmp_path):
    sample_order = ["S1", "S2"]

    def fake_union_headers(paths, sample_order=None):
        names = sample_order if sample_order is not None else ["S1"]
        return _header_with_formats(names)

    monkeypatch.setattr(merging, "union_headers", fake_union_headers)
    monkeypatch.setattr(merging, "apply_metadata_to_header", lambda header, **_: header)
    monkeypatch.setattr(merging, "preprocess_vcf", lambda path: path)
    monkeypatch.setattr(merging, "log_message", lambda *args, **kwargs: None)
    monkeypatch.setattr(merging, "_filter_vcf_records", lambda *args, **kwargs: None)
    monkeypatch.setattr(merging, "_ensure_info_header_lines", lambda header: None)

    reader_instances = {}

    class _Reader:
        def __init__(self, path):
            self.header = _header_with_formats(["S1"])
            self._records = [_StubRecord("S1", "0/1")]

        def __iter__(self):
            for record in self._records:
                yield record.copy()

        def close(self):
            return None

    class _Writer:
        def __init__(self, path, header):
            self.path = Path(path)
            self.header = header
            self.records = []
            self._handle = self.path.open("w", encoding="utf-8")

        def write_record(self, record):
            self.records.append(record.copy())
            self._handle.write("record\n")

        def close(self):
            self._handle.close()
            reader_instances["writer"] = self

    def fake_reader_from_path(path):
        instance = _Reader(path)
        return instance

    class _ReaderFactory:
        @staticmethod
        def from_path(path):
            return fake_reader_from_path(path)

    class _WriterFactory:
        @staticmethod
        def from_path(path, header):
            writer = _Writer(path, header)
            reader_instances["writer"] = writer
            return writer

    monkeypatch.setattr(merging.vcfpy, "Call", _StubCall)
    monkeypatch.setattr(merging.vcfpy, "Reader", _ReaderFactory)
    monkeypatch.setattr(merging.vcfpy, "Writer", _WriterFactory)
    monkeypatch.setattr(
        merging.vcfpy,
        "header",
        SimpleNamespace(FormatHeaderLine=object, InfoHeaderLine=object),
    )

    def fake_run(cmd, *args, **kwargs):
        if cmd[:2] == ["bcftools", "merge"]:
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:2] == ["bcftools", "+fill-tags"]:
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:2] == ["bcftools", "view"]:
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:1] == ["bgzip"]:
            target = Path(cmd[-1])
            gz_path = Path(str(target) + ".gz")
            gz_path.write_text("stub")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:1] == ["tabix"]:
            Path(cmd[-1] + ".tbi").write_text("stub-index")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(merging, "subprocess", SimpleNamespace(run=fake_run))

    output_dir = tmp_path
    input_path = output_dir / "sample1.vcf"
    input_path.write_text("placeholder")

    result_path = merging.merge_vcfs([str(input_path)], str(output_dir), sample_order=sample_order)

    writer = reader_instances["writer"]
    assert writer.records, "Expected merge to emit at least one record"
    merged_record = writer.records[0]

    assert set(merged_record.call_for_sample) == {"S1", "S2"}
    assert merged_record.call_for_sample["S2"].data["GT"] == "./."
    assert merged_record.call_for_sample["S1"].data["GT"] == "0/1"
    assert result_path.endswith(".gz")


def test_append_metadata_to_merged_vcf_strips_sample_columns(tmp_path):
    header_checks: list[Path] = []
    patcher = MonkeyPatch()
    _VCF_STORAGE.clear()
    fake_vcfpy = _FakeVcfpyModule()
    patcher.setattr(metadata_module, "vcfpy", fake_vcfpy, raising=False)
    patcher.setattr(metadata_module, "VCFPY_AVAILABLE", True, raising=False)

    def _writer_enter(self):
        return self

    def _writer_exit(self, exc_type, exc, tb):
        self.close()
        return False

    patcher.setattr(fake_vcfpy.Writer, "__enter__", _writer_enter, raising=False)
    patcher.setattr(fake_vcfpy.Writer, "__exit__", _writer_exit, raising=False)

    original_writer_from_path = fake_vcfpy.Writer.from_path.__func__

    def _patched_writer_from_path(cls, path, header, *args, **kwargs):
        base_columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
        if (
            hasattr(header, "samples")
            and hasattr(header.samples, "names")
            and not header.samples.names
            and hasattr(header, "_columns_line")
        ):
            header._columns_line = base_columns
        return original_writer_from_path(cls, path, header, *args, **kwargs)

    patcher.setattr(
        fake_vcfpy.Writer,
        "from_path",
        classmethod(_patched_writer_from_path),
        raising=False,
    )

    class _Record:
        def __init__(
            self,
            CHROM,
            POS,
            ID,
            REF,
            ALT,
            QUAL,
            FILTER,
            INFO,
            FORMAT,
            calls,
        ):
            self.CHROM = CHROM
            self.POS = POS
            self.ID = ID
            self.REF = REF
            if isinstance(ALT, (list, tuple)):
                self.ALT = list(ALT)
            elif ALT in (None, ".", ""):
                self.ALT = []
                
            else:
                self.ALT = [ALT]
            self.QUAL = QUAL
            if isinstance(FILTER, (list, tuple)):
                self.FILTER = list(FILTER)
            elif FILTER in (None, "", "."):
                self.FILTER = []
            else:
                self.FILTER = [FILTER]
            self.INFO = INFO
            self.FORMAT = FORMAT
            self.calls = calls

        def copy(self):
            return _Record(
                self.CHROM,
                self.POS,
                self.ID,
                self.REF,
                list(self.ALT),
                self.QUAL,
                list(self.FILTER),
                copy.deepcopy(self.INFO),
                self.FORMAT,
                [call for call in self.calls],
            )

        def to_line(self):
            info_parts = []
            for key, value in (self.INFO or {}).items():
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
                str(self.QUAL) if self.QUAL not in {None, ""} else ".",
                ";".join(self.FILTER) if self.FILTER else ".",
                info_field,
            ]

            if self.FORMAT:
                fields.append(":".join(self.FORMAT))
                for call in self.calls:
                    if hasattr(call, "to_string"):
                        fields.append(call.to_string(self.FORMAT))
                    else:
                        fields.append(".")

            return "\t".join(fields)

    patcher.setattr(fake_vcfpy, "Record", _Record, raising=False)

    for cls_name in ("InfoHeaderLine", "SampleHeaderLine"):
        cls = getattr(fake_vcfpy.header, cls_name)

        def _from_mapping(inner_cls, mapping, _cls=cls):
            return _cls(mapping)

        patcher.setattr(cls, "from_mapping", classmethod(_from_mapping), raising=False)

    _configure_fake_bcftools(patcher, metadata_module, header_checks)

    merged_vcf = tmp_path / "merged.vcf"
    merged_vcf.write_text("placeholder")

    try:
        final_vcf = metadata_module.append_metadata_to_merged_vcf(str(merged_vcf))
    finally:
        patcher.undo()

    contents = Path(final_vcf).read_text()
    header_lines = [line for line in contents.splitlines() if line.startswith("##")]
    assert any(line.startswith("##fileformat=VCFv4.2") for line in header_lines)
    assert all(not line.startswith("##FORMAT=") for line in header_lines)

    columns_line = next(line for line in contents.splitlines() if line.startswith("#CHROM"))
    assert columns_line.split("\t") == [
        "#CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
    ]

    data_lines = [line for line in contents.splitlines() if line and not line.startswith("#")]
    assert data_lines == [SAMPLE_BODY_TRIMMED.strip()]
    assert all(len(line.split("\t")) == 8 for line in data_lines)
