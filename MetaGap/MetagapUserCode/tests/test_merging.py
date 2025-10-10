"""Unit tests covering merging helpers and anonymization workflows."""

from __future__ import annotations

import shutil
from typing import Callable
import datetime
import gzip
import logging
from collections import OrderedDict
from pathlib import Path
import time
import sys
from types import SimpleNamespace
from typing import Sequence

import pytest
from _pytest.monkeypatch import MonkeyPatch

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from MetagapUserCode.merge_vcf import merging
from MetagapUserCode.merge_vcf import metadata as metadata_module
from MetagapUserCode.tests.test_append_metadata_cli_smoke import (
    SAMPLE_BODY_TRIMMED,
    _configure_fake_bcftools,
)
from MetagapUserCode.tests.test_large_merge_workflow import (
    _FakeCall,
    _FakeRecord,
    _FakeVcfpyModule,
    _HeaderDefinition,
    _VCF_STORAGE,
    _generate_demo_vcf,
)


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

    def __init__(self, sample: str, genotype, *, alts=None, format_keys=None):
        self.CHROM = "1"
        self.POS = 100
        self.ID = "."
        self.REF = "A"
        self.ALT = list(alts) if alts is not None else ["G"]
        self.QUAL = "."
        self.FILTER = []
        self.INFO = {}
        self.FORMAT = list(format_keys) if format_keys is not None else ["GT"]
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


def _serialize_record(record) -> str:
    alt = ",".join(record.ALT) if getattr(record, "ALT", None) else "."
    filt = ";".join(record.FILTER) if getattr(record, "FILTER", None) else "."
    qual = getattr(record, "QUAL", None)
    if qual in {None, "", "."}:
        qual_field = "."
    else:
        qual_field = str(qual)

    fields = [
        str(record.CHROM),
        str(record.POS),
        getattr(record, "ID", ".") or ".",
        getattr(record, "REF", "."),
        alt if alt else ".",
        qual_field,
        filt if filt else ".",
        ".",
    ]

    fmt_keys = list(getattr(record, "FORMAT", []) or [])
    if fmt_keys:
        fields.append(":".join(fmt_keys))
        for call in getattr(record, "calls", []):
            values = []
            data = getattr(call, "data", {}) or {}
            for key in fmt_keys:
                value = data.get(key)
                if isinstance(value, list):
                    if not value:
                        values.append(".")
                    else:
                        serialized = ",".join(
                            "." if entry in {None, "", "."} else str(entry)
                            for entry in value
                        )
                        values.append(serialized)
                elif value in {None, "", "."}:
                    values.append(".")
                else:
                    values.append(str(value))
            fields.append(":".join(values))

    return "\t".join(fields)

class _MultiSampleRecord:
    """Record supporting multiple samples for _merge_colliding_records tests."""

    def __init__(self, format_keys, sample_data, *, chrom="1", pos=100, ref="A", alt=None):
        self.CHROM = chrom
        self.POS = pos
        self.ID = "."
        self.REF = ref
        self.ALT = list(alt or ["G"])
        self.QUAL = "."
        self.FILTER = []
        self.INFO = {}
        self.FORMAT = list(format_keys)
        self.calls = [_StubCall(sample, dict(data)) for sample, data in sample_data.items()]
        self.call_for_sample = {call.sample: call for call in self.calls}

    def copy(self) -> "_MultiSampleRecord":
        duplicate = _MultiSampleRecord.__new__(_MultiSampleRecord)
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

    def __deepcopy__(self, memo) -> "_MultiSampleRecord":  # pragma: no cover - delegation
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


def test_merge_colliding_records_recomputes_ac_an_af_from_calls():
    header = _header_with_formats(["S1", "S2"])
    r1 = _StubRecord("S1", "0/1")
    r1.INFO.update({"AC": [5], "AN": 8, "AF": [0.6]})
    r2 = _StubRecord("S2", "1/1")
    r2.INFO.update({"AC": [2], "AN": 4, "AF": [0.5]})

    merged = merging._merge_colliding_records([(r1, 0), (r2, 1)], header, ["S1", "S2"])
    merging._recompute_ac_an_af(merged)

    assert merged.INFO["AC"] == [3]
    assert merged.INFO["AN"] == 4
    assert merged.INFO["AF"] == pytest.approx([0.75])


def test_recompute_ac_an_af_defaults_zero_frequencies_when_an_missing():
    record = _StubRecord("S1", "./.")
    record.INFO.update({"AC": [7], "AN": 10, "AF": [0.7]})

    merging._recompute_ac_an_af(record)

    assert record.INFO["AC"] == [0]
    assert record.INFO["AN"] == 0
    assert record.INFO["AF"] == [0.0]
def test_merge_colliding_records_inserts_missing_sample_defaults(monkeypatch):
    header = _header_with_formats(["X"])

    def _format_info(key):
        mappings = {
            "GT": SimpleNamespace(number="1"),
            "DP": SimpleNamespace(number="1"),
            "GQ": SimpleNamespace(number="1"),
            "AD": SimpleNamespace(number="R"),
        }
        return mappings.get(key)

    header.get_format_field_info = _format_info

    rec1 = _StubRecord("X", "0/1")
    rec1.FORMAT.extend(["DP", "GQ"])
    rec1.call_for_sample["X"].data["DP"] = 7
    rec1.call_for_sample["X"].data["GQ"] = 45

    rec2 = _StubRecord("X", "0/0")
    rec2.FORMAT.append("AD")
    rec2.call_for_sample["X"].data["AD"] = [10, 5]

    monkeypatch.setattr(merging.vcfpy, "Call", _StubCall)

    merged = merging._merge_colliding_records(
        [(rec1, 0), (rec2, 1)],
def test_merge_colliding_records_pads_missing_samples():
    header = _header_with_formats(["X", "Y"])

    def _format_info(key):
        if key == "GT":
            return SimpleNamespace(number="1")
        if key == "DP":
            return SimpleNamespace(number="1")
        if key == "AD":
            return SimpleNamespace(number="R")
        return None

    header.get_format_field_info = _format_info

    left = _MultiSampleRecord(["GT", "DP"], {"X": {"GT": "0/1", "DP": 12}})
    right = _MultiSampleRecord(["GT", "AD"], {"X": {"GT": "0/0", "AD": [3, 5]}})

    merged = merging._merge_colliding_records(
        [(left, 0), (right, 1)],
        header,
        ["X", "Y"],
    )

    assert list(merged.call_for_sample) == ["X", "Y"]
    missing = merged.call_for_sample["Y"].data
    assert missing["GT"] == "./."
    assert missing["DP"] is None
    assert missing["GQ"] is None
    assert missing["AD"] == []


def test_merge_vcfs_pads_missing_samples(monkeypatch, tmp_path):
    sample_order = ["S1", "S2"]

    def _with_format_info(header):
        def _format_info(key):
            if key == "GT":
                return SimpleNamespace(number="1")
            if key == "DP":
                return SimpleNamespace(number="1")
            if key == "AD":
                return SimpleNamespace(number="R")
            return None

        header.get_format_field_info = _format_info
        return header

    def fake_union_headers(paths, sample_order=None):
        names = sample_order if sample_order is not None else ["S1"]
        return _with_format_info(_header_with_formats(names))

    monkeypatch.setattr(merging, "union_headers", fake_union_headers)
    monkeypatch.setattr(merging, "apply_metadata_to_header", lambda header, **_: header)
    monkeypatch.setattr(merging, "preprocess_vcf", lambda path: path)
    monkeypatch.setattr(merging, "log_message", lambda *args, **kwargs: None)
    monkeypatch.setattr(merging, "_filter_vcf_records", lambda *args, **kwargs: None)
    monkeypatch.setattr(merging, "_ensure_info_header_lines", lambda header: None)

    reader_registry: dict[str, Callable[[], object]] = {}

    def _register_reader(path: Path, factory):
        reader_registry[str(path)] = factory

    format_header = _with_format_info(_header_with_formats(sample_order))

    records_by_path: dict[str, list[_MultiSampleRecord]] = {}
    headers_by_path: dict[str, SimpleNamespace] = {}

    class _Reader:
        def __init__(self, header, records):
            self.header = header
            self._records = [rec.copy() for rec in records]

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

            def _factory():
                return _Reader(self.header.copy(), self.records)

            _register_reader(self.path, _factory)

            self._handle.write("##fileformat=VCFv4.2\n")
            columns = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
            sample_names = getattr(getattr(self.header, "samples", None), "names", []) or []
            if sample_names:
                columns += "\tFORMAT"
                for name in sample_names:
                    columns += f"\t{name}"
            self._handle.write(columns + "\n")

        def write_record(self, record):
            self.records.append(record.copy())
            self._handle.write(_serialize_record(record) + "\n")

        def close(self):
            if self._handle:
                self._handle.close()
                self._handle = None

            def _factory():
                return _Reader(self.header.copy(), self.records)

            _register_reader(self.path, _factory)

    class _ReaderFactory:
        @staticmethod
        def from_path(path):
            return fake_reader_from_path(path)
    class _WriterFactory:
        @staticmethod
        def from_path(path, header):
            writer = _Writer(path, header)
            if "merged_writer" not in reader_registry:
                reader_registry["merged_writer"] = lambda: writer
            return writer

    def _register_input(path: Path, records: list[_MultiSampleRecord], samples=None):
        sample_names = list(samples) if samples is not None else ["S1"]
        headers_by_path[str(path)] = _with_format_info(_header_with_formats(sample_names))
        records_by_path[str(path)] = [record.copy() for record in records]

    monkeypatch.setattr(merging.vcfpy, "Call", _StubCall)
    monkeypatch.setattr(merging.vcfpy, "Reader", _ReaderFactory)
    monkeypatch.setattr(merging.vcfpy, "Writer", _WriterFactory)
    monkeypatch.setattr(
        merging.vcfpy,
        "header",
        SimpleNamespace(FormatHeaderLine=object, InfoHeaderLine=object),
    )

    class _PysamStub:
        @staticmethod
        def tabix_compress(src, dest, force=True):
            with open(src, "rt", encoding="utf-8") as src_handle:
                with gzip.open(dest, "wt", encoding="utf-8") as dest_handle:
                    dest_handle.write(src_handle.read())

        @staticmethod
        def tabix_index(path, preset="vcf", force=True):
            Path(path + ".tbi").write_text("stub-index", encoding="utf-8")

    monkeypatch.setattr(merging, "pysam", _PysamStub())

    output_dir = tmp_path
    input_path_1 = output_dir / "sample1.vcf"
    input_path_2 = output_dir / "sample2.vcf"
    input_path_1.write_text("placeholder")
    input_path_2.write_text("placeholder")

    record_a = _MultiSampleRecord(["GT", "DP"], {"S1": {"GT": "0/1", "DP": 9}})
    record_b = _MultiSampleRecord(["GT", "AD"], {"S1": {"GT": "0/0", "AD": [4, 6]}})

    _register_input(input_path_1, [record_a], samples=["S1"])
    _register_input(input_path_2, [record_b], samples=["S1"])

    def _input_reader_factory():
        header = _header_with_formats(["S1"])
        record = _StubRecord("S1", "0/1")
        return _Reader(header, [record])

    _register_reader(input_path, _input_reader_factory)

    result_path = merging.merge_vcfs([str(input_path)], str(output_dir), sample_order=sample_order)

    writer = reader_registry["merged_writer"]()
    assert writer.records, "Expected merge to emit at least one record"
    merged_record = writer.records[0]

    assert set(merged_record.call_for_sample) == {"S1", "S2"}
    s1_call = merged_record.call_for_sample["S1"].data
    s2_call = merged_record.call_for_sample["S2"].data

    assert s2_call["GT"] == "./."
    assert s2_call["DP"] is None
    assert s2_call["AD"] == []
    assert s1_call["GT"] == "0/0"
    assert s1_call["DP"] is None
    assert s1_call["AD"] == [4, 6]
    assert merged_record.FORMAT == ["GT", "DP", "AD"]
    assert result_path.endswith(".gz")

    gz_path = Path(result_path)
    assert gz_path.exists()
    with gzip.open(gz_path, "rt", encoding="utf-8") as handle:
        lines = [line.strip() for line in handle.readlines() if line.strip()]

    assert any(line.startswith("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO") for line in lines)
    data_rows = [line for line in lines if not line.startswith("#")]
    assert data_rows, "Expected data rows in compressed output"
    assert all(len(row.split("\t")) == 8 for row in data_rows)


def test_filter_vcf_records_applies_thresholds_and_filters(monkeypatch, tmp_path, caplog):
    class _Variant:
        def __init__(self, qual, an, filters):
            self.QUAL = qual
            self.INFO = {"AN": an} if an is not None else {}
            self.FILTER = list(filters)

    records = [
        _Variant(qual=60, an=120, filters=[]),
        _Variant(qual=10, an=200, filters=["PASS"]),
        _Variant(qual=80, an=20, filters=["PASS"]),
        _Variant(qual=70, an=200, filters=["q10"]),
        _Variant(qual=90, an=200, filters=["q10", "LowQual"]),
    ]

    input_path = tmp_path / "variants.vcf"
    input_path.write_text("placeholder")

    reader_instances = []
    writer_instances = []

    class _FakeReader:
        def __init__(self, path):
            assert str(path) == str(input_path)
            self.header = SimpleNamespace()
            self.closed = False

        def __iter__(self):
            for record in records:
                yield record

        def close(self):
            self.closed = True

    class _ReaderFactory:
        @staticmethod
        def from_path(path):
            reader = _FakeReader(path)
            reader_instances.append(reader)
            return reader

    class _FakeWriter:
        def __init__(self, path, header):
            assert path.endswith(".filtered")
            self.path = Path(path)
            self.header = header
            self.records = []
            self.closed = False
            self._handle = self.path.open("w", encoding="utf-8")

        def write_record(self, record):
            self.records.append(record)
            self._handle.write("record\n")

        def close(self):
            if not self.closed:
                self._handle.close()
                self.closed = True

    class _WriterFactory:
        @staticmethod
        def from_path(path, header):
            writer = _FakeWriter(path, header)
            writer_instances.append(writer)
            return writer

    monkeypatch.setattr(merging.vcfpy, "Reader", _ReaderFactory)
    monkeypatch.setattr(merging.vcfpy, "Writer", _WriterFactory)

    with caplog.at_level(logging.INFO, logger="vcf_merger"):
        merging._filter_vcf_records(
            str(input_path),
            qual_threshold=20,
            an_threshold=50,
            allowed_filter_values=("PASS", "q10"),
            verbose=False,
        )

    reader = reader_instances[0]
    writer = writer_instances[0]

    assert reader.closed
    assert writer.closed
    assert writer.records == [records[0], records[3]]

    final_lines = input_path.read_text().splitlines()
    assert len(final_lines) == 2

    summary_messages = [
        rec.message for rec in caplog.records if "Applied variant filter" in rec.message
    ]
    assert summary_messages == ["Applied variant filter: kept 2 of 5."]


def test_merge_colliding_records_rewrites_genotypes_for_reordered_alts():
    header = _header_with_formats(["S1", "S2", "S3"])

    record_a = _StubRecord("S1", "1|0", alts=["G", "T"])
    record_b = _StubRecord("S2", "1|0", alts=["T", "G"])
    record_c = _StubRecord("S3", "0/1", alts=["T", "G"])

    merged = merging._merge_colliding_records(
        [(record_a, 0), (record_b, 1), (record_c, 2)],
        header,
        ["S1", "S2", "S3"],
    )

    assert [alt for alt in merged.ALT] == ["G", "T"]
    assert merged.call_for_sample["S1"].data["GT"] == "1|0"
    assert merged.call_for_sample["S2"].data["GT"] == "2|0"
    assert merged.call_for_sample["S3"].data["GT"] == "0/2"


def test_filter_vcf_records_applies_thresholds(tmp_path, monkeypatch, caplog):
    class _Record:
        def __init__(self, qual, an, filters):
            self.QUAL = qual
            self.INFO = {} if an is None else {"AN": an}
            self.FILTER = filters

    records = [
        _Record(50, 120, ["PASS"]),
        _Record(35, 80, ["q10"]),
        _Record(10, 80, ["q10"]),
        _Record(40, 40, ["PASS"]),
        _Record(45, 90, ["LowQual"]),
        _Record(55, [150], []),
        _Record(60, None, ["PASS"]),
    ]

    written_records = []

    class _StubReader:
        def __init__(self, entries):
            self.header = SimpleNamespace()
            self._entries = entries

        def __iter__(self):
            for entry in self._entries:
                yield entry

        def close(self):
            return None

    class _StubWriter:
        def __init__(self, path, header):
            self.path = Path(path)
            self.header = header
            self._handle = self.path.open("w", encoding="utf-8")

        def write_record(self, record):
            written_records.append(record)
            self._handle.write("record\n")

        def close(self):
            self._handle.close()

    class _ReaderFactory:
        @staticmethod
        def from_path(path):
            return _StubReader(records)

    class _WriterFactory:
        @staticmethod
        def from_path(path, header):
            return _StubWriter(path, header)

    monkeypatch.setattr(merging.vcfpy, "Reader", _ReaderFactory)
    monkeypatch.setattr(merging.vcfpy, "Writer", _WriterFactory)

    input_path = tmp_path / "input.vcf"
    input_path.write_text("original\n", encoding="utf-8")

    with caplog.at_level(logging.INFO, logger="vcf_merger"):
        merging._filter_vcf_records(
            str(input_path),
            qual_threshold=20,
            an_threshold=50,
            allowed_filter_values=("PASS", "q10"),
            verbose=False,
        )

    assert written_records == [records[0], records[1], records[5]]
    assert input_path.read_text(encoding="utf-8").count("record\n") == 3

    messages = [rec.message for rec in caplog.records if "Applied variant filter" in rec.message]
    assert messages == ["Applied variant filter: kept 3 of 7."]


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

def test_append_metadata_cleanup_on_compression_failure(monkeypatch, tmp_path):
    import gzip
    import pytest
    from pathlib import Path
    from types import SimpleNamespace
    from collections import OrderedDict

    merged_vcf = tmp_path / "merged.vcf.gz"
    input_lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t10\t.\tA\tG\t99\tPASS\t.",
    ]
    with gzip.open(merged_vcf, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(input_lines) + "\n")

    record = SimpleNamespace(
        CHROM="1",
        POS=10,
        ID=".",
        REF="A",
        ALT=["G"],
        QUAL=99,
        FILTER=["PASS"],
        INFO=OrderedDict(),
        calls=[SimpleNamespace(data={"GT": "0/1"})],
    )

    class _StubReader:
        def __init__(self, entries):
            self._entries = entries
            self.header = SimpleNamespace()

        def __iter__(self):
            return iter(self._entries)

        def close(self):
            pass

    class _ReaderFactory:
        @classmethod
        def from_stream(cls, stream):
            return _StubReader([record])

    fake_vcfpy = SimpleNamespace(Reader=_ReaderFactory)
    monkeypatch.setattr(metadata_module, "vcfpy", fake_vcfpy, raising=False)
    monkeypatch.setattr(metadata_module, "VCFPY_AVAILABLE", True, raising=False)

    class _FailingPysam:
        @staticmethod
        def bgzip(*_args, **_kwargs):
            raise OSError("bgzip failed")

        @staticmethod
        def tabix_index(*_args, **_kwargs):
            raise AssertionError("tabix_index should not be called")

    monkeypatch.setattr(metadata_module, "pysam", _FailingPysam, raising=False)
    monkeypatch.setattr(metadata_module, "PYSAM_AVAILABLE", True, raising=False)
    monkeypatch.setattr(metadata_module, "log_message", lambda *a, **k: None, raising=False)

    joint_temp = Path(f"{merged_vcf}.joint.temp.vcf")
    filtered_temp = Path(f"{merged_vcf}.filtered.temp.vcf")
    joint_temp.write_text("joint", encoding="utf-8")
    filtered_temp.write_text("filtered", encoding="utf-8")

    captured = {}

    def _fake_handle(message, **_kwargs):
        captured["message"] = message
        raise RuntimeError("critical failure")

    monkeypatch.setattr(metadata_module, "handle_critical_error", _fake_handle, raising=False)

    with pytest.raises(RuntimeError):
        metadata_module.append_metadata_to_merged_vcf(str(merged_vcf))

    assert "Failed to compress or index anonymized VCF" in captured["message"]
    assert not joint_temp.exists()
    assert not filtered_temp.exists()

    final_plain = tmp_path / "merged.anonymized.vcf"
    final_gzip = Path(str(final_plain) + ".gz")
    assert not final_plain.exists()
    assert not final_gzip.exists()


def test_merge_vcfs_handles_many_shards_with_cleanup(tmp_path):
    import time
    from pathlib import Path

    output_dir = tmp_path / "output"
    output_dir.mkdir()
    shard_dir = tmp_path / "shards"
    shard_dir.mkdir()

    def _write_shard(path: Path, sample_name: str, genotype: str) -> None:
        ac_value = "2" if genotype == "1/1" else "1" if genotype == "0/1" else "0"
        info_field = f"AC={ac_value};AN=60;AF=0.05"
        lines = [
            "##fileformat=VCFv4.2",
            "##reference=GRCh38",
            '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
            '##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">',
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">',
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
            f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}",
            f"1\t100\t.\tA\tG\t60\tPASS\t{info_field}\tGT\t{genotype}",
        ]
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    shard_paths = []
    for index in range(55):
        genotype = "0/1" if index % 3 else "0/0"
        shard_path = shard_dir / f"sample_{index:02d}.vcf"
        _write_shard(shard_path, f"Sample{index:02d}", genotype)
        shard_paths.append(str(shard_path))

    try:
        import psutil
        process = psutil.Process()
        start_mem = process.memory_info().rss
    except Exception:  # pragma: no cover
        process = None
        start_mem = None

    start_time = time.monotonic()
    merged_path = merging.merge_vcfs(shard_paths, str(output_dir), verbose=True)
    duration = time.monotonic() - start_time

    assert duration < 20
    if process is not None and start_mem is not None:
        end_mem = process.memory_info().rss
        assert end_mem - start_mem < 128 * 1024 * 1024

    merged_file = Path(merged_path)
    assert merged_file.exists()
    index_file = Path(str(merged_file) + ".tbi")
    assert index_file.exists()

    remaining = {entry.name for entry in output_dir.iterdir()}
    assert remaining == {merged_file.name, index_file.name}


def test_merge_vcfs_emits_anonymized_gzip(tmp_path, monkeypatch):
    import gzip, shutil
    from pathlib import Path
    from types import SimpleNamespace

    shard_template = """##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Allele number">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
1\t{pos}\t.\tA\tC\t99\tPASS\tAC=1;AN=2;AF=0.5\tGT\t0/1
"""

    shard_paths = []
    for idx, pos in enumerate((100, 200), start=1):
        p = tmp_path / f"shard{idx}.vcf"
        p.write_text(shard_template.format(sample=f"S{idx}", pos=pos), encoding="utf-8")
        shard_paths.append(str(p))

    compressed = {}

    def fake_tabix_compress(src, dest, force=True):
        with open(src, "rb") as s, gzip.open(dest, "wb") as d:
            shutil.copyfileobj(s, d)
        compressed["source"] = Path(src)
        compressed["dest"] = Path(dest)

    def fake_tabix_index(path, preset="vcf", force=True):
        idx = Path(f"{path}.tbi}")
        idx.write_text("stub-index", encoding="utf-8")
        compressed["index"] = idx

    monkeypatch.setattr(
        merging,
        "pysam",
        SimpleNamespace(tabix_compress=fake_tabix_compress, tabix_index=fake_tabix_index),
    )

    result_path = Path(
        merging.merge_vcfs(
            shard_paths,
            str(tmp_path),
            qual_threshold=None,
            an_threshold=None,
        )
    )

    assert result_path.suffixes[-2:] == [".vcf", ".gz"]
    assert "dest" in compressed

    with gzip.open(result_path, "rt", encoding="utf-8") as fh:
        lines = [ln.rstrip("\n") for ln in fh if ln.strip()]

    header_lines = [ln for ln in lines if ln.startswith("##")]
    assert header_lines
    assert all(not ln.startswith("##FORMAT=") for ln in header_lines)

    columns_line = next(ln for ln in lines if ln.startswith("#CHROM"))
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

    data_lines = [ln for ln in lines if not ln.startswith("#")]
    assert data_lines
    assert all(len(ln.split("\t")) == 8 for ln in data_lines)

    base_vcf_path = Path(str(result_path)[:-3])
    assert not base_vcf_path.exists()
