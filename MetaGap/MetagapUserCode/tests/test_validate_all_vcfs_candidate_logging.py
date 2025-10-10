"""Tests covering logging for candidate VCF discovery in ``validate_all_vcfs``."""

from __future__ import annotations

import gzip
import sys
from types import SimpleNamespace
from pathlib import Path

import importlib
import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

validation = importlib.import_module("merge_vcf.validation")

_VCF_TEXT = "\n".join(
    [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##contig=<ID=1,length=1000>",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        "1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
    ]
) + "\n"


class _StubHeaderLine:
    def __init__(self, key: str, mapping: dict[str, str] | None = None, value: str | None = None):
        self.key = key
        self.value = value
        self.mapping = mapping or {}
        self.id = self.mapping.get("ID")


class _StubHeader:
    def __init__(self) -> None:
        self.samples = SimpleNamespace(names=["SAMPLE"])
        self.lines = [
            SimpleNamespace(key="fileformat", value="VCFv4.2"),
            SimpleNamespace(key="reference", value="GRCh38"),
            _StubHeaderLine("contig", {"ID": "1", "length": "1000"}),
            _StubHeaderLine("FILTER", {"ID": "PASS", "Description": "All filters passed"}),
            _StubHeaderLine(
                "INFO",
                {
                    "ID": "DP",
                    "Number": "1",
                    "Type": "Integer",
                    "Description": "Total Depth",
                },
            ),
            _StubHeaderLine(
                "FORMAT",
                {
                    "ID": "GT",
                    "Number": "1",
                    "Type": "String",
                    "Description": "Genotype",
                },
            ),
        ]

    def get_lines(self, key: str):
        key_upper = key.upper()
        return [line for line in self.lines if getattr(line, "key", "").upper() == key_upper]


class _StubReader:
    def __init__(self, path: str) -> None:
        self.path = path
        self.header = _StubHeader()
        self._records = iter(
            [
                SimpleNamespace(CHROM="1", POS=10),
                SimpleNamespace(CHROM="1", POS=20),
            ]
        )

    def __iter__(self):
        return self

    def __next__(self):
        return next(self._records)

    def close(self) -> None:
        pass


def _install_reader_stub(monkeypatch, module) -> None:
    def _fake_from_path(cls, path, *_, **__):
        return _StubReader(path)

    monkeypatch.setattr(module.vcfpy.Reader, "from_path", classmethod(_fake_from_path))


def _write_vcf(path: Path) -> None:
    if path.suffix == ".gz":
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(_VCF_TEXT)
    else:
        path.write_text(_VCF_TEXT)


def _prepare_successful_validation(monkeypatch, module, prepared_inputs):
    monkeypatch.setattr(module, "discover_and_prepare_inputs", lambda *_, **__: prepared_inputs)
    monkeypatch.setattr(module, "preprocess_vcf", lambda file_path, **__: str(file_path))
    monkeypatch.setattr(module, "validate_vcf", lambda *_, **__: True)
    _install_reader_stub(monkeypatch, module)


@pytest.mark.skipif(not getattr(validation, "VCFPY_AVAILABLE", False), reason="vcfpy dependency is required for validation tests")
def test_validate_all_vcfs_logs_counts_with_only_gz(tmp_path, monkeypatch):
    module = validation

    shard_path = tmp_path / "only.vcf.gz"
    _write_vcf(shard_path)

    prepared = [
        module.PreparedVCFInput(
            source_path=str(shard_path),
            compressed_path=str(shard_path),
            index_path=str(shard_path) + ".tbi",
        )
    ]

    _prepare_successful_validation(monkeypatch, module, prepared)

    logged_messages: list[str] = []
    original_log_message = module.log_message

    def _capture(message, *args, **kwargs):
        logged_messages.append(message)
        return original_log_message(message, *args, **kwargs)

    monkeypatch.setattr(module, "log_message", _capture)

    valid_files, samples = module.validate_all_vcfs(
        str(tmp_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert valid_files == [str(shard_path)]
    assert samples == ["SAMPLE"]

    assert any("Discovered 1 VCF shard(s) prior to validation." in msg for msg in logged_messages)
    assert any("Validation summary: 1 of 1 shard(s) passed." in msg for msg in logged_messages)


@pytest.mark.skipif(not getattr(validation, "VCFPY_AVAILABLE", False), reason="vcfpy dependency is required for validation tests")
def test_validate_all_vcfs_logs_counts_for_nested_vcfs(tmp_path, monkeypatch):
    module = validation

    nested_dir = tmp_path / "nested"
    nested_dir.mkdir()
    shard_path = nested_dir / "shard.vcf"
    _write_vcf(shard_path)

    prepared = [
        module.PreparedVCFInput(
            source_path=str(shard_path),
            compressed_path=str(shard_path),
            index_path=str(shard_path) + ".tbi",
        )
    ]

    _prepare_successful_validation(monkeypatch, module, prepared)

    logged_messages: list[str] = []
    original_log_message = module.log_message

    def _capture(message, *args, **kwargs):
        logged_messages.append(message)
        return original_log_message(message, *args, **kwargs)

    monkeypatch.setattr(module, "log_message", _capture)

    valid_files, samples = module.validate_all_vcfs(
        str(tmp_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert valid_files == [str(shard_path)]
    assert samples == ["SAMPLE"]

    assert any("Discovered 1 VCF shard(s) prior to validation." in msg for msg in logged_messages)
    assert any("Validation summary: 1 of 1 shard(s) passed." in msg for msg in logged_messages)
