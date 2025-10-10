"""Regression tests for sample blocks that imply mixed ploidy depth encodings.

These cases rely on vcfpy's conversion warnings to model inputs where genotype
fields carry more alleles than declared. The tests ensure ``validate_vcf`` and
``validate_all_vcfs`` treat the shards as invalid by routing the warning through
``handle_non_critical_error`` rather than allowing silent coercion.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

@pytest.fixture(autouse=True)
def _skip_without_vcfpy(merge_script_module):
    if not getattr(merge_script_module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for mixed ploidy tests")


def _copy_fixture(tmp_path: Path, name: str) -> Path:
    source = Path(__file__).with_name("data") / name
    destination = tmp_path / name
    destination.write_text(source.read_text())
    return destination


def _filter_integer_conversion():
    warnings.filterwarnings("error", r".*cannot be converted to Integer.*")


def test_validate_vcf_rejects_mixed_depth_encoding(
    tmp_path, monkeypatch, merge_script_module
):
    module = merge_script_module

    shard = _copy_fixture(tmp_path, "mixed_ploidy_depth.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    with warnings.catch_warnings():
        _filter_integer_conversion()
        is_valid = module.validate_vcf(
            str(shard),
            ref_genome="GRCh38",
            vcf_version="VCFv4.2",
        )

    assert not is_valid
    assert any("cannot be converted to Integer" in message for message in messages)


def test_validate_all_vcfs_skips_mixed_ploidy_shard(
    tmp_path, monkeypatch, merge_script_module
):
    module = merge_script_module

    good_text = "\n".join(
        [
            "##fileformat=VCFv4.2",
            "##reference=GRCh38",
            "##contig=<ID=1,length=1000>",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
            "1\t100\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
        ]
    )
    good_shard = tmp_path / "good.vcf"
    good_shard.write_text(good_text + "\n")

    _copy_fixture(tmp_path, "mixed_ploidy_depth.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    with warnings.catch_warnings():
        _filter_integer_conversion()
        valid_files, samples = module.validate_all_vcfs(
            str(tmp_path),
            ref_genome="GRCh38",
            vcf_version="VCFv4.2",
        )

    assert valid_files == [str(good_shard)]
    assert samples == ["SAMPLE"]
    assert any("cannot be converted to Integer" in message for message in messages)
