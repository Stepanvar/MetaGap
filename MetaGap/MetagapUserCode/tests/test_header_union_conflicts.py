"""Tests for header union edge cases where duplicate definitions must be rejected.

The scenarios here focus on conflicting header declarations that would make a
merged cohort unsafe. By elevating vcfpy's duplicate header warnings to errors we
verify that ``validate_vcf`` and ``validate_all_vcfs`` surface the problem via
``handle_non_critical_error`` rather than allowing the shard to pass silently.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

@pytest.fixture(autouse=True)
def _skip_without_vcfpy(merge_script_module):
    if not getattr(merge_script_module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for header conflict tests")


def _fixture_path(name: str) -> Path:
    return Path(__file__).with_name("data") / name


def _copy_fixture(tmp_path: Path, name: str) -> Path:
    target = tmp_path / name
    target.write_text(_fixture_path(name).read_text())
    return target


def test_duplicate_format_ids_fail_validation(tmp_path, monkeypatch, merge_script_module):
    module = merge_script_module

    shard = _copy_fixture(tmp_path, "duplicate_format_conflict.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    with warnings.catch_warnings():
        warnings.simplefilter(
            "error", category=module.vcfpy.header.DuplicateHeaderLineWarning
        )
        is_valid = module.validate_vcf(
            str(shard),
            ref_genome="GRCh38",
            vcf_version="VCFv4.2",
        )

    assert not is_valid
    assert any("FORMAT header" in message for message in messages)


def test_validate_all_vcfs_skips_duplicate_format_shards(
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

    _copy_fixture(tmp_path, "duplicate_format_conflict.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    with warnings.catch_warnings():
        warnings.simplefilter(
            "error", category=module.vcfpy.header.DuplicateHeaderLineWarning
        )
        valid_files, samples = module.validate_all_vcfs(
            str(tmp_path),
            ref_genome="GRCh38",
            vcf_version="VCFv4.2",
        )

    assert valid_files == [str(good_shard)]
    assert samples == ["SAMPLE"]
    assert any("FORMAT header" in message for message in messages)


def test_inconsistent_contig_definitions_are_reported(
    tmp_path, monkeypatch, merge_script_module
):
    module = merge_script_module

    shard = _copy_fixture(tmp_path, "inconsistent_contigs.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    with warnings.catch_warnings():
        warnings.simplefilter(
            "error", category=module.vcfpy.header.DuplicateHeaderLineWarning
        )
        is_valid = module.validate_vcf(
            str(shard),
            ref_genome="GRCh38",
            vcf_version="VCFv4.2",
        )

    assert not is_valid
    assert any("contig header" in message for message in messages)
