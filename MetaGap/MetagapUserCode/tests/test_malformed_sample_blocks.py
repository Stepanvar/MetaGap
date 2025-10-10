"""Coverage for malformed ``##SAMPLE`` metadata blocks encountered during validation.

These tests document that missing sample identifiers trigger a non-critical
warning and prevent the offending shard from participating in downstream merge
workflows.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[2] / "MetagapUserCode" / "test_merge_vcf.py"


def load_user_module():
    spec = importlib.util.spec_from_file_location("user_test_merge_vcf", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)  # type: ignore[call-arg]
    return module


_module = load_user_module()
pytestmark = pytest.mark.skipif(
    not getattr(_module, "VCFPY_AVAILABLE", True),
    reason="vcfpy dependency is required for sample metadata tests",
)


def _copy_fixture(tmp_path: Path, name: str) -> Path:
    source = Path(__file__).with_name("data") / name
    destination = tmp_path / name
    destination.write_text(source.read_text())
    return destination


def test_validate_vcf_flags_missing_sample_id(tmp_path, monkeypatch):
    module = load_user_module()

    shard = _copy_fixture(tmp_path, "malformed_sample_metadata.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    is_valid = module.validate_vcf(
        str(shard),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert not is_valid
    assert any("Missing key \"ID\"" in message for message in messages)


def test_validate_all_vcfs_skips_malformed_sample_metadata(tmp_path, monkeypatch):
    module = load_user_module()

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

    _copy_fixture(tmp_path, "malformed_sample_metadata.vcf")

    messages: list[str] = []
    monkeypatch.setattr(module, "handle_non_critical_error", messages.append)

    valid_files, samples = module.validate_all_vcfs(
        str(tmp_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert valid_files == [str(good_shard)]
    assert samples == ["SAMPLE"]
    assert any("Missing key \"ID\"" in message for message in messages)
