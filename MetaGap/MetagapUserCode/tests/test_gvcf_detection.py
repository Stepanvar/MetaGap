from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def _skip_without_vcfpy(merge_script_module):
    """Skip these tests if the target module reports vcfpy is unavailable."""
    if not getattr(merge_script_module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for gVCF detection tests")


@pytest.fixture
def gvcf_path(tmp_path):
    """Copy the sample gVCF into a temp dir so it can be safely modified per test."""
    fixture_path = Path(__file__).with_name("data") / "sample.g.vcf"
    destination = tmp_path / "sample.g.vcf"
    destination.write_text(fixture_path.read_text())
    return destination


def test_validate_vcf_rejects_gvcf_header(monkeypatch, gvcf_path, merge_script_module):
    messages = []

    def capture_warning(message):
        messages.append(message)

    monkeypatch.setattr(merge_script_module, "handle_non_critical_error", capture_warning)

    is_valid = merge_script_module.validate_vcf(
        str(gvcf_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert not is_valid
    assert any("gVCF" in message for message in messages)


def test_validate_vcf_rejects_non_ref_alt(monkeypatch, tmp_path, merge_script_module):
    fixture_path = Path(__file__).with_name("data") / "sample.g.vcf"
    alt_only_path = tmp_path / "sample_without_gvcf_header.vcf"
    contents = fixture_path.read_text().replace("##GVCFBlock=1-100\n", "")
    alt_only_path.write_text(contents)

    messages = []

    def capture_warning(message):
        messages.append(message)

    monkeypatch.setattr(merge_script_module, "handle_non_critical_error", capture_warning)

    is_valid = merge_script_module.validate_vcf(
        str(alt_only_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert not is_valid
    assert any("<NON_REF>" in message for message in messages)
