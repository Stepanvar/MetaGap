from pathlib import Path
import importlib.util

import pytest


_MODULE_PATH = Path(__file__).resolve().parents[1] / "test_merge_vcf.py"
_SPEC = importlib.util.spec_from_file_location("test_merge_vcf_module", _MODULE_PATH)
test_merge_vcf = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None  # for mypy/type checkers
_SPEC.loader.exec_module(test_merge_vcf)

pytestmark = pytest.mark.skipif(
    not getattr(test_merge_vcf, "VCFPY_AVAILABLE", True),
    reason="vcfpy dependency is required for gVCF detection tests",
)


@pytest.fixture
def gvcf_path(tmp_path):
    fixture_path = Path(__file__).with_name("data") / "sample.g.vcf"
    destination = tmp_path / "sample.g.vcf"
    destination.write_text(fixture_path.read_text())
    return destination


def test_validate_vcf_rejects_gvcf_header(monkeypatch, gvcf_path):
    messages = []

    def capture_warning(message):
        messages.append(message)

    monkeypatch.setattr(test_merge_vcf, "handle_non_critical_error", capture_warning)

    is_valid = test_merge_vcf.validate_vcf(
        str(gvcf_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert not is_valid
    assert any("gVCF" in message for message in messages)


def test_validate_vcf_rejects_non_ref_alt(monkeypatch, tmp_path):
    fixture_path = Path(__file__).with_name("data") / "sample.g.vcf"
    alt_only_path = tmp_path / "sample_without_gvcf_header.vcf"
    contents = fixture_path.read_text().replace("##GVCFBlock=1-100\n", "")
    alt_only_path.write_text(contents)

    messages = []

    def capture_warning(message):
        messages.append(message)

    monkeypatch.setattr(test_merge_vcf, "handle_non_critical_error", capture_warning)

    is_valid = test_merge_vcf.validate_vcf(
        str(alt_only_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert not is_valid
    assert any("<NON_REF>" in message for message in messages)
