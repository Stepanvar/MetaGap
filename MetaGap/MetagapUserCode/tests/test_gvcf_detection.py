from pathlib import Path
import importlib.util
from functools import lru_cache

import pytest


_MODULE_PATH = Path(__file__).resolve().parents[1] / "test_merge_vcf.py"
_DATA_DIR = Path(__file__).with_name("data")
_BASE_GVCF_CONTENT = (_DATA_DIR / "sample.g.vcf").read_text()


@lru_cache(maxsize=1)
def _load_test_merge_vcf_module():
    spec = importlib.util.spec_from_file_location("test_merge_vcf_module", _MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def test_merge_vcf_module():
    if not _MODULE_PATH.exists():
        pytest.skip("test_merge_vcf.py is required for gVCF detection tests")

    module = _load_test_merge_vcf_module()
    if not getattr(module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for gVCF detection tests")
    return module


@pytest.mark.parametrize(
    "filename, contents, expected_keyword",
    [
        pytest.param(
            "sample.g.vcf",
            _BASE_GVCF_CONTENT,
            "gVCF",
            id="gvcf-header",
        ),
        pytest.param(
            "sample_without_gvcf_header.vcf",
            _BASE_GVCF_CONTENT.replace("##GVCFBlock=1-100\n", ""),
            "<NON_REF>",
            id="non-ref-alt",
        ),
    ],
)
def test_validate_vcf_rejects_gvcf_variants(
    filename,
    contents,
    expected_keyword,
    tmp_path,
    monkeypatch,
    test_merge_vcf_module,
):
    shard_path = tmp_path / filename
    shard_path.write_text(contents)

    messages = []
    monkeypatch.setattr(test_merge_vcf_module, "handle_non_critical_error", messages.append)

    is_valid = test_merge_vcf_module.validate_vcf(
        str(shard_path),
        ref_genome="GRCh38",
        vcf_version="VCFv4.2",
    )

    assert is_valid is False
    assert any(expected_keyword in message for message in messages)
