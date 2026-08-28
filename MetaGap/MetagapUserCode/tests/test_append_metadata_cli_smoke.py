import gzip
import io
import sys
from pathlib import Path

import pytest

SAMPLE_HEADER = """##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##SAMPLE=<ID=S1,Description="Primary sample">
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""
SAMPLE_BODY_WITH_FORMAT = "1\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=2;AF=0.5\tGT\t0/1\n"
SAMPLE_BODY_TRIMMED = "1\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=2;AF=0.5\n"
SAMPLE_VCF = SAMPLE_HEADER + SAMPLE_BODY_WITH_FORMAT


def _write_merged_vcf(path: Path) -> None:
    """Write a realistic merged VCF for the native (vcfpy/pysam) finalize path.

    The merge pipeline no longer shells out to bcftools/bgzip/tabix -- and
    ``merge_vcf/__init__.py`` statically forbids it -- so these tests drive the
    real implementation instead of faking a subprocess that is never spawned.
    """

    if str(path).endswith(".gz"):
        with gzip.open(path, "wt", encoding="utf-8") as handle:
            handle.write(SAMPLE_VCF)
    else:
        path.write_text(SAMPLE_VCF)


def _metadata_module(module):
    return sys.modules[module.append_metadata_to_merged_vcf.__module__]


def _spy_header_validation(monkeypatch, module, header_checks):
    """Record which artifact the anonymized-header validation runs against."""

    metadata_module = _metadata_module(module)
    original = metadata_module._validate_anonymized_vcf_header

    def spy(final_vcf_path, *args, **kwargs):
        header_checks.append(Path(final_vcf_path))
        return original(final_vcf_path, *args, **kwargs)

    monkeypatch.setattr(metadata_module, "_validate_anonymized_vcf_header", spy)


def test_append_metadata_runs_header_validation_for_gzip(
    monkeypatch, tmp_path, merge_script_module
):
    module = merge_script_module
    header_checks = []
    _spy_header_validation(monkeypatch, module, header_checks)

    merged_vcf = tmp_path / "merged.vcf.gz"
    _write_merged_vcf(merged_vcf)

    final_vcf = module.append_metadata_to_merged_vcf(str(merged_vcf))

    assert final_vcf.endswith(".vcf.gz")
    assert Path(final_vcf) == header_checks[0]

    with gzip.open(final_vcf, "rt", encoding="utf-8") as handle:
        contents = handle.read()
    assert "##fileformat=VCFv4.2" in contents
    assert "##FORMAT=" not in contents
    assert "##SAMPLE=" not in contents
    assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in contents


def test_append_metadata_runs_header_validation_for_plain_vcf(
    monkeypatch, tmp_path, merge_script_module
):
    module = merge_script_module
    header_checks = []
    _spy_header_validation(monkeypatch, module, header_checks)

    merged_vcf = tmp_path / "merged.vcf"
    _write_merged_vcf(merged_vcf)

    final_vcf = module.append_metadata_to_merged_vcf(str(merged_vcf))

    assert final_vcf.endswith(".vcf")
    assert Path(final_vcf) == header_checks[0]
    contents = Path(final_vcf).read_text()
    assert "##fileformat=VCFv4.2" in contents
    assert "##FORMAT=" not in contents
    assert "##SAMPLE=" not in contents


def test_append_metadata_compresses_and_indexes(tmp_path, merge_script_module):
    """Was ``test_append_metadata_triggers_bgzip_and_tabix``.

    Same requirement -- the gzip artifact is BGZF-compressed, tabix-indexed, and
    the plain intermediate is removed -- asserted against the pysam calls that
    replaced the bgzip/tabix subprocesses.
    """

    module = merge_script_module

    merged_vcf = tmp_path / "merged.vcf.gz"
    _write_merged_vcf(merged_vcf)

    final_vcf = module.append_metadata_to_merged_vcf(str(merged_vcf))

    final_path = Path(final_vcf)
    assert final_path.suffix == ".gz"
    assert final_path.exists()
    assert final_path.with_suffix(final_path.suffix + ".tbi").exists()

    plain_path = Path(str(final_vcf)[:-3])
    assert not plain_path.exists()


def test_validate_header_failure_raises(monkeypatch, tmp_path, merge_script_module):
    module = merge_script_module
    metadata_module = _metadata_module(module)

    captured = {}

    def fake_handle_error(message, **_kwargs):
        captured["message"] = message
        raise RuntimeError("critical failure")

    monkeypatch.setattr(metadata_module, "handle_critical_error", fake_handle_error)

    unreadable = tmp_path / "broken.vcf.gz"
    unreadable.write_bytes(b"this is not a VCF")

    with pytest.raises(RuntimeError):
        metadata_module._validate_anonymized_vcf_header(str(unreadable))

    assert "vcfpy" in captured["message"]
