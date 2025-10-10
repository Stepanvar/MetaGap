import importlib
import pathlib

import pytest


@pytest.fixture
def vcf_header() -> str:
    return "##fileformat=VCFv4.2\n"


def _load_preprocess_vcf():
    module = importlib.import_module("MetagapUserCode.merge_vcf.merging")
    return getattr(module, "preprocess_vcf")


def test_preprocess_vcf_no_changes(
    tmp_path: pathlib.Path, vcf_header: str, cli_module
) -> None:
    vcf_path = tmp_path / "sample.vcf"
    vcf_path.write_text(
        vcf_header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "chr1\t1\t.\tA\tT\t.\tPASS\t.\n",
        encoding="utf-8",
    )

    preprocess_vcf = _load_preprocess_vcf()

    result = preprocess_vcf(str(vcf_path))

    assert result == str(vcf_path)
    assert vcf_path.read_text(encoding="utf-8") == (
        vcf_header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "chr1\t1\t.\tA\tT\t.\tPASS\t.\n"
    )


def test_preprocess_vcf_requires_normalization(
    tmp_path: pathlib.Path, vcf_header: str, cli_module
) -> None:
    vcf_path = tmp_path / "normalize_me.vcf"
    vcf_path.write_text(
        vcf_header
        + "#CHROM POS ID REF ALT QUAL FILTER INFO\n"
        + "chr1 1 . A T . PASS .\n"
        + "chr2 2 rs1 G C . q10 NS=1\n",
        encoding="utf-8",
    )

    preprocess_vcf = _load_preprocess_vcf()

    result = preprocess_vcf(str(vcf_path), chunk_size=1)

    assert result != str(vcf_path)
    normalized_text = pathlib.Path(result).read_text(encoding="utf-8")
    assert normalized_text == (
        vcf_header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "chr1\t1\t.\tA\tT\t.\tPASS\t.\n"
        + "chr2\t2\trs1\tG\tC\t.\tq10\tNS=1\n"
    )
    # Original input remains untouched.
    assert vcf_path.read_text(encoding="utf-8") == (
        vcf_header
        + "#CHROM POS ID REF ALT QUAL FILTER INFO\n"
        + "chr1 1 . A T . PASS .\n"
        + "chr2 2 rs1 G C . q10 NS=1\n"
    )


def test_preprocess_vcf_converts_spaces_to_tabs(
    tmp_path: pathlib.Path, vcf_header: str, cli_module
) -> None:
    vcf_path = tmp_path / "mixed_spaces.vcf"
    original_text = (
        vcf_header
        + "#CHROM POS  ID   REF ALT  QUAL FILTER INFO\n"
        + "chr1  1 .  A T  . PASS .\n"
        + "chr2   2 rs1 G  C .  q10 NS=1\n"
    )
    vcf_path.write_text(original_text, encoding="utf-8")

    preprocess_vcf = _load_preprocess_vcf()

    result = preprocess_vcf(str(vcf_path), chunk_size=1)
    result_path = pathlib.Path(result)

    try:
        assert result_path.name.endswith(".tmp")

        normalized_text = result_path.read_text(encoding="utf-8")

        lines = normalized_text.strip().splitlines()
        # Skip comment metadata lines; ensure tab separation for header and records.
        content_lines = [line for line in lines if not line.startswith("##")]
        assert content_lines, "Expected header and records in normalized output"
        for line in content_lines:
            # Each VCF line should have eight columns separated by tabs, thus seven tab characters.
            assert line.count("\t") == 7

        # Original input should remain unchanged and retain spaces.
        assert vcf_path.read_text(encoding="utf-8") == original_text
    finally:
        if result_path.exists() and result_path != vcf_path:
            result_path.unlink()
