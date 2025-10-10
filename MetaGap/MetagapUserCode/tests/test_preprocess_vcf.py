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


def test_preprocess_vcf_normalizes_spaces_and_preserves_original(
    tmp_path: pathlib.Path, vcf_header: str, cli_module
) -> None:
    vcf_path = tmp_path / "spaces.vcf"
    original_content = (
        vcf_header
        + "#CHROM POS  ID   REF ALT  QUAL FILTER  INFO\n"
        + "chr1  1 .  A T  . PASS .\n"
        + "chr2   2 rs1  G   C .  q10 NS=1\n"
    )
    vcf_path.write_text(original_content, encoding="utf-8")

    preprocess_vcf = _load_preprocess_vcf()

    result = preprocess_vcf(str(vcf_path), chunk_size=1)

    try:
        assert result.endswith(".tmp")
        normalized_path = pathlib.Path(result)
        normalized_text = normalized_path.read_text(encoding="utf-8")
        assert normalized_text == (
            vcf_header
            + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            + "chr1\t1\t.\tA\tT\t.\tPASS\t.\n"
            + "chr2\t2\trs1\tG\tC\t.\tq10\tNS=1\n"
        )

        # Original input remains untouched and retains spaces.
        assert vcf_path.read_text(encoding="utf-8") == original_content
    finally:
        if result != str(vcf_path):
            normalized = pathlib.Path(result)
            if normalized.exists():
                normalized.unlink()
