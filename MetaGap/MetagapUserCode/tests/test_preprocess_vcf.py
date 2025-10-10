import gzip
import importlib
import pathlib

import pytest


@pytest.fixture
def vcf_header() -> str:
    return "##fileformat=VCFv4.2\n"


@pytest.fixture(params=(".vcf", ".vcf.gz"))
def pre_tabbed_vcf(tmp_path: pathlib.Path, vcf_header: str, request) -> tuple[pathlib.Path, str]:
    extension = request.param
    vcf_path = tmp_path / f"pre_tabbed{extension}"
    text = (
        vcf_header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "chr1\t1\t.\tA\tT\t.\tPASS\t.\n"
    )

    if extension.endswith(".gz"):
        with gzip.open(vcf_path, "wt", encoding="utf-8") as handle:
            handle.write(text)
    else:
        vcf_path.write_text(text, encoding="utf-8")

    return vcf_path, text


def _load_preprocess_vcf():
    module = importlib.import_module("MetagapUserCode.merge_vcf.merging")
    return getattr(module, "preprocess_vcf")


@pytest.fixture(params=[".vcf", ".vcf.gz"])
def pre_tabbed_vcf(
    tmp_path: pathlib.Path, vcf_header: str, request: pytest.FixtureRequest
) -> tuple[pathlib.Path, str]:
    content = (
        vcf_header
        + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        + "chr1\t1\t.\tA\tT\t.\tPASS\t.\n"
    )
    suffix = request.param
    vcf_path = tmp_path / f"pre_tabbed{suffix}"

    if suffix.endswith(".gz"):
        with gzip.open(vcf_path, "wt", encoding="utf-8") as fh:
            fh.write(content)
    else:
        vcf_path.write_text(content, encoding="utf-8")

    return vcf_path, content


def test_preprocess_vcf_preserves_pre_tabbed(
    pre_tabbed_vcf: tuple[pathlib.Path, str], tmp_path: pathlib.Path, cli_module
) -> None:
    preprocess_vcf = _load_preprocess_vcf()

    vcf_path, expected_content = pre_tabbed_vcf

    result = preprocess_vcf(str(vcf_path))
    assert result == str(vcf_path)

    if vcf_path.suffix == ".gz":
        with gzip.open(vcf_path, "rt", encoding="utf-8") as fh:
            assert fh.read() == expected_content
    else:
        assert vcf_path.read_text(encoding="utf-8") == expected_content

    # no temp files left behind
    assert not list(tmp_path.glob("*.tmp"))


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
