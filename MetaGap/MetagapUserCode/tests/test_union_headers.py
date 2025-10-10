"""Targeted regression tests for ``union_headers`` conflict handling."""

from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def _skip_without_vcfpy(merge_script_module):
    if not getattr(merge_script_module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for union header tests")


def _merge_conflict_error_cls(module):
    error_cls = getattr(module, "MergeConflictError", None)
    if error_cls is not None:
        return error_cls

    logging_utils = getattr(module, "logging_utils", None)
    if logging_utils is not None:
        error_cls = getattr(logging_utils, "MergeConflictError", None)
        if error_cls is not None:
            return error_cls

    merge_pkg = getattr(module, "merge_vcf", None)
    if merge_pkg is not None:
        error_cls = getattr(getattr(merge_pkg, "logging_utils", None), "MergeConflictError", None)
        if error_cls is not None:
            return error_cls

    raise AssertionError("MergeConflictError could not be located on merge module")


def _write_header_only_vcf(path: Path, *extra_lines: str) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        *extra_lines,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    path.write_text("\n".join(lines) + "\n")


@pytest.mark.parametrize(
    "first_info,second_info,expected",
    [
        (
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##INFO=<ID=DP,Number=.,Type=Integer,Description="Depth">',
            "Field 'Number' for INFO 'DP'",
        ),
        (
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##INFO=<ID=DP,Number=1,Type=Float,Description="Depth">',
            "Field 'Type' for INFO 'DP'",
        ),
    ],
)
def test_union_headers_conflicting_info_definitions(
    tmp_path: Path,
    merge_script_module,
    first_info: str,
    second_info: str,
    expected: str,
):
    module = merge_script_module

    vcf_one = tmp_path / "info_one.vcf"
    vcf_two = tmp_path / "info_two.vcf"

    _write_header_only_vcf(vcf_one, first_info)
    _write_header_only_vcf(vcf_two, second_info)

    with pytest.raises(_merge_conflict_error_cls(module)) as excinfo:
        module.union_headers([str(vcf_one), str(vcf_two)])

    assert expected in str(excinfo.value)


def test_union_headers_conflicting_filter_definitions(tmp_path: Path, merge_script_module):
    module = merge_script_module

    vcf_one = tmp_path / "filter_one.vcf"
    vcf_two = tmp_path / "filter_two.vcf"

    _write_header_only_vcf(
        vcf_one,
        '##FILTER=<ID=lowq,Description="Fails quality">',
    )
    _write_header_only_vcf(
        vcf_two,
        '##FILTER=<ID=lowq,Description="Low quality threshold">',
    )

    with pytest.raises(_merge_conflict_error_cls(module)) as excinfo:
        module.union_headers([str(vcf_one), str(vcf_two)])

    assert "Description for FILTER 'lowq'" in str(excinfo.value)
