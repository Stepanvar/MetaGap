"""Regression tests for informative union_headers conflict errors."""

from __future__ import annotations

from pathlib import Path

import pytest


@pytest.fixture(autouse=True)
def _skip_without_vcfpy(merge_script_module):
    if not getattr(merge_script_module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for header union conflict tests")


def _write_header_only_vcf(path: Path, extra_header_lines: list[str]) -> None:
    header_lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        *extra_header_lines,
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
    ]
    path.write_text("\n".join(header_lines) + "\n")


@pytest.mark.parametrize("merge_script_module", ["cli"], indirect=True)
@pytest.mark.parametrize(
    "first_info_line, second_info_line, field_name, info_id",
    [
        (
            '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">',
            '##INFO=<ID=DP,Number=.,Type=Integer,Description="Depth">',
            "Number",
            "DP",
        ),
        (
            '##INFO=<ID=AD,Number=2,Type=Integer,Description="Allelic depths">',
            '##INFO=<ID=AD,Number=2,Type=Float,Description="Allelic depths">',
            "Type",
            "AD",
        ),
    ],
)
def test_union_headers_conflicting_info_fields(
    tmp_path: Path,
    merge_script_module,
    first_info_line: str,
    second_info_line: str,
    field_name: str,
    info_id: str,
) -> None:
    module = merge_script_module

    vcf_one = tmp_path / "one.vcf"
    vcf_two = tmp_path / "two.vcf"

    _write_header_only_vcf(vcf_one, [first_info_line])
    _write_header_only_vcf(vcf_two, [second_info_line])

    with pytest.raises(module.MergeConflictError) as excinfo:
        module.union_headers([str(vcf_one), str(vcf_two)])

    message = str(excinfo.value)
    assert f"Field '{field_name}'" in message
    assert f"INFO '{info_id}'" in message


@pytest.mark.parametrize("merge_script_module", ["cli"], indirect=True)
def test_union_headers_conflicting_filter_descriptions(tmp_path: Path, merge_script_module) -> None:
    module = merge_script_module

    vcf_one = tmp_path / "one.vcf"
    vcf_two = tmp_path / "two.vcf"

    _write_header_only_vcf(
        vcf_one,
        ['##FILTER=<ID=LowQ,Description="Low quality">'],
    )
    _write_header_only_vcf(
        vcf_two,
        ['##FILTER=<ID=LowQ,Description="Stricter low quality">'],
    )

    with pytest.raises(module.MergeConflictError) as excinfo:
        module.union_headers([str(vcf_one), str(vcf_two)])

    message = str(excinfo.value)
    assert "FILTER" in message
    assert "Description" in message
    assert "LowQ" in message
