"""Tests for detailed INFO/FORMAT header conflict reporting."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

validation = importlib.import_module("merge_vcf.validation")


@pytest.mark.skipif(
    not getattr(validation, "VCFPY_AVAILABLE", False),
    reason="vcfpy dependency is required for header conflict tests",
)
def test_validate_all_vcfs_reports_number_type_conflicts(tmp_path):
    vcf_one = tmp_path / "one.vcf"
    vcf_two = tmp_path / "two.vcf"

    common_lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##contig=<ID=1,length=1000>",
    ]

    vcf_one.write_text(
        "\n".join(
            common_lines
            + [
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1",
                "1\t10\t.\tA\tC\t.\tPASS\t.\tAD\t10,5",
            ]
        )
        + "\n"
    )

    vcf_two.write_text(
        "\n".join(
            common_lines
            + [
                '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Allelic depths">',
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS2",
                "1\t20\t.\tG\tT\t.\tPASS\t.\tAD\t15",
            ]
        )
        + "\n"
    )

    messages: list[str] = []

    original_handle_critical_error = validation.handle_critical_error
    try:
        def fake_handle_critical_error(message: str) -> None:
            messages.append(message)
            raise SystemExit(message)

        validation.handle_critical_error = fake_handle_critical_error

        with pytest.raises(SystemExit) as excinfo:
            validation.validate_all_vcfs(
                str(tmp_path),
                ref_genome="GRCh38",
                vcf_version="VCFv4.2",
            )
    finally:
        validation.handle_critical_error = original_handle_critical_error

    assert messages, "Expected the validation to report a critical error"
    message = messages[0]
    assert excinfo.value.code == message
    assert "FORMAT header definitions conflict across shards" in message
    assert "AD" in message
    assert str(vcf_one) in message
    assert str(vcf_two) in message
    assert "Number='R', Type='Integer'" in message
    assert "Number=1" in message and "Type='Integer'" in message
