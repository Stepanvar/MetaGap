"""Regression test ensuring :func:`validate_all_vcfs` imports ``glob``."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path
from unittest.mock import patch

import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))

validation = importlib.import_module("merge_vcf.validation")


def test_validate_all_vcfs_handles_glob_import(tmp_path: Path):
    """Ensure ``validate_all_vcfs`` runs without raising ``NameError``."""

    module = validation

    if not getattr(module, "VCFPY_AVAILABLE", False):
        pytest.skip("vcfpy dependency is required for validation tests")

    shard_path = tmp_path / "shard.vcf"
    shard_path.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##reference=GRCh38",
                "##contig=<ID=1,length=1000>",
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
                "1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t0/1",
            ]
        )
        + "\n"
    )

    prepared = module.PreparedVCFInput(
        source_path=str(shard_path),
        compressed_path=str(shard_path),
        index_path=str(shard_path) + ".tbi",
    )

    def fake_preprocess(file_path: str, *, chunk_size: int = 1024) -> str:
        return str(file_path)

    valid_files = samples = None

    with (
        patch.object(
            module,
            "discover_and_prepare_inputs",
            return_value=[prepared],
        ),
        patch.object(module, "preprocess_vcf", side_effect=fake_preprocess),
    ):
        try:
            valid_files, samples = module.validate_all_vcfs(
                str(tmp_path),
                ref_genome="GRCh38",
                vcf_version="VCFv4.2",
            )
        except NameError as exc:  # pragma: no cover - explicit regression guard
            pytest.fail(f"validate_all_vcfs raised unexpected NameError: {exc}")

    assert valid_files == [str(shard_path)]
    assert samples == ["SAMPLE"]
