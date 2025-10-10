import importlib.util
import os
import sys
from pathlib import Path
from types import ModuleType, SimpleNamespace

import pytest


@pytest.fixture
def cli_module(monkeypatch):
    base_dir = Path(__file__).resolve().parents[1]

    package_module = ModuleType("MetagapUserCode")
    package_module.__path__ = [str(base_dir)]
    monkeypatch.setitem(sys.modules, "MetagapUserCode", package_module)

    merge_pkg_name = "MetagapUserCode.merge_vcf"
    merge_pkg_path = base_dir / "merge_vcf" / "__init__.py"
    merge_spec = importlib.util.spec_from_file_location(
        merge_pkg_name,
        merge_pkg_path,
        submodule_search_locations=[str(base_dir / "merge_vcf")],
    )
    merge_pkg = importlib.util.module_from_spec(merge_spec)
    assert merge_spec.loader is not None
    monkeypatch.setitem(sys.modules, merge_pkg_name, merge_pkg)
    merge_spec.loader.exec_module(merge_pkg)

    module_name = "MetagapUserCode.merge_vcf.cli"
    module_path = base_dir / "merge_vcf" / "cli.py"
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    monkeypatch.setitem(sys.modules, module_name, module)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def gvcf_inputs(tmp_path):
    input_dir = tmp_path / "input"
    output_dir = tmp_path / "output"
    input_dir.mkdir()
    output_dir.mkdir()
    fixture_path = Path(__file__).with_name("data") / "sample.g.vcf"
    gvcf_path = input_dir / "sample.g.vcf"
    gvcf_path.write_text(fixture_path.read_text())
    return gvcf_path, input_dir, output_dir


def _build_args(input_dir, output_dir, allow_gvcf):
    return SimpleNamespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        ref="GRCh38",
        vcf_version="4.2",
        meta=[],
        sample_metadata_entries=[],
        header_metadata_lines=[],
        verbose=False,
        allow_gvcf=allow_gvcf,
    )


def test_cli_accepts_gvcf_when_flag_enabled(monkeypatch, cli_module, gvcf_inputs):
    gvcf_path, input_dir, output_dir = gvcf_inputs
    recorded = {}

    def fake_validate_all_vcfs(
        input_dir,
        ref_genome,
        vcf_version,
        verbose=False,
        allow_gvcf=False,
    ):
        assert allow_gvcf is True
        discovered = sorted(Path(input_dir).glob("*.vcf"))
        recorded["validated_files"] = [str(path) for path in discovered]
        return recorded["validated_files"], ["SAMPLE"]

    def fake_merge_vcfs(valid_files, output_dir, verbose, **kwargs):
        recorded["valid_files"] = list(valid_files)
        merged_path = os.path.join(output_dir, "merged.vcf")
        Path(merged_path).write_text(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "1\t100\t.\tA\tT\t.\tPASS\t.\n"
        )
        return merged_path

    monkeypatch.setattr(cli_module, "validate_all_vcfs", fake_validate_all_vcfs)
    monkeypatch.setattr(cli_module, "merge_vcfs", fake_merge_vcfs)
    monkeypatch.setattr(
        cli_module,
        "append_metadata_to_merged_vcf",
        lambda merged_vcf, **kwargs: merged_vcf,
    )
    monkeypatch.setattr(
        cli_module,
        "validate_merged_vcf",
        lambda *args, **kwargs: None,
    )
    monkeypatch.setattr(
        cli_module,
        "summarize_produced_vcfs",
        lambda output_dir, final_vcf: (output_dir, os.path.basename(final_vcf), 1),
    )

    args = _build_args(input_dir, output_dir, allow_gvcf=True)
    monkeypatch.setattr(cli_module, "parse_arguments", lambda: args)

    cli_module.main()

    assert recorded["validated_files"] == [str(gvcf_path)]
    assert recorded["valid_files"] == recorded["validated_files"]


def test_cli_rejects_gvcf_without_flag(monkeypatch, cli_module, gvcf_inputs, capsys):
    _, input_dir, output_dir = gvcf_inputs

    monkeypatch.setattr(
        cli_module,
        "merge_vcfs",
        lambda *args, **kwargs: pytest.fail("merge_vcfs should not run when gVCFs are rejected"),
    )
    monkeypatch.setattr(
        cli_module,
        "append_metadata_to_merged_vcf",
        lambda *args, **kwargs: pytest.fail(
            "append_metadata_to_merged_vcf should not run when gVCFs are rejected"
        ),
    )

    args = _build_args(input_dir, output_dir, allow_gvcf=False)
    monkeypatch.setattr(cli_module, "parse_arguments", lambda: args)

    with pytest.raises(SystemExit):
        cli_module.main()

    captured = capsys.readouterr()
    assert "Use --allow-gvcf" in captured.out
