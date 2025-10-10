import os
from pathlib import Path
from types import SimpleNamespace

import pytest


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


@pytest.mark.parametrize("merge_script_module", ["cli"], indirect=True)
def test_cli_accepts_gvcf_when_flag_enabled(
    monkeypatch, merge_script_module, gvcf_inputs
):
    cli_module = merge_script_module
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


@pytest.mark.parametrize("merge_script_module", ["cli"], indirect=True)
def test_cli_rejects_gvcf_without_flag(
    monkeypatch, merge_script_module, gvcf_inputs, capsys
):
    cli_module = merge_script_module
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
