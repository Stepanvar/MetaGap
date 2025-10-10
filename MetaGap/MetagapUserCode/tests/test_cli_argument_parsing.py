import sys

import pytest


def test_parse_arguments_requires_input_dir(monkeypatch, cli_module):
    monkeypatch.setattr(sys, "argv", ["merge_vcf"])

    with pytest.raises(SystemExit):
        cli_module.parse_arguments()


def test_parse_arguments_rejects_unknown_metadata_tokens(monkeypatch, cli_module, tmp_path):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()

    argv = [
        "merge_vcf",
        "--input-dir",
        str(input_dir),
        "--meta",
        "--UNEXPECTED=Value",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(SystemExit):
        cli_module.parse_arguments()


def test_parse_arguments_supports_optional_flags(monkeypatch, cli_module, tmp_path):
    input_dir = tmp_path / "inputs"
    output_dir = tmp_path / "outputs"
    template_path = tmp_path / "metadata_template.txt"

    input_dir.mkdir()
    output_dir.mkdir()
    template_path.write_text("##dummy=entry\n")

    argv = [
        "merge_vcf",
        "--input-dir",
        str(input_dir),
        "--output",
        str(output_dir),
        "--ref",
        "GRCh38",
        "--vcf-version",
        "4.2",
        "--allow-gvcf",
        "--meta",
        "ID=SAMPLE1",
        "--meta",
        "Sex=Female",
        "--sample-metadata",
        "ID=SAMPLE1",
        "--sample-metadata",
        "Status=Affected",
        "--header-metadata",
        "INFO=DP=1",
        "--metadata-template",
        str(template_path),
        "--verbose",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    args = cli_module.parse_arguments()

    assert args.input_dir == str(input_dir)
    assert args.output_dir == str(output_dir)
    assert args.ref == "GRCh38"
    assert args.vcf_version == "4.2"
    assert args.allow_gvcf is True
    assert args.meta == ["ID=SAMPLE1", "Sex=Female"]
    assert args.sample_metadata_entries == ["ID=SAMPLE1", "Status=Affected"]
    assert args.header_metadata_lines == ["INFO=DP=1"]
    assert args.metadata_template_path == str(template_path)
    assert args.verbose is True
