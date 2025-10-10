import importlib
import importlib.util
import os
import sys
from pathlib import Path
from types import SimpleNamespace

_fixture_spec = importlib.util.spec_from_file_location(
    "cli_fixture_module", Path(__file__).with_name("test_cli_gvcf_flag.py")
)
_fixture_module = importlib.util.module_from_spec(_fixture_spec)
assert _fixture_spec.loader is not None
_fixture_spec.loader.exec_module(_fixture_module)
cli_module = _fixture_module.cli_module  # noqa: F401  (re-exported fixture)


def test_cli_metadata_integration(tmp_path, monkeypatch, cli_module):
    metadata_module = importlib.import_module("MetagapUserCode.merge_vcf.metadata")

    if not hasattr(metadata_module, "template_sample_mapping"):
        monkeypatch.setattr(metadata_module, "template_sample_mapping", None, raising=False)
    if not hasattr(metadata_module, "template_simple_lines"):
        monkeypatch.setattr(metadata_module, "template_simple_lines", [], raising=False)
    if not hasattr(metadata_module, "template_header_lines"):
        monkeypatch.setattr(metadata_module, "template_header_lines", [], raising=False)
    if not hasattr(metadata_module, "template_serialized_sample"):
        monkeypatch.setattr(metadata_module, "template_serialized_sample", None, raising=False)

    class _SampleHeaderLine:
        def __init__(self, mapping):
            self.mapping = dict(mapping)
            self.key = "SAMPLE"

        @classmethod
        def from_mapping(cls, mapping):
            return cls(mapping)

    class _SimpleHeaderLine:
        def __init__(self, key, value):
            self.key = key
            self.value = value

    stub = SimpleNamespace(
        SampleHeaderLine=_SampleHeaderLine,
        SimpleHeaderLine=_SimpleHeaderLine,
    )

    monkeypatch.setattr(metadata_module, "vcfpy", stub, raising=False)
    monkeypatch.setattr(metadata_module, "VCFPY_AVAILABLE", True, raising=False)

    input_dir = tmp_path / "inputs"
    output_dir = tmp_path / "outputs"
    input_dir.mkdir()
    output_dir.mkdir()

    stub_vcf = input_dir / "sample.vcf"
    stub_vcf.write_text(
        "\n".join(
            [
                "##fileformat=VCFv4.2",
                "##reference=GRCh38",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
                "1\t100\t.\tA\tT\t.\tPASS\t.",
            ]
        )
        + "\n"
    )

    captured = {}

    def fake_validate_all_vcfs(input_dir, ref_genome, vcf_version, verbose=False, allow_gvcf=False):
        captured["validated"] = {
            "input_dir": input_dir,
            "ref_genome": ref_genome,
            "vcf_version": vcf_version,
            "allow_gvcf": allow_gvcf,
        }
        return [str(stub_vcf)], ["SAMPLE1"]

    def fake_merge_vcfs(valid_files, output_dir, verbose, **kwargs):
        captured["merge"] = {
            "valid_files": list(valid_files),
            "output_dir": output_dir,
            "kwargs": kwargs,
        }
        merged_path = os.path.join(output_dir, "merged_output.vcf")
        Path(merged_path).write_text("##fileformat=VCFv4.2\n")
        return merged_path

    def fake_append_metadata(
        merged_vcf,
        sample_metadata_entries=None,
        header_metadata_lines=None,
        serialized_sample_line=None,
        verbose=False,
    ):
        header_metadata_lines = list(header_metadata_lines or [])
        final_path = os.path.join(output_dir, "final_output.vcf")
        header_lines = ["##fileformat=VCFv4.2"] + header_metadata_lines
        if serialized_sample_line:
            header_lines.append(serialized_sample_line)
        Path(final_path).write_text("\n".join(header_lines) + "\n")

        captured["append"] = {
            "merged_vcf": merged_vcf,
            "sample_metadata_entries": list(sample_metadata_entries or []),
            "header_metadata_lines": header_metadata_lines,
            "serialized_sample_line": serialized_sample_line,
        }
        captured["final_path"] = final_path
        captured["header_text"] = Path(final_path).read_text()
        return final_path

    def fake_validate_merged_vcf(final_vcf, verbose=False):
        captured["validated_final"] = final_vcf

    monkeypatch.setattr(cli_module, "validate_all_vcfs", fake_validate_all_vcfs)
    monkeypatch.setattr(cli_module, "merge_vcfs", fake_merge_vcfs)
    monkeypatch.setattr(cli_module, "append_metadata_to_merged_vcf", fake_append_metadata)
    monkeypatch.setattr(cli_module, "validate_merged_vcf", fake_validate_merged_vcf)
    monkeypatch.setattr(
        cli_module,
        "summarize_produced_vcfs",
        lambda output_dir, final_vcf: (output_dir, os.path.basename(final_vcf), 1),
    )

    metadata_args = [
        "--meta",
        "ID=SAMPLE001",
        "--meta",
        "Description=Study Cohort",
        "--sample-metadata",
        "Site=MGH",
        "--sample-metadata",
        "Investigator=Doe,John",
        "--header-metadata",
        "analysis=pilot",
    ]

    argv = [
        "merge-vcf",
        "--input-dir",
        str(input_dir),
        "--output-dir",
        str(output_dir),
        "--ref",
        "GRCh38",
        "--vcf-version",
        "4.2",
        *metadata_args,
    ]

    monkeypatch.setattr(sys, "argv", argv)

    cli_module.main()

    header_text = captured["header_text"]
    sample_line = captured["append"]["serialized_sample_line"]

    assert captured["final_path"].endswith("final_output.vcf")
    assert sample_line.startswith("##SAMPLE=<")
    for fragment in [
        "ID=SAMPLE001",
        'Description="Study Cohort"',
        "Site=MGH",
        'Investigator="Doe,John"',
    ]:
        assert fragment in sample_line
        assert fragment in header_text

    assert "##analysis=pilot" in header_text
    assert captured["validated_final"] == captured["final_path"]
