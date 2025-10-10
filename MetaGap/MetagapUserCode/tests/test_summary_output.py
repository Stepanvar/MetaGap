import importlib.util
import sys
import types
from collections import OrderedDict
from types import SimpleNamespace
from pathlib import Path


MODULE_PATH = Path(__file__).resolve().parents[2] / "MetagapUserCode" / "test_merge_vcf.py"


def _install_vcfpy_stub():
    if "vcfpy" in sys.modules:
        return

    stub = types.ModuleType("vcfpy")

    class _BaseHeaderLine:
        def __init__(self, key=None, value=None, mapping=None):
            self.key = key
            self.value = value
            self.mapping = dict(mapping or {})
            self.id = self.mapping.get("ID")

    class SampleHeaderLine(_BaseHeaderLine):
        @classmethod
        def from_mapping(cls, mapping):
            return cls(key="SAMPLE", mapping=mapping)

    class SimpleHeaderLine(_BaseHeaderLine):
        def __init__(self, key, value):
            super().__init__(key=key, value=value)

    header_namespace = types.SimpleNamespace(
        InfoHeaderLine=type("InfoHeaderLine", (), {}),
        FilterHeaderLine=type("FilterHeaderLine", (), {}),
        ContigHeaderLine=type("ContigHeaderLine", (), {}),
        SampleHeaderLine=SampleHeaderLine,
        FormatHeaderLine=type("FormatHeaderLine", (), {}),
    )

    class Reader:
        def __init__(self):
            self.header = types.SimpleNamespace(
                lines=[], samples=types.SimpleNamespace(names=[])
            )

        @classmethod
        def from_path(cls, *args, **kwargs):
            return cls()

        @classmethod
        def from_stream(cls, *args, **kwargs):
            return cls()

        def __iter__(self):
            return iter([])

        def close(self):
            pass

    class Writer:
        @classmethod
        def from_path(cls, *args, **kwargs):
            return cls()

        def write_record(self, *args, **kwargs):
            pass

        def close(self):
            pass

    class Call:
        def __init__(self, name, data):
            self.name = name
            self.data = data

    stub.SampleHeaderLine = SampleHeaderLine
    stub.SimpleHeaderLine = SimpleHeaderLine
    stub.Reader = Reader
    stub.Writer = Writer
    stub.Call = Call
    stub.header = header_namespace

    sys.modules["vcfpy"] = stub


def load_user_module():
    _install_vcfpy_stub()
    spec = importlib.util.spec_from_file_location("user_test_merge_vcf", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_summarize_produced_vcfs_prefers_sample_files(tmp_path):
    module = load_user_module()

    output_dir = tmp_path / "outputs"
    output_dir.mkdir()
    (output_dir / "SAMPLE_0002.vcf").touch()
    (output_dir / "SAMPLE_0001.vcf.gz").touch()

    fallback = output_dir / "merged_output.vcf.gz"
    fallback.touch()

    directory, filename, count = module.summarize_produced_vcfs(
        str(output_dir), str(fallback)
    )

    assert Path(directory) == output_dir.resolve()
    assert filename == "SAMPLE_0001.vcf.gz"
    assert count == 2


def test_main_emits_compact_summary(tmp_path, monkeypatch, capsys):
    module = load_user_module()

    input_dir = tmp_path / "inputs"
    output_dir = tmp_path / "outputs"
    input_dir.mkdir()
    output_dir.mkdir()

    # Create files representing per-sample cohort shards.
    (output_dir / "SAMPLE_0001.vcf.gz").touch()
    (output_dir / "SAMPLE_0002.vcf").touch()

    final_vcf = output_dir / "merged_output.vcf.gz"
    final_vcf.touch()

    args = SimpleNamespace(
        input_dir=str(input_dir),
        output_dir=str(output_dir),
        ref="GRCh38",
        vcf_version="4.2",
        meta=[],
        sample_metadata_entries=[],
        header_metadata_lines=[],
        verbose=False,
        allow_gvcf=False,
    )

    monkeypatch.setattr(module, "parse_arguments", lambda: args)
    monkeypatch.setattr(
        module,
        "validate_all_vcfs",
        lambda *unused, **unused_kwargs: ([str(input_dir / "SAMPLE_0001.vcf")], ["SAMPLE_0001"]),
    )
    monkeypatch.setattr(module, "merge_vcfs", lambda *unused, **unused_kwargs: str(final_vcf))
    monkeypatch.setattr(
        module,
        "append_metadata_to_merged_vcf",
        lambda *unused, **unused_kwargs: str(final_vcf),
    )
    monkeypatch.setattr(module, "validate_merged_vcf", lambda *args, **kwargs: None)

    module.main()

    captured = capsys.readouterr()
    summary_line = captured.out.strip().splitlines()[-1]
    expected_path = output_dir.resolve() / "SAMPLE_0001.vcf.gz"
    assert summary_line == f"Wrote: {expected_path} x 2"


def test_parse_metadata_arguments_merges_template(tmp_path):
    module = load_user_module()

    template_path = tmp_path / "template_header.txt"
    template_path.write_text(
        "\n".join(
            [
                "##reference=GRCh38",
                "##SAMPLE=<ID=TEMPLATE123,Description=Template Provided>",
                "##contact=ResearchLab",
                "",
            ]
        )
    )

    args = SimpleNamespace(
        sample_metadata_entries=[],
        header_metadata_lines=[],
        metadata_template_path=str(template_path),
    )

    (
        sample_header_line,
        simple_header_lines,
        sample_mapping,
        sanitized_header_lines,
        serialized_sample_line,
    ) = module.parse_metadata_arguments(args, verbose=False)

    assert sample_header_line is not None
    assert serialized_sample_line and serialized_sample_line.startswith("##SAMPLE=<")
    assert sample_mapping == OrderedDict(
        [("ID", "TEMPLATE123"), ("Description", "Template Provided")]
    )

    assert "##reference=GRCh38" in sanitized_header_lines
    assert "##contact=ResearchLab" in sanitized_header_lines
    simple_pairs = {(line.key, line.value) for line in simple_header_lines}
    assert ("reference", "GRCh38") in simple_pairs
    assert ("contact", "ResearchLab") in simple_pairs


def test_parse_metadata_arguments_template_allows_cli_overrides(tmp_path):
    module = load_user_module()

    template_path = tmp_path / "template_header.txt"
    template_path.write_text(
        "\n".join(
            [
                "##project=Genome",
                "##SAMPLE=<ID=TEMPLATE123,Description=Template Provided>",
            ]
        )
    )

    args = SimpleNamespace(
        sample_metadata_entries=["Description=CLI Override", "NewField=Extra"],
        header_metadata_lines=["##analysis=complete"],
        metadata_template_path=str(template_path),
    )

    (
        _sample_header_line,
        _simple_header_lines,
        sample_mapping,
        sanitized_header_lines,
        serialized_sample_line,
    ) = module.parse_metadata_arguments(args, verbose=False)

    assert serialized_sample_line and "NewField=Extra" in serialized_sample_line
    assert sample_mapping["ID"] == "TEMPLATE123"
    assert sample_mapping["Description"] == "CLI Override"
    assert sample_mapping["NewField"] == "Extra"
    assert "##project=Genome" in sanitized_header_lines
    assert "##analysis=complete" in sanitized_header_lines
