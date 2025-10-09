import importlib.util
import sys
import types
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
