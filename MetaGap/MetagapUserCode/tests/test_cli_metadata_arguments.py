import sys
from pathlib import Path
import importlib
import types

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))


def _install_vcfpy_stub():
    if "vcfpy" in sys.modules:
        return

    stub = types.ModuleType("vcfpy")

    class SampleHeaderLine:
        def __init__(self, mapping=None):
            self.key = "SAMPLE"
            self.mapping = dict(mapping or {})

        @classmethod
        def from_mapping(cls, mapping):
            return cls(mapping=mapping)

    class SimpleHeaderLine:
        def __init__(self, key, value):
            self.key = key
            self.value = value

    stub.SampleHeaderLine = SampleHeaderLine
    stub.SimpleHeaderLine = SimpleHeaderLine
    stub.header = types.SimpleNamespace(InfoHeaderLine=type("InfoHeaderLine", (), {}))

    sys.modules["vcfpy"] = stub


_install_vcfpy_stub()

cli = importlib.import_module("MetagapUserCode.merge_vcf.cli")
metadata = importlib.import_module("MetagapUserCode.merge_vcf.metadata")


def test_parse_arguments_honors_legacy_meta_flag(monkeypatch):
    argv = [
        "merge_vcf",
        "--input-dir",
        "/tmp/input",
        "--meta",
        "ID=legacy",
        "--meta",
        "Sex=Female",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    args = cli.parse_arguments()
    (
        sample_header_line,
        simple_header_lines,
        sample_metadata_entries,
        sanitized_header_lines,
        serialized_sample_line,
    ) = metadata.parse_metadata_arguments(args, verbose=False)

    assert serialized_sample_line == "##SAMPLE=<ID=legacy,Sex=Female>"
    assert list(sample_metadata_entries.items()) == [
        ("ID", "legacy"),
        ("Sex", "Female"),
    ]
    assert sample_header_line is not None
    assert simple_header_lines == []
    assert sanitized_header_lines == []
