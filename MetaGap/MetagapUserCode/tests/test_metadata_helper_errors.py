import importlib
import sys
import types
from pathlib import Path

import pytest


PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))


@pytest.fixture
def metadata_module(monkeypatch):
    """Load the metadata module with stubbed dependencies for testing."""

    sys.modules.pop("vcfpy", None)
    sys.modules.pop("MetagapUserCode.merge_vcf.metadata", None)

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

    module = importlib.import_module("MetagapUserCode.merge_vcf.metadata")

    monkeypatch.setattr(module, "template_sample_mapping", None, raising=False)
    monkeypatch.setattr(module, "template_simple_lines", [], raising=False)
    monkeypatch.setattr(module, "template_header_lines", [], raising=False)
    monkeypatch.setattr(module, "template_serialized_sample", None, raising=False)

    return module


def test_load_metadata_template_missing_file(metadata_module, tmp_path):
    missing_path = tmp_path / "missing.txt"

    with pytest.raises(metadata_module.ValidationError) as excinfo:
        metadata_module._load_metadata_template(str(missing_path))

    assert "does not exist" in str(excinfo.value)


def test_load_metadata_template_invalid_sample_line(metadata_module, tmp_path):
    template_path = tmp_path / "template.txt"
    template_path.write_text("##SAMPLE=<Invalid>")

    with pytest.raises(metadata_module.ValidationError) as excinfo:
        metadata_module._load_metadata_template(str(template_path))

    assert "Invalid SAMPLE" in str(excinfo.value)


def test_load_metadata_lines_requires_path(metadata_module):
    with pytest.raises(metadata_module.ValidationError) as excinfo:
        metadata_module.load_metadata_lines("")

    assert "must be provided" in str(excinfo.value)


def test_load_metadata_lines_rejects_invalid_line(metadata_module, tmp_path):
    metadata_path = tmp_path / "metadata.txt"
    metadata_path.write_text("Invalid line without equals\n")

    with pytest.raises(metadata_module.ValidationError) as excinfo:
        metadata_module.load_metadata_lines(str(metadata_path))

    assert "'##key=value'" in str(excinfo.value)
