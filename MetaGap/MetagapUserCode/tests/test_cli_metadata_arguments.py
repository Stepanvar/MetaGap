import importlib
import sys
import types
from collections import OrderedDict
from pathlib import Path
from types import SimpleNamespace

import pytest

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
if str(PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(PACKAGE_ROOT))


@pytest.fixture
def metadata_module(monkeypatch):
    """Return the metadata module with a stubbed ``vcfpy`` dependency."""

    sys.modules.pop("vcfpy", None)
    sys.modules.pop("MetagapUserCode.merge_vcf.metadata", None)

    stub = types.ModuleType("vcfpy")

    class SampleHeaderLine:
        def __init__(self, mapping=None):
            self.key = "SAMPLE"
            self.mapping = OrderedDict(mapping or {})

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


def test_parse_metadata_arguments_combines_multiple_sources(metadata_module):
    module = metadata_module

    module.template_sample_mapping = OrderedDict(
        [("ID", "TEMPLATE"), ("Description", "Template Provided")]
    )
    module.template_serialized_sample = (
        "##SAMPLE=<ID=TEMPLATE,Description=Template Provided>"
    )
    module.template_simple_lines = [
        module.vcfpy.SimpleHeaderLine("reference", "GRCh37")
    ]
    module.template_header_lines = ["##reference=GRCh37"]

    args = SimpleNamespace(
        sample_metadata_entries=["", "ID=cli-id", "Description=CLI Provided", "  "],
        meta=["Status=Legacy", "ID=Legacy"],
        header_metadata_lines=[
            "##contact=ResearchLab",
            "  ",
            "Quality=High",
            "##reference=GRCh37",
        ],
    )

    (
        sample_header_line,
        simple_header_lines,
        sample_mapping,
        sanitized_header_lines,
        serialized_sample_line,
    ) = module.parse_metadata_arguments(args, verbose=False)

    assert isinstance(sample_header_line, module.vcfpy.SampleHeaderLine)
    assert list(sample_header_line.mapping.items()) == [
        ("ID", "Legacy"),
        ("Description", "CLI Provided"),
        ("Status", "Legacy"),
    ]
    assert [
        (line.key, line.value) for line in simple_header_lines
    ] == [
        ("reference", "GRCh37"),
        ("contact", "ResearchLab"),
        ("Quality", "High"),
    ]
    assert list(sample_mapping.items()) == [
        ("ID", "Legacy"),
        ("Description", "CLI Provided"),
        ("Status", "Legacy"),
    ]
    assert sanitized_header_lines == [
        "##reference=GRCh37",
        "##contact=ResearchLab",
        "##Quality=High",
    ]
    assert serialized_sample_line == (
        '##SAMPLE=<ID=Legacy,Description="CLI Provided",Status=Legacy>'
    )


def test_parse_metadata_arguments_requires_id(metadata_module):
    module = metadata_module

    args = SimpleNamespace(
        sample_metadata_entries=["Description=Missing"],
        meta=[],
        header_metadata_lines=[],
    )

    with pytest.raises(module.ValidationError) as excinfo:
        module.parse_metadata_arguments(args, verbose=False)

    assert "non-empty ID" in str(excinfo.value)


def test_parse_metadata_arguments_rejects_malformed_header_line(metadata_module):
    module = metadata_module

    args = SimpleNamespace(
        sample_metadata_entries=["ID=S1"],
        meta=[],
        header_metadata_lines=["Invalid header"],
    )

    with pytest.raises(module.ValidationError) as excinfo:
        module.parse_metadata_arguments(args, verbose=False)

    assert "##key=value" in str(excinfo.value)
