"""Shared pytest fixtures for the MetagapUserCode test suite."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

import pytest


@pytest.fixture
def cli_module(monkeypatch):
    """Return the ``MetagapUserCode.merge_vcf.cli`` module for testing."""

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
