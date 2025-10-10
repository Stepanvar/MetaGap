"""Shared pytest fixtures for the MetagapUserCode test suite."""

import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Dict

import pytest
from pytest import MonkeyPatch

# Cache loaded modules across tests (especially useful for session-scoped fixtures)
_MODULE_CACHE: Dict[str, ModuleType] = {}


def _load_script_module(base_dir: Path) -> ModuleType:
    """
    Load a standalone user script (e.g., test_merge_vcf.py) as a module.

    This uses importlib to import the file directly by path so tests can call its
    functions without needing it to be on sys.path.
    """
    module_path = base_dir / "test_merge_vcf.py"
    spec = importlib.util.spec_from_file_location("user_test_merge_vcf", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_cli_module(base_dir: Path, monkeypatch: pytest.MonkeyPatch | None = None) -> ModuleType:
    """
    Load the package-style CLI: MetagapUserCode.merge_vcf.cli

    We create a synthetic top-level package (MetagapUserCode) rooted at base_dir so that
    relative imports within the package resolve as expected during tests.
    """
    package_name = "MetagapUserCode"
    package_module = ModuleType(package_name)
    package_module.__path__ = [str(base_dir)]
    if monkeypatch is not None:
        monkeypatch.setitem(sys.modules, package_name, package_module)
    else:
        sys.modules.setdefault(package_name, package_module)

    # Load the package 'MetagapUserCode.merge_vcf'
    merge_pkg_name = f"{package_name}.merge_vcf"
    merge_pkg_path = base_dir / "merge_vcf" / "__init__.py"
    merge_spec = importlib.util.spec_from_file_location(
        merge_pkg_name,
        merge_pkg_path,
        submodule_search_locations=[str(base_dir / "merge_vcf")],
    )
    merge_module = importlib.util.module_from_spec(merge_spec)
    assert merge_spec is not None and merge_spec.loader is not None
    if monkeypatch is not None:
        monkeypatch.setitem(sys.modules, merge_pkg_name, merge_module)
    else:
        sys.modules[merge_pkg_name] = merge_module
    merge_spec.loader.exec_module(merge_module)

    # Load the CLI module 'MetagapUserCode.merge_vcf.cli'
    cli_name = f"{merge_pkg_name}.cli"
    cli_path = base_dir / "merge_vcf" / "cli.py"
    cli_spec = importlib.util.spec_from_file_location(cli_name, cli_path)
    cli_module = importlib.util.module_from_spec(cli_spec)
    assert cli_spec is not None and cli_spec.loader is not None
    if monkeypatch is not None:
        monkeypatch.setitem(sys.modules, cli_name, cli_module)
    else:
        sys.modules[cli_name] = cli_module
    cli_spec.loader.exec_module(cli_module)

    return cli_module


@pytest.fixture
def cli_module(monkeypatch):
    """Return the ``MetagapUserCode.merge_vcf.cli`` module for testing."""
    base_dir = Path(__file__).resolve().parents[1]
    return _load_cli_module(base_dir, monkeypatch=monkeypatch)


@pytest.fixture(scope="session")
def merge_script_module(request) -> ModuleType:
    """
    Return either the user's standalone script module or the CLI module.

    Usage:
        @pytest.mark.parametrize("merge_script_module", ["script"], indirect=True)
        def test_something(merge_script_module): ...

        @pytest.mark.parametrize("merge_script_module", ["cli"], indirect=True)
        def test_cli(merge_script_module): ...
    """
    variant = getattr(request, "param", "script")  # "script" or "cli"

    cached = _MODULE_CACHE.get(variant)
    if cached is not None:
        return cached

    base_dir = Path(__file__).resolve().parents[1]
    script_path = base_dir / "test_merge_vcf.py"

    if variant == "cli" or not script_path.exists():
        mp = MonkeyPatch()
        request.addfinalizer(mp.undo)
        module = _load_cli_module(base_dir, monkeypatch=mp)
    else:
        module = _load_script_module(base_dir)

    _MODULE_CACHE[variant] = module
    return module
