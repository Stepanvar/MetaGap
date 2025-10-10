import importlib.util
import sys
from pathlib import Path
from types import ModuleType
from typing import Dict

import pytest


_MODULE_CACHE: Dict[str, ModuleType] = {}


def _load_script_module(base_dir: Path) -> ModuleType:
    module_path = base_dir / "test_merge_vcf.py"
    spec = importlib.util.spec_from_file_location("user_test_merge_vcf", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _load_cli_module(base_dir: Path) -> ModuleType:
    package_name = "MetagapUserCode"
    package_module = ModuleType(package_name)
    package_module.__path__ = [str(base_dir)]
    sys.modules.setdefault(package_name, package_module)

    merge_pkg_name = f"{package_name}.merge_vcf"
    merge_pkg_path = base_dir / "merge_vcf" / "__init__.py"
    merge_spec = importlib.util.spec_from_file_location(
        merge_pkg_name,
        merge_pkg_path,
        submodule_search_locations=[str(base_dir / "merge_vcf")],
    )
    merge_module = importlib.util.module_from_spec(merge_spec)
    assert merge_spec.loader is not None
    sys.modules[merge_pkg_name] = merge_module
    merge_spec.loader.exec_module(merge_module)

    cli_name = f"{merge_pkg_name}.cli"
    cli_path = base_dir / "merge_vcf" / "cli.py"
    cli_spec = importlib.util.spec_from_file_location(cli_name, cli_path)
    cli_module = importlib.util.module_from_spec(cli_spec)
    assert cli_spec.loader is not None
    sys.modules[cli_name] = cli_module
    cli_spec.loader.exec_module(cli_module)
    return cli_module


@pytest.fixture(scope="session")
def merge_script_module(request) -> ModuleType:
    """Return the user merge script or CLI module, importing it once per test session."""

    variant = getattr(request, "param", "script")
    cached = _MODULE_CACHE.get(variant)
    if cached is not None:
        return cached

    base_dir = Path(__file__).resolve().parents[1]
    if variant == "cli":
        module = _load_cli_module(base_dir)
    else:
        module = _load_script_module(base_dir)

    _MODULE_CACHE[variant] = module
    return module
