"""MetaGap VCF merging utilities with built-in safety checks.

This package exposes the Python helpers used by MetaGap to merge and validate
Variant Call Format (VCF) files. Importing the package immediately verifies
that the required runtime dependencies :mod:`vcfpy` and :mod:`pysam` are
available so that later operations can rely on them without deferred import
errors. In addition, the package performs an abstract syntax tree (AST) audit
of its modules to ensure subprocess calls never invoke external utilities such
as ``bcftools``, ``bgzip``, or ``tabix``. This guard enforces the use of the
native Python implementations bundled with MetaGap when handling VCF data.
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Iterable


def _import_dependency(name: str):
    try:
        module = __import__(name)
    except ImportError as exc:  # pragma: no cover - exercised when dependency missing
        raise ModuleNotFoundError(
            f"The '{name}' package is required for MetaGap VCF utilities. "
            "Please install it with 'pip install {name}'."
        ) from exc
    return module


vcfpy = _import_dependency("vcfpy")
pysam = _import_dependency("pysam")

VCFPY_AVAILABLE = True
PYSAM_AVAILABLE = True


def _read_module_source(relative_path: str) -> str:
    package_dir = Path(__file__).resolve().parent
    module_path = package_dir / relative_path
    if not module_path.exists():  # pragma: no cover - defensive
        return ""
    return module_path.read_text(encoding="utf-8")


def _iter_forbidden_commands(tree: ast.AST, forbidden: Iterable[str]):
    forbidden_set = set(forbidden)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        func = node.func
        if not isinstance(func, ast.Attribute):
            continue
        if not isinstance(func.value, ast.Name):
            continue
        if func.value.id != "subprocess":
            continue
        if not node.args:
            continue
        first = node.args[0]
        values: list[str] = []
        if isinstance(first, ast.List):
            for elt in first.elts:
                if isinstance(elt, ast.Constant) and isinstance(elt.value, str):
                    values.append(elt.value)
        elif isinstance(first, ast.Constant) and isinstance(first.value, str):
            values.append(first.value)
        if forbidden_set.intersection(values):
            yield values


def _ensure_no_forbidden_subprocess_usage() -> None:
    forbidden = {"bcftools", "bgzip", "tabix"}
    for module_name in ("merging.py", "metadata.py"):
        source = _read_module_source(module_name)
        if not source:
            continue
        try:
            tree = ast.parse(source)
        except SyntaxError:  # pragma: no cover - defensive
            continue
        for values in _iter_forbidden_commands(tree, forbidden):
            joined = ", ".join(values)
            raise RuntimeError(
                "Forbidden subprocess command detected in MetaGap tooling: "
                f"{joined}. Please use the native Python implementations."
            )


_ensure_no_forbidden_subprocess_usage()


__all__ = ["vcfpy", "pysam", "VCFPY_AVAILABLE", "PYSAM_AVAILABLE"]

