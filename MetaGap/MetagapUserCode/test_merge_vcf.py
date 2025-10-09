#!/usr/bin/env python3
"""Backward compatible wrapper for the :mod:`merge_vcf` package."""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure the repository root is on ``sys.path`` so the ``merge_vcf`` package
# can be imported when this file is executed directly.
_THIS_FILE = Path(__file__).resolve()
_REPO_ROOT = _THIS_FILE.parents[2]
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

from merge_vcf.main import *  # noqa: F401,F403 - re-export legacy symbols
from merge_vcf import main as _main


if __name__ == "__main__":  # pragma: no cover - CLI entry point
    _main()
