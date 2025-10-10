"""Convenience entry point for running the MetaGap VCF merging CLI."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional, Sequence

_PACKAGE_ROOT = Path(__file__).resolve().parents[1]
if str(_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(_PACKAGE_ROOT))

from MetaGap.MetagapUserCode.merge_vcf.cli import main as _cli_main


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Dispatch to the production CLI that powers the merge workflow.

    The function mirrors :func:`MetaGap.MetagapUserCode.merge_vcf.cli.main` while
    accepting an optional *argv* sequence so callers can invoke the workflow
    programmatically. When *argv* is ``None`` the arguments are read directly
    from :data:`sys.argv`.
    """

    if argv is None:
        _cli_main()
        return

    original_argv = sys.argv
    try:
        sys.argv = [original_argv[0], *argv]
        _cli_main()
    finally:
        sys.argv = original_argv


__all__ = ["main"]


if __name__ == "__main__":  # pragma: no cover - manual execution helper
    main()
