"""Compatibility shim that forwards to the canonical MetaGap workflow."""

from MetaGap.MetagapUserCode.merge_vcf.cli import main as _workflow_main


def main() -> None:
    """Delegate execution to the shared MetaGap merge VCF workflow."""

    _workflow_main()


if __name__ == "__main__":  # pragma: no cover - manual execution helper
    main()
