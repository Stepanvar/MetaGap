"""Command-line entry point for the merge_vcf package."""

from MetaGap.MetagapUserCode.merge_vcf.cli import main as _workflow_main


def main() -> None:
    """Execute the merge_vcf command-line interface."""

    _workflow_main()


if __name__ == "__main__":  # pragma: no cover - entry point
    main()
