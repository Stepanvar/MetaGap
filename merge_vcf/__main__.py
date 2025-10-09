"""Command-line entry point for the merge_vcf package."""

from . import main as _main


def main() -> None:
    """Execute the merge_vcf command-line interface."""
    _main()


if __name__ == "__main__":  # pragma: no cover - entry point
    main()
