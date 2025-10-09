"""Command line interface for the VCF merging workflow."""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import os
import sys
from typing import Any


def parse_arguments(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse command line arguments for the VCF merging workflow."""

    parser = argparse.ArgumentParser(
        description="Consolidate multiple VCF files into a single merged VCF file."
    )
    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing the VCF files to merge.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        "--output",
        dest="output_dir",
        help="Directory for the merged VCF. Defaults to the input directory when omitted.",
    )
    parser.add_argument(
        "--ref",
        help="Expected reference genome build. When omitted the script attempts to auto-detect it.",
    )
    parser.add_argument(
        "--vcf-version",
        help="Expected VCF version (e.g., 4.2). When omitted the script attempts to auto-detect it.",
    )
    parser.add_argument(
        "--meta",
        action="append",
        default=[],
        metavar="KEY=VALUE",
        help=(
            "Sample metadata in KEY=VALUE form. Repeat for multiple keys. "
            "An ID entry is required when metadata is provided."
        ),
    )
    parser.add_argument(
        "--allow-gvcf",
        action="store_true",
        help="Allow input files that contain gVCF annotations such as <NON_REF> ALT alleles.",
    )
    parser.add_argument(
        "--sample-metadata",
        dest="sample_metadata_entries",
        action="append",
        metavar="KEY=VALUE",
        default=[],
        help="Add a key/value pair to the ##SAMPLE metadata line. Provide at least ID=... to emit the line.",
    )
    parser.add_argument(
        "--header-metadata",
        dest="header_metadata_lines",
        action="append",
        metavar="LINE",
        default=[],
        help="Add an arbitrary metadata header line (##key=value). The '##' prefix is optional.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Enable verbose console logging in addition to the log file.",
    )
    return parser.parse_args(argv)


def _load_workflow_module():
    candidates = ("MetagapUserCode.test_merge_vcf", "test_merge_vcf")
    for name in candidates:
        module = sys.modules.get(name)
        if module is not None:
            return module

    workflow_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_merge_vcf.py")
    normalized_workflow_path = os.path.normpath(os.path.abspath(workflow_path))
    for module in list(sys.modules.values()):
        module_path = getattr(module, "__file__", None)
        if module_path and os.path.normpath(os.path.abspath(module_path)) == normalized_workflow_path:
            return module

    for name in candidates:
        try:
            return importlib.import_module(name)
        except ModuleNotFoundError:
            continue

    spec = importlib.util.spec_from_file_location("test_merge_vcf", workflow_path)
    if spec is None or spec.loader is None:
        raise ImportError("Unable to locate workflow module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sys.modules.setdefault("test_merge_vcf", module)
    return module


def main(args: Any | None = None, workflow_module: Any | None = None):
    """Entry point for the command line interface."""

    if workflow_module is None:
        workflow_module = _load_workflow_module()
    if args is None:
        parse_args = getattr(workflow_module, "parse_arguments", parse_arguments)
        args = parse_args()
    return workflow_module.run_workflow(args)


if __name__ == "__main__":  # pragma: no cover - manual execution helper
    main()
