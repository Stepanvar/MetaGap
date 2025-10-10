import importlib.util
import datetime as _dt
import gzip
import sys
from pathlib import Path
from types import ModuleType

import pytest


def _write_sample_vcf(path: Path, alt_order: list[str]) -> None:
    lines = [
        "##fileformat=VCFv4.2",
        "##reference=GRCh38",
        "##contig=<ID=chr1,length=100>",
        "##contig=<ID=chr2,length=200>",
        f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "chr1\t1\t.\tA\tC\t50\tPASS\t.",
        "chr2\t5\t.\tT\t{alt}\t50\tPASS\t.".format(alt=",".join(alt_order)),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _load_cli_module(patcher: pytest.MonkeyPatch):
    base_dir = Path(__file__).resolve().parents[1]

    package_name = "MetagapUserCode"
    package_module = ModuleType(package_name)
    package_module.__path__ = [str(base_dir)]
    patcher.setitem(sys.modules, package_name, package_module)

    merge_pkg_name = f"{package_name}.merge_vcf"
    merge_pkg_path = base_dir / "merge_vcf" / "__init__.py"
    merge_spec = importlib.util.spec_from_file_location(
        merge_pkg_name,
        merge_pkg_path,
        submodule_search_locations=[str(base_dir / "merge_vcf")],
    )
    merge_module = importlib.util.module_from_spec(merge_spec)
    assert merge_spec is not None and merge_spec.loader is not None
    patcher.setitem(sys.modules, merge_pkg_name, merge_module)
    merge_spec.loader.exec_module(merge_module)

    cli_name = f"{merge_pkg_name}.cli"
    cli_path = base_dir / "merge_vcf" / "cli.py"
    cli_spec = importlib.util.spec_from_file_location(cli_name, cli_path)
    cli_module = importlib.util.module_from_spec(cli_spec)
    assert cli_spec is not None and cli_spec.loader is not None
    patcher.setitem(sys.modules, cli_name, cli_module)
    cli_spec.loader.exec_module(cli_module)

    return cli_module


class _SequenceDateTime(_dt.datetime):
    _counter = 0

    @classmethod
    def now(cls):
        base = _dt.datetime(2024, 1, 1, 0, 0, 0)
        value = base + _dt.timedelta(seconds=cls._counter)
        cls._counter += 1
        return value


def _install_fake_subprocess(patcher, module):
    def fake_run(cmd, *args, **kwargs):
        if cmd[:2] == ["bcftools", "merge"]:
            output_path = Path(cmd[cmd.index("-o") + 1])
            output_path.write_text("", encoding="utf-8")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:1] == ["bgzip"]:
            target = Path(cmd[-1])
            gz_path = target.with_suffix(target.suffix + ".gz") if target.suffix else Path(
                str(target) + ".gz"
            )
            with open(target, "rb") as src, gzip.open(gz_path, "wb") as dest:
                dest.write(src.read())
            target.unlink()
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        if cmd[:1] == ["tabix"]:
            index_path = Path(cmd[-1] + ".tbi")
            index_path.write_bytes(b"")
            return SimpleNamespace(returncode=0, stdout="", stderr="")
        raise AssertionError(f"Unexpected command: {cmd}")

    from types import SimpleNamespace
    
    merging_module = sys.modules[module.merge_vcfs.__module__]
    patcher.setattr(merging_module.subprocess, "run", fake_run)


def test_merge_outputs_are_deterministic(tmp_path):
    sample_a = tmp_path / "sample_a.vcf"
    sample_b = tmp_path / "sample_b.vcf"
    _write_sample_vcf(sample_a, ["G", "T"])
    _write_sample_vcf(sample_b, ["T", "G"])

    patcher = pytest.MonkeyPatch()
    first_contents = ""
    second_contents = ""
    try:
        module = _load_cli_module(patcher)

        if not getattr(module, "VCFPY_AVAILABLE", True):
            pytest.skip("vcfpy dependency is required for merge ordering tests")

        _install_fake_subprocess(patcher, module)
        patcher.setattr(module.datetime, "datetime", _SequenceDateTime)
        _SequenceDateTime._counter = 0

        output_dir = tmp_path / "out"
        output_dir.mkdir()

        first_path = Path(
            module.merge_vcfs(
                [str(sample_a), str(sample_b)],
                output_dir=str(output_dir),
                verbose=False,
                sample_order=[],
                an_threshold=None,
            )
        )

        second_path = Path(
            module.merge_vcfs(
                [str(sample_b), str(sample_a)],
                output_dir=str(output_dir),
                verbose=False,
                sample_order=[],
                an_threshold=None,
            )
        )

        with gzip.open(first_path, "rt", encoding="utf-8") as handle:
            first_contents = handle.read()
        with gzip.open(second_path, "rt", encoding="utf-8") as handle:
            second_contents = handle.read()

    finally:
        patcher.undo()

    assert first_contents == second_contents

    contig_lines = [
        line for line in first_contents.splitlines() if line.startswith("##contig=")
    ]
    assert contig_lines == [
        "##contig=<ID=chr1,length=100>",
        "##contig=<ID=chr2,length=200>",
    ]

    data_lines = [line for line in first_contents.splitlines() if not line.startswith("#")]
    assert any(line.split("\t")[0:2] == ["chr1", "1"] for line in data_lines)
    chr2_line = next(line for line in data_lines if line.startswith("chr2"))
    fields = chr2_line.split("\t")
    assert fields[4] == "G,T"
