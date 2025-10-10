"""Integration-style test covering merge cleanup for numerous shards."""

from __future__ import annotations

import gzip
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from MetagapUserCode.merge_vcf import merging
from MetagapUserCode.tests.test_large_merge_workflow import (
    _FakeVcfpyModule,
    _VCF_STORAGE,
    _generate_demo_vcf,
)


@pytest.mark.usefixtures("monkeypatch")
def test_merge_handles_many_shards_without_leaving_artifacts(tmp_path, monkeypatch):
    """Merge 50+ shards while ensuring cleanup and sane resource usage."""

    # Ensure a clean slate for the in-memory VCF stub used by merge helpers.
    _VCF_STORAGE.clear()

    # Patch the vcfpy dependency with the rich stub from the large workflow tests.
    fake_vcfpy = _FakeVcfpyModule()

    class _FilterHeaderLine:
        def __init__(self, mapping):
            self.key = "FILTER"
            self.mapping = dict(mapping)
            self.id = self.mapping.get("ID")

        def copy(self):  # pragma: no cover - simple data container
            return _FilterHeaderLine(self.mapping)

        def to_line(self):  # pragma: no cover - serialization helper
            parts = ",".join(f"{k}={v}" for k, v in self.mapping.items())
            return f"##FILTER=<{parts}>"

    fake_vcfpy.header.FilterHeaderLine = _FilterHeaderLine

    class _ContigHeaderLine:
        def __init__(self, mapping):
            self.key = "CONTIG"
            self.mapping = dict(mapping)
            self.id = self.mapping.get("ID")

        def copy(self):  # pragma: no cover - simple data container
            return _ContigHeaderLine(self.mapping)

        def to_line(self):  # pragma: no cover - serialization helper
            parts = ",".join(f"{k}={v}" for k, v in self.mapping.items())
            return f"##contig=<{parts}>"

    fake_vcfpy.header.ContigHeaderLine = _ContigHeaderLine

    original_reader_cls = fake_vcfpy.Reader
    original_call_cls = fake_vcfpy.Call
    original_copy = original_call_cls.copy

    def _copy_with_sample(self):  # pragma: no cover - helper to retain sample attribute
        duplicate = original_copy(self)
        setattr(duplicate, "sample", getattr(duplicate, "name", None))
        return duplicate

    original_call_cls.copy = _copy_with_sample

    class _Reader(original_reader_cls):  # type: ignore[misc]
        def __init__(self, path):
            super().__init__(path)
            for rec in getattr(self, "_records", []):
                for call in getattr(rec, "calls", []):
                    if not hasattr(call, "sample"):
                        setattr(call, "sample", getattr(call, "name", None))

        @classmethod
        def from_path(cls, path, *args, **kwargs):
            return cls(path)

    fake_vcfpy.Reader = _Reader

    original_writer_cls = fake_vcfpy.Writer

    class _Writer(original_writer_cls):  # type: ignore[misc]
        def write_record(self, record):
            filt = getattr(record, "FILTER", None)
            if isinstance(filt, list):
                record.FILTER = ";".join(filt) if filt else "PASS"
            super().write_record(record)

    fake_vcfpy.Writer = _Writer

    class _Call(fake_vcfpy.Call):  # type: ignore[misc]
        def __init__(self, name, data):
            super().__init__(name, data)
            self.sample = self.name

        def copy(self):  # pragma: no cover - simple delegation
            duplicate = super().copy()
            setattr(duplicate, "sample", getattr(duplicate, "name", None))
            return duplicate

    fake_vcfpy.Call = _Call
    monkeypatch.setattr(merging, "vcfpy", fake_vcfpy, raising=False)
    monkeypatch.setattr(merging, "VCFPY_AVAILABLE", True, raising=False)

    # Patch metadata/header helpers so the stub header survives untouched.
    monkeypatch.setattr(merging, "apply_metadata_to_header", lambda header, **_: header)
    monkeypatch.setattr(merging, "_ensure_info_header_lines", lambda header: None)

    # Capture verbose logging to assert progress reporting occurred.
    captured_logs: list[str] = []

    def _capture(message: str, verbose: bool = False):
        if verbose:
            captured_logs.append(str(message))

    monkeypatch.setattr(merging, "log_message", _capture)

    # Provide lightweight pysam replacements for compression and indexing.
    def _fake_tabix_compress(src: str, dest: str, force: bool = True):
        with open(src, "rb") as read_handle, gzip.open(dest, "wb") as write_handle:
            write_handle.write(read_handle.read())

    def _fake_tabix_index(path: str, preset: str = "vcf", force: bool = True):
        index_path = Path(path + ".tbi")
        index_path.write_bytes(b"")

    monkeypatch.setattr(
        merging,
        "pysam",
        SimpleNamespace(tabix_compress=_fake_tabix_compress, tabix_index=_fake_tabix_index),
        raising=False,
    )

    # Generate a cohort of small yet valid shards.
    input_dir = tmp_path / "inputs"
    output_dir = tmp_path / "outputs"
    input_dir.mkdir()
    output_dir.mkdir()

    shard_paths: list[str] = []
    for index in range(55):
        sample_name = f"SAMPLE_{index + 1:03d}"
        genotype = "0/1" if index % 2 == 0 else "1/1"
        shard_path = input_dir / f"{sample_name}.vcf"
        _generate_demo_vcf(shard_path, sample_name, genotype)
        shard_paths.append(str(shard_path))

    # Track runtime and resident memory before kicking off the merge.
    start = time.monotonic()
    rss_before = None
    try:
        import psutil  # type: ignore

        rss_before = psutil.Process().memory_info().rss
    except Exception:  # pragma: no cover - psutil may be absent
        psutil = None  # type: ignore

    merged_path = merging.merge_vcfs(
        shard_paths,
        str(output_dir),
        verbose=True,
        qual_threshold=None,
        an_threshold=None,
    )

    duration = time.monotonic() - start
    rss_after = None
    if "psutil" in locals() and psutil is not None:  # pragma: no branch
        rss_after = psutil.Process().memory_info().rss

    # The merge completed successfully and produced compressed + indexed output.
    final_path = Path(merged_path)
    assert final_path.exists()
    assert final_path.suffix == ".gz"

    produced = sorted(p.name for p in output_dir.iterdir())
    expected = sorted([final_path.name, final_path.name + ".tbi"])
    assert produced == expected

    # Intermediate artifacts should have been removed by cleanup logic.
    assert not list(output_dir.glob("*.anon"))
    assert not list(output_dir.glob("*.filtered"))
    assert not list(output_dir.glob("*.tmp"))

    # Verbose logging captured progress updates, proving verbose mode is honoured.
    assert any("Discovered" in entry for entry in captured_logs)
    assert any("Compressing" in entry for entry in captured_logs)

    # Resource usage stayed within reasonable bounds for the synthetic dataset.
    assert duration < 10, f"merge_vcfs unexpectedly slow: {duration:.2f}s"
    if rss_before is not None and rss_after is not None:
        assert rss_after - rss_before < 200 * 1024 * 1024

