import gzip
import importlib
import sys
from pathlib import Path
import types

import pytest


_BASE_DIR = Path(__file__).resolve().parents[1]
if str(_BASE_DIR) not in sys.path:
    sys.path.insert(0, str(_BASE_DIR))

merge_pkg = types.ModuleType("merge_vcf")
merge_pkg.__path__ = [str(_BASE_DIR / "merge_vcf")]
try:
    merge_pkg.vcfpy = importlib.import_module("vcfpy")
    merge_pkg.VCFPY_AVAILABLE = True
except ModuleNotFoundError:  # pragma: no cover - dependency missing
    merge_pkg.vcfpy = None
    merge_pkg.VCFPY_AVAILABLE = False
try:
    merge_pkg.pysam = importlib.import_module("pysam")
    merge_pkg.PYSAM_AVAILABLE = True
except ModuleNotFoundError:  # pragma: no cover - dependency missing
    merge_pkg.pysam = None
    merge_pkg.PYSAM_AVAILABLE = False
sys.modules.setdefault("merge_vcf", merge_pkg)

logging_spec = importlib.util.spec_from_file_location(
    "merge_vcf.logging_utils", _BASE_DIR / "merge_vcf" / "logging_utils.py"
)
logging_module = importlib.util.module_from_spec(logging_spec)
sys.modules[logging_spec.name] = logging_module
logging_spec.loader.exec_module(logging_module)

merging_stub = types.ModuleType("merge_vcf.merging")
def _stub_preprocess_vcf(*_args, **_kwargs):  # pragma: no cover - helper for imports
    raise NotImplementedError("preprocess_vcf is not available in tests")


merging_stub.preprocess_vcf = _stub_preprocess_vcf
sys.modules.setdefault("merge_vcf.merging", merging_stub)

validation_module = importlib.import_module("merge_vcf.validation")
metadata_module = importlib.import_module("merge_vcf.metadata")


def _parse_info_field(info_field):
    entries = {}
    for part in info_field.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        entries[key] = value
    return entries

def test_validate_merged_vcf_missing_info_is_tolerated(tmp_path, capsys):
    module = validation_module

    vcf_content = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\trs1\tA\tC\t.\tPASS\t.
"""
    vcf_path = tmp_path / "missing_info.vcf"
    vcf_path.write_text(vcf_content)

    module.validate_merged_vcf(str(vcf_path))

    captured = capsys.readouterr()
    assert "Validation completed successfully" in captured.out


def test_validate_merged_vcf_supports_bgzipped_input(tmp_path, capsys):
    module = validation_module

    vcf_content = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t2000\trs2\tC\tT\t.\tPASS\tAC=1
"""

    gz_path = tmp_path / "anonymized.vcf.gz"

    with gzip.open(gz_path, "wt") as handle:
        handle.write(vcf_content)

    module.validate_merged_vcf(str(gz_path))

    captured = capsys.readouterr()
    assert "Validation completed successfully" in captured.out


def test_validate_merged_vcf_reports_missing_required_info_field(tmp_path):
    module = validation_module

    if not getattr(module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for merged VCF validation tests")

    vcf_content = """##fileformat=VCFv4.2
##reference=GRCh38
##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t3000\trs3\tG\tA\t.\tPASS\tAC=1
"""

    vcf_path = tmp_path / "missing_required_info.vcf"
    vcf_path.write_text(vcf_content)

    messages: list[str] = []
    original_handler = module.handle_non_critical_error
    try:
        module.handle_non_critical_error = messages.append  # type: ignore[assignment]
        module.validate_merged_vcf(str(vcf_path))
    finally:
        module.handle_non_critical_error = original_handler

    assert any("missing required INFO fields: NS" in message for message in messages)


@pytest.mark.parametrize("use_vcfpy", [True, False])
def test_recalculate_cohort_info_tags_populates_ac_an_af(
    tmp_path, monkeypatch, use_vcfpy
):
    module = metadata_module

    if not hasattr(module, "recalculate_cohort_info_tags"):
        pytest.skip("recalculate_cohort_info_tags is not available in the metadata module")

    original_flag = getattr(module, "VCFPY_AVAILABLE", False)
    if use_vcfpy:
        if not original_flag:
            pytest.skip("vcfpy dependency is required to exercise the vcfpy code path")
    else:
        monkeypatch.setattr(module, "VCFPY_AVAILABLE", False)

    vcf_content = """##fileformat=VCFv4.2
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency for each ALT allele">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
1\t100\t.\tA\tG\t.\tPASS\tAC=.;AN=.;AF=.\tGT\t0/1\t1/1
"""

    vcf_path = tmp_path / "needs_recalc.vcf"
    vcf_path.write_text(vcf_content)

    module.recalculate_cohort_info_tags(str(vcf_path), verbose=True)

    lines = [line.strip() for line in vcf_path.read_text().splitlines() if line.strip()]
    data_line = next(line for line in lines if not line.startswith("##") and not line.startswith("#CHROM"))
    fields = data_line.split("\t")
    info_entries = _parse_info_field(fields[7])

    assert info_entries["AN"] == "4"
    assert [int(value) for value in info_entries["AC"].split(",")] == [3]
    assert abs(float(info_entries["AF"].split(",")[0]) - 0.75) < 1e-6
