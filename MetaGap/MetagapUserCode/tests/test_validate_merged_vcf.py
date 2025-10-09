import gzip
import importlib.util
from pathlib import Path


def _parse_info_field(info_field):
    entries = {}
    for part in info_field.split(";"):
        if "=" not in part:
            continue
        key, value = part.split("=", 1)
        entries[key] = value
    return entries

MODULE_PATH = Path(__file__).resolve().parents[2] / "MetagapUserCode" / "test_merge_vcf.py"


def load_user_module():
    spec = importlib.util.spec_from_file_location("user_test_merge_vcf", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_validate_merged_vcf_missing_info_is_tolerated(tmp_path, capsys):
    module = load_user_module()

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
    module = load_user_module()

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


def test_recalculate_cohort_info_tags_populates_ac_an_af(tmp_path):
    module = load_user_module()

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
