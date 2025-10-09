import importlib.util
from pathlib import Path

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
