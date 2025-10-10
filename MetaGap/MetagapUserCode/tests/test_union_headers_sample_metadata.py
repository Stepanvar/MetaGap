from pathlib import Path

import pytest

@pytest.fixture(autouse=True)
def _skip_without_vcfpy(merge_script_module):
    if not getattr(merge_script_module, "VCFPY_AVAILABLE", True):
        pytest.skip("vcfpy dependency is required for header union tests")


def _write_vcf(path, sample_line):
    content = """##fileformat=VCFv4.2
##reference=GRCh38
{sample}
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
1\t1000\trs1\tA\tC\t.\tPASS\t.\tGT\t0/1
""".format(sample=sample_line)
    path.write_text(content)


def _write_vcf_with_samples(path, samples):
    sample_columns = "\t".join(samples)
    sample_calls = "\t".join("0/1" for _ in samples)
    content = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_columns}
1\t1000\trs1\tA\tC\t.\tPASS\t.\tGT\t{sample_calls}
""".format(sample_columns=sample_columns, sample_calls=sample_calls)
    path.write_text(content)


def test_union_headers_merges_sample_metadata(tmp_path, merge_script_module):
    module = merge_script_module

    sample_one = (
        '##SAMPLE=<ID=S1,Description="Primary sample",Meta="{\\"foo\\": \\"bar,baz\\"}">' 
    )
    sample_two = '##SAMPLE=<ID=S1,Extra=42>'

    vcf_one = tmp_path / "one.vcf"
    vcf_two = tmp_path / "two.vcf"
    _write_vcf(vcf_one, sample_one)
    _write_vcf(vcf_two, sample_two)

    header = module.union_headers([str(vcf_one), str(vcf_two)])

    sample_lines = [
        line for line in header.lines if getattr(line, "key", None) == "SAMPLE"
    ]
    assert len(sample_lines) == 1

    mapping = sample_lines[0].mapping
    assert mapping["ID"] == "S1"
    assert mapping["Description"] == "Primary sample"
    assert mapping["Meta"] == '{"foo": "bar,baz"}'
    assert mapping["Extra"] == "42"


def test_union_headers_preserves_sample_order(tmp_path, merge_script_module):
    module = merge_script_module

    vcf_one = tmp_path / "one.vcf"
    vcf_two = tmp_path / "two.vcf"

    _write_vcf_with_samples(vcf_one, ["Alpha", "Beta"])
    _write_vcf_with_samples(vcf_two, ["Beta", "Gamma", "Alpha"])

    header = module.union_headers([str(vcf_one), str(vcf_two)])

    assert header.samples.names == ["Alpha", "Beta", "Gamma"]
