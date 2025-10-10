import gzip
import io
from pathlib import Path

import pytest

SAMPLE_HEADER = """##fileformat=VCFv4.2
##reference=GRCh38
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
"""
SAMPLE_BODY_WITH_FORMAT = "1\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=2;AF=0.5\tGT\t0/1\n"
SAMPLE_BODY_TRIMMED = "1\t100\t.\tA\tG\t.\tPASS\tAC=1;AN=2;AF=0.5\n"
SAMPLE_VCF = SAMPLE_HEADER + SAMPLE_BODY_WITH_FORMAT


def _configure_fake_bcftools(monkeypatch, module, header_checks, call_log=None):
    def fake_run(cmd, check=True, stdout=None, stdin=None, **kwargs):
        if call_log is not None:
            call_log.append(list(cmd))
        if cmd[:2] == ["bcftools", "+fill-tags"]:
            output_path = Path(cmd[cmd.index("-o") + 1])
            output_path.write_text(SAMPLE_VCF)
            return None

        if cmd[:2] == ["bcftools", "view"] and "-i" in cmd:
            output_path = Path(cmd[cmd.index("-o") + 1])
            output_path.write_text(SAMPLE_VCF)
            return None

        if cmd[:3] == ["bcftools", "view", "-h"] and stdout not in (
            None,
            module.subprocess.PIPE,
        ):
            stdout.write(SAMPLE_HEADER)
            return None

        if cmd[:3] == ["bcftools", "view", "-h"] and stdout == module.subprocess.PIPE:
            header_checks.append(Path(cmd[-1]))
            return None

        if cmd and cmd[0] == "cut":
            if stdout is None:
                raise AssertionError("cut command requires stdout to be provided")
            stdout.write(SAMPLE_BODY_TRIMMED)
            return None

        if cmd and cmd[0] == "bgzip":
            if len(cmd) < 2:
                raise AssertionError("bgzip command missing target path")
            target = Path(cmd[-1])
            compressed = target.with_name(target.name + ".gz")
            data = target.read_bytes()
            with gzip.open(compressed, "wb") as handle:
                handle.write(data)
            return None

        if cmd and cmd[0] == "tabix":
            if "-f" not in cmd:
                raise AssertionError("tabix command missing -f flag")
            target = Path(cmd[-1])
            index_path = target.with_suffix(target.suffix + ".tbi")
            index_path.write_text("index")
            return None

        raise AssertionError(f"Unexpected command: {cmd}")

    class FakePopen:
        def __init__(self, cmd, stdout=None, **kwargs):
            if cmd[:3] != ["bcftools", "view", "-H"]:
                raise AssertionError(f"Unexpected Popen command: {cmd}")
            if stdout != module.subprocess.PIPE:
                raise AssertionError("Popen must be called with stdout=PIPE")
            self.stdout = io.BytesIO(SAMPLE_BODY_WITH_FORMAT.encode("utf-8"))

        def wait(self):
            return 0

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.setattr(module.subprocess, "Popen", FakePopen)


def test_append_metadata_runs_header_validation_for_gzip(
    monkeypatch, tmp_path, merge_script_module
):
    module = merge_script_module
    header_checks = []
    _configure_fake_bcftools(monkeypatch, module, header_checks)

    merged_vcf = tmp_path / "merged.vcf.gz"
    merged_vcf.write_text("placeholder")

    final_vcf = module.append_metadata_to_merged_vcf(str(merged_vcf))

    assert final_vcf.endswith(".vcf.gz")
    assert Path(final_vcf) == header_checks[0]

    with gzip.open(final_vcf, "rt", encoding="utf-8") as handle:
        contents = handle.read()
    assert "##fileformat=VCFv4.2" in contents
    assert "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" in contents


def test_append_metadata_runs_header_validation_for_plain_vcf(
    monkeypatch, tmp_path, merge_script_module
):
    module = merge_script_module
    header_checks = []
    _configure_fake_bcftools(monkeypatch, module, header_checks)

    merged_vcf = tmp_path / "merged.vcf"
    merged_vcf.write_text("placeholder")

    final_vcf = module.append_metadata_to_merged_vcf(str(merged_vcf))

    assert final_vcf.endswith(".vcf")
    assert Path(final_vcf) == header_checks[0]
    assert "##fileformat=VCFv4.2" in Path(final_vcf).read_text()


def test_append_metadata_triggers_bgzip_and_tabix(
    monkeypatch, tmp_path, merge_script_module
):
    module = merge_script_module
    header_checks = []
    call_log = []
    _configure_fake_bcftools(monkeypatch, module, header_checks, call_log=call_log)

    merged_vcf = tmp_path / "merged.vcf.gz"
    merged_vcf.write_text("placeholder")

    final_vcf = module.append_metadata_to_merged_vcf(str(merged_vcf))

    assert any(cmd[0] == "bgzip" for cmd in call_log)
    assert any(cmd[0] == "tabix" for cmd in call_log)

    final_path = Path(final_vcf)
    assert final_path.exists()
    assert final_path.with_suffix(final_path.suffix + ".tbi").exists()

    plain_path = Path(str(final_vcf)[:-3])
    assert not plain_path.exists()


def test_validate_header_failure_raises(monkeypatch, merge_script_module):
    module = merge_script_module

    def fake_run(cmd, check=True, **kwargs):
        raise module.subprocess.CalledProcessError(1, cmd)

    captured = {}

    def fake_handle_error(message, **_kwargs):
        captured["message"] = message
        raise RuntimeError("critical failure")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    monkeypatch.setattr(module, "handle_critical_error", fake_handle_error)

    with pytest.raises(RuntimeError):
        module._validate_anonymized_vcf_header("output.vcf.gz")

    assert "bcftools" in captured["message"]
