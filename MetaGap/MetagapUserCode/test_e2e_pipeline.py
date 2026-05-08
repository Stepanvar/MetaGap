#!/usr/bin/env python3
"""End-to-end pipeline test: generate multiVCF → merge → import → verify."""
from __future__ import annotations

import gzip
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

# ── paths ────────────────────────────────────────────────────────────────────
REPO_ROOT = Path(__file__).resolve().parent.parent          # MetaGap/MetaGap/
SCRIPT_DIR = Path(__file__).resolve().parent                # MetagapUserCode/
WES_VCF = SCRIPT_DIR / "ZD8U010E-GED-E080TSO-1715980255-BWA-DPV-hs37d5.filter.WES.vcf.gz"
METADATA_TXT = SCRIPT_DIR / "metadata.txt"
WORK_DIR = Path(tempfile.mkdtemp(prefix="metagap_e2e_"))
MULTIVCF_PLAIN = WORK_DIR / "multivcf.vcf"
MULTIVCF_BGZ = WORK_DIR / "multivcf.vcf.gz"
MERGED_DIR = WORK_DIR / "merged"

N_SAMPLES = 15
SAMPLE_NAMES = [f"SAMPLE_{i:03d}" for i in range(1, N_SAMPLES + 1)]

FAILURES: list[str] = []
PASSES: list[str] = []


def fail(msg: str) -> None:
    FAILURES.append(msg)
    print(f"  FAIL: {msg}")


def ok(msg: str) -> None:
    PASSES.append(msg)
    print(f"  OK:   {msg}")


# ── Phase 1: Generate MultiVCF ───────────────────────────────────────────────

def generate_multivcf() -> bool:
    """Read WES VCF, clone sample column ×15, write plain-text + BGZF multiVCF.

    Adds ##reference=GRCh38 if absent (required by merge_vcf validation).
    Returns True on success.
    """
    print("\n=== Phase 1: Generate MultiVCF ===")
    if not WES_VCF.exists():
        fail(f"WES VCF not found: {WES_VCF}")
        return False

    try:
        with gzip.open(WES_VCF, "rt", encoding="utf-8") as fh:
            lines = fh.readlines()
    except Exception as exc:
        fail(f"Cannot read WES VCF: {exc}")
        return False

    header: list[str] = [l for l in lines if l.startswith("#")]
    data: list[str] = [l for l in lines if not l.startswith("#")]

    # Ensure ##reference= present (merge_vcf validation requires it)
    has_reference = any(l.lower().startswith("##reference=") for l in header)
    if not has_reference:
        insert_at = next(
            (i for i, l in enumerate(header) if l.lower().startswith("##fileformat=")), 0
        )
        header.insert(insert_at + 1, "##reference=GRCh38\n")
        ok("Injected missing ##reference=GRCh38 into header")

    # Rewrite #CHROM line with 15 sample columns
    chrom_idx = next(i for i, l in enumerate(header) if l.startswith("#CHROM"))
    chrom_fields = header[chrom_idx].rstrip("\n").split("\t")
    fixed_cols = chrom_fields[:9]  # CHROM..FORMAT
    new_chrom = "\t".join(fixed_cols + SAMPLE_NAMES) + "\n"
    header[chrom_idx] = new_chrom

    # Clone sample column N_SAMPLES times per data row
    new_data: list[str] = []
    skipped = 0
    for line in data:
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 10:
            skipped += 1
            continue
        sample_col = fields[9]
        new_data.append("\t".join(fields[:9] + [sample_col] * N_SAMPLES) + "\n")

    if skipped:
        fail(f"{skipped} data lines had < 10 fields and were skipped")

    MULTIVCF_PLAIN.write_text("".join(header + new_data), encoding="utf-8")
    ok(f"Plain multiVCF written: {len(new_data)} variants × {N_SAMPLES} samples")

    # BGZF-compress + tabix-index
    try:
        import pysam
        pysam.tabix_compress(str(MULTIVCF_PLAIN), str(MULTIVCF_BGZ), force=True)
        pysam.tabix_index(str(MULTIVCF_BGZ), preset="vcf", force=True)
        ok(f"BGZF compressed and indexed: {MULTIVCF_BGZ}")
    except Exception as exc:
        fail(f"BGZF compress/index failed: {exc}")
        return False

    return True


# ── Phase 2: Run merge_vcf ───────────────────────────────────────────────────

def run_merge_vcf() -> Path | None:
    """Run merge_vcf on the multiVCF with metadata template.

    Returns path to merged .vcf.gz on success, None on failure.
    """
    print("\n=== Phase 2: Run merge_vcf ===")
    if not MULTIVCF_BGZ.exists():
        fail("MultiVCF bgz not found — Phase 1 must succeed first")
        return None

    MERGED_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, "-m", "MetagapUserCode.merge_vcf",
        str(MULTIVCF_BGZ),
        "-m", str(METADATA_TXT),
        "-o", str(MERGED_DIR),
        "--an-threshold", "-1",   # disable: default 50 silently drops all with <25 samples
        "--qual-threshold", "-1", # disable: keep all variants for import testing
        "-v",
    ]

    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=str(REPO_ROOT),
            timeout=300,
        )
    except subprocess.TimeoutExpired:
        fail("merge_vcf timed out after 300s")
        return None
    except Exception as exc:
        fail(f"merge_vcf subprocess error: {exc}")
        return None

    print(f"  [merge_vcf exit code: {result.returncode}]")
    if result.stdout.strip():
        for line in result.stdout.strip().splitlines()[-20:]:
            print(f"  stdout: {line}")
    if result.stderr.strip():
        for line in result.stderr.strip().splitlines()[-30:]:
            print(f"  stderr: {line}")

    if result.returncode != 0:
        fail(f"merge_vcf exited {result.returncode}")
    else:
        ok("merge_vcf exited 0")

    # Locate merged output
    candidates = sorted(MERGED_DIR.glob("*.vcf.gz"))
    if not candidates:
        fail(f"No .vcf.gz found in {MERGED_DIR}")
        return None

    merged = candidates[0]
    ok(f"Merged output: {merged} ({merged.stat().st_size // 1024} KB)")

    # Spot-check merged VCF (cohort_final is anonymized — 0 samples by design)
    try:
        import pysam
        with pysam.VariantFile(str(merged)) as vcf:
            n_records = sum(1 for _ in vcf.fetch())
            n_samples = len(vcf.header.samples)
        ok(f"Merged VCF: {n_records} records, {n_samples} samples (anonymized → samples stripped)")
        if n_records == 0:
            fail("Merged VCF has 0 records — all variants were filtered out")
    except Exception as exc:
        fail(f"Cannot open merged VCF with pysam: {exc}")

    return merged


# ── Phase 3: Import to Django ────────────────────────────────────────────────

def setup_django() -> bool:
    """Configure Django for in-process use."""
    print("\n=== Phase 3: Django Setup ===")
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "MetaGap.settings")
    os.environ["DEBUG"] = "1"
    sys.path.insert(0, str(REPO_ROOT))
    try:
        import django
        django.setup()
        ok("Django configured")
        return True
    except Exception as exc:
        fail(f"django.setup() failed: {exc}")
        return False


def import_vcf(merged_vcf: Path) -> Any | None:
    """Create test user + org, POST merged VCF to /import/, return SampleGroup or None."""
    print("\n=== Phase 3: Import VCF ===")
    try:
        from django.contrib.auth.models import User
        from django.core.files.uploadedfile import SimpleUploadedFile
        from django.test import Client
        from django.test.utils import setup_test_environment
        setup_test_environment()
    except Exception as exc:
        fail(f"Django import error: {exc}")
        return None

    # Create test user (delete previous run's user if exists)
    username = "e2e_pipeline_user"
    try:
        from django.contrib.auth.models import User
        User.objects.filter(username=username).delete()
        user = User.objects.create_user(username=username, password="e2e_pass_123")
        ok(f"Test user created: {username}")
    except Exception as exc:
        fail(f"Could not create test user: {exc}")
        return None

    # OrganizationProfile auto-created by signal
    try:
        org = user.organization_profile
        ok(f"OrganizationProfile found: pk={org.pk}")
    except Exception as exc:
        fail(f"OrganizationProfile missing after user create: {exc}")
        return None

    # Read merged VCF bytes
    try:
        vcf_bytes = merged_vcf.read_bytes()
        ok(f"Read merged VCF: {len(vcf_bytes) // 1024} KB")
    except Exception as exc:
        fail(f"Cannot read merged VCF file: {exc}")
        return None

    # POST to /en/import/ (all app URLs wrapped in i18n_patterns → locale prefix required)
    client = Client()
    client.force_login(user)
    uploaded = SimpleUploadedFile(
        "merged.vcf.gz", vcf_bytes, content_type="application/gzip"
    )
    try:
        response = client.post("/en/import/", {"data_file": uploaded}, follow=True)
    except Exception as exc:
        fail(f"POST /import/ raised exception: {exc}")
        return None

    print(f"  [POST /en/import/ → HTTP {response.status_code}]")
    if response.status_code != 200:
        fail(f"POST /en/import/ returned HTTP {response.status_code}")
    else:
        ok("POST /en/import/ → HTTP 200")

    # Flash messages — check context (consumed by template on redirect) + session storage
    try:
        ctx_msgs = list(response.context.get("messages", []) if response.context else [])
    except Exception:
        ctx_msgs = []
    from django.contrib.messages import get_messages
    session_msgs = list(get_messages(response.wsgi_request))
    all_msgs = ctx_msgs or session_msgs
    if not all_msgs:
        # Last resort: peek at response HTML for known error/success markers
        html = response.content.decode("utf-8", errors="replace")
        if "alert-success" in html or "successfully" in html.lower():
            ok("Import success detected in HTML")
        elif "alert-danger" in html or "error" in html.lower():
            snippet = html[max(0, html.lower().find("error")-50):html.lower().find("error")+200]
            fail(f"Import error detected in HTML: ...{snippet}...")
        else:
            print("  [flash] no flash messages captured")
    for msg in all_msgs:
        tag = str(msg.tags)
        print(f"  [flash {tag}] {msg}")
        if "error" in tag:
            fail(f"Import flash error: {msg}")
        elif "success" in tag:
            ok(f"Import flash success: {msg}")

    # Retrieve created SampleGroup
    try:
        from app.models import SampleGroup
        sg = SampleGroup.objects.filter(created_by=org).order_by("-pk").first()
        if sg is None:
            fail("No SampleGroup found in DB after import")
            return None
        ok(f"SampleGroup created: pk={sg.pk}, name={sg.name!r}")
        return sg
    except Exception as exc:
        fail(f"SampleGroup query failed: {exc}")
        return None


# ── Phase 4a: Verify Storage ─────────────────────────────────────────────────

def verify_storage(sg: Any, expected_variants: int) -> None:
    """Check DB state against expectations."""
    print("\n=== Phase 4a: Verify Storage ===")

    af_count = sg.allele_frequencies.count()
    print(f"  AlleleFrequency rows: {af_count} (expected ≈{expected_variants})")
    if af_count == 0:
        fail("AlleleFrequency count is 0 — no variants imported")
    elif af_count != expected_variants:
        fail(
            f"AlleleFrequency count mismatch: {af_count} stored vs {expected_variants} in source VCF"
        )
    else:
        ok(f"AlleleFrequency count matches: {af_count}")

    # Metadata FK fields
    checks: list[tuple[str, Any]] = [
        ("reference_genome_build", sg.reference_genome_build),
        ("sequencing_platform", sg.sequencing_platform),
        ("bioinfo_alignment", sg.bioinfo_alignment),
        ("bioinfo_variant_calling", sg.bioinfo_variant_calling),
        ("sample_origin", sg.sample_origin),
        ("library_construction", sg.library_construction),
    ]
    for field_name, value in checks:
        if value:
            ok(f"FK populated: {field_name} = {value}")
        else:
            fail(f"FK not populated: {field_name} is None")

    # chrom field lengths — regression for C4 fix (max_length=64)
    from app.models import AlleleFrequency
    chroms = list(
        sg.allele_frequencies.values_list("chrom", flat=True).distinct()
    )
    for chrom in chroms:
        if len(chrom) > 10:
            ok(f"Long chrom stored correctly (len={len(chrom)}): {chrom}")
    ok(f"Distinct chroms in DB: {sorted(chroms)[:10]}")

    # additional_metadata
    if sg.additional_metadata:
        ok(f"additional_metadata present: {list(sg.additional_metadata.keys())[:5]}")
    else:
        fail("additional_metadata is empty — metadata header may not have parsed")


# ── Phase 4b: Verify HTTP Views ──────────────────────────────────────────────

def verify_views(sg: Any) -> None:
    """GET detail + search + export views, check HTTP 200 and content."""
    print("\n=== Phase 4b: Verify HTTP Views ===")
    from django.contrib.auth.models import User
    from django.test import Client

    user = User.objects.filter(username="e2e_pipeline_user").first()
    if user is None:
        fail("Test user missing — cannot verify views")
        return

    client = Client()
    client.force_login(user)

    # Detail view
    detail_url = f"/en/profile/sample-groups/{sg.pk}/"
    try:
        resp = client.get(detail_url)
        if resp.status_code == 200:
            ok(f"GET {detail_url} → 200")
            content = resp.content.decode("utf-8", errors="replace")
            if sg.name in content:
                ok(f"Sample group name '{sg.name}' found in detail page")
            else:
                fail(f"Sample group name '{sg.name}' NOT found in detail page HTML")
            if "Sequencing" in content or "sequencing" in content:
                ok("Sequencing metadata section rendered")
            else:
                fail("Sequencing metadata section missing from detail page")
        else:
            fail(f"GET {detail_url} → HTTP {resp.status_code}")
    except Exception as exc:
        fail(f"GET {detail_url} raised: {exc}")

    # Search view — query for first variant's position
    try:
        from app.models import AlleleFrequency
        first_af = sg.allele_frequencies.order_by("pk").first()
        if first_af:
            search_url = (
                f"/en/search/?chrom={first_af.chrom}"
                f"&pos_min={first_af.pos}&pos_max={first_af.pos}"
            )
            resp = client.get(search_url)
            if resp.status_code == 200:
                ok(f"GET {search_url} → 200")
                content = resp.content.decode("utf-8", errors="replace")
                if str(first_af.pos) in content:
                    ok(f"Variant pos {first_af.pos} found in search results")
                else:
                    fail(f"Variant pos {first_af.pos} NOT found in search HTML")
            else:
                fail(f"GET {search_url} → HTTP {resp.status_code}")
        else:
            fail("No AlleleFrequency rows to search for")
    except Exception as exc:
        fail(f"Search view raised: {exc}")

    # Export endpoint
    export_url = f"/en/profile/sample-groups/{sg.pk}/export/?format=csv"
    try:
        resp = client.get(export_url)
        if resp.status_code == 200:
            ok(f"GET {export_url} → 200 ({len(resp.content)} bytes)")
        else:
            fail(f"GET {export_url} → HTTP {resp.status_code}")
    except Exception as exc:
        fail(f"GET {export_url} raised: {exc}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main() -> None:
    print(f"Work dir: {WORK_DIR}")
    print(f"WES source: {WES_VCF}")

    # Phase 1
    phase1_ok = generate_multivcf()

    # Phase 2
    merged_vcf: Path | None = None
    if phase1_ok:
        merged_vcf = run_merge_vcf()

    # Count source variants for comparison
    expected_variants = 0
    try:
        with gzip.open(WES_VCF, "rt") as fh:
            expected_variants = sum(1 for l in fh if not l.startswith("#") and l.strip())
    except Exception:
        pass

    # Phase 3
    sg = None
    if merged_vcf:
        if setup_django():
            sg = import_vcf(merged_vcf)

    # Phase 4
    if sg is not None:
        verify_storage(sg, expected_variants)
        verify_views(sg)

    # ── Final report ─────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print(f"RESULTS: {len(PASSES)} passed, {len(FAILURES)} failed")
    print("=" * 60)
    if FAILURES:
        print("\nFAILURES:")
        for i, f in enumerate(FAILURES, 1):
            print(f"  {i}. {f}")
    else:
        print("\nAll checks passed!")
    print()


if __name__ == "__main__":
    main()
