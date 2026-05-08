# E2E Pipeline Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the full MetaGap pipeline end-to-end — generate a synthetic 15-sample VCF from the real WES file, merge it with metadata via `merge_vcf`, import to the Django app, and report every failure found.

**Architecture:** Single self-contained script (`test_e2e_pipeline.py`) runs four sequential phases; each phase catches its own errors so later phases still run. A `FAILURES` list accumulates every problem; a final report is printed. Uses Django's test client (no real HTTP server) for the import and verification steps.

**Tech Stack:** Python 3.11, pysam (bgzf + tabix), Django 5.1 test client, `MetagapUserCode.merge_vcf` CLI (subprocess), pytest-less (plain script).

---

## File Map

| File | Action | Purpose |
|------|--------|---------|
| `MetaGap/MetagapUserCode/test_e2e_pipeline.py` | **Create** | Main test script — all four phases |
| `MetaGap/MetagapUserCode/metadata.txt` | Read-only | Metadata template for merge_vcf `-m` flag |
| `MetaGap/MetagapUserCode/ZD8U010E-GED-E080TSO-1715980255-BWA-DPV-hs37d5.filter.WES.vcf.gz` | Read-only | Real WES VCF source |

---

## Task 1: Generate 15-Sample MultiVCF from Real WES File

**Files:**
- Create: `MetaGap/MetagapUserCode/test_e2e_pipeline.py` (phase 1 section)

- [ ] **Step 1.1: Create the script with imports and failure tracker**

Create `MetaGap/MetagapUserCode/test_e2e_pipeline.py`:

```python
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
```

- [ ] **Step 1.2: Add `generate_multivcf()` function**

Append to `test_e2e_pipeline.py`:

```python
def generate_multivcf() -> bool:
    """Read WES VCF, clone sample column ×15, write plain-text multiVCF.
    
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
        # Insert after ##fileformat line
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
```

- [ ] **Step 1.3: Smoke-test phase 1 in isolation**

```bash
cd /home/szuev/MetaGap/MetaGap
DEBUG=1 ~/miniconda3/envs/metagap-py311/bin/python - <<'EOF'
import sys; sys.path.insert(0, ".")
exec(open("MetagapUserCode/test_e2e_pipeline.py").read())
generate_multivcf()
EOF
```

Expected output:
```
=== Phase 1: Generate MultiVCF ===
  OK:   Injected missing ##reference=GRCh38 into header
  OK:   Plain multiVCF written: 8326 variants × 15 samples
  OK:   BGZF compressed and indexed: /tmp/metagap_e2e_XXX/multivcf.vcf.gz
```

If `pysam` import fails, `BGZF compress/index failed` will appear — expected failure to capture.

---

## Task 2: Run merge_vcf to Produce exVCF

**Files:**
- Modify: `MetaGap/MetagapUserCode/test_e2e_pipeline.py` (append phase 2)

- [ ] **Step 2.1: Add `run_merge_vcf()` function**

Append to `test_e2e_pipeline.py`:

```python
def run_merge_vcf() -> Path | None:
    """Run merge_vcf on the multiVCF with metadata template.
    
    Returns path to merged .vcf.gz on success, None on failure.
    Captures stdout/stderr and logs any warnings/errors as failures.
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
        # Still try to find output — merge may produce partial output
    else:
        ok("merge_vcf exited 0")

    # Locate merged output (merge_vcf writes merged.vcf.gz or merged_<timestamp>.vcf.gz)
    candidates = sorted(MERGED_DIR.glob("*.vcf.gz"))
    if not candidates:
        fail(f"No .vcf.gz found in {MERGED_DIR}")
        return None

    merged = candidates[0]
    ok(f"Merged output: {merged} ({merged.stat().st_size // 1024} KB)")

    # Spot-check: count variant records in merged output
    try:
        import pysam
        with pysam.VariantFile(str(merged)) as vcf:
            n_records = sum(1 for _ in vcf.fetch())
            n_samples = len(vcf.header.samples)
        ok(f"Merged VCF: {n_records} records, {n_samples} samples")
        if n_samples != N_SAMPLES:
            fail(f"Expected {N_SAMPLES} samples in merged VCF, got {n_samples}")
    except Exception as exc:
        fail(f"Cannot open merged VCF with pysam: {exc}")

    return merged
```

- [ ] **Step 2.2: Smoke-test phase 2**

```bash
cd /home/szuev/MetaGap/MetaGap
DEBUG=1 ~/miniconda3/envs/metagap-py311/bin/python - <<'EOF'
import sys; sys.path.insert(0, ".")
exec(open("MetagapUserCode/test_e2e_pipeline.py").read())
generate_multivcf()
merged = run_merge_vcf()
print("merged path:", merged)
EOF
```

Expected: either `OK: merge_vcf exited 0` + merged path, or one or more `FAIL:` lines describing what went wrong (missing `##reference` enforcement, AC/AN/AF issues, etc.). All failures are captured — the script continues regardless.

---

## Task 3: Import exVCF to Django via Test Client

**Files:**
- Modify: `MetaGap/MetagapUserCode/test_e2e_pipeline.py` (append phase 3)

- [ ] **Step 3.1: Add Django setup block**

Append to `test_e2e_pipeline.py` (before import function, so Django is configured once):

```python
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
```

- [ ] **Step 3.2: Add `import_vcf()` function**

Append to `test_e2e_pipeline.py`:

```python
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

    # Create (or reuse) test user
    username = "e2e_pipeline_user"
    try:
        from django.contrib.auth.models import User
        User.objects.filter(username=username).delete()
        user = User.objects.create_user(username=username, password="e2e_pass_123")
        ok(f"Test user created: {username}")
    except Exception as exc:
        fail(f"Could not create test user: {exc}")
        return None

    # OrganizationProfile should be auto-created by signal
    try:
        org = user.organization_profile
        ok(f"OrganizationProfile found: pk={org.pk}")
    except Exception as exc:
        fail(f"OrganizationProfile missing after user create: {exc}")
        return None

    # Read merged VCF
    try:
        vcf_bytes = merged_vcf.read_bytes()
        ok(f"Read merged VCF: {len(vcf_bytes) // 1024} KB")
    except Exception as exc:
        fail(f"Cannot read merged VCF file: {exc}")
        return None

    # POST to /import/
    client = Client()
    client.force_login(user)
    uploaded = SimpleUploadedFile(
        "merged.vcf.gz", vcf_bytes, content_type="application/gzip"
    )
    try:
        response = client.post("/import/", {"data_file": uploaded}, follow=True)
    except Exception as exc:
        fail(f"POST /import/ raised exception: {exc}")
        return None

    print(f"  [POST /import/ → HTTP {response.status_code}]")
    if response.status_code != 200:
        fail(f"POST /import/ returned HTTP {response.status_code}")
    else:
        ok("POST /import/ → HTTP 200")

    # Extract flash messages from response
    from django.contrib.messages import get_messages
    msgs = list(get_messages(response.wsgi_request))
    for msg in msgs:
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
```

- [ ] **Step 3.3: Smoke-test import phase (requires phase 1+2 to succeed)**

```bash
cd /home/szuev/MetaGap/MetaGap
DEBUG=1 ~/miniconda3/envs/metagap-py311/bin/python - <<'EOF'
import sys; sys.path.insert(0, ".")
exec(open("MetagapUserCode/test_e2e_pipeline.py").read())
generate_multivcf()
merged = run_merge_vcf()
if merged:
    setup_django()
    sg = import_vcf(merged)
    print("SampleGroup:", sg)
EOF
```

Expected: `OK: Import flash success: Imported ...` or one or more `FAIL:` entries if import pipeline chokes (pysam parse error, metadata mapping failure, DB constraint, etc.).

---

## Task 4: Verify Storage and HTTP Views, Print Report

**Files:**
- Modify: `MetaGap/MetagapUserCode/test_e2e_pipeline.py` (append phase 4 + `main()`)

- [ ] **Step 4.1: Add `verify_storage()` function**

Append to `test_e2e_pipeline.py`:

```python
def verify_storage(sg: Any, expected_variants: int) -> None:
    """Check DB state against expectations."""
    print("\n=== Phase 4a: Verify Storage ===")
    from app.models import AlleleFrequency, SampleGroup

    # Variant count
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

    # chrom field lengths (regression for C4 fix — max_length=64)
    long_chroms = (
        sg.allele_frequencies
        .values_list("chrom", flat=True)
        .distinct()
    )
    for chrom in long_chroms:
        if len(chrom) > 10:
            ok(f"Long chrom stored correctly (len={len(chrom)}): {chrom}")

    # additional_metadata catch-all
    if sg.additional_metadata:
        ok(f"additional_metadata present: {list(sg.additional_metadata.keys())[:5]}")
    else:
        fail("additional_metadata is empty — metadata header may not have parsed")
```

- [ ] **Step 4.2: Add `verify_views()` function**

Append to `test_e2e_pipeline.py`:

```python
def verify_views(sg: Any) -> None:
    """GET detail + search views, check HTTP 200 and content."""
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
    detail_url = f"/profile/sample-groups/{sg.pk}/"
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
                ok("Sequencing metadata section rendered in detail page")
            else:
                fail("Sequencing metadata section missing from detail page")
        else:
            fail(f"GET {detail_url} → HTTP {resp.status_code}")
    except Exception as exc:
        fail(f"GET {detail_url} raised: {exc}")

    # Search view — use first variant's chrom as query
    try:
        from app.models import AlleleFrequency
        first_af = sg.allele_frequencies.order_by("pk").first()
        if first_af:
            search_url = f"/search/?chrom={first_af.chrom}&pos_min={first_af.pos}&pos_max={first_af.pos}"
            resp = client.get(search_url)
            if resp.status_code == 200:
                ok(f"GET {search_url} → 200")
                content = resp.content.decode("utf-8", errors="replace")
                if str(first_af.pos) in content:
                    ok(f"Variant pos {first_af.pos} found in search results")
                else:
                    fail(f"Variant pos {first_af.pos} NOT found in search results HTML")
            else:
                fail(f"GET {search_url} → HTTP {resp.status_code}")
        else:
            fail("No AlleleFrequency rows to search for")
    except Exception as exc:
        fail(f"Search view raised: {exc}")

    # Export endpoint
    export_url = f"/profile/sample-groups/{sg.pk}/export/?format=csv"
    try:
        resp = client.get(export_url)
        if resp.status_code == 200:
            ok(f"GET {export_url} → 200 ({len(resp.content)} bytes)")
        else:
            fail(f"GET {export_url} → HTTP {resp.status_code}")
    except Exception as exc:
        fail(f"GET {export_url} raised: {exc}")
```

- [ ] **Step 4.3: Add `main()` and final report**

Append to `test_e2e_pipeline.py`:

```python
def main() -> None:
    print(f"Work dir: {WORK_DIR}")
    print(f"WES source: {WES_VCF}")

    # Phase 1
    phase1_ok = generate_multivcf()

    # Phase 2
    merged_vcf: Path | None = None
    if phase1_ok:
        merged_vcf = run_merge_vcf()

    # Count source variants for later comparison
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
```

- [ ] **Step 4.4: Run the full pipeline**

```bash
cd /home/szuev/MetaGap/MetaGap
DEBUG=1 ~/miniconda3/envs/metagap-py311/bin/python MetagapUserCode/test_e2e_pipeline.py 2>&1 | tee /tmp/e2e_report.txt
```

Expected: full output ending with a RESULTS section. All `FAIL:` lines are the issues to fix in the next session. Do not fix anything now — just collect the report.

- [ ] **Step 4.5: Commit the test script and plan**

```bash
cd /home/szuev/MetaGap
git add MetaGap/MetagapUserCode/test_e2e_pipeline.py docs/superpowers/plans/2026-05-08-e2e-pipeline-test.md
git commit -m "test: add e2e pipeline test script for multiVCF → merge → import → verify"
```

---

## Self-Review Checklist

- [x] **Spec coverage:** Phase 1 (generate multiVCF) ✓ Phase 2 (merge_vcf) ✓ Phase 3 (import) ✓ Phase 4 (verify storage + views) ✓ Failure report ✓
- [x] **Placeholders:** None — all code blocks are complete
- [x] **Type consistency:** `FAILURES`/`PASSES` lists used consistently; `fail()`/`ok()` helpers used in all phases; `sg` passed from phase 3 → phase 4
- [x] **Edge cases covered:** phase 2 skips if phase 1 failed; phase 3 skips if phase 2 failed; phase 4 skips if phase 3 failed — no uncaught crash cascade
- [x] **Known expected failures documented:** `##reference=GRCh38` injected in phase 1 (validation requirement); merge_vcf may still fail on missing AC/AN/AF headers in WES input — that failure will surface in report
