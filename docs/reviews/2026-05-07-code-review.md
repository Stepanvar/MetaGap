# MetaGap Code Review — 2026-05-07

**Scope:** Full app — `app/` (Django 5.1 web) + `MetagapUserCode/merge_vcf/` (CLI tool)
**Test baseline:** 154 tests pass, 64s, 76% total coverage
**Branch:** master (14 commits ahead of origin)

---

## Recon Summary

Two independent subsystems share the repo:

- **Django web app** (`MetaGap/app/`, `MetaGap/MetaGap/`): models + VCF import pipeline + views/forms/filters. Entry: `VCFImporter.import_file` → `VCFMetadataParser` → `VCFDatabaseWriter`. Config-driven via `app/config/metadata_fields.yaml`.
- **CLI VCF merger** (`MetagapUserCode/merge_vcf/`): standalone pysam+vcfpy tool. No Django imports; fully independent.

No cross-imports between subsystems. `app/` uses only Django models and its own services.

Hotspots by line count: `vcf_database.py` (1379), `views.py` (860), `forms.py` (773), `vcf_metadata.py` (676), `models.py` (671).

---

## Issues by Severity

### CORRECTNESS

#### C1 — Dead code in `VCFImporter.import_file`: unreachable `sample_group.delete()` branch
**File:** `app/services/vcf_importer.py:110-112`
```python
cleanup_pysam_sample_group()        # sets pysam_sample_group=None, deletes if set
if sample_group is not None:        # sample_group is ALWAYS None here (never assigned from pysam at this branch)
    sample_group.delete()           # UNREACHABLE
```
`sample_group` is only assigned inside `_import_with_pysam`, which raised an exception to get here. The second `if` check dead-letter wastes reader attention and could mask future refactors that do set `sample_group` earlier.

#### C2 — Broken `get_or_create_sample_group()` in vcf_database.py — always returns None
**File:** `app/services/vcf_database.py:188-217`
```python
def get_or_create_sample_group(...):
    ...
    sample_group = SampleGroup.objects.create(
        name=sample_group_name,
        description=sample_group_description,   # field doesn't exist on SampleGroup
        organization_profile=user_organization_profile,  # field is `created_by`
    )
```
`SampleGroup` has no `description` or `organization_profile` field. `objects.create()` raises `TypeError`, caught by the bare `except Exception: return None`. Function is exported in `__all__` and labeled "for testing" but is silently broken. Any caller gets `None` with no diagnostic.

#### C3 — Double queryset evaluation in `SampleGroupTableView.get_context_data`
**File:** `app/views.py:113-121`
```python
def get_context_data(self, **kwargs):
    context = super().get_context_data(**kwargs)   # calls get_queryset() internally via ListView
    ...
    table = table_class(self.get_queryset())        # get_queryset() called AGAIN
    RequestConfig(...).configure(table)
    context["table"] = table
    return context
```
`ListView.get_context_data` already calls `get_queryset()` (which has 13-relation `select_related` + `prefetch_related`). Calling it a second time to build the table doubles DB I/O. The table's data is unused by the template in favor of the `table` object anyway — the two querysets drift.

**Fix:** use `context["object_list"]` (already populated) instead of `self.get_queryset()`.

#### C4 — `AlleleFrequency.chrom` max_length=10 truncates modern assembly chr names
**File:** `app/models.py:647`
```python
chrom = models.CharField(max_length=10)
```
hg38 alt-contig names reach 26+ chars (e.g. `chr1_GL383518v1_alt`). Importing such a VCF hits `DataError` from the DB engine. The pipeline catches `DataError` and wraps as `ValidationError` with a generic message (`"One or more metadata values are out of range"`), losing the actual field context.

#### C5 — Syntax error in MetagapUserCode test file blocks pytest collection entirely
**File:** `MetagapUserCode/tests/test_merging.py:331`
```python
merged = merging._merge_colliding_records(
    [(rec1, 0), (rec2, 1)],
def test_merge_colliding_records_pads_missing_samples():   # parenthesis never closed
```
Unclosed `(` on line 331. `pytest` cannot collect `test_merging.py`. Because it interrupts collection with `--interrupt-on-collection-errors` behavior, it can block the full pytest run depending on config.

---

### SECURITY

#### S1 — Weak SECRET_KEY and DEBUG defaults in settings.py
**File:** `MetaGap/MetaGap/settings.py:30-31`
```python
DEBUG = os.getenv("DEBUG", "1") == "1"
SECRET_KEY = os.getenv("SECRET_KEY", "dev-secret")
```
Default `DEBUG=True` exposes full tracebacks in prod if env var is unset. `SECRET_KEY="dev-secret"` is public knowledge — Django uses it to sign session cookies and CSRF tokens. Any prod deployment without explicit env vars is trivially exploitable (session forgery, CSRF bypass).

**Fix:** no default for `SECRET_KEY` (raise `ImproperlyConfigured`); default `DEBUG=False`.

#### S2 — Broad `ALLOWED_HOSTS` wildcard includes GitHub Codespaces by default
**File:** `MetaGap/MetaGap/settings.py:32-34`
```python
ALLOWED_HOSTS = _split_env_list("ALLOWED_HOSTS", "localhost,127.0.0.1,0.0.0.0,*.github.dev")
```
`0.0.0.0` and `*.github.dev` in the default list are dev conveniences, but if this config reaches a prod server without overriding `ALLOWED_HOSTS`, any request with a matching Host header is accepted.

**Fix:** remove `0.0.0.0` and `*.github.dev` from default; keep only `localhost,127.0.0.1`.

#### S3 — Temp VCF upload not cleaned up on `ValidationError` path in `ImportDataView`
**File:** `app/views.py:757-806`
The `finally:` block in `form_valid` only runs after the `with transaction.atomic()` block. If `ValidationError` is raised inside the importer, it's caught and re-raised as `ImporterValidationError`. The `default_storage.delete(temp_path)` in `finally:` does run (finally runs regardless of exception). This is actually correct — but the `except Exception` clause at line 778 logs `exc_info=exc` which passes the exception object directly to `logger.exception` (which already captures `exc_info` from the current exception). Minor logging redundancy but not a bug.

No critical issue here — noted for awareness only.

---

### COVERAGE

#### COV1 — Three test files use `pytest.mark.django_db` but `pytest-django` is NOT installed
**Files:** `app/tests/test_filters.py:20`, `app/tests/test_vcf_database.py:28`, `app/tests/test_views_additional.py:35`

`pytestmark = pytest.mark.django_db` requires `pytest-django` to enforce DB transaction isolation and Django setup. Without it:
- `manage.py test` ignores all pytest-style test files → these 422 pytest tests are **never run** by CI.
- Running `pytest` directly raises `PytestUnknownMarkWarning` — the mark is silently no-op, meaning tests run without proper Django test DB setup and may silently share state or fail differently than expected.

Coverage numbers for `app/filters.py` (66%), `app/services/vcf_database.py` (79%), and `app/views.py` (88%) are misleading: lines covered by skipped pytest tests aren't counted. Real coverage in untested pytest-paths is lower.

**Fix:** add `pytest-django` to `requirements.txt`; add `conftest.py` at project root with `django_settings` fixture, or ensure `pytest.ini`/`setup.cfg` has `DJANGO_SETTINGS_MODULE`.

#### COV2 — `app/filters.py` at 66% is the lowest-covered production module
**File:** `app/filters.py` (173 stmts, 59 missed)
Primary search feature. 59 uncovered statements in filter logic that interprets range queries on INFO fields. The coverage hole predates the pytest-django gap; even with pytest-django installed, range filter branches are not fully exercised.

#### COV3 — `test_import_data_view.py` is an empty stub
**File:** `app/tests/test_import_data_view.py`
Contains only `from __future__ import annotations` (3 lines). No tests. `ImportDataView` HTTP-level behavior (form rendering, error message wiring, file upload rejection) is not covered by the Django test runner.

---

### STYLE / ARCHITECTURE

#### ST1 — Module-level side-effect in `logging_utils.py`: opens log file on import
**File:** `MetagapUserCode/merge_vcf/logging_utils.py:38-43`
```python
fh = logging.FileHandler(LOG_FILE)   # opens "script_execution.log" in CWD at import time
logger.addHandler(fh)
...
configure_logging()   # called at module bottom (line 166)
```
Importing this package creates `script_execution.log` in whatever the current working directory is. Tests, Django imports, or any tool that imports `merge_vcf` will litter the filesystem. `configure_logging()` is idempotent and re-entrant, but the module-level handler addition before it is not guarded.

**Fix:** remove the top-level `fh = logging.FileHandler(...)` + `logger.addHandler(fh)` block (lines 38-43); rely solely on `configure_logging()` at line 166.

#### ST2 — Dual source of truth: module-level dicts duplicate YAML config
**File:** `app/services/vcf_metadata.py:183-366`
`METADATA_SECTION_MAP`, `METADATA_MODEL_MAP`, `METADATA_FIELD_ALIASES`, `SECTION_PRIMARY_FIELD` are defined both as Python dicts in `vcf_metadata.py` AND in `metadata_fields.yaml`. `VCFMetadataParser` uses the YAML-loaded `MetadataConfiguration`. Only `METADATA_SECTION_MAP` is actually imported elsewhere (`vcf_file_utils.py:17`). All others are unused dead exports.

**Fix:** remove the unused module-level dicts; import `METADATA_SECTION_MAP` equivalent from `load_metadata_configuration()` in `vcf_file_utils.py`.

#### ST3 — `vcf_database.py` is a 1379-line God object
`VCFDatabaseWriter` conflates sample-group creation logic (metadata extraction/coercion/consuming, 800+ lines) with low-level variant record persistence (`create_allele_frequency`, `create_info_instance`, `create_format_instance`). These responsibilities could be split into `MetadataBinder` and `VariantWriter` classes for readability and testability.

Not blocking, but worth tracking.

#### ST4 — Mixed indentation in `settings.py`
**File:** `MetaGap/MetaGap/settings.py:75-84`
MIDDLEWARE list uses tabs; rest of file uses 4-space indentation. Python accepts both since Python 3.x doesn't mix them in the same block, but it breaks style checks and editors configured for spaces.

#### ST5 — `django_tables2` deprecated accessor syntax
**File:** `app/tables.py:88`
```
DeprecationWarning: Use '__' to separate path components, not '.' in accessor 'format.genotype'
```
Two columns use dot notation (`format.genotype`, `format.payload`). Will break in django_tables2 v3.

**Fix:** `accessor="format__genotype"` and `accessor="format__payload"`.

---

### NITS

#### N1 — `export_sample_group_variants` uses `getattr(request.user, "organization_profile")` to access RelatedObjectDoesNotExist
**File:** `app/views.py:243-245`
```python
try:
    organization_profile = getattr(request.user, "organization_profile")
except ObjectDoesNotExist:
    organization_profile = None
```
`getattr` without a default won't raise `ObjectDoesNotExist` — it returns the related manager descriptor. `OrganizationProfile.DoesNotExist` is a subclass of `ObjectDoesNotExist` but is raised only when accessing `.organization_profile` through direct attribute access, not through `getattr`. The `try/except` may silently not catch if the accessor is lazy. Should use `getattr(request.user, "organization_profile", None)` + guard.

Actually the `RelatedObjectDoesNotExist` IS raised on attribute access when the OneToOne has no row — `getattr` without default propagates it. The `except ObjectDoesNotExist` does catch it correctly. Minor style concern only: `getattr(user, "organization_profile")` is idiomatic Django but is cleaner written as `user.organization_profile`.

#### N2 — `_po_count.py` is untracked stray at repo root
**File:** `MetaGap/_po_count.py` (untracked, `git status`)
Appears to be a helper script. Should be moved to `MetagapUserCode/` or committed/deleted.

#### N3 — `docs/overnight/` untracked directory with 7 files
`git status` shows 7 untracked `docs/overnight/*.md` files. Either add to `.gitignore` or commit.

---

## Top 5 Escalated Issues — Draft Patches

### PATCH 1 — Fix weak SECRET_KEY/DEBUG defaults (S1)

```diff
--- a/MetaGap/MetaGap/settings.py
+++ b/MetaGap/MetaGap/settings.py
@@ -1,6 +1,6 @@
+from django.core.exceptions import ImproperlyConfigured
+
 DEBUG = os.getenv("DEBUG", "0") == "1"
-SECRET_KEY = os.getenv("SECRET_KEY", "dev-secret")
+_secret_key = os.getenv("SECRET_KEY")
+if not _secret_key and not DEBUG:
+    raise ImproperlyConfigured(
+        "SECRET_KEY environment variable must be set in production (DEBUG=0)."
+    )
+SECRET_KEY = _secret_key or "dev-only-insecure-key-change-me"
 ALLOWED_HOSTS = _split_env_list(
-    "ALLOWED_HOSTS", "localhost,127.0.0.1,0.0.0.0,*.github.dev"
+    "ALLOWED_HOSTS", "localhost,127.0.0.1"
 )
```

### PATCH 2 — Fix broken `get_or_create_sample_group` (C2)

Simplest fix: remove or mark as unsupported. It's labeled "for testing" but is broken and exported.

```diff
--- a/MetaGap/MetaGap/app/services/vcf_database.py
+++ b/MetaGap/MetaGap/app/services/vcf_database.py
@@ -188,29 +188,0 @@
-def get_or_create_sample_group(
-    sample_group_name: str,
-    sample_group_description: str,
-    user_organization_profile: Any,
-    metadata: Dict[str, Any],
-) -> Optional[SampleGroup]:
-    """Get or create a SampleGroup with the given metadata.
-    ...
-    Note:
-        This is a simplified helper for testing. Production code should use
-        VCFDatabaseWriter.create_sample_group() instead.
-    """
-    try:
-        sample_group = SampleGroup.objects.create(
-            name=sample_group_name,
-            description=sample_group_description,
-            organization_profile=user_organization_profile,
-        )
-        return sample_group
-    except Exception:
-        return None

 # Also remove from __all__:
-    "get_or_create_sample_group",
```

### PATCH 3 — Fix double queryset in `SampleGroupTableView` (C3)

```diff
--- a/MetaGap/MetaGap/app/views.py
+++ b/MetaGap/MetaGap/app/views.py
@@ -113,8 +113,8 @@ class SampleGroupTableView(ListView):
     def get_context_data(self, **kwargs):
         context = super().get_context_data(**kwargs)
         table_class = create_dynamic_table(
             SampleGroup, table_name="SampleGroupTable", include_related=True
         )
-        table = table_class(self.get_queryset())
+        table = table_class(context["object_list"])
         RequestConfig(self.request, paginate={"per_page": self.paginate_by}).configure(table)
         context["table"] = table
         return context
```

### PATCH 4 — Fix `AlleleFrequency.chrom` max_length (C4)

```diff
--- a/MetaGap/MetaGap/app/models.py
+++ b/MetaGap/MetaGap/app/models.py
@@ -647,1 +647,1 @@
-    chrom = models.CharField(max_length=10)
+    chrom = models.CharField(max_length=64)

# Also fix ref/alt while here:
-    ref = models.CharField(max_length=50)
-    alt = models.CharField(max_length=50)
+    ref = models.CharField(max_length=255)
+    alt = models.CharField(max_length=255)
```
Requires a new migration: `python manage.py makemigrations app`.

### PATCH 5 — Fix syntax error in `test_merging.py` (C5)

```diff
--- a/MetaGap/MetaGap/MetagapUserCode/tests/test_merging.py
+++ b/MetaGap/MetaGap/MetagapUserCode/tests/test_merging.py
@@ -329,5 +329,5 @@
     monkeypatch.setattr(merging.vcfpy, "Call", _StubCall)
 
     merged = merging._merge_colliding_records(
-        [(rec1, 0), (rec2, 1)],          # parenthesis opened but never closed
+        [(rec1, 0), (rec2, 1)],
+    )
 def test_merge_colliding_records_pads_missing_samples():
```

Note: line 333 starts a new function without closing the prior call. The closing `)` must be inserted before the `def`. Exact fix depends on the intended arguments — the function signature must be inspected to determine what follows `[(rec1, 0), (rec2, 1)]`.

---

## Additional Recommendations (Not Patched Here)

- **Install `pytest-django`**: add to `requirements.txt`, configure `DJANGO_SETTINGS_MODULE=MetaGap.settings` in `pytest.ini`. Without this, 422 pytest tests are invisible to CI.
- **Remove module-level logging side-effect** (`logging_utils.py:38-43`): move handler setup into `configure_logging()` only.
- **Fix django_tables2 deprecated accessors** (`tables.py:88`): `format.genotype` → `format__genotype`.
- **Delete or move `MetaGap/_po_count.py`** to `MetagapUserCode/` or add to `.gitignore`.
- **Fill `test_import_data_view.py`** stub with actual HTTP-level import tests.
- **Remove duplicate module-level dicts** from `vcf_metadata.py` (ST2) — after confirming `vcf_file_utils.py` loads `METADATA_SECTION_MAP` from `load_metadata_configuration()` instead.

---

## Artifacts

- `docs/reviews/2026-05-07-code-review.md` — this report
