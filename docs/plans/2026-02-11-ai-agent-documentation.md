# AI Agent Documentation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Create structured reference documentation for AI agents in `.claude/docs/` covering project overview, models, workflows, and architecture.

**Architecture:** Four focused markdown files replace the monolithic CLAUDE.md approach. Each file documents a distinct layer of the system with precise file paths, class names, and relationships that AI agents need to navigate the codebase effectively.

**Tech Stack:** Django 5.1, Python 3.11

---

## Task 1: Create project-overview.md

**Files:**
- Create: `.claude/docs/project-overview.md`

**Step 1: Create directory structure**

```bash
mkdir -p .claude/docs
```

**Step 2: Write project-overview.md**

Create `.claude/docs/project-overview.md`:

```markdown
# MetaGap Project Overview

## What is MetaGap?

MetaGap is a Django 5.1 web application for managing sample group metadata and allele frequency information from genomic sequencing workflows. It ingests VCF (Variant Call Format) files, extracts rich metadata from custom header sections, and provides search/filtering capabilities for variant data.

**Core capabilities:**
- VCF file import with metadata extraction from custom headers
- Sample group management with 13 metadata model types
- Allele frequency search with INFO field range filters
- CSV/TSV export of variant data
- Multi-language support (English/Russian)

## Tech Stack

**Framework:** Django 5.1 (Python 3.11)

**Database:** SQLite (default), PostgreSQL/MySQL via `DATABASE_URL` env var

**Key dependencies:**
- `pysam` / `vcfpy` — VCF parsing
- `django-filter` — Search filtering
- `django-tables2` — Dynamic table generation
- `django-bootstrap5` — UI components
- `dj-database-url` — Database configuration
- `whitenoise` — Static file serving
- `gunicorn` — WSGI server (production)

**Frontend:** Bootstrap 5, DataTables.js, vanilla JavaScript

## Project Structure

```
MetaGap/                          # Django project root (working directory)
├── manage.py
├── requirements.txt
├── MetaGap/                      # Settings package
│   ├── settings.py               # All configuration
│   ├── urls.py                   # Root routing (i18n_patterns wrapper)
│   └── context_processors.py     # SITE_NAME/SITE_COLORS
├── app/                          # Main Django application
│   ├── models.py                 # 13 data models
│   ├── views.py                  # All views (CBVs)
│   ├── forms.py                  # BootstrapFormMixin + form classes
│   ├── filters.py                # django-filter search configs
│   ├── tables.py                 # django-tables2 table builders
│   ├── urls.py                   # App URL patterns
│   ├── signals.py                # Auto-create OrganizationProfile
│   ├── services/                 # VCF import pipeline
│   │   ├── vcf_importer.py       # Orchestrator
│   │   ├── vcf_metadata.py       # YAML config loader + parser
│   │   ├── vcf_database.py       # Database writer
│   │   ├── vcf_file_utils.py     # Text fallback parser
│   │   └── import_exceptions.py  # Custom exceptions
│   ├── config/
│   │   └── metadata_fields.yaml  # Metadata mapping config
│   ├── tests/                    # 20+ test modules, 320+ tests
│   └── fixtures/demo_data.json   # Demo users + sample groups
└── MetagapUserCode/              # Standalone VCF merge CLI (separate from Django)
    ├── merge_vcf/cli.py
    └── tests/
```

## Environment Variables

Configured in `MetaGap/MetaGap/settings.py`:

| Variable | Default | Purpose |
|----------|---------|---------|
| `DEBUG` | `1` | Enable debug mode |
| `SECRET_KEY` | `dev-secret` | Django secret key |
| `DATABASE_URL` | `sqlite:///db.sqlite3` | Database connection |
| `ALLOWED_HOSTS` | `localhost,127.0.0.1,*.github.dev` | Allowed host headers |
| `CSRF_TRUSTED_ORIGINS` | `https://*.github.dev,...` | CSRF origins |
| `LOG_LEVEL` | `INFO` | Logging verbosity |
| `SITE_NAME` | `MetaGaP` | Branding name |
| `SITE_COLOR_PRIMARY` | `#1f6feb` | Primary theme color |
| `SITE_COLOR_SECONDARY` | `#2ea043` | Secondary theme color |
| `METAGAP_MAX_UPLOAD_SIZE_MB` | `100` | Max VCF upload size |

## Common Commands

```bash
# Run development server
python MetaGap/manage.py runserver

# Run all tests
python MetaGap/manage.py test

# Run specific test module
python MetaGap/manage.py test app.tests.test_vcf_importer

# Run tests with coverage
cd MetaGap && coverage run manage.py test && coverage report

# Run migrations
python MetaGap/manage.py migrate

# Load demo data (users: demo_org/demo_lab, password: metagap-demo)
python MetaGap/manage.py loaddata app/fixtures/demo_data.json

# Compile translations
cd MetaGap && DJANGO_SETTINGS_MODULE=MetaGap.settings django-admin compilemessages

# Collect static files (production)
python MetaGap/manage.py collectstatic
```

## Sub-Projects

**1. Django Web Application** (main project)
- Location: `MetaGap/app/`
- Purpose: VCF import, metadata management, search, export
- Entry point: `MetaGap/manage.py runserver`

**2. VCF Merge CLI** (standalone tool)
- Location: `MetaGap/MetagapUserCode/merge_vcf/`
- Purpose: Merge multiple VCF files with metadata conflict resolution
- Entry point: `python -m MetaGap.MetagapUserCode.merge_vcf.cli`
- **Not integrated into Django app** — separate test suite, logging config

## Testing

- **Framework:** pytest + Django TestCase
- **Test count:** 320+ tests across 20 modules
- **Fixtures:** `app/tests/conftest.py` (pytest), `app/fixtures/demo_data.json` (Django)
- **Coverage targets (CI):** `app/` ≥ 70%, `MetagapUserCode/` ≥ 75%

## CI/CD

GitHub Actions (`.github/workflows/ci.yml`) runs on push to `main` and PRs:
1. Python 3.11 setup
2. Install dependencies + coverage
3. Run `coverage run manage.py test`
4. Enforce coverage thresholds
5. Scan logs for "Unhandled metadata key" warnings → CI annotations
6. Upload test artifacts
```

**Step 3: Commit**

```bash
git add .claude/docs/project-overview.md
git commit -m "docs: add project overview for AI agents"
```

---

## Task 2: Create models.md

**Files:**
- Create: `.claude/docs/models.md`

**Step 1: Write models.md**

Create `.claude/docs/models.md`:

```markdown
# MetaGap Data Models

All models defined in `MetaGap/app/models.py`.

## 1. User & Organization

### OrganizationProfile

OneToOne relationship with Django `User`. Required for VCF imports.

**Fields:**
- `user`: OneToOneField(User, CASCADE, related_name="organization_profile")
- `organization_name`: CharField(255, blank/null)

**Business logic:**
- Auto-created via `post_save` signal in `app/signals.py:create_user_profile`
- Signal checks `kwargs.get("raw", False)` to skip fixture loading
- `VCFImporter.import_file()` raises `ImporterConfigurationError` if missing

**Related by:**
- `SampleGroup.created_by` → OrganizationProfile (CASCADE)

---

## 2. SampleGroup (Central Hub)

### SampleGroup

The main entity. Links all metadata models together.

**Fields:**
- `name`: CharField(255)
- `created_by`: ForeignKey(OrganizationProfile, CASCADE)
- `doi`: CharField(100, blank/null)
- `source_lab`: CharField(255, blank/null)
- `contact_email`: EmailField(blank/null)
- `contact_phone`: CharField(50, blank/null)
- `total_samples`: IntegerField(blank/null)
- `inclusion_criteria`: TextField(blank/null)
- `exclusion_criteria`: TextField(blank/null)
- `comments`: TextField(blank/null)
- `additional_metadata`: JSONField(blank/null) — catches unmapped VCF headers
- `sequencing_platform`: CharField(32, choices=SequencingPlatform.choices, blank/null)

**Metadata relationships (all SET_NULL):**
- `reference_genome_build`: ForeignKey(ReferenceGenomeBuild)
- `genome_complexity`: ForeignKey(GenomeComplexity)
- `sample_origin`: ForeignKey(SampleOrigin)
- `material_type`: ForeignKey(MaterialType)
- `library_construction`: ForeignKey(LibraryConstruction)
- `input_quality`: OneToOneField(InputQuality)

**Sequencing platform relationships (all SET_NULL):**
- `illumina_seq`: ForeignKey(IlluminaSeq)
- `ont_seq`: ForeignKey(OntSeq)
- `pacbio_seq`: ForeignKey(PacBioSeq)
- `iontorrent_seq`: ForeignKey(IonTorrentSeq)

**Bioinformatics relationships (all SET_NULL):**
- `bioinfo_alignment`: ForeignKey(BioinfoAlignment)
- `bioinfo_variant_calling`: ForeignKey(BioinfoVariantCalling)
- `bioinfo_post_proc`: ForeignKey(BioinfoPostProc)

**Reverse relationships:**
- `allele_frequencies` ← AlleleFrequency (CASCADE)

**Business logic:**
- `delete()` method (lines 595-636): Deletes unshared metadata instances
  - Iterates through all FK/OneToOne fields
  - Checks if metadata instance is shared by other SampleGroups
  - Only deletes if not shared (protects reusable metadata)
  - Always deletes OneToOneField instances (InputQuality)
- `get_active_sequencing_platform()`: Returns (label, instance) tuple
  - Checks `additional_metadata["active_sequencing_platform"]` first
  - Falls back to first non-null platform FK
- `_normalize_platform_key()`: Normalizes platform identifiers (snake_case, lowercase)

**Class attributes:**
- `SEQUENCING_PLATFORM_FIELDS`: Tuple of (field_name, label) pairs
- `SequencingPlatform`: TextChoices enum
- `PLATFORM_FIELD_MAP`: Maps enum values to field names

---

## 3. Metadata Models

All models use SET_NULL on delete (except InputQuality which uses CASCADE via OneToOne).

### ReferenceGenomeBuild

**Fields:**
- `build_name`: CharField(100)
- `build_version`: CharField(100, blank/null)
- `additional_info`: JSONField(blank/null)

### GenomeComplexity

**Fields:**
- `size`: CharField(50, blank/null)
- `ploidy`: CharField(50, blank/null)
- `gc_content`: CharField(50, blank/null)

### SampleOrigin

**Fields:**
- `tissue`: CharField(100, blank/null)
- `collection_method`: CharField(100, blank/null)
- `storage_conditions`: CharField(100, blank/null)
- `time_stored`: CharField(100, blank/null)

### MaterialType

**Fields:**
- `material_type`: CharField(50, choices=MATERIAL_CHOICES, blank/null)
  - Choices: DNA, RNA, cDNA
- `integrity_number`: CharField(50, blank/null)

### LibraryConstruction

**Fields:**
- `kit`: CharField(100, blank/null)
- `fragmentation`: CharField(100, blank/null)
- `adapter_ligation_efficiency`: CharField(50, blank/null)
- `pcr_cycles`: IntegerField(blank/null)

### InputQuality

**Fields:**
- `a260_a280`: FloatField(blank/null)
- `a260_a230`: FloatField(blank/null)
- `dna_concentration`: FloatField(blank/null)
- `rna_concentration`: FloatField(blank/null)
- `notes`: TextField(blank/null)
- `additional_metrics`: JSONField(blank/null)

**Business logic:**
- OneToOne with SampleGroup (always deleted with group)

---

## 4. Sequencing Platform Models

All inherit from abstract `SequencingInstrument` base class.

### SequencingInstrument (abstract)

**Fields:**
- `instrument`: CharField(100)
- `flow_cell`: CharField(100, blank/null)
- `additional_info`: JSONField(blank/null)

### IlluminaSeq

**Additional fields:**
- `channel_method`: CharField(50, blank/null)
- `cluster_density`: CharField(50, blank/null)
- `qc_software`: CharField(100, blank/null)

### OntSeq (Oxford Nanopore)

**Additional fields:**
- `flow_cell_version`: CharField(50, blank/null)
- `pore_type`: CharField(50, blank/null)
- `bias_voltage`: CharField(50, blank/null)

### PacBioSeq

**Additional fields:**
- `smrt_cell_type`: CharField(50, blank/null)
- `zmw_density`: CharField(50, blank/null)

### IonTorrentSeq

**Additional fields:**
- `chip_type`: CharField(50, blank/null)
- `ph_calibration`: CharField(50, blank/null)
- `flow_order`: CharField(50, blank/null)
- `ion_sphere_metrics`: CharField(50, blank/null)

---

## 5. Bioinformatics Models

### BioinfoAlignment

**Fields:**
- `tool`: CharField(100, blank/null)
- `params`: CharField(255, blank/null)
- `ref_genome_version`: CharField(100, blank/null)
- `recalibration_settings`: CharField(255, blank/null)

### BioinfoVariantCalling

**Fields:**
- `tool`: CharField(100, blank/null)
- `version`: CharField(50, blank/null)
- `filtering_thresholds`: CharField(255, blank/null)
- `duplicate_handling`: CharField(50, blank/null)
- `mq`: CharField(50, blank/null)

### BioinfoPostProc

**Fields:**
- `normalization`: CharField(100, blank/null)
- `harmonization`: CharField(100, blank/null)

---

## 6. Variant Data Models

### AlleleFrequency

Stores variant records from VCF files.

**Fields:**
- `sample_group`: ForeignKey(SampleGroup, CASCADE, related_name="allele_frequencies")
- `chrom`: CharField(10)
- `pos`: IntegerField
- `variant_id`: CharField(100, blank/null)
- `ref`: CharField(50)
- `alt`: CharField(50)
- `qual`: FloatField(blank/null)
- `filter`: CharField(50, blank/null)
- `info`: ForeignKey(Info, CASCADE, blank/null)
- `format`: ForeignKey(Format, CASCADE, blank/null)
- `comments`: TextField(blank/null)

**Constraints:**
- Unique: (sample_group, chrom, pos, ref, alt)

**Indexes:**
- (sample_group, chrom, pos)
- variant_id

### Info

Stores VCF INFO field data.

**Fields (VCF standard):**
- `aa`, `ac`, `af`, `an`, `bq`, `cigar`, `db`, `dp`, `end`, `h2`, `h3`, `mq`, `mq0`, `qd`, `fs`, `sor`, `ns`, `sb`: CharField(50, blank/null)
- `additional`: JSONField(blank/null) — unmapped INFO fields

**Business logic:**
- `save()` calls `_normalize_placeholder_values()` before persist
- `_normalize_placeholder_values()`:
  - Converts placeholder strings (".", "") to None for numeric fields
  - Introspects model to find IntegerField/FloatField/DecimalField
  - Also normalizes numeric values in `additional` JSONField
  - Removes None entries from `additional`

**Constants:**
- `PLACEHOLDER_STRINGS = {".", ""}`
- `NUMERIC_FIELD_TYPES = (IntegerField, FloatField, DecimalField)`

### Format

Stores VCF FORMAT field data.

**Fields:**
- `genotype`: CharField(50, blank/null)
- `payload`: JSONField(blank/null)

**Properties:**
- `fields`: Returns `payload["fields"]` dict or empty dict
- `additional`: Returns `payload["additional"]` dict or empty dict

---

## Utility Functions

### _format_attributes(*attributes)

**Location:** `models.py:110-118`

**Purpose:** Format model attributes for `__str__()` methods

**Logic:**
- Takes tuple pairs: (label, value)
- Filters out None/"" values
- Returns comma-separated "Label: Value" string
- Returns "Not provided" if all empty
```

**Step 2: Commit**

```bash
git add .claude/docs/models.md
git commit -m "docs: add data models reference for AI agents"
```

---

## Task 3: Create workflows.md

**Files:**
- Create: `.claude/docs/workflows.md`

**Step 1: Write workflows.md**

Create `.claude/docs/workflows.md`:

```markdown
# MetaGap User Workflows

User-facing flows through the Django application. All views in `MetaGap/app/views.py`.

## 1. User Registration

**Flow:**
1. User visits `/signup/` → `UserRegistrationView` (CreateView)
2. Form: `CustomUserCreationForm` (`app/forms.py:30-60`)
   - Validates unique email
   - Validates password strength (Django validators)
3. User created → `post_save` signal fires (`app/signals.py:8-15`)
4. `OrganizationProfile.objects.create(user=instance)` auto-creates profile
5. Redirect to `/login/`

**Views:**
- `UserRegistrationView` (lines 280-290)
  - Template: `signup.html`
  - Success URL: `reverse_lazy("login")`

**Forms:**
- `CustomUserCreationForm` (lines 30-60)
  - Fields: username, email, password1, password2
  - Clean method: checks `User.objects.filter(email__iexact=email).exists()`

**Signals:**
- `create_user_profile` (signals.py:8-15)
  - Skips if `kwargs.get("raw", False)` (fixture loading)
  - Creates OrganizationProfile only if `created=True`

---

## 2. VCF Import

**Flow:**
1. User visits `/import/` (LoginRequired) → `ImportDataView` (FormView)
2. Form: `ImportDataForm` validates file extension (.vcf, .vcf.gz)
3. File saved to `MEDIA_ROOT` via `default_storage.save()`
4. `VCFImporter(user).import_file(file_path)` orchestrates:
   - **Step A:** Validate `user.organization_profile` exists
     - Raises `ImporterConfigurationError` if missing
   - **Step B:** Try pysam parser
     - Calls `_import_with_pysam()` which:
       - Opens VCF with `pysam.VariantFile()`
       - Extracts metadata headers (`##SAMPLE_GROUP`, `##REFERENCE_GENOME_BUILD`, etc.)
       - Calls `VCFMetadataParser.parse()` to normalize headers
       - Calls `VCFDatabaseWriter.create_sample_group()` to create models
       - Iterates variants, calls `VCFDatabaseWriter.write_allele_frequencies()`
   - **Step C:** Fallback on pysam failure
     - If gzipped, inflates to temp file and retries pysam
     - If still fails, warns and falls back to `parse_vcf_text_fallback()`
     - Text parser manually extracts headers and variant rows
   - **Step D:** Cleanup
     - Removes temp files
     - Deletes partial SampleGroup on error
5. Success: redirect to `/profile/`
6. Error: display `ImporterValidationError` message to user

**Views:**
- `ImportDataView` (lines 300-380)
  - Template: `import_data.html`
  - Form: `ImportDataForm`
  - `form_valid()`: saves file, calls importer, redirects

**Services:**
- `VCFImporter` (`services/vcf_importer.py:39-180`)
  - `__init__(user)`: stores user, creates parser/writer
  - `import_file(file_path)`: main entry point
  - `_import_with_pysam()`: pysam parsing
  - `_inflate_gzip_to_temp()`: decompress .gz files
  - `_extract_metadata_text_fallback()`: manual header extraction

- `VCFMetadataParser` (`services/vcf_metadata.py:150-400`)
  - `parse(vcf_metadata)`: normalizes VCF headers via YAML config
  - Loads `app/config/metadata_fields.yaml`:
    - `metadata_section_map`: VCF section → snake_case identifier
    - `metadata_models`: identifier → Django model class path
    - `metadata_field_aliases`: canonical field + legacy aliases
    - `section_primary_field`: fallback field per section
  - Logs warnings for unmapped keys
  - Stores unmapped sections in `additional_metadata` JSONField

- `VCFDatabaseWriter` (`services/vcf_database.py:50-500`)
  - `create_sample_group(metadata, org_profile)`: creates models
  - `write_allele_frequencies(variants, sample_group)`: bulk creates
  - For each variant:
    - Creates `Info` instance (normalizes placeholders)
    - Creates `Format` instance (stores genotype + JSON)
    - Creates `AlleleFrequency` instance
    - Enforces unique constraint on (sample_group, chrom, pos, ref, alt)

**Exceptions:**
- `ImporterConfigurationError` — missing org profile, bad config
- `ImporterValidationError` — VCF parse error, invalid format
- `ImporterError` — generic fallback

**Error handling:**
- `try/except` wraps entire import in `transaction.atomic()`
- Partial SampleGroup deleted on any exception
- Temp files cleaned up in `finally` block

---

## 3. Search

**Flow:**
1. User visits `/` → `HomePageView` (TemplateView)
   - Template: `home.html`
   - Context: `SearchForm` (empty form)
2. User submits search form → POST to `/search/`
3. `SearchResultsView` (FilterView + SingleTableMixin)
   - Filter: `AlleleFrequencySearchFilter` (`app/filters.py:30-150`)
   - Table: `build_allele_frequency_table()` (`app/tables.py:200-300`)
4. Filtering:
   - Query: universal search across chrom/pos/ref/alt
   - Chromosome: exact match
   - Position: exact, min, max
   - Reference/alternate: exact match
   - QUAL: min, max
   - Filter pass: boolean
   - INFO fields: AF, DP, MQ, QD, FS, SOR (min/max ranges)
   - Related models: source_lab, tissue, alignment tool, variant caller
5. Table generation:
   - `create_dynamic_table()` introspects `AlleleFrequency` model
   - Includes related model columns (SampleGroup, Info)
   - Applies CSS classes (`text-end` for numeric columns)
   - Configures column ordering and exclusions
6. Pagination:
   - DataTables.js client-side (10 rows/page)
   - Django pagination disabled (all results sent to JS)

**Views:**
- `HomePageView` (lines 100-110)
  - Template: `home.html`
  - Context: `{"form": SearchForm()}`

- `SearchResultsView` (lines 120-180)
  - Inherits: `SingleTableMixin`, `FilterView`
  - `filterset_class = AlleleFrequencySearchFilter`
  - `get_table()`: calls `build_allele_frequency_table()`
  - `get_queryset()`: `select_related("sample_group", "info", "format")`

**Filters:**
- `AlleleFrequencySearchFilter` (filters.py:30-150)
  - Uses django-filter `FilterSet`
  - Field types: CharFilter, NumberFilter, RangeFilter, BooleanFilter
  - `filter_query()` method: combines Q objects for universal search

**Tables:**
- `build_allele_frequency_table()` (tables.py:200-300)
  - Returns django-tables2 `Table` class dynamically
  - Calls `create_dynamic_table(AlleleFrequency, ...)`
  - Priority columns: chrom, pos, ref, alt, qual, filter
  - Excludes: id, comments, format

---

## 4. Sample Group CRUD

All views enforce ownership via `OrganizationSampleGroupMixin` (`app/mixins.py`).

### Detail View

**Flow:**
1. User visits `/profile/sample-groups/<pk>/` → `SampleGroupDetailView`
2. Ownership check: `get_queryset()` filters by `created_by=user.organization_profile`
3. Query optimization: `select_related()` all FK/OneToOne metadata fields
4. Template displays all metadata + variant count

**Views:**
- `SampleGroupDetailView` (lines 400-420)
  - Inherits: `LoginRequiredMixin`, `OrganizationSampleGroupMixin`, `DetailView`
  - Template: `sample_group_detail.html`
  - Context: sample group + related metadata

### Update View

**Flow:**
1. User visits `/profile/sample-groups/<pk>/edit/` → `SampleGroupUpdateView`
2. Ownership check via mixin
3. Form: `SampleGroupForm` (BootstrapFormMixin + ModelForm)
4. On submit: saves changes, redirects to detail view

**Views:**
- `SampleGroupUpdateView` (lines 430-450)
  - Inherits: `LoginRequiredMixin`, `OrganizationSampleGroupMixin`, `UpdateView`
  - Form: `SampleGroupForm`
  - Success URL: `sample_group_detail`

### Delete View

**Flow:**
1. User visits `/profile/sample-groups/<pk>/delete/` → `SampleGroupDeleteView`
2. Ownership check via mixin
3. Confirmation page
4. On confirm: `SampleGroup.delete()` called
   - Custom delete logic (models.py:595-636):
     - Iterates FK/OneToOne metadata fields
     - Checks if metadata shared by other groups
     - Deletes unshared metadata
     - Always deletes OneToOne (InputQuality)
5. Redirect to profile

**Views:**
- `SampleGroupDeleteView` (lines 460-480)
  - Inherits: `LoginRequiredMixin`, `OrganizationSampleGroupMixin`, `DeleteView`
  - Template: `sample_group_confirm_delete.html`
  - Success URL: `profile`

### Profile List

**Flow:**
1. User visits `/profile/` → `ProfileView`
2. Query: `SampleGroup.objects.filter(created_by=user.organization_profile)`
3. Optimization: `select_related()` all metadata, `prefetch_related("allele_frequencies")`
4. Template displays table of owned groups

**Views:**
- `ProfileView` (lines 350-380)
  - Inherits: `LoginRequiredMixin`, `TemplateView`
  - Template: `profile.html`
  - Context: `{"sample_groups": queryset}`

**Mixins:**
- `OrganizationSampleGroupMixin` (mixins.py:10-30)
  - Overrides `get_queryset()`: filters by `created_by`
  - Raises `Http404` if not owner

---

## 5. Export

**Flow:**
1. User visits `/profile/sample-groups/<pk>/export/` → `export_sample_group_variants`
2. Ownership check: `get_object_or_404(SampleGroup, pk=pk, created_by=org_profile)`
3. Query: `AlleleFrequency.objects.filter(sample_group=group).select_related("info", "format")`
4. Format: CSV or TSV (defaults to CSV)
5. Streaming response:
   - Headers: chrom, pos, variant_id, ref, alt, qual, filter, INFO fields, FORMAT fields
   - Rows: one per variant
   - Filename: `{slugify(group_name)}-variants.csv`

**Views:**
- `export_sample_group_variants` (lines 500-600)
  - Function view (not CBV)
  - Decorator: `@login_required`
  - Ownership check: manual `get_object_or_404()` with `created_by` filter
  - Returns: `HttpResponse` with `Content-Type: text/csv`
  - Uses Python `csv.writer()` for streaming

**Query optimization:**
- `select_related("info", "format")` to avoid N+1 queries
- Streams results (no pagination, all rows exported)

**Filename format:**
- `slugify(sample_group.name)` + `-variants.` + `csv|tsv`
- Example: `lung-tumor-cohort-variants.csv`
```

**Step 2: Commit**

```bash
git add .claude/docs/workflows.md
git commit -m "docs: add user workflows reference for AI agents"
```

---

## Task 4: Create architecture.md

**Files:**
- Create: `.claude/docs/architecture.md`

**Step 1: Write architecture.md**

Create `.claude/docs/architecture.md`:

```markdown
# MetaGap Architecture

## Request Lifecycle

```
HTTP Request
    ↓
Root URL Router (MetaGap/urls.py)
    ↓ i18n_patterns() wrapper
    ↓
Middleware Stack
    ├── SecurityMiddleware
    ├── WhiteNoiseMiddleware (static files)
    ├── SessionMiddleware
    ├── LocaleMiddleware (i18n)
    ├── CommonMiddleware
    ├── CsrfViewMiddleware
    ├── AuthenticationMiddleware
    ├── MessageMiddleware
    └── XFrameOptionsMiddleware
    ↓
App URL Router (app/urls.py)
    ↓
View Layer (CBVs + mixins)
    ↓ (only for VCF import)
    ↓
Service Layer (services/ package)
    ↓
Model Layer (ORM)
    ↓
Template Rendering (Django templates)
    ↓ django-bootstrap5 + django-tables2
    ↓
HTTP Response
```

**What doesn't exist:**
- No REST API / DRF serializers
- No Celery / async task queue
- No Redis / caching layer
- No GraphQL endpoint
- No WebSockets

## Layer Map

### 1. Configuration Layer

**Files:**
- `MetaGap/settings.py` — all Django settings
- `MetaGap/context_processors.py` — SITE_NAME/SITE_COLORS injection
- `app/config/metadata_fields.yaml` — VCF import config

**Responsibilities:**
- Environment variable parsing (`os.getenv()`)
- Database configuration (`dj_database_url`)
- Middleware ordering
- Installed apps
- Logging setup
- I18n settings (LANGUAGES, LOCALE_PATHS)
- Static/media file paths

**Dependencies:**
- Python stdlib (`os`, `pathlib`)
- `dj_database_url` — DB URL parsing
- `yaml` — metadata config loading

### 2. Routing Layer

**Files:**
- `MetaGap/urls.py` — root URL config
- `app/urls.py` — app URL patterns

**Root patterns:**
```python
urlpatterns = [
    path("i18n/", include("django.conf.urls.i18n")),
]
urlpatterns += i18n_patterns(
    path("", include("app.urls")),
    path("admin/", admin.site.urls),
)
```

**App patterns (app/urls.py):**
- `/` → HomePageView
- `/dashboard/` → DashboardView
- `/search/` → SearchResultsView
- `/import/` → ImportDataView
- `/profile/` → ProfileView
- `/profile/sample-groups/<pk>/` → SampleGroupDetailView
- `/profile/sample-groups/<pk>/edit/` → SampleGroupUpdateView
- `/profile/sample-groups/<pk>/delete/` → SampleGroupDeleteView
- `/profile/sample-groups/<pk>/export/` → export_sample_group_variants
- `/signup/` → UserRegistrationView
- `/login/` → LoginView
- `/logout/` → LogoutView
- `/contact/` → ContactView
- `/about/` → AboutView

**i18n patterns:**
- All app URLs wrapped in `i18n_patterns()`
- Language prefix: `/<lang>/...` (e.g., `/en/profile/`, `/ru/profile/`)
- Language selection: POST to `/i18n/setlang/` with `language` param

### 3. View Layer

**Files:**
- `app/views.py` — all views (800+ lines)
- `app/mixins.py` — OrganizationSampleGroupMixin
- `app/forms.py` — all forms

**Pattern:** Class-based views (CBVs) with mixins

**Common mixins:**
- `LoginRequiredMixin` — enforces authentication
- `OrganizationSampleGroupMixin` — enforces ownership
- `SingleTableMixin` — django-tables2 integration
- `FilterView` — django-filter integration

**View types:**
- `TemplateView` — home, dashboard, profile, about, contact
- `FormView` — import data
- `CreateView` — user registration
- `DetailView` — sample group detail
- `UpdateView` — sample group edit, profile edit
- `DeleteView` — sample group delete, account delete
- `FilterView` — search results
- Function view — export (streaming response)

**Query optimization patterns:**
- `select_related()` for FK/OneToOne (avoid N+1)
- `prefetch_related()` for reverse FK (allele_frequencies)
- Applied in `get_queryset()` methods

### 4. Service Layer

**Files:**
- `app/services/vcf_importer.py` — orchestrator
- `app/services/vcf_metadata.py` — YAML config + parser
- `app/services/vcf_database.py` — database writer
- `app/services/vcf_file_utils.py` — text fallback parser
- `app/services/import_exceptions.py` — custom exceptions

**Purpose:** Only exists for VCF import. All other operations use views → models directly.

**Flow:**
```
ImportDataView
    ↓
VCFImporter.import_file()
    ↓ Try pysam
    ↓
VCFMetadataParser.parse()
    ↓ Loads metadata_fields.yaml
    ↓ Normalizes headers
    ↓
VCFDatabaseWriter.create_sample_group()
    ↓ Creates model instances
    ↓
VCFDatabaseWriter.write_allele_frequencies()
    ↓ Bulk creates variants
```

**Error handling:**
- Wraps entire import in `transaction.atomic()`
- Deletes partial SampleGroup on error
- Cleans up temp files in `finally` block

**Dependencies:**
- `pysam` — VCF parsing (C extension)
- `vcfpy` — alternative VCF parser (pure Python)
- `PyYAML` — config loading

### 5. Model Layer

**Files:**
- `app/models.py` — all 13 models
- `app/signals.py` — post_save signal for OrganizationProfile
- `app/admin.py` — Django admin registration

**Patterns:**
- Abstract base class: `SequencingInstrument`
- Custom `delete()` method: `SampleGroup` (protects shared metadata)
- Custom `save()` method: `Info` (normalizes placeholders)
- Signal: auto-create `OrganizationProfile` on User creation

**Relationships:**
- SampleGroup → OrganizationProfile (CASCADE)
- SampleGroup → metadata models (SET_NULL)
- SampleGroup ← AlleleFrequency (CASCADE)
- AlleleFrequency → Info (CASCADE)
- AlleleFrequency → Format (CASCADE)

**Constraints:**
- Unique: (sample_group, chrom, pos, ref, alt) on AlleleFrequency
- Indexes: (sample_group, chrom, pos), variant_id

### 6. Presentation Layer

**Files:**
- `app/forms.py` — BootstrapFormMixin + all forms
- `app/filters.py` — django-filter FilterSet classes
- `app/tables.py` — django-tables2 Table builders
- `app/templates/*.html` — Django templates

**BootstrapFormMixin:**
- Auto-applies Bootstrap CSS classes to widgets
- Sets data attributes for JS interop
- Applied to all forms via inheritance

**Dynamic table generation:**
- `create_dynamic_table()` introspects model fields
- Builds django-tables2 `Table` class at runtime
- Supports column exclusion, priority ordering
- Applies CSS classes based on field type

**Client-side:**
- Bootstrap 5 components
- DataTables.js for pagination
- Vanilla JavaScript (no jQuery dependency)

### 7. Static Files

**Middleware:** `whitenoise.middleware.WhiteNoiseMiddleware`

**Storage:** `whitenoise.storage.CompressedManifestStaticFilesStorage`

**Process:**
1. Development: served from `app/static/` directories
2. Production: `python manage.py collectstatic` → `MetaGap/staticfiles/`
3. WhiteNoise compresses and adds cache headers
4. Manifest file tracks versioned filenames

**Media files:**
- Location: `MetaGap/media/`
- Purpose: VCF uploads, exports
- Served by Django in dev, web server in production

### 8. Internationalization

**Middleware:** `django.middleware.locale.LocaleMiddleware`

**Languages:** English (en), Russian (ru)

**Process:**
1. Developer edits `.po` message catalogs in `app/locale/`
2. `django-admin compilemessages` generates `.mo` binaries
3. Auto-compilation on app startup via `AppConfig.ready()` (app/apps.py)
4. Language selection: POST to `/i18n/setlang/` sets `metagap_language` cookie

**Auto-compilation:**
```python
# app/apps.py
def ready(self):
    ensure_compiled_catalogs(Path(...) / 'locale')
    import app.signals
```

## Component Dependencies

### Core Django Components

**settings.py depends on:**
- `os`, `pathlib` (stdlib)
- `dj_database_url` (DB config)
- `django.utils.translation` (i18n)

**models.py depends on:**
- Django ORM (`django.db.models`)
- `django.contrib.auth.models.User`
- No business logic imports (self-contained)

**views.py depends on:**
- `app/models.py` (SampleGroup, AlleleFrequency)
- `app/forms.py` (all form classes)
- `app/filters.py` (AlleleFrequencySearchFilter)
- `app/tables.py` (build_allele_frequency_table)
- `app/mixins.py` (OrganizationSampleGroupMixin)
- `app/services/vcf_importer.py` (VCFImporter)
- `django_filters.views.FilterView`
- `django_tables2.views.SingleTableMixin`

**forms.py depends on:**
- Django forms (`django.forms`)
- `app/models.py` (SampleGroup, User)
- Bootstrap 5 (CSS classes only, no Python dep)

**filters.py depends on:**
- `django_filters` (FilterSet, filters)
- `app/models.py` (AlleleFrequency, SampleGroup)

**tables.py depends on:**
- `django_tables2` (Table, Column)
- Django ORM (introspection via `_meta`)

### Service Layer Dependencies

**vcf_importer.py depends on:**
- `pysam` (VCF parsing)
- `app/models.py` (SampleGroup, OrganizationProfile)
- `app/services/vcf_metadata.py` (VCFMetadataParser)
- `app/services/vcf_database.py` (VCFDatabaseWriter)
- `app/services/vcf_file_utils.py` (text fallback)
- `app/services/import_exceptions.py` (custom exceptions)

**vcf_metadata.py depends on:**
- `pysam` (VariantFile)
- `PyYAML` (config loading)
- `app/models.py` (all metadata models)
- `app/config/metadata_fields.yaml` (YAML config)

**vcf_database.py depends on:**
- `app/models.py` (SampleGroup, AlleleFrequency, Info, Format)
- Django ORM (get_or_create, bulk_create)

### Third-Party Integrations

| Library | Used By | Purpose |
|---------|---------|---------|
| `pysam` | services/vcf_importer.py, services/vcf_metadata.py | VCF file parsing (C extension wrapping htslib) |
| `vcfpy` | services/vcf_file_utils.py | Alternative VCF parser (pure Python fallback) |
| `PyYAML` | services/vcf_metadata.py | Load metadata_fields.yaml config |
| `django-filter` | filters.py, views.py | Search filtering with FilterSet classes |
| `django-tables2` | tables.py, views.py | Dynamic table generation |
| `django-bootstrap5` | forms.py, templates | Bootstrap 5 form rendering |
| `dj-database-url` | settings.py | Parse DATABASE_URL env var |
| `whitenoise` | settings.py (middleware) | Static file serving with compression |
| `gunicorn` | Production only | WSGI server (not imported in code) |

### Circular Dependency Prevention

**No circular imports:**
- Models import nothing from app/ (except in TYPE_CHECKING)
- Signals import models via `apps.get_model()` to avoid circular imports
- Services import models directly (one-way dependency)
- Views import everything (top of dependency tree)

**Import order:**
1. Stdlib
2. Django
3. Third-party
4. Local app (models → forms → filters → tables → views)

## Deployment Architecture

**Development:**
- Django dev server (`manage.py runserver`)
- SQLite database (`db.sqlite3`)
- Static files served by Django
- DEBUG=1

**Production:**
- Gunicorn WSGI server (systemd service: `gunicorn.service`)
- PostgreSQL/MySQL via DATABASE_URL
- WhiteNoise for static files
- Web server (nginx/Apache) for media files
- DEBUG=0
- LOG_LEVEL=INFO (full tracebacks logged)

**Environment config:**
- All settings via env vars (no code changes needed)
- `.env` file loaded by Django (not committed to git)
```

**Step 2: Commit**

```bash
git add .claude/docs/architecture.md
git commit -m "docs: add architecture reference for AI agents"
```

---

## Task 5: Final verification and commit

**Files:**
- Read: `.claude/docs/*.md`

**Step 1: Verify all files created**

```bash
ls -lh .claude/docs/
```

Expected output:
```
project-overview.md
models.md
workflows.md
architecture.md
```

**Step 2: Verify files are well-formed markdown**

```bash
head -20 .claude/docs/project-overview.md
head -20 .claude/docs/models.md
head -20 .claude/docs/workflows.md
head -20 .claude/docs/architecture.md
```

**Step 3: Final commit**

```bash
git add .claude/docs/
git commit -m "docs: complete AI agent documentation structure

- project-overview.md: tech stack, commands, environment
- models.md: all 13 models with fields, relationships, business logic
- workflows.md: 5 user-facing flows with code paths
- architecture.md: request lifecycle, component dependencies, integrations

Co-Authored-By: Claude Sonnet 4.5 <noreply@anthropic.com>"
```

**Step 4: Verify commit**

```bash
git log -1 --stat
```

Expected: Shows 4 files added in `.claude/docs/`
