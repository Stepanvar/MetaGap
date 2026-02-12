# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MetaGap is a Django 5.1 web application for managing sample group metadata and allele frequency information from sequencing workflows. It provides VCF file import/export, metadata management, and search capabilities for genomic variant data. Python 3.11.

## Commands

All commands run from `MetaGap/` subdirectory (where `manage.py` lives):

```bash
# Install dependencies
pip install -r MetaGap/requirements.txt

# Run development server
python MetaGap/manage.py runserver

# Run all tests
python MetaGap/manage.py test

# Run a single test module
python MetaGap/manage.py test app.tests.test_vcf_importer

# Run a specific test class
python MetaGap/manage.py test app.tests.test_import_helpers.MetadataConfigurationTests

# Run tests with coverage
cd MetaGap && coverage run manage.py test && coverage report

# Run migrations
python MetaGap/manage.py migrate

# Load demo data (creates demo_org/demo_lab users, password: metagap-demo)
python MetaGap/manage.py loaddata app/fixtures/demo_data.json

# Collect static files (production)
python MetaGap/manage.py collectstatic

# Compile translation message catalogs
cd MetaGap && DJANGO_SETTINGS_MODULE=MetaGap.settings django-admin compilemessages
```

## Architecture

### Directory Layout

```
/home/szuev/MetaGap/              # Repository root
├── CLAUDE.md                     # This file - project guidance
├── README.md
├── docs/                         # AI agent documentation and planning
│   └── plans/
└── MetaGap/                      # Django project root (working directory for manage.py)
    ├── manage.py
    ├── requirements.txt
    ├── db.sqlite3                # Development database
    ├── MetaGap/                  # Django settings package
    │   ├── settings.py           # Configuration (env-var driven)
    │   ├── urls.py               # Root URL routing (i18n_patterns wrapping app.urls)
    │   ├── context_processors.py # SITE_NAME/SITE_COLORS branding
    │   └── locale/               # Project-level translations (en/ru)
    ├── app/                      # Main Django application
    │   ├── models.py             # All data models (see below)
    │   ├── views.py              # All views (CBVs)
    │   ├── forms.py              # BootstrapFormMixin + all forms
    │   ├── filters.py            # django-filter definitions for search
    │   ├── tables.py             # django-tables2 dynamic table generation
    │   ├── signals.py            # Auto-creates OrganizationProfile on User creation
    │   ├── services/             # VCF import pipeline
    │   │   ├── vcf_importer.py   # Orchestrator (entry point: VCFImporter.import_file)
    │   │   ├── vcf_metadata.py   # YAML config loader + metadata parser
    │   │   ├── vcf_database.py   # Database persistence (VCFDatabaseWriter)
    │   │   ├── vcf_file_utils.py # Fallback VCF text parser (when pysam fails)
    │   │   └── import_exceptions.py # ImporterConfigurationError, ImporterValidationError
    │   ├── config/
    │   │   └── metadata_fields.yaml # Metadata section maps, model bindings, field aliases
    │   ├── locale/               # App-level translations (en/ru)
    │   ├── templates/            # Django templates
    │   ├── static/               # CSS, JS, images
    │   ├── tests/                # 20+ test modules, 320+ tests
    │   └── fixtures/demo_data.json
    └── MetagapUserCode/          # **IMPORTANT:** All project-related code that is NOT part of the Django web app belongs here
        ├── merge_vcf/            # Standalone VCF merge CLI tool
        │   ├── cli.py            # Entry point: python -m MetaGap.MetagapUserCode.merge_vcf.cli
        │   └── ...
        ├── demo_vcfs/            # Sample VCF files for testing
        ├── tests/                # Tests for standalone tools
        └── ...                   # Other utilities (shell scripts, test data, etc.)
```

**Structure Rules:**
- Django web app code: `/home/szuev/MetaGap/MetaGap/{app,MetaGap}/`
- Non-Django utilities/tools: `/home/szuev/MetaGap/MetaGap/MetagapUserCode/` ONLY
- No loose files in repository root except documentation (CLAUDE.md, README.md)
- Test VCF files belong in `MetagapUserCode/demo_vcfs/`, not in root

### Data Model (app/models.py)

**SampleGroup** is the central hub. It links to:
- **Organization**: `OrganizationProfile` (OneToOne with User, auto-created via signal)
- **Metadata models** (all FK with SET_NULL): `ReferenceGenomeBuild`, `GenomeComplexity`, `SampleOrigin`, `MaterialType`, `LibraryConstruction`, `InputQuality`
- **Sequencing platform models** (abstract base `SequencingInstrument`): `IlluminaSeq`, `OntSeq`, `PacBioSeq`, `IonTorrentSeq`
- **Bioinformatics models**: `BioinfoAlignment`, `BioinfoVariantCalling`, `BioinfoPostProc`
- **Variant data** (reverse FK): `AlleleFrequency` → `Info` + `Format` (CASCADE)

`SampleGroup.delete()` cascades but selectively protects shared metadata instances. AlleleFrequency has a unique constraint on (sample_group, chrom, pos, ref, alt).

### VCF Import Pipeline

`VCFImporter.import_file()` → parses VCF headers via `VCFMetadataParser` (driven by `metadata_fields.yaml`) → creates model instances via `VCFDatabaseWriter` → stores variants as `AlleleFrequency` records. Falls back to manual text parsing if pysam fails.

The YAML config (`app/config/metadata_fields.yaml`) has four mappings:
- `metadata_section_map` — VCF header sections → snake_case identifiers
- `metadata_models` — identifiers → Django model dotted paths
- `metadata_field_aliases` — canonical field names + legacy aliases
- `section_primary_field` — fallback field per section

Unmapped VCF header sections are stored in `SampleGroup.additional_metadata`.

### Adding/Modifying Metadata Fields

1. Update `metadata_fields.yaml` with new section/field/aliases
2. Add/modify Django model field + create migration
3. Run `python manage.py test app.tests.test_import_helpers.MetadataConfigurationTests` to validate config
4. Update `tables.py` friendly_names if UI display needed
5. Add test coverage in `test_import_helpers.py` and `test_vcf_importer.py`

### Key Patterns

- **BootstrapFormMixin** (`forms.py`): auto-applies Bootstrap CSS classes and data attributes to form widgets
- **Dynamic tables** (`tables.py`): `create_dynamic_table()` introspects model fields to build django-tables2 tables at runtime
- **Search** (`filters.py`): `AlleleFrequencySearchFilter` supports range queries on INFO fields and related model filters
- **I18n**: English/Russian, `.po` catalogs auto-compiled on app startup via `AppConfig.ready()`
- **Access control**: `OrganizationSampleGroupMixin` enforces ownership; `LoginRequiredMixin` for authenticated views

## CI/CD

GitHub Actions (`.github/workflows/ci.yml`) runs on push to `main` and PRs:
- Python 3.11, runs `coverage run manage.py test`
- Coverage thresholds: `app/` ≥ 70%, `MetagapUserCode/` ≥ 75%
- Scans test output for "Unhandled metadata key" warnings → CI annotations

## Environment Variables

Key settings configured via env vars (see `MetaGap/MetaGap/settings.py`):
- `DEBUG` (default: `1`), `SECRET_KEY`, `DATABASE_URL` (default: sqlite)
- `ALLOWED_HOSTS`, `CSRF_TRUSTED_ORIGINS` (comma-separated)
- `LOG_LEVEL` (default: `INFO`), `SITE_NAME`, `SITE_COLOR_PRIMARY`, `SITE_COLOR_SECONDARY`
- `METAGAP_MAX_UPLOAD_SIZE_MB`
