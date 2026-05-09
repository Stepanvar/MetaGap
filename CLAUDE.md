# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MetaGap is a Django 5.1 web application for managing sample group metadata and allele frequency information from sequencing workflows. It provides VCF file import/export, metadata management, and search capabilities for genomic variant data. Python 3.11.

## Commands

All commands run from repo root `~/MetaGap/` (commands below prefix `MetaGap/` for the Django project dir):

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

Full detail in `.claude/docs/` (architecture.md, models.md, workflows.md). Quick map:

- `MetaGap/MetaGap/` — Django settings package (`settings.py`, `urls.py`, `locale/`)
- `MetaGap/app/` — main Django app (models, views, forms, filters, tables, services/, tests/)
- `MetaGap/app/services/` — VCF import pipeline (entry: `VCFImporter.import_file`)
- `MetaGap/app/config/metadata_fields.yaml` — YAML driving VCF → model field mapping
- `MetaGap/MetagapUserCode/` — **all non-Django code here** (CLI tools, demo VCFs, standalone tests)

**Rule:** Django web code → `{app,MetaGap}/`. Non-Django → `MetagapUserCode/` ONLY. No loose files at repo root.

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
