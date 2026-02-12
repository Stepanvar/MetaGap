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
/home/szuev/MetaGap/              # Repository root
├── CLAUDE.md                     # Project guidance for AI agents
├── README.md                     # User documentation
├── docs/                         # AI agent docs and planning
│   └── plans/
├── .claude/                      # Claude Code configuration
│   └── docs/                     # Detailed architecture references
│       ├── project-overview.md   # This file
│       ├── models.md             # Data model documentation
│       ├── workflows.md          # User workflows
│       └── architecture.md       # System architecture
└── MetaGap/                      # Django project root (working directory)
    ├── manage.py
    ├── requirements.txt
    ├── db.sqlite3                # Development database
    ├── MetaGap/                  # Settings package
    │   ├── settings.py           # All configuration
    │   ├── urls.py               # Root routing (i18n_patterns wrapper)
    │   ├── context_processors.py # SITE_NAME/SITE_COLORS
    │   └── locale/               # Project-level translations (en/ru)
    ├── app/                      # Main Django application
    │   ├── models.py             # 13 data models
    │   ├── views.py              # All views (CBVs)
    │   ├── forms.py              # BootstrapFormMixin + form classes
    │   ├── filters.py            # django-filter search configs
    │   ├── tables.py             # django-tables2 table builders
    │   ├── urls.py               # App URL patterns
    │   ├── signals.py            # Auto-create OrganizationProfile
    │   ├── services/             # VCF import pipeline
    │   │   ├── vcf_importer.py   # Orchestrator
    │   │   ├── vcf_metadata.py   # YAML config loader + parser
    │   │   ├── vcf_database.py   # Database writer
    │   │   ├── vcf_file_utils.py # Text fallback parser
    │   │   └── import_exceptions.py # Custom exceptions
    │   ├── config/
    │   │   └── metadata_fields.yaml # Metadata mapping config
    │   ├── locale/               # App-level translations (en/ru)
    │   ├── templates/            # Django templates
    │   ├── static/               # CSS, JS, images
    │   ├── tests/                # 20+ test modules, 320+ tests
    │   └── fixtures/demo_data.json # Demo users + sample groups
    └── MetagapUserCode/          # **IMPORTANT:** All non-Django project code
        ├── merge_vcf/            # Standalone VCF merge CLI tool
        │   ├── cli.py            # Entry point
        │   ├── __init__.py
        │   └── __main__.py
        ├── demo_vcfs/            # Sample VCF files for testing
        ├── tests/                # Tests for standalone tools
        └── ...                   # Other utilities (shell scripts, test data)
```

**Structure Rules:**
- **Django web app code:** `/home/szuev/MetaGap/MetaGap/{app,MetaGap}/` ONLY
- **Non-Django utilities/tools:** `/home/szuev/MetaGap/MetaGap/MetagapUserCode/` ONLY
- **Documentation:** `/home/szuev/MetaGap/{CLAUDE.md,README.md,docs/,.claude/}`
- **No loose files in repository root** except documentation
- **Test VCF files:** Belong in `MetagapUserCode/demo_vcfs/`, NOT in repository root

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
