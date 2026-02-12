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
