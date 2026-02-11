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
