# MetaGap

MetaGap is a Django application for managing sample group metadata and allele frequency information collected from sequencing workflows.

## Demo data

A curated fixture is included to help you explore the interface without entering data manually. It creates:

- Two demo organisation accounts (`demo_org` and `demo_lab`) with the shared password `metagap-demo`.
- Rich sequencing metadata for two sample groups ("Lung Tumor Cohort" and "Rare Disease Trio").
- Linked reference information covering library preparation, sequencing platforms, and representative allele frequency records.

### Loading the fixture

```bash
python manage.py migrate
python manage.py loaddata app/fixtures/demo_data.json
```

After loading the data you can sign in with either demo account and visit the profile page to view the populated sample group table or run a search to explore the aggregated metadata.

## Development notes

- Install dependencies with `pip install -r requirements.txt`.
- Run checks with `python manage.py check` and execute the test suite using `python manage.py test`.
- Static files are served from the `static/` directory; new assets should be collected with `python manage.py collectstatic` in production environments.
- Set `LOG_LEVEL` (default: `INFO`) to tune runtime verbosity.  Even when `DEBUG=0`
  the configured logging ensures full tracebacks for request errors reach the
  server logs so developers can diagnose issues without exposing Django's debug
  pages to end users.

## Metadata configuration workflow

MetaGap's VCF importer is driven by a configuration module at
`MetaGap/app/services/vcf_metadata.py`.  It maps VCF header sections to Django
models (`METADATA_MODEL_MAP`), enumerates supported aliases
(`METADATA_FIELD_ALIASES`), defines section-level fallbacks, and normalizes any
input keys to snake_case before they are persisted.【F:MetaGap/app/services/vcf_metadata.py†L30-L200】【F:MetaGap/app/services/vcf_metadata.py†L320-L397】

### Alias mapping guidelines

- Prefer descriptive canonical field names and list common legacy spellings as
  aliases; the importer strips punctuation and lower-cases everything before
  comparison, so aliases should focus on semantic variants rather than casing
  changes.【F:MetaGap/app/services/vcf_metadata.py†L320-L368】【F:MetaGap/app/services/vcf_metadata.py†L384-L397】
- When retiring or renaming a field, keep the previous spelling in the alias
  list until all upstream feeds are migrated so older VCFs continue to ingest
  cleanly.【F:MetaGap/app/services/vcf_metadata.py†L78-L190】

### Adding or renaming metadata fields

1. Extend the relevant section in `vcf_metadata.py` with the new field and any
   aliases, and update `SECTION_PRIMARY_FIELD` if the field should become the
   section's display default.【F:MetaGap/app/services/vcf_metadata.py†L78-L200】
2. Ensure the corresponding Django model exposes the new attribute and, when
   necessary, add a migration so it is stored alongside imported data.【F:MetaGap/app/models.py†L360-L458】
3. Review the importer/database integration tests to keep expectations aligned:
   `test_import_helpers.py` exercises alias resolution while
   `test_vcf_importer.py` verifies the end-to-end import payload.  Update or add
   fixtures so the new field is covered by these guard rails.【F:MetaGap/app/tests/test_import_helpers.py†L134-L210】【F:MetaGap/app/tests/test_vcf_importer.py†L242-L320】
4. If the field should appear prominently in search results or detail tables,
   update the `friendly_names` mapping and curated ordering in
   `MetaGap/app/tables.py` so the UI renders a readable label and groups the
   column with related metadata.【F:MetaGap/app/tables.py†L160-L260】

### Removing metadata fields

- Remove the field and its aliases from `vcf_metadata.py`, drop or adjust the
  associated model field and migration, and prune any tests that depended on the
  retired metadata.  When backwards compatibility is required, leave an alias in
  place and route the value into `SampleGroup.additional_metadata` instead of
  deleting it outright.【F:MetaGap/app/services/vcf_metadata.py†L78-L200】【F:MetaGap/app/tests/test_vcf_importer.py†L301-L320】

### Configuration loading and validation

`VCFImporter` wires the configuration into the runtime importer by instantiating
`VCFMetadataParser` and `VCFDatabaseWriter` for each import.  This loader applies
the dictionaries above before any rows are written, ensuring updates to the
configuration module immediately affect new uploads.【F:MetaGap/app/services/vcf_importer.py†L31-L80】【F:MetaGap/app/services/vcf_metadata.py†L30-L200】【F:MetaGap/app/services/vcf_database.py†L1-L104】

After making configuration changes, run `python manage.py test` so the importer,
table rendering, and integration suites confirm the workflow remains consistent
with the documented expectations.【F:MetaGap/app/tests/test_import_helpers.py†L134-L210】【F:MetaGap/app/tests/test_tables.py†L109-L202】
### Maintaining VCF metadata configuration

Metadata section mappings, model bindings, and field aliases are now stored in
`MetaGap/app/config/metadata_fields.yaml`.  The loader in
`app/services/vcf_metadata.py` validates the file, caches the parsed result, and
exposes it to the importer/database helpers.  The YAML file contains four
top-level mappings:

* `metadata_section_map` – normalises header section names (e.g. `SAMPLE_GROUP`
  → `sample_group`).
* `metadata_models` – maps section identifiers to the Django model class that
  persists the parsed data.  Each value must be a dotted import path such as
  `app.models.IlluminaSeq`.
* `metadata_field_aliases` – lists the recognised aliases for each model field
  within a section.  Provide every alternate header key that should populate a
  field.  Empty lists are allowed when no aliases exist.
* `section_primary_field` – designates the field that should mirror the overall
  section value when a header omits structured keys.

To add a new metadata field or section:

1. Update `metadata_fields.yaml` with the new section key, model mapping, and
   alias list.  Keep keys upper-cased in `metadata_section_map` to match VCF
   header conventions and ensure the model path is importable.
2. Run `python manage.py test app.tests.test_import_helpers.MetadataConfigurationTests`
   to confirm the configuration parses successfully.  The tests also surface
   meaningful errors if the YAML is malformed or the model reference cannot be
   imported.
3. If a new Django model or field is required, add it to `app/models.py` and
   create a migration before adjusting the configuration.

The importer and database writer use the cached configuration automatically, so
no further code changes are needed once the YAML file and models are in sync.

### Logging format for the VCF merger

The `MetaGap.MetagapUserCode.merge_vcf` workflow logs to both `script_execution.log`
and the console using the format:

```
YYYY-MM-DD HH:MM:SS | LEVEL | vcf_merger | module | message
```

The shared formatter ensures severity and module details are visible across file and
console handlers so operational logs can be correlated easily.

### Continuous integration expectations

- Every pull request and push to `main` runs the Django test suite (`python manage.py test`) under coverage on GitHub Actions.
- The workflow archives the raw test log and coverage summary so failures can be triaged without re-running the suite locally.
- Coverage is enforced for critical apps:
  - `MetaGap/app`: minimum 80% line coverage.
  - `MetaGap/MetagapUserCode`: minimum 75% line coverage.
- Test output is scanned for `Unhandled metadata key` warnings; when present they surface as CI annotations so developers can address missing metadata handlers before merging.
- New features should maintain or improve coverage. If code changes introduce new modules, add tests to keep the coverage thresholds passing.

## Merging VCF files

The production CLI for the MetaGap VCF merging workflow lives in
`MetaGap.MetagapUserCode.merge_vcf.cli:main`. Invoke it either directly
(`python -m MetaGap.MetagapUserCode.merge_vcf.cli`) or via the
compatibility shim that ships with the package (`python -m merge_vcf`).
Both routes execute the same tested workflow that powers
`MetagapUserCode.test_merge_vcf`.

### Logging configuration

Shared helpers under ``MetaGap.MetagapUserCode.merge_vcf.logging_utils``
expose a ``configure_logging`` function so downstream tooling can reuse the
workflow's log formatting while directing output to a custom destination.
For example, to record only warnings to a dedicated file:

```python
from MetaGap.MetagapUserCode.merge_vcf.logging_utils import configure_logging

configure_logging(log_level="WARNING", log_file="/tmp/merge.log")
```

The helper clears existing handlers before applying the new configuration, so
it can be invoked multiple times without creating duplicate console or file
output streams.
