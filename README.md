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
