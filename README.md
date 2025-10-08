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
