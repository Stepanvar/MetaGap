# QA Dataset Manifest for Import/Search/Export Checks

Use this manifest with `docs/current-website-check-tasks.md` sections 6-9. It gives testers reproducible files, expected sample-group names, metadata checks, and row-level expectations without requiring knowledge of the importer internals.

## General tester setup

1. Start from a clean test database or a dedicated QA account so imported sample-group names do not collide with earlier runs.
2. Log in as a non-staff tester with a completed organization profile.
3. Open `/import/` and upload one dataset at a time.
4. After each successful import, open `/profile/`, then the imported sample-group detail page, then `/search/`, then `/profile/sample-groups/<id>/export/`.
5. For search checks, use exact row fields from the tables below: variant ID, chromosome, position, sample-group name, and numeric filters such as `DP ≥`, `AC ≥`, or `AN ≥` where available.
6. For export checks, download CSV and TSV if both are exposed. Confirm at minimum these columns are present when data exists: `chrom`, `pos`, `ref`, `alt`, `variant_id`, `info_dp`, and any dataset-specific INFO columns listed below.
7. For compressed formats generated during QA, keep the generated output under `/tmp/metagap-qa/` and delete it after the run.

## Required metadata baseline

For this QA checklist, a dataset is considered complete when the imported sample-group detail page exposes enough metadata for a tester to identify the source and VCF structure:

- sample group name;
- reference genome build or `##reference` value;
- contig metadata or chromosome naming convention;
- INFO field definitions needed by search/export checks, especially `DP`, `AC`, `AN`, and `AF` when present;
- FORMAT field definitions and sample identifier when genotype columns are present.

If the current application accepts a file that is missing this baseline, record the result as a product/validation issue instead of assuming the tester made a mistake.

## Positive fixtures

### P1: Small compressed single-sample VCF

| Field | Value |
| --- | --- |
| Path/source | `MetaGap/MetagapUserCode/demo_vcfs/1.vcf.gz` |
| Format | `.vcf.gz` |
| Size | 439 bytes |
| Expected sample-group name | `1.vcf` (fallback from the uploaded filename) |
| Expected sample identifier | `SAMPLE_001` |
| Expected key metadata fields | `fileformat=VCFv4.2`; `reference=GRCh37`; `contig=20`; INFO definitions `NS`, `DP`; FORMAT definitions `GT`, `DP`, `GQ`; sample-level FORMAT payload includes `GT`, `DP`, `GQ`, and `sample_id=SAMPLE_001`. |

Expected allele-frequency rows for search/export:

| Search cue | chrom | pos | variant_id | ref | alt | filter | info_dp | expected genotype | expected AF note |
| --- | --- | ---: | --- | --- | --- | --- | ---: | --- | --- |
| `rs6054257` or `14370` | `20` | 14370 | `rs6054257` | `G` | `C` | `PASS` | 13 | `0/1` | No `info_af` field in this file; manual genotype AF is 0.5. |
| `rs6054258` or `14420` | `20` | 14420 | `rs6054258` | `T` | `G` | `PASS` | 16 | `1/1` | No `info_af` field in this file; manual genotype AF is 1.0. |
| `rs6054259` or `14470` | `20` | 14470 | `rs6054259` | `A` | `T` | `PASS` | 19 | `./.` | No-call genotype; use row presence, not AF, as the pass criterion. |

### P2: Second small compressed VCF for authorization/data isolation

| Field | Value |
| --- | --- |
| Path/source | `MetaGap/MetagapUserCode/demo_vcfs/2.vcf.gz` |
| Format | `.vcf.gz` |
| Size | 436 bytes |
| Expected sample-group name | `2.vcf` (fallback from the uploaded filename) |
| Expected sample identifier | `SAMPLE_002` |
| Expected key metadata fields | `fileformat=VCFv4.2`; `reference=GRCh37`; `contig=20`; INFO definitions `NS`, `DP`; FORMAT definitions `GT`, `DP`, `GQ`; sample-level FORMAT payload includes `GT`, `DP`, `GQ`, and `sample_id=SAMPLE_002`. |

Expected allele-frequency rows for search/export:

| Search cue | chrom | pos | variant_id | ref | alt | filter | info_dp | expected genotype | expected AF note |
| --- | --- | ---: | --- | --- | --- | --- | ---: | --- | --- |
| `rs6054257` or `14370` | `20` | 14370 | `rs6054257` | `G` | `T` | `PASS` | 18 | `1/1` | No `info_af` field in this file; manual genotype AF is 1.0. |
| `rs6054258` or `14420` | `20` | 14420 | `rs6054258` | `T` | `A` | `PASS` | 9 | `./.` | No-call genotype; use row presence, not AF, as the pass criterion. |
| `rs6054259` or `14470` | `20` | 14470 | `rs6054259` | `A` | `C` | `PASS` | 12 | `0/0` | No `info_af` field in this file; manual genotype AF is 0.0. |

### P3: Site-only 1000 Genomes Y-chromosome VCF

| Field | Value |
| --- | --- |
| Path/source | `MetaGap/MetagapUserCode/trio.2010_06.ychr.sites.vcf.gz` |
| Format | `.vcf.gz` |
| Size | 12,937 bytes |
| Expected sample-group name | `trio.2010_06.ychr.sites.vcf` (fallback from the uploaded filename) |
| Expected sample identifier | None; this is a site-only VCF with no genotype sample columns. |
| Expected key metadata fields | `fileformat=VCFv4.0`; `fileDate=20100610`; `source=glfTools v3`; `reference=1000GenomesPilot-NCBI36`; `phasing=NA`; INFO definitions `NS`, `DP`, `DB`, `H2`, `AC`, `AN`; FILTER definition `NUYR`; FORMAT definitions `GT`, `GQ`, `DP` in header only. |

Expected allele-frequency rows for search/export:

| Search cue | chrom | pos | variant_id | ref | alt | filter | info_ac | info_an | info_dp | expected AF note |
| --- | --- | ---: | --- | --- | --- | --- | ---: | ---: | ---: | --- |
| `rs2058276` or `2728456` | `Y` | 2728456 | `rs2058276` | `T` | `C` | blank/`.` | 2 | 2 | 182 | Computed QA AF is 1.0; `info_af` is blank unless the app derives it. |
| `2734240` | `Y` | 2734240 | blank | `G` | `A` | blank/`.` | 1 | 2 | 196 | Computed QA AF is 0.5; `info_af` is blank unless the app derives it. |
| `2743242` | `Y` | 2743242 | blank | `C` | `T` | blank/`.` | 1 | 2 | 275 | Computed QA AF is 0.5; `info_af` is blank unless the app derives it. |
| `2746727` | `Y` | 2746727 | blank | `A` | `G` | blank/`.` | 2 | 2 | 179 | Computed QA AF is 1.0; `info_af` is blank unless the app derives it. |
| `rs2075640` or `2782506` | `Y` | 2782506 | `rs2075640` | `A` | `G` | blank/`.` | 1 | 2 | 254 | Computed QA AF is 0.5; `info_af` is blank unless the app derives it. |

### P4: DeepVariant WES VCF for larger import/performance smoke checks

| Field | Value |
| --- | --- |
| Path/source | `MetaGap/MetagapUserCode/ZD8U010E-GED-E080TSO-1715980255-BWA-DPV-hs37d5.filter.WES.vcf.gz` |
| Format | `.vcf.gz` |
| Size | 174,206 bytes |
| Expected sample-group name | `ZD8U010E-GED-E080TSO-1715980255-BWA-DPV-hs37d5.filter.WES.vcf` (fallback from the uploaded filename) |
| Expected sample identifier | `ZD8U010E` |
| Expected key metadata fields | `fileformat=VCFv4.2`; `DeepVariant_version=1.6.0`; contigs `chr1` through `chr22`, `chrX`, `chrY`, `chrMT`; FILTER definitions `PASS`, `RefCall`, `LowQual`, `NoCall`; FORMAT definitions `GT`, `GQ`, `DP`, `MIN_DP`, `AD`, `VAF`, `PL`, `MED_DP`; sample-level FORMAT payload includes `GT`, `GQ`, `DP`, `AD`, `VAF`, `PL`, and `sample_id=ZD8U010E`. |

Expected allele-frequency rows for search/export:

| Search cue | chrom | pos | variant_id | ref | alt | filter | expected genotype | expected FORMAT DP | expected FORMAT VAF | expected AF note |
| --- | --- | ---: | --- | --- | --- | --- | --- | ---: | ---: | --- |
| `977330` | `chr1` | 977330 | blank | `T` | `C` | `PASS` | `1/1` | 169 | 0.988166 | VAF is stored in FORMAT payload; export may not include it unless FORMAT export is implemented. |
| `981931` | `chr1` | 981931 | blank | `A` | `G` | `PASS` | `1/1` | 72 | 0.986111 | VAF is stored in FORMAT payload; export may not include it unless FORMAT export is implemented. |
| `982994` | `chr1` | 982994 | blank | `T` | `C` | `PASS` | `1/1` | 143 | 1.0 | VAF is stored in FORMAT payload; export may not include it unless FORMAT export is implemented. |
| `984302` | `chr1` | 984302 | blank | `T` | `C` | `PASS` | `1/1` | 58 | 1.0 | VAF is stored in FORMAT payload; export may not include it unless FORMAT export is implemented. |
| `987200` | `chr1` | 987200 | blank | `C` | `T` | `PASS` | `1/1` | 114 | 0.991228 | VAF is stored in FORMAT payload; export may not include it unless FORMAT export is implemented. |

## Generated positive-format variants

Use these only when the UI specifically needs `.vcf`, `.vcf.bgz`, or `.bcf` coverage and the repository does not already contain that exact format.

| Generated fixture | Source | Format | Generation command | Expected size | Expected sample group | Expected rows |
| --- | --- | --- | --- | --- | --- | --- |
| `/tmp/metagap-qa/1.vcf` | `MetaGap/MetagapUserCode/demo_vcfs/1.vcf.gz` | `.vcf` | `mkdir -p /tmp/metagap-qa && gzip -cd MetaGap/MetagapUserCode/demo_vcfs/1.vcf.gz > /tmp/metagap-qa/1.vcf` | Check with `wc -c /tmp/metagap-qa/1.vcf` after generation. | `1` if uploaded as `1.vcf`; rows match P1. | Rows match P1. |
| `/tmp/metagap-qa/1.vcf.bgz` | `/tmp/metagap-qa/1.vcf` | `.vcf.bgz` | Preferred: `bgzip -c /tmp/metagap-qa/1.vcf > /tmp/metagap-qa/1.vcf.bgz`; fallback extension-only smoke check: `cp MetaGap/MetagapUserCode/demo_vcfs/1.vcf.gz /tmp/metagap-qa/1.vcf.bgz`. | Check with `wc -c /tmp/metagap-qa/1.vcf.bgz` after generation. | `1.vcf` if uploaded as `1.vcf.bgz`; rows match P1. | Rows match P1. |
| `/tmp/metagap-qa/1.bcf` | `/tmp/metagap-qa/1.vcf` | `.bcf` | `bcftools view -Ob -o /tmp/metagap-qa/1.bcf /tmp/metagap-qa/1.vcf` | Check with `wc -c /tmp/metagap-qa/1.bcf` after generation. | `1` if uploaded as `1.bcf`; rows match P1 if BCF import is working. | Rows match P1. |

If `bgzip` or `bcftools` is not installed in the QA environment, mark only the generated-format check as blocked and continue with repository fixtures P1-P4.

## Negative fixtures

| ID | Path/source | Format | Size | How to use | Expected behavior |
| --- | --- | --- | ---: | --- | --- |
| N1 | `docs/qa-fixtures/unsupported-extension.txt` | unsupported `.txt` | 104 bytes | Upload from `/import/`. | The form rejects before import and displays an unsupported-file-type validation message. No sample group or allele rows are created. |
| N2 | `docs/qa-fixtures/malformed.vcf` | `.vcf` | 115 bytes | Upload from `/import/`. | The import fails safely because `POS=ABC` is not an integer. The user sees a safe validation error; logs contain the failure without a traceback escaping to the browser. |
| N3 | Generated oversized VCF | `.vcf` | Larger than configured `METAGAP_MAX_UPLOAD_SIZE_MB` | Generate with the command below, then upload. | The form rejects by size limit before parsing. No sample group or allele rows are created. |
| N4 | `docs/qa-fixtures/missing-required-metadata.vcf` | `.vcf` | 108 bytes | Upload from `/import/`. | The expected QA baseline is rejection or a clear warning because reference/contig/sample-group metadata is missing. If the app imports it using fallback metadata, record that as a validation gap. |

Oversized fixture generation command:

```bash
mkdir -p /tmp/metagap-qa
python - <<'PY'
from pathlib import Path
limit_mb = 51  # Use one MB above the default 50 MB limit; adjust if METAGAP_MAX_UPLOAD_SIZE_MB differs.
out = Path('/tmp/metagap-qa/oversized.vcf')
header = '##fileformat=VCFv4.2\n##reference=GRCh38\n##contig=<ID=1,length=100000000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'
row = '1\t1\trsOversized\tA\tC\t50\tPASS\tAC=1;AN=2;DP=10\n'
out.write_text(header)
with out.open('a') as handle:
    while out.stat().st_size <= limit_mb * 1024 * 1024:
        handle.write(row)
print(out, out.stat().st_size)
PY
```

## Reproducible import/search/export checklist

For each positive fixture:

1. Upload the file at `/import/`.
2. Confirm redirect to `/profile/` and a success message.
3. Confirm a sample group with the expected name appears under the current tester's organization.
4. Open the sample group detail page and verify expected metadata fields from this manifest.
5. Open `/search/` and run these checks:
   - `query=<variant_id>` for rows with non-blank IDs;
   - `chrom=<chrom>`, `pos_min=<pos>`, and `pos_max=<pos>` for rows with blank IDs;
   - `sample_group_name=<expected sample-group name>`;
   - `dp_min=<info_dp>` when the expected row has INFO DP;
   - `ac_min=<info_ac>` and `an_min=<info_an>` when the expected row has AC/AN.
6. Confirm each search result includes only rows that match the current tester's expected visibility policy.
7. Export CSV from `/profile/sample-groups/<id>/export/` and TSV from `/profile/sample-groups/<id>/export/?format=tsv` if available.
8. Compare exported rows against the tables above. Blank values are acceptable where the source file lacks that field; unexpected values or missing rows are failures.
9. Delete the imported sample group or reset the database before repeating the same fixture under the same account.

For authorization/data-isolation checks:

1. User A imports P1.
2. User B imports P2.
3. User A must see P1 in profile and must not see P2 in profile.
4. User B must see P2 in profile and must not see P1 in profile.
5. Direct access to the other user's detail, edit, delete, and export URLs must be denied or return not found.
6. Public search visibility must match the documented product policy; if cross-organization rows are visible, record the policy and confirm export remains owner-only.
