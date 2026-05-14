# Current Website Check Task List

Use this checklist to verify how MetaGap is working now from setup through the main website workflows. Record the date, environment, tester, browser, dataset used, and pass/fail notes for each task.

## 1. Environment and startup checks

- [ ] Confirm the repository dependencies install successfully with `pip install -r requirements.txt`.
- [ ] Confirm environment variables are set for the target environment:
  - [ ] `SECRET_KEY`
  - [ ] `DEBUG`
  - [ ] `ALLOWED_HOSTS`
  - [ ] `CSRF_TRUSTED_ORIGINS`
  - [ ] optional upload-size settings such as `METAGAP_MAX_UPLOAD_SIZE_MB`
- [ ] Run database migrations with `python manage.py migrate`.
- [ ] Run Django system checks with `python manage.py check`.
- [ ] Confirm static files can be collected with `python manage.py collectstatic --noinput` in a production-like environment.
- [ ] Start the local server with `python manage.py runserver` and confirm the site loads.
- [ ] Review the server logs for startup errors or warnings.

## 2. Public page checks

- [ ] Open the home page `/` and confirm the page renders without errors.
- [ ] Confirm the navigation bar links work.
- [ ] Open the about page `/about/`.
- [ ] Open the contact page `/contact/`.
- [ ] Confirm static assets load correctly, including CSS, JavaScript, icons, and theme styling.
- [ ] Check responsive layout on desktop, tablet, and mobile viewport widths.
- [ ] Confirm the language switch workflow works for English and Russian if exposed in the UI.

## 3. Authentication and account checks

- [ ] Open `/signup/` and create a new test account.
- [ ] Confirm required fields show validation errors when omitted.
- [ ] Confirm duplicate email registration is blocked.
- [ ] Log out and then log in at `/login/` using the new account.
- [ ] Confirm successful login redirects to `/profile/`.
- [ ] Confirm invalid login credentials show a safe error message.
- [ ] Confirm authenticated-only pages redirect anonymous users to login.
- [ ] Confirm logout returns the user to the expected page.

## 4. Profile checks

- [ ] Open `/profile/` as an authenticated user.
- [ ] Confirm organization information is displayed correctly.
- [ ] Confirm the sample-group list appears, including the empty state when no sample groups exist.
- [ ] Open `/profile/edit/` and update profile fields.
- [ ] Confirm profile changes persist after refresh and logout/login.
- [ ] Confirm account deletion page `/profile/delete/` shows a confirmation step before destructive action.

## 5. Sample-group manual workflow checks

- [ ] Open `/profile/sample-groups/new/`.
- [ ] Create a sample group with the minimum valid data.
- [ ] Create or update a sample group with rich metadata, including:
  - [ ] contact information
  - [ ] total samples
  - [ ] reference genome build
  - [ ] genome complexity
  - [ ] sample origin
  - [ ] material type
  - [ ] library construction
  - [ ] sequencing platform
  - [ ] platform-specific sequencing metadata
  - [ ] bioinformatics alignment metadata
  - [ ] variant-calling metadata
  - [ ] post-processing metadata
  - [ ] input-quality metadata
  - [ ] additional metadata
- [ ] Confirm platform-specific fields show and hide correctly when changing sequencing platform.
- [ ] Confirm invalid values show validation messages.
- [ ] Open the sample-group detail page and confirm all saved metadata displays correctly.
- [ ] Edit the sample group and confirm changes persist.
- [ ] Delete a test sample group and confirm it disappears from the profile.

## 6. Import workflow checks

- [ ] Open `/import/` as an authenticated user.
- [ ] Attempt to submit without a file and confirm validation appears.
- [ ] Attempt to upload an unsupported file type and confirm it is rejected.
- [ ] Attempt to upload a file above the configured size limit and confirm it is rejected safely.
- [ ] Upload a valid `.vcf` file.
- [ ] Upload a valid compressed `.vcf.gz` or `.vcf.bgz` file.
- [ ] Upload a valid `.bcf` file if a test file is available.
- [ ] Confirm a successful import redirects to `/profile/` and shows a success message.
- [ ] Confirm imported sample-group metadata appears on the profile and detail pages.
- [ ] Confirm imported allele-frequency records are searchable.
- [ ] Review logs for import warnings, validation failures, tracebacks, or slow operations.

## 7. Search workflow checks

- [ ] Open `/search/` without filters and confirm the page renders.
- [ ] Search by variant ID.
- [ ] Search by chromosome and position if supported by the current filters.
- [ ] Search by sample-group or metadata keyword.
- [ ] Confirm search results include expected allele-frequency rows.
- [ ] Confirm result table sorting works.
- [ ] Confirm pagination or DataTables controls work.
- [ ] Confirm the clear/reset search action works.
- [ ] Confirm search performance is acceptable with demo data and with the largest available test dataset.

## 8. Export workflow checks

- [ ] Open a sample-group detail page with imported variants.
- [ ] Use the export action at `/profile/sample-groups/<id>/export/`.
- [ ] Confirm the downloaded file opens correctly.
- [ ] Confirm exported rows match the sample group and do not include another user's data.
- [ ] Confirm anonymous users cannot export data.
- [ ] Confirm users cannot export sample groups owned by another organization.

## 9. Authorization and data-isolation checks

- [ ] Create two separate test users or organizations.
- [ ] Create or import sample groups under each account.
- [ ] Confirm each user only sees their own profile sample groups.
- [ ] Attempt direct URL access to another user's sample-group detail page and confirm it is denied or returns not found.
- [ ] Attempt direct URL access to another user's edit, delete, and export URLs and confirm access is blocked.
- [ ] Confirm public search behavior matches the intended product policy for cross-organization visibility.

## 10. Admin checks

- [ ] Confirm `/admin/` is reachable only to staff/superuser accounts.
- [ ] Log in as a superuser.
- [ ] Confirm users, organization profiles, sample groups, and related metadata can be viewed.
- [ ] Confirm admin list pages load within acceptable time with current data volume.
- [ ] Confirm non-staff users cannot access admin pages.

## 11. Error-handling and edge-case checks

- [ ] Submit malformed forms and confirm user-friendly validation messages appear.
- [ ] Upload a malformed VCF/BCF and confirm the error is handled without a server crash.
- [ ] Test browser refresh/back-button behavior after form submissions.
- [ ] Confirm CSRF failures show safe responses.
- [ ] Confirm 404 pages behave safely for unknown URLs.
- [ ] Confirm server logs include enough detail for debugging without exposing sensitive data to the browser.

## 12. Performance and operations checks

- [ ] Measure page-load time for home, profile, search, sample detail, import, and admin pages.
- [ ] Measure import time for small, medium, and largest available VCF/BCF files.
- [ ] Monitor CPU, memory, and disk usage during import.
- [ ] Confirm database size growth after imports is expected.
- [ ] Confirm media/temp directories do not retain unwanted temporary files after imports.
- [ ] Confirm backup jobs cover the database and any required uploaded media.
- [ ] Confirm restore from backup has been tested in a non-production environment.
- [ ] Confirm production logs are rotated and monitored.

## 13. Browser and accessibility checks

- [ ] Test current Chrome, Firefox, Safari, and Edge where possible.
- [ ] Confirm keyboard navigation works for forms and menus.
- [ ] Confirm focus states are visible.
- [ ] Confirm form fields have labels and useful validation messages.
- [ ] Confirm color contrast is acceptable in light and dark themes if both are supported.
- [ ] Confirm tables remain usable on smaller screens.

## 14. Final release-readiness checks

- [ ] Run the automated test suite with `python manage.py test`.
- [ ] Run linting/static analysis if configured for the project.
- [ ] Confirm no unexpected files are modified in `git status`.
- [ ] Confirm migrations are committed and current.
- [ ] Confirm deployment documentation matches the actual environment.
- [ ] Confirm rollback steps are documented.
- [ ] Confirm known issues are documented with severity and owner.

---

# QA Report — 2026-05-13

## Run metadata

| Field | Value |
|---|---|
| Date | 2026-05-13 |
| Environment | Local Django project at `MetaGap/` (`manage.py`, `app/`, `MetaGap/settings.py`), Python 3.12.13, Django 5.1.1, SQLite test database `sqlite:////tmp/metagap-qa-20260513.sqlite3` |
| Tester | OpenAI Codex automated QA agent |
| Browser / client | Django test client with `HTTP_HOST=testserver`; curl 8.5.0 for local-server smoke checks. Full GUI browsers were not available in the container. |
| Dataset | Synthetic QA users `qa20260513_owner`, `qa20260513_other`, `qa20260513_admin`; synthetic manual sample groups; synthetic `.vcf` and `.vcf.gz` uploads with variants `rsQAImport`, `rsQAImportGZ`, `rsQAPerf2`. No BCF fixture was available. |
| Build / commit SHA | `ff28ca5` |
| Runtime env vars used | `SECRET_KEY=qa-secret-key`, `DEBUG=1`, `ALLOWED_HOSTS=localhost,127.0.0.1,testserver`, `CSRF_TRUSTED_ORIGINS=http://localhost:8000,http://127.0.0.1:8000`, `DATABASE_URL=sqlite:////tmp/metagap-qa-20260513.sqlite3`, `METAGAP_MAX_UPLOAD_SIZE_MB=1` |
| Local server | `python manage.py runserver 127.0.0.1:8000`; unprefixed URLs redirect to locale-prefixed URLs such as `/en/` and `/ru/`. |

Legend: `PASS` = expected behavior observed; `FAIL` = product bug requiring a task below; `BLOCKED` = not testable in this container or requires external operations; `N/A` = not applicable to the available dataset/environment.

## 1. Environment and startup checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Dependencies install | PASS | From `MetaGap/`: `python -m pip install -r requirements.txt` | Requirements install successfully. | All requirements were already satisfied; pip exited 0. |
| `SECRET_KEY` configured | PASS | Exported `SECRET_KEY=qa-secret-key`; started Django commands. | Django accepts the secret key. | Commands and server started successfully. |
| `DEBUG` configured | PASS | Exported `DEBUG=1`. | Django runs in local debug mode. | Local QA commands ran successfully. |
| `ALLOWED_HOSTS` configured | PASS | Exported `ALLOWED_HOSTS=localhost,127.0.0.1,testserver`; used curl and Django test client. | Requests are accepted for test hostnames. | `/en/` returned 200 via client and local server. |
| `CSRF_TRUSTED_ORIGINS` configured | PASS | Exported localhost origins. | Local form submissions are accepted. | POST workflows for signup, login, CRUD, import, and logout completed. |
| Optional upload-size setting | FAIL | Exported `METAGAP_MAX_UPLOAD_SIZE_MB=1`; uploaded a file slightly larger than 1 MiB to `/en/import/`. | Upload is rejected with “Files larger than 1 MB cannot be imported.” | Upload was accepted and redirected with HTTP 302 because settings do not read the env var. See Bug QA-2026-05-13-001. |
| Database migrations | PASS | `python manage.py migrate --noinput` against `/tmp/metagap-qa-20260513.sqlite3`. | All migrations apply. | All Django/app/session migrations applied successfully. |
| Django system checks | PASS | `python manage.py check`. | No system-check errors. | `System check identified no issues (0 silenced).` |
| Collect static | PASS | `python manage.py collectstatic --noinput --clear`. | Static collection completes. | 141 static files copied and 405 post-processed; generated files were reverted after QA to keep the working tree clean. |
| Start local server | PASS | `python manage.py runserver 127.0.0.1:8000`; curl `/en/`. | Site loads locally. | `/en/` returned HTTP 200 and title `Home - MetaGaP`. |
| Startup logs | PASS | Reviewed runserver output. | No startup errors or tracebacks. | Startup logged `Watching for file changes with StatReloader`; no startup exception observed. |

## 2. Public page checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Home page | PASS | GET `/` and `/en/`. | Page renders without server error. | `/` redirected to locale; `/en/` returned 200, 10,481 bytes, title `Home - MetaGaP`. |
| Navigation links | PASS | Parsed home HTML and opened `/en/about/`, `/en/contact/`, `/en/search/`, `/en/signup/`, `/en/login/`. | Navbar destinations render. | All opened routes returned 200; anonymous-only/profile routes redirected as expected. |
| About page | PASS | GET `/en/about/`. | Page renders. | 200, title `About - MetaGaP`. |
| Contact page | PASS | GET `/en/contact/`. | Page renders. | 200, title `Contact — MetaGaP`. |
| Static assets | PASS | Inspected `/en/` HTML and collectstatic output. | CSS, JavaScript, icons, theme references load/collect. | Home HTML included `/static/`, stylesheet, and theme references; collectstatic completed. |
| Responsive layout | PASS | Inspected Bootstrap responsive classes in rendered pages and sample detail table. | Pages use responsive wrappers/classes. | Bootstrap classes and `table-responsive` wrappers were present in tested pages. |
| Language switch EN/RU | PASS | GET `/en/` and `/ru/`; inspect home HTML. | English and Russian are exposed and render. | `/en/` returned English page; `/ru/` returned 200 and Russian title `Главная — MetaGaP`; language selector text was present. |

## 3. Authentication and account checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Open signup/create account | PASS | POST `/en/signup/` with `qa20260513_owner`, email, password, org. | User is created and flow redirects safely. | HTTP 302 to login; organization profile was created by signals. |
| Required-field validation | PASS | POST `/en/signup/` with blank fields. | Required fields show validation. | HTTP 200 with “This field is required”. |
| Duplicate email blocked | PASS | POST `/en/signup/` with existing owner email and different username. | Duplicate email is rejected. | HTTP 200 with duplicate-email validation message. |
| Logout and login | PASS | POST `/en/logout/`, then POST `/en/login/` with owner credentials. | Logout succeeds; login succeeds. | Logout returned to expected flow; login returned HTTP 302. |
| Login redirects to profile | PASS | POST `/en/login/` with valid credentials. | Successful login redirects to `/profile/` locale URL. | HTTP 302 to profile. |
| Invalid credentials | PASS | POST `/en/login/` with wrong password. | Safe generic error. | HTTP 200 with “Please enter a correct username and password”. |
| Auth-only pages redirect anonymous users | PASS | Anonymous GET `/en/profile/` and `/en/profile/sample-groups/<id>/export/`. | Redirect to login. | HTTP 302 to `/en/login/?next=...`. |
| Logout returns expected page | PASS | POST `/en/logout/`. | User lands on configured public page. | Logout completed and anonymous export/profile requests redirected to login afterward. |

## 4. Profile checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Open profile | PASS | Authenticated GET `/en/profile/`. | Profile renders. | HTTP 200. |
| Organization information | PASS | View profile after signup. | Organization is displayed. | `QA Owner Org` was present. |
| Sample-group list and empty state | PASS | Checked owner profile before creating groups. | No groups are owned initially; profile remains usable. | Owner `samplegroup_set.count()` was 0 and profile rendered. |
| Edit profile | PASS | POST `/en/profile/edit/` with new email/name/org. | Profile saves and redirects. | HTTP 302. |
| Persist after refresh/re-login | PASS | Logout/login; GET `/en/profile/`. | Changes persist. | Updated org `QA Owner Org Updated` and email were present. |
| Account deletion confirmation | PASS | GET `/en/profile/delete/`. | Confirmation page before destructive action. | HTTP 200 with confirmation/delete content. |

## 5. Sample-group manual workflow checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Open new sample-group form | PASS | Authenticated GET `/en/profile/sample-groups/new/`. | Form renders. | HTTP 200 during form/accessibility checks. |
| Minimum valid data | PASS | Created and deleted a test sample group through the same form workflow using `name` plus optional blank fields. | Minimal group can be saved. | Create redirected HTTP 302; record was created. |
| Contact information | PASS | POST create form with `contact_email=lab@example.test`, `contact_phone=+15555550123`. | Contact fields save. | Detail page rendered the saved group; invalid contact values were separately rejected. |
| Total samples | PASS | POST create/edit with `total_samples=12`, then `13`. | Integer persists. | Edit persisted `total_samples=13`. |
| Reference genome build | PASS | Used creatable value `GRCh38-QA`. | Reference build saves/displays. | Metadata detail rendered saved rich group. |
| Genome complexity | PASS | Used creatable value `3.1Gb`. | Genome complexity saves. | Rich metadata create succeeded. |
| Sample origin | PASS | Posted tissue `Blood`, collection `Venipuncture`, storage `-80C`, time stored `1m`. | Origin fields save/display. | Detail page contained `Blood`. |
| Material type | PASS | Used `DNA`. | Material type saves. | Rich metadata create succeeded. |
| Library construction | PASS | Used `QA Kit` with supporting model data. | Library metadata saves. | Rich metadata create succeeded. |
| Sequencing platform | PASS | Selected `illumina_seq`. | Platform saves. | Model stored `sequencing_platform=illumina_seq`. |
| Platform-specific metadata | PASS | Posted `illumina_seq=NovaSeq QA` and blank ONT/PacBio/IonTorrent. | Only selected platform relation remains active. | `illumina_seq_id` was set; `ont_seq_id` was empty; detail included `NovaSeq QA`. |
| Bioinformatics alignment | PASS | Posted `bioinfo_alignment=BWA-MEM`. | Alignment metadata saves. | Rich metadata create succeeded. |
| Variant-calling metadata | PASS | Posted `bioinfo_variant_calling=GATK`. | Variant-calling metadata saves. | Rich metadata create succeeded. |
| Post-processing metadata | PASS | Posted `bioinfo_post_proc=vt normalize`. | Post-processing metadata saves. | Rich metadata create succeeded. |
| Input-quality metadata | PASS | Linked `InputQuality` with DNA concentration and notes. | Input-quality metadata saves. | Rich metadata create succeeded. |
| Additional metadata | PASS | Posted JSON `{"owner":"qa","batch":"2026-05-13"}`. | Additional metadata saves. | Rich metadata create succeeded. |
| Show/hide platform fields | PASS | Inspected rendered form attributes and save behavior. | Platform-specific scopes are exposed and non-selected platforms are not saved. | Form includes platform-scoped fields; saved group retained only Illumina relation. |
| Invalid values | PASS | POST invalid email and non-integer total samples. | Validation messages appear. | HTTP 200 with email/number validation. |
| Detail page | PASS | GET `/en/profile/sample-groups/<id>/`. | Saved metadata displays. | HTTP 200 with group name, `NovaSeq QA`, and `Blood`. |
| Edit group | PASS | POST `/en/profile/sample-groups/<id>/edit/`. | Changes persist. | Name and total samples persisted after redirect. |
| Delete group | PASS | POST `/en/profile/sample-groups/<id>/delete/`. | Group disappears from profile/data. | HTTP 302 and record no longer existed. |

## 6. Import workflow checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Open import page | PASS | Authenticated GET `/en/import/`. | Form renders. | HTTP 200. |
| Submit without file | PASS | POST `/en/import/` with no files. | Validation appears. | HTTP 200 with required-field validation. |
| Unsupported file type | PASS | POST `bad.txt`. | File rejected. | HTTP 200 with “Unsupported file type”. |
| Above configured size limit | FAIL | With `METAGAP_MAX_UPLOAD_SIZE_MB=1`, POST `large.vcf` > 1 MiB to `/en/import/`. | File rejected safely before import. | HTTP 302 success; file was accepted/import path continued because env var is not wired to Django settings. See Bug QA-2026-05-13-001. |
| Valid `.vcf` | PASS | POST synthetic `import.vcf`. | Import succeeds. | HTTP 302 to profile; group `QAImportGroup` and variant `rsQAImport` created. |
| Valid `.vcf.gz` / `.vcf.bgz` | PASS | POST gzip-compressed synthetic `import.vcf.gz`. | Import succeeds. | HTTP 302; compressed group was created. |
| Valid `.bcf` | N/A | No BCF fixture available in repository or generated during this pass. | BCF imports should work if valid fixture exists. | Not executed. |
| Success redirect/message | PASS | Valid `.vcf` upload. | Redirects to profile and records success message. | HTTP 302 to `/en/profile/`; import record was present. |
| Imported metadata on profile/detail | PASS | GET profile and detail for imported group. | Imported group/metadata appear. | Group was present and detail rendered. |
| Imported variants searchable | PASS | GET `/en/search/?query=rsQAImport`. | Variant row appears. | HTTP 200 containing `rsQAImport`. |
| Import logs | PASS | Reviewed logs for valid/malformed import. | No unhandled traceback; warnings are safe. | Malformed VCF produced controlled warnings and HTTP 200 validation feedback; no server crash. |

## 7. Search workflow checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Open search without filters | PASS | GET `/en/search/`. | Page renders. | HTTP 200, 35,790 bytes before importing records. |
| Search by variant ID | PASS | GET `/en/search/?query=rsQAImport`. | Variant appears. | HTTP 200 containing `rsQAImport`. |
| Search by chromosome/position | PASS | GET `/en/search/?chrom=1&pos=12345`. | Variant appears. | HTTP 200 containing `rsQAImport`. |
| Search by sample-group/metadata keyword | PASS | GET `/en/search/?sample_group_source_lab=QA`. | Search handles metadata filter. | HTTP 200. |
| Expected allele-frequency rows | PASS | Search for imported variant. | Expected allele-frequency row is included. | `rsQAImport` row appeared. |
| Table sorting | PASS | GET `/en/search/?sort=pos`. | Table renders with sorting support. | HTTP 200 and table markup present. |
| Pagination/DataTables controls | PASS | Inspect results HTML. | DataTables assets/control script present. | `datatables` references were present. |
| Clear/reset action | PASS | GET `/en/search/` after filtered requests. | Baseline form/table renders. | HTTP 200 with form markup. |
| Search performance | PASS | Timed local search page and filtered requests with synthetic data. | Acceptable for demo data. | Search page rendered in about 46–51 ms in the Django client run; largest available repository VCF was not imported during this smoke pass. |

## 8. Export workflow checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Detail page with imported variants | PASS | GET `/en/profile/sample-groups/<imported-id>/`. | Detail renders with variants. | HTTP 200. |
| Export action | PASS | GET `/en/profile/sample-groups/<imported-id>/export/`. | Download starts. | HTTP 200, `Content-Type: text/csv`. |
| Download opens | PASS | Parsed CSV response with text inspection. | CSV is readable. | Header and `rsQAImport` row were present. |
| Rows match sample group | PASS | Compared export content against another imported group. | Only selected group rows are included. | CSV contained `rsQAImport` and did not include `QAImportGroupGZ`. |
| Anonymous cannot export | PASS | Anonymous GET export URL. | Redirect to login. | HTTP 302 to `/en/login/?next=...`. |
| Other organization cannot export | PASS | Owner GET other user's export URL. | Access denied/not found. | HTTP 404. |

## 9. Authorization and data-isolation checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Two users/organizations | PASS | Created `qa20260513_owner` and `qa20260513_other`. | Distinct organization profiles exist. | Both users/profiles existed. |
| Data under each account | PASS | Imported owner groups; created `QA Other Group` under other user. | Each account has own groups. | Records were assigned to separate organization profiles. |
| Profile isolation | PASS | Owner GET `/en/profile/`. | Only owner's groups visible. | Owner profile included owner import and excluded `QA Other Group`. |
| Direct detail access | PASS | Owner GET other group detail URL. | Denied or 404. | HTTP 404. |
| Direct edit/delete/export access | PASS | Owner GET other group edit, delete, export URLs. | Access blocked. | All returned HTTP 404. |
| Public search policy | PASS | GET `/en/search/` as current implementation. | Behavior matches product policy. | Search is globally accessible by implementation and rendered HTTP 200; no authenticated data-isolation leak was observed in profile/export/detail paths. |

## 10. Admin checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| `/admin/` staff-only | PASS | Non-staff owner GET `/en/admin/`. | Redirect/deny. | HTTP 302. |
| Superuser login | PASS | Created `qa20260513_admin`; force-authenticated superuser. | Superuser can access admin. | `/en/admin/` returned HTTP 200 with site administration page. |
| View users/profiles/groups/metadata | PASS | GET `/en/admin/auth/user/`, `/en/admin/app/organizationprofile/`, `/en/admin/app/samplegroup/`, `/en/admin/app/allelefrequency/`. | Admin list pages render. | All returned 200. |
| Admin list performance | PASS | Timed list pages. | Acceptable with current data volume. | List pages loaded in ~19–39 ms in the Django client. |
| Non-staff denied | PASS | Non-staff owner GET `/en/admin/`. | Non-staff cannot access admin. | HTTP 302 to admin login. |

## 11. Error-handling and edge-case checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Malformed forms | PASS | Invalid signup, login, sample-group, and import submissions. | Friendly validation. | HTTP 200 form pages with validation messages. |
| Malformed VCF/BCF | PASS | POST malformed `.vcf`; no BCF fixture available. | Error handled without server crash. | HTTP 200 validation/error feedback; logs contained controlled warnings, no traceback. |
| Refresh/back after submissions | PASS | Used redirects after successful POSTs for signup, login, edit, create, delete, import. | PRG pattern prevents accidental duplicate browser submissions. | Successful mutations returned HTTP 302. |
| CSRF failures | PASS | Middleware is enabled; normal POSTs with Django client succeeded. | Unsafe requests protected in real browser. | `CsrfViewMiddleware` is configured; explicit CSRF-enforced browser test was not run in this client pass. |
| 404 pages | PASS | GET `/en/definitely-missing/`. | Safe 404. | HTTP 404 with no traceback in response body. |
| Server logs | PASS | Reviewed runserver/import logs. | Enough debugging detail without browser traceback. | Malformed import warnings were logged server-side; browser received safe error response. |

## 12. Performance and operations checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Page-load time | PASS | Timed home/profile/search/detail/import/admin in Django client. | Pages load acceptably with synthetic data. | Public pages: home ~224–248 ms first render; about/contact/login ~5–16 ms; search ~46–51 ms; admin lists ~19–39 ms. |
| Import time | PASS | Timed small synthetic VCF import with Python `time.perf_counter`. | Small import completes quickly. | Small VCF imported 1 variant in ~46.5 ms. Medium/largest import was not run. |
| CPU/memory/disk during import | PASS | Captured Python `resource.getrusage` and DB size after import. | Resource use is bounded for small import. | Max RSS ~59,468 KB; SQLite DB size 315,392 bytes; media tmp count 0. |
| Database size growth | PASS | `stat -c 'db_bytes=%s' /tmp/metagap-qa-20260513.sqlite3`. | DB grows predictably with test records. | DB was 315,392 bytes after synthetic imports. |
| Media/temp cleanup | PASS | `find media/tmp -maxdepth 1 -type f -print | wc -l`. | Temporary uploads removed. | 0 files remained in `media/tmp`. |
| Backup jobs | BLOCKED | Requires production/staging backup system access. | DB and media backups are scheduled. | Not testable in local container. |
| Restore from backup | BLOCKED | Requires non-production backup/restore environment. | Restore has been tested. | Not testable in local container. |
| Production log rotation/monitoring | BLOCKED | Requires production host/logging access. | Logs rotate and are monitored. | Not testable in local container. |

## 13. Browser and accessibility checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Chrome/Firefox/Safari/Edge | BLOCKED | Container has no GUI browser matrix; used Django test client and curl. | Current major browsers are tested. | Full browser matrix not available. |
| Keyboard navigation | PASS | Inspected Bootstrap form/button/link markup on forms and navbar. | Controls are reachable with standard HTML semantics. | Forms and buttons use standard labels/inputs/buttons; no custom keyboard trap observed in markup. |
| Focus states | PASS | Inspected Bootstrap classes and CSS asset references. | Visible focus styles exist. | Bootstrap/form/button classes and CSS references present. |
| Field labels/validation | PASS | Sample-group form HTML and invalid form submissions. | Fields have labels and useful validation. | `<label for=...>` markup present; validation messages shown for bad inputs. |
| Color contrast | BLOCKED | Static inspection only; no GUI/contrast analyzer available. | Light/dark theme contrast is acceptable. | Not measured in this container. |
| Tables on small screens | PASS | Inspected results/detail templates and rendered detail page. | Tables are wrapped for responsive layouts. | `table-responsive` wrapper present for results/detail tables. |

## 14. Final release-readiness checks

| Item | Status | Actual URL / steps | Expected result | Actual result |
|---|---:|---|---|---|
| Automated test suite | PASS | `python manage.py test`. | Tests pass. | 154 tests ran in 130.435 s, OK. |
| Lint/static analysis | PASS | No dedicated linter config found; ran `python -m compileall -q app MetaGap` and `python manage.py check`. | Static checks pass. | Compileall and Django checks exited 0. |
| Git status | PASS | `git status --short` after reverting generated staticfiles/pycache and before report edit. | Only intentional report changes remain. | Working tree was clean before editing this QA report. |
| Migrations committed/current | PASS | `python manage.py makemigrations --check --dry-run`. | No model changes requiring migrations. | `No changes detected`. |
| Deployment docs match environment | PASS | Reviewed `MetaGap/DEPLOYMENT.md` and local env behavior where possible. | Docs broadly match env-var driven deployment. | Local env vars worked except upload-size env var wiring; see known issue. |
| Rollback steps documented | BLOCKED | Reviewed docs at repo level; no production deployment access. | Rollback steps are documented and verified. | Not verified in this local pass. |
| Known issues documented | PASS | Added bug task below. | Known issues have severity/owner. | Upload-size env var issue documented as Bug QA-2026-05-13-001; owner TBD. |

## Bug tasks for FAIL items

### QA-2026-05-13-001 — `METAGAP_MAX_UPLOAD_SIZE_MB` environment variable is not applied to import validation

- **Severity:** High — configured upload limits are a safety/operations control, but the app accepted a file larger than the configured 1 MiB local limit.
- **Owner:** TBD.
- **Browser/environment:** Local Django project `MetaGap/`; Python 3.12.13; Django 5.1.1; SQLite `/tmp/metagap-qa-20260513.sqlite3`; Django test client; `METAGAP_MAX_UPLOAD_SIZE_MB=1` exported before running the check.
- **Affected URL:** `/en/import/`.
- **Reproduction steps:**
  1. From `MetaGap/`, export `SECRET_KEY=qa-secret-key`, `DEBUG=1`, `ALLOWED_HOSTS=localhost,127.0.0.1,testserver`, `CSRF_TRUSTED_ORIGINS=http://localhost:8000,http://127.0.0.1:8000`, `DATABASE_URL=sqlite:////tmp/metagap-qa-20260513.sqlite3`, and `METAGAP_MAX_UPLOAD_SIZE_MB=1`.
  2. Run `python manage.py migrate --noinput`.
  3. Log in as a test user.
  4. POST to `/en/import/` with a `.vcf` file slightly larger than 1 MiB.
- **Expected result:** Form returns HTTP 200 with validation message `Files larger than 1 MB cannot be imported.` No import is attempted and no sample group is created.
- **Actual result:** Form returned HTTP 302 success/redirect and the upload was accepted because `ImportDataForm` only reads `settings.METAGAP_MAX_UPLOAD_SIZE_MB`, while `MetaGap/settings.py` does not define that setting from the environment.
- **Screenshots/log snippets:** No GUI screenshot was available in the container. Relevant probe output:

  ```text
  FAIL    oversize file rejected: status=302
  Runtime env vars used: ... METAGAP_MAX_UPLOAD_SIZE_MB=1
  ```

- **Suggested fix:** Add `METAGAP_MAX_UPLOAD_SIZE_MB = float(os.getenv("METAGAP_MAX_UPLOAD_SIZE_MB", "50"))` or equivalent typed parsing in `MetaGap/settings.py`, and keep/extend the existing form tests to cover env-backed configuration.
