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

Use the reproducible dataset definitions in [QA Dataset Manifest for Import/Search/Export Checks](qa-dataset-manifest.md) before running this section. Record the positive fixture ID (P1-P4 or a generated positive-format variant) and any negative fixture ID (N1-N4) used for each result.

- [ ] Open `/import/` as an authenticated user.
- [ ] Attempt to submit without a file and confirm validation appears.
- [ ] Attempt to upload the unsupported-extension negative fixture N1 and confirm it is rejected.
- [ ] Generate or upload the oversized negative fixture N3 and confirm it is rejected safely.
- [ ] Upload a valid `.vcf` file using the generated positive-format variant from the manifest.
- [ ] Upload a valid compressed `.vcf.gz` fixture P1, P2, P3, or P4, and upload a generated `.vcf.bgz` variant if bgzip coverage is required.
- [ ] Upload a valid generated `.bcf` file if `bcftools` is available in the QA environment.
- [ ] Confirm a successful import redirects to `/profile/` and shows a success message.
- [ ] Confirm imported sample-group metadata appears on the profile and detail pages.
- [ ] Confirm imported allele-frequency records are searchable by comparing results with the expected row tables in the manifest.
- [ ] Review logs for import warnings, validation failures, tracebacks, or slow operations.

## 7. Search workflow checks

Use the expected allele-frequency rows in [QA Dataset Manifest for Import/Search/Export Checks](qa-dataset-manifest.md) for all row-level assertions. Testers should search by visible UI filters only, not by database queries or internal importer behavior.

- [ ] Open `/search/` without filters and confirm the page renders.
- [ ] Search by variant ID.
- [ ] Search by chromosome and position if supported by the current filters.
- [ ] Search by sample-group or metadata keyword.
- [ ] Confirm search results include the expected allele-frequency rows from the manifest for the imported fixture.
- [ ] Confirm result table sorting works.
- [ ] Confirm pagination or DataTables controls work.
- [ ] Confirm the clear/reset search action works.
- [ ] Confirm search performance is acceptable with demo data and fixture P4, the largest repository fixture in the manifest.

## 8. Export workflow checks

Use the export column and row expectations in [QA Dataset Manifest for Import/Search/Export Checks](qa-dataset-manifest.md). The tester should be able to validate output by opening the downloaded CSV/TSV and comparing it to the manifest tables.

- [ ] Open a sample-group detail page with imported variants.
- [ ] Use the export action at `/profile/sample-groups/<id>/export/`.
- [ ] Confirm the downloaded file opens correctly.
- [ ] Confirm exported rows match the manifest expectations for the sample group and do not include another user's data.
- [ ] Confirm anonymous users cannot export data.
- [ ] Confirm users cannot export sample groups owned by another organization.

## 9. Authorization and data-isolation checks

Use fixtures P1 and P2 from [QA Dataset Manifest for Import/Search/Export Checks](qa-dataset-manifest.md) so each test user imports a small dataset with overlapping positions but different alternate alleles and sample identifiers.

- [ ] Create two separate test users or organizations.
- [ ] Import fixture P1 under user/organization A and fixture P2 under user/organization B.
- [ ] Confirm each user only sees their own profile sample groups.
- [ ] Attempt direct URL access to another user's sample-group detail page and confirm it is denied or returns not found.
- [ ] Attempt direct URL access to another user's edit, delete, and export URLs and confirm access is blocked.
- [ ] Confirm public search behavior matches the intended product policy for cross-organization visibility, using the P1/P2 row expectations to identify data leakage or intended cross-organization visibility.

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
