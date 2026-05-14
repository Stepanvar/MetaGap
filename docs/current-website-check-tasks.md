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

Product/security policy: `/search/` is an authenticated organization-private catalogue.
Search results must be scoped to allele-frequency records whose `sample_group.created_by`
matches the current user's organization profile; anonymous users must be redirected to
login and must not see private records.

- [ ] Open `/search/` without filters as an authenticated user and confirm the page renders.
- [ ] Confirm anonymous access to `/search/` redirects to login without rendering allele-frequency rows.
- [ ] Search by variant ID.
- [ ] Search by chromosome and position if supported by the current filters.
- [ ] Search by sample-group or metadata keyword within the current organization.
- [ ] Confirm search results include expected allele-frequency rows for the current organization.
- [ ] Confirm search results exclude allele-frequency rows from other organizations.
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

Product/security policy: cross-organization public search is not allowed. `/search/` is
private to the authenticated user's organization and must filter allele-frequency records
through `sample_group__created_by`; this intentionally prevents user A from seeing user
B's variants and prevents anonymous users from seeing private allele-frequency records.

- [ ] Create two separate test users or organizations.
- [ ] Create or import sample groups under each account.
- [ ] Confirm each user only sees their own profile sample groups.
- [ ] Attempt direct URL access to another user's sample-group detail page and confirm it is denied or returns not found.
- [ ] Attempt direct URL access to another user's edit, delete, and export URLs and confirm access is blocked.
- [ ] Confirm `/search/` redirects anonymous users to login without rendering private allele-frequency records.
- [ ] Confirm `/search/` for user A excludes user B's variants, including keyword, variant-ID, chromosome/position, and metadata searches.

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
