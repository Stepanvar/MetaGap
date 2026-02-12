# Project Structure Cleanup - Summary
**Date:** 2026-02-12
**Status:** ✅ COMPLETED

## Objective
Clarify and enforce project structure rules: Django web app code vs non-Django utilities.

## User Requirement
> "All code that is related to my web app but not should be inside it is in the /home/szuev/MetaGap/MetaGap/MetagapUserCode and only in it."

## Actions Taken

### 1. Removed Loose Files from Repository Root ✅
Deleted 3 test VCF files that were incorrectly placed in `/home/szuev/MetaGap/`:
- `cohort_final.vcf.gz`
- `cohort_final(1).vcf.gz`
- `common_all_20180531_papu.vcf.gz`

**Correct location:** Test VCF files belong in `/home/szuev/MetaGap/MetaGap/MetagapUserCode/demo_vcfs/`

### 2. Updated CLAUDE.md ✅
Added explicit structure rules and full directory tree showing:
- Repository root structure
- Django project location
- MetagapUserCode location
- locale/ folders for translations
- static/ and templates/ for frontend

**Key addition:**
```markdown
**Structure Rules:**
- Django web app code: /home/szuev/MetaGap/MetaGap/{app,MetaGap}/ ONLY
- Non-Django utilities/tools: /home/szuev/MetaGap/MetaGap/MetagapUserCode/ ONLY
- No loose files in repository root except documentation (CLAUDE.md, README.md)
- Test VCF files belong in MetagapUserCode/demo_vcfs/, not in root
```

### 3. Updated .claude/docs/project-overview.md ✅
Expanded directory tree to show full structure with:
- Repository root (`/home/szuev/MetaGap/`)
- Documentation folders (`docs/`, `.claude/docs/`)
- Django project root (`MetaGap/`)
- All subdirectories properly organized

Added explicit structure rules section matching CLAUDE.md.

### 4. Verified merge_vcf Location ✅
Confirmed `merge_vcf/` is correctly located in:
```
/home/szuev/MetaGap/MetaGap/MetagapUserCode/merge_vcf/
```

Contains:
- `cli.py` - Entry point
- `__init__.py`
- `__main__.py`

Entry point documented as: `python -m MetaGap.MetagapUserCode.merge_vcf.cli`

### 5. Verified presentation Folder Removed ✅
User confirmed they removed the `presentation/` folder (was a failed test, not needed for MVP).

## Final Structure

```
/home/szuev/MetaGap/                    # Repository root
├── CLAUDE.md                           # AI agent guidance (updated ✅)
├── README.md                           # User documentation
├── LICENSE.txt
├── docs/                               # Planning and documentation
│   └── plans/
├── .claude/                            # Claude Code config
│   └── docs/                           # Architecture references (updated ✅)
│       ├── project-overview.md
│       ├── models.md
│       ├── workflows.md
│       └── architecture.md
└── MetaGap/                            # Django project root
    ├── manage.py
    ├── requirements.txt
    ├── db.sqlite3
    ├── MetaGap/                        # Django settings package
    │   ├── settings.py
    │   ├── urls.py
    │   ├── context_processors.py
    │   └── locale/                     # Project-level translations
    ├── app/                            # Main Django application
    │   ├── models.py
    │   ├── views.py
    │   ├── forms.py
    │   ├── filters.py
    │   ├── tables.py
    │   ├── urls.py
    │   ├── signals.py
    │   ├── services/                   # VCF import pipeline
    │   ├── config/
    │   ├── locale/                     # App-level translations
    │   ├── templates/
    │   ├── static/
    │   ├── tests/
    │   └── fixtures/
    └── MetagapUserCode/                # ✅ All non-Django utilities
        ├── merge_vcf/                  # Standalone VCF merge CLI
        │   ├── cli.py
        │   ├── __init__.py
        │   └── __main__.py
        ├── demo_vcfs/                  # Sample VCF files
        ├── tests/                      # Tests for utilities
        ├── environment.yml             # Conda environment
        ├── joint_call.sh               # Shell scripts
        ├── mk_demo_vcf.sh
        ├── requirements.txt            # Utility dependencies
        └── ...                         # Other test data/scripts
```

## Structure Enforcement

### ✅ Django Web App (ONLY in these locations)
- `/home/szuev/MetaGap/MetaGap/app/` - Main Django app
- `/home/szuev/MetaGap/MetaGap/MetaGap/` - Settings package

### ✅ Non-Django Utilities (ONLY in this location)
- `/home/szuev/MetaGap/MetaGap/MetagapUserCode/` - CLI tools, scripts, test data

### ✅ Documentation (Repository root only)
- `/home/szuev/MetaGap/CLAUDE.md`
- `/home/szuev/MetaGap/README.md`
- `/home/szuev/MetaGap/docs/`
- `/home/szuev/MetaGap/.claude/`

### ❌ Prohibited
- No loose VCF files in repository root
- No test data in `/home/szuev/MetaGap/` (use `MetagapUserCode/demo_vcfs/`)
- No utility scripts in Django app folders
- No presentation/demo folders (were removed)

## Git Status

**Commit:** `dd9ff5f6` - "docs: refine project structure and update documentation"

**Changes:**
- Modified: `CLAUDE.md`
- Modified: `.claude/docs/project-overview.md`
- Deleted: `cohort_final.vcf.gz`
- Deleted: `cohort_final(1).vcf.gz`
- Deleted: `common_all_20180531_papu.vcf.gz`

**Pushed to:** `origin/master` ✅

## Next Steps

Structure is now clean and properly documented. Ready to proceed with **Phase 2: Add SampleGroupCreateView**.
