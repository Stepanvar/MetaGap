# Project Overview

## Project Summary
MetaGap is a Django-based bioinformatics application designed for managing genetic variant data from VCF files, calculating and storing allele frequencies, and organizing sample groups with comprehensive metadata. The system enables organizations to manage and analyze genetic variant data with rich contextual information about sample groups, sequencing workflows, and bioinformatics processing pipelines. It supports ingestion of VCF files and extraction of structured metadata for variant analysis and cohort management.

## Tech Stack
- **Backend**: Python 3.x, Django 5.1.1
- **Database**: SQLite (with dj-database-url for flexible deployment)
- **Bioinformatics Libraries**: Pysam 0.22.1, VCFPy 0.14.2
- **Web Framework**: Django with Bootstrap 5, django-tables2, django-filter
- **Data Processing**: Pandas (via VCF processing), NumPy
- **Serialization**: PyYAML 6.0.2 for configuration files
- **Deployment**: Gunicorn, Whitenoise for static files

## Key Directories
- `app/services`: Contains VCF ingestion logic (`vcf_importer.py`), metadata parsing (`vcf_metadata.py`), and database writing utilities
- `app/models`: Defines the database schema including `SampleGroup`, `AlleleFrequency`, and related metadata models for sequencing, bioinformatics, and sample information
- `MetagapUserCode`: Houses custom scripts for VCF merging/joint-calling in the `merge_vcf` package with CLI tools, validation, and metadata handling
- `app/fixtures`: Stores demo data (`demo_data.json`) for easy exploration of the interface with pre-populated sample groups and metadata