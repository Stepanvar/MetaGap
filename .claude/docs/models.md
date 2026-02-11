# MetaGap Data Models

All models defined in `MetaGap/app/models.py`.

## 1. User & Organization

### OrganizationProfile

OneToOne relationship with Django `User`. Required for VCF imports.

**Fields:**
- `user`: OneToOneField(User, CASCADE, related_name="organization_profile")
- `organization_name`: CharField(255, blank/null)

**Business logic:**
- Auto-created via `post_save` signal in `app/signals.py:create_user_profile`
- Signal checks `kwargs.get("raw", False)` to skip fixture loading
- `VCFImporter.import_file()` raises `ImporterConfigurationError` if missing

**Related by:**
- `SampleGroup.created_by` → OrganizationProfile (CASCADE)

---

## 2. SampleGroup (Central Hub)

### SampleGroup

The main entity. Links all metadata models together.

**Fields:**
- `name`: CharField(255)
- `created_by`: ForeignKey(OrganizationProfile, CASCADE)
- `doi`: CharField(100, blank/null)
- `source_lab`: CharField(255, blank/null)
- `contact_email`: EmailField(blank/null)
- `contact_phone`: CharField(50, blank/null)
- `total_samples`: IntegerField(blank/null)
- `inclusion_criteria`: TextField(blank/null)
- `exclusion_criteria`: TextField(blank/null)
- `comments`: TextField(blank/null)
- `additional_metadata`: JSONField(blank/null) — catches unmapped VCF headers
- `sequencing_platform`: CharField(32, choices=SequencingPlatform.choices, blank/null)

**Metadata relationships (all SET_NULL):**
- `reference_genome_build`: ForeignKey(ReferenceGenomeBuild)
- `genome_complexity`: ForeignKey(GenomeComplexity)
- `sample_origin`: ForeignKey(SampleOrigin)
- `material_type`: ForeignKey(MaterialType)
- `library_construction`: ForeignKey(LibraryConstruction)
- `input_quality`: OneToOneField(InputQuality)

**Sequencing platform relationships (all SET_NULL):**
- `illumina_seq`: ForeignKey(IlluminaSeq)
- `ont_seq`: ForeignKey(OntSeq)
- `pacbio_seq`: ForeignKey(PacBioSeq)
- `iontorrent_seq`: ForeignKey(IonTorrentSeq)

**Bioinformatics relationships (all SET_NULL):**
- `bioinfo_alignment`: ForeignKey(BioinfoAlignment)
- `bioinfo_variant_calling`: ForeignKey(BioinfoVariantCalling)
- `bioinfo_post_proc`: ForeignKey(BioinfoPostProc)

**Reverse relationships:**
- `allele_frequencies` ← AlleleFrequency (CASCADE)

**Business logic:**
- `delete()` method (lines 595-636): Deletes unshared metadata instances
  - Iterates through all FK/OneToOne fields
  - Checks if metadata instance is shared by other SampleGroups
  - Only deletes if not shared (protects reusable metadata)
  - Always deletes OneToOneField instances (InputQuality)
- `get_active_sequencing_platform()`: Returns (label, instance) tuple
  - Checks `additional_metadata["active_sequencing_platform"]` first
  - Falls back to first non-null platform FK
- `_normalize_platform_key()`: Normalizes platform identifiers (snake_case, lowercase)

**Class attributes:**
- `SEQUENCING_PLATFORM_FIELDS`: Tuple of (field_name, label) pairs
- `SequencingPlatform`: TextChoices enum
- `PLATFORM_FIELD_MAP`: Maps enum values to field names

---

## 3. Metadata Models

All models use SET_NULL on delete (except InputQuality which uses CASCADE via OneToOne).

### ReferenceGenomeBuild

**Fields:**
- `build_name`: CharField(100)
- `build_version`: CharField(100, blank/null)
- `additional_info`: JSONField(blank/null)

### GenomeComplexity

**Fields:**
- `size`: CharField(50, blank/null)
- `ploidy`: CharField(50, blank/null)
- `gc_content`: CharField(50, blank/null)

### SampleOrigin

**Fields:**
- `tissue`: CharField(100, blank/null)
- `collection_method`: CharField(100, blank/null)
- `storage_conditions`: CharField(100, blank/null)
- `time_stored`: CharField(100, blank/null)

### MaterialType

**Fields:**
- `material_type`: CharField(50, choices=MATERIAL_CHOICES, blank/null)
  - Choices: DNA, RNA, cDNA
- `integrity_number`: CharField(50, blank/null)

### LibraryConstruction

**Fields:**
- `kit`: CharField(100, blank/null)
- `fragmentation`: CharField(100, blank/null)
- `adapter_ligation_efficiency`: CharField(50, blank/null)
- `pcr_cycles`: IntegerField(blank/null)

### InputQuality

**Fields:**
- `a260_a280`: FloatField(blank/null)
- `a260_a230`: FloatField(blank/null)
- `dna_concentration`: FloatField(blank/null)
- `rna_concentration`: FloatField(blank/null)
- `notes`: TextField(blank/null)
- `additional_metrics`: JSONField(blank/null)

**Business logic:**
- OneToOne with SampleGroup (always deleted with group)

---

## 4. Sequencing Platform Models

All inherit from abstract `SequencingInstrument` base class.

### SequencingInstrument (abstract)

**Fields:**
- `instrument`: CharField(100)
- `flow_cell`: CharField(100, blank/null)
- `additional_info`: JSONField(blank/null)

### IlluminaSeq

**Additional fields:**
- `channel_method`: CharField(50, blank/null)
- `cluster_density`: CharField(50, blank/null)
- `qc_software`: CharField(100, blank/null)

### OntSeq (Oxford Nanopore)

**Additional fields:**
- `flow_cell_version`: CharField(50, blank/null)
- `pore_type`: CharField(50, blank/null)
- `bias_voltage`: CharField(50, blank/null)

### PacBioSeq

**Additional fields:**
- `smrt_cell_type`: CharField(50, blank/null)
- `zmw_density`: CharField(50, blank/null)

### IonTorrentSeq

**Additional fields:**
- `chip_type`: CharField(50, blank/null)
- `ph_calibration`: CharField(50, blank/null)
- `flow_order`: CharField(50, blank/null)
- `ion_sphere_metrics`: CharField(50, blank/null)

---

## 5. Bioinformatics Models

### BioinfoAlignment

**Fields:**
- `tool`: CharField(100, blank/null)
- `params`: CharField(255, blank/null)
- `ref_genome_version`: CharField(100, blank/null)
- `recalibration_settings`: CharField(255, blank/null)

### BioinfoVariantCalling

**Fields:**
- `tool`: CharField(100, blank/null)
- `version`: CharField(50, blank/null)
- `filtering_thresholds`: CharField(255, blank/null)
- `duplicate_handling`: CharField(50, blank/null)
- `mq`: CharField(50, blank/null)

### BioinfoPostProc

**Fields:**
- `normalization`: CharField(100, blank/null)
- `harmonization`: CharField(100, blank/null)

---

## 6. Variant Data Models

### AlleleFrequency

Stores variant records from VCF files.

**Fields:**
- `sample_group`: ForeignKey(SampleGroup, CASCADE, related_name="allele_frequencies")
- `chrom`: CharField(10)
- `pos`: IntegerField
- `variant_id`: CharField(100, blank/null)
- `ref`: CharField(50)
- `alt`: CharField(50)
- `qual`: FloatField(blank/null)
- `filter`: CharField(50, blank/null)
- `info`: ForeignKey(Info, CASCADE, blank/null)
- `format`: ForeignKey(Format, CASCADE, blank/null)
- `comments`: TextField(blank/null)

**Constraints:**
- Unique: (sample_group, chrom, pos, ref, alt)

**Indexes:**
- (sample_group, chrom, pos)
- variant_id

### Info

Stores VCF INFO field data.

**Fields (VCF standard):**
- `aa`, `ac`, `af`, `an`, `bq`, `cigar`, `db`, `dp`, `end`, `h2`, `h3`, `mq`, `mq0`, `qd`, `fs`, `sor`, `ns`, `sb`: CharField(50, blank/null)
- `additional`: JSONField(blank/null) — unmapped INFO fields

**Business logic:**
- `save()` calls `_normalize_placeholder_values()` before persist
- `_normalize_placeholder_values()`:
  - Converts placeholder strings (".", "") to None for numeric fields
  - Introspects model to find IntegerField/FloatField/DecimalField
  - Also normalizes numeric values in `additional` JSONField
  - Removes None entries from `additional`

**Constants:**
- `PLACEHOLDER_STRINGS = {".", ""}`
- `NUMERIC_FIELD_TYPES = (IntegerField, FloatField, DecimalField)`

### Format

Stores VCF FORMAT field data.

**Fields:**
- `genotype`: CharField(50, blank/null)
- `payload`: JSONField(blank/null)

**Properties:**
- `fields`: Returns `payload["fields"]` dict or empty dict
- `additional`: Returns `payload["additional"]` dict or empty dict

---

## Utility Functions

### _format_attributes(*attributes)

**Location:** `models.py:110-118`

**Purpose:** Format model attributes for `__str__()` methods

**Logic:**
- Takes tuple pairs: (label, value)
- Filters out None/"" values
- Returns comma-separated "Label: Value" string
- Returns "Not provided" if all empty
