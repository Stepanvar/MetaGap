# app/models.py

from django.db import models
from django.contrib.auth.models import User

class SampleGroup(models.Model):
    """
    Model representing a group of samples along with their metadata.
    """
    name = models.CharField(max_length=255)
    created_by = models.ForeignKey(User, on_delete=models.CASCADE)
    doi = models.CharField(max_length=100, blank=True, null=True)

    # Metadata fields from VCF headers
    source_lab = models.CharField(max_length=255, blank=True, null=True)
    contact_email = models.EmailField(blank=True, null=True)
    contact_phone = models.CharField(max_length=20, blank=True, null=True)
    device_type = models.CharField(max_length=100, blank=True, null=True)
    device_id = models.CharField(max_length=100, blank=True, null=True)
    tissue = models.CharField(max_length=100, blank=True, null=True)
    species = models.CharField(max_length=100, blank=True, null=True)
    storage_temperature = models.CharField(max_length=50, blank=True, null=True)
    storage_duration = models.CharField(max_length=50, blank=True, null=True)
    absence_of_relatives = models.BooleanField(default=False)
    preparation_method = models.CharField(max_length=255, blank=True, null=True)
    preparation_kit = models.CharField(max_length=255, blank=True, null=True)
    is_transcriptome = models.BooleanField(default=False)
    sequencing_type = models.CharField(max_length=100, blank=True, null=True)
    analysis_software = models.CharField(max_length=100, blank=True, null=True)
    analysis_version = models.CharField(max_length=50, blank=True, null=True)
    total_samples = models.IntegerField(blank=True, null=True)
    phenotype = models.CharField(max_length=255, blank=True, null=True)
    inclusion_criteria = models.TextField(blank=True, null=True)
    exclusion_criteria = models.TextField(blank=True, null=True)
    ancestry = models.CharField(max_length=255, blank=True, null=True)

    # Comments
    comments = models.TextField(blank=True, null=True)

    def __str__(self):
        return self.name


class AlleleFrequency(models.Model):
    """
    Model representing allele frequency data along with variant information.
    """
    sample_group = models.ForeignKey(SampleGroup, on_delete=models.CASCADE)
    chrom = models.CharField(max_length=10)
    pos = models.IntegerField()
    variant_id = models.CharField(max_length=100, blank=True, null=True)
    ref = models.CharField(max_length=50)
    alt = models.CharField(max_length=50)
    qual = models.FloatField(blank=True, null=True)
    filter = models.CharField(max_length=50, blank=True, null=True)

    # INFO field as a JSON field to store AF, DP, MQ, etc.
    info = models.JSONField()

    # Additional comments
    comments = models.TextField(blank=True, null=True)

    def __str__(self):
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"


class QualityControlMetrics(models.Model):
    """
    Model representing quality control metrics associated with a sample group.
    """
    sample_group = models.OneToOneField(SampleGroup, on_delete=models.CASCADE)

    # Sample Collection and Preservation
    dna_integrity = models.FloatField(blank=True, null=True)
    rna_integrity = models.FloatField(blank=True, null=True)
    contamination_level = models.CharField(max_length=50, blank=True, null=True)
    storage_conditions = models.CharField(max_length=255, blank=True, null=True)

    # DNA/RNA Extraction and Quantification
    extraction_yield = models.CharField(max_length=50, blank=True, null=True)
    sample_purity = models.CharField(max_length=50, blank=True, null=True)
    library_yield = models.CharField(max_length=50, blank=True, null=True)

    # Library Preparation QC
    fragment_size_distribution = models.CharField(max_length=50, blank=True, null=True)
    pcr_duplication_rate = models.CharField(max_length=50, blank=True, null=True)
    phred_quality_score = models.CharField(max_length=50, blank=True, null=True)
    percent_passing_filter = models.CharField(max_length=50, blank=True, null=True)

    # Sequencing QC
    coverage_depth = models.CharField(max_length=50, blank=True, null=True)
    base_calling_accuracy = models.CharField(max_length=50, blank=True, null=True)
    per_base_sequence_quality = models.CharField(max_length=50, blank=True, null=True)
    per_sequence_gc_content = models.CharField(max_length=50, blank=True, null=True)

    # Alignment and BAM File QC
    mapping_quality = models.CharField(max_length=50, blank=True, null=True)
    alignment_percentage = models.CharField(max_length=50, blank=True, null=True)
    insert_size_distribution = models.CharField(max_length=50, blank=True, null=True)

    # Variant Calling and gVCF QC
    variant_call_quality = models.CharField(max_length=50, blank=True, null=True)
    depth_of_coverage = models.CharField(max_length=50, blank=True, null=True)
    genotype_quality = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"QC Metrics for {self.sample_group.name}"
