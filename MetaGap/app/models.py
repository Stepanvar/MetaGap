# app/models.py

from django.conf import settings
from django.contrib.auth.models import User
from django.db import models
from .signals import create_user_profile, save_user_profile
class OrganizationProfile(models.Model):
    user = models.OneToOneField(User, on_delete=models.CASCADE, related_name="organization_profile")
    organization_name = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self):
        return f"{self.user.username}'s Organization Profile"

class ReferenceGenomeBuild(models.Model):
    """
    Model representing reference genome build metadata.
    """
    build_name = models.CharField(max_length=100)
    build_version = models.CharField(max_length=100, blank=True, null=True)
    additional_info = models.JSONField(blank=True, null=True)

    def __str__(self):
        return self.build_name


class GenomeComplexity(models.Model):
    """
    Model representing genome complexity information.
    """
    size = models.CharField(max_length=50, blank=True, null=True)
    ploidy = models.CharField(max_length=50, blank=True, null=True)
    gc_content = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"Size: {self.size}, Ploidy: {self.ploidy}, GC: {self.gc_content}"


class SampleOrigin(models.Model):
    """
    Model representing sample origin and collection information.
    """
    tissue = models.CharField(max_length=100, blank=True, null=True)
    collection_method = models.CharField(max_length=100, blank=True, null=True)
    storage_conditions = models.CharField(max_length=100, blank=True, null=True)
    time_stored = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return f"Tissue: {self.tissue}, Collection Method: {self.collection_method}"


class MaterialType(models.Model):
    """
    Model representing material type and integrity information.
    """
    material_type = models.CharField(max_length=50, choices=[('DNA', 'DNA'), ('RNA', 'RNA'), ('cDNA', 'cDNA')], blank=True, null=True)
    integrity_number = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"Type: {self.material_type}, Integrity Number: {self.integrity_number}"


class LibraryConstruction(models.Model):
    """
    Model representing library construction details.
    """
    kit = models.CharField(max_length=100, blank=True, null=True)
    fragmentation = models.CharField(max_length=100, blank=True, null=True)
    adapter_ligation_efficiency = models.CharField(max_length=50, blank=True, null=True)
    pcr_cycles = models.IntegerField(blank=True, null=True)

    def __str__(self):
        return f"Kit: {self.kit}, Fragmentation: {self.fragmentation}"


class SequencingInstrument(models.Model):
    """
    Base model representing sequencing instrument-related details.
    """
    instrument = models.CharField(max_length=100)
    flow_cell = models.CharField(max_length=100, blank=True, null=True)
    additional_info = models.JSONField(blank=True, null=True)

    class Meta:
        abstract = True


class IlluminaSeq(SequencingInstrument):
    """
    Model representing Illumina sequencing specifics.
    """
    channel_method = models.CharField(max_length=50, blank=True, null=True)
    cluster_density = models.CharField(max_length=50, blank=True, null=True)
    qc_software = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return f"Illumina - {self.instrument}, Flow Cell: {self.flow_cell}"


class OntSeq(SequencingInstrument):
    """
    Model representing Oxford Nanopore sequencing specifics.
    """
    flow_cell_version = models.CharField(max_length=50, blank=True, null=True)
    pore_type = models.CharField(max_length=50, blank=True, null=True)
    bias_voltage = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"ONT - {self.instrument}, Flow Cell Version: {self.flow_cell_version}"


class PacBioSeq(SequencingInstrument):
    """
    Model representing PacBio sequencing specifics.
    """
    smrt_cell_type = models.CharField(max_length=50, blank=True, null=True)
    zmw_density = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"PacBio - {self.instrument}, SMRT Cell Type: {self.smrt_cell_type}"


class IonTorrentSeq(SequencingInstrument):
    """
    Model representing Ion Torrent sequencing specifics.
    """
    chip_type = models.CharField(max_length=50, blank=True, null=True)
    ph_calibration = models.CharField(max_length=50, blank=True, null=True)
    flow_order = models.CharField(max_length=50, blank=True, null=True)
    ion_sphere_metrics = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"Ion Torrent - {self.instrument}, Chip Type: {self.chip_type}"


class PlatformIndependent(models.Model):
    """
    Model for platform-independent sequencing details.
    """
    pooling = models.CharField(max_length=100, blank=True, null=True)
    sequencing_kit = models.CharField(max_length=100, blank=True, null=True)
    base_calling_alg = models.CharField(max_length=100, blank=True, null=True)
    q30 = models.CharField(max_length=50, blank=True, null=True)
    normalized_coverage = models.CharField(max_length=100, blank=True, null=True)
    run_specific_calibration = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return f"Pooling: {self.pooling}, Sequencing Kit: {self.sequencing_kit}"


class BioinfoAlignment(models.Model):
    """
    Model for bioinformatics alignment settings.
    """
    software = models.CharField(max_length=100, blank=True, null=True)
    params = models.CharField(max_length=255, blank=True, null=True)
    ref_genome_version = models.CharField(max_length=100, blank=True, null=True)
    recalibration_settings = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self):
        return f"Software: {self.software}, Ref Genome Version: {self.ref_genome_version}"


class BioinfoVariantCalling(models.Model):
    """
    Model for bioinformatics variant calling settings.
    """
    tool = models.CharField(max_length=100, blank=True, null=True)
    version = models.CharField(max_length=50, blank=True, null=True)
    filtering_thresholds = models.CharField(max_length=255, blank=True, null=True)
    duplicate_handling = models.CharField(max_length=50, blank=True, null=True)
    mq = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self):
        return f"Tool: {self.tool}, Version: {self.version}"


class BioinfoPostProc(models.Model):
    """
    Model for bioinformatics post-processing settings.
    """
    normalization = models.CharField(max_length=100, blank=True, null=True)
    harmonization = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self):
        return f"Normalization: {self.normalization}, Harmonization: {self.harmonization}"

class Info(models.Model):
    """
    Model representing INFO fields in a VCF file.
    """
    aa = models.CharField(max_length=50, blank=True, null=True)
    ac = models.CharField(max_length=50, blank=True, null=True)
    af = models.CharField(max_length=50, blank=True, null=True)
    an = models.CharField(max_length=50, blank=True, null=True)
    bq = models.CharField(max_length=50, blank=True, null=True)
    cigar = models.CharField(max_length=50, blank=True, null=True)
    db = models.CharField(max_length=50, blank=True, null=True)
    dp = models.CharField(max_length=50, blank=True, null=True)
    end = models.CharField(max_length=50, blank=True, null=True)
    h2 = models.CharField(max_length=50, blank=True, null=True)
    h3 = models.CharField(max_length=50, blank=True, null=True)
    mq = models.CharField(max_length=50, blank=True, null=True)
    mq0 = models.CharField(max_length=50, blank=True, null=True)
    ns = models.CharField(max_length=50, blank=True, null=True)
    sb = models.CharField(max_length=50, blank=True, null=True)
    additional = models.JSONField(blank=True, null=True)  # Stores additional INFO fields

    def __str__(self):
        return f"Info for {self.pk}"

class Format(models.Model):
    """
    Model representing FORMAT fields in a VCF file.
    """
    ad = models.CharField(max_length=50, blank=True, null=True)
    adf = models.CharField(max_length=50, blank=True, null=True)
    adr = models.CharField(max_length=50, blank=True, null=True)
    dp = models.CharField(max_length=50, blank=True, null=True)
    ec = models.CharField(max_length=50, blank=True, null=True)
    ft = models.CharField(max_length=50, blank=True, null=True)
    gl = models.CharField(max_length=50, blank=True, null=True)
    gp = models.CharField(max_length=50, blank=True, null=True)
    gq = models.CharField(max_length=50, blank=True, null=True)
    gt = models.CharField(max_length=50, blank=True, null=True)
    hq = models.CharField(max_length=50, blank=True, null=True)
    mq = models.CharField(max_length=50, blank=True, null=True)
    pl = models.CharField(max_length=50, blank=True, null=True)
    pq = models.CharField(max_length=50, blank=True, null=True)
    ps = models.CharField(max_length=50, blank=True, null=True)
    additional = models.JSONField(blank=True, null=True)  # Stores additional FORMAT fields

    def __str__(self):
        return f"Format for {self.pk}"


class AlleleFrequency(models.Model):
    """
    Model representing allele frequency data along with variant information.
    """
    chrom = models.CharField(max_length=10)
    pos = models.IntegerField()
    variant_id = models.CharField(max_length=100, blank=True, null=True)
    ref = models.CharField(max_length=50)
    alt = models.CharField(max_length=50)
    qual = models.FloatField(blank=True, null=True)
    filter = models.CharField(max_length=50, blank=True, null=True)

    info = models.ForeignKey(Info, on_delete=models.CASCADE, blank=True, null=True)
    format = models.ForeignKey(Format, on_delete=models.CASCADE, blank=True, null=True)


    # Additional comments
    comments = models.TextField(blank=True, null=True)

    def __str__(self):
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"


class SampleGroup(models.Model):
    """
    Model representing a group of samples along with their metadata.
    """
    name = models.CharField(max_length=255)
    created_by = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name="sample_groups",
    )
    tissue = models.CharField(max_length=100, blank=True, null=True)
    collection_method = models.CharField(max_length=100, blank=True, null=True)
    storage_conditions = models.CharField(max_length=100, blank=True, null=True)
    comments = models.TextField(blank=True, null=True)

    # Foreign keys to other models
    allele_frequency = models.ForeignKey(AlleleFrequency, on_delete=models.SET_NULL, null=True, blank=True)
    reference_genome_build = models.ForeignKey(ReferenceGenomeBuild, on_delete=models.SET_NULL, null=True, blank=True)
    genome_complexity = models.ForeignKey(GenomeComplexity, on_delete=models.SET_NULL, null=True, blank=True)
    sample_origin = models.ForeignKey(SampleOrigin, on_delete=models.SET_NULL, null=True, blank=True)
    material_type = models.ForeignKey(MaterialType, on_delete=models.SET_NULL, null=True, blank=True)
    library_construction = models.ForeignKey(LibraryConstruction, on_delete=models.SET_NULL, null=True, blank=True)
    illumina_seq = models.ForeignKey(IlluminaSeq, on_delete=models.SET_NULL, null=True, blank=True)
    ont_seq = models.ForeignKey(OntSeq, on_delete=models.SET_NULL, null=True, blank=True)
    pacbio_seq = models.ForeignKey(PacBioSeq, on_delete=models.SET_NULL, null=True, blank=True)
    iontorrent_seq = models.ForeignKey(IonTorrentSeq, on_delete=models.SET_NULL, null=True, blank=True)
    platform_independent = models.ForeignKey(PlatformIndependent, on_delete=models.SET_NULL, null=True, blank=True)
    bioinfo_alignment = models.ForeignKey(BioinfoAlignment, on_delete=models.SET_NULL, null=True, blank=True)
    bioinfo_variant_calling = models.ForeignKey(BioinfoVariantCalling, on_delete=models.SET_NULL, null=True, blank=True)
    bioinfo_post_proc = models.ForeignKey(BioinfoPostProc, on_delete=models.SET_NULL, null=True, blank=True)

    def __str__(self):
        return self.name
