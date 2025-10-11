# app/models.py

from typing import Any, Dict

from django.conf import settings
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.db import models


class OrganizationProfile(models.Model):
    user = models.OneToOneField(
        User,
        on_delete=models.CASCADE,
        related_name="organization_profile",
    )
    organization_name = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self) -> str:
        return f"{self.user.username}'s Organization Profile"


class ReferenceGenomeBuild(models.Model):
    """Reference genome build metadata."""

    build_name = models.CharField(max_length=100)
    build_version = models.CharField(max_length=100, blank=True, null=True)
    additional_info = models.JSONField(blank=True, null=True)

    def __str__(self) -> str:
        return self.build_name


class GenomeComplexity(models.Model):
    """Genome complexity information."""

    size = models.CharField(max_length=50, blank=True, null=True)
    ploidy = models.CharField(max_length=50, blank=True, null=True)
    gc_content = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self) -> str:
        return f"Size: {self.size}, Ploidy: {self.ploidy}, GC: {self.gc_content}"


class SampleOrigin(models.Model):
    """Sample origin and collection information."""

    tissue = models.CharField(max_length=100, blank=True, null=True)
    collection_method = models.CharField(max_length=100, blank=True, null=True)
    storage_conditions = models.CharField(max_length=100, blank=True, null=True)
    time_stored = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self) -> str:
        return f"Tissue: {self.tissue}, Collection Method: {self.collection_method}"


class MaterialType(models.Model):
    """Material type and integrity information."""

    MATERIAL_CHOICES = [
        ("DNA", "DNA"),
        ("RNA", "RNA"),
        ("cDNA", "cDNA"),
    ]

    material_type = models.CharField(
        max_length=50,
        choices=MATERIAL_CHOICES,
        blank=True,
        null=True,
    )
    integrity_number = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self) -> str:
        return f"Type: {self.material_type}, Integrity Number: {self.integrity_number}"


class LibraryConstruction(models.Model):
    """Library construction details."""

    kit = models.CharField(max_length=100, blank=True, null=True)
    fragmentation = models.CharField(max_length=100, blank=True, null=True)
    adapter_ligation_efficiency = models.CharField(max_length=50, blank=True, null=True)
    pcr_cycles = models.IntegerField(blank=True, null=True)

    def __str__(self) -> str:
        return f"Kit: {self.kit}, Fragmentation: {self.fragmentation}"


def _format_attributes(*attributes: tuple[str, Any]) -> str:
    """Return a comma separated list of populated attribute labels."""

    parts = []
    for label, value in attributes:
        if value not in (None, ""):
            parts.append(f"{label}: {value}" if label else str(value))

    return ", ".join(parts) if parts else "Not provided"


class SequencingInstrument(models.Model):
    """Base model representing sequencing instrument-related details."""

    instrument = models.CharField(max_length=100)
    flow_cell = models.CharField(max_length=100, blank=True, null=True)
    additional_info = models.JSONField(blank=True, null=True)

    class Meta:
        abstract = True


class IlluminaSeq(SequencingInstrument):
    """Illumina sequencing specifics."""

    channel_method = models.CharField(max_length=50, blank=True, null=True)
    cluster_density = models.CharField(max_length=50, blank=True, null=True)
    qc_software = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("", self.instrument),
            ("Flow Cell", self.flow_cell),
            ("Channel method", self.channel_method),
            ("Cluster density", self.cluster_density),
            ("QC software", self.qc_software),
        )


class OntSeq(SequencingInstrument):
    """Oxford Nanopore sequencing specifics."""

    flow_cell_version = models.CharField(max_length=50, blank=True, null=True)
    pore_type = models.CharField(max_length=50, blank=True, null=True)
    bias_voltage = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("", self.instrument),
            ("Flow Cell Version", self.flow_cell_version),
            ("Pore type", self.pore_type),
            ("Bias voltage", self.bias_voltage),
        )


class PacBioSeq(SequencingInstrument):
    """PacBio sequencing specifics."""

    smrt_cell_type = models.CharField(max_length=50, blank=True, null=True)
    zmw_density = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("", self.instrument),
            ("SMRT Cell Type", self.smrt_cell_type),
            ("ZMW density", self.zmw_density),
        )


class IonTorrentSeq(SequencingInstrument):
    """Ion Torrent sequencing specifics."""

    chip_type = models.CharField(max_length=50, blank=True, null=True)
    ph_calibration = models.CharField(max_length=50, blank=True, null=True)
    flow_order = models.CharField(max_length=50, blank=True, null=True)
    ion_sphere_metrics = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("", self.instrument),
            ("Chip Type", self.chip_type),
            ("pH calibration", self.ph_calibration),
            ("Flow order", self.flow_order),
            ("Ion sphere metrics", self.ion_sphere_metrics),
        )


class BioinfoAlignment(models.Model):
    """Bioinformatics alignment settings."""

    tool = models.CharField(max_length=100, blank=True, null=True)
    params = models.CharField(max_length=255, blank=True, null=True)
    ref_genome_version = models.CharField(max_length=100, blank=True, null=True)
    recalibration_settings = models.CharField(max_length=255, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("Tool", self.tool or "Unknown"),
            ("Parameters", self.params),
            ("Ref Genome Version", self.ref_genome_version),
            ("Recalibration settings", self.recalibration_settings),
        )


class BioinfoVariantCalling(models.Model):
    """Bioinformatics variant calling settings."""

    tool = models.CharField(max_length=100, blank=True, null=True)
    version = models.CharField(max_length=50, blank=True, null=True)
    filtering_thresholds = models.CharField(max_length=255, blank=True, null=True)
    duplicate_handling = models.CharField(max_length=50, blank=True, null=True)
    mq = models.CharField(max_length=50, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("Tool", self.tool),
            ("Version", self.version),
            ("Filtering thresholds", self.filtering_thresholds),
            ("Duplicate handling", self.duplicate_handling),
            ("MQ", self.mq),
        )


class BioinfoPostProc(models.Model):
    """Bioinformatics post-processing settings."""

    normalization = models.CharField(max_length=100, blank=True, null=True)
    harmonization = models.CharField(max_length=100, blank=True, null=True)

    def __str__(self) -> str:
        return _format_attributes(
            ("Normalization", self.normalization),
            ("Harmonization", self.harmonization),
        )


class InputQuality(models.Model):
    """Quality metrics of the input biological material."""

    a260_a280 = models.FloatField(blank=True, null=True)
    a260_a230 = models.FloatField(blank=True, null=True)
    dna_concentration = models.FloatField(blank=True, null=True)
    rna_concentration = models.FloatField(blank=True, null=True)
    notes = models.TextField(blank=True, null=True)
    additional_metrics = models.JSONField(blank=True, null=True)

    def __str__(self) -> str:
        return f"Input quality #{self.pk}"


class Info(models.Model):
    """Representation of INFO fields in a VCF file."""

    PLACEHOLDER_STRINGS = {".", ""}
    NUMERIC_FIELD_TYPES = (
        models.IntegerField,
        models.FloatField,
        models.DecimalField,
    )

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
    qd = models.CharField(max_length=50, blank=True, null=True)
    fs = models.CharField(max_length=50, blank=True, null=True)
    sor = models.CharField(max_length=50, blank=True, null=True)
    ns = models.CharField(max_length=50, blank=True, null=True)
    sb = models.CharField(max_length=50, blank=True, null=True)
    additional = models.JSONField(blank=True, null=True)

    def __str__(self) -> str:
        return f"Info for {self.pk}"

    def save(self, *args: Any, **kwargs: Any) -> None:
        """Normalize placeholder numeric values before persisting."""

        self._normalize_placeholder_values()
        super().save(*args, **kwargs)

    def _normalize_placeholder_values(self) -> None:
        numeric_field_names = self._normalize_numeric_fields()
        self._normalize_additional_numeric_entries(numeric_field_names)

    def _normalize_numeric_fields(self) -> set[str]:
        numeric_fields: set[str] = set()
        for field in self._meta.concrete_fields:
            if field.primary_key:
                continue
            if not isinstance(field, self.NUMERIC_FIELD_TYPES):
                continue

            numeric_fields.add(field.name.lower())
            raw_value = getattr(self, field.attname)
            normalized_value = self._coerce_placeholder_value(raw_value)
            if normalized_value is not None:
                try:
                    normalized_value = field.to_python(normalized_value)
                except (TypeError, ValueError, ValidationError):
                    normalized_value = None
            setattr(self, field.attname, normalized_value)

        return numeric_fields

    def _normalize_additional_numeric_entries(self, numeric_field_names: set[str]) -> None:
        if not isinstance(self.additional, dict):
            return

        cleaned_additional: Dict[Any, Any] = {}
        modified = False

        for key, value in self.additional.items():
            lookup_key = str(key).lower()
            if lookup_key in numeric_field_names:
                coerced_value = self._coerce_placeholder_value(value)
                if coerced_value is None:
                    modified = True
                    continue
                if coerced_value != value:
                    modified = True
                    value = coerced_value

            cleaned_additional[key] = value

        if cleaned_additional:
            if modified or cleaned_additional != self.additional:
                self.additional = cleaned_additional
        elif self.additional:
            self.additional = None

    @classmethod
    def _coerce_placeholder_value(cls, value: Any) -> Any:
        if isinstance(value, str) and value.strip() in cls.PLACEHOLDER_STRINGS:
            return None
        return value


class Format(models.Model):
    """Serialized representation of FORMAT fields in a VCF file."""

    genotype = models.CharField(max_length=50, blank=True, null=True)
    payload = models.JSONField(blank=True, null=True)

    def __str__(self) -> str:
        return f"Format for {self.pk}"

    @property
    def fields(self) -> Dict[str, Any]:
        """Return structured FORMAT key/value pairs if available."""
        payload: Dict[str, Any] = self.payload or {}
        return payload.get("fields") or {}

    @property
    def additional(self) -> Dict[str, Any]:
        """Return unmapped FORMAT key/value pairs."""

        payload = self.payload or {}
        return payload.get("additional") or {}


class SampleGroup(models.Model):
    """A group of samples along with the associated metadata."""

    class SequencingPlatform(models.TextChoices):
        ILLUMINA = "illumina", "Illumina"
        ONT = "ont", "ONT"
        PACBIO = "pacbio", "PacBio"
        ION_TORRENT = "ion_torrent", "Ion Torrent"

    PLATFORM_FIELD_MAP = {
        SequencingPlatform.ILLUMINA: "illumina_seq",
        SequencingPlatform.ONT: "ont_seq",
        SequencingPlatform.PACBIO: "pacbio_seq",
        SequencingPlatform.ION_TORRENT: "iontorrent_seq",
    }

    name = models.CharField(max_length=255)
    created_by = models.ForeignKey(OrganizationProfile, on_delete=models.CASCADE)
    doi = models.CharField(max_length=100, blank=True, null=True)
    source_lab = models.CharField(max_length=255, blank=True, null=True)
    contact_email = models.EmailField(blank=True, null=True)
    contact_phone = models.CharField(max_length=50, blank=True, null=True)
    total_samples = models.IntegerField(blank=True, null=True)
    inclusion_criteria = models.TextField(blank=True, null=True)
    exclusion_criteria = models.TextField(blank=True, null=True)
    input_quality = models.OneToOneField(
        InputQuality,
        on_delete=models.SET_NULL,
        blank=True,
        null=True,
        related_name="sample_group",
    )
    comments = models.TextField(blank=True, null=True)
    additional_metadata = models.JSONField(blank=True, null=True)

    reference_genome_build = models.ForeignKey(
        ReferenceGenomeBuild,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    genome_complexity = models.ForeignKey(
        GenomeComplexity,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    sample_origin = models.ForeignKey(
        SampleOrigin,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    material_type = models.ForeignKey(
        MaterialType,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    library_construction = models.ForeignKey(
        LibraryConstruction,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    sequencing_platform = models.CharField(
        max_length=32,
        choices=SequencingPlatform.choices,
        blank=True,
        null=True,
    )
    illumina_seq = models.ForeignKey(
        IlluminaSeq,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    ont_seq = models.ForeignKey(
        OntSeq,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    pacbio_seq = models.ForeignKey(
        PacBioSeq,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    iontorrent_seq = models.ForeignKey(
        IonTorrentSeq,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    bioinfo_alignment = models.ForeignKey(
        BioinfoAlignment,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    bioinfo_variant_calling = models.ForeignKey(
        BioinfoVariantCalling,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    bioinfo_post_proc = models.ForeignKey(
        BioinfoPostProc,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )

    def __str__(self) -> str:
        return self.name

    @property
    def active_sequencing_platform(self) -> SequencingPlatform | None:
        """Return the active sequencing platform choice if configured."""

        value = self.sequencing_platform
        if not value:
            return None
        try:
            return self.SequencingPlatform(value)
        except ValueError:
            return None

    def get_active_platform_instance(self) -> models.Model | None:
        """Return the related sequencing instrument instance for the active platform."""

        platform = self.active_sequencing_platform
        if platform is None:
            return None

        field_name = self.PLATFORM_FIELD_MAP.get(platform)
        if not field_name:
            return None
        return getattr(self, field_name, None)

    def delete(self, using: Any | None = None, keep_parents: bool = False) -> None:
        """Delete the sample group along with unshared related metadata."""

        related_field_names = [
            "reference_genome_build",
            "genome_complexity",
            "sample_origin",
            "material_type",
            "library_construction",
            "illumina_seq",
            "ont_seq",
            "pacbio_seq",
            "iontorrent_seq",
            "bioinfo_alignment",
            "bioinfo_variant_calling",
            "bioinfo_post_proc",
            "input_quality",
        ]

        related_instances: list[models.Model] = []
        for field_name in related_field_names:
            instance = getattr(self, field_name, None)
            if instance is None:
                continue

            field = self._meta.get_field(field_name)
            if isinstance(field, models.OneToOneField):
                related_instances.append(instance)
                continue

            shared = (
                SampleGroup.objects.filter(**{field_name: instance})
                .exclude(pk=self.pk)
                .exists()
            )
            if not shared:
                related_instances.append(instance)

        super().delete(using=using, keep_parents=keep_parents)

        for instance in related_instances:
            instance.delete()


class AlleleFrequency(models.Model):
    """Allele frequency data along with variant information."""

    sample_group = models.ForeignKey(
        SampleGroup,
        on_delete=models.CASCADE,
        related_name="allele_frequencies",
    )
    chrom = models.CharField(max_length=10)
    pos = models.IntegerField()
    variant_id = models.CharField(max_length=100, blank=True, null=True)
    ref = models.CharField(max_length=50)
    alt = models.CharField(max_length=50)
    qual = models.FloatField(blank=True, null=True)
    filter = models.CharField(max_length=50, blank=True, null=True)
    info = models.ForeignKey(Info, on_delete=models.CASCADE, blank=True, null=True)
    format = models.ForeignKey(Format, on_delete=models.CASCADE, blank=True, null=True)
    comments = models.TextField(blank=True, null=True)

    class Meta:
        indexes = [
            models.Index(fields=["sample_group", "chrom", "pos"]),
            models.Index(fields=["variant_id"]),
        ]
        constraints = [
            models.UniqueConstraint(
                fields=["sample_group", "chrom", "pos", "ref", "alt"],
                name="allele_frequency_unique_variant",
            )
        ]

    def __str__(self) -> str:
        return f"{self.chrom}:{self.pos} {self.ref}>{self.alt}"
