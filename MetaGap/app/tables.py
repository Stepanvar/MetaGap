from typing import Dict, Iterable, List, Sequence

import django_tables2 as tables
from django.db import models

from .models import AlleleFrequency


NUMERIC_FIELD_TYPES = (
    models.AutoField,
    models.BigAutoField,
    models.SmallIntegerField,
    models.IntegerField,
    models.BigIntegerField,
    models.PositiveIntegerField,
    models.PositiveSmallIntegerField,
    models.FloatField,
    models.DecimalField,
    models.DurationField,
)


def _is_numeric_field(field: models.Field) -> bool:
    """Return True when the provided model field stores numeric values."""
    return isinstance(field, NUMERIC_FIELD_TYPES)


def _column_attrs(is_numeric: bool) -> Dict[str, Dict[str, str]]:
    """Generate django-tables2 column attrs with shared styling rules."""
    attrs: Dict[str, Dict[str, str]] = {"th": {"class": "fw-bold"}}
    if is_numeric:
        attrs["td"] = {"class": "text-end"}
    return attrs

def create_dynamic_table(
    primary_model,
    table_name: str = "DynamicTable",
    include_related: bool = True,
    priority_fields: Sequence[str] | None = None,
    exclude_fields: Iterable[str] | None = None,
):
    """
    Dynamically creates a django_tables2.Table subclass for the given primary_model.
    If include_related is True, for each forward relation (ForeignKey/OneToOneField),
    it adds columns for all of the related model's non-auto-created fields using a
    'relation__field' naming scheme.

    :param primary_model: The primary Django model class.
    :param table_name: Name of the generated table class.
    :param include_related: Whether to include fields from related models.
    :param priority_fields: Optional iterable of column names to position at the
        front of the table in the order provided.
    :param exclude_fields: Optional iterable of column names to omit from the
        generated table entirely.
    :return: A dynamically created subclass of tables.Table.
    """
    columns: Dict[str, tables.Column] = {}
    field_order: List[str] = []
    priority_list: List[str] = list(priority_fields or [])
    exclude_set = set(exclude_fields or [])

    def register_column(name: str, column: tables.Column) -> None:
        if name in exclude_set:
            return
        columns[name] = column
        field_order.append(name)

    # Loop over all fields in the primary model.
    for field in primary_model._meta.get_fields():
        # Skip auto-created reverse relations.
        if field.auto_created:
            continue

        if isinstance(field, models.Field):
            is_numeric = _is_numeric_field(field)
            # For related fields (ForeignKey, OneToOneField) include their own fields.
            if field.is_relation and field.related_model is not None:
                if include_related and (field.many_to_one or field.one_to_one):
                    for rel_field in field.related_model._meta.get_fields():
                        if not rel_field.auto_created and isinstance(rel_field, models.Field):
                            col_name = f"{field.name}__{rel_field.name}"
                            rel_is_numeric = _is_numeric_field(rel_field)
                            register_column(
                                col_name,
                                tables.Column(
                                    accessor=f"{field.name}.{rel_field.name}",
                                    attrs=_column_attrs(rel_is_numeric),
                                ),
                            )
                else:
                    # If not expanding related fields, add the field as-is.
                    register_column(
                        field.name, tables.Column(attrs=_column_attrs(is_numeric))
                    )
            else:
                # Standard (non-relation) field.
                register_column(
                    field.name, tables.Column(attrs=_column_attrs(is_numeric))
                )

    seen = set()
    ordered_meta_fields: List[str] = []
    for field_name in priority_list:
        if field_name in columns and field_name not in seen:
            ordered_meta_fields.append(field_name)
            seen.add(field_name)

    for field_name in field_order:
        if field_name not in seen:
            ordered_meta_fields.append(field_name)
            seen.add(field_name)

    # Create a dynamic Meta class.
    Meta = type("Meta", (), {
        "model": primary_model,
        "fields": ordered_meta_fields,
        "attrs": {"class": "table table-hover table-sm"}
    })

    # Attach Meta to the columns dict.
    columns["Meta"] = Meta

    # Create and return the table class.
    return type(table_name, (tables.Table,), columns)


DEFAULT_ALLELE_PRIORITY_FIELDS: Sequence[str] = (
    "chrom",
    "pos",
    "variant_id",
    "ref",
    "alt",
    "qual",
    "filter",
    "info__af",
    "info__ac",
    "info__an",
    "info__dp",
    "info__mq",
)

DEFAULT_ALLELE_EXCLUDE_FIELDS: Sequence[str] = (
    "info__additional",
    "format__payload",
    "sample_group__created_by",
)


def build_allele_frequency_table(
    *,
    priority_extra: Sequence[str] | None = None,
    exclude_extra: Iterable[str] | None = None,
    table_name: str = "AlleleFrequencyTable",
    include_related: bool = True,
):
    """Create an :class:`~django_tables2.Table` for :class:`AlleleFrequency` records.

    Parameters allow callers to extend the shared priority and exclusion lists while
    keeping their defaults aligned across the application.
    """

    priority_fields: List[str] = list(DEFAULT_ALLELE_PRIORITY_FIELDS)
    if priority_extra:
        priority_fields.extend(priority_extra)

    exclude_fields = list(DEFAULT_ALLELE_EXCLUDE_FIELDS)
    if exclude_extra:
        exclude_fields.extend(exclude_extra)

    table_class = create_dynamic_table(
        AlleleFrequency,
        table_name=table_name,
        include_related=include_related,
        priority_fields=priority_fields,
        exclude_fields=exclude_fields,
    )

    friendly_names = {
        "chrom": "Chromosome",
        "pos": "Position (bp)",
        "ref": "Reference allele",
        "alt": "Alternate allele",
        "qual": "QUAL",
        "filter": "Filter status",
        "info__af": "Allele frequency (AF)",
        "info__ac": "Allele count (AC)",
        "info__an": "Alleles analysed (AN)",
        "info__dp": "Read depth (DP)",
        "info__mq": "Mapping quality (MQ)",
        "info__qd": "Quality by depth (QD)",
        "info__fs": "Fisher strand (FS)",
        "info__sor": "Strand odds ratio (SOR)",
        "info__ns": "Samples counted (NS)",
        "info__sb": "Strand bias (SB)",
        "variant_id": "Variant ID",
        "comments": "Variant comments",
        "sample_group__name": "Sample group",
        "sample_group__source_lab": "Source laboratory",
        "sample_group__doi": "DOI",
        "sample_group__contact_email": "Contact email",
        "sample_group__contact_phone": "Contact phone",
        "sample_group__total_samples": "Total samples",
        "sample_group__inclusion_criteria": "Inclusion criteria",
        "sample_group__exclusion_criteria": "Exclusion criteria",
        "sample_group__comments": "Group comments",
        "sample_group__reference_genome_build__build_name": "Reference build",
        "sample_group__reference_genome_build__build_version": "Reference build version",
        "sample_group__genome_complexity__size": "Genome size",
        "sample_group__genome_complexity__ploidy": "Genome ploidy",
        "sample_group__genome_complexity__gc_content": "Genome GC content",
        "sample_group__sample_origin__tissue": "Sample tissue",
        "sample_group__sample_origin__collection_method": "Collection method",
        "sample_group__sample_origin__storage_conditions": "Storage conditions",
        "sample_group__sample_origin__time_stored": "Time stored",
        "sample_group__material_type__material_type": "Material type",
        "sample_group__material_type__integrity_number": "Integrity number",
        "sample_group__library_construction__kit": "Library kit",
        "sample_group__library_construction__fragmentation": "Fragmentation",
        "sample_group__library_construction__adapter_ligation_efficiency": "Adapter ligation efficiency",
        "sample_group__library_construction__pcr_cycles": "PCR cycles",
        "sample_group__illumina_seq__instrument": "Illumina instrument",
        "sample_group__illumina_seq__flow_cell": "Illumina flow cell",
        "sample_group__illumina_seq__channel_method": "Illumina channel method",
        "sample_group__illumina_seq__cluster_density": "Illumina cluster density",
        "sample_group__illumina_seq__qc_software": "Illumina QC software",
        "sample_group__ont_seq__instrument": "ONT instrument",
        "sample_group__ont_seq__flow_cell": "ONT flow cell",
        "sample_group__ont_seq__flow_cell_version": "ONT flow cell version",
        "sample_group__ont_seq__pore_type": "ONT pore type",
        "sample_group__ont_seq__bias_voltage": "ONT bias voltage",
        "sample_group__pacbio_seq__instrument": "PacBio instrument",
        "sample_group__pacbio_seq__flow_cell": "PacBio flow cell",
        "sample_group__pacbio_seq__smrt_cell_type": "PacBio SMRT cell type",
        "sample_group__pacbio_seq__zmw_density": "PacBio ZMW density",
        "sample_group__iontorrent_seq__instrument": "Ion Torrent instrument",
        "sample_group__iontorrent_seq__flow_cell": "Ion Torrent flow cell",
        "sample_group__iontorrent_seq__chip_type": "Ion Torrent chip type",
        "sample_group__iontorrent_seq__ph_calibration": "Ion Torrent pH calibration",
        "sample_group__iontorrent_seq__flow_order": "Ion Torrent flow order",
        "sample_group__iontorrent_seq__ion_sphere_metrics": "Ion Sphere metrics",
        "sample_group__bioinfo_alignment__tool": "Alignment tool",
        "sample_group__bioinfo_alignment__params": "Alignment parameters",
        "sample_group__bioinfo_alignment__ref_genome_version": "Reference genome version",
        "sample_group__bioinfo_alignment__recalibration_settings": "Recalibration settings",
        "sample_group__bioinfo_variant_calling__tool": "Variant calling tool",
        "sample_group__bioinfo_variant_calling__version": "Variant caller version",
        "sample_group__bioinfo_variant_calling__filtering_thresholds": "Variant filtering thresholds",
        "sample_group__bioinfo_variant_calling__duplicate_handling": "Duplicate handling",
        "sample_group__bioinfo_variant_calling__mq": "Variant caller MQ",
        "sample_group__bioinfo_post_proc__normalization": "Post-processing normalization",
        "sample_group__bioinfo_post_proc__harmonization": "Post-processing harmonization",
        "sample_group__input_quality__a260_a280": "Input A260/A280",
        "sample_group__input_quality__a260_a230": "Input A260/A230",
        "sample_group__input_quality__dna_concentration": "DNA concentration",
        "sample_group__input_quality__rna_concentration": "RNA concentration",
        "sample_group__input_quality__notes": "Input quality notes",
    }

    curated_sequence: List[str] = []

    def _append_sequence(name: str) -> None:
        if name in table_class.base_columns and name not in curated_sequence:
            curated_sequence.append(name)

    for field_name in priority_fields:
        _append_sequence(field_name)

    for field_name in [
        "info__qd",
        "info__fs",
        "info__sor",
        "info__ns",
        "info__sb",
        "comments",
    ]:
        _append_sequence(field_name)

    for field_name in [
        "sample_group__name",
        "sample_group__source_lab",
        "sample_group__doi",
        "sample_group__contact_email",
        "sample_group__contact_phone",
        "sample_group__total_samples",
        "sample_group__reference_genome_build__build_name",
        "sample_group__reference_genome_build__build_version",
        "sample_group__genome_complexity__size",
        "sample_group__genome_complexity__ploidy",
        "sample_group__genome_complexity__gc_content",
        "sample_group__sample_origin__tissue",
        "sample_group__sample_origin__collection_method",
        "sample_group__sample_origin__storage_conditions",
        "sample_group__sample_origin__time_stored",
        "sample_group__material_type__material_type",
        "sample_group__material_type__integrity_number",
        "sample_group__library_construction__kit",
        "sample_group__library_construction__fragmentation",
        "sample_group__library_construction__adapter_ligation_efficiency",
        "sample_group__library_construction__pcr_cycles",
        "sample_group__illumina_seq__instrument",
        "sample_group__illumina_seq__flow_cell",
        "sample_group__illumina_seq__channel_method",
        "sample_group__illumina_seq__cluster_density",
        "sample_group__illumina_seq__qc_software",
        "sample_group__ont_seq__instrument",
        "sample_group__ont_seq__flow_cell",
        "sample_group__ont_seq__flow_cell_version",
        "sample_group__ont_seq__pore_type",
        "sample_group__ont_seq__bias_voltage",
        "sample_group__pacbio_seq__instrument",
        "sample_group__pacbio_seq__flow_cell",
        "sample_group__pacbio_seq__smrt_cell_type",
        "sample_group__pacbio_seq__zmw_density",
        "sample_group__iontorrent_seq__instrument",
        "sample_group__iontorrent_seq__flow_cell",
        "sample_group__iontorrent_seq__chip_type",
        "sample_group__iontorrent_seq__ph_calibration",
        "sample_group__iontorrent_seq__flow_order",
        "sample_group__iontorrent_seq__ion_sphere_metrics",
        "sample_group__bioinfo_alignment__tool",
        "sample_group__bioinfo_alignment__params",
        "sample_group__bioinfo_alignment__ref_genome_version",
        "sample_group__bioinfo_alignment__recalibration_settings",
        "sample_group__bioinfo_variant_calling__tool",
        "sample_group__bioinfo_variant_calling__version",
        "sample_group__bioinfo_variant_calling__filtering_thresholds",
        "sample_group__bioinfo_variant_calling__duplicate_handling",
        "sample_group__bioinfo_variant_calling__mq",
        "sample_group__bioinfo_post_proc__normalization",
        "sample_group__bioinfo_post_proc__harmonization",
        "sample_group__input_quality__a260_a280",
        "sample_group__input_quality__a260_a230",
        "sample_group__input_quality__dna_concentration",
        "sample_group__input_quality__rna_concentration",
        "sample_group__input_quality__notes",
    ]:
        _append_sequence(field_name)

    if curated_sequence:
        table_class.Meta.sequence = tuple(curated_sequence + ["..."])

    highlight_columns = {"chrom", "pos", "ref", "alt", "info__af"}

    for column_name, column in table_class.base_columns.items():
        if column_name in friendly_names:
            column.verbose_name = friendly_names[column_name]

        attrs = column.attrs or {}
        th_classes = attrs.setdefault("th", {}).get("class", "")
        th_class_list = [cls for cls in th_classes.split() if cls and cls != "fw-bold"]

        if column_name in highlight_columns:
            th_class_list.extend(["text-uppercase", "small", "fw-semibold"])
        else:
            th_class_list.append("fw-semibold")

        deduped = []
        for cls in th_class_list:
            if cls not in deduped:
                deduped.append(cls)
        attrs["th"]["class"] = " ".join(deduped)

        column.attrs = attrs

    return table_class
