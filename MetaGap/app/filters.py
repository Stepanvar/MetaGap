# filters.py
from django import forms
import django_filters
from django.db.models import (
    CharField,
    F,
    FloatField,
    Q,
    Value,
)
from django.db.models.functions import Cast, NullIf
from django.db.models.fields.json import KeyTextTransform

from .models import AlleleFrequency, SampleGroup


class AlleleFrequencyFilter(django_filters.FilterSet):
    query = django_filters.CharFilter(method="universal_search", label="Search")
    chrom = django_filters.CharFilter(method="filter_chrom", label="Chrom")
    pos_min = django_filters.NumberFilter(
        field_name="pos", lookup_expr="gte", label="Position (min)"
    )
    pos_max = django_filters.NumberFilter(
        field_name="pos", lookup_expr="lte", label="Position (max)"
    )
    pos = django_filters.NumberFilter(
        field_name="pos", lookup_expr="exact", label="Position"
    )
    ref = django_filters.CharFilter(
        field_name="ref", lookup_expr="iexact", label="Reference"
    )
    alt = django_filters.CharFilter(
        field_name="alt", lookup_expr="iexact", label="Alternate"
    )
    filter_pass = django_filters.BooleanFilter(
        method="filter_pass_status", label="Filter is PASS"
    )
    qual_min = django_filters.NumberFilter(
        field_name="qual", lookup_expr="gte", label="QUAL (min)"
    )
    qual_max = django_filters.NumberFilter(
        field_name="qual", lookup_expr="lte", label="QUAL (max)"
    )
    af_min = django_filters.NumberFilter(
        method="filter_info_numeric", label="AF (min)"
    )
    af_max = django_filters.NumberFilter(
        method="filter_info_numeric", label="AF (max)"
    )
    dp_min = django_filters.NumberFilter(
        method="filter_info_numeric", label="DP (min)"
    )
    dp_max = django_filters.NumberFilter(
        method="filter_info_numeric", label="DP (max)"
    )
    mq_min = django_filters.NumberFilter(
        method="filter_info_numeric", label="MQ (min)"
    )
    mq_max = django_filters.NumberFilter(
        method="filter_info_numeric", label="MQ (max)"
    )
    qd_min = django_filters.NumberFilter(
        method="filter_info_numeric", label="QD (min)"
    )
    qd_max = django_filters.NumberFilter(
        method="filter_info_numeric", label="QD (max)"
    )
    fs_min = django_filters.NumberFilter(
        method="filter_info_numeric", label="FS (min)"
    )
    fs_max = django_filters.NumberFilter(
        method="filter_info_numeric", label="FS (max)"
    )
    sor_min = django_filters.NumberFilter(
        method="filter_info_numeric", label="SOR (min)"
    )
    sor_max = django_filters.NumberFilter(
        method="filter_info_numeric", label="SOR (max)"
    )

    sample_group_source_lab = django_filters.CharFilter(
        field_name="sample_group__source_lab",
        lookup_expr="icontains",
        label="Source lab",
    )
    sample_group_sample_origin_tissue = django_filters.CharFilter(
        field_name="sample_group__sample_origin__tissue",
        lookup_expr="icontains",
        label="Sample tissue",
    )
    sample_group_bioinfo_variant_calling_tool = django_filters.CharFilter(
        field_name="sample_group__bioinfo_variant_calling__tool",
        lookup_expr="icontains",
        label="Variant calling tool",
    )
    sample_group_bioinfo_variant_calling_version = django_filters.CharFilter(
        field_name="sample_group__bioinfo_variant_calling__version",
        lookup_expr="icontains",
        label="Variant calling version",
    )
    sample_group_bioinfo_alignment_tool = django_filters.CharFilter(
        field_name="sample_group__bioinfo_alignment__tool",
        lookup_expr="icontains",
        label="Alignment tool",
    )
    sample_group_bioinfo_alignment_ref_genome = django_filters.CharFilter(
        field_name="sample_group__bioinfo_alignment__ref_genome_version",
        lookup_expr="icontains",
        label="Alignment reference genome",
    )
    sample_group_bioinfo_alignment_recalibration = django_filters.CharFilter(
        field_name="sample_group__bioinfo_alignment__recalibration_settings",
        lookup_expr="icontains",
        label="Alignment recalibration",
    )

    field_order = [
        "query",
        "chrom",
        "pos",
        "pos_min",
        "pos_max",
        "ref",
        "alt",
        "filter_pass",
        "qual_min",
        "qual_max",
        "af_min",
        "af_max",
        "dp_min",
        "dp_max",
        "mq_min",
        "mq_max",
        "qd_min",
        "qd_max",
        "fs_min",
        "fs_max",
        "sor_min",
        "sor_max",
        "sample_group_source_lab",
        "sample_group_sample_origin_tissue",
        "sample_group_bioinfo_variant_calling_tool",
        "sample_group_bioinfo_variant_calling_version",
        "sample_group_bioinfo_alignment_tool",
        "sample_group_bioinfo_alignment_ref_genome",
        "sample_group_bioinfo_alignment_recalibration",
    ]

    NUMERIC_INFO_EXPRESSIONS = {
        "af": F("info__af"),
        "dp": F("info__dp"),
        "mq": F("info__mq"),
        "qd": KeyTextTransform("QD", "info__additional"),
        "fs": KeyTextTransform("FS", "info__additional"),
        "sor": KeyTextTransform("SOR", "info__additional"),
    }

    class Meta:
        model = AlleleFrequency
        fields = []

    def universal_search(self, queryset, name, value):
        if value:
            value = value.strip()
            filters_q = Q(chrom__icontains=value) | Q(variant_id__icontains=value)
            if value.isdigit():
                filters_q |= Q(pos=int(value))
            return queryset.filter(filters_q)
        return queryset

    @staticmethod
    def _normalize_chrom(value):
        if value is None:
            return None
        normalized = str(value).strip()
        if not normalized:
            return None
        return normalized.upper()

    def filter_chrom(self, queryset, name, value):
        normalized = self._normalize_chrom(value)
        if not normalized:
            return queryset
        return queryset.filter(chrom__iexact=normalized)

    def filter_pass_status(self, queryset, name, value):
        if value is None:
            return queryset

        if value:
            return queryset.filter(filter__iexact="PASS")
        return queryset.exclude(filter__iexact="PASS")

    def filter_info_numeric(self, queryset, name, value):
        if value is None:
            return queryset

        lookup_map = {"min": "gte", "max": "lte"}
        parts = name.rsplit("_", 1)
        if len(parts) != 2 or parts[1] not in lookup_map:
            return queryset

        field_key, bound = parts
        lookup = lookup_map[bound]
        alias = f"numeric_{field_key}"
        expression = self.NUMERIC_INFO_EXPRESSIONS.get(field_key)
        if expression is None:
            return queryset
        cleaned_expression = Cast(
            NullIf(
                NullIf(expression, Value(".", output_field=CharField())),
                Value("", output_field=CharField()),
            ),
            FloatField(),
        )
        queryset = queryset.annotate(**{alias: cleaned_expression})
        return queryset.filter(**{f"{alias}__{lookup}": value})


class AlleleFrequencySearchFilter(django_filters.FilterSet):
    """Comprehensive filter set for refining allele frequency results."""

    query = django_filters.CharFilter(
        method="filter_query",
        label="Keyword",
    )
    chrom = django_filters.CharFilter(
        field_name="chrom",
        lookup_expr="icontains",
        label="Chromosome",
    )
    pos_min = django_filters.NumberFilter(
        field_name="pos",
        lookup_expr="gte",
        label="Position ≥",
    )
    pos_max = django_filters.NumberFilter(
        field_name="pos",
        lookup_expr="lte",
        label="Position ≤",
    )
    ref = django_filters.CharFilter(
        field_name="ref",
        lookup_expr="iexact",
        label="Reference",
    )
    alt = django_filters.CharFilter(
        field_name="alt",
        lookup_expr="iexact",
        label="Alternate",
    )
    pass_only = django_filters.BooleanFilter(
        method="filter_pass_only",
        label="PASS only",
    )
    qual_min = django_filters.NumberFilter(
        field_name="qual",
        lookup_expr="gte",
        label="Quality ≥",
    )
    af_min = django_filters.NumberFilter(
        field_name="info__af",
        method="filter_info_float",
        label="AF ≥",
    )
    ac_min = django_filters.NumberFilter(
        field_name="info__ac",
        method="filter_info_int",
        label="AC ≥",
    )
    an_min = django_filters.NumberFilter(
        field_name="info__an",
        method="filter_info_int",
        label="AN ≥",
    )
    dp_min = django_filters.NumberFilter(
        field_name="info__dp",
        method="filter_info_int",
        label="DP ≥",
    )
    mq_min = django_filters.NumberFilter(
        field_name="info__mq",
        method="filter_info_float",
        label="MQ ≥",
    )
    sample_group_name = django_filters.CharFilter(
        field_name="sample_group__name",
        lookup_expr="icontains",
        label="Sample group",
    )
    sample_group_source_lab = django_filters.CharFilter(
        field_name="sample_group__source_lab",
        lookup_expr="icontains",
        label="Source lab",
    )
    sample_origin_tissue = django_filters.CharFilter(
        field_name="sample_group__sample_origin__tissue",
        lookup_expr="icontains",
        label="Tissue",
    )
    variant_calling_tool = django_filters.CharFilter(
        field_name="sample_group__bioinfo_variant_calling__tool",
        lookup_expr="icontains",
        label="Variant caller",
    )

    class Meta:
        model = AlleleFrequency
        fields = []

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._apply_bootstrap_widgets()

    # ------------------------------------------------------------------
    # Filter helpers
    # ------------------------------------------------------------------
    def filter_query(self, queryset, name, value):
        if not value:
            return queryset

        value = value.strip()
        if not value:
            return queryset

        search_filter = (
            Q(chrom__icontains=value)
            | Q(variant_id__icontains=value)
            | Q(ref__icontains=value)
            | Q(alt__icontains=value)
            | Q(sample_group__name__icontains=value)
            | Q(sample_group__source_lab__icontains=value)
            | Q(sample_group__sample_origin__tissue__icontains=value)
        )

        if value.isdigit():
            search_filter |= Q(pos=int(value))

        return queryset.filter(search_filter)

    def filter_pass_only(self, queryset, name, value):
        truthy_values = {True, "True", "true", "1", 1, "on", "ON", "On"}
        if value in truthy_values:
            return queryset.filter(filter__iexact="PASS")
        return queryset

    def filter_info_float(self, queryset, name, value):
        if value in (None, ""):
            return queryset

        filter_instance = self.filters.get(name)
        field_name = getattr(filter_instance, "field_name", name)
        annotation = f"{name}_cast"
        return queryset.annotate(
            **{annotation: Cast(F(field_name), FloatField())}
        ).filter(**{f"{annotation}__gte": value})

    def filter_info_int(self, queryset, name, value):
        if value in (None, ""):
            return queryset

        filter_instance = self.filters.get(name)
        field_name = getattr(filter_instance, "field_name", name)
        annotation = f"{name}_cast"
        return queryset.annotate(
            **{annotation: Cast(F(field_name), IntegerField())}
        ).filter(**{f"{annotation}__gte": value})

    # ------------------------------------------------------------------
    # Presentation helpers
    # ------------------------------------------------------------------
    def _apply_bootstrap_widgets(self) -> None:
        form = self.form
        for name, field in form.fields.items():
            widget = field.widget
            if isinstance(widget, forms.CheckboxInput):
                existing = widget.attrs.get("class", "")
                widget.attrs["class"] = self._merge_css_classes(
                    existing, "form-check-input"
                )
            else:
                existing = widget.attrs.get("class", "")
                widget.attrs["class"] = self._merge_css_classes(
                    existing, "form-control"
                )

    @staticmethod
    def _merge_css_classes(existing: str, new_class: str) -> str:
        classes = [cls for cls in existing.split() if cls]
        if new_class and new_class not in classes:
            classes.append(new_class)
        return " ".join(classes)


class SampleGroupFilter(django_filters.FilterSet):
    query = django_filters.CharFilter(method='universal_search', label='Search')

    class Meta:
        model = SampleGroup
        fields = []

    def universal_search(self, queryset, name, value):
        if not value:
            return queryset

        value = value.strip()
        if not value:
            return queryset

        search_filter = (
            Q(name__icontains=value)
            | Q(doi__icontains=value)
            | Q(source_lab__icontains=value)
            | Q(contact_email__icontains=value)
            | Q(contact_phone__icontains=value)
            | Q(inclusion_criteria__icontains=value)
            | Q(exclusion_criteria__icontains=value)
            | Q(comments__icontains=value)
            | Q(created_by__organization_name__icontains=value)
            | Q(reference_genome_build__build_name__icontains=value)
            | Q(reference_genome_build__build_version__icontains=value)
            | Q(genome_complexity__size__icontains=value)
            | Q(genome_complexity__ploidy__icontains=value)
            | Q(genome_complexity__gc_content__icontains=value)
            | Q(sample_origin__tissue__icontains=value)
            | Q(sample_origin__collection_method__icontains=value)
            | Q(sample_origin__storage_conditions__icontains=value)
            | Q(sample_origin__time_stored__icontains=value)
            | Q(material_type__material_type__icontains=value)
            | Q(material_type__integrity_number__icontains=value)
            | Q(library_construction__kit__icontains=value)
            | Q(library_construction__fragmentation__icontains=value)
            | Q(library_construction__adapter_ligation_efficiency__icontains=value)
            | Q(illumina_seq__instrument__icontains=value)
            | Q(illumina_seq__flow_cell__icontains=value)
            | Q(illumina_seq__channel_method__icontains=value)
            | Q(illumina_seq__cluster_density__icontains=value)
            | Q(illumina_seq__qc_software__icontains=value)
            | Q(ont_seq__instrument__icontains=value)
            | Q(ont_seq__flow_cell__icontains=value)
            | Q(ont_seq__flow_cell_version__icontains=value)
            | Q(ont_seq__pore_type__icontains=value)
            | Q(pacbio_seq__instrument__icontains=value)
            | Q(pacbio_seq__flow_cell__icontains=value)
            | Q(pacbio_seq__smrt_cell_type__icontains=value)
            | Q(iontorrent_seq__instrument__icontains=value)
            | Q(iontorrent_seq__flow_cell__icontains=value)
            | Q(iontorrent_seq__chip_type__icontains=value)
            | Q(bioinfo_alignment__tool__icontains=value)
            | Q(bioinfo_alignment__params__icontains=value)
            | Q(bioinfo_alignment__ref_genome_version__icontains=value)
            | Q(bioinfo_alignment__recalibration_settings__icontains=value)
            | Q(bioinfo_variant_calling__tool__icontains=value)
            | Q(bioinfo_variant_calling__version__icontains=value)
            | Q(bioinfo_variant_calling__filtering_thresholds__icontains=value)
            | Q(bioinfo_variant_calling__duplicate_handling__icontains=value)
            | Q(bioinfo_variant_calling__mq__icontains=value)
            | Q(bioinfo_post_proc__normalization__icontains=value)
            | Q(bioinfo_post_proc__harmonization__icontains=value)
            | Q(allele_frequencies__chrom__icontains=value)
            | Q(allele_frequencies__variant_id__icontains=value)
            | Q(allele_frequencies__ref__icontains=value)
            | Q(allele_frequencies__alt__icontains=value)
            | Q(allele_frequencies__filter__icontains=value)
            | Q(allele_frequencies__comments__icontains=value)
        )

        if value.isdigit():
            numeric_value = int(value)
            search_filter |= Q(total_samples=numeric_value)
            search_filter |= Q(library_construction__pcr_cycles=numeric_value)
            search_filter |= Q(allele_frequencies__pos=numeric_value)

        return queryset.filter(search_filter).distinct()
