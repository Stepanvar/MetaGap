# filters.py
from django import forms
import django_filters
from django.db.models import F, Q
from django.db.models.functions import Cast
from django.db.models import FloatField, IntegerField

from .models import AlleleFrequency, SampleGroup


class AlleleFrequencyFilter(django_filters.FilterSet):
    query = django_filters.CharFilter(method="universal_search", label="Search")

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
