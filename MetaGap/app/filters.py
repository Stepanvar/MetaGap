# filters.py
import django_filters
from django.db.models import Q
from .models import AlleleFrequency, SampleGroup


class AlleleFrequencyFilter(django_filters.FilterSet):
    query = django_filters.CharFilter(method='universal_search', label='Search')

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
            | Q(tissue__icontains=value)
            | Q(collection_method__icontains=value)
            | Q(storage_conditions__icontains=value)
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
            | Q(platform_independent__pooling__icontains=value)
            | Q(platform_independent__instrument__icontains=value)
            | Q(platform_independent__sequencing_kit__icontains=value)
            | Q(platform_independent__base_calling_alg__icontains=value)
            | Q(platform_independent__q30__icontains=value)
            | Q(platform_independent__normalized_coverage__icontains=value)
            | Q(platform_independent__run_specific_calibration__icontains=value)
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
