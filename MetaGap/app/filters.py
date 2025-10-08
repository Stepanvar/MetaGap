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
        if value:
            value = value.strip()
            return queryset.filter(
                Q(name__icontains=value)
                | Q(tissue__icontains=value)
                | Q(collection_method__icontains=value)
                | Q(storage_conditions__icontains=value)
                | Q(comments__icontains=value)
                | Q(sample_origin__tissue__icontains=value)
                | Q(sample_origin__collection_method__icontains=value)
                | Q(sample_origin__storage_conditions__icontains=value)
                | Q(sample_origin__time_stored__icontains=value)
            )
        return queryset
