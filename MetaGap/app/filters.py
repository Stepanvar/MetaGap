# filters.py
import django_filters
from .models import AlleleFrequency, SampleGroup

class AlleleFrequencyFilter(django_filters.FilterSet):
    query = django_filters.CharFilter(method='universal_search', label='Search')

    class Meta:
        model = AlleleFrequency
        fields = []

    def universal_search(self, queryset, name, value):
        if value:
            return queryset.filter(
                Q(chromosome__icontains=value) |
                Q(position__icontains=value) |
                Q(allele_id__icontains=value)
            )
        return queryset

class SampleGroupFilter(django_filters.FilterSet):
    query = django_filters.CharFilter(method='universal_search', label='Search')

    class Meta:
        model = SampleGroup
        fields = []

    def universal_search(self, queryset, name, value):
        if value:
            return queryset.filter(
                Q(name__icontains=value) |
                Q(tissue__icontains=value) |
                Q(collection_method__icontains=value) |
                Q(storage_conditions__icontains=value) |
                Q(time_stored__icontains=value) |
                Q(comments__icontains=value)
            )
        return queryset
