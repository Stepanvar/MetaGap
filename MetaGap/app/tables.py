import django_tables2 as tables
from .models import AlleleFrequency, SampleGroup

class AlleleFrequencyTable(tables.Table):
    class Meta:
        model = AlleleFrequency
        template_name = 'django_tables2/bootstrap5.html'
        fields = ('chrom', 'pos', 'ref', 'alt', 'qual', 'filter')

class SampleGroupTable(tables.Table):
    class Meta:
        model = SampleGroup
        template_name = 'django_tables2/bootstrap5.html'
        fields = ('name', 'doi', 'created_by')
