from django.contrib import admin
from .models import SampleGroup, AlleleFrequency, QualityControlMetrics

admin.site.register(SampleGroup)
admin.site.register(AlleleFrequency)
admin.site.register(QualityControlMetrics)
