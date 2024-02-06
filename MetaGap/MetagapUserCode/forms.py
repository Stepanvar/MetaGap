from django import forms
from .models import PopulationFrequency, Genotype, Common


class PopulationFrequencyForm(forms.ModelForm):
    class Meta:
        model = PopulationFrequency
        fields = "__all__"


class GenotypeForm(forms.ModelForm):
    class Meta:
        model = Genotype
        fields = "__all__"


class CommonForm(forms.ModelForm):
    class Meta:
        model = Common
        fields = "__all__"
