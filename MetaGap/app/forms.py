from django import forms
from django.contrib.auth.forms import UserCreationForm
from django.contrib.auth.models import User
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


class SearchForm(forms.Form):
	query = forms.CharField(
		required=False,
		label="",
		widget=forms.TextInput(
			attrs={"class": "form-control", "placeholder": "Search"}
		),
	)


class CustomUserCreationForm(UserCreationForm):
	email = forms.EmailField(
		required=True,
		widget=forms.EmailInput(
			attrs={"class": "form-control", "placeholder": "Email Address"}
		),
	)
	username = forms.CharField(
		required=True,
		widget=forms.TextInput(
			attrs={"class": "form-control", "placeholder": "Username"}
		),
	)
	password1 = forms.CharField(
		label="Password",
		widget=forms.PasswordInput(
			attrs={"class": "form-control", "placeholder": "Password"}
		),
	)
	password2 = forms.CharField(
		label="Confirm Password",
		widget=forms.PasswordInput(
			attrs={"class": "form-control", "placeholder": "Confirm Password"}
		),
	)

	class Meta:
		model = User
		fields = ("username", "email", "password1", "password2")

	def save(self, commit=True):
		user = super().save(commit=False)
		user.email = self.cleaned_data["email"]
		if commit:
			user.save()
		return user
