# app/forms.py

from django import forms
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.contrib.auth.models import User
from .models import SampleGroup

class SearchForm(forms.Form):
    query = forms.CharField(
        required=False,
        label='',
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Search...'})
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

class SampleGroupForm(forms.ModelForm):
    class Meta:
        model = SampleGroup
        exclude = ['created_by']  # Exclude fields that will be set in the view
        widgets = {
            'name': forms.TextInput(attrs={'class': 'form-control'}),
            # Add widgets for other fields as needed
        }

class ImportDataForm(forms.Form):
    data_file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput(attrs={'class': 'form-control'}),
        help_text='Upload a VCF file.'
    )

class EditProfileForm(UserChangeForm):
    password = None  # Exclude password field

    class Meta:
        model = User
        fields = ('username', 'email', 'first_name', 'last_name')
        widgets = {
            'username': forms.TextInput(attrs={'class': 'form-control'}),
            # Add widgets for other fields
        }

class DeleteAccountForm(forms.Form):
    confirm = forms.BooleanField(
        required=True,
        label="I confirm that I want to delete my account.",
    )

class ImportDataForm(forms.Form):
    data_file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput(attrs={'class': 'form-control'}),
        help_text='Upload a VCF file.'
    )