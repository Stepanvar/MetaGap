from django import forms
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.contrib.auth.models import User

from .models import OrganizationProfile, SampleGroup

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
    organization_name = forms.CharField(
        required=False,
        label="Organization Name",
        widget=forms.TextInput(
            attrs={"class": "form-control", "placeholder": "Organization"}
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
        organization_name = self.cleaned_data.get("organization_name", "").strip()

        if commit:
            user.save()
            self.save_m2m()
            OrganizationProfile.objects.update_or_create(
                user=user,
                defaults={"organization_name": organization_name or None},
            )
        else:
            # Store the organization name so it can be accessed after saving the instance.
            self._pending_organization_name = organization_name
        return user

    def save_m2m(self):
        super().save_m2m()
        # Ensure that the organization profile is saved if save() was called with commit=False.
        if hasattr(self, "_pending_organization_name") and self.instance.pk:
            organization_name = self._pending_organization_name.strip()
            OrganizationProfile.objects.update_or_create(
                user=self.instance,
                defaults={"organization_name": organization_name or None},
            )
            del self._pending_organization_name

class SampleGroupForm(forms.ModelForm):
    """ModelForm that assigns ``created_by`` using the authenticated user."""

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)

    class Meta:
        model = SampleGroup
        exclude = ["created_by"]  # Exclude fields that will be set in the view
        widgets = {
            "name": forms.TextInput(attrs={"class": "form-control"}),
            # Add widgets for other fields as needed
        }

    def save(self, commit=True):
        sample_group = super().save(commit=False)

        if sample_group.created_by_id is None:
            if self.user is None:
                raise ValueError("SampleGroupForm.save() requires a user when creating a sample group.")
            try:
                organization_profile = self.user.organization_profile
            except OrganizationProfile.DoesNotExist as exc:
                raise ValueError(
                    "SampleGroupForm.save() requires the user to have an organization profile."
                ) from exc
            sample_group.created_by = organization_profile

        if commit:
            sample_group.save()
            self.save_m2m()

        return sample_group

class ImportDataForm(forms.Form):
    data_file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput(attrs={'class': 'form-control'}),
        help_text='Upload a VCF file.'
    )

class EditProfileForm(UserChangeForm):
    password = None  # Exclude password field
    organization_name = forms.CharField(
        required=False,
        label="Organization Name",
        widget=forms.TextInput(
            attrs={"class": "form-control", "placeholder": "Organization"}
        ),
    )

    class Meta:
        model = User
        fields = ("username", "email", "first_name", "last_name")
        widgets = {
            "username": forms.TextInput(attrs={"class": "form-control"}),
            "first_name": forms.TextInput(attrs={"class": "form-control"}),
            "last_name": forms.TextInput(attrs={"class": "form-control"}),
            "email": forms.EmailInput(attrs={"class": "form-control"}),
        }

    field_order = ("username", "email", "first_name", "last_name", "organization_name")

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        organization_profile = getattr(self.instance, "organization_profile", None)
        if organization_profile:
            self.fields["organization_name"].initial = organization_profile.organization_name

    def save(self, commit=True):
        user = super().save(commit=False)
        organization_name = self.cleaned_data.get("organization_name", "").strip()
        if commit:
            user.save()
            self.save_m2m()
            if self.instance.pk:
                OrganizationProfile.objects.update_or_create(
                    user=self.instance,
                    defaults={"organization_name": organization_name or None},
                )
        else:
            self._pending_organization_name = organization_name
        return user

    def save_m2m(self):
        super().save_m2m()
        if hasattr(self, "_pending_organization_name") and self.instance.pk:
            organization_name = self._pending_organization_name.strip()
            OrganizationProfile.objects.update_or_create(
                user=self.instance,
                defaults={"organization_name": organization_name or None},
            )
            del self._pending_organization_name

class DeleteAccountForm(forms.Form):
    confirm = forms.BooleanField(
        required=True,
        label="I confirm that I want to delete my account.",
        widget=forms.CheckboxInput(attrs={"class": "form-check-input"}),
    )
