from django import forms
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.contrib.auth.models import User

from .models import OrganizationProfile, SampleGroup, SampleOrigin


class BootstrapFormMixin:
    """Apply Bootstrap-friendly defaults to Django form widgets."""

    bootstrap_input_class = "form-control"
    bootstrap_select_class = "form-select"
    bootstrap_check_class = "form-check-input"
    bootstrap_file_class = "form-control"

    text_input_widgets = (
        forms.TextInput,
        forms.EmailInput,
        forms.NumberInput,
        forms.URLInput,
        forms.PasswordInput,
        forms.Textarea,
        forms.DateInput,
        forms.DateTimeInput,
        forms.TimeInput,
    )
    select_widgets = (forms.Select, forms.SelectMultiple)
    checkbox_widgets = (forms.CheckboxInput, forms.CheckboxSelectMultiple)
    file_widgets = (forms.FileInput, forms.ClearableFileInput)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._apply_bootstrap_metadata()

    def _apply_bootstrap_metadata(self) -> None:
        for name, field in self.fields.items():
            widget = field.widget
            if isinstance(widget, self.checkbox_widgets):
                default_class = self.bootstrap_check_class
            elif isinstance(widget, self.select_widgets):
                default_class = self.bootstrap_select_class
            elif isinstance(widget, self.file_widgets):
                default_class = self.bootstrap_file_class
            elif isinstance(widget, self.text_input_widgets) or getattr(widget, "input_type", None):
                default_class = self.bootstrap_input_class
            else:
                default_class = self.bootstrap_input_class

            existing_classes = widget.attrs.get("class", "")
            widget.attrs["class"] = self._merge_css_classes(existing_classes, default_class)
            widget.attrs.setdefault("data-field-name", name)
            widget.attrs.setdefault("data-label", field.label or "")
            widget.attrs.setdefault("data-required", "true" if field.required else "false")
            widget.attrs.setdefault("data_widget_class_name", widget.__class__.__name__)
            widget.attrs.setdefault(
                "data_widget_input_type",
                getattr(widget, "input_type", "") or "",
            )
            if field.help_text:
                widget.attrs.setdefault("data-help-text", str(field.help_text))

    @staticmethod
    def _merge_css_classes(existing: str, new_class: str) -> str:
        classes = [cls for cls in existing.split() if cls]
        if new_class and new_class not in classes:
            classes.append(new_class)
        return " ".join(classes)


class _OrganizationProfileFormMixin:
    """Shared helpers for deferring :class:`OrganizationProfile` updates."""

    _pending_organization_name_attr = "_pending_organization_name"

    def _store_pending_profile_update(self, organization_name: str) -> None:
        """Defer profile creation until ``save_m2m`` runs."""

        setattr(self, self._pending_organization_name_attr, organization_name)
        original_save_m2m = getattr(self, "save_m2m", None)

        def deferred_save_m2m():
            if callable(original_save_m2m):
                original_save_m2m()
            self._apply_pending_organization_update()

        self.save_m2m = deferred_save_m2m

    def _update_organization_profile(self, organization_name: str) -> None:
        if self.instance.pk:
            OrganizationProfile.objects.update_or_create(
                user=self.instance,
                defaults={"organization_name": (organization_name or None)},
            )

    def _apply_pending_organization_update(self) -> None:
        if hasattr(self, self._pending_organization_name_attr):
            organization_name = getattr(self, self._pending_organization_name_attr, "").strip()
            self._update_organization_profile(organization_name)
            self._clear_pending_organization_name()

    def _clear_pending_organization_name(self) -> None:
        if hasattr(self, self._pending_organization_name_attr):
            delattr(self, self._pending_organization_name_attr)

class SearchForm(BootstrapFormMixin, forms.Form):
    query = forms.CharField(
        required=False,
        label='',
        widget=forms.TextInput(attrs={'class': 'form-control', 'placeholder': 'Search...'})
    )

class CustomUserCreationForm(BootstrapFormMixin, _OrganizationProfileFormMixin, UserCreationForm):
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
            self._update_organization_profile(organization_name)
            self._clear_pending_organization_name()
        else:
            self._store_pending_profile_update(organization_name)
        return user

class SampleGroupForm(BootstrapFormMixin, forms.ModelForm):
    """ModelForm that assigns ``created_by`` using the authenticated user."""

    sample_origin_tissue = forms.CharField(
        required=False,
        label="Tissue",
        widget=forms.TextInput(attrs={"class": "form-control"}),
    )
    sample_origin_collection_method = forms.CharField(
        required=False,
        label="Collection method",
        widget=forms.TextInput(attrs={"class": "form-control"}),
    )
    sample_origin_storage_conditions = forms.CharField(
        required=False,
        label="Storage conditions",
        widget=forms.TextInput(attrs={"class": "form-control"}),
    )
    sample_origin_time_stored = forms.CharField(
        required=False,
        label="Time stored",
        widget=forms.TextInput(attrs={"class": "form-control"}),
    )

    def __init__(self, *args, **kwargs):
        self.user = kwargs.pop("user", None)
        super().__init__(*args, **kwargs)

        origin = getattr(self.instance, "sample_origin", None)
        if origin is not None:
            self.fields["sample_origin_tissue"].initial = origin.tissue
            self.fields["sample_origin_collection_method"].initial = (
                origin.collection_method
            )
            self.fields["sample_origin_storage_conditions"].initial = (
                origin.storage_conditions
            )
            self.fields["sample_origin_time_stored"].initial = origin.time_stored

    class Meta:
        model = SampleGroup
        exclude = ["created_by", "sample_origin"]  # Managed manually via form fields
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

        origin_payload = {
            "tissue": self._normalize_origin_value(
                self.cleaned_data.get("sample_origin_tissue")
            ),
            "collection_method": self._normalize_origin_value(
                self.cleaned_data.get("sample_origin_collection_method")
            ),
            "storage_conditions": self._normalize_origin_value(
                self.cleaned_data.get("sample_origin_storage_conditions")
            ),
            "time_stored": self._normalize_origin_value(
                self.cleaned_data.get("sample_origin_time_stored")
            ),
        }

        existing_origin = sample_group.sample_origin
        has_origin_updates = any(value is not None for value in origin_payload.values())

        origin_instance = None
        if existing_origin is not None:
            origin_instance = existing_origin
        elif has_origin_updates:
            origin_instance = SampleOrigin()

        if origin_instance is not None:
            for field_name, value in origin_payload.items():
                if value is not None or hasattr(origin_instance, field_name):
                    setattr(origin_instance, field_name, value)

            if not commit:
                sample_group.sample_origin = origin_instance

        if commit:
            if origin_instance is not None:
                origin_instance.save()
                sample_group.sample_origin = origin_instance
            sample_group.save()
            self.save_m2m()

        return sample_group

    @staticmethod
    def _normalize_origin_value(value):
        if value is None:
            return None
        if isinstance(value, str):
            stripped = value.strip()
            return stripped or None
        return value

class ImportDataForm(BootstrapFormMixin, forms.Form):
    data_file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput(attrs={'class': 'form-control'}),
        help_text='Upload a VCF file.'
    )

class EditProfileForm(BootstrapFormMixin, _OrganizationProfileFormMixin, UserChangeForm):
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
            self._update_organization_profile(organization_name)
            self._clear_pending_organization_name()
        else:
            self._store_pending_profile_update(organization_name)
        return user

class DeleteAccountForm(BootstrapFormMixin, forms.Form):
    confirm = forms.BooleanField(
        required=True,
        label="I confirm that I want to delete my account.",
        widget=forms.CheckboxInput(attrs={"class": "form-check-input"}),
    )
