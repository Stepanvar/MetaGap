from django import forms
from django.contrib.auth.forms import UserChangeForm, UserCreationForm
from django.contrib.auth.models import User
from django.core.exceptions import ValidationError
from django.utils.html import format_html, format_html_join

from django.conf import settings

from .models import (
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    GenomeComplexity,
    IlluminaSeq,
    IonTorrentSeq,
    LibraryConstruction,
    MaterialType,
    OntSeq,
    OrganizationProfile,
    PacBioSeq,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)


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
    """Lightweight search form that mirrors the available filter parameters."""

    advanced_field_names = (
        "chrom",
        "pos_min",
        "pos_max",
        "ref",
        "alt",
        "qual_min",
        "af_min",
        "ac_min",
        "an_min",
        "dp_min",
        "mq_min",
        "sample_group_name",
        "sample_group_source_lab",
        "sample_origin_tissue",
        "variant_calling_tool",
    )

    query = forms.CharField(
        required=False,
        label="Keyword",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Variant ID, coordinate, or keyword",
            }
        ),
    )
    chrom = forms.CharField(
        required=False,
        label="Chromosome",
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g. 1",
            }
        ),
    )
    pos_min = forms.IntegerField(
        required=False,
        label="Position ≥",
        widget=forms.NumberInput(
            attrs={
                "placeholder": "Start position",
                "min": 0,
            }
        ),
    )
    pos_max = forms.IntegerField(
        required=False,
        label="Position ≤",
        widget=forms.NumberInput(
            attrs={
                "placeholder": "End position",
                "min": 0,
            }
        ),
    )
    ref = forms.CharField(
        required=False,
        label="Reference",
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g. A",
            }
        ),
    )
    alt = forms.CharField(
        required=False,
        label="Alternate",
        widget=forms.TextInput(
            attrs={
                "placeholder": "e.g. T",
            }
        ),
    )
    pass_only = forms.BooleanField(
        required=False,
        label="PASS only",
    )
    qual_min = forms.FloatField(
        required=False,
        label="Quality ≥",
        widget=forms.NumberInput(
            attrs={
                "step": "any",
                "min": 0,
            }
        ),
    )
    af_min = forms.FloatField(
        required=False,
        label="AF ≥",
        widget=forms.NumberInput(
            attrs={
                "step": "any",
                "min": 0,
                "max": 1,
            }
        ),
    )
    ac_min = forms.IntegerField(
        required=False,
        label="AC ≥",
        widget=forms.NumberInput(
            attrs={
                "min": 0,
            }
        ),
    )
    an_min = forms.IntegerField(
        required=False,
        label="AN ≥",
        widget=forms.NumberInput(
            attrs={
                "min": 0,
            }
        ),
    )
    dp_min = forms.IntegerField(
        required=False,
        label="DP ≥",
        widget=forms.NumberInput(
            attrs={
                "min": 0,
            }
        ),
    )
    mq_min = forms.FloatField(
        required=False,
        label="MQ ≥",
        widget=forms.NumberInput(
            attrs={
                "step": "any",
                "min": 0,
            }
        ),
    )
    sample_group_name = forms.CharField(
        required=False,
        label="Sample group",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Name contains",
            }
        ),
    )
    sample_group_source_lab = forms.CharField(
        required=False,
        label="Source lab",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Lab contains",
            }
        ),
    )
    sample_origin_tissue = forms.CharField(
        required=False,
        label="Tissue",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Tissue contains",
            }
        ),
    )
    variant_calling_tool = forms.CharField(
        required=False,
        label="Variant caller",
        widget=forms.TextInput(
            attrs={
                "placeholder": "Tool contains",
            }
        ),
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

    def clean_email(self):
        email = (self.cleaned_data.get("email") or "").strip()
        if not email:
            raise ValidationError("Please provide an email address so we can contact you.")
        if User.objects.filter(email__iexact=email).exists():
            raise ValidationError("An account with this email address already exists.")
        return email

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

    class DatalistTextInput(forms.TextInput):
        """Text input that renders an associated ``<datalist>`` element."""

        def __init__(self, *args, datalist=None, datalist_id=None, **kwargs):
            super().__init__(*args, **kwargs)
            self.datalist = [] if datalist is None else list(dict.fromkeys(filter(None, datalist)))
            self.datalist_id = datalist_id

        def render(self, name, value, attrs=None, renderer=None):
            attrs = attrs or {}
            datalist_id = self.datalist_id or attrs.get("list") or f"{attrs.get('id', name)}_options"
            attrs["list"] = datalist_id
            input_html = super().render(name, value, attrs, renderer)
            options_html = format_html_join(
                "",
                "<option value=\"{}\"></option>",
                ((option,) for option in self.datalist),
            )
            datalist_html = format_html("<datalist id=\"{}\">{}</datalist>", datalist_id, options_html)
            return format_html("{}{}", input_html, datalist_html)

    class CreatableModelChoiceField(forms.ModelChoiceField):
        """ModelChoiceField that creates related objects when no match exists."""

        def __init__(self, *args, create_field: str, **kwargs):
            self.create_field = create_field
            kwargs.setdefault("to_field_name", create_field)
            super().__init__(*args, **kwargs)

        def prepare_value(self, value):
            if isinstance(value, self.queryset.model):
                return getattr(value, self.create_field, super().prepare_value(value))
            return super().prepare_value(value)

        def clean(self, value):
            if isinstance(value, str):
                value = value.strip()

            try:
                return super().clean(value)
            except ValidationError:
                if not value:
                    raise
                if isinstance(value, str):
                    existing = self.queryset.filter(**{self.create_field: value}).first()
                    if existing is not None:
                        return existing
                    return self.queryset.model.objects.create(**{self.create_field: value})
                raise

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

        self._enhance_creatable_metadata_fields()

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

        self._apply_bootstrap_metadata()

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

    def clean_total_samples(self):
        total_samples = self.cleaned_data.get("total_samples")
        if total_samples is None:
            return total_samples
        if total_samples <= 0:
            raise ValidationError("Total samples must be a positive integer.")
        return total_samples

    def clean_contact_phone(self):
        phone = self.cleaned_data.get("contact_phone")
        if not phone:
            return phone

        cleaned = phone.strip()
        if not cleaned:
            return None

        allowed_chars = set("0123456789+-.() extEXT ")
        if any(char not in allowed_chars for char in cleaned):
            raise ValidationError(
                "Use digits and common separators (spaces, dashes, parentheses) for the phone number."
            )
        if sum(char.isdigit() for char in cleaned) < 6:
            raise ValidationError("Enter a phone number with at least six digits.")
        return cleaned

    CREATABLE_METADATA_FIELDS = {
        "reference_genome_build": (ReferenceGenomeBuild, "build_name"),
        "genome_complexity": (GenomeComplexity, "size"),
        "material_type": (MaterialType, "material_type"),
        "library_construction": (LibraryConstruction, "kit"),
        "illumina_seq": (IlluminaSeq, "instrument"),
        "ont_seq": (OntSeq, "instrument"),
        "pacbio_seq": (PacBioSeq, "instrument"),
        "iontorrent_seq": (IonTorrentSeq, "instrument"),
        "bioinfo_alignment": (BioinfoAlignment, "tool"),
        "bioinfo_variant_calling": (BioinfoVariantCalling, "tool"),
        "bioinfo_post_proc": (BioinfoPostProc, "normalization"),
    }

    def _enhance_creatable_metadata_fields(self) -> None:
        for field_name, (model, lookup_field) in self.CREATABLE_METADATA_FIELDS.items():
            if field_name not in self.fields:
                continue

            original_field = self.fields[field_name]
            queryset = getattr(original_field, "queryset", model.objects.all())
            datalist_options = self._build_metadata_datalist(field_name, queryset, lookup_field)
            widget = self.DatalistTextInput(
                datalist=datalist_options,
                attrs={"class": self.bootstrap_input_class, "placeholder": "Enter or select a value"},
            )

            help_text = original_field.help_text or ""
            helper_suffix = "Start typing to select an existing value or enter a new one."
            help_text = f"{help_text} {helper_suffix}".strip()

            self.fields[field_name] = self.CreatableModelChoiceField(
                queryset=queryset,
                create_field=lookup_field,
                required=original_field.required,
                empty_label=getattr(original_field, "empty_label", None),
                label=original_field.label,
                widget=widget,
                help_text=help_text,
            )

            if self.instance.pk:
                related = getattr(self.instance, field_name, None)
                if related is not None:
                    initial_value = getattr(related, lookup_field, None)
                    if initial_value is not None:
                        self.fields[field_name].initial = initial_value

    def _build_metadata_datalist(self, field_name, queryset, lookup_field):
        values = []
        for obj in queryset:
            value = getattr(obj, lookup_field, None)
            if value:
                values.append(value)

        if field_name == "material_type":
            values.extend(choice for choice, _ in getattr(MaterialType, "MATERIAL_CHOICES", []))

        seen = {}
        return [seen.setdefault(val, val) for val in values if val and val not in seen]

class ImportDataForm(BootstrapFormMixin, forms.Form):
    """Validate uploads for the variant import workflow."""

    DEFAULT_MAX_SIZE_MB = 50
    SUPPORTED_EXTENSIONS = {".vcf", ".vcf.gz", ".vcf.bgz", ".bcf"}
    SUPPORTED_CONTENT_TYPES = {
        "text/plain",
        "text/vcf",
        "application/octet-stream",
        "application/gzip",
        "application/x-gzip",
        "application/bcf",
    }

    data_file = forms.FileField(
        required=True,
        widget=forms.ClearableFileInput(attrs={"class": "form-control"}),
        help_text=(
            "Upload a VCF (.vcf, .vcf.gz) or BCF file containing the variants you want to ingest."
        ),
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        widget = self.fields["data_file"].widget
        widget.attrs.setdefault("accept", ",".join(sorted(self.SUPPORTED_EXTENSIONS)))

        configured_limit = getattr(settings, "METAGAP_MAX_UPLOAD_SIZE_MB", None)
        self.max_upload_size_mb = configured_limit or self.DEFAULT_MAX_SIZE_MB
        self.max_upload_size_bytes = int(self.max_upload_size_mb * 1024 * 1024)
        help_text = (
            "Files up to %(limit)s&nbsp;MB are accepted. Ensure the header includes metadata"
            " sections for your sample group." % {"limit": self.max_upload_size_mb}
        )
        self.fields["data_file"].help_text = format_html("{}<br>{}", self.fields["data_file"].help_text, help_text)

    def clean_data_file(self):
        uploaded = self.cleaned_data.get("data_file")
        if not uploaded:
            raise ValidationError("Select a VCF or BCF file to import.")

        name = uploaded.name or ""
        normalized_name = name.lower()
        is_supported = any(normalized_name.endswith(ext) for ext in self.SUPPORTED_EXTENSIONS)
        if not is_supported:
            raise ValidationError(
                "Unsupported file type. Please upload a file ending in %s." % ", ".join(sorted(self.SUPPORTED_EXTENSIONS))
            )

        file_size = getattr(uploaded, "size", None)
        if file_size is None:
            file_size = uploaded.file.size if hasattr(uploaded, "file") else None
        if file_size is not None and file_size > self.max_upload_size_bytes:
            raise ValidationError(
                f"Files larger than {self.max_upload_size_mb} MB cannot be imported."
            )

        content_type = getattr(uploaded, "content_type", "")
        if content_type and content_type not in self.SUPPORTED_CONTENT_TYPES:
            raise ValidationError(
                "The uploaded file does not look like a VCF/BCF file. Please check the format and try again."
            )

        return uploaded

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

    def clean_email(self):
        email = (self.cleaned_data.get("email") or "").strip()
        if not email:
            raise ValidationError("Please provide a valid email address.")
        queryset = User.objects.filter(email__iexact=email)
        if self.instance.pk:
            queryset = queryset.exclude(pk=self.instance.pk)
        if queryset.exists():
            raise ValidationError("Another user is already using this email address.")
        return email

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
