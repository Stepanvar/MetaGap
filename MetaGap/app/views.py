"""Application views."""

import os
from typing import Any, Dict, Iterable, Optional, Tuple

import pysam
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.files.storage import default_storage
from django.urls import reverse_lazy
from django.views.generic import CreateView, FormView, ListView, TemplateView
from django_filters.views import FilterView
from django_tables2 import RequestConfig
from django_tables2.views import SingleTableMixin

from .forms import (
    CustomUserCreationForm,
    DeleteAccountForm,
    EditProfileForm,
    ImportDataForm,
    SearchForm,
)
from .models import AlleleFrequency, Format, Info, SampleGroup
from .tables import create_dynamic_table

class SampleGroupTableView(ListView):
    """Display a django-tables2 listing of sample groups."""

    model = SampleGroup
    template_name = "sample_group_table.html"
    context_object_name = "sample_groups"
    paginate_by = 10

    def get_queryset(self):
        return super().get_queryset().select_related("allele_frequency")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        table_class = create_dynamic_table(
            SampleGroup, table_name="SampleGroupTable", include_related=True
        )
        table = table_class(self.get_queryset())
        RequestConfig(self.request, paginate={"per_page": self.paginate_by}).configure(table)
        context["table"] = table
        return context

class SearchResultsView(SingleTableMixin, FilterView):
    template_name = 'results.html'
    model = SampleGroup
    context_object_name = 'sample_groups'
    paginate_by = 10

    # Dynamically create a table class that includes all related fields.
    table_class = create_dynamic_table(SampleGroup, table_name="CombinedTable", include_related=True)

    def get_queryset(self):
        # Adjust the queryset if needed (for example, using select_related/prefetch_related)
        return SampleGroup.objects.all()


class HomePageView(TemplateView):
    """Landing page that exposes the primary search form."""

    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["form"] = SearchForm()
        return context


class ProfileView(LoginRequiredMixin, TemplateView):
    """Display the user's organisation profile and related import tools."""

    template_name = "profile.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        organization_profile = getattr(self.request.user, "organization_profile", None)
        sample_groups = (
            SampleGroup.objects.filter(created_by=organization_profile).order_by("name")
            if organization_profile
            else SampleGroup.objects.none()
        )
        context.update(
            {
                "organization_profile": organization_profile,
                "sample_groups": sample_groups,
                "import_form": ImportDataForm(),
                "import_form_action": reverse_lazy("import_data"),
                "import_form_enctype": "multipart/form-data",
            }
        )
        return context


class EditProfileView(LoginRequiredMixin, FormView):
    form_class = EditProfileForm
    template_name = "edit_profile.html"
    success_url = reverse_lazy("profile")

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["instance"] = self.request.user
        return kwargs

    def form_valid(self, form):
        form.save()
        messages.success(self.request, "Your profile has been updated.")
        return super().form_valid(form)


class DeleteAccountView(LoginRequiredMixin, FormView):
    form_class = DeleteAccountForm
    template_name = "confirm_delete.html"
    success_url = reverse_lazy("home")

    def form_valid(self, form):
        user = self.request.user
        logout(self.request)
        user.delete()
        messages.success(self.request, "Your account has been deleted.")
        return super().form_valid(form)

class UserRegistrationView(CreateView):
    form_class = CustomUserCreationForm
    template_name = "signup.html"
    success_url = reverse_lazy("login")


class ContactView(TemplateView):
    template_name = "contact.html"


class AboutView(TemplateView):
    template_name = "about.html"

class ImportDataView(LoginRequiredMixin, FormView):
    """Handle ingestion of VCF uploads into the relational schema."""

    template_name = "import_data.html"
    form_class = ImportDataForm
    success_url = reverse_lazy("profile")

    INFO_FIELDS = {
        field.name for field in Info._meta.fields if field.name not in {"id", "additional"}
    }
    FORMAT_FIELDS = {
        field.name for field in Format._meta.fields if field.name not in {"id", "additional"}
    }

    def form_valid(self, form):
        data_file = form.cleaned_data["data_file"]
        temp_path = default_storage.save(f"tmp/{data_file.name}", data_file)
        full_path = os.path.join(settings.MEDIA_ROOT, temp_path)

        try:
            created_group = self.parse_vcf_file(full_path)
            messages.success(
                self.request,
                f"Imported {created_group.name} successfully.",
            )
        except Exception as exc:  # pragma: no cover - defensive feedback channel
            messages.error(self.request, f"An error occurred: {exc}")
        finally:
            default_storage.delete(temp_path)

        return super().form_valid(form)

    def parse_vcf_file(self, file_path: str) -> SampleGroup:
        organization_profile = getattr(self.request.user, "organization_profile", None)
        if organization_profile is None:
            raise ValueError("The current user does not have an organisation profile.")

        with pysam.VariantFile(file_path) as vcf_in:
            metadata = self.extract_sample_group_metadata(vcf_in)
            sample_group = SampleGroup.objects.create(
                name=metadata.get("name")
                or os.path.splitext(os.path.basename(file_path))[0],
                created_by=organization_profile,
                comments=metadata.get("description"),
            )

            created_alleles = []
            for record in vcf_in.fetch():
                info_instance = self._create_info_instance(record.info)
                format_instance, format_sample = self._create_format_instance(record.samples)

                allele = AlleleFrequency.objects.create(
                    chrom=record.chrom,
                    pos=record.pos,
                    variant_id=record.id,
                    ref=record.ref,
                    alt=self._serialize_alt(record.alts),
                    qual=record.qual,
                    filter=self._serialize_filter(record.filter),
                    info=info_instance,
                    format=format_instance,
                )

                if format_instance and format_sample:
                    format_instance.additional = {
                        **(format_instance.additional or {}),
                        "sample_id": format_sample,
                    }
                    format_instance.save(update_fields=["additional"])

                created_alleles.append(allele)

        if created_alleles and sample_group.allele_frequency is None:
            sample_group.allele_frequency = created_alleles[0]
            sample_group.save(update_fields=["allele_frequency"])

        return sample_group

    def extract_sample_group_metadata(self, vcf_in: pysam.VariantFile) -> Dict[str, Any]:
        metadata: Dict[str, Any] = {}
        for record in vcf_in.header.records:
            if record.key == "SAMPLE":
                metadata.setdefault("name", record.get("ID"))
                for key, value in record.items():
                    if key == "ID":
                        continue
                    metadata[key.lower()] = value
        return metadata

    def _create_info_instance(self, info: Any) -> Optional[Info]:
        info_dict = dict(info)
        structured: Dict[str, Any] = {}
        additional: Dict[str, Any] = {}

        for key, value in info_dict.items():
            normalized = key.lower()
            target = structured if normalized in self.INFO_FIELDS else additional
            target[normalized] = self._stringify(value)

        if not structured and not additional:
            return None

        return Info.objects.create(**structured, additional=additional or None)

    def _create_format_instance(
        self, samples: Any
    ) -> Tuple[Optional[Format], Optional[str]]:
        if not samples:
            return None, None

        sample_name, sample_data = next(iter(samples.items()))
        structured: Dict[str, Any] = {}
        additional: Dict[str, Any] = {}

        for key in sample_data.keys():
            normalized = key.lower()
            if normalized == "gt":
                serialized = self._serialize_genotype(sample_data, key)
            else:
                serialized = self._stringify(sample_data[key])
            if normalized in self.FORMAT_FIELDS:
                structured[normalized] = serialized
            else:
                additional[normalized] = serialized

        payload = additional or None
        if not structured and payload is None:
            return None, sample_name

        format_instance = Format.objects.create(**structured, additional=payload)
        return format_instance, sample_name

    @staticmethod
    def _serialize_alt(alts: Optional[Iterable[str]]) -> str:
        return ",".join(alts or [])

    @staticmethod
    def _serialize_filter(filter_field: Any) -> Optional[str]:
        if not filter_field:
            return None
        values = filter_field.keys() if hasattr(filter_field, "keys") else filter_field
        serialized = [str(value) for value in values]
        return ",".join(serialized) if serialized else None

    @staticmethod
    def _serialize_genotype(sample_data: Any, key: str) -> Optional[str]:
        value = sample_data[key]
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            separator = "|" if getattr(sample_data, "phased", False) else "/"
            return separator.join(str(item) for item in value)
        return str(value)

    @staticmethod
    def _stringify(value: Any) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            return ",".join(str(item) for item in value)
        return str(value)
