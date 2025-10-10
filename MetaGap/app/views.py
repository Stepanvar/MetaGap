"""Application views."""

import csv
import logging
import os
from typing import Any, Dict, Iterable, Optional, Tuple
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.files.storage import default_storage
from django.db import models as django_models
from django.http import Http404, HttpResponse
from django.shortcuts import get_object_or_404
from django.urls import reverse_lazy
from django.utils.text import slugify
from django.views.generic import (
    CreateView,
    DeleteView,
    DetailView,
    FormView,
    ListView,
    TemplateView,
    UpdateView,
)
from django_filters.views import FilterView
from django_tables2 import RequestConfig
from django_tables2.views import SingleTableMixin

from .filters import AlleleFrequencySearchFilter, SampleGroupFilter
from .forms import (
    CustomUserCreationForm,
    DeleteAccountForm,
    EditProfileForm,
    ImportDataForm,
    SampleGroupForm,
    SearchForm,
)
from .models import (
    AlleleFrequency,
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    Format,
    GenomeComplexity,
    IlluminaSeq,
    InputQuality,
    IonTorrentSeq,
    Info,
    LibraryConstruction,
    MaterialType,
    OntSeq,
    PacBioSeq,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)
from .services.vcf_importer import VCFImporter
from .tables import create_dynamic_table
from .mixins import OrganizationSampleGroupMixin


logger = logging.getLogger(__name__)

class SampleGroupTableView(ListView):
    """Display a django-tables2 listing of sample groups."""

    model = SampleGroup
    template_name = "sample_group_table.html"
    context_object_name = "sample_groups"
    paginate_by = 10
    ordering = ("name",)

    def get_queryset(self):
        # Optimize query by selecting related metadata and prefetching variants
        queryset = (
            super()
            .get_queryset()
            .select_related(
                "reference_genome_build",
                "genome_complexity",
                "sample_origin",
                "material_type",
                "library_construction",
                "illumina_seq",
                "ont_seq",
                "pacbio_seq",
                "iontorrent_seq",
                "bioinfo_alignment",
                "bioinfo_variant_calling",
                "bioinfo_post_proc",
                "input_quality",
            )
            .prefetch_related("allele_frequencies")
        )
        return queryset.order_by("name", "pk")

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
    template_name = "results.html"
    model = AlleleFrequency
    context_object_name = "allele_frequencies"
    paginate_by = 10
    filterset_class = AlleleFrequencySearchFilter
    ordering = ("chrom", "pos", "ref", "alt", "pk")

    # Dynamically create a table class prioritising variant descriptors first.
    table_class = create_dynamic_table(
        AlleleFrequency,
        table_name="AlleleFrequencyTable",
        include_related=True,
        priority_fields=(
            "chrom",
            "pos",
            "ref",
            "alt",
            "qual",
            "filter",
            "info__af",
            "info__ac",
            "info__an",
            "info__dp",
            "info__mq",
        ),
        exclude_fields=("info__additional", "format__payload"),
    )

    def get_queryset(self):
        base_queryset = (
            AlleleFrequency.objects.select_related(
                "info",
                "format",
                "sample_group",
                "sample_group__reference_genome_build",
                "sample_group__genome_complexity",
                "sample_group__sample_origin",
                "sample_group__material_type",
                "sample_group__library_construction",
                "sample_group__illumina_seq",
                "sample_group__ont_seq",
                "sample_group__pacbio_seq",
                "sample_group__iontorrent_seq",
                "sample_group__bioinfo_alignment",
                "sample_group__bioinfo_variant_calling",
                "sample_group__bioinfo_post_proc",
                "sample_group__input_quality",
                "sample_group__created_by",
            )
        )

        # Instantiate the filter set manually so that we can reuse it in the template context.
        self.filterset = self.filterset_class(
            self.request.GET or None, queryset=base_queryset.order_by(*self.ordering)
        )
        return self.filterset.qs.order_by(*self.ordering).distinct()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault("filter", getattr(self, "filterset", None))
        if hasattr(self, "filterset"):
            context["filter_form"] = self.filterset.form
        return context


class HomePageView(TemplateView):
    """Landing page that exposes the primary search form."""

    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["form"] = SearchForm()
        return context


class ProfileView(LoginRequiredMixin, OrganizationSampleGroupMixin, TemplateView):
    """Display the user's organisation profile and related import tools."""

    template_name = "profile.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        organization_profile = self.get_organization_profile()
        sample_groups = self.get_owned_sample_groups().order_by("name")

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


def _serialize_info(info: Optional[Info]) -> str:
    """Serialise an ``Info`` instance into a flattened key/value string."""

    if not info:
        return ""

    serialized: Dict[str, Any] = {}
    for field in info._meta.concrete_fields:
        if field.name == "id":
            continue
        value = getattr(info, field.name)
        if value in (None, "", {}):
            continue
        if field.name == "additional":
            if isinstance(value, dict):
                for key, additional_value in value.items():
                    if additional_value in (None, "", {}):
                        continue
                    serialized[str(key).upper()] = additional_value
            continue
        serialized[field.name.upper()] = value

    return ";".join(f"{key}={value}" for key, value in serialized.items())


@login_required
def export_sample_group_variants(request, pk: int) -> HttpResponse:
    """Export allele frequency data for a user's sample group."""

    sample_group = get_object_or_404(
        SampleGroup.objects.select_related("created_by"), pk=pk
    )
    organization_profile = getattr(request.user, "organization_profile", None)

    if organization_profile is None or sample_group.created_by != organization_profile:
        raise Http404("Sample group not found.")

    export_format = request.GET.get("format", "csv").lower()
    if export_format == "tsv":
        delimiter = "\t"
        file_extension = "tsv"
        content_type = "text/tab-separated-values"
    else:
        delimiter = ","
        file_extension = "csv"
        content_type = "text/csv"

    filename = slugify(sample_group.name) or f"sample-group-{sample_group.pk}"
    response = HttpResponse(content_type=content_type)
    response["Content-Disposition"] = (
        f"attachment; filename=\"{filename}.{file_extension}\""
    )

    allele_frequencies_qs = sample_group.allele_frequencies.select_related("info").order_by(
        "chrom", "pos", "pk"
    )
    allele_frequencies = list(allele_frequencies_qs)

    info_field_columns = [
        f"info_{field.name.lower()}"
        for field in Info._meta.concrete_fields
        if field.name not in {"id", "additional"}
    ]

    additional_info_keys = set()
    for allele in allele_frequencies:
        info = allele.info
        if info and isinstance(info.additional, dict):
            additional_info_keys.update(str(key) for key in info.additional.keys())

    additional_info_columns = [
        f"info_{key.lower()}" for key in sorted(additional_info_keys)
    ]

    fieldnames = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "variant_id",
        *info_field_columns,
        *additional_info_columns,
    ]

    writer = csv.DictWriter(response, fieldnames=fieldnames, delimiter=delimiter)
    writer.writeheader()

    for allele in allele_frequencies:
        info = allele.info
        row: Dict[str, Any] = {
            "chrom": allele.chrom,
            "pos": allele.pos,
            "ref": allele.ref,
            "alt": allele.alt,
            "variant_id": allele.variant_id or "",
        }

        if info:
            for field in Info._meta.concrete_fields:
                if field.name in {"id", "additional"}:
                    continue
                column = f"info_{field.name.lower()}"
                value = getattr(info, field.name)
                if value not in (None, ""):
                    row[column] = value

            if isinstance(info.additional, dict):
                for key, value in info.additional.items():
                    column = f"info_{str(key).lower()}"
                    if value not in (None, ""):
                        row[column] = value

        for column in fieldnames:
            row.setdefault(column, "")

        writer.writerow(row)

    return response


class DashboardView(LoginRequiredMixin, TemplateView):
    """Authenticated dashboard summarising recent datasets and activity."""

    template_name = "dashboard.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        organization_profile = getattr(self.request.user, "organization_profile", None)

        dataset_queryset = (
            SampleGroup.objects.select_related("created_by", "created_by__user")
            .order_by("-pk")
        )
        action_queryset = (
            AlleleFrequency.objects.select_related(
                "sample_group",
                "sample_group__created_by",
                "sample_group__created_by__user",
            )
            .order_by("-pk")
        )

        if organization_profile is not None:
            dataset_queryset = dataset_queryset.filter(created_by=organization_profile)
            action_queryset = action_queryset.filter(
                sample_group__created_by=organization_profile
            )

        context.update(
            {
                "recent_datasets": list(dataset_queryset[:6]),
                "recent_actions": list(action_queryset[:6]),
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

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault(
            "form_extra_actions",
            [
                {
                    "url": reverse_lazy("profile"),
                    "class": "btn btn-secondary",
                    "label": "Cancel",
                    "text": "Cancel",
                }
            ],
        )
        return context

class UserRegistrationView(CreateView):
    form_class = CustomUserCreationForm
    template_name = "signup.html"
    success_url = reverse_lazy("login")


class ContactView(TemplateView):
    template_name = "contact.html"


class AboutView(TemplateView):
    template_name = "about.html"


class SampleGroupUpdateView(
    LoginRequiredMixin, OrganizationSampleGroupMixin, UpdateView
):
    """Allow organisation members to edit their sample group metadata."""

    model = SampleGroup
    form_class = SampleGroupForm
    template_name = "sample_group_form.html"
    success_url = reverse_lazy("profile")

    def get_queryset(self):
        """Restrict editing to groups owned by the user's organisation."""

        return self.get_owned_sample_groups()

    def get_form_kwargs(self) -> Dict[str, Any]:
        kwargs = super().get_form_kwargs()
        kwargs["user"] = self.request.user
        return kwargs

    def form_valid(self, form):
        messages.success(self.request, "Sample group metadata updated successfully.")
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault(
            "form_extra_actions",
            [
                {
                    "url": reverse_lazy("profile"),
                    "class": "btn btn-outline-secondary",
                    "label": "Cancel",
                    "text": "Cancel",
                }
            ],
        )
        return context


class SampleGroupDeleteView(
    LoginRequiredMixin, OrganizationSampleGroupMixin, DeleteView
):
    """Provide a confirmation flow for removing imported sample groups."""

    model = SampleGroup
    template_name = "sample_group_confirm_delete.html"
    success_url = reverse_lazy("profile")

    def get_queryset(self):
        return self.get_owned_sample_groups()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault(
            "cancel_url",
            reverse_lazy("sample_group_detail", args=[self.object.pk]),
        )
        return context

    def delete(self, request, *args, **kwargs):
        self.object = self.get_object()
        group_name = self.object.name
        response = super().delete(request, *args, **kwargs)
        messages.success(request, f"Deleted {group_name} successfully.")
        return response


class SampleGroupDetailView(
    LoginRequiredMixin, OrganizationSampleGroupMixin, DetailView
):
    """Display an individual sample group's metadata and variant catalogue."""

    model = SampleGroup
    template_name = "sample_group_detail.html"
    context_object_name = "sample_group"

    def get_queryset(self):
        """Limit access to groups owned by the requesting organisation."""

        return (
            self.get_owned_sample_groups()
            .select_related(
                "created_by",
                "created_by__user",
                "reference_genome_build",
                "genome_complexity",
                "sample_origin",
                "material_type",
                "library_construction",
                "illumina_seq",
                "ont_seq",
                "pacbio_seq",
                "iontorrent_seq",
                "bioinfo_alignment",
                "bioinfo_variant_calling",
                "bioinfo_post_proc",
                "input_quality",
            )
            .prefetch_related(
                "allele_frequencies",
                "allele_frequencies__info",
                "allele_frequencies__format",
            )
        )

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        sample_group_field_names = [
            field.name
            for field in SampleGroup._meta.get_fields()
            if isinstance(field, django_models.Field) and not field.auto_created
        ]
        exclude_columns = [
            "sample_group",
            "sample_group__id",
            *[f"sample_group__{name}" for name in sample_group_field_names],
            "info__id",
            "info__additional",
            "format__id",
            "format__payload",
        ]
        priority_columns = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "qual",
            "filter",
            "info__af",
            "info__ac",
            "info__an",
            "info__dp",
            "info__mq",
            "variant_id",
            "format__genotype",
        ]
        table_class = create_dynamic_table(
            AlleleFrequency,
            table_name="AlleleFrequencyTable",
            include_related=True,
            priority_fields=priority_columns,
            exclude_fields=exclude_columns,
        )
        allele_qs = self.object.allele_frequencies.all()
        table = table_class(allele_qs)
        RequestConfig(self.request, paginate={"per_page": 25}).configure(table)
        sample_group = self.object

        def build_section(title: str, items: Iterable[Tuple[str, Any, Optional[str]]]):
            section: Dict[str, Any] = {"title": title, "items": []}
            for key, value, item_type in items:
                label = key
                value_key = slugify(key).replace("-", "_")
                section[value_key] = value
                item_context = {"label": label, "value": value}
                if item_type:
                    item_context["type"] = item_type
                section["items"].append(item_context)
            return section

        metadata_sections = [
            build_section(
                "Summary",
                [
                    ("Name", sample_group.name, None),
                    ("DOI", sample_group.doi, None),
                    ("Source lab", sample_group.source_lab, None),
                    ("Total samples", sample_group.total_samples, None),
                    ("Contact email", sample_group.contact_email, "email"),
                    ("Contact phone", sample_group.contact_phone, None),
                    (
                        "Created by",
                        getattr(sample_group.created_by, "organization_name", None)
                        or sample_group.created_by,
                        None,
                    ),
                ],
            ),
            build_section(
                "Criteria",
                [
                    ("Inclusion", sample_group.inclusion_criteria, None),
                    ("Exclusion", sample_group.exclusion_criteria, None),
                    ("Comments", sample_group.comments, None),
                ],
            ),
            build_section(
                "Sample & Material",
                [
                    (
                        "Reference genome build",
                        getattr(
                            sample_group.reference_genome_build,
                            "build_name",
                            sample_group.reference_genome_build,
                        ),
                        None,
                    ),
                    (
                        "Genome complexity",
                        getattr(
                            sample_group.genome_complexity,
                            "size",
                            sample_group.genome_complexity,
                        ),
                        None,
                    ),
                    ("Sample origin", sample_group.sample_origin, None),
                    ("Material type", sample_group.material_type, None),
                    ("Library construction", sample_group.library_construction, None),
                ],
            ),
            build_section(
                "Sequencing & Bioinformatics",
                [
                    ("Illumina", sample_group.illumina_seq, None),
                    ("Oxford Nanopore", sample_group.ont_seq, None),
                    ("PacBio", sample_group.pacbio_seq, None),
                    ("Ion Torrent", sample_group.iontorrent_seq, None),
                    ("Alignment", sample_group.bioinfo_alignment, None),
                    ("Variant calling", sample_group.bioinfo_variant_calling, None),
                    ("Post-processing", sample_group.bioinfo_post_proc, None),
                ],
            ),
            build_section(
                "Input Quality",
                [
                    ("A260/A280", getattr(sample_group.input_quality, "a260_a280", None), None),
                    ("A260/A230", getattr(sample_group.input_quality, "a260_a230", None), None),
                    (
                        "DNA concentration",
                        getattr(sample_group.input_quality, "dna_concentration", None),
                        None,
                    ),
                    (
                        "RNA concentration",
                        getattr(sample_group.input_quality, "rna_concentration", None),
                        None,
                    ),
                    ("Notes", getattr(sample_group.input_quality, "notes", None), None),
                ],
            ),
        ]

        additional_metadata = sample_group.additional_metadata
        if isinstance(additional_metadata, dict) and additional_metadata:
            metadata_sections.append(
                build_section(
                    "Custom metadata",
                    [
                        (key, value, None)
                        for key, value in additional_metadata.items()
                    ],
                )
            )

        context["metadata_sections"] = metadata_sections
        context["allele_frequency_table"] = table
        context["variant_table"] = table
        return context


class ImportDataView(LoginRequiredMixin, OrganizationSampleGroupMixin, FormView):
    """Handle ingestion of VCF uploads into the relational schema."""

    template_name = "import_data.html"
    form_class = ImportDataForm
    success_url = reverse_lazy("profile")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        sample_groups = self.get_owned_sample_groups().order_by("name")
        context.setdefault("sample_groups", sample_groups)
        return context

    def form_valid(self, form):
        data_file = form.cleaned_data["data_file"]
        temp_path = default_storage.save(f"tmp/{data_file.name}", data_file)
        full_path = os.path.join(settings.MEDIA_ROOT, temp_path)

        importer = VCFImporter(self.request.user)
        try:
            created_group = importer.import_file(full_path)
        except Exception as exc:  # pragma: no cover - defensive feedback channel
            messages.error(self.request, f"An error occurred: {exc}")
        else:
            messages.success(
                self.request,
                f"Imported {created_group.name} successfully.",
            )
            for warning in importer.warnings:
                messages.warning(self.request, warning)
        finally:
            default_storage.delete(temp_path)

        return super().form_valid(form)

