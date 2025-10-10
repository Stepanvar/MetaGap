"""Views related to sample group management and display."""

from __future__ import annotations

import csv
from typing import Any, Dict, Iterable, Optional, Tuple

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.db import models as django_models
from django.http import Http404, HttpResponse
from django.shortcuts import get_object_or_404
from django.urls import reverse_lazy
from django.utils.text import slugify
from django.views.generic import DeleteView, DetailView, ListView, UpdateView
from django_tables2 import RequestConfig

from ..forms import SampleGroupForm
from ..mixins import OrganizationSampleGroupMixin
from ..models import Info, SampleGroup
from ..tables import build_allele_frequency_table, create_dynamic_table


class SampleGroupTableView(ListView):
    """Display a django-tables2 listing of sample groups."""

    model = SampleGroup
    template_name = "sample_group_table.html"
    context_object_name = "sample_groups"
    paginate_by = 10
    ordering = ("name",)

    def get_queryset(self):
        """Optimise the queryset for related metadata and variant prefetching."""

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
            SampleGroup,
            table_name="SampleGroupTable",
            include_related=True,
        )
        table = table_class(self.get_queryset())
        RequestConfig(self.request, paginate={"per_page": self.paginate_by}).configure(table)
        context["table"] = table
        return context


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
            "format__id",
        ]
        table_class = build_allele_frequency_table(
            priority_extra=("variant_id", "format__genotype"),
            exclude_extra=exclude_columns,
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
                    ("Created", sample_group.created_at, "datetime"),
                    ("Last updated", sample_group.updated_at, "datetime"),
                    ("Comments", sample_group.comments, None),
                ],
            ),
            build_section(
                "Reference & complexity",
                [
                    (
                        "Reference genome build",
                        getattr(sample_group.reference_genome_build, "build_name", None),
                        None,
                    ),
                    (
                        "Reference genome version",
                        getattr(sample_group.reference_genome_build, "build_version", None),
                        None,
                    ),
                    (
                        "Genome size",
                        getattr(sample_group.genome_complexity, "size", None),
                        None,
                    ),
                    (
                        "Ploidy",
                        getattr(sample_group.genome_complexity, "ploidy", None),
                        None,
                    ),
                    (
                        "GC content",
                        getattr(sample_group.genome_complexity, "gc_content", None),
                        None,
                    ),
                ],
            ),
            build_section(
                "Sample origin",
                [
                    (
                        "Tissue",
                        getattr(sample_group.sample_origin, "tissue", None),
                        None,
                    ),
                    (
                        "Collection method",
                        getattr(sample_group.sample_origin, "collection_method", None),
                        None,
                    ),
                    (
                        "Storage conditions",
                        getattr(sample_group.sample_origin, "storage_conditions", None),
                        None,
                    ),
                    (
                        "Time stored",
                        getattr(sample_group.sample_origin, "time_stored", None),
                        None,
                    ),
                ],
            ),
            build_section(
                "Material type",
                [
                    (
                        "Material",
                        getattr(sample_group.material_type, "material_type", None),
                        None,
                    ),
                    (
                        "Integrity",
                        getattr(sample_group.material_type, "integrity_number", None),
                        None,
                    ),
                ],
            ),
            build_section(
                "Library construction",
                [
                    ("Kit", getattr(sample_group.library_construction, "kit", None), None),
                    (
                        "Fragmentation",
                        getattr(sample_group.library_construction, "fragmentation", None),
                        None,
                    ),
                    (
                        "Adapter ligation efficiency",
                        getattr(
                            sample_group.library_construction,
                            "adapter_ligation_efficiency",
                            None,
                        ),
                        None,
                    ),
                    (
                        "PCR cycles",
                        getattr(sample_group.library_construction, "pcr_cycles", None),
                        None,
                    ),
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
                    (
                        "A260/A280",
                        getattr(sample_group.input_quality, "a260_a280", None),
                        None,
                    ),
                    (
                        "A260/A230",
                        getattr(sample_group.input_quality, "a260_a230", None),
                        None,
                    ),
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
                    (
                        "Notes",
                        getattr(sample_group.input_quality, "notes", None),
                        None,
                    ),
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

__all__ = [
    "SampleGroupTableView",
    "SampleGroupDetailView",
    "SampleGroupUpdateView",
    "SampleGroupDeleteView",
    "export_sample_group_variants",
]
