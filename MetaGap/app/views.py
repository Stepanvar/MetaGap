"""Application views."""

import csv
import json
import os
from typing import Any, Dict, Iterable, Optional, Tuple

import pysam
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
from django.views.generic import CreateView, FormView, ListView, TemplateView, UpdateView
from django_filters.views import FilterView
from django_tables2 import RequestConfig
from django_tables2.views import SingleTableMixin

from .filters import SampleGroupFilter
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
    PlatformIndependent,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)
from .tables import create_dynamic_table

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
                "platform_independent",
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
    model = SampleGroup
    context_object_name = "sample_groups"
    paginate_by = 10
    filterset_class = SampleGroupFilter
    ordering = ("name", "pk")

    # Dynamically create a table class that includes all related fields.
    table_class = create_dynamic_table(
        SampleGroup, table_name="CombinedTable", include_related=True
    )

    def get_queryset(self):
        base_queryset = SampleGroup.objects.select_related(
            "reference_genome_build",
            "genome_complexity",
            "sample_origin",
            "material_type",
            "library_construction",
            "illumina_seq",
            "ont_seq",
            "pacbio_seq",
            "iontorrent_seq",
            "platform_independent",
            "bioinfo_alignment",
            "bioinfo_variant_calling",
            "bioinfo_post_proc",
            "created_by",
        ).prefetch_related("allele_frequencies").order_by(*self.ordering)

        # Instantiate the filter set manually so that we can reuse it in the template context.
        self.filterset = self.filterset_class(
            self.request.GET or None, queryset=base_queryset
        )
        return self.filterset.qs.order_by(*self.ordering).distinct()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault("filter", getattr(self, "filterset", None))
        context["form"] = SearchForm(self.request.GET or None)
        return context


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
        user = self.request.user
        organization_profile = getattr(user, "organization_profile", None)

        if organization_profile is None:
            sample_groups = SampleGroup.objects.none()
        else:
            sample_groups = SampleGroup.objects.filter(
                created_by=organization_profile
            ).order_by("name")

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

    writer = csv.writer(response, delimiter=delimiter)
    writer.writerow(["Chrom", "Pos", "Ref", "Alt", "AF", "INFO"])

    allele_frequencies = sample_group.allele_frequencies.select_related("info").order_by(
        "chrom", "pos", "pk"
    )

    for allele in allele_frequencies:
        info = allele.info
        writer.writerow(
            [
                allele.chrom,
                allele.pos,
                allele.ref,
                allele.alt,
                info.af if info and info.af else "",
                _serialize_info(info),
            ]
        )

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


class SampleGroupUpdateView(LoginRequiredMixin, UpdateView):
    """Allow organisation members to edit their sample group metadata."""

    model = SampleGroup
    form_class = SampleGroupForm
    template_name = "sample_group_form.html"
    success_url = reverse_lazy("profile")

    def get_queryset(self):
        """Restrict editing to groups owned by the user's organisation."""

        organization_profile = getattr(self.request.user, "organization_profile", None)
        if organization_profile is None:
            return SampleGroup.objects.none()
        return SampleGroup.objects.filter(created_by=organization_profile)

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
                }
            ],
        )
        return context


class ImportDataView(LoginRequiredMixin, FormView):
    """Handle ingestion of VCF uploads into the relational schema."""

    template_name = "import_data.html"
    form_class = ImportDataForm
    success_url = reverse_lazy("profile")

    METADATA_SECTION_MAP = {
        "SAMPLE_GROUP": "sample_group",
        "SAMPLEGROUP": "sample_group",
        "GROUP": "sample_group",
        "REFERENCE_GENOME_BUILD": "reference_genome_build",
        "REFERENCE_GENOME": "reference_genome_build",
        "REFERENCE": "reference_genome_build",
        "GENOME_COMPLEXITY": "genome_complexity",
        "SAMPLE_ORIGIN": "sample_origin",
        "ORIGIN": "sample_origin",
        "MATERIAL_TYPE": "material_type",
        "LIBRARY_CONSTRUCTION": "library_construction",
        "LIBRARY_PREP": "library_construction",
        "ILLUMINA_SEQ": "illumina_seq",
        "ONT_SEQ": "ont_seq",
        "PACBIO_SEQ": "pacbio_seq",
        "IONTORRENT_SEQ": "iontorrent_seq",
        "ION_TORRENT_SEQ": "iontorrent_seq",
        "PLATFORM_INDEPENDENT": "platform_independent",
        "BIOINFO_ALIGNMENT": "bioinfo_alignment",
        "ALIGNMENT": "bioinfo_alignment",
        "BIOINFO_VARIANT_CALLING": "bioinfo_variant_calling",
        "VARIANT_CALLING": "bioinfo_variant_calling",
        "BIOINFO_POSTPROC": "bioinfo_post_proc",
        "BIOINFO_POST_PROC": "bioinfo_post_proc",
        "BIOINFO_POST_PROCESSING": "bioinfo_post_proc",
        "INPUT_QUALITY": "input_quality",
    }

    METADATA_MODEL_MAP = {
        "reference_genome_build": ReferenceGenomeBuild,
        "genome_complexity": GenomeComplexity,
        "sample_origin": SampleOrigin,
        "material_type": MaterialType,
        "library_construction": LibraryConstruction,
        "illumina_seq": IlluminaSeq,
        "ont_seq": OntSeq,
        "pacbio_seq": PacBioSeq,
        "iontorrent_seq": IonTorrentSeq,
        "platform_independent": PlatformIndependent,
        "bioinfo_alignment": BioinfoAlignment,
        "bioinfo_variant_calling": BioinfoVariantCalling,
        "bioinfo_post_proc": BioinfoPostProc,
        "input_quality": InputQuality,
    }

    METADATA_FIELD_ALIASES = {
        "sample_group": {
            "name": ["group", "group_name", "dataset", "id"],
            "doi": ["dataset_doi", "group_doi"],
            "source_lab": ["lab", "lab_name", "source"],
            "contact_email": ["email", "lab_email", "contact"],
            "contact_phone": ["phone", "lab_phone"],
            "total_samples": ["samples", "sample_count", "n_samples"],
            "inclusion_criteria": ["inclusion", "inclusioncriteria"],
            "exclusion_criteria": ["exclusion", "exclusioncriteria"],
            "comments": ["description", "notes"],
        },
        "input_quality": {
            "a260_a280": ["a260_280", "ratio_a260_a280"],
            "a260_a230": ["a260_230", "ratio_a260_a230"],
            "dna_concentration": ["dna_conc", "dna_concentration_ng_ul", "concentration"],
            "rna_concentration": ["rna_conc", "rna_concentration_ng_ul"],
            "notes": ["note", "comment", "comments"],
            "additional_metrics": ["metrics", "additional_metrics"],
        },
        "reference_genome_build": {
            "build_name": ["name", "reference", "build"],
            "build_version": ["version", "build_version"],
            "additional_info": ["additional", "additional_info"],
        },
        "genome_complexity": {
            "size": ["genome_size", "size_bp"],
            "ploidy": ["ploidy_level"],
            "gc_content": ["gc", "gc_percent"],
        },
        "sample_origin": {
            "tissue": [
                "tissue_type",
                "sample_group_tissue",
                "samplegroup_tissue",
            ],
            "collection_method": [
                "collection",
                "method",
                "sample_group_collection_method",
                "samplegroup_collection_method",
            ],
            "storage_conditions": [
                "storage",
                "storage_conditions",
                "sample_group_storage_conditions",
                "samplegroup_storage_conditions",
            ],
            "time_stored": ["time", "storage_time"],
        },
        "material_type": {
            "material_type": ["type"],
            "integrity_number": ["rin", "din", "integrity"],
        },
        "library_construction": {
            "kit": ["library_kit", "kit_name"],
            "fragmentation": ["fragmentation_method"],
            "adapter_ligation_efficiency": ["adapter_efficiency"],
            "pcr_cycles": ["pcr", "pcr_cycles"],
        },
        "illumina_seq": {
            "instrument": ["machine", "instrument"],
            "flow_cell": ["flowcell", "flow_cell_id"],
            "channel_method": ["channel"],
            "cluster_density": ["cluster"],
            "qc_software": ["qc", "software"],
        },
        "ont_seq": {
            "instrument": ["machine", "instrument"],
            "flow_cell": ["flowcell", "flow_cell_id"],
            "flow_cell_version": ["flowcell_version"],
            "pore_type": ["pore"],
            "bias_voltage": ["bias"],
        },
        "pacbio_seq": {
            "instrument": ["machine", "instrument"],
            "flow_cell": ["flowcell", "flow_cell_id"],
            "smrt_cell_type": ["smrt_cell", "cell_type"],
            "zmw_density": ["zmw"],
        },
        "iontorrent_seq": {
            "instrument": ["machine", "instrument"],
            "flow_cell": ["flowcell", "flow_cell_id"],
            "chip_type": ["chip"],
            "ph_calibration": ["ph"],
            "flow_order": ["floworder"],
            "ion_sphere_metrics": ["ionsphere", "sphere_metrics"],
        },
        "platform_independent": {
            "pooling": ["pool"],
            "sequencing_kit": ["seq_kit", "kit"],
            "base_calling_alg": ["basecalling", "base_calling"],
            "q30": ["q30_rate"],
            "normalized_coverage": ["coverage"],
            "run_specific_calibration": ["calibration"],
        },
        "bioinfo_alignment": {
            "tool": ["aligner", "software", "tool"],
            "params": ["parameters", "params"],
            "ref_genome_version": ["reference_version"],
            "recalibration_settings": ["recalibration", "recal_settings"],
        },
        "bioinfo_variant_calling": {
            "tool": ["caller", "tool"],
            "version": ["tool_version", "version"],
            "filtering_thresholds": ["filters", "thresholds"],
            "duplicate_handling": ["duplicates"],
            "mq": ["mapping_quality"],
        },
        "bioinfo_post_proc": {
            "normalization": ["normalisation", "normalization_method"],
            "harmonization": ["harmonisation", "harmonization_method"],
        },
    }

    SECTION_PRIMARY_FIELD = {
        "reference_genome_build": "build_name",
        "genome_complexity": "size",
        "sample_origin": "tissue",
        "material_type": "material_type",
        "library_construction": "kit",
        "illumina_seq": "instrument",
        "ont_seq": "instrument",
        "pacbio_seq": "instrument",
        "iontorrent_seq": "instrument",
        "platform_independent": "instrument",
        "bioinfo_alignment": "tool",
        "bioinfo_variant_calling": "tool",
        "bioinfo_post_proc": "normalization",
    }

    INFO_FIELD_MAP = {
        "aa": "aa",
        "ac": "ac",
        "af": "af",
        "an": "an",
        "bq": "bq",
        "cigar": "cigar",
        "db": "db",
        "dp": "dp",
        "end": "end",
        "h2": "h2",
        "h3": "h3",
        "mq": "mq",
        "mq0": "mq0",
        "ns": "ns",
        "sb": "sb",
    }

    FORMAT_FIELD_MAP = {
        "ad": "ad",
        "adf": "adf",
        "adr": "adr",
        "dp": "dp",
        "ec": "ec",
        "ft": "ft",
        "gl": "gl",
        "gp": "gp",
        "gq": "gq",
        "gt": "gt",
        "hq": "hq",
        "mq": "mq",
        "pl": "pl",
        "pq": "pq",
        "ps": "ps",
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

        metadata: Dict[str, Any] = {}
        sample_group: Optional[SampleGroup] = None
        try:
            with pysam.VariantFile(file_path) as vcf_in:
                metadata = self.extract_sample_group_metadata(vcf_in)
                sample_group = self._create_sample_group(
                    metadata, file_path, organization_profile
                )

                self._populate_sample_group_from_pysam(vcf_in, sample_group)
        except (OSError, ValueError) as exc:  # pragma: no cover - defensive fallback
            messages.warning(
                self.request,
                f"Could not parse VCF metadata with pysam: {exc}. Falling back to a text parser.",
            )
            if sample_group is not None:
                sample_group.delete()
            metadata = self._extract_metadata_text_fallback(file_path)
            sample_group = self._create_sample_group(
                metadata, file_path, organization_profile
            )
            self._parse_vcf_text_fallback(file_path, sample_group)

        assert sample_group is not None
        return sample_group

    def _populate_sample_group_from_pysam(
        self, vcf_in: pysam.VariantFile, sample_group: SampleGroup
    ) -> None:
        for record in vcf_in.fetch():
            info_instance = self._create_info_instance(record.info)
            format_instance, format_sample = self._create_format_instance(record.samples)

            self._create_allele_frequency(
                sample_group,
                chrom=record.chrom,
                pos=record.pos,
                variant_id=record.id,
                ref=record.ref,
                alt=self._serialize_alt(record.alts),
                qual=record.qual,
                filter_value=self._serialize_filter(record.filter),
                info=info_instance,
                format_instance=format_instance,
                format_sample=format_sample,
            )

        return None

    def _create_sample_group(
        self,
        metadata: Dict[str, Any],
        file_path: str,
        organization_profile: Any,
    ) -> SampleGroup:
        group_data, _, _ = self._extract_section_data(
            metadata, "sample_group", SampleGroup
        )

        fallback_name = os.path.splitext(os.path.basename(file_path))[0]
        name = group_data.pop("name", None) or metadata.get("name") or fallback_name
        comments = (
            group_data.pop("comments", None)
            or metadata.get("comments")
            or metadata.get("description")
        )

        sample_group = SampleGroup.objects.create(
            name=name,
            created_by=organization_profile,
            comments=comments,
            **group_data,
        )

        update_fields: list[str] = []
        for section, model_cls in self.METADATA_MODEL_MAP.items():
            section_data, section_consumed, additional = self._extract_section_data(
                metadata, section, model_cls
            )

            if not section_data and additional is None:
                continue

            payload = {key: value for key, value in section_data.items() if value is not None}
            if additional is not None:
                additional_field = self._resolve_additional_field(model_cls)
                if additional_field:
                    existing = payload.get(additional_field)
                    if isinstance(existing, dict) and isinstance(additional, dict):
                        payload[additional_field] = {**existing, **additional}
                    elif existing is None:
                        payload[additional_field] = additional
                    else:
                        payload[additional_field] = additional

            if not payload:
                continue

            instance = model_cls.objects.create(**payload)
            setattr(sample_group, section, instance)
            update_fields.append(section)

        if update_fields:
            sample_group.save(update_fields=update_fields)

        return sample_group

    def _extract_section_data(
        self,
        metadata: Dict[str, Any],
        section: str,
        model_cls: Any,
    ) -> Tuple[Dict[str, Any], set[str], Optional[Dict[str, Any]]]:
        alias_map = self.METADATA_FIELD_ALIASES.get(section, {})
        section_data: Dict[str, Any] = {}
        consumed: set[str] = set()

        for field_name, aliases in alias_map.items():
            model_field = self._get_model_field(model_cls, field_name)
            if model_field is None:
                continue
            key = self._find_metadata_key(metadata, section, field_name, aliases)
            if key is None:
                continue
            raw_value = metadata[key]
            section_data[field_name] = self._coerce_model_value(model_field, raw_value)
            consumed.add(key)

        primary_field = self.SECTION_PRIMARY_FIELD.get(section)
        if (
            primary_field
            and primary_field not in section_data
            and section in metadata
        ):
            model_field = self._get_model_field(model_cls, primary_field)
            if model_field is not None:
                section_data[primary_field] = self._coerce_model_value(
                    model_field, metadata[section]
                )
                consumed.add(section)

        additional = self._build_additional_payload(metadata, section, consumed)
        return section_data, consumed, additional

    @staticmethod
    def _get_model_field(model_cls: Any, field_name: str) -> Optional[django_models.Field]:
        try:
            return model_cls._meta.get_field(field_name)
        except django_models.FieldDoesNotExist:  # pragma: no cover - defensive
            return None

    def _find_metadata_key(
        self,
        metadata: Dict[str, Any],
        section: str,
        field_name: str,
        aliases: Iterable[str],
    ) -> Optional[str]:
        ordered_aliases = [field_name.lower()]
        for alias in aliases:
            normalized_alias = alias.lower()
            if normalized_alias not in ordered_aliases:
                ordered_aliases.append(normalized_alias)

        normalized_variants: list[str] = []
        for candidate in ordered_aliases:
            for variant in (
                candidate,
                candidate.replace(" ", "_"),
                candidate.replace("-", "_"),
            ):
                if variant not in normalized_variants:
                    normalized_variants.append(variant)

        parts = section.split("_")
        section_variants_set: set[str] = set()

        def build_section_variants(index: int, current: str) -> None:
            current += parts[index]
            if index == len(parts) - 1:
                section_variants_set.add(current)
                return
            build_section_variants(index + 1, current + "_")
            build_section_variants(index + 1, current)

        if parts:
            build_section_variants(0, "")
        else:  # pragma: no cover - defensive
            section_variants_set.add(section)

        section_variants = list(section_variants_set)
        section_variants.sort(key=len)

        prefixes: list[str] = []
        for variant in section_variants:
            for suffix in ("_", "."):
                candidate_prefix = f"{variant}{suffix}"
                if candidate_prefix not in prefixes:
                    prefixes.append(candidate_prefix)
        prefixes.append("")

        for candidate in normalized_variants:
            for prefix in prefixes:
                key = f"{prefix}{candidate}" if prefix else candidate
                if key in metadata:
                    return key
        return None

    def _coerce_model_value(
        self, field: django_models.Field, value: Any
    ) -> Optional[Any]:
        if value is None:
            return None
        if isinstance(value, str):
            stripped = value.strip()
            if stripped == "":
                return None
        else:
            stripped = value

        if isinstance(field, django_models.IntegerField):
            try:
                return int(stripped)
            except (TypeError, ValueError):
                return None
        if isinstance(field, django_models.FloatField):
            try:
                return float(stripped)
            except (TypeError, ValueError):
                return None
        if isinstance(field, django_models.JSONField):
            if isinstance(stripped, str):
                try:
                    return json.loads(stripped)
                except json.JSONDecodeError:
                    pass
            return stripped
        return stripped

    def _build_additional_payload(
        self,
        metadata: Dict[str, Any],
        section: str,
        consumed: Iterable[str],
    ) -> Optional[Dict[str, Any]]:
        consumed_set = set(consumed)
        prefixes = (f"{section}_", f"{section}.")
        additional: Dict[str, Any] = {}

        for key, value in metadata.items():
            if key in consumed_set:
                continue
            for prefix in prefixes:
                if key.startswith(prefix):
                    trimmed = key[len(prefix) :]
                    additional[trimmed] = self._coerce_additional_value(value)
                    break

        return additional or None

    @staticmethod
    def _coerce_additional_value(value: Any) -> Any:
        if isinstance(value, str):
            stripped = value.strip()
            if stripped == "":
                return None
            try:
                return json.loads(stripped)
            except json.JSONDecodeError:
                pass
            try:
                return int(stripped)
            except ValueError:
                try:
                    return float(stripped)
                except ValueError:
                    return stripped
        return value

    @staticmethod
    def _resolve_additional_field(model_cls: Any) -> Optional[str]:
        for candidate in ("additional_info", "additional_metrics", "additional"):
            if hasattr(model_cls, candidate):
                return candidate
        return None

    def _extract_metadata_text_fallback(self, file_path: str) -> Dict[str, Any]:
        metadata: Dict[str, Any] = {}
        with open(file_path, "r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("##SAMPLE"):
                    start = stripped.find("<")
                    end = stripped.rfind(">")
                    if start == -1 or end == -1 or end <= start:
                        continue
                    content = stripped[start + 1 : end]
                    for item in content.split(","):
                        if "=" not in item:
                            continue
                        key, value = item.split("=", 1)
                        metadata[key.lower()] = value
                    metadata.setdefault("name", metadata.get("id"))
                if stripped.startswith("#CHROM"):
                    break
        return metadata

    def _create_allele_frequency(
        self,
        sample_group: SampleGroup,
        *,
        chrom: str,
        pos: int,
        variant_id: Optional[str],
        ref: str,
        alt: str,
        qual: Optional[float],
        filter_value: Optional[str],
        info: Optional[Info],
        format_instance: Optional[Format],
        format_sample: Optional[str],
    ) -> AlleleFrequency:
        allele = AlleleFrequency.objects.create(
            sample_group=sample_group,
            chrom=chrom,
            pos=pos,
            variant_id=variant_id,
            ref=ref,
            alt=alt,
            qual=qual,
            filter=filter_value,
            info=info,
            format=format_instance,
        )

        if format_instance and format_sample:
            payload: Dict[str, Any] = dict(format_instance.payload or {})
            additional = dict(payload.get("additional") or {})
            if additional.get("sample_id") != format_sample:
                additional["sample_id"] = format_sample
                if additional:
                    payload["additional"] = additional
                elif "additional" in payload:
                    payload.pop("additional")
                format_instance.payload = payload or None
                format_instance.save(update_fields=["payload"])

        return allele

    def _parse_vcf_text_fallback(
        self, file_path: str, sample_group: SampleGroup
    ) -> None:
        header_sample: Optional[str] = None

        with open(file_path, "r", encoding="utf-8") as handle:
            for line in handle:
                stripped = line.strip()
                if not stripped:
                    continue
                if stripped.startswith("#CHROM"):
                    columns = stripped.lstrip("#").split("\t")
                    if len(columns) > 9:
                        header_sample = columns[9]
                    continue
                if stripped.startswith("#"):
                    continue

                fields = stripped.split("\t")
                if len(fields) < 8:
                    continue

                chrom, pos, variant_id, ref, alt, qual, filter_value, info_field = fields[:8]
                info_mapping: Dict[str, Any] = {}
                for entry in info_field.split(";"):
                    if not entry:
                        continue
                    if "=" in entry:
                        key, value = entry.split("=", 1)
                        info_mapping[key] = value
                    else:
                        info_mapping[entry] = True

                info_instance = self._create_info_instance(info_mapping)

                format_instance: Optional[Format] = None
                sample_identifier: Optional[str] = None
                if len(fields) > 9:
                    format_keys = fields[8].split(":") if len(fields) > 8 else []
                    format_values = fields[9].split(":")
                    sample_identifier = header_sample or "Sample"
                    samples = (
                        {sample_identifier: dict(zip(format_keys, format_values))}
                        if format_keys
                        else {}
                    )
                    if samples:
                        format_instance, sample_identifier = self._create_format_instance(samples)

                self._create_allele_frequency(
                    sample_group,
                    chrom=chrom,
                    pos=int(pos),
                    variant_id=None if variant_id == "." else variant_id,
                    ref=ref,
                    alt=self._serialize_alt((alt,)),
                    qual=None if qual in {".", ""} else float(qual),
                    filter_value=self._serialize_filter(filter_value),
                    info=info_instance,
                    format_instance=format_instance,
                    format_sample=sample_identifier,
                )

        return None

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
            mapped_field = self.INFO_FIELD_MAP.get(normalized)
            if mapped_field:
                structured[mapped_field] = self._stringify(value)
            else:
                additional[normalized] = self._stringify(value)
            target = (
                structured if normalized in self.INFO_FIELD_MAP else additional
            )
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

            mapped_field = self.FORMAT_FIELD_MAP.get(normalized)
            if mapped_field:
                structured[mapped_field] = serialized
            else:
                additional[normalized] = serialized

        payload: Dict[str, Any] = {}
        if structured:
            payload["fields"] = structured
        if additional:
            payload["additional"] = additional

        if not payload:
            return None, sample_name

        format_instance = Format.objects.create(
            genotype=structured.get("gt"),
            payload=payload,
        )
        return format_instance, sample_name

    @staticmethod
    def _serialize_alt(alts: Optional[Iterable[str]]) -> str:
        return ",".join(alts or [])

    @staticmethod
    def _serialize_filter(filter_field: Any) -> Optional[str]:
        if not filter_field:
            return None
        if isinstance(filter_field, str):
            return filter_field
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
