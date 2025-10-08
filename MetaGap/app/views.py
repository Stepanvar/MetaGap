"""Application views."""

import os
from typing import Any, Dict, Iterable, Optional, Tuple
from pathlib import Path

import pysam
from django.conf import settings
from django.contrib import messages
from django.contrib.auth import logout
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.files.storage import default_storage
from django.db import models as django_models
from django.db import transaction
from django.urls import reverse_lazy
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

    def get_queryset(self):
        # Optimize query by selecting related metadata and prefetching variants
        return (
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
        ).prefetch_related("allele_frequencies")

        # Instantiate the filter set manually so that we can reuse it in the template context.
        self.filterset = self.filterset_class(
            self.request.GET or None, queryset=base_queryset
        )
        return self.filterset.qs.distinct()

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
            "tissue": ["tissue_type"],
            "collection_method": ["collection", "method"],
            "storage_conditions": ["storage", "storage_conditions"],
            "comments": ["description", "notes"],
        },
        "input_quality": {
            "a260_a280": ["a260_280", "ratio_a260_a280"],
            "a260_a230": ["a260_230", "ratio_a260_a230"],
            "dna_concentration": ["dna_conc", "dna_concentration_ng_ul", "concentration"],
            "rna_concentration": ["rna_conc", "rna_concentration_ng_ul"],
            "notes": ["note", "comment", "comments"],
        },
        "reference_genome_build": {
            "build_name": ["name", "reference", "build"],
            "build_version": ["version", "build_version"],
        },
        "genome_complexity": {
            "size": ["genome_size", "size_bp"],
            "ploidy": ["ploidy_level"],
            "gc_content": ["gc", "gc_percent"],
        },
        "sample_origin": {
            "tissue": ["tissue_type"],
            "collection_method": ["collection", "method"],
            "storage_conditions": ["storage", "storage_conditions"],
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
            "software": ["aligner", "software"],
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
        try:
            with pysam.VariantFile(file_path) as vcf_in:
                metadata = self.extract_sample_group_metadata(vcf_in)
        except (OSError, ValueError) as exc:  # pragma: no cover - defensive fallback
            messages.warning(
                self.request,
                f"Could not parse VCF metadata with pysam: {exc}. Falling back to a text parser.",
            )
            metadata = self._extract_metadata_text_fallback(file_path)

            created_alleles = []
            for record in vcf_in.fetch():
                info_instance = self._create_info_instance(record.info)
                format_instance, format_sample = self._create_format_instance(record.samples)

                allele = AlleleFrequency.objects.create(
                    sample_group=sample_group,
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
        sample_group = SampleGroup.objects.create(
            name=metadata.get("name")
            or os.path.splitext(os.path.basename(file_path))[0],
            created_by=organization_profile,
            comments=metadata.get("description"),
        )

        self._parse_vcf_text_fallback(file_path, sample_group)

        return sample_group

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
            format_instance.additional = {
                **(format_instance.additional or {}),
                "sample_id": format_sample,
            }
            format_instance.save(update_fields=["additional"])

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
            if normalized in self.FORMAT_FIELD_MAP:
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
