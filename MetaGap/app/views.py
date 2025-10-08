# views.py

import os
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
from django.views.generic import CreateView, FormView, ListView, TemplateView
from django_tables2 import RequestConfig
from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin

from .forms import (
    CustomUserCreationForm,
    DeleteAccountForm,
    EditProfileForm,
    ImportDataForm,
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
        # Get the dynamically created table class
        SampleGroupTable = create_dynamic_table(
            SampleGroup, table_name="SampleGroupTable", include_related=True
        )
        table = SampleGroupTable(self.get_queryset())
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
    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["form"] = SearchForm()
        return context
    def get(self, request, *args, **kwargs):
        form = SearchForm()
        return self.render_to_response({'form': form})

class ProfileView(LoginRequiredMixin, TemplateView):
    template_name = "profile.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        organization_profile = getattr(self.request.user, "organization_profile", None)
        if organization_profile:
            sample_groups = SampleGroup.objects.filter(
                created_by=organization_profile
            ).order_by("name")
        else:
            sample_groups = SampleGroup.objects.none()

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
        # Save the uploaded file temporarily
        file_path = default_storage.save("tmp/" + data_file.name, data_file)
        full_path = os.path.join(settings.MEDIA_ROOT, file_path)
        try:
            self.parse_vcf_file(full_path, data_file.name, self.request.user)
            messages.success(self.request, "Data imported successfully.")
        except Exception as e:
            messages.error(self.request, f"An error occurred: {e}")
        finally:
            # Delete the temporary file
            default_storage.delete(file_path)
        return super().form_valid(form)

    def parse_vcf_file(self, file_path, original_filename, user):
        organization_profile = getattr(user, "organization_profile", None)
        if organization_profile is None:
            raise ValueError("The current user does not have an organization profile.")

        vcf_in = pysam.VariantFile(file_path)  # auto-detect input format (VCF or BCF)
        try:
            metadata = self.extract_sample_group_metadata(vcf_in)
            default_name = Path(original_filename).stem

            with transaction.atomic():
                sample_group = self.create_sample_group_with_metadata(
                    metadata,
                    organization_profile,
                    default_name,
                )

                if hasattr(vcf_in, "reset"):
                    vcf_in.reset()
                allele_frequency_objects = self._create_allele_frequencies(vcf_in, sample_group)
                if allele_frequency_objects:
                    AlleleFrequency.objects.bulk_create(allele_frequency_objects)
        finally:
            vcf_in.close()

    def extract_sample_group_metadata(self, vcf_in):
        metadata = {}
        header = vcf_in.header
        for record in header.records:
            section_key, values = self._parse_header_record(record)
            if not section_key or not values:
                continue
            section_name = self.METADATA_SECTION_MAP.get(section_key)
            if not section_name:
                continue
            section = metadata.setdefault(section_name, {})
            section.update(values)
        return metadata

    def create_sample_group_with_metadata(self, metadata, organization_profile, default_name):
        sample_group_metadata = metadata.get("sample_group", {})
        sample_group_kwargs = self._prepare_sample_group_kwargs(
            sample_group_metadata,
            default_name,
            organization_profile,
        )

        input_quality_obj = None
        input_quality_metadata = metadata.get("input_quality")
        if input_quality_metadata:
            input_quality_kwargs = self._prepare_section_kwargs(
                "input_quality",
                input_quality_metadata,
                InputQuality,
            )
            if input_quality_kwargs:
                input_quality_obj = InputQuality.objects.create(**input_quality_kwargs)
        sample_group_kwargs["input_quality"] = input_quality_obj

        metadata_objects = self._build_metadata_objects(metadata)
        sample_group_kwargs.update(metadata_objects)

        return SampleGroup.objects.create(**sample_group_kwargs)

    def _prepare_sample_group_kwargs(self, metadata, default_name, organization_profile):
        normalized = {self._normalize_header_key(k): v for k, v in metadata.items()}
        aliases = self.METADATA_FIELD_ALIASES.get("sample_group", {})
        fields = [
            "name",
            "doi",
            "source_lab",
            "contact_email",
            "contact_phone",
            "total_samples",
            "inclusion_criteria",
            "exclusion_criteria",
            "tissue",
            "collection_method",
            "storage_conditions",
            "comments",
        ]
        kwargs = {}
        for field in fields:
            value = normalized.get(field)
            if value in (None, ""):
                for alias in aliases.get(field, []):
                    alias_value = normalized.get(alias)
                    if alias_value not in (None, ""):
                        value = alias_value
                        break
            if field == "total_samples":
                value = self._to_int(value)
            elif field in {"contact_phone", "contact_email", "doi", "source_lab", "tissue", "collection_method", "storage_conditions"}:
                value = self._clean_string(value)
            elif field in {"inclusion_criteria", "exclusion_criteria", "comments"} and value is not None:
                value = str(value)
            elif field == "name":
                value = self._clean_string(value)
            kwargs[field] = value

        kwargs["name"] = kwargs.get("name") or default_name or "Sample Group"
        kwargs["created_by"] = organization_profile
        return kwargs

    def _build_metadata_objects(self, metadata):
        created = {}
        for section, model in self.METADATA_MODEL_MAP.items():
            section_data = metadata.get(section)
            if not section_data:
                continue
            kwargs = self._prepare_section_kwargs(section, section_data, model)
            if not kwargs:
                continue
            created[section] = model.objects.create(**kwargs)
        return created

    def _prepare_section_kwargs(self, section_key, data, model):
        normalized = {self._normalize_header_key(k): v for k, v in data.items()}
        aliases = self.METADATA_FIELD_ALIASES.get(section_key, {})
        kwargs = {}
        consumed_keys = set()
        json_fields = {
            field.name
            for field in model._meta.get_fields()
            if isinstance(field, django_models.Field)
            and not field.auto_created
            and field.get_internal_type() == "JSONField"
        }

        for field in model._meta.get_fields():
            if not isinstance(field, django_models.Field) or field.auto_created:
                continue
            if field.name in {"id", "sample_group"}:
                continue
            if isinstance(field, (django_models.ForeignKey, django_models.OneToOneField)):
                continue

            value = None
            consumed_key = None
            raw_value = normalized.get(field.name)
            if raw_value not in (None, ""):
                value = raw_value
                consumed_key = field.name
            else:
                for alias in aliases.get(field.name, []):
                    alias_value = normalized.get(alias)
                    if alias_value not in (None, ""):
                        value = alias_value
                        consumed_key = alias
                        break

            if consumed_key is None:
                continue

            cast_value = self._cast_field_value(field, value)
            if cast_value is None and field.get_internal_type() != "JSONField":
                consumed_keys.add(consumed_key)
                continue
            kwargs[field.name] = cast_value
            consumed_keys.add(consumed_key)

        extras = {
            key: value
            for key, value in normalized.items()
            if key not in consumed_keys and value not in (None, "")
        }
        if extras and json_fields:
            json_field_name = next(iter(json_fields))
            kwargs.setdefault(json_field_name, extras)

        return kwargs

    def _cast_field_value(self, field, value):
        if value in (None, ""):
            return None

        internal_type = field.get_internal_type()
        if internal_type in {"IntegerField", "BigIntegerField", "PositiveIntegerField", "SmallIntegerField"}:
            return self._to_int(value)
        if internal_type in {"FloatField", "DecimalField"}:
            try:
                return float(value)
            except (TypeError, ValueError):
                return None
        if internal_type == "BooleanField":
            if isinstance(value, bool):
                return value
            value_str = str(value).strip().lower()
            if value_str in {"true", "1", "yes"}:
                return True
            if value_str in {"false", "0", "no"}:
                return False
            return None
        if internal_type == "JSONField":
            return value if isinstance(value, (dict, list)) else {"value": value}
        return str(value)

    def _to_int(self, value):
        if value in (None, ""):
            return None
        try:
            if isinstance(value, str):
                value = value.strip()
            return int(float(value))
        except (TypeError, ValueError):
            return None

    def _clean_string(self, value):
        if value in (None, ""):
            return None
        if isinstance(value, bool):
            return "true" if value else "false"
        return str(value).strip() or None

    def _parse_header_record(self, record):
        record_text = str(record).strip()
        if not record_text.startswith("##"):
            return None, {}
        record_text = record_text[2:]
        if "=" not in record_text:
            return None, {}
        key, remainder = record_text.split("=", 1)
        key = key.strip().upper()
        remainder = remainder.strip()
        if remainder.startswith("<") and remainder.endswith(">"):
            remainder = remainder[1:-1]
            data = {}
            for part in self._split_metadata_values(remainder):
                if "=" not in part:
                    continue
                part_key, part_value = part.split("=", 1)
                data[self._normalize_header_key(part_key)] = self._clean_header_value(part_value)
            return key, data
        cleaned_value = self._clean_header_value(remainder)
        return key, {self._normalize_header_key(key): cleaned_value}

    def _split_metadata_values(self, value):
        parts = []
        current = []
        in_quotes = False
        for char in value:
            if char == '"':
                in_quotes = not in_quotes
            if char == "," and not in_quotes:
                part = "".join(current).strip()
                if part:
                    parts.append(part)
                current = []
            else:
                current.append(char)
        if current:
            part = "".join(current).strip()
            if part:
                parts.append(part)
        return parts

    def _clean_header_value(self, value):
        value = value.strip()
        if value.startswith("<") and value.endswith(">"):
            value = value[1:-1]
        if value.startswith('"') and value.endswith('"'):
            value = value[1:-1]
        return value.strip()

    def _normalize_header_key(self, key):
        key = key.replace("-", "_").replace("/", "_").replace(" ", "_")
        normalized = []
        for char in key:
            if char.isalnum() or char == "_":
                normalized.append(char.lower())
            else:
                normalized.append("_")
        normalized_key = "".join(normalized).strip("_")
        while "__" in normalized_key:
            normalized_key = normalized_key.replace("__", "_")
        return normalized_key

    def _create_allele_frequencies(self, vcf_in, sample_group):
        allele_frequency_objects = []
        for record in vcf_in:
            info_obj = self._create_info_object(record)
            format_obj = self._create_format_object(record)

            alt_value = ",".join(record.alts) if record.alts else ""
            allele_frequency = AlleleFrequency(
                sample_group=sample_group,
                chrom=record.chrom,
                pos=record.pos,
                variant_id=record.id or None,
                ref=record.ref,
                alt=alt_value,
                qual=record.qual if record.qual is not None else None,
                filter=self._stringify_filter(record),
                info=info_obj,
                format=format_obj,
            )
            allele_frequency_objects.append(allele_frequency)
        return allele_frequency_objects

    def _create_info_object(self, record):
        info_data = {}
        additional = {}
        for key, value in record.info.items():
            normalized_key = self._normalize_header_key(key)
            string_value = self._stringify_variant_value(value)
            if normalized_key in self.INFO_FIELD_MAP:
                info_data[self.INFO_FIELD_MAP[normalized_key]] = string_value
            else:
                additional[normalized_key] = self._normalize_variant_value(value)

        if additional:
            info_data["additional"] = additional

        if not info_data:
            return None
        return Info.objects.create(**info_data)

    def _create_format_object(self, record):
        if not record.samples:
            return None
        sample = next(iter(record.samples.values()))
        format_data = {}
        additional = {}
        for key, value in sample.items():
            normalized_key = self._normalize_header_key(key)
            string_value = self._stringify_variant_value(value)
            if normalized_key in self.FORMAT_FIELD_MAP:
                format_data[self.FORMAT_FIELD_MAP[normalized_key]] = string_value
            else:
                additional[normalized_key] = self._normalize_variant_value(value)

        if additional:
            format_data["additional"] = additional

        if not format_data:
            return None
        return Format.objects.create(**format_data)

    def _normalize_variant_value(self, value):
        if isinstance(value, bytes):
            return value.decode(errors="ignore")
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        if isinstance(value, (list, tuple)):
            return [self._normalize_variant_value(v) for v in value]
        return str(value)

    def _stringify_variant_value(self, value):
        normalized = self._normalize_variant_value(value)
        if normalized is None:
            return None
        if isinstance(normalized, list):
            return ",".join(str(v) for v in normalized)
        if isinstance(normalized, bool):
            return "true" if normalized else "false"
        return str(normalized)

    def _stringify_filter(self, record):
        keys_method = getattr(record.filter, "keys", None)
        if callable(keys_method):
            keys = list(keys_method())
        else:
            keys = list(record.filter)
        if not keys:
            return None
        return ",".join(str(key) for key in keys)
