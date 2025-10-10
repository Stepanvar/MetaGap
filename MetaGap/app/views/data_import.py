"""Views that encapsulate VCF data ingestion."""

from __future__ import annotations

import json
import logging
import os
import re
from decimal import Decimal, InvalidOperation
from typing import Any, Dict, Iterable, Optional, Tuple

import pysam
from django.conf import settings
from django.contrib import messages
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.files.storage import default_storage
from django.db import models as django_models
from django.urls import reverse_lazy
from django.views.generic import FormView

from ..forms import ImportDataForm
from ..mixins import OrganizationSampleGroupMixin
from ..models import (
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
from ..services.vcf_importer import VCFImporter

logger = logging.getLogger(__name__)

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

    METADATA_SECTION_MAP = {
        "SAMPLE_GROUP": "sample_group",
        "SAMPLEGROUP": "sample_group",
        "GROUP": "sample_group",
        "SAMPLE": "sample_group",
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
        "BIOINFO_ALIGNMENT": "bioinfo_alignment",
        "ALIGNMENT": "bioinfo_alignment",
        "BIOINFO_VARIANT_CALLING": "bioinfo_variant_calling",
        "VARIANT_CALLING": "bioinfo_variant_calling",
        "BIOINFO_POSTPROC": "bioinfo_post_proc",
        "BIOINFO_POST_PROC": "bioinfo_post_proc",
        "BIOINFO_POST_PROCESSING": "bioinfo_post_proc",
        "INPUT_QUALITY": "input_quality",
        "PLATFORM_INDEPENDENT": "platform_independent",
        "PLATFORM-INDEPENDENT": "platform_independent",
        "PLATFORMINDEPENDENT": "platform_independent",
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
        "bioinfo_alignment": BioinfoAlignment,
        "bioinfo_variant_calling": BioinfoVariantCalling,
        "bioinfo_post_proc": BioinfoPostProc,
        "input_quality": InputQuality,
    }

    METADATA_FIELD_ALIASES = {
        "sample_group": {
            "name": ["group", "group_name", "dataset", "id"],
            "doi": ["dataset_doi", "group_doi"],
            "source_lab": ["lab", "lab_name", "source", "center"],
            "contact_email": ["email", "lab_email", "contact"],
            "contact_phone": ["phone", "lab_phone"],
            "total_samples": [
                "samples",
                "sample_count",
                "n_samples",
                "n",
                "num_samples",
                "number_of_samples",
            ],
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
        "bioinfo_alignment": "tool",
        "bioinfo_variant_calling": "tool",
        "bioinfo_post_proc": "normalization",
    }

    INFO_FIELD_STRING = "string"
    INFO_FIELD_INT = "int"
    INFO_FIELD_FLOAT = "float"
    INFO_PLACEHOLDER_VALUES = {".", ""}

    INFO_FIELD_MAP = {
        "aa": ("aa", INFO_FIELD_STRING),
        "ac": ("ac", INFO_FIELD_STRING),
        "af": ("af", INFO_FIELD_STRING),
        "an": ("an", INFO_FIELD_STRING),
        "bq": ("bq", INFO_FIELD_STRING),
        "cigar": ("cigar", INFO_FIELD_STRING),
        "db": ("db", INFO_FIELD_STRING),
        "dp": ("dp", INFO_FIELD_STRING),
        "end": ("end", INFO_FIELD_STRING),
        "h2": ("h2", INFO_FIELD_STRING),
        "h3": ("h3", INFO_FIELD_STRING),
        "mq": ("mq", INFO_FIELD_STRING),
        "mq0": ("mq0", INFO_FIELD_STRING),
        "ns": ("ns", INFO_FIELD_STRING),
        "qd": ("qd", INFO_FIELD_STRING),
        "fs": ("fs", INFO_FIELD_STRING),
        "sor": ("sor", INFO_FIELD_STRING),
        "sb": ("sb", INFO_FIELD_STRING),
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
        group_data, group_consumed, group_additional = self._extract_section_data(
            metadata, "sample_group", SampleGroup
        )

        consumed_keys = set(group_consumed)
        additional_metadata: Dict[str, Any] = {}
        if group_additional:
            additional_metadata.update(group_additional)

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

            consumed_keys.update(section_consumed)
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
            if section not in update_fields:
                update_fields.append(section)

        for key, value in metadata.items():
            if key in consumed_keys:
                continue
            coerced = self._coerce_additional_value(value)
            if key not in additional_metadata:
                additional_metadata[key] = coerced

        additional_payload = additional_metadata or None
        if additional_payload != getattr(sample_group, "additional_metadata", None):
            sample_group.additional_metadata = additional_payload
            if "additional_metadata" not in update_fields:
                update_fields.append("additional_metadata")

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

        additional, additional_consumed = self._build_additional_payload(
            metadata, section, consumed
        )
        consumed.update(additional_consumed)
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
            for suffix in ("_", ".", "-"):
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
    ) -> Tuple[Optional[Dict[str, Any]], set[str]]:
        consumed_set = set(consumed)
        prefixes = (f"{section}_", f"{section}.", f"{section}-")
        additional: Dict[str, Any] = {}
        additional_consumed: set[str] = set()

        for key, value in metadata.items():
            if key in consumed_set:
                continue
            for prefix in prefixes:
                if key.startswith(prefix):
                    trimmed = key[len(prefix) :]
                    additional[trimmed] = self._coerce_additional_value(value)
                    additional_consumed.add(key)
                    break

        payload = additional or None
        return payload, additional_consumed

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

    @staticmethod
    def _split_sample_attributes(content: str) -> Iterable[str]:
        items: list[str] = []
        current: list[str] = []
        quote_char: Optional[str] = None
        escape = False
        bracket_stack: list[str] = []
        opening = {"{": "}", "[": "]", "(": ")"}
        closing = {value: key for key, value in opening.items()}

        for char in content:
            if quote_char:
                current.append(char)
                if escape:
                    escape = False
                    continue
                if char == "\\":
                    escape = True
                    continue
                if char == quote_char:
                    quote_char = None
                continue

            if char in {'"', "'"}:
                quote_char = char
                current.append(char)
                continue

            if char in opening:
                bracket_stack.append(char)
                current.append(char)
                continue

            if char in closing:
                if bracket_stack and bracket_stack[-1] == closing[char]:
                    bracket_stack.pop()
                current.append(char)
                continue

            if char == "," and not bracket_stack:
                item = "".join(current).strip()
                if item:
                    items.append(item)
                current = []
                continue

            current.append(char)

        tail = "".join(current).strip()
        if tail:
            items.append(tail)

        return items

    @staticmethod
    def _normalize_metadata_value(value: str) -> str:
        stripped = value.strip()
        if len(stripped) >= 2 and stripped[0] == stripped[-1] and stripped[0] in {'"', "'"}:
            stripped = stripped[1:-1]
        if "\\" in stripped:
            try:
                stripped = bytes(stripped, "utf-8").decode("unicode_escape")
            except UnicodeDecodeError:  # pragma: no cover - defensive decoding
                pass
        return stripped

    @staticmethod
    def _normalize_metadata_key(key: Any) -> str:
        normalized = re.sub(r"[^0-9a-z]+", "", str(key).lower())
        return normalized

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
                    for item in self._split_sample_attributes(content):
                        if "=" not in item:
                            continue
                        key, value = item.split("=", 1)
                        metadata[key.lower()] = self._normalize_metadata_value(value)
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
            key = (record.key or "").upper()
            if key == "SEQUENCING_PLATFORM":
                self._ingest_sequencing_platform_record(metadata, record)
                continue

            section = self.METADATA_SECTION_MAP.get(key)
            if not section:
                continue

            items = self._collect_record_items(record)
            self._process_metadata_section(metadata, section, items)

        if "name" not in metadata and "sample_group_name" in metadata:
            metadata["name"] = metadata["sample_group_name"]
        return metadata

    def _ingest_sequencing_platform_record(
        self, metadata: Dict[str, Any], record: pysam.libcbcf.VariantHeaderRecord
    ) -> None:
        items = self._collect_record_items(record)
        platform_value = None
        for candidate in ("platform", "Platform"):
            if candidate in items:
                platform_value = items[candidate]
                break

        section = self._determine_platform_section(platform_value)
        if not section:
            logger.warning(
                "Unsupported sequencing platform %r; storing raw metadata only.",
                platform_value,
            )
            section = "platform_independent"

        self._process_metadata_section(metadata, section, items)

    def _determine_platform_section(self, platform_value: Any) -> Optional[str]:
        if not platform_value:
            return None

        normalized = str(platform_value).strip().lower()
        if "illumina" in normalized:
            return "illumina_seq"
        if "nanopore" in normalized or "ont" in normalized or "oxford" in normalized:
            return "ont_seq"
        if "pacbio" in normalized or "sequel" in normalized or "revio" in normalized:
            return "pacbio_seq"
        if "ion" in normalized:
            return "iontorrent_seq"
        return None

    def _collect_record_items(
        self, record: pysam.libcbcf.VariantHeaderRecord
    ) -> Dict[str, Any]:
        collected: Dict[str, Any] = {}
        for key, value in record.items():
            normalized_key = str(key)
            if isinstance(value, str):
                normalized_value = self._normalize_metadata_value(value)
            else:
                normalized_value = value
            collected[normalized_key] = normalized_value
        return collected

    def _process_metadata_section(
        self, metadata: Dict[str, Any], section: str, items: Dict[str, Any]) -> None:
        def normalize_alias_key(candidate: str) -> str:
            stripped = candidate.rstrip("?!.,;:")
            return stripped or candidate

        # alias lookup including punctuation-stripped variants
        alias_lookup: Dict[str, str] = {}
        for field_name, aliases in self.METADATA_FIELD_ALIASES.get(section, {}).items():
            canonical = field_name.lower()
            for candidate in [canonical, *aliases]:
                norm = str(candidate).lower()
                alias_lookup[norm] = field_name
                stripped = normalize_alias_key(norm)
                if stripped != norm:
                    alias_lookup[stripped] = field_name

        recognized: Dict[str, Any] = {}
        leftovers: Dict[str, Any] = {}

        for raw_key, value in items.items():
            nk = str(raw_key).lower()
            sk = normalize_alias_key(nk)
            field = alias_lookup.get(nk) or alias_lookup.get(sk)
            if field:
                recognized[field] = value
            else:
                leftovers[sk] = value

        for field_name, value in recognized.items():
            if section == "sample_group":
                metadata[field_name] = value
                if field_name == "name" and value is not None:
                    metadata.setdefault("name", value)
                elif field_name == "comments" and value is not None:
                    metadata.setdefault("comments", value)
            else:
                metadata_key = f"{section}_{field_name}"
                metadata[metadata_key] = value

        custom_keys: list[str] = []
        for raw_key, value in leftovers.items():
            if not raw_key:
                continue

            normalized_custom_key = self._normalize_metadata_key(raw_key)
            if not normalized_custom_key:
                continue

            metadata_key = f"{section}_{normalized_custom_key}"
            suffix = 2
            while metadata_key in metadata:
                metadata_key = f"{section}_{normalized_custom_key}_{suffix}"
                suffix += 1

            metadata[metadata_key] = value
            custom_keys.append(str(raw_key))

        if custom_keys:
            joined_keys = ", ".join(sorted(custom_keys))
            logger.info(
                "Preserved custom metadata keys (%s) for section '%s'.",
                joined_keys,
                section,
            )

    def _create_info_instance(self, info: Any) -> Optional[Info]:
        info_dict = dict(info)
        structured: Dict[str, Any] = {}
        additional: Dict[str, Any] = {}

        for key, value in info_dict.items():
            normalized = key.lower()
            mapped_field = self.INFO_FIELD_MAP.get(normalized)
            if mapped_field:
                field_name, field_type = mapped_field
                structured[field_name] = self._coerce_info_value(value, field_type)
            else:
                additional_value = self._normalize_additional_info_value(value)
                if additional_value is not None:
                    additional[normalized] = additional_value

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

    @classmethod
    def _coerce_info_value(cls, value: Any, field_type: str) -> Any:
        normalized = cls._normalize_info_scalar(value)
        return cls._stringify(normalized)

    @classmethod
    def _normalize_info_scalar(cls, value: Any) -> Any:
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            normalized_items = []
            for item in value:
                normalized_item = cls._normalize_info_scalar(item)
                if normalized_item is None:
                    continue
                normalized_items.append(normalized_item)
            return normalized_items or None
        if isinstance(value, str):
            stripped = value.strip()
            if not stripped or stripped in cls.INFO_PLACEHOLDER_VALUES:
                return None
            return stripped
        return value

    @staticmethod
    def _coerce_int(value: Any) -> Optional[int]:
        if value is None:
            return None
        if isinstance(value, bool):
            return int(value)
        try:
            decimal_value = Decimal(str(value))
        except (InvalidOperation, TypeError, ValueError):
            return None

        try:
            integral_value = decimal_value.to_integral_value()
        except InvalidOperation:
            return None

        if integral_value != decimal_value:
            return None
        return int(integral_value)

    @staticmethod
    def _coerce_float(value: Any) -> Optional[float]:
        if value is None:
            return None
        if isinstance(value, bool):
            return float(value)
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @classmethod
    def _normalize_additional_info_value(cls, value: Any) -> Any:
        if isinstance(value, (list, tuple)):
            normalized_list = []
            for item in value:
                normalized_item = cls._normalize_additional_info_value(item)
                if normalized_item is None:
                    continue
                normalized_list.append(normalized_item)
            return normalized_list or None
        if value in (None, "", "."):
            return None
        if isinstance(value, (int, float, bool)):
            return value
        return str(value)

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

__all__ = ["ImportDataView"]

