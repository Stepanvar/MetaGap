"""Metadata parsing utilities for VCF imports."""

from __future__ import annotations

import logging
import re
from typing import Any, Dict, Iterable, Optional

import pysam

from ..models import (
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    GenomeComplexity,
    IlluminaSeq,
    InputQuality,
    IonTorrentSeq,
    LibraryConstruction,
    MaterialType,
    OntSeq,
    PacBioSeq,
    ReferenceGenomeBuild,
    SampleOrigin,
)

logger = logging.getLogger(__name__)


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


class VCFMetadataParser:
    """Parse metadata embedded in a VCF header."""

    def __init__(self, warnings: Optional[list[str]] = None) -> None:
        self.warnings = warnings if warnings is not None else []

    def extract_sample_group_metadata(self, vcf_in: pysam.VariantFile) -> Dict[str, Any]:
        metadata: Dict[str, Any] = {}
        for record in vcf_in.header.records:
            items = self._collect_record_items(record)
            self.ingest_metadata_items(metadata, record.key, items)

        if "name" not in metadata and "sample_group_name" in metadata:
            metadata["name"] = metadata["sample_group_name"]
        return metadata

    def ingest_metadata_items(
        self,
        metadata: Dict[str, Any],
        key: Optional[str],
        items: Dict[str, Any],
    ) -> None:
        normalized_key = (key or "").upper()
        if normalized_key == "SEQUENCING_PLATFORM":
            self._ingest_sequencing_platform_items(metadata, items)
            return

        section = METADATA_SECTION_MAP.get(normalized_key)
        if not section:
            if self._is_metadata_section_candidate(normalized_key):
                log_message = (
                    f"Unsupported metadata section '{key}' encountered in the VCF header."
                )
                logger.warning("%s", log_message)
                warning = f"Unsupported metadata section '{key}'"
                self.warnings.append(warning)
                additional = metadata.setdefault("additional_metadata", {})
                additional_key = normalize_metadata_key(key)
                additional[additional_key] = items
            return

        self._process_metadata_section(metadata, section, items)

    def _ingest_sequencing_platform_items(
        self, metadata: Dict[str, Any], items: Dict[str, Any]
    ) -> None:
        platform_value = None
        for candidate in ("platform", "Platform"):
            if candidate in items:
                platform_value = items[candidate]
                break

        section = self._determine_platform_section(platform_value)
        if not section:
            warning = (
                "Unsupported sequencing platform "
                f"{platform_value!r}; storing raw metadata only."
            )
            logger.warning("%s", warning)
            self.warnings.append(warning)
            section = "platform_independent"

        self._process_metadata_section(metadata, section, items)

    @staticmethod
    def _is_metadata_section_candidate(key: str) -> bool:
        normalized = (key or "").upper()
        metadata_prefixes = (
            "SAMPLE_",
            "SAMPLEGROUP",
            "SAMPLE-GROUP",
            "GROUP_",
        )
        return any(normalized.startswith(prefix) for prefix in metadata_prefixes)

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
        if "iontorrent" in normalized or "ion" in normalized:
            return "iontorrent_seq"
        return None

    def _collect_record_items(
        self, record: pysam.libcbcf.VariantHeaderRecord
    ) -> Dict[str, Any]:
        collected: Dict[str, Any] = {}
        for key, value in record.items():
            normalized_key = str(key)
            if isinstance(value, str):
                normalized_value = normalize_metadata_value(value)
            else:
                normalized_value = value
            collected[normalized_key] = normalized_value
        return collected

    def _process_metadata_section(
        self, metadata: Dict[str, Any], section: str, items: Dict[str, Any]
    ) -> None:
        def normalize_alias_key(candidate: str) -> str:
            stripped = candidate.rstrip("?!.,;:")
            return normalize_metadata_key(stripped)

        alias_map = METADATA_FIELD_ALIASES.get(section, {})
        normalized_section = normalize_metadata_key(section)

        for key, value in items.items():
            normalized_key = normalize_metadata_key(key)
            value_key = normalized_key

            for field_name, aliases in alias_map.items():
                alias_candidates = [normalize_alias_key(alias) for alias in aliases]
                if normalized_key in alias_candidates:
                    normalized_key = normalize_metadata_key(field_name)
                    break

            qualified_key = f"{normalized_section}_{value_key}"
            metadata[normalized_key] = value
            metadata[qualified_key] = value

        if section not in metadata:
            primary_aliases: Iterable[str] = alias_map.get(section, [])
            for alias in primary_aliases:
                alias_key = normalize_metadata_key(alias)
                if alias_key in metadata:
                    metadata[section] = metadata[alias_key]
                    break


def normalize_metadata_value(value: str) -> str:
    stripped = value.strip()
    if len(stripped) >= 2 and stripped[0] == stripped[-1] and stripped[0] in {'"', "'"}:
        stripped = stripped[1:-1]
    if "\\" in stripped:
        try:
            stripped = bytes(stripped, "utf-8").decode("unicode_escape")
        except UnicodeDecodeError:
            pass
    return stripped


def normalize_metadata_key(key: Any) -> str:
    """Normalize metadata keys to snake_case-friendly strings.

    The importer expects metadata keys to be lower case and delimited by
    underscores.  Previously, normalization removed non-alphanumeric
    characters entirely, which collapsed delimiters and caused section
    qualified keys (e.g., ``sample_group_contact_email``) to lose their
    underscores.  This new implementation replaces runs of
    non-alphanumeric characters with single underscores and trims leading or
    trailing underscores to maintain consistent snake_case keys.
    """

    normalized = re.sub(r"[^0-9a-z]+", "_", str(key).lower())
    return normalized.strip("_")


__all__ = [
    "METADATA_FIELD_ALIASES",
    "METADATA_MODEL_MAP",
    "METADATA_SECTION_MAP",
    "SECTION_PRIMARY_FIELD",
    "VCFMetadataParser",
    "normalize_metadata_key",
    "normalize_metadata_value",
]
