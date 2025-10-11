"""Metadata parsing utilities for VCF imports."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from importlib import import_module
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

import pysam
import yaml

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
    SampleGroup,
    SampleOrigin,
)

logger = logging.getLogger(__name__)


class MetadataConfigurationError(RuntimeError):
    """Raised when the metadata configuration file cannot be loaded."""


@dataclass(frozen=True)
class MetadataConfiguration:
    """Container for metadata configuration entries."""

    section_map: Mapping[str, str]
    models: Mapping[str, Any]
    field_aliases: Mapping[str, Mapping[str, Iterable[str]]]
    section_primary_field: Mapping[str, str]


_CONFIG_CACHE: dict[Path, MetadataConfiguration] = {}


def _default_config_path() -> Path:
    return Path(__file__).resolve().parent.parent / "config" / "metadata_fields.yaml"


def load_metadata_configuration(config_path: Optional[Path | str] = None) -> MetadataConfiguration:
    """Load the metadata configuration from disk.

    The configuration is cached in memory so repeated lookups avoid touching the
    filesystem.  Callers may supply an alternate ``config_path`` for testing or
    advanced scenarios.
    """

    path = Path(config_path) if config_path is not None else _default_config_path()
    path = path.resolve()

    cached = _CONFIG_CACHE.get(path)
    if cached is not None:
        return cached

    if not path.exists():
        raise MetadataConfigurationError(
            f"Metadata configuration file not found: {path}"
        )

    try:
        with path.open("r", encoding="utf-8") as handle:
            raw_config = yaml.safe_load(handle) or {}
    except yaml.YAMLError as exc:  # pragma: no cover - rare parse errors
        raise MetadataConfigurationError(
            f"Failed to parse metadata configuration: {exc}"
        ) from exc

    if not isinstance(raw_config, dict):
        raise MetadataConfigurationError(
            "Metadata configuration must be a mapping at the top level."
        )

    section_map = _coerce_str_dict(
        raw_config.get("metadata_section_map", {}),
        "metadata_section_map",
    )
    models = _load_model_map(raw_config.get("metadata_models", {}))
    field_aliases = _load_field_aliases(raw_config.get("metadata_field_aliases", {}))
    section_primary_field = _coerce_str_dict(
        raw_config.get("section_primary_field", {}),
        "section_primary_field",
    )

    configuration = MetadataConfiguration(
        section_map=section_map,
        models=models,
        field_aliases=field_aliases,
        section_primary_field=section_primary_field,
    )

    _CONFIG_CACHE[path] = configuration
    return configuration


def _coerce_str_dict(raw: Any, field_name: str) -> Dict[str, str]:
    if not isinstance(raw, dict):
        raise MetadataConfigurationError(
            f"'{field_name}' must be defined as a mapping in the metadata configuration."
        )

    coerced: Dict[str, str] = {}
    for key, value in raw.items():
        if not isinstance(key, str) or not isinstance(value, str):
            raise MetadataConfigurationError(
                f"Entries in '{field_name}' must map strings to strings."
            )
        coerced[key] = value
    return coerced


def _load_model_map(raw: Any) -> Dict[str, Any]:
    mapping = _coerce_str_dict(raw, "metadata_models")
    resolved: Dict[str, Any] = {}
    for section, import_path in mapping.items():
        module_path, _, attribute = import_path.rpartition(".")
        if not module_path or not attribute:
            raise MetadataConfigurationError(
                f"Invalid model import path '{import_path}' for section '{section}'."
            )
        try:
            module = import_module(module_path)
        except ModuleNotFoundError as exc:
            raise MetadataConfigurationError(
                f"Could not import module '{module_path}' for metadata model '{section}'."
            ) from exc
        try:
            resolved_model = getattr(module, attribute)
        except AttributeError as exc:
            raise MetadataConfigurationError(
                f"Module '{module_path}' does not define '{attribute}' for metadata section '{section}'."
            ) from exc
        resolved[section] = resolved_model
    return resolved


def _load_field_aliases(raw: Any) -> Dict[str, Dict[str, list[str]]]:
    if not isinstance(raw, dict):
        raise MetadataConfigurationError(
            "'metadata_field_aliases' must be defined as a mapping of sections."
        )

    field_aliases: Dict[str, Dict[str, list[str]]] = {}
    for section, section_aliases in raw.items():
        if not isinstance(section, str) or not isinstance(section_aliases, dict):
            raise MetadataConfigurationError(
                "Each entry in 'metadata_field_aliases' must map a section name to a mapping of fields."
            )
        aliases_for_section: Dict[str, list[str]] = {}
        for field_name, aliases in section_aliases.items():
            if not isinstance(field_name, str):
                raise MetadataConfigurationError(
                    "Field names in 'metadata_field_aliases' must be strings."
                )
            if aliases is None:
                aliases_list: list[str] = []
            elif isinstance(aliases, (list, tuple)):
                aliases_list = [str(alias) for alias in aliases]
            else:
                raise MetadataConfigurationError(
                    f"Aliases for field '{field_name}' in section '{section}' must be a list."
                )
            aliases_for_section[field_name] = aliases_list
        field_aliases[section] = aliases_for_section
    return field_aliases
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

PLATFORM_SECTION_MAP = SampleGroup.PLATFORM_FIELD_MAP
PLATFORM_SECTION_FIELDS = set(PLATFORM_SECTION_MAP.values())
SECTION_PLATFORM_MAP = {
    field_name: platform
    for platform, field_name in PLATFORM_SECTION_MAP.items()
}

METADATA_FIELD_ALIASES = {
    "sample_group": {
        "name": ["group", "group_name", "dataset", "id"],
        "doi": ["dataset_doi", "group_doi"],
        "source_lab": ["lab", "lab_name", "source", "center"],
        "contact_email": ["email", "lab_email", "contact"],
        "contact_phone": ["phone", "lab_phone"],
        "sequencing_platform": ["platform", "sequencing_platform"],
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
        "build_name": ["name", "reference", "build", "buildname"],
        "build_version": ["version", "build_version", "buildversion"],
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
        "pcr_cycles": ["pcr", "pcr_cycles", "pcrcycle", "pcrcycles", "pcrcyles"],
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

    def __init__(
        self,
        warnings: Optional[list[str]] = None,
        *,
        configuration: Optional[MetadataConfiguration] = None,
    ) -> None:
        self.warnings = warnings if warnings is not None else []
        self.configuration = configuration or load_metadata_configuration()

    def extract_sample_group_metadata(self, vcf_in: pysam.VariantFile) -> Dict[str, Any]:
        metadata: Dict[str, Any] = {}
        self._consumed_keys.clear()
        for record in vcf_in.header.records:
            items = self._collect_record_items(record)
            self.ingest_metadata_items(metadata, record.key, items)

        if "name" not in metadata and "sample_group_name" in metadata:
            metadata["name"] = metadata["sample_group_name"]

        if self._consumed_keys:
            metadata["_consumed_keys"] = sorted(self._consumed_keys)
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

        section = self.configuration.section_map.get(normalized_key)
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
        platform_choice = None
        if section:
            platform_choice = SECTION_PLATFORM_MAP.get(section)

        if not section:
            warning = (
                "Unsupported sequencing platform "
                f"{platform_value!r}; storing raw metadata only."
            )
            logger.warning("%s", warning)
            self.warnings.append(warning)
            section = "platform_independent"

        if platform_choice is not None:
            metadata["sequencing_platform"] = platform_choice.value
            self._record_section_keys("sample_group", {"sequencing_platform"})
            self._active_platform_section = section
            self._prune_inactive_platform_sections(metadata)

        restrict = section in PLATFORM_SECTION_FIELDS
        self._process_metadata_section(
            metadata, section, items, restrict_to_section=restrict
        )

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
        self,
        metadata: Dict[str, Any],
        section: str,
        items: Dict[str, Any],
        *,
        restrict_to_section: bool = False,
    ) -> None:
        if (
            section in PLATFORM_SECTION_FIELDS
            and self._active_platform_section is not None
            and section != self._active_platform_section
        ):
            return

        def normalize_alias_key(candidate: str) -> str:
            stripped = candidate.rstrip("?!.,;:")
            return normalize_metadata_key(stripped)

        alias_map = self.configuration.field_aliases.get(section, {})
        normalized_section = normalize_metadata_key(section)
        section_values: Dict[str, Any] = {}
        section_keys: set[str] = set()

        if restrict_to_section and alias_map:
            alias_keys_to_remove = {
                normalize_metadata_key(field_name)
                for field_name in alias_map.keys()
            }
            for aliases in alias_map.values():
                for alias in aliases:
                    alias_keys_to_remove.add(normalize_alias_key(alias))

            for alias_key in alias_keys_to_remove:
                metadata.pop(alias_key, None)

        for key, value in items.items():
            normalized_key = normalize_metadata_key(key)

            for field_name, aliases in alias_map.items():
                canonical_field = normalize_metadata_key(field_name)
                alias_candidates = {canonical_field, canonical_field.replace("_", "")}
                for alias in aliases:
                    normalized_alias = normalize_alias_key(alias)
                    if normalized_alias:
                        alias_candidates.add(normalized_alias)
                        alias_candidates.add(normalized_alias.replace("_", ""))
                if normalized_key in alias_candidates:
                    normalized_key = canonical_field
                    break

            section_values[normalized_key] = value
            if normalized_key.startswith(f"{normalized_section}_"):
                qualified_key = normalized_key
            else:
                qualified_key = f"{normalized_section}_{normalized_key}"
            if not restrict_to_section:
                metadata[normalized_key] = value
                section_keys.add(normalized_key)
            metadata[qualified_key] = value
            section_keys.add(qualified_key)

            if section == "sample_group" and normalized_key in {
                "name",
                "comments",
                "description",
            }:
                self._consumed_keys.add(normalized_key)

        primary_field = self.configuration.section_primary_field.get(section)
        if primary_field:
            normalized_primary = normalize_metadata_key(primary_field)
            if normalized_primary in section_values:
                metadata[section] = section_values[normalized_primary]
                section_keys.add(section)
        elif not restrict_to_section and section not in metadata:
            primary_aliases: Iterable[str] = alias_map.get(section, [])
            for alias in primary_aliases:
                alias_key = normalize_metadata_key(alias)
                if alias_key in metadata:
                    metadata[section] = metadata[alias_key]
                    section_keys.add(section)
                    break

        if section_keys:
            self._record_section_keys(section, section_keys)

    def _record_section_keys(self, section: str, keys: set[str]) -> None:
        if not keys:
            return
        recorded = self._section_keys.setdefault(section, set())
        recorded.update(keys)

    def _prune_inactive_platform_sections(self, metadata: Dict[str, Any]) -> None:
        for section in PLATFORM_SECTION_FIELDS:
            if section == self._active_platform_section:
                continue
            self._prune_section(metadata, section)

    def _prune_section(self, metadata: Dict[str, Any], section: str) -> None:
        keys = self._section_keys.pop(section, set())
        for key in keys:
            metadata.pop(key, None)


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
    "MetadataConfiguration",
    "MetadataConfigurationError",
    "VCFMetadataParser",
    "load_metadata_configuration",
    "normalize_metadata_key",
    "normalize_metadata_value",
]
