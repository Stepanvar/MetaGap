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
        if not section:
            warning = (
                "Unsupported sequencing platform "
                f"{platform_value!r}; storing raw metadata only."
            )
            logger.warning("%s", warning)
            self.warnings.append(warning)
            section = "platform_independent"

        restrict = section != "platform_independent"
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
        def normalize_alias_key(candidate: str) -> str:
            stripped = candidate.rstrip("?!.,;:")
            return normalize_metadata_key(stripped)

        alias_map = self.configuration.field_aliases.get(section, {})
        normalized_section = normalize_metadata_key(section)
        section_values: Dict[str, Any] = {}

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
            value_key = normalized_key

            for field_name, aliases in alias_map.items():
                alias_candidates = [normalize_alias_key(alias) for alias in aliases]
                if normalized_key in alias_candidates:
                    normalized_key = normalize_metadata_key(field_name)
                    break

            section_values[normalized_key] = value
            qualified_key = f"{normalized_section}_{value_key}"
            if not restrict_to_section:
                metadata[normalized_key] = value
            metadata[qualified_key] = value

        primary_field = self.configuration.section_primary_field.get(section)
        if primary_field:
            normalized_primary = normalize_metadata_key(primary_field)
            if normalized_primary in section_values:
                metadata[section] = section_values[normalized_primary]
        elif not restrict_to_section and section not in metadata:
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
    "MetadataConfiguration",
    "MetadataConfigurationError",
    "VCFMetadataParser",
    "load_metadata_configuration",
    "normalize_metadata_key",
    "normalize_metadata_value",
]
