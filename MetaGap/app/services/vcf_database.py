"""Database persistence helpers for the VCF importer."""

from __future__ import annotations

import json
import os
import re
from decimal import Decimal, InvalidOperation
from typing import Any, Dict, Iterable, Optional, Tuple

from django.core.exceptions import ValidationError
from django.db import DataError, IntegrityError, models as django_models

from ..models import AlleleFrequency, Format, Info, SampleGroup
from .vcf_metadata import load_metadata_configuration, normalize_metadata_key


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


class VCFDatabaseWriter:
    """Handle creation of database records from parsed VCF content."""

    def __init__(self) -> None:
        configuration = load_metadata_configuration()
        self.metadata_field_aliases = configuration.field_aliases
        self.metadata_model_map = configuration.models
        self.section_primary_field = configuration.section_primary_field

    # ------------------------------------------------------------------
    # Sample group handling
    # ------------------------------------------------------------------
    def create_sample_group(
        self,
        metadata: Dict[str, Any],
        file_path: str,
        organization_profile: Any,
    ) -> SampleGroup:
        parser_consumed_raw = metadata.pop("_consumed_keys", None)
        parser_consumed = set(parser_consumed_raw or [])

        group_data, group_consumed, group_additional = self._extract_section_data(
            metadata, "sample_group", SampleGroup
        )

        consumed_keys = set(group_consumed)
        consumed_keys.update(parser_consumed)
        section_overrides = metadata.pop("_section_additional_fields", None)
        if isinstance(section_overrides, dict):
            consumed_keys.add("_section_additional_fields")
        additional_metadata: Dict[str, Any] = {}
        parser_additional = metadata.pop("additional_metadata", None)
        if isinstance(parser_additional, dict):
            for key, value in parser_additional.items():
                additional_metadata[key] = self._coerce_additional_value(value)
            consumed_keys.add("additional_metadata")
        if group_additional:
            additional_metadata.update(group_additional)

        self._consume_additional_metadata_buffer(
            metadata, additional_metadata, consumed_keys
        )

        fallback_name = os.path.splitext(os.path.basename(file_path))[0]
        name_candidate = group_data.pop("name", None)
        metadata_name = metadata.get("name")
        if name_candidate:
            name = name_candidate
        elif metadata_name:
            name = metadata_name
            consumed_keys.add("name")
        else:
            name = fallback_name
            if "name" in metadata:
                consumed_keys.add("name")

        comments_candidate = group_data.pop("comments", None)
        metadata_comments = metadata.get("comments")
        metadata_description = metadata.get("description")
        if comments_candidate:
            comments = comments_candidate
        elif metadata_comments:
            comments = metadata_comments
            consumed_keys.add("comments")
        elif metadata_description:
            comments = metadata_description
            consumed_keys.add("description")
        else:
            comments = None
            if "comments" in metadata:
                consumed_keys.add("comments")
            if "description" in metadata:
                consumed_keys.add("description")

        try:
            sample_group = SampleGroup.objects.create(
                name=name,
                created_by=organization_profile,
                comments=comments,
                **group_data,
            )

            update_fields: list[str] = []
            for section, model_cls in self.metadata_model_map.items():
                section_data, section_consumed, additional = self._extract_section_data(
                    metadata, section, model_cls, skip_keys=consumed_keys
                )

                consumed_keys.update(section_consumed)
                self._consume_additional_metadata_buffer(
                    metadata, additional_metadata, consumed_keys
                )
                if not section_data and additional is None:
                    continue

                payload = {
                    key: value for key, value in section_data.items() if value is not None
                }
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

            self._consume_additional_metadata_buffer(
                metadata, additional_metadata, consumed_keys
            )

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
        except IntegrityError as exc:
            raise ValidationError(
                f"A dataset named '{name}' already exists. Please choose a different name."
            ) from exc
        except DataError as exc:
            raise ValidationError(
                "One or more metadata values are out of range. Please review your dataset metadata."
            ) from exc

    def _extract_section_data(
        self,
        metadata: Dict[str, Any],
        section: str,
        model_cls: Any,
        skip_keys: Optional[Iterable[str]] = None,
    ) -> Tuple[Dict[str, Any], set[str], Optional[Dict[str, Any]]]:
        alias_map = self.metadata_field_aliases.get(section, {})
        section_data: Dict[str, Any] = {}
        consumed: set[str] = set()
        skip_set = set(skip_keys or [])
        normalized_section = normalize_metadata_key(section)
        drop_fields, instruction_consumed = self._collect_section_additional_fields(
            metadata, section, alias_map, skip_set
        )
        if drop_fields:
            section_overrides = metadata.setdefault("_section_additional_fields", {})
            if isinstance(section_overrides, dict):
                section_overrides[section] = sorted(drop_fields)
        if instruction_consumed:
            consumed.update(instruction_consumed)

        for field_name, aliases in alias_map.items():
            model_field = self._get_model_field(model_cls, field_name)
            if model_field is None:
                continue
            normalized_field = normalize_metadata_key(field_name)
            if normalized_field in drop_fields:
                key = self._find_metadata_key(
                    metadata,
                    section,
                    field_name,
                    aliases,
                    skip_keys=skip_set | consumed,
                )
                if key is None:
                    continue
                raw_value = metadata[key]
                self._store_dropped_field_value(
                    metadata, section, normalized_field, raw_value
                )

                consumed.add(key)
                consumed.add(field_name)
                consumed.add(normalized_field)

                normalized_key = normalize_metadata_key(key)
                for prefix in (
                    f"{normalized_section}_",
                    f"{normalized_section}.",
                    f"{normalized_section}-",
                ):
                    if normalized_key.startswith(prefix):
                        suffix = normalized_key[len(prefix) :]
                        if suffix:
                            consumed.add(suffix)
                        break

                equivalent_keys = self._resolve_equivalent_metadata_keys(
                    metadata, section, field_name, aliases, key
                )
                if equivalent_keys:
                    consumed.update(equivalent_keys)
                continue
            key = self._find_metadata_key(
                metadata,
                section,
                field_name,
                aliases,
                skip_keys=skip_set | consumed,
            )
            if key is None:
                continue
            raw_value = metadata[key]
            if (
                section == "sample_group"
                and field_name == "sequencing_platform"
                and str(raw_value) not in SampleGroup.SequencingPlatform.values
            ):
                consumed.add(key)
                continue
            section_data[field_name] = self._coerce_model_value(
                model_field, raw_value
            )

            # Consume the original key and normalized variants
            consumed.add(key)
            consumed.add(field_name)  # keep for reverse lookups if you rely on it elsewhere
            consumed.add(normalize_metadata_key(field_name))

            normalized_key = normalize_metadata_key(key)
            for prefix in (
                f"{normalized_section}_",
                f"{normalized_section}.",
                f"{normalized_section}-",
            ):
                if normalized_key.startswith(prefix):
                    suffix = normalized_key[len(prefix) :]
                    if suffix:
                        consumed.add(suffix)
                    break

            # Also consume equivalent/aliased keys (from codex branch)
            equivalent_keys = self._resolve_equivalent_metadata_keys(
                metadata, section, field_name, aliases, key
            )
            if equivalent_keys:
                consumed.update(equivalent_keys)

        primary_field = self.section_primary_field.get(section)
        if (
            primary_field
            and primary_field not in section_data
            and section in metadata
            and section not in skip_set
        ):
            model_field = self._get_model_field(model_cls, primary_field)
            if model_field is not None:
                section_data[primary_field] = self._coerce_model_value(
                    model_field, metadata[section]
                )
                consumed.add(section)

        additional, additional_consumed = self._build_additional_payload(
            metadata, section, consumed | skip_set
        )
        consumed.update(additional_consumed)
        return section_data, consumed, additional

    def _collect_section_additional_fields(
        self,
        metadata: Dict[str, Any],
        section: str,
        alias_map: Dict[str, Iterable[str]],
        skip_set: set[str],
    ) -> Tuple[set[str], set[str]]:
        normalized_section = normalize_metadata_key(section)
        candidate_keys = [
            "additional_fields",
            "additionalfields",
            f"{section}_additional_fields",
            f"{section}_additionalfields",
            f"{normalized_section}_additional_fields",
            f"{normalized_section}_additionalfields",
        ]
        lookup = self._build_normalized_lookup(metadata)
        drop_fields: set[str] = set()
        consumed: set[str] = set()
        seen_keys: set[str] = set()

        for candidate in candidate_keys:
            normalized_candidate = normalize_metadata_key(candidate)
            if not normalized_candidate:
                continue
            for key in lookup.get(normalized_candidate, []):
                if key in skip_set or key in seen_keys:
                    continue
                seen_keys.add(key)
                consumed.add(key)
                raw_value = metadata.get(key)
                entries = self._parse_additional_field_entries(raw_value)
                for entry in entries:
                    canonical = self._resolve_additional_field_name(
                        entry, alias_map, normalized_section
                    )
                    if canonical:
                        drop_fields.add(canonical)

        if (
            section == "library_construction"
            and self._is_metadata_file_style(metadata)
        ):
            metadata_keys = self._build_normalized_lookup(metadata)
            target_keys = {
                normalize_metadata_key("pcr_cycles"),
                normalize_metadata_key(f"{section}_pcr_cycles"),
            }
            if any(metadata_keys.get(target) for target in target_keys):
                drop_fields.add("pcr_cycles")

        return drop_fields, consumed

    def _consume_additional_metadata_buffer(
        self,
        metadata: Dict[str, Any],
        additional_metadata: Dict[str, Any],
        consumed_keys: set[str],
    ) -> None:
        buffer = metadata.pop("additional_metadata", None)
        if not isinstance(buffer, dict):
            return

        for key, value in buffer.items():
            coerced = self._coerce_additional_value(value)
            existing = additional_metadata.get(key)
            if isinstance(existing, dict) and isinstance(coerced, dict):
                additional_metadata[key] = {**existing, **coerced}
            else:
                additional_metadata[key] = coerced

        consumed_keys.add("additional_metadata")

    def _is_metadata_file_style(self, metadata: Dict[str, Any]) -> bool:
        sentinel_keys = {
            "sample_group_source_lab",
            "sample_group_contact_email",
            "sample_group_contact_phone",
            "sample_group_total_samples",
        }
        metadata_keys = {normalize_metadata_key(key) for key in metadata.keys()}
        return any(key in metadata_keys for key in sentinel_keys)

    def _parse_additional_field_entries(self, raw_value: Any) -> set[str]:
        if raw_value is None:
            return set()

        entries: Iterable[Any]
        if isinstance(raw_value, str):
            candidate = raw_value.strip()
            if not candidate:
                return set()
            if (
                len(candidate) >= 2
                and candidate[0] in {"[", "{"}
                and candidate[-1] in {"]", "}"}
            ):
                try:
                    loaded = json.loads(candidate)
                except (TypeError, ValueError):
                    loaded = None
                if isinstance(loaded, dict):
                    entries = list(loaded.keys())
                elif isinstance(loaded, (list, tuple, set)):
                    entries = list(loaded)
                else:
                    entries = [candidate]
            else:
                parts = [segment.strip() for segment in re.split(r"[|,;]", candidate)]
                entries = [segment for segment in parts if segment]
        elif isinstance(raw_value, dict):
            entries = list(raw_value.keys())
        elif isinstance(raw_value, (list, tuple, set)):
            entries = list(raw_value)
        else:
            entries = [raw_value]

        normalized: set[str] = set()
        for entry in entries:
            normalized_entry = normalize_metadata_key(entry)
            if normalized_entry:
                normalized.add(normalized_entry)
        return normalized

    def _resolve_additional_field_name(
        self,
        normalized_entry: str,
        alias_map: Dict[str, Iterable[str]],
        normalized_section: str,
    ) -> Optional[str]:
        if not normalized_entry:
            return None

        candidates: set[str] = {normalized_entry, normalized_entry.replace("_", "")}
        collapsed_section = normalized_section.replace("_", "")

        for prefix in (normalized_section, collapsed_section):
            if not prefix:
                continue
            if normalized_entry.startswith(prefix):
                remainder = normalized_entry[len(prefix) :].lstrip("_.-")
                if remainder:
                    candidates.add(remainder)
                    candidates.add(remainder.replace("_", ""))

        expanded_candidates = set(candidates)
        for candidate in candidates:
            collapsed = candidate.replace("_", "")
            if collapsed:
                expanded_candidates.add(collapsed)

        for field_name, aliases in alias_map.items():
            canonical = normalize_metadata_key(field_name)
            alias_candidates = {canonical, canonical.replace("_", "")}
            for alias in aliases:
                normalized_alias = normalize_metadata_key(alias)
                if normalized_alias:
                    alias_candidates.add(normalized_alias)
                    alias_candidates.add(normalized_alias.replace("_", ""))
            if expanded_candidates & alias_candidates:
                return canonical

        return None

    def _store_dropped_field_value(
        self,
        metadata: Dict[str, Any],
        section: str,
        canonical_field: str,
        value: Any,
    ) -> None:
        normalized_section = normalize_metadata_key(section)
        if canonical_field.startswith(f"{normalized_section}_"):
            qualified_key = canonical_field
        else:
            qualified_key = f"{normalized_section}_{canonical_field}"

        additional = metadata.setdefault("additional_metadata", {})
        if not isinstance(additional, dict):
            additional = {}
            metadata["additional_metadata"] = additional

        additional.setdefault(qualified_key, value)

    @staticmethod
    def _get_model_field(model_cls: Any, field_name: str) -> Optional[django_models.Field]:
        try:
            return model_cls._meta.get_field(field_name)
        except django_models.FieldDoesNotExist:
            return None

    def _find_metadata_key(
        self,
        metadata: Dict[str, Any],
        section: str,
        field_name: str,
        aliases: Iterable[str],
        *,
        skip_keys: Optional[Iterable[str]] = None,
    ) -> Optional[str]:
        skip_exact: set[str] = set()
        skip_normalized: set[str] = set()

        for entry in skip_keys or []:
            if entry is None:
                continue
            text = str(entry)
            if not text:
                continue
            skip_exact.add(text)
            normalized = normalize_metadata_key(text)
            if normalized:
                skip_normalized.add(normalized)
                collapsed = normalized.replace("_", "")
                if collapsed:
                    skip_normalized.add(collapsed)

        combined_skip = skip_exact | skip_normalized

        candidate_order = self._enumerate_metadata_candidates(
            section,
            field_name,
            aliases,
            skip_keys=combined_skip,
        )

        for candidate in candidate_order:
            if candidate in skip_exact:
                continue
            if candidate in metadata:
                return candidate

        normalized_lookup = self._build_normalized_lookup(metadata)

        for candidate in candidate_order:
            normalized_candidate = normalize_metadata_key(candidate)
            if not normalized_candidate:
                continue
            if normalized_candidate in skip_normalized:
                continue

            mapped_keys = normalized_lookup.get(normalized_candidate)
            if mapped_keys:
                for mapped_key in mapped_keys:
                    if mapped_key in skip_exact:
                        continue
                    normalized_mapped = normalize_metadata_key(mapped_key)
                    if normalized_mapped and normalized_mapped in skip_normalized:
                        continue
                    return mapped_key

            collapsed_candidate = normalized_candidate.replace("_", "")
            if collapsed_candidate != normalized_candidate:
                if collapsed_candidate in skip_normalized:
                    continue
                mapped_keys = normalized_lookup.get(collapsed_candidate)
                if mapped_keys:
                    for mapped_key in mapped_keys:
                        if mapped_key in skip_exact:
                            continue
                        normalized_mapped = normalize_metadata_key(mapped_key)
                        if normalized_mapped and normalized_mapped in skip_normalized:
                            continue
                        return mapped_key

        return None

    def _enumerate_metadata_candidates(
        self,
        section: str,
        field_name: str,
        aliases: Iterable[str],
        skip_keys: Optional[Iterable[str]] = None,
    ) -> list[str]:
        skip_set: set[str] = set()
        if skip_keys:
            for entry in skip_keys:
                if entry is None:
                    continue
                text = str(entry)
                if not text:
                    continue
                skip_set.add(text)

        def _dedupe_append(collection: list[str], candidate: str) -> None:
            if candidate and candidate not in collection and candidate not in skip_set:
                collection.append(candidate)

        normalized_section = normalize_metadata_key(section)
        normalized_field = normalize_metadata_key(field_name)

        alias_candidates: list[str] = []
        for candidate in (str(field_name), normalized_field):
            _dedupe_append(alias_candidates, candidate)

        for alias in aliases:
            alias_text = str(alias)
            normalized_alias = normalize_metadata_key(alias_text)
            _dedupe_append(alias_candidates, alias_text)
            _dedupe_append(alias_candidates, normalized_alias)

        candidate_order: list[str] = []

        section_variants: list[str] = []
        for section_variant in (section, normalized_section):
            if section_variant:
                _dedupe_append(section_variants, section_variant)

        collapsed_variants: list[str] = []
        for section_variant in section_variants:
            collapsed = section_variant.replace("_", "")
            if collapsed and collapsed != section_variant:
                _dedupe_append(collapsed_variants, collapsed)

        section_variants.extend(collapsed_variants)

        for section_variant in section_variants:
            for alias_candidate in alias_candidates:
                _dedupe_append(
                    candidate_order, f"{section_variant}_{alias_candidate}"
                )

        for alias_candidate in alias_candidates:
            _dedupe_append(candidate_order, alias_candidate)

        return candidate_order

    @staticmethod
    def _build_normalized_lookup(metadata: Dict[str, Any]) -> Dict[str, list[str]]:
        lookup: Dict[str, list[str]] = {}

        def _append_lookup(key: str, original_key: str) -> None:
            if not key:
                return
            bucket = lookup.setdefault(key, [])
            if original_key not in bucket:
                bucket.append(original_key)

        for original_key in metadata.keys():
            normalized_key = normalize_metadata_key(original_key)
            _append_lookup(normalized_key, original_key)
            collapsed_key = normalized_key.replace("_", "")
            if collapsed_key != normalized_key:
                _append_lookup(collapsed_key, original_key)

        return lookup

    def _resolve_equivalent_metadata_keys(
        self,
        metadata: Dict[str, Any],
        section: str,
        field_name: str,
        aliases: Iterable[str],
        matched_key: str,
    ) -> set[str]:
        equivalents: set[str] = {matched_key}
        candidate_order = self._enumerate_metadata_candidates(section, field_name, aliases)
        normalized_lookup = self._build_normalized_lookup(metadata)

        for candidate in candidate_order:
            if candidate in metadata:
                equivalents.add(candidate)

        for candidate in candidate_order:
            normalized_candidate = normalize_metadata_key(candidate)
            if not normalized_candidate:
                continue

            for mapped_key in normalized_lookup.get(normalized_candidate, []):
                equivalents.add(mapped_key)

            collapsed_candidate = normalized_candidate.replace("_", "")
            if collapsed_candidate != normalized_candidate:
                for mapped_key in normalized_lookup.get(collapsed_candidate, []):
                    equivalents.add(mapped_key)

        return equivalents

    def _coerce_model_value(
        self, field: django_models.Field, value: Any
    ) -> Any:
        if isinstance(field, (django_models.CharField, django_models.TextField)):
            return str(value)
        if isinstance(field, django_models.BooleanField):
            if isinstance(value, str):
                normalized = value.strip().lower()
                return normalized in {"true", "1", "yes"}
            return bool(value)
        if isinstance(field, django_models.IntegerField):
            try:
                return int(value)
            except (TypeError, ValueError):
                return None
        if isinstance(field, django_models.FloatField):
            try:
                return float(value)
            except (TypeError, ValueError):
                return None
        if isinstance(field, django_models.JSONField):
            if isinstance(value, str):
                try:
                    return json.loads(value)
                except json.JSONDecodeError:
                    pass
            return value
        return value

    def _build_additional_payload(
        self,
        metadata: Dict[str, Any],
        section: str,
        consumed: Iterable[str],
    ) -> Tuple[Optional[Dict[str, Any]], set[str]]:
        consumed_set = set(consumed)
        normalized_section = normalize_metadata_key(section)
        section_variants = {section, normalized_section}
        prefixes: list[tuple[str, str]] = []
        for variant in filter(None, section_variants):
            lower_variant = variant.lower()
            prefixes.extend(
                [
                    (f"{variant}_", f"{lower_variant}_"),
                    (f"{variant}.", f"{lower_variant}."),
                    (f"{variant}-", f"{lower_variant}-"),
                ]
            )

        normalized_consumed = {
            normalize_metadata_key(entry) for entry in consumed_set if entry
        }
        consumed_suffixes: set[str] = set()
        for entry in consumed_set:
            normalized_entry = normalize_metadata_key(entry)
            for prefix in (
                f"{normalized_section}_",
                f"{normalized_section}.",
                f"{normalized_section}-",
            ):
                if normalized_entry.startswith(prefix):
                    suffix = normalized_entry[len(prefix) :]
                    if suffix:
                        consumed_suffixes.add(suffix)
                    break

        additional: Dict[str, Any] = {}
        additional_consumed: set[str] = set()

        for key, value in metadata.items():
            if key in consumed_set:
                continue
            normalized_key = normalize_metadata_key(key)
            if normalized_key in normalized_consumed or normalized_key in consumed_suffixes:
                additional_consumed.add(key)
                continue
            for prefix, prefix_lower in prefixes:
                if key.lower().startswith(prefix_lower):
                    trimmed = key[len(prefix) :]
                    if not trimmed:
                        additional_consumed.add(key)
                        break
                    normalized_trimmed = normalize_metadata_key(trimmed)
                    if (
                        normalized_trimmed in normalized_consumed
                        or normalized_trimmed in consumed_suffixes
                    ):
                        additional_consumed.add(key)
                        break
                    if trimmed not in additional:
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

    # ------------------------------------------------------------------
    # Variant record creation
    # ------------------------------------------------------------------
    def create_info_instance(self, info: Any) -> Optional[Info]:
        if not info:
            return None

        mapped_fields: Dict[str, Any] = {}
        additional: Dict[str, Any] = {}

        for key, value in dict(info or {}).items():
            normalized = normalize_metadata_key(key)
            if normalized in INFO_FIELD_MAP:
                mapped_key, field_type = INFO_FIELD_MAP[normalized]
                coerced = self._coerce_info_value(value, field_type)
                mapped_fields[mapped_key] = coerced
            else:
                normalized_value = self._normalize_additional_info_value(value)
                if normalized_value is not None:
                    additional[normalized] = normalized_value

        additional_payload = additional or None

        if not mapped_fields and additional_payload is None:
            return None

        return Info.objects.create(**mapped_fields, additional=additional_payload)

    def create_format_instance(
        self, samples: Any
    ) -> Tuple[Optional[Format], Optional[str]]:
        if not samples:
            return None, None

        sample_name: Optional[str] = None
        sample_data: Optional[Any] = None
        for name, data in samples.items():
            sample_name = name
            sample_data = data
            break

        if not sample_data:
            return None, sample_name

        structured: Dict[str, Any] = {}
        additional: Dict[str, Any] = {}

        for key in getattr(sample_data, "keys", lambda: sample_data.keys())():
            normalized = normalize_metadata_key(key)
            if normalized == "gt":
                serialized = self.serialize_genotype(sample_data, key)
            else:
                serialized = self.stringify(sample_data[key])

            mapped_field = FORMAT_FIELD_MAP.get(normalized)
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
        return cls.stringify(normalized)

    @classmethod
    def _normalize_info_scalar(cls, value: Any) -> Any:
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            normalized_items = [
                cls._normalize_info_scalar(item)
                for item in value
                if item not in (None, "")
            ]
            flattened = [item for item in normalized_items if item not in (None, "")]
            if not flattened:
                return None
            return flattened
        if isinstance(value, str):
            stripped = value.strip()
            if not stripped or stripped in INFO_PLACEHOLDER_VALUES:
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
            normalized_list = [
                cls._normalize_additional_info_value(item)
                for item in value
                if item not in (None, "")
            ]
            return normalized_list or None
        if value in (None, "", "."):
            return None
        if isinstance(value, (int, float, bool)):
            return value
        return str(value)

    @staticmethod
    def serialize_alt(alts: Optional[Iterable[str]]) -> str:
        return ",".join(alts or [])

    @staticmethod
    def serialize_filter(filter_field: Any) -> Optional[str]:
        if not filter_field:
            return None
        if isinstance(filter_field, str):
            return filter_field
        values = filter_field.keys() if hasattr(filter_field, "keys") else filter_field
        serialized = [str(value) for value in values]
        return ",".join(serialized) if serialized else None

    @staticmethod
    def serialize_genotype(sample_data: Any, key: str) -> Optional[str]:
        value = sample_data[key]
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            separator = "|" if getattr(sample_data, "phased", False) else "/"
            return separator.join(str(item) for item in value)
        return str(value)

    @staticmethod
    def stringify(value: Any) -> Optional[str]:
        if value is None:
            return None
        if isinstance(value, (list, tuple)):
            return ",".join(str(item) for item in value)
        return str(value)

    def create_allele_frequency(
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
        variant_label = f"{chrom}:{pos} {ref}>{alt}"

        try:
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
        except IntegrityError as exc:
            raise ValidationError(
                f"The variant {variant_label} is already present in this dataset."
            ) from exc
        except DataError as exc:
            raise ValidationError(
                f"Variant {variant_label} contains out-of-range or invalid values."
            ) from exc


__all__ = [
    "FORMAT_FIELD_MAP",
    "INFO_FIELD_MAP",
    "INFO_FIELD_FLOAT",
    "INFO_FIELD_INT",
    "INFO_FIELD_STRING",
    "VCFDatabaseWriter",
]
