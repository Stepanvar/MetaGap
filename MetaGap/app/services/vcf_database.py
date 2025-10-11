"""Database persistence helpers for the VCF importer."""

from __future__ import annotations

import json
import os
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
        group_data, group_consumed, group_additional = self._extract_section_data(
            metadata, "sample_group", SampleGroup
        )

        consumed_keys = set(group_consumed)
        additional_metadata: Dict[str, Any] = {}
        parser_additional = metadata.get("additional_metadata")
        if isinstance(parser_additional, dict):
            for key, value in parser_additional.items():
                additional_metadata[key] = self._coerce_additional_value(value)
            consumed_keys.add("additional_metadata")
        if group_additional:
            additional_metadata.update(group_additional)

        fallback_name = os.path.splitext(os.path.basename(file_path))[0]
        name = group_data.pop("name", None) or metadata.get("name") or fallback_name
        comments = (
            group_data.pop("comments", None)
            or metadata.get("comments")
            or metadata.get("description")
        )

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
                    metadata, section, model_cls
                )

                consumed_keys.update(section_consumed)
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
    ) -> Tuple[Dict[str, Any], set[str], Optional[Dict[str, Any]]]:
        alias_map = self.metadata_field_aliases.get(section, {})
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

        primary_field = self.section_primary_field.get(section)
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
        except django_models.FieldDoesNotExist:
            return None

    def _find_metadata_key(
        self,
        metadata: Dict[str, Any],
        section: str,
        field_name: str,
        aliases: Iterable[str],
    ) -> Optional[str]:
        def _dedupe_append(collection: list[str], candidate: str) -> None:
            if candidate and candidate not in collection:
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

        for candidate in candidate_order:
            if candidate in metadata:
                return candidate

        normalized_lookup: Dict[str, str] = {}

        def _register_lookup(normalized_key: str, original_key: str) -> None:
            if normalized_key and normalized_key not in normalized_lookup:
                normalized_lookup[normalized_key] = original_key

        for key in metadata.keys():
            normalized_key = normalize_metadata_key(key)
            _register_lookup(normalized_key, key)
            collapsed_key = normalized_key.replace("_", "")
            _register_lookup(collapsed_key, key)

        for candidate in candidate_order:
            normalized_candidate = normalize_metadata_key(candidate)
            mapped_key = normalized_lookup.get(normalized_candidate)
            if mapped_key is not None:
                return mapped_key

            collapsed_candidate = normalized_candidate.replace("_", "")
            mapped_key = normalized_lookup.get(collapsed_candidate)
            if mapped_key is not None:
                return mapped_key

        return None

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
