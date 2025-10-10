from typing import Any, Dict, Iterable, List, Sequence

import django_tables2 as tables
from django.db import models
from django.utils.html import format_html

from .models import AlleleFrequency


NUMERIC_FIELD_TYPES = (
    models.AutoField,
    models.BigAutoField,
    models.SmallIntegerField,
    models.IntegerField,
    models.BigIntegerField,
    models.PositiveIntegerField,
    models.PositiveSmallIntegerField,
    models.FloatField,
    models.DecimalField,
    models.DurationField,
)


PRIORITY_COLUMN_DETAILS: Dict[str, Dict[str, Any]] = {
    "variant_id": {
        "abbr": "ID",
        "label": "Variant ID",
        "description": "Identifier assigned to the variant in the VCF (ID field).",
        "icon": "fa-fingerprint",
        "tooltip": "Unique identifier sourced from the VCF ID column.",
        "unit": None,
        "badge_class": "bg-primary text-white",
        "attrs": {
            "th": {"class": "priority-metric-header"},
            "td": {"class": "font-monospace text-nowrap"},
        },
    },
    "info__af": {
        "abbr": "AF",
        "label": "Allele frequency",
        "description": "Proportion of alternate allele observations",
        "icon": "fa-chart-column",
        "tooltip": "Allele frequency (AF) expressed as a proportion between 0 and 1.",
        "unit": "0–1",
        "badge_class": "bg-primary text-white",
        "attrs": {
            "th": {"class": "priority-metric-header"},
            "td": {"class": "priority-metric-cell"},
        },
    },
    "info__ac": {
        "abbr": "AC",
        "label": "Allele count",
        "description": "Total count of alternate alleles observed",
        "icon": "fa-hashtag",
        "tooltip": "Allele count (AC) representing the total number of alternate alleles observed.",
        "unit": "count",
        "badge_class": "bg-secondary text-white",
        "attrs": {
            "th": {"class": "priority-metric-header"},
            "td": {"class": "priority-metric-cell"},
        },
    },
    "info__dp": {
        "abbr": "DP",
        "label": "Depth",
        "description": "Total read depth at the locus",
        "icon": "fa-layer-group",
        "tooltip": "Read depth (DP) indicating the total number of reads covering the locus.",
        "unit": "reads",
        "badge_class": "bg-info text-dark",
        "attrs": {
            "th": {"class": "priority-metric-header"},
            "td": {"class": "priority-metric-cell"},
        },
    },
    "info__mq": {
        "abbr": "MQ",
        "label": "Mapping quality",
        "description": "Phred-scaled mapping quality",
        "icon": "fa-wave-square",
        "tooltip": "Mapping quality (MQ) summarising alignment confidence on a Phred scale.",
        "unit": "Phred",
        "badge_class": "bg-warning text-dark",
        "attrs": {
            "th": {"class": "priority-metric-header"},
            "td": {"class": "priority-metric-cell"},
        },
    },
}


def _is_numeric_field(field: models.Field) -> bool:
    """Return True when the provided model field stores numeric values."""
    return isinstance(field, NUMERIC_FIELD_TYPES)


def _metric_header(name: str, details: Dict[str, Any]) -> str:
    """Render a compact badge + description header for metric columns."""
    abbr = details.get("abbr", name.upper())
    description = details.get("description", "")
    unit = details.get("unit")
    icon_class = details.get("icon")

    description_text = description
    if unit:
        description_text = f"{description} ({unit})"

    icon_markup = ""
    if icon_class:
        icon_markup = format_html(
            '<i class="fa-solid {} text-primary"></i>', icon_class
        )

    return format_html(
        """
        <span class="d-inline-flex align-items-center gap-2">
            {icon}
            <span class="badge rounded-pill {badge_class} fw-semibold">{abbr}</span>
            <span class="text-muted small">{description}</span>
        </span>
        """,
        icon=icon_markup,
        badge_class=details.get("badge_class", "bg-primary text-white"),
        abbr=abbr,
        description=description_text,
    )


def _column_attrs(is_numeric: bool) -> Dict[str, Dict[str, str]]:
    """Generate django-tables2 column attrs with shared styling rules."""
    attrs: Dict[str, Dict[str, str]] = {"th": {"class": "fw-bold text-nowrap"}}
    if is_numeric:
        attrs["td"] = {"class": "text-end"}
    return attrs


def _merge_class(existing: str | None, new: str | None) -> str | None:
    """Combine CSS classes while preserving uniqueness."""

    if not existing:
        return new
    if not new:
        return existing
    combined = existing.split()
    for cls in new.split():
        if cls not in combined:
            combined.append(cls)
    return " ".join(combined)


def _merge_column_attrs(
    base: Dict[str, Dict[str, str]],
    extra: Dict[str, Dict[str, str]] | None,
) -> Dict[str, Dict[str, str]]:
    """Merge base column attributes with column-specific overrides."""

    if not extra:
        return base

    merged: Dict[str, Dict[str, str]] = {}
    for key, value in base.items():
        merged[key] = dict(value)

    for cell_type, attrs in extra.items():
        target = merged.setdefault(cell_type, {})
        for attr_name, attr_value in attrs.items():
            if attr_name == "class":
                target["class"] = _merge_class(target.get("class"), attr_value)
            else:
                target[attr_name] = attr_value
    return merged


def _build_column(
    name: str,
    *,
    accessor: str | None = None,
    is_numeric: bool = False,
) -> tables.Column:
    """Create a tables.Column respecting priority metric configuration."""

    base_attrs = _column_attrs(is_numeric)
    column_config = PRIORITY_COLUMN_DETAILS.get(name, {})
    attrs = _merge_column_attrs(base_attrs, column_config.get("attrs"))

    column_kwargs: Dict[str, Any] = {"attrs": attrs}
    if accessor is not None:
        column_kwargs["accessor"] = accessor

    tooltip = column_config.get("tooltip")
    if tooltip:
        attrs.setdefault("th", {})
        attrs["th"]["title"] = tooltip
        attrs["th"].setdefault("data-bs-toggle", "tooltip")

    if name in PRIORITY_COLUMN_DETAILS:
        column_kwargs["verbose_name"] = _metric_header(name, column_config)

    for key in ("verbose_name", "orderable", "linkify"):
        if key in column_config:
            column_kwargs[key] = column_config[key]

    return tables.Column(**column_kwargs)

def create_dynamic_table(
    primary_model,
    table_name: str = "DynamicTable",
    include_related: bool = True,
    priority_fields: Sequence[str] | None = None,
    exclude_fields: Iterable[str] | None = None,
):
    """
    Dynamically creates a django_tables2.Table subclass for the given primary_model.
    If include_related is True, for each forward relation (ForeignKey/OneToOneField),
    it adds columns for all of the related model's non-auto-created fields using a
    'relation__field' naming scheme.

    :param primary_model: The primary Django model class.
    :param table_name: Name of the generated table class.
    :param include_related: Whether to include fields from related models.
    :param priority_fields: Optional iterable of column names to position at the
        front of the table in the order provided.
    :param exclude_fields: Optional iterable of column names to omit from the
        generated table entirely.
    :return: A dynamically created subclass of tables.Table.
    """
    columns: Dict[str, tables.Column] = {}
    field_order: List[str] = []
    priority_list: List[str] = list(priority_fields or [])
    exclude_set = set(exclude_fields or [])

    def register_column(name: str, column: tables.Column) -> None:
        if name in exclude_set:
            return
        columns[name] = column
        field_order.append(name)

    # Loop over all fields in the primary model.
    for field in primary_model._meta.get_fields():
        # Skip auto-created reverse relations.
        if field.auto_created:
            continue

        if isinstance(field, models.Field):
            is_numeric = _is_numeric_field(field)
            # For related fields (ForeignKey, OneToOneField) include their own fields.
            if field.is_relation and field.related_model is not None:
                if include_related and (field.many_to_one or field.one_to_one):
                    for rel_field in field.related_model._meta.get_fields():
                        if not rel_field.auto_created and isinstance(rel_field, models.Field):
                            col_name = f"{field.name}__{rel_field.name}"
                            rel_is_numeric = _is_numeric_field(rel_field)
                            register_column(
                                col_name,
                                _build_column(
                                    col_name,
                                    accessor=f"{field.name}.{rel_field.name}",
                                    is_numeric=rel_is_numeric,
                                ),
                            )
                else:
                    # If not expanding related fields, add the field as-is.
                    register_column(
                        field.name,
                        _build_column(field.name, is_numeric=is_numeric),
                    )
            else:
                # Standard (non-relation) field.
                register_column(
                    field.name, _build_column(field.name, is_numeric=is_numeric)
                )

    seen = set()
    ordered_meta_fields: List[str] = []
    for field_name in priority_list:
        if field_name in columns and field_name not in seen:
            ordered_meta_fields.append(field_name)
            seen.add(field_name)

    for field_name in field_order:
        if field_name not in seen:
            ordered_meta_fields.append(field_name)
            seen.add(field_name)

    # Create a dynamic Meta class.
    Meta = type(
        "Meta",
        (),
        {
            "model": primary_model,
            "fields": ordered_meta_fields,
            "attrs": {
                "class": "table table-hover table-striped table-sm align-middle",
            },
        },
    )

    # Attach Meta to the columns dict.
    columns["Meta"] = Meta

    # Create and return the table class.
    return type(table_name, (tables.Table,), columns)


DEFAULT_ALLELE_PRIORITY_FIELDS: Sequence[str] = (
    "chrom",
    "pos",
    "ref",
    "alt",
    "qual",
    "filter",
    "info__af",
    "info__ac",
    "info__an",
    "info__dp",
    "info__mq",
)

DEFAULT_ALLELE_EXCLUDE_FIELDS: Sequence[str] = (
    "info__additional",
    "format__payload",
)


def build_allele_frequency_table(
    *,
    priority_extra: Sequence[str] | None = None,
    exclude_extra: Iterable[str] | None = None,
    table_name: str = "AlleleFrequencyTable",
    include_related: bool = True,
):
    """Create an :class:`~django_tables2.Table` for :class:`AlleleFrequency` records.

    Parameters allow callers to extend the shared priority and exclusion lists while
    keeping their defaults aligned across the application.
    """

    priority_fields: List[str] = list(DEFAULT_ALLELE_PRIORITY_FIELDS)
    if priority_extra:
        priority_fields.extend(priority_extra)

    exclude_fields = list(DEFAULT_ALLELE_EXCLUDE_FIELDS)
    if exclude_extra:
        exclude_fields.extend(exclude_extra)

    return create_dynamic_table(
        AlleleFrequency,
        table_name=table_name,
        include_related=include_related,
        priority_fields=priority_fields,
        exclude_fields=exclude_fields,
    )


PRIORITY_METRIC_SUMMARY: Sequence[Dict[str, Any]] = tuple(
    {
        "key": column,
        "abbr": details.get("abbr", column.upper()),
        "label": details.get("label", column.replace("_", " ").title()),
        "description": details.get("description", ""),
        "icon": details.get("icon"),
        "tooltip": details.get("tooltip"),
        "unit": details.get("unit"),
        "badge_class": details.get("badge_class", "bg-primary text-white"),
    }
    for column, details in PRIORITY_COLUMN_DETAILS.items()
)
