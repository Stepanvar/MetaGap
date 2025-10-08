from typing import Dict

import django_tables2 as tables
from django.db import models


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


def _is_numeric_field(field: models.Field) -> bool:
    """Return True when the provided model field stores numeric values."""
    return isinstance(field, NUMERIC_FIELD_TYPES)


def _column_attrs(is_numeric: bool) -> Dict[str, Dict[str, str]]:
    """Generate django-tables2 column attrs with shared styling rules."""
    attrs: Dict[str, Dict[str, str]] = {"th": {"class": "fw-bold"}}
    if is_numeric:
        attrs["td"] = {"class": "text-end"}
    return attrs

def create_dynamic_table(primary_model, table_name="DynamicTable", include_related=True):
    """
    Dynamically creates a django_tables2.Table subclass for the given primary_model.
    If include_related is True, for each forward relation (ForeignKey/OneToOneField),
    it adds columns for all of the related model's non-auto-created fields using a
    'relation__field' naming scheme.

    :param primary_model: The primary Django model class.
    :param table_name: Name of the generated table class.
    :param include_related: Whether to include fields from related models.
    :return: A dynamically created subclass of tables.Table.
    """
    columns = {}
    meta_fields = []

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
                            columns[col_name] = tables.Column(
                                accessor=f"{field.name}.{rel_field.name}",
                                attrs=_column_attrs(rel_is_numeric),
                            )
                            meta_fields.append(col_name)
                else:
                    # If not expanding related fields, add the field as-is.
                    columns[field.name] = tables.Column(attrs=_column_attrs(is_numeric))
                    meta_fields.append(field.name)
            else:
                # Standard (non-relation) field.
                columns[field.name] = tables.Column(attrs=_column_attrs(is_numeric))
                meta_fields.append(field.name)

    # Create a dynamic Meta class.
    Meta = type("Meta", (), {
        "model": primary_model,
        "fields": meta_fields,
        "attrs": {"class": "table table-hover table-sm"}
    })

    # Attach Meta to the columns dict.
    columns["Meta"] = Meta

    # Create and return the table class.
    return type(table_name, (tables.Table,), columns)
