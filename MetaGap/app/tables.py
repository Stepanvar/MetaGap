import django_tables2 as tables
from django.db import models

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
            # For related fields (ForeignKey, OneToOneField) include their own fields.
            if field.is_relation and field.related_model is not None:
                if include_related and (field.many_to_one or field.one_to_one):
                    for rel_field in field.related_model._meta.get_fields():
                        if not rel_field.auto_created and isinstance(rel_field, models.Field):
                            col_name = f"{field.name}__{rel_field.name}"
                            columns[col_name] = tables.Column(accessor=f"{field.name}.{rel_field.name}")
                            meta_fields.append(col_name)
                else:
                    # If not expanding related fields, add the field as-is.
                    columns[field.name] = tables.Column()
                    meta_fields.append(field.name)
            else:
                # Standard (non-relation) field.
                columns[field.name] = tables.Column()
                meta_fields.append(field.name)

    # Create a dynamic Meta class.
    Meta = type("Meta", (), {
        "model": primary_model,
        "fields": meta_fields,
        "attrs": {"class": "table table-striped"}
    })

    # Attach Meta to the columns dict.
    columns["Meta"] = Meta

    # Create and return the table class.
    return type(table_name, (tables.Table,), columns)
