from __future__ import annotations

from django.db import models
from django.utils.translation import gettext_lazy as _


class PublishedQuerySet(models.QuerySet):
    def published(self) -> "PublishedQuerySet":
        return self.filter(is_published=True)


class ProductQuerySet(PublishedQuerySet):
    def featured(self) -> "ProductQuerySet":
        return self.filter(is_published=True, is_featured=True)


class Category(models.Model):
    name = models.CharField(max_length=255)
    slug = models.SlugField(unique=True)
    description = models.TextField(blank=True)
    images = models.JSONField(blank=True, null=True, help_text=_("Optional list of image data."))
    metadata = models.JSONField(blank=True, null=True)
    is_published = models.BooleanField(default=False)
    is_featured = models.BooleanField(default=False)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = PublishedQuerySet.as_manager()

    class Meta:
        ordering = ("name",)
        verbose_name_plural = "categories"

    def __str__(self) -> str:
        return self.name


class Product(models.Model):
    name = models.CharField(max_length=255)
    slug = models.SlugField(unique=True)
    description = models.TextField(blank=True)
    detail_html = models.TextField(blank=True, help_text=_("Rich product content or HTML."))
    images = models.JSONField(blank=True, null=True, help_text=_("List of image URLs or metadata."))
    metadata = models.JSONField(blank=True, null=True)
    is_published = models.BooleanField(default=False)
    is_featured = models.BooleanField(default=False)
    categories = models.ManyToManyField(Category, related_name="products", blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = ProductQuerySet.as_manager()

    class Meta:
        ordering = ("name",)

    def __str__(self) -> str:
        return self.name


class ProductDetail(models.Model):
    product = models.OneToOneField(
        Product,
        on_delete=models.CASCADE,
        related_name="detail",
    )
    headline = models.CharField(max_length=255, blank=True)
    detail_html = models.TextField(blank=True)
    images = models.JSONField(blank=True, null=True, help_text=_("Additional gallery or media."))
    metadata = models.JSONField(blank=True, null=True)
    is_published = models.BooleanField(default=False)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)

    objects = PublishedQuerySet.as_manager()

    class Meta:
        ordering = ("product__name",)
        verbose_name = "Product detail"
        verbose_name_plural = "Product details"

    def __str__(self) -> str:
        return self.headline or str(self.product)
