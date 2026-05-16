from django.contrib import admin

from .models import Category, Product, ProductDetail


class ProductDetailInline(admin.StackedInline):
    model = ProductDetail
    extra = 0
    classes = ["collapse"]


@admin.register(Category)
class CategoryAdmin(admin.ModelAdmin):
    list_display = ("name", "is_published", "is_featured", "created_at", "updated_at")
    list_filter = ("is_published", "is_featured")
    search_fields = ("name", "slug", "description")
    prepopulated_fields = {"slug": ("name",)}
    readonly_fields = ("created_at", "updated_at")


@admin.register(Product)
class ProductAdmin(admin.ModelAdmin):
    list_display = ("name", "is_published", "is_featured", "created_at", "updated_at")
    list_filter = ("is_published", "is_featured", "categories")
    search_fields = ("name", "slug", "description")
    prepopulated_fields = {"slug": ("name",)}
    filter_horizontal = ("categories",)
    readonly_fields = ("created_at", "updated_at")
    inlines = [ProductDetailInline]


@admin.register(ProductDetail)
class ProductDetailAdmin(admin.ModelAdmin):
    list_display = ("product", "headline", "is_published", "created_at", "updated_at")
    list_filter = ("is_published",)
    search_fields = ("headline", "product__name")
    readonly_fields = ("created_at", "updated_at")
