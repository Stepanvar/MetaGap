from django.contrib import admin

from .models import NewsArticle, NewsCategory


@admin.register(NewsCategory)
class NewsCategoryAdmin(admin.ModelAdmin):
    list_display = ("name", "slug")
    prepopulated_fields = {"slug": ("name",)}
    search_fields = ("name", "slug")
    ordering = ("name",)


@admin.register(NewsArticle)
class NewsArticleAdmin(admin.ModelAdmin):
    list_display = ("title", "author", "published_at", "featured", "updated_at")
    list_filter = ("featured", "published_at", "categories")
    search_fields = ("title", "slug", "body", "author__username", "author__email")
    prepopulated_fields = {"slug": ("title",)}
    ordering = ("-published_at", "-updated_at")
    filter_horizontal = ("categories",)
    list_select_related = ("author",)
