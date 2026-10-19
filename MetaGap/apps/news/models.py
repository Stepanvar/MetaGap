from django.conf import settings
from django.db import models
from django.utils import timezone
from django.utils.translation import gettext_lazy as _


class NewsCategory(models.Model):
    """A taxonomy category for grouping news articles."""

    name = models.CharField(_("name"), max_length=200, unique=True)
    slug = models.SlugField(_("slug"), max_length=200, unique=True)

    class Meta:
        verbose_name = _("news category")
        verbose_name_plural = _("news categories")
        ordering = ("name",)

    def __str__(self) -> str:
        return self.name


class NewsArticleQuerySet(models.QuerySet):
    def published(self) -> "NewsArticleQuerySet":
        """Return articles with a publication date not in the future."""

        return self.filter(published_at__isnull=False, published_at__lte=timezone.now())

    def featured(self) -> "NewsArticleQuerySet":
        """Return articles flagged as featured."""

        return self.filter(featured=True)


class NewsArticleManager(models.Manager):
    def get_queryset(self) -> NewsArticleQuerySet:  # type: ignore[override]
        return NewsArticleQuerySet(self.model, using=self._db)

    def published(self) -> NewsArticleQuerySet:
        return self.get_queryset().published()

    def featured(self) -> NewsArticleQuerySet:
        return self.get_queryset().featured()


class NewsArticle(models.Model):
    """A news article for the MetaGap platform."""

    title = models.CharField(_("title"), max_length=255)
    slug = models.SlugField(_("slug"), max_length=255, unique=True)
    body = models.TextField(_("body"))
    author = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.CASCADE,
        related_name="news_articles",
        verbose_name=_("author"),
    )
    categories = models.ManyToManyField(
        NewsCategory,
        related_name="articles",
        verbose_name=_("categories"),
        blank=True,
    )
    hero_image = models.ImageField(
        _("hero image"), upload_to="news/hero_images/", blank=True, null=True
    )
    metadata = models.JSONField(_("metadata"), default=dict, blank=True)
    featured = models.BooleanField(_("featured"), default=False)
    published_at = models.DateTimeField(_("published at"), blank=True, null=True)
    updated_at = models.DateTimeField(_("updated at"), auto_now=True)

    objects = NewsArticleManager()

    class Meta:
        verbose_name = _("news article")
        verbose_name_plural = _("news articles")
        ordering = ("-published_at", "-updated_at")

    def __str__(self) -> str:
        return self.title
