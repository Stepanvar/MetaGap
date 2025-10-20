"""Custom context processors for global template variables."""

from __future__ import annotations

from django.conf import settings


def branding(request):
    """Expose project branding information to every template."""

    return {
        "SITE_NAME": getattr(settings, "SITE_NAME", "MetaGaP"),
        "SITE_COLORS": getattr(settings, "SITE_COLORS", {}),
    }
