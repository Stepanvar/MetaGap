"""Tests for context processors."""

from django.conf import settings
from django.test import SimpleTestCase, override_settings

from MetaGap.context_processors import branding


class BrandingContextProcessorTests(SimpleTestCase):
    """Tests for the branding context processor."""

    @override_settings(SITE_NAME="Custom", SITE_COLORS={"primary": "#fff"})
    def test_branding_uses_custom_settings(self):
        """The context processor should reflect custom branding settings."""

        self.assertEqual(
            branding(request=None),
            {"SITE_NAME": "Custom", "SITE_COLORS": {"primary": "#fff"}},
        )

    def test_branding_defaults_when_settings_missing(self):
        """The context processor should fall back to default branding values."""

        with override_settings():
            for attr in ("SITE_NAME", "SITE_COLORS"):
                if hasattr(settings, attr):
                    delattr(settings, attr)

            self.assertEqual(
                branding(request=None),
                {"SITE_NAME": "MetaGaP", "SITE_COLORS": {}},
            )
