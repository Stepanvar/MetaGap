"""Simple content pages."""

from django.views.generic import TemplateView


class ContactView(TemplateView):
    template_name = "contact.html"


class AboutView(TemplateView):
    template_name = "about.html"


__all__ = ["ContactView", "AboutView"]
