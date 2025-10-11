"""Reusable helpers for organization-aware views.

This module bundles small mixins that encapsulate common organization
resolution logic for class-based views.  In particular,
``OrganizationSampleGroupMixin`` exposes convenience methods to fetch the
requesting user's organization profile and to scope ``SampleGroup`` queryset
lookups to records created by that organization.  Views that include the mixin
can therefore rely on it for access control, ensuring that a user only sees
``SampleGroup`` instances owned by their organization without having to repeat
filtering logic.
"""

from __future__ import annotations

from django.db.models import QuerySet

from .models import SampleGroup


class OrganizationSampleGroupMixin:
    """Provide helpers to access organization-specific sample groups."""

    def get_organization_profile(self):
        """Return the requesting user's organization profile if available."""

        request = getattr(self, "request", None)
        if request is None:
            return None
        user = getattr(request, "user", None)
        if user is None:
            return None
        return getattr(user, "organization_profile", None)

    def get_owned_sample_groups(self) -> QuerySet[SampleGroup]:
        """Return sample groups owned by the requesting organization."""

        organization_profile = self.get_organization_profile()
        if organization_profile is None:
            return SampleGroup.objects.none()
        return SampleGroup.objects.filter(created_by=organization_profile)
