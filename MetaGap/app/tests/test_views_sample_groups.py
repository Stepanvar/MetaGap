"""Tests focused on core sample group CRUD views and mixins."""

from __future__ import annotations

from types import SimpleNamespace
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import RequestFactory, TestCase
from django.urls import NoReverseMatch, reverse

from ..mixins import OrganizationSampleGroupMixin
from ..models import (
    InputQuality,
    ReferenceGenomeBuild,
    SampleGroup,
)
from ..views import (
    SampleGroupDeleteView,
    SampleGroupDetailView,
    SampleGroupUpdateView,
)
from .sample_group_test_mixins import SampleGroupViewMatrixMixin


class OrganizationSampleGroupMixinTests(TestCase):
    """Validate helper behaviour for organization-aware mixins."""

    def setUp(self) -> None:
        super().setUp()
        self.factory = RequestFactory()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="mixin_user",
            password="secure-pass",
            email="mixin@example.com",
        )
        self.other_user = User.objects.create_user(
            username="mixin_other",
            password="secure-pass",
            email="other@example.com",
        )

    def _build_view(self, request):
        class DummyView(OrganizationSampleGroupMixin):
            def __init__(self, request):
                self.request = request

        return DummyView(request)

    def test_get_organization_profile_returns_user_profile(self) -> None:
        request = self.factory.get("/profile")
        request.user = self.user

        mixin = self._build_view(request)

        self.assertEqual(
            mixin.get_organization_profile(), self.user.organization_profile
        )

    def test_get_owned_sample_groups_returns_empty_without_profile(self) -> None:
        request = self.factory.get("/profile")
        request.user = SimpleNamespace()

        mixin = self._build_view(request)

        self.assertEqual(list(mixin.get_owned_sample_groups()), [])

    def test_get_owned_sample_groups_filters_by_profile(self) -> None:
        SampleGroup.objects.create(
            name="Owned", created_by=self.user.organization_profile
        )
        SampleGroup.objects.create(
            name="Other", created_by=self.other_user.organization_profile
        )

        request = self.factory.get("/profile")
        request.user = self.user

        mixin = self._build_view(request)

        queryset = mixin.get_owned_sample_groups().order_by("name")
        self.assertQuerySetEqual(
            queryset,
            SampleGroup.objects.filter(
                created_by=self.user.organization_profile
            ).order_by("name"),
            transform=lambda group: group,
        )


class SampleGroupDetailViewTests(SampleGroupViewMatrixMixin, TestCase):
    """Exercise the sample group detail view and its access controls."""

    def detail_url(self) -> str:
        for name in (
            "sample_group_detail",
            "profile_sample_group_detail",
            "sample-group-detail",
        ):
            try:
                return reverse(name, args=[self.sample_group.pk])
            except NoReverseMatch:
                continue
        self.fail("Sample group detail route is not configured")

    def test_owner_can_view_detail(self) -> None:
        self.client.force_login(self.owner)

        response = self.client.get(self.detail_url())

        self.assertEqual(response.status_code, 200)
        context_object = response.context.get("sample_group") or response.context.get(
            "object"
        )
        self.assertEqual(context_object, self.sample_group)

    def test_forbids_access_for_non_owner(self) -> None:
        self.client.force_login(self.other_user)

        response = self.client.get(self.detail_url())

        self.assertIn(response.status_code, {403, 404})

    def test_metadata_sections_and_variant_table_present(self) -> None:
        self.client.force_login(self.owner)

        response = self.client.get(self.detail_url())

        metadata_sections = response.context.get("metadata_sections")
        if metadata_sections is None:
            metadata_sections = response.context.get("metadata")
        self.assertIsNotNone(metadata_sections)

        if isinstance(metadata_sections, dict):
            sections = metadata_sections.values()
        else:
            sections = metadata_sections

        metadata_values = {
            str(value)
            for section in sections
            if isinstance(section, dict)
            for value in section.values()
            if value not in (None, "", [])
        }

        for expected_value in {
            self.sample_group.name,
            self.sample_group.source_lab,
            self.sample_group.contact_email,
            self.sample_group.reference_genome_build.build_name,
            self.sample_group.genome_complexity.size,
        }:
            self.assertIn(expected_value, metadata_values)

        variant_table = response.context.get("variant_table")
        if variant_table is None:
            variant_table = response.context.get("table")
        self.assertIsNotNone(variant_table)

        header_names = [column.name for column in variant_table.columns]
        expected_prefix = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "qual",
            "filter",
            "info__af",
            "info__ac",
            "info__an",
            "info__dp",
            "info__mq",
            "variant_id",
        ]
        self.assertGreaterEqual(len(header_names), len(expected_prefix))
        self.assertEqual(header_names[: len(expected_prefix)], expected_prefix)
        self.assertNotIn("sample_group", header_names)
        self.assertNotIn("format__payload", header_names)

        table_html = variant_table.as_html(response.wsgi_request)
        self.assertIn(self.allele.variant_id, table_html)
        self.assertIn(str(self.allele.pos), table_html)
        self.assertIn(self.allele.ref, table_html)
        self.assertIn(self.allele.alt, table_html)
        self.assertIn(self.allele.info.af, table_html)
        self.assertIn(self.allele.info.ac, table_html)
        self.assertIn(self.allele.info.an, table_html)
        self.assertIn(self.allele.info.dp, table_html)
        self.assertIn(self.allele.info.mq, table_html)

    def test_detail_view_get_queryset_uses_mixin_helper(self) -> None:
        request = RequestFactory().get(self.detail_url())
        request.user = self.owner
        view = SampleGroupDetailView()
        view.setup(request, pk=self.sample_group.pk)

        base_queryset = SampleGroup.objects.filter(pk=self.sample_group.pk)

        with patch.object(
            SampleGroupDetailView, "get_owned_sample_groups", return_value=base_queryset
        ) as mock_groups:
            result = view.get_queryset()

        mock_groups.assert_called_once_with()
        expected = (
            base_queryset.select_related(
                "created_by",
                "created_by__user",
                "reference_genome_build",
                "genome_complexity",
                "sample_origin",
                "material_type",
                "library_construction",
                "illumina_seq",
                "ont_seq",
                "pacbio_seq",
                "iontorrent_seq",
                "bioinfo_alignment",
                "bioinfo_variant_calling",
                "bioinfo_post_proc",
                "input_quality",
            ).prefetch_related(
                "allele_frequencies",
                "allele_frequencies__info",
                "allele_frequencies__format",
            )
        )
        self.assertEqual(str(result.query), str(expected.query))


class SampleGroupUpdateViewTests(TestCase):
    """Ensure sample group metadata can be edited securely."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="editor",
            email="editor@example.com",
            password="pass12345",
        )
        self.other_user = User.objects.create_user(
            username="intruder",
            email="intruder@example.com",
            password="pass12345",
        )
        self.sample_group = SampleGroup.objects.create(
            name="Editable Group",
            created_by=self.user.organization_profile,
            contact_email="old@example.com",
        )

    def test_login_required(self) -> None:
        response = self.client.get(
            reverse("sample_group_edit", args=[self.sample_group.pk])
        )

        self.assertEqual(response.status_code, 302)
        self.assertIn(reverse("login"), response.url)

    def test_user_can_update_owned_group(self) -> None:
        self.client.force_login(self.user)

        response = self.client.post(
            reverse("sample_group_edit", args=[self.sample_group.pk]),
            {
                "name": "Updated Group",
                "contact_email": "new@example.com",
                "total_samples": 42,
            },
            follow=True,
        )

        self.assertRedirects(response, reverse("profile"))
        self.sample_group.refresh_from_db()
        self.assertEqual(self.sample_group.name, "Updated Group")
        self.assertEqual(self.sample_group.contact_email, "new@example.com")
        self.assertEqual(self.sample_group.total_samples, 42)

    def test_user_cannot_edit_other_organisations_group(self) -> None:
        other_group = SampleGroup.objects.create(
            name="Locked Group",
            created_by=self.other_user.organization_profile,
        )

        self.client.force_login(self.user)
        response = self.client.get(reverse("sample_group_edit", args=[other_group.pk]))

        self.assertEqual(response.status_code, 404)

    def test_update_view_get_queryset_uses_mixin_helper(self) -> None:
        request = RequestFactory().get("/sample-group/edit/")
        request.user = self.user
        view = SampleGroupUpdateView()
        view.setup(request, pk=self.sample_group.pk)

        queryset = SampleGroup.objects.filter(
            created_by=self.user.organization_profile
        )

        with patch.object(
            SampleGroupUpdateView, "get_owned_sample_groups", return_value=queryset
        ) as mock_groups:
            result = view.get_queryset()

        mock_groups.assert_called_once_with()
        self.assertIs(result, queryset)


class SampleGroupDeleteViewTests(TestCase):
    """Ensure the imported sample groups can be safely removed."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="deleter",
            email="deleter@example.com",
            password="pass12345",
        )
        self.other_user = User.objects.create_user(
            username="outsider",
            email="outsider@example.com",
            password="pass12345",
        )

        self.input_quality = InputQuality.objects.create(a260_a280=1.9)
        self.reference = ReferenceGenomeBuild.objects.create(
            build_name="GRCh38",
            build_version="v2",
        )
        self.sample_group = SampleGroup.objects.create(
            name="Deletable Group",
            created_by=self.user.organization_profile,
            contact_email="delete@example.com",
            input_quality=self.input_quality,
            reference_genome_build=self.reference,
        )

    def delete_url(self) -> str:
        return reverse("sample_group_delete", args=[self.sample_group.pk])

    def test_login_required(self) -> None:
        response = self.client.get(self.delete_url())
        login_url = reverse("login")
        self.assertEqual(response.status_code, 302)
        self.assertTrue(response.url.startswith(f"{login_url}?next="))

    def test_owner_sees_confirmation(self) -> None:
        self.client.force_login(self.user)
        response = self.client.get(self.delete_url())
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Delete Sample Group")
        self.assertContains(response, self.sample_group.name)

    def test_owner_can_delete_sample_group(self) -> None:
        self.client.force_login(self.user)
        response = self.client.post(self.delete_url(), follow=True)
        self.assertRedirects(response, reverse("profile"))
        self.assertFalse(SampleGroup.objects.filter(pk=self.sample_group.pk).exists())
        self.assertFalse(
            ReferenceGenomeBuild.objects.filter(pk=self.reference.pk).exists()
        )
        self.assertFalse(
            InputQuality.objects.filter(pk=self.input_quality.pk).exists()
        )

    def test_non_owner_cannot_delete_sample_group(self) -> None:
        self.client.force_login(self.other_user)
        response = self.client.post(self.delete_url())
        self.assertEqual(response.status_code, 404)
        self.assertTrue(SampleGroup.objects.filter(pk=self.sample_group.pk).exists())

    def test_shared_metadata_is_preserved(self) -> None:
        shared_reference = ReferenceGenomeBuild.objects.create(
            build_name="SharedRef",
            build_version="v1",
        )
        sibling = SampleGroup.objects.create(
            name="Sibling Group",
            created_by=self.user.organization_profile,
            reference_genome_build=shared_reference,
        )
        self.sample_group.reference_genome_build = shared_reference
        self.sample_group.save(update_fields=["reference_genome_build"])

        self.client.force_login(self.user)
        response = self.client.post(self.delete_url())
        self.assertRedirects(response, reverse("profile"))

        sibling.refresh_from_db()
        self.assertTrue(
            ReferenceGenomeBuild.objects.filter(pk=shared_reference.pk).exists()
        )

    def test_delete_view_get_queryset_uses_mixin_helper(self) -> None:
        request = RequestFactory().get(self.delete_url())
        request.user = self.user
        view = SampleGroupDeleteView()
        view.setup(request, pk=self.sample_group.pk)

        queryset = SampleGroup.objects.filter(
            created_by=self.user.organization_profile
        )

        with patch.object(
            SampleGroupDeleteView, "get_owned_sample_groups", return_value=queryset
        ) as mock_groups:
            result = view.get_queryset()

        mock_groups.assert_called_once_with()
        self.assertIs(result, queryset)
