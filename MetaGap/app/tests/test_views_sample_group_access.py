"""Tests covering access control and mutation behaviour for sample group views."""

from __future__ import annotations

import os

from unittest.mock import patch

from django.conf import settings
from django.contrib.messages import get_messages
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import RequestFactory, TestCase
from django.urls import reverse

from ..forms import ImportDataForm, SampleGroupForm
from ..models import SampleGroup
from ..views import SampleGroupDeleteView
from .sample_group_test_mixins import SampleGroupViewMatrixMixin


class SampleGroupViewAccessMatrixTests(SampleGroupViewMatrixMixin, TestCase):
    """Matrix covering anonymous and non-owner access to protected views."""

    def test_anonymous_users_are_redirected(self) -> None:
        login_url = reverse("login")
        for url in (
            self.get_detail_url(),
            self.get_update_url(),
            self.get_delete_url(),
            self.get_import_url(),
        ):
            response = self.client.get(url)
            self.assertEqual(response.status_code, 302)
            self.assertTrue(
                response.url.startswith(login_url),
                msg=f"Expected redirect to login for {url}, got {response.url}",
            )

    def test_non_owners_receive_not_found_for_object_views(self) -> None:
        self.client.force_login(self.intruder)
        for url in (
            self.get_detail_url(),
            self.get_update_url(),
            self.get_delete_url(),
        ):
            response = self.client.get(url)
            self.assertIn(response.status_code, {403, 404})

    def test_owner_detail_view_renders_metadata(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(self.get_detail_url())

        self.assertEqual(response.status_code, 200)
        context_object = response.context.get("sample_group") or response.context.get(
            "object"
        )
        self.assertEqual(context_object, self.owned_group)
        self.assertIn(self.owned_allele.variant_id, response.content.decode())

    def test_import_view_uses_correct_form_setup(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(self.get_import_url())

        self.assertEqual(response.status_code, 200)
        self.assertIsInstance(response.context["form"], ImportDataForm)
        self.assertEqual(response.context.get("form_action"), reverse("import_data"))
        self.assertEqual(
            response.context.get("form_enctype"),
            "multipart/form-data",
        )


class SampleGroupMutationSecurityTests(SampleGroupViewMatrixMixin, TestCase):
    """Exercise mutation endpoints with ownership and form validation checks."""

    def test_owner_update_get_includes_csrf_and_bound_form(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(self.get_update_url())

        self.assertEqual(response.status_code, 200)
        self.assertIn("csrfmiddlewaretoken", response.content.decode())
        self.assertIsInstance(response.context["form"], SampleGroupForm)
        self.assertEqual(
            response.context["form"].instance.pk,
            self.owned_group.pk,
        )

    def test_owner_successful_update_post_sets_message(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.post(
            self.get_update_url(),
            {"name": "Matrix Updated Cohort"},
            follow=True,
        )

        self.assertRedirects(response, reverse("profile"))
        self.owned_group.refresh_from_db()
        self.assertEqual(self.owned_group.name, "Matrix Updated Cohort")

        messages = [message.message for message in get_messages(response.wsgi_request)]
        self.assertIn("Sample group metadata updated successfully.", messages)

    def test_owner_invalid_update_post_shows_errors(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.post(self.get_update_url(), {"name": ""})

        self.assertEqual(response.status_code, 200)
        form = response.context["form"]
        self.assertTrue(form.errors)
        self.assertIn("name", form.errors)
        messages = [message.message for message in get_messages(response.wsgi_request)]
        self.assertEqual(messages, [])

    def test_non_owner_update_post_is_blocked(self) -> None:
        original_name = self.owned_group.name
        self.client.force_login(self.intruder)
        response = self.client.post(
            self.get_update_url(),
            {"name": "Hacked Cohort"},
        )

        self.assertEqual(response.status_code, 404)
        self.owned_group.refresh_from_db()
        self.assertEqual(self.owned_group.name, original_name)

    def test_owner_delete_get_includes_csrf(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(self.get_delete_url())

        self.assertEqual(response.status_code, 200)
        self.assertIn("csrfmiddlewaretoken", response.content.decode())
        self.assertEqual(response.context["cancel_url"], self.get_detail_url())

    def test_owner_delete_post_cleans_up_and_redirects(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.post(self.get_delete_url(), follow=True)

        self.assertRedirects(response, reverse("profile"))
        self.assertFalse(SampleGroup.objects.filter(pk=self.owned_group.pk).exists())

    def test_delete_view_emits_success_message(self) -> None:
        request = RequestFactory().post(self.get_delete_url())
        request.user = self.owner
        view = SampleGroupDeleteView()
        view.setup(request, pk=self.owned_group.pk)

        with patch("app.views.messages.success") as mock_success:
            response = view.delete(request, pk=self.owned_group.pk)

        self.assertEqual(response.status_code, 302)
        mock_success.assert_called_once_with(
            request, "Deleted Matrix Owner Cohort successfully."
        )

    def test_non_owner_delete_post_is_blocked(self) -> None:
        self.client.force_login(self.intruder)
        response = self.client.post(self.get_delete_url())

        self.assertEqual(response.status_code, 404)
        self.assertTrue(SampleGroup.objects.filter(pk=self.owned_group.pk).exists())

    def test_import_get_filters_sample_groups_by_owner(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.get(self.get_import_url())

        self.assertEqual(response.status_code, 200)
        sample_groups = response.context["sample_groups"].order_by("name")
        self.assertQuerySetEqual(
            sample_groups,
            SampleGroup.objects.filter(
                created_by=self.owner.organization_profile
            ).order_by("name"),
            transform=lambda group: group,
        )
        self.assertNotIn(self.intruder_group, list(sample_groups))
        self.assertIn("csrfmiddlewaretoken", response.content.decode())

    def test_import_get_excludes_foreign_groups_for_intruder(self) -> None:
        self.client.force_login(self.intruder)
        response = self.client.get(self.get_import_url())

        self.assertEqual(response.status_code, 200)
        sample_groups = list(response.context["sample_groups"].order_by("name"))
        self.assertEqual(sample_groups, [self.intruder_group])
        self.assertNotIn(self.owned_group, sample_groups)

    def test_import_post_requires_login(self) -> None:
        response = self.client.post(self.get_import_url())

        login_url = reverse("login")
        self.assertEqual(response.status_code, 302)
        self.assertTrue(response.url.startswith(login_url))

    def test_import_post_without_file_shows_form_errors(self) -> None:
        self.client.force_login(self.owner)
        response = self.client.post(self.get_import_url(), {})

        self.assertEqual(response.status_code, 200)
        form = response.context["form"]
        self.assertIn("data_file", form.errors)

    def test_import_post_with_file_triggers_importer(self) -> None:
        self.client.force_login(self.owner)
        upload = SimpleUploadedFile("matrix.vcf", b"Fake VCF content", content_type="text/vcf")

        temp_path = "tmp/matrix.vcf"
        expected_full_path = os.path.join(settings.MEDIA_ROOT, temp_path)

        with patch("app.views.default_storage") as mock_storage, patch(
            "app.views.VCFImporter"
        ) as mock_importer_cls:
            mock_storage.save.return_value = temp_path
            mock_importer = mock_importer_cls.return_value
            mock_importer.import_file.return_value = self.owned_group
            mock_importer.warnings = ["Matrix warning"]

            response = self.client.post(
                self.get_import_url(),
                {"data_file": upload},
                follow=True,
            )

        self.assertRedirects(response, reverse("profile"))
        mock_storage.save.assert_called_once()
        mock_importer.import_file.assert_called_once_with(expected_full_path)
        mock_storage.delete.assert_called_once_with(temp_path)

        messages = [message.message for message in get_messages(response.wsgi_request)]
        self.assertIn("Imported Matrix Owner Cohort successfully.", messages)
        self.assertIn("Matrix warning", messages)
