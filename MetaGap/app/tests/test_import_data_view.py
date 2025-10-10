"""Tests for the VCF import form view."""

from __future__ import annotations

from unittest import mock

from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse


class ImportDataViewTests(TestCase):
    """Ensure the import view delegates to the VCF importer service."""

    def setUp(self) -> None:
        super().setUp()
        self.user = get_user_model().objects.create_user(
            username="vcf_user",
            email="vcf@example.com",
            password="import-pass",
        )
        self.client.login(username="vcf_user", password="import-pass")

    def test_success_flow_uses_service_and_displays_messages(self) -> None:
        """Successful uploads call the importer and surface warnings."""

        importer_instance = mock.Mock()
        importer_instance.import_file.return_value = mock.Mock(name="SampleGroupMock")
        importer_instance.import_file.return_value.name = "Imported Group"
        importer_instance.warnings = ["Fallback message"]

        with mock.patch("app.views.VCFImporter", return_value=importer_instance) as importer_cls:
            upload = SimpleUploadedFile(
                "import.vcf",
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                content_type="text/vcf",
            )
            response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        importer_cls.assert_called_once_with(self.user)
        importer_instance.import_file.assert_called_once()

        messages = [message.message for message in response.wsgi_request._messages]
        self.assertIn("Imported Imported Group successfully.", messages)
        self.assertIn("Fallback message", messages)

    def test_error_flow_surfaces_exception_message(self) -> None:
        """Errors from the importer are shown to the user."""

        importer_instance = mock.Mock()
        importer_instance.import_file.side_effect = ValueError("boom")
        importer_instance.warnings = []

        with mock.patch("app.views.VCFImporter", return_value=importer_instance):
            upload = SimpleUploadedFile(
                "import.vcf",
                b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                content_type="text/vcf",
            )
            response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        importer_instance.import_file.assert_called_once()

        messages = [message.message for message in response.wsgi_request._messages]
        self.assertIn("An error occurred: boom", messages)
