"""Tests for forms declared in :mod:`app.forms`."""

from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import SimpleTestCase

from app import forms


class ImportDataFormTests(SimpleTestCase):
    """Validate the behaviour of :class:`app.forms.ImportDataForm`."""

    def test_valid_file_populates_cleaned_data(self):
        form = forms.ImportDataForm(
            files={"data_file": SimpleUploadedFile("variants.vcf", b"##fileformat=VCF\n")}
        )

        self.assertTrue(form.is_valid())
        self.assertIn("data_file", form.cleaned_data)

    def test_missing_file_fails_validation(self):
        form = forms.ImportDataForm(data={})

        self.assertFalse(form.is_valid())
        self.assertIn("data_file", form.errors)
