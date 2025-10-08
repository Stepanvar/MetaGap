"""Tests for forms declared in :mod:`app.forms`."""

from django.core.files.uploadedfile import SimpleUploadedFile
from django.contrib.auth.models import User
from django.test import SimpleTestCase, TestCase

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


class SampleGroupFormTests(TestCase):
    """Tests for :class:`app.forms.SampleGroupForm`."""

    def test_save_with_commit_true_assigns_organization_profile(self):
        user = User.objects.create_user("creator", "creator@example.com", "secret")
        form = forms.SampleGroupForm(data={"name": "Test Group"}, user=user)

        self.assertTrue(form.is_valid())

        sample_group = form.save()

        self.assertIsNotNone(sample_group.pk)
        self.assertEqual(sample_group.created_by, user.organization_profile)

    def test_save_with_commit_false_assigns_organization_profile(self):
        user = User.objects.create_user("creator", "creator@example.com", "secret")
        form = forms.SampleGroupForm(data={"name": "Deferred Group"}, user=user)

        self.assertTrue(form.is_valid())

        sample_group = form.save(commit=False)

        self.assertIsNone(sample_group.pk)
        self.assertEqual(sample_group.created_by, user.organization_profile)

    def test_missing_user_raises_value_error(self):
        form = forms.SampleGroupForm(data={"name": "Nameless"})

        self.assertTrue(form.is_valid())

        with self.assertRaisesMessage(
            ValueError, "SampleGroupForm.save() requires a user when creating a sample group."
        ):
            form.save()
