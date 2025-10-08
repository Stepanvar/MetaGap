"""Tests for forms declared in :mod:`app.forms`."""

from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import SimpleTestCase, TestCase, override_settings
from app import forms
from app.models import OrganizationProfile


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


@override_settings(MIGRATION_MODULES={"app": None})
class CustomUserCreationFormTests(TestCase):
    """Behaviour specific to :class:`app.forms.CustomUserCreationForm`."""

    def test_commit_false_defers_profile_creation_until_save_m2m(self):
        form = forms.CustomUserCreationForm(
            data={
                "username": "new-user",
                "email": "user@example.com",
                "password1": "complex-pass-123",
                "password2": "complex-pass-123",
                "organization_name": "New Org",
            }
        )

        self.assertTrue(form.is_valid())

        user = form.save(commit=False)

        # Saving with commit=False should not touch the database yet.
        self.assertFalse(User.objects.filter(username="new-user").exists())
        self.assertFalse(OrganizationProfile.objects.exists())

        # Once the caller saves the user and invokes ``save_m2m`` the profile appears.
        user.save()
        form.save_m2m()

        profile = OrganizationProfile.objects.get(user=user)
        self.assertEqual(profile.organization_name, "New Org")
        self.assertFalse(hasattr(form, "_pending_organization_name"))


@override_settings(MIGRATION_MODULES={"app": None})
class EditProfileFormTests(TestCase):
    """Behaviour specific to :class:`app.forms.EditProfileForm`."""

    def setUp(self):
        self.user = User.objects.create_user(
            username="existing", email="before@example.com", password="secret-pass"
        )
        self.profile = self.user.organization_profile
        self.profile.organization_name = "Original Org"
        self.profile.save()

    def test_commit_false_defers_profile_update_until_save_m2m(self):
        form = forms.EditProfileForm(
            data={
                "username": "existing",
                "email": "after@example.com",
                "first_name": "Updated",
                "last_name": "User",
                "organization_name": "Updated Org",
            },
            instance=self.user,
        )

        self.assertTrue(form.is_valid())

        user = form.save(commit=False)

        # The organization profile should still reflect the previous value at this point.
        self.profile.refresh_from_db()
        self.assertEqual(self.profile.organization_name, "Original Org")

        user.save()
        form.save_m2m()

        self.profile.refresh_from_db()
        self.assertEqual(self.profile.organization_name, "Updated Org")
        self.assertFalse(hasattr(form, "_pending_organization_name"))
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
