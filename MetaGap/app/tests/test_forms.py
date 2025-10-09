"""Tests for forms declared in :mod:`app.forms`."""

from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import SimpleTestCase, TestCase, override_settings
from app import forms
from app.models import (
    BioinfoAlignment,
    MaterialType,
    OrganizationProfile,
    SampleGroup,
    SampleOrigin,
)


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

    def test_sample_origin_fields_create_related_instance(self):
        user = User.objects.create_user("creator", "creator@example.com", "secret")
        form = forms.SampleGroupForm(
            data={
                "name": "Origin Group",
                "sample_origin_tissue": "Brain",
                "sample_origin_collection_method": "Stereotactic biopsy",
                "sample_origin_storage_conditions": "-80C",
                "sample_origin_time_stored": "3 months",
            },
            user=user,
        )

        self.assertTrue(form.is_valid())

        sample_group = form.save()

        self.assertIsNotNone(sample_group.sample_origin)
        self.assertEqual(sample_group.sample_origin.tissue, "Brain")
        self.assertEqual(
            sample_group.sample_origin.collection_method, "Stereotactic biopsy"
        )
        self.assertEqual(sample_group.sample_origin.storage_conditions, "-80C")
        self.assertEqual(sample_group.sample_origin.time_stored, "3 months")

    def test_sample_origin_fields_update_existing_instance(self):
        user = User.objects.create_user("creator", "creator@example.com", "secret")
        origin = SampleOrigin.objects.create(
            tissue="Kidney",
            collection_method="Biopsy",
            storage_conditions="Chilled",
            time_stored="1 week",
        )
        sample_group = SampleGroup.objects.create(
            name="Kidney Group",
            created_by=user.organization_profile,
            sample_origin=origin,
        )

        form = forms.SampleGroupForm(
            data={
                "name": "Kidney Group",
                "sample_origin_tissue": "Liver",
                "sample_origin_collection_method": "Surgical",
                "sample_origin_storage_conditions": "Ambient",
                "sample_origin_time_stored": "",
            },
            instance=sample_group,
            user=user,
        )

        self.assertTrue(form.is_valid())

        updated_group = form.save()

        self.assertEqual(updated_group.pk, sample_group.pk)
        self.assertIsNotNone(updated_group.sample_origin)
        updated_group.sample_origin.refresh_from_db()
        self.assertEqual(updated_group.sample_origin.tissue, "Liver")
        self.assertEqual(updated_group.sample_origin.collection_method, "Surgical")
        self.assertEqual(updated_group.sample_origin.storage_conditions, "Ambient")
        self.assertIsNone(updated_group.sample_origin.time_stored)

    def test_creatable_metadata_fields_allow_new_entries(self):
        user = User.objects.create_user("creator", "creator@example.com", "secret")
        form = forms.SampleGroupForm(
            data={
                "name": "Metadata Rich Group",
                "material_type": "DNA",
                "bioinfo_alignment": "BWA-MEM",
            },
            user=user,
        )

        self.assertTrue(form.is_valid())

        sample_group = form.save()

        self.assertEqual(MaterialType.objects.count(), 1)
        self.assertIsNotNone(sample_group.material_type)
        self.assertEqual(sample_group.material_type.material_type, "DNA")

        self.assertEqual(BioinfoAlignment.objects.count(), 1)
        self.assertIsNotNone(sample_group.bioinfo_alignment)
        self.assertEqual(sample_group.bioinfo_alignment.tool, "BWA-MEM")
