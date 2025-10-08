"""Focused regression tests for the primary application views."""

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from ..forms import ImportDataForm, SearchForm
from ..models import AlleleFrequency, SampleGroup


class HomePageViewTests(TestCase):
    """Ensure the landing page renders the expected search form."""

    def test_search_form_present_in_context(self) -> None:
        response = self.client.get(reverse("home"))

        self.assertEqual(response.status_code, 200)
        self.assertIn("form", response.context)
        self.assertIsInstance(response.context["form"], SearchForm)


class ProfileViewTests(TestCase):
    """Exercise the profile dashboard and its supporting context."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="profile_user",
            password="password123",
            email="user@example.com",
        )
        self.other_user = User.objects.create_user(
            username="other_user",
            password="password123",
            email="other@example.com",
        )

    def test_profile_lists_user_groups_and_import_form(self) -> None:
        SampleGroup.objects.create(
            name="VCF Cohort",
            created_by=self.user.organization_profile,
        )

        self.client.login(username="profile_user", password="password123")
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["organization_profile"], self.user.organization_profile)
        self.assertIsInstance(response.context["import_form"], ImportDataForm)
        self.assertEqual(response.context["import_form_action"], reverse("import_data"))
        self.assertEqual(response.context["import_form_enctype"], "multipart/form-data")

    def test_profile_context_includes_only_owned_sample_groups(self) -> None:
        SampleGroup.objects.create(
            name="Alpha",
            created_by=self.user.organization_profile,
        )
        SampleGroup.objects.create(
            name="Beta",
            created_by=self.user.organization_profile,
        )
        SampleGroup.objects.create(
            name="Gamma",
            created_by=self.other_user.organization_profile,
        )

        self.client.force_login(self.user)
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)

        expected_queryset = SampleGroup.objects.filter(
            created_by=self.user.organization_profile
        ).order_by("name")
        self.assertQuerySetEqual(
            response.context["sample_groups"],
            SampleGroup.objects.filter(
                created_by=self.user.organization_profile
            ).order_by("name"),
            transform=lambda group: group,
        )
        self.assertIsInstance(response.context["import_form"], ImportDataForm)
        self.assertEqual(response.context["import_form_action"], reverse("import_data"))
        self.assertEqual(
            response.context["import_form_enctype"],
            "multipart/form-data",
        )


class ImportDataViewTests(TestCase):
    """Validate the VCF import workflow end to end."""

    VCF_CONTENT = """##fileformat=VCFv4.2
##SAMPLE=<ID=GroupA,Description=Imported group>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t1234\trsTest\tA\tT\t99\tPASS\tAF=0.5;CLNSIG=Pathogenic\tGT:GQ\t0/1:99
"""

    def setUp(self) -> None:
        super().setUp()
        self.user = get_user_model().objects.create_user(
            username="vcf_user",
            email="vcf@example.com",
            password="import-pass",
        )

    def test_import_creates_sample_group_and_variant(self) -> None:
        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "import.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        sample_group = SampleGroup.objects.get(name="GroupA")
        self.assertEqual(sample_group.allele_frequencies.count(), 1)
        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertEqual(allele.info.af, "0.5")
        self.assertEqual(allele.info.additional["clnsig"], "Pathogenic")
        self.assertEqual(allele.format.gt, "0/1")
        self.assertEqual(allele.format.gq, "99")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

        self.assertGreaterEqual(AlleleFrequency.objects.count(), 1)


class DeleteAccountViewTests(TestCase):
    """Exercise the account deletion endpoint."""

    def test_confirmed_post_deletes_user_and_logs_out(self) -> None:
        User = get_user_model()
        user = User.objects.create_user(
            username="delete_me",
            password="secure-pass",
            email="deleteme@example.com",
        )

        self.client.force_login(user)

        response = self.client.post(
            reverse("delete_account"),
            {"confirm": True},
            follow=True,
        )

        self.assertRedirects(response, reverse("home"))
        self.assertFalse(User.objects.filter(pk=user.pk).exists())
        self.assertNotIn("_auth_user_id", self.client.session)

        messages = list(get_messages(response.wsgi_request))
        self.assertIn(
            "Your account has been deleted.",
            [message.message for message in messages],
        )


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
