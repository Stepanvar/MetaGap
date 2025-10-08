"""Focused regression tests for the primary application views."""

from app.models import SampleGroup

from __future__ import annotations

from django.contrib.auth import get_user_model
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
        self.user = get_user_model().objects.create_user(
            username="importer",
            email="importer@example.com",
            password="super-secret",
        )

    def test_profile_lists_user_groups_and_import_form(self) -> None:
        SampleGroup.objects.create(
            name="VCF Cohort",
            created_by=self.user.organization_profile,
        )

        self.client.login(username="importer", password="super-secret")
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)


class ProfileViewTests(TestCase):
    def setUp(self):
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="profile_user", password="password123", email="user@example.com"
        )
        self.other_user = User.objects.create_user(
            username="other_user", password="password123", email="other@example.com"
        )

    def test_profile_context_includes_owned_sample_groups(self):
        SampleGroup.objects.create(name="Alpha", created_by=self.user)
        SampleGroup.objects.create(name="Beta", created_by=self.user)
        SampleGroup.objects.create(name="Gamma", created_by=self.other_user)

        self.client.force_login(self.user)
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)

        sample_groups = list(response.context["sample_groups"])
        self.assertEqual([group.name for group in sample_groups], ["Alpha", "Beta"])

        self.assertEqual(
            response.context["organization_profile"],
            self.user.organization_profile,
        )
        self.assertQuerySetEqual(
            response.context["sample_groups"],
            SampleGroup.objects.filter(created_by=self.user.organization_profile).order_by("name"),
            transform=lambda group: group,
        )
        self.assertIsInstance(response.context["import_form"], ImportDataForm)
        self.assertEqual(response.context["import_form_action"], reverse("import_data"))
        self.assertEqual(response.context["import_form_enctype"], "multipart/form-data")


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
