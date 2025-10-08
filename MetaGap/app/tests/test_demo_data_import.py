"""Tests validating that the demo data fixture imports as expected."""

from django.apps import apps
from django.core.management import call_command
from django.test import TestCase


class DemoDataImportTests(TestCase):
    """Ensure the demo fixture loads cleanly and links related models."""

    def setUp(self):
        super().setUp()
        # Load the fixture explicitly to mimic the documented workflow.
        call_command("loaddata", "demo_data.json", verbosity=0)

    def test_users_and_profiles_created(self):
        User = apps.get_model("auth", "User")
        OrganizationProfile = apps.get_model("app", "OrganizationProfile")

        self.assertEqual(User.objects.count(), 3)
        self.assertEqual(OrganizationProfile.objects.count(), 3)

        partner_profile = OrganizationProfile.objects.get(user__username="demo_partner")
        self.assertEqual(partner_profile.organization_name, "Microbiome Discovery Lab")

    def test_sample_groups_cross_linked(self):
        SampleGroup = apps.get_model("app", "SampleGroup")

        self.assertEqual(SampleGroup.objects.count(), 3)

        microbiome_group = SampleGroup.objects.get(name="Microbiome Case Study")
        self.assertEqual(microbiome_group.reference_genome_build.build_name, "T2T-CHM13")
        self.assertEqual(
            microbiome_group.allele_frequencies.get().variant_id,
            "bft_v1",
        )
        self.assertEqual(microbiome_group.bioinfo_alignment.software, "Bowtie2")

    def test_variant_metadata_round_trip(self):
        AlleleFrequency = apps.get_model("app", "AlleleFrequency")

        allele = AlleleFrequency.objects.get(variant_id="rs121913529")
        self.assertEqual(allele.info.additional["clinvar_significance"], "Pathogenic")

        microbiome_allele = AlleleFrequency.objects.get(variant_id="bft_v1")
        self.assertEqual(
            microbiome_allele.format.additional["sample_id"],
            "Microbiome_Case01",
        )
