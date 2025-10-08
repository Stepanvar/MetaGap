"""Integration tests covering the VCF import view."""

from __future__ import annotations

from unittest import mock

from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from ..models import AlleleFrequency, SampleGroup


class ImportDataViewIntegrationTests(TestCase):
    """Ensure VCF uploads produce related allele frequency records."""

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

    def test_import_links_alleles_to_sample_group(self) -> None:
        """Posting a small VCF links the created alleles to the sample group."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "import.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))

        sample_group = SampleGroup.objects.get(name="GroupA")
        alleles = list(sample_group.allele_frequencies.order_by("id"))

        self.assertEqual(len(alleles), 1)
        allele = alleles[0]

        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertEqual(allele.info.af, "0.5")
        self.assertEqual(allele.info.additional["clnsig"], "Pathogenic")
        self.assertEqual(allele.format.gt, "0/1")
        self.assertEqual(allele.format.gq, "99")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

        self.assertEqual(
            list(sample_group.allele_frequencies.all()),
            alleles,
        )
        self.assertGreaterEqual(AlleleFrequency.objects.count(), 1)

    @mock.patch("app.views.pysam.VariantFile", side_effect=OSError("boom"))
    def test_import_falls_back_to_text_parser(self, mocked_variant_file: mock.Mock) -> None:
        """If pysam fails to open the VCF, the text fallback still imports data."""

        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "import.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        mocked_variant_file.assert_called()

        sample_group = SampleGroup.objects.get(name="GroupA")
        self.assertEqual(sample_group.comments, "Imported group")

        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertEqual(allele.info.af, "0.5")
        self.assertEqual(allele.format.gt, "0/1")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")
