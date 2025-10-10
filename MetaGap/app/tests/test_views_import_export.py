"""Tests for sample group import and export view workflows."""

from __future__ import annotations

import csv
import io

from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import NoReverseMatch, reverse

from ..models import (
    AlleleFrequency,
    Format,
    Info,
    SampleGroup,
)
from ..views import ImportDataView
from .sample_group_test_mixins import SampleGroupTestDataMixin


class SampleGroupExportViewTests(SampleGroupTestDataMixin, TestCase):
    """Validate CSV exports for sample groups and their access controls."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.owner = User.objects.create_user(
            username="export_owner",
            password="export-pass",
            email="export-owner@example.com",
        )
        self.other_user = User.objects.create_user(
            username="export_other",
            password="export-pass",
            email="export-other@example.com",
        )
        self.sample_group, self.allele = self.create_sample_group_with_variant(self.owner)

    def export_url(self) -> str:
        for name in (
            "sample_group_export",
            "profile_sample_group_export",
            "sample-group-export",
        ):
            try:
                return reverse(name, args=[self.sample_group.pk])
            except NoReverseMatch:
                continue
        self.fail("Sample group export route is not configured")

    def test_owner_receives_csv_with_variant_rows(self) -> None:
        self.client.force_login(self.owner)

        response = self.client.get(self.export_url())

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "text/csv")

        content_stream = io.StringIO(response.content.decode())
        reader = csv.DictReader(content_stream)
        self.assertIsNotNone(reader.fieldnames)

        normalized_headers = {header.lower() for header in reader.fieldnames}
        for expected_header in {
            "chrom",
            "pos",
            "ref",
            "alt",
            "variant_id",
            "info_af",
            "info_clinvar_significance",
        }:
            self.assertIn(expected_header, normalized_headers)

        rows = list(reader)
        self.assertEqual(len(rows), 1)

        row = rows[0]
        self.assertEqual(row.get("chrom") or row.get("Chrom"), self.allele.chrom)
        self.assertEqual(row.get("pos") or row.get("Pos"), str(self.allele.pos))
        self.assertEqual(row.get("ref") or row.get("Ref"), self.allele.ref)
        self.assertEqual(row.get("alt") or row.get("Alt"), self.allele.alt)
        self.assertEqual(
            row.get("variant_id") or row.get("Variant_ID"), self.allele.variant_id
        )

        self.assertEqual(row.get("info_af"), str(self.allele.info.af))
        self.assertEqual(row.get("info_clinvar_significance"), "Pathogenic")

    def test_non_owner_cannot_export(self) -> None:
        self.client.force_login(self.other_user)

        response = self.client.get(self.export_url())

        self.assertIn(response.status_code, {403, 404})


class ImportDataViewTests(TestCase):
    """Validate the VCF import workflow end to end."""

    VCF_CONTENT = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
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
        self.assertAlmostEqual(float(allele.info.af), 0.5)
        self.assertEqual(allele.info.additional["clnsig"], "Pathogenic")
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.fields["gq"], "99")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

        self.assertGreaterEqual(AlleleFrequency.objects.count(), 1)


class ImportDataPageTests(TestCase):
    """Verify the import page lists previously uploaded sample groups."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="importviewer",
            email="importviewer@example.com",
            password="pass12345",
        )
        self.sample_group = SampleGroup.objects.create(
            name="Listed Group",
            created_by=self.user.organization_profile,
        )

    def test_import_page_displays_existing_groups(self) -> None:
        self.client.force_login(self.user)
        response = self.client.get(reverse("import_data"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, self.sample_group.name)
        self.assertContains(
            response,
            reverse("sample_group_detail", args=[self.sample_group.pk]),
        )
        self.assertContains(
            response,
            reverse("sample_group_delete", args=[self.sample_group.pk]),
        )

    def test_import_page_uses_mixin_helper(self) -> None:
        self.client.force_login(self.user)
        queryset = SampleGroup.objects.filter(
            created_by=self.user.organization_profile
        )

        with patch.object(
            ImportDataView, "get_owned_sample_groups", return_value=queryset
        ) as mock_groups:
            response = self.client.get(reverse("import_data"))

        mock_groups.assert_called_once_with()
        self.assertEqual(
            list(response.context["sample_groups"]),
            list(queryset.order_by("name")),
        )
