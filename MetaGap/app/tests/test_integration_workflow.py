"""High level integration scenarios spanning import, search, and export."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from ..models import AlleleFrequency, SampleGroup
from ..services.vcf_importer import VCFImporter


class UserWorkflowIntegrationTests(TestCase):
    """Simulate a representative MetaGap user workflow."""

    WORKFLOW_VCF = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##SAMPLE=<ID=WorkflowGroup,Description=Integration test cohort>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t555\trsWorkflow\tA\tG\t99\tPASS\tAF=0.42\tGT:GQ\t0/1:85
"""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="workflow_user",
            email="workflow@example.com",
            password="workflow-pass",
        )

    def _import_vcf(self) -> SampleGroup:
        importer = VCFImporter(self.user)
        with tempfile.TemporaryDirectory(prefix="workflow-vcf") as temp_dir:
            path = Path(temp_dir) / "workflow_import.vcf"
            path.write_text(self.WORKFLOW_VCF)
            sample_group = importer.import_file(str(path))
        # importer.import_file deletes intermediate data only when errors occur;
        # ensure we return the created group for assertions below.
        return sample_group

    def test_user_can_import_search_and_export_variants(self) -> None:
        sample_group = self._import_vcf()
        allele = sample_group.allele_frequencies.get()

        self.client.force_login(self.user)

        profile_response = self.client.get(reverse("profile"))
        self.assertEqual(profile_response.status_code, 200)
        self.assertContains(profile_response, sample_group.name)

        search_response = self.client.get(reverse("search_results"), {"query": "rsWorkflow"})
        self.assertEqual(search_response.status_code, 200)
        table = search_response.context["table"]
        records = [row.record for row in table.rows]
        self.assertEqual(records, [allele])

        detail_response = self.client.get(
            reverse("sample_group_detail", args=[sample_group.pk])
        )
        self.assertEqual(detail_response.status_code, 200)
        self.assertEqual(
            detail_response.context.get("sample_group")
            or detail_response.context.get("object"),
            sample_group,
        )

        export_response = self.client.get(
            reverse("sample_group_export", args=[sample_group.pk])
        )
        self.assertEqual(export_response.status_code, 200)
        self.assertEqual(export_response["Content-Type"], "text/csv")

        reader = csv.DictReader(export_response.content.decode().splitlines())
        rows = list(reader)
        self.assertEqual(len(rows), 1)
        exported_row = rows[0]
        self.assertEqual(exported_row["chrom"], allele.chrom)
        self.assertEqual(exported_row["pos"], str(allele.pos))
        self.assertEqual(exported_row["ref"], allele.ref)
        self.assertEqual(exported_row["alt"], allele.alt)
        self.assertEqual(exported_row["variant_id"], allele.variant_id)

        self.assertTrue(
            SampleGroup.objects.filter(pk=sample_group.pk).exists(),
            "Sample group should persist after export",
        )
        self.assertTrue(
            AlleleFrequency.objects.filter(pk=allele.pk).exists(),
            "Variant should persist after integration workflow",
        )
