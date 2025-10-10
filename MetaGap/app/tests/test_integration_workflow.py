"""High level integration scenarios spanning import, search, and export."""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path
from unittest import mock

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.contrib.auth.forms import PasswordResetForm
from django.core import mail
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from ..forms import CustomUserCreationForm
from ..models import AlleleFrequency, SampleGroup
from ..services.vcf_importer import VCFImporter
from .test_vcf_importer import VCFImporterTests as VCFImporterFixtures


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

    SECONDARY_WORKFLOW_VCF = """##fileformat=VCFv4.2
##contig=<ID=2>
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">
##SAMPLE=<ID=WorkflowGroupB,Description=Second integration cohort>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample002
2\t777\trsWorkflowTwo\tC\tT\t85\tPASS\tAF=0.58\tGT:GQ\t0/1:92
"""

    def setUp(self) -> None:
        super().setUp()
        self.user = self._create_user(
            username="workflow_user",
            email="workflow@example.com",
        )
        mail.outbox.clear()

    def _create_user(
        self, *, username: str, email: str, password: str = "workflow-pass"
    ):
        User = get_user_model()
        return User.objects.create_user(
            username=username,
            email=email,
            password=password,
        )

    def _import_vcf_content(self, content: str, filename: str) -> SampleGroup:
        importer = VCFImporter(self.user)
        with tempfile.TemporaryDirectory(prefix="workflow-vcf") as temp_dir:
            path = Path(temp_dir) / filename
            path.write_text(content)
            sample_group = importer.import_file(str(path))
        # importer.import_file deletes intermediate data only when errors occur;
        # ensure we return the created group for assertions below.
        return sample_group

    def _import_vcf(self) -> SampleGroup:
        return self._import_vcf_content(self.WORKFLOW_VCF, "workflow_import.vcf")

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

    def test_import_with_warning_feedback_during_onboarding(self) -> None:
        registration_form = CustomUserCreationForm(
            data={
                "username": "onboard_user",
                "email": "onboard@example.com",
                "password1": "ComplexPass123",
                "password2": "ComplexPass123",
                "organization_name": "Onboard Labs",
            }
        )
        self.assertTrue(registration_form.is_valid(), registration_form.errors)
        onboard_user = registration_form.save()
        self.addCleanup(onboard_user.delete)

        reset_form = PasswordResetForm(data={"email": onboard_user.email})
        self.assertTrue(reset_form.is_valid(), reset_form.errors)
        with self.settings(ROOT_URLCONF="app.tests.urls"):
            reset_form.save(domain_override="testserver")
        self.assertEqual(len(mail.outbox), 1)

        self.client.force_login(onboard_user)
        uploaded_file = SimpleUploadedFile(
            "malformed_header.vcf",
            VCFImporterFixtures.MALFORMED_HEADER_VCF.encode(),
            content_type="text/vcf",
        )

        with mock.patch(
            "app.services.vcf_importer.pysam.VariantFile",
            side_effect=ValueError("malformed header"),
        ):
            response = self.client.post(
                reverse("import_data"),
                {"data_file": uploaded_file},
                follow=True,
            )

        self.assertEqual(response.status_code, 200)
        messages = list(response.context["messages"])
        self.assertTrue(
            any("Imported" in message.message for message in messages),
            "Successful import should trigger a confirmation message.",
        )
        self.assertTrue(
            any("Falling back to a text parser" in message.message for message in messages),
            "Importer warnings should surface to the user during onboarding.",
        )

        sample_group = SampleGroup.objects.get(name="malformed_header")
        self.assertEqual(sample_group.created_by, onboard_user.organization_profile)
        self.assertEqual(sample_group.allele_frequencies.count(), 1)

    def test_search_filters_across_multiple_sample_groups(self) -> None:
        first_group = self._import_vcf()
        first_group.name = "Cohort Alpha"
        first_group.source_lab = "Alpha Lab"
        first_group.save()

        second_group = self._import_vcf_content(
            self.SECONDARY_WORKFLOW_VCF,
            "workflow_import_second.vcf",
        )
        second_group.name = "Cohort Beta"
        second_group.source_lab = "Beta Lab"
        second_group.save()

        allele_one = first_group.allele_frequencies.get()
        allele_two = second_group.allele_frequencies.get()

        response_all = self.client.get(reverse("search_results"), {"query": "rsWorkflow"})
        self.assertEqual(response_all.status_code, 200)
        records_all = [row.record for row in response_all.context["table"].rows]
        self.assertCountEqual(records_all, [allele_one, allele_two])

        response_filtered = self.client.get(
            reverse("search_results"),
            {"sample_group_source_lab": "Beta Lab"},
        )
        self.assertEqual(response_filtered.status_code, 200)
        filtered_records = [
            row.record for row in response_filtered.context["table"].rows
        ]
        self.assertEqual(filtered_records, [allele_two])

    def test_export_respects_query_parameters(self) -> None:
        sample_group = self._import_vcf()
        allele = sample_group.allele_frequencies.get()

        self.client.force_login(self.user)
        export_response = self.client.get(
            reverse("sample_group_export", args=[sample_group.pk]),
            {"format": "tsv"},
        )

        self.assertEqual(export_response.status_code, 200)
        self.assertEqual(
            export_response["Content-Type"],
            "text/tab-separated-values",
        )

        lines = export_response.content.decode().splitlines()
        self.assertGreaterEqual(len(lines), 2)
        _header, row = lines[0], lines[1]
        self.assertIn("\t", row)
        values = row.split("\t")
        self.assertIn(str(allele.pos), values)
        self.assertIn(allele.variant_id, values)

    def test_sample_group_delete_respects_ownership(self) -> None:
        sample_group = SampleGroup.objects.create(
            name="Workflow Delete Cohort",
            created_by=self.user.organization_profile,
            contact_email="owner@example.com",
        )
        intruder = self._create_user(
            username="workflow_intruder",
            email="intruder@example.com",
        )
        delete_url = reverse("sample_group_delete", args=[sample_group.pk])

        with mock.patch("app.views.messages.success") as mock_success:
            self.client.force_login(intruder)
            intruder_response = self.client.post(delete_url)
            self.assertEqual(intruder_response.status_code, 404)
            self.assertTrue(SampleGroup.objects.filter(pk=sample_group.pk).exists())
            intruder_messages = list(get_messages(intruder_response.wsgi_request))
            self.assertFalse(intruder_messages)
            mock_success.assert_not_called()

            self.client.logout()
            self.client.force_login(self.user)
            owner_response = self.client.delete(delete_url)

        self.assertEqual(owner_response.status_code, 302)
        self.assertEqual(owner_response["Location"], reverse("profile"))
        self.assertFalse(SampleGroup.objects.filter(pk=sample_group.pk).exists())

        mock_success.assert_called_once_with(
            mock.ANY,
            f"Deleted {sample_group.name} successfully.",
        )

        follow_response = self.client.get(reverse("profile"))
        self.assertEqual(follow_response.status_code, 200)
