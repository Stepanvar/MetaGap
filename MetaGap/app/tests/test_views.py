"""Focused regression tests for the primary application views."""

from __future__ import annotations

import csv
import io

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import RequestFactory, TestCase
from django.urls import NoReverseMatch, reverse

from ..filters import SampleGroupFilter
from ..forms import ImportDataForm, SearchForm
from ..models import (
    AlleleFrequency,
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    Format,
    GenomeComplexity,
    IlluminaSeq,
    Info,
    InputQuality,
    IonTorrentSeq,
    LibraryConstruction,
    MaterialType,
    OntSeq,
    PacBioSeq,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)
from ..views import EditProfileView


class SampleGroupTestDataMixin:
    """Utility helpers for constructing richly populated sample groups."""

    def create_sample_group_with_variant(self, owner):
        reference = ReferenceGenomeBuild.objects.create(
            build_name="GRCh38",
            build_version="v1",
        )
        genome_complexity = GenomeComplexity.objects.create(
            size="3.2Gb",
            ploidy="Diploid",
            gc_content="41%",
        )
        sample_origin = SampleOrigin.objects.create(
            tissue="Lung",
            collection_method="Biopsy",
            storage_conditions="-80C",
        )
        material_type = MaterialType.objects.create(
            material_type="DNA",
            integrity_number="9.8",
        )
        library_construction = LibraryConstruction.objects.create(
            kit="MetaPrep",
            fragmentation="Acoustic",
            adapter_ligation_efficiency="95%",
        )
        illumina_seq = IlluminaSeq.objects.create(
            instrument="NovaSeq 6000",
            flow_cell="S4",
        )
        ont_seq = OntSeq.objects.create(
            instrument="PromethION",
            flow_cell_version="R10.4",
        )
        pacbio_seq = PacBioSeq.objects.create(
            instrument="Sequel II",
            smrt_cell_type="8M",
        )
        iontorrent_seq = IonTorrentSeq.objects.create(
            instrument="Ion S5",
            chip_type="530",
        )
        bioinfo_alignment = BioinfoAlignment.objects.create(
            tool="BWA",
            ref_genome_version="GRCh38",
        )
        bioinfo_variant_calling = BioinfoVariantCalling.objects.create(
            tool="GATK",
            version="4.2",
        )
        bioinfo_post_proc = BioinfoPostProc.objects.create(
            normalization="bcftools",
        )
        input_quality = InputQuality.objects.create(
            a260_a280=1.8,
            dna_concentration=15.2,
        )

        sample_group = SampleGroup.objects.create(
            name="Detail Cohort",
            created_by=owner.organization_profile,
            doi="10.1000/detail-cohort",
            source_lab="MetaLab",
            contact_email="lab@example.com",
            contact_phone="123-456-7890",
            total_samples=5,
            inclusion_criteria="Adults",
            exclusion_criteria="Under 18",
            comments="Rich metadata snapshot",
            reference_genome_build=reference,
            genome_complexity=genome_complexity,
            sample_origin=sample_origin,
            material_type=material_type,
            library_construction=library_construction,
            illumina_seq=illumina_seq,
            ont_seq=ont_seq,
            pacbio_seq=pacbio_seq,
            iontorrent_seq=iontorrent_seq,
            bioinfo_alignment=bioinfo_alignment,
            bioinfo_variant_calling=bioinfo_variant_calling,
            bioinfo_post_proc=bioinfo_post_proc,
            input_quality=input_quality,
        )

        info = Info.objects.create(
            af="0.5",
            additional={"clinvar_significance": "Pathogenic"},
        )
        fmt = Format.objects.create(
            genotype="0/1",
            payload={
                "fields": {"gq": "99"},
                "additional": {"sample_id": "Sample001"},
            },
        )

        allele = AlleleFrequency.objects.create(
            sample_group=sample_group,
            chrom="1",
            pos=123456,
            variant_id="rsDetail1",
            ref="A",
            alt="T",
            qual=42.0,
            filter="PASS",
            info=info,
            format=fmt,
        )

        return sample_group, allele


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


class EditProfileViewTests(TestCase):
    """Validate the edit profile workflow and form configuration."""

    def setUp(self) -> None:
        super().setUp()
        self.factory = RequestFactory()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="edit_user",
            password="test-pass-123",
            email="initial@example.com",
            first_name="Initial",
            last_name="User",
        )

    def test_form_kwargs_uses_authenticated_user_instance(self) -> None:
        request = self.factory.get(reverse("edit_profile"))
        request.user = self.user

        view = EditProfileView()
        view.setup(request)

        form_kwargs = view.get_form_kwargs()

        self.assertIs(form_kwargs["instance"], self.user)

    def test_successful_profile_update_redirects_and_persists(self) -> None:
        self.client.force_login(self.user)

        response = self.client.post(
            reverse("edit_profile"),
            {
                "username": "edit_user",
                "email": "updated@example.com",
                "first_name": "Updated",
                "last_name": "Name",
                "organization_name": "Updated Org",
            },
        )

        self.assertRedirects(response, reverse("profile"))

        self.user.refresh_from_db()
        self.assertEqual(self.user.email, "updated@example.com")
        self.assertEqual(self.user.first_name, "Updated")
        self.assertEqual(self.user.last_name, "Name")
        self.assertEqual(
            self.user.organization_profile.organization_name,
            "Updated Org",
        )

        messages = list(get_messages(response.wsgi_request))
        self.assertIn(
            "Your profile has been updated.",
            [message.message for message in messages],
        )


class SearchResultsViewTests(TestCase):
    """Validate search behaviour, filter context, and table population."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="search_user",
            password="search-pass",
            email="search@example.com",
        )

        self.kidney_origin = SampleOrigin.objects.create(
            tissue="Kidney",
            collection_method="Biopsy",
            storage_conditions="Cryogenic",
        )
        self.liver_origin = SampleOrigin.objects.create(
            tissue="Liver",
            collection_method="Surgical",
            storage_conditions="Room Temperature",
        )

        self.kidney_group = SampleGroup.objects.create(
            name="Kidney Cohort",
            sample_origin=self.kidney_origin,
            created_by=self.user.organization_profile,
        )
        self.liver_group = SampleGroup.objects.create(
            name="Liver Cohort",
            sample_origin=self.liver_origin,
            created_by=self.user.organization_profile,
        )

    def test_search_query_filters_expected_sample_group(self) -> None:
        response = self.client.get(reverse("search_results"), {"query": "Kidney"})

        self.assertEqual(response.status_code, 200)

        table = response.context["table"]
        self.assertEqual([row.record for row in table.rows], [self.kidney_group])

        sample_filter = response.context["filter"]
        self.assertIsInstance(sample_filter, SampleGroupFilter)
        self.assertEqual(sample_filter.data.get("query"), "Kidney")
        self.assertQuerySetEqual(
            sample_filter.qs,
            [self.kidney_group],
            transform=lambda group: group,
        )

        form = response.context["form"]
        self.assertIsInstance(form, SearchForm)
        self.assertEqual(form.data.get("query"), "Kidney")

    def test_empty_query_returns_all_records(self) -> None:
        response = self.client.get(reverse("search_results"), {"query": ""})

        self.assertEqual(response.status_code, 200)

        table = response.context["table"]
        self.assertCountEqual(
            [row.record for row in table.rows],
            [self.kidney_group, self.liver_group],
        )

        sample_filter = response.context["filter"]
        self.assertIsInstance(sample_filter, SampleGroupFilter)
        self.assertEqual(sample_filter.data.get("query"), "")
        self.assertQuerySetEqual(
            sample_filter.qs.order_by("pk"),
            SampleGroup.objects.order_by("pk"),
            transform=lambda group: group,
        )

        form = response.context["form"]
        self.assertIsInstance(form, SearchForm)
        self.assertEqual(form.data.get("query"), "")


class DashboardViewTests(TestCase):
    """Confirm the dashboard only surfaces recent data from the logged-in user."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()

        self.user_one = User.objects.create_user(
            username="dashboard_user_one",
            password="dashboard-pass-one",
            email="one@example.com",
        )
        self.user_two = User.objects.create_user(
            username="dashboard_user_two",
            password="dashboard-pass-two",
            email="two@example.com",
        )

        # Create more than six datasets for the first user to exercise the cap.
        self.user_one_groups = [
            SampleGroup.objects.create(
                name=f"User One Group {index}",
                created_by=self.user_one.organization_profile,
            )
            for index in range(8)
        ]

        self.user_two_groups = [
            SampleGroup.objects.create(
                name=f"User Two Group {index}",
                created_by=self.user_two.organization_profile,
            )
            for index in range(2)
        ]

        # Populate allele frequencies for each organization.
        self.user_one_actions = [
            AlleleFrequency.objects.create(
                sample_group=group,
                chrom="1",
                pos=index + 1,
                ref="A",
                alt="T",
            )
            for index, group in enumerate(self.user_one_groups)
        ]

        self.user_two_actions = [
            AlleleFrequency.objects.create(
                sample_group=group,
                chrom="2",
                pos=index + 10,
                ref="C",
                alt="G",
            )
            for index, group in enumerate(self.user_two_groups)
        ]

    def test_dashboard_limits_and_filters_datasets_and_actions(self) -> None:
        self.client.force_login(self.user_one)

        response = self.client.get(reverse("dashboard"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "dashboard.html")

        recent_datasets = response.context["recent_datasets"]
        expected_datasets = list(
            SampleGroup.objects.filter(
                created_by=self.user_one.organization_profile
            ).order_by("-pk")[:6]
        )

        self.assertEqual(recent_datasets, expected_datasets)
        self.assertTrue(all(dataset.created_by == self.user_one.organization_profile for dataset in recent_datasets))
        self.assertLessEqual(len(recent_datasets), 6)

        recent_actions = response.context["recent_actions"]
        expected_actions = list(
            AlleleFrequency.objects.filter(
                sample_group__created_by=self.user_one.organization_profile
            ).order_by("-pk")[:6]
        )

        self.assertEqual(recent_actions, expected_actions)
        self.assertTrue(
            all(
                action.sample_group.created_by == self.user_one.organization_profile
                for action in recent_actions
            )
        )
        self.assertLessEqual(len(recent_actions), 6)

    def test_dashboard_respects_organization_for_different_user(self) -> None:
        self.client.force_login(self.user_two)

        response = self.client.get(reverse("dashboard"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "dashboard.html")

        self.assertEqual(
            response.context["recent_datasets"],
            list(
                SampleGroup.objects.filter(
                    created_by=self.user_two.organization_profile
                ).order_by("-pk")[:6]
            ),
        )
        self.assertEqual(
            response.context["recent_actions"],
            list(
                AlleleFrequency.objects.filter(
                    sample_group__created_by=self.user_two.organization_profile
                ).order_by("-pk")[:6]
            ),
        )


class SampleGroupDetailViewTests(SampleGroupTestDataMixin, TestCase):
    """Exercise the sample group detail view and its access controls."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.owner = User.objects.create_user(
            username="detail_owner",
            password="detail-pass",
            email="owner@example.com",
        )
        self.other_user = User.objects.create_user(
            username="detail_other",
            password="detail-pass",
            email="other@example.com",
        )
        self.sample_group, self.allele = self.create_sample_group_with_variant(self.owner)

    def detail_url(self) -> str:
        for name in (
            "sample_group_detail",
            "profile_sample_group_detail",
            "sample-group-detail",
        ):
            try:
                return reverse(name, args=[self.sample_group.pk])
            except NoReverseMatch:
                continue
        self.fail("Sample group detail route is not configured")

    def test_owner_can_view_detail(self) -> None:
        self.client.force_login(self.owner)

        response = self.client.get(self.detail_url())

        self.assertEqual(response.status_code, 200)
        context_object = response.context.get("sample_group") or response.context.get(
            "object"
        )
        self.assertEqual(context_object, self.sample_group)

    def test_forbids_access_for_non_owner(self) -> None:
        self.client.force_login(self.other_user)

        response = self.client.get(self.detail_url())

        self.assertIn(response.status_code, {403, 404})

    def test_metadata_sections_and_variant_table_present(self) -> None:
        self.client.force_login(self.owner)

        response = self.client.get(self.detail_url())

        metadata_sections = response.context.get("metadata_sections")
        if metadata_sections is None:
            metadata_sections = response.context.get("metadata")
        self.assertIsNotNone(metadata_sections)

        if isinstance(metadata_sections, dict):
            sections = metadata_sections.values()
        else:
            sections = metadata_sections

        metadata_values = {
            str(value)
            for section in sections
            if isinstance(section, dict)
            for value in section.values()
            if value not in (None, "", [])
        }

        for expected_value in {
            self.sample_group.name,
            self.sample_group.source_lab,
            self.sample_group.contact_email,
            self.sample_group.reference_genome_build.build_name,
            self.sample_group.genome_complexity.size,
        }:
            self.assertIn(expected_value, metadata_values)

        variant_table = response.context.get("variant_table")
        if variant_table is None:
            variant_table = response.context.get("table")
        self.assertIsNotNone(variant_table)

        header_names = [column.name for column in variant_table.columns]
        expected_prefix = [
            "chrom",
            "pos",
            "variant_id",
            "ref",
            "alt",
            "qual",
            "filter",
        ]
        self.assertGreaterEqual(len(header_names), len(expected_prefix))
        self.assertEqual(header_names[: len(expected_prefix)], expected_prefix)

        table_html = variant_table.as_html(response.wsgi_request)
        self.assertIn(self.allele.variant_id, table_html)
        self.assertIn(str(self.allele.pos), table_html)
        self.assertIn(self.allele.ref, table_html)
        self.assertIn(self.allele.alt, table_html)


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
        self.assertEqual(allele.info.af, "0.5")
        self.assertEqual(allele.info.additional["clnsig"], "Pathogenic")
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.fields["gq"], "99")
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


class SampleGroupDeleteViewTests(TestCase):
    """Ensure the imported sample groups can be safely removed."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="deleter",
            email="deleter@example.com",
            password="pass12345",
        )
        self.other_user = User.objects.create_user(
            username="outsider",
            email="outsider@example.com",
            password="pass12345",
        )

        self.input_quality = InputQuality.objects.create(a260_a280=1.9)
        self.reference = ReferenceGenomeBuild.objects.create(
            build_name="GRCh38",
            build_version="v2",
        )
        self.sample_group = SampleGroup.objects.create(
            name="Deletable Group",
            created_by=self.user.organization_profile,
            contact_email="delete@example.com",
            input_quality=self.input_quality,
            reference_genome_build=self.reference,
        )

    def delete_url(self) -> str:
        return reverse("sample_group_delete", args=[self.sample_group.pk])

    def test_login_required(self) -> None:
        response = self.client.get(self.delete_url())
        login_url = reverse("login")
        self.assertEqual(response.status_code, 302)
        self.assertTrue(response.url.startswith(f"{login_url}?next="))

    def test_owner_sees_confirmation(self) -> None:
        self.client.force_login(self.user)
        response = self.client.get(self.delete_url())
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Delete Sample Group")
        self.assertContains(response, self.sample_group.name)

    def test_owner_can_delete_sample_group(self) -> None:
        self.client.force_login(self.user)
        response = self.client.post(self.delete_url(), follow=True)
        self.assertRedirects(response, reverse("profile"))
        self.assertFalse(SampleGroup.objects.filter(pk=self.sample_group.pk).exists())
        self.assertFalse(
            ReferenceGenomeBuild.objects.filter(pk=self.reference.pk).exists()
        )
        self.assertFalse(
            InputQuality.objects.filter(pk=self.input_quality.pk).exists()
        )

    def test_non_owner_cannot_delete_sample_group(self) -> None:
        self.client.force_login(self.other_user)
        response = self.client.post(self.delete_url())
        self.assertEqual(response.status_code, 404)
        self.assertTrue(SampleGroup.objects.filter(pk=self.sample_group.pk).exists())

    def test_shared_metadata_is_preserved(self) -> None:
        shared_reference = ReferenceGenomeBuild.objects.create(
            build_name="SharedRef",
            build_version="v1",
        )
        sibling = SampleGroup.objects.create(
            name="Sibling Group",
            created_by=self.user.organization_profile,
            reference_genome_build=shared_reference,
        )
        self.sample_group.reference_genome_build = shared_reference
        self.sample_group.save(update_fields=["reference_genome_build"])

        self.client.force_login(self.user)
        response = self.client.post(self.delete_url())
        self.assertRedirects(response, reverse("profile"))

        sibling.refresh_from_db()
        self.assertTrue(
            ReferenceGenomeBuild.objects.filter(pk=shared_reference.pk).exists()
        )


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
