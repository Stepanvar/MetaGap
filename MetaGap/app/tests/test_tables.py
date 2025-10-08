"""Tests for dynamic table generation and rendering."""

from __future__ import annotations

import django_tables2 as tables
from django.contrib.auth import get_user_model
from django.contrib.auth.models import AnonymousUser
from django.test import RequestFactory, TestCase

from ..models import (
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    GenomeComplexity,
    IlluminaSeq,
    InputQuality,
    IonTorrentSeq,
    LibraryConstruction,
    MaterialType,
    OntSeq,
    PacBioSeq,
    PlatformIndependent,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)
from ..tables import create_dynamic_table
from ..views import SampleGroupTableView


class CreateDynamicTableTests(TestCase):
    """Validate the behaviour of the dynamic table factory."""

    def test_includes_related_fields_when_requested(self) -> None:
        table_class = create_dynamic_table(
            SampleGroup,
            table_name="SampleGroupRelatedTable",
            include_related=True,
        )

        self.assertTrue(issubclass(table_class, tables.Table))

        self.assertIn("name", table_class.Meta.fields)
        self.assertIn("reference_genome_build__build_name", table_class.Meta.fields)
        self.assertIn("created_by__organization_name", table_class.Meta.fields)
        self.assertNotIn("reference_genome_build", table_class.Meta.fields)

    def test_excludes_related_fields_when_disabled(self) -> None:
        table_class = create_dynamic_table(
            SampleGroup,
            table_name="SampleGroupFlatTable",
            include_related=False,
        )

        self.assertTrue(issubclass(table_class, tables.Table))

        self.assertIn("name", table_class.Meta.fields)
        self.assertIn("reference_genome_build", table_class.Meta.fields)
        self.assertNotIn("reference_genome_build__build_name", table_class.Meta.fields)


class SampleGroupTableViewTests(TestCase):
    """Ensure the SampleGroup table view renders safely with dynamic tables."""

    def setUp(self) -> None:
        super().setUp()
        self.factory = RequestFactory()
        user = get_user_model().objects.create_user(
            username="table_user",
            email="table_user@example.com",
            password="password123",
        )

        self.reference_genome_build = ReferenceGenomeBuild.objects.create(
            build_name="GRCh38",
            build_version="v1",
        )
        self.genome_complexity = GenomeComplexity.objects.create()
        self.sample_origin = SampleOrigin.objects.create()
        self.material_type = MaterialType.objects.create(material_type="DNA")
        self.library_construction = LibraryConstruction.objects.create(kit="Kit A")
        self.illumina_seq = IlluminaSeq.objects.create(instrument="NovaSeq")
        self.ont_seq = OntSeq.objects.create(instrument="PromethION")
        self.pacbio_seq = PacBioSeq.objects.create(instrument="Sequel")
        self.iontorrent_seq = IonTorrentSeq.objects.create(instrument="Ion S5")
        self.platform_independent = PlatformIndependent.objects.create(instrument="Generic")
        self.bioinfo_alignment = BioinfoAlignment.objects.create(tool="BWA")
        self.bioinfo_variant_calling = BioinfoVariantCalling.objects.create(tool="GATK")
        self.bioinfo_post_proc = BioinfoPostProc.objects.create(normalization="bcftools")
        self.input_quality = InputQuality.objects.create(a260_a280=1.8)

        self.sample_group = SampleGroup.objects.create(
            name="Rendered Cohort",
            created_by=user.organization_profile,
            reference_genome_build=self.reference_genome_build,
            genome_complexity=self.genome_complexity,
            sample_origin=self.sample_origin,
            material_type=self.material_type,
            library_construction=self.library_construction,
            illumina_seq=self.illumina_seq,
            ont_seq=self.ont_seq,
            pacbio_seq=self.pacbio_seq,
            iontorrent_seq=self.iontorrent_seq,
            platform_independent=self.platform_independent,
            bioinfo_alignment=self.bioinfo_alignment,
            bioinfo_variant_calling=self.bioinfo_variant_calling,
            bioinfo_post_proc=self.bioinfo_post_proc,
            input_quality=self.input_quality,
        )

    def test_view_renders_dynamic_table_html(self) -> None:
        request = self.factory.get("/sample-groups/")
        request.user = AnonymousUser()

        response = SampleGroupTableView.as_view()(request)
        response.render()

        self.assertEqual(response.status_code, 200)

        table = response.context_data["table"]
        html = table.as_html(request)

        self.assertIn(self.sample_group.name, html)
        self.assertIn(self.reference_genome_build.build_name, html)
