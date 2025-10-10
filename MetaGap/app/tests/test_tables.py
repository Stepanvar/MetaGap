"""Tests for dynamic table generation and rendering."""

from __future__ import annotations

import django_tables2 as tables
from django.contrib.auth import get_user_model
from django.contrib.auth.models import AnonymousUser
from django.test import RequestFactory, TestCase

from ..models import (
    AlleleFrequency,
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
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)
from ..tables import build_allele_frequency_table, create_dynamic_table
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


class AlleleFrequencyTableBuilderTests(TestCase):
    """Ensure the allele frequency table helper enforces shared defaults."""

    def test_defaults_prioritise_core_columns(self) -> None:
        table_class = build_allele_frequency_table()

        self.assertTrue(issubclass(table_class, tables.Table))

        expected_prefix = [
            "chrom",
            "pos",
            "ref",
            "alt",
            "qual",
            "filter",
            "info__af",
            "info__ac",
            "info__an",
            "info__dp",
            "info__mq",
        ]
        self.assertEqual(list(table_class.Meta.fields[: len(expected_prefix)]), expected_prefix)
        self.assertNotIn("info__additional", table_class.Meta.fields)
        self.assertNotIn("format__payload", table_class.Meta.fields)

    def test_allows_view_specific_overrides(self) -> None:
        table_class = build_allele_frequency_table(
            priority_extra=("variant_id",),
            exclude_extra=("info__id",),
        )

        fields = list(table_class.Meta.fields)
        self.assertEqual(
            fields[:12],
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "qual",
                "filter",
                "info__af",
                "info__ac",
                "info__an",
                "info__dp",
                "info__mq",
                "variant_id",
            ],
        )
        self.assertNotIn("info__id", fields)

    def test_sample_group_header_uses_friendly_label(self) -> None:
        user = get_user_model().objects.create_user(
            username="friendly_user",
            email="friendly@example.com",
            password="password123",
        )
        sample_group = SampleGroup.objects.create(
            name="Friendly Cohort",
            created_by=user.organization_profile,
            source_lab="Genome Lab",
        )
        AlleleFrequency.objects.create(
            sample_group=sample_group,
            chrom="1",
            pos=12345,
            ref="A",
            alt="G",
        )

        table_class = build_allele_frequency_table()
        table = table_class(AlleleFrequency.objects.all())

        factory = RequestFactory()
        request = factory.get("/allele-table/")
        html = table.as_html(request)

        self.assertIn(">Sample group<", html)
