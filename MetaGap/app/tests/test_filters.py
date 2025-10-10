"""Tests for the filter classes defined in :mod:`app.filters`."""

from django.contrib.auth.models import User
from django.test import TestCase

from app import filters
from app.models import (
    AlleleFrequency,
    BioinfoAlignment,
    BioinfoVariantCalling,
    Info,
    SampleGroup,
    SampleOrigin,
)


class AlleleFrequencyFilterTests(TestCase):
    """Validate the behaviour of :class:`app.filters.AlleleFrequencyFilter`."""

    def setUp(self):
        self.user = User.objects.create(username="filter-user")
        self.sample_origin = SampleOrigin.objects.create(tissue="Liver")
        self.other_origin = SampleOrigin.objects.create(tissue="Blood")
        self.alignment_primary = BioinfoAlignment.objects.create(
            tool="AlignerA",
            params="--fast",
            ref_genome_version="GRCh38",
            recalibration_settings="BQSR",
        )
        self.alignment_secondary = BioinfoAlignment.objects.create(
            tool="AlignerB",
            params="--sensitive",
            ref_genome_version="GRCh37",
            recalibration_settings="None",
        )
        self.variant_primary = BioinfoVariantCalling.objects.create(
            tool="CallerA",
            version="1.0",
            filtering_thresholds="QD>10",
            duplicate_handling="Remove",
            mq="60",
        )
        self.variant_secondary = BioinfoVariantCalling.objects.create(
            tool="CallerB",
            version="2.1",
            filtering_thresholds="QD>15",
            duplicate_handling="Mark",
            mq="50",
        )

        self.sample_group = SampleGroup.objects.create(
            name="Population A",
            created_by=self.user.organization_profile,
            source_lab="North Lab",
            sample_origin=self.sample_origin,
            bioinfo_alignment=self.alignment_primary,
            bioinfo_variant_calling=self.variant_primary,
        )
        self.other_group = SampleGroup.objects.create(
            name="Population B",
            created_by=self.user.organization_profile,
            source_lab="South Lab",
            sample_origin=self.other_origin,
            bioinfo_alignment=self.alignment_secondary,
            bioinfo_variant_calling=self.variant_secondary,
        )

        self.info_chr1 = Info.objects.create(
            af="0.12",
            dp="25",
            mq="55",
            additional={"QD": "15.0", "FS": "0.4", "SOR": "1.5"},
        )
        self.info_chr2 = Info.objects.create(
            af=None,
            dp=None,
            mq=None,
            additional=None,
        )
        self.info_chr3 = Info.objects.create(
            af="0.45",
            dp="100",
            mq="70",
            additional={"QD": "20.5", "FS": "3.1", "SOR": "0.9"},
        )

        self.record_chr1 = AlleleFrequency.objects.create(
            sample_group=self.sample_group,
            chrom="chr1",
            pos=12345,
            variant_id="rs123",
            ref="A",
            alt="T",
            qual=75.0,
            filter="PASS",
            info=self.info_chr1,
        )
        self.record_chr2 = AlleleFrequency.objects.create(
            sample_group=self.sample_group,
            chrom="chr2",
            pos=54321,
            variant_id="rs543",
            ref="G",
            alt="C",
            qual=15.0,
            filter="LowQual",
            info=self.info_chr2,
        )
        self.record_chr3 = AlleleFrequency.objects.create(
            sample_group=self.other_group,
            chrom="chr3",
            pos=88888,
            variant_id="rs789",
            ref="C",
            alt="G",
            qual=120.0,
            filter="PASS",
            info=self.info_chr3,
        )

    def test_search_matches_chrom_field(self):
        qs = filters.AlleleFrequencyFilter(
            {"query": "chr1"}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr1])

    def test_search_matches_position_field(self):
        qs = filters.AlleleFrequencyFilter(
            {"query": "12345"}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr1])

    def test_search_matches_variant_id_field(self):
        qs = filters.AlleleFrequencyFilter(
            {"query": "rs543"}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr2])

    def test_chrom_filter_matches_exact_value(self):
        qs = filters.AlleleFrequencyFilter(
            {"chrom": "chr3"}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr3])

    def test_position_min_and_max_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"pos_min": 20000, "pos_max": 60000},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs), [self.record_chr2])

    def test_reference_and_alternate_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"ref": "G", "alt": "C"}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr2])

    def test_filter_pass_true_includes_only_pass_variants(self):
        qs = filters.AlleleFrequencyFilter(
            {"filter_pass": True}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr1, self.record_chr3])

    def test_filter_pass_false_excludes_pass_variants(self):
        qs = filters.AlleleFrequencyFilter(
            {"filter_pass": False}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr2])

    def test_qual_range_filters(self):
        qs_min = filters.AlleleFrequencyFilter(
            {"qual_min": 100}, queryset=AlleleFrequency.objects.all()
        ).qs
        qs_max = filters.AlleleFrequencyFilter(
            {"qual_max": 50}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs_min), [self.record_chr3])
        self.assertEqual(list(qs_max), [self.record_chr2])

    def test_af_range_filters_handle_placeholder_values(self):
        placeholder_info = Info.objects.create(
            af=".",
            dp=".",
            mq=".",
            additional={"QD": ".", "FS": ".", "SOR": "."},
        )
        placeholder_record = AlleleFrequency.objects.create(
            sample_group=self.sample_group,
            chrom="chrX",
            pos=11111,
            variant_id="rsPlaceholder",
            ref="A",
            alt="G",
            qual=5.0,
            filter="PASS",
            info=placeholder_info,
        )

        qs = filters.AlleleFrequencyFilter(
            {"af_min": 0.1, "af_max": 0.2}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr1])
        self.assertNotIn(self.record_chr2, qs)
        self.assertNotIn(placeholder_record, qs)

    def test_search_filter_sanitizes_non_numeric_placeholders(self):
        placeholder_info = Info.objects.create(
            af="N/A",
            dp="n/a",
            mq="NA",
            additional=None,
        )
        placeholder_record = AlleleFrequency.objects.create(
            sample_group=self.sample_group,
            chrom="chrY",
            pos=22222,
            variant_id="rsInvalid",
            ref="T",
            alt="C",
            qual=12.0,
            filter="PASS",
            info=placeholder_info,
        )

        qs = filters.AlleleFrequencySearchFilter(
            {"af_min": 0.2, "dp_min": 20, "mq_min": 50},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs), [self.record_chr3])
        self.assertNotIn(placeholder_record, qs)

    def test_dp_range_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"dp_min": 50}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr3])

    def test_mq_range_filters(self):
        qs_max = filters.AlleleFrequencyFilter(
            {"mq_max": 60}, queryset=AlleleFrequency.objects.all()
        ).qs
        qs_min = filters.AlleleFrequencyFilter(
            {"mq_min": 60}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs_max), [self.record_chr1])
        self.assertEqual(list(qs_min), [self.record_chr3])

    def test_qd_range_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"qd_min": 18, "qd_max": 19}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [])

        qs_max_only = filters.AlleleFrequencyFilter(
            {"qd_max": 16}, queryset=AlleleFrequency.objects.all()
        ).qs
        qs_min_only = filters.AlleleFrequencyFilter(
            {"qd_min": 18}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs_max_only), [self.record_chr1])
        self.assertEqual(list(qs_min_only), [self.record_chr3])

    def test_fs_range_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"fs_max": 1.0}, queryset=AlleleFrequency.objects.all()
        ).qs
        qs_min = filters.AlleleFrequencyFilter(
            {"fs_min": 2.0}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr1])
        self.assertEqual(list(qs_min), [self.record_chr3])

    def test_sor_range_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"sor_min": 1.0}, queryset=AlleleFrequency.objects.all()
        ).qs
        qs_max = filters.AlleleFrequencyFilter(
            {"sor_max": 1.0}, queryset=AlleleFrequency.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.record_chr1])
        self.assertEqual(list(qs_max), [self.record_chr3])

    def test_combined_numeric_and_string_filters(self):
        qs = filters.AlleleFrequencyFilter(
            {"chrom": "chr1", "af_min": 0.1, "qd_max": 16},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs), [self.record_chr1])

    def test_metadata_filters_for_source_lab(self):
        qs = filters.AlleleFrequencyFilter(
            {"sample_group_source_lab": "South"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs), [self.record_chr3])

    def test_metadata_filters_for_sample_origin(self):
        qs = filters.AlleleFrequencyFilter(
            {"sample_group_sample_origin_tissue": "Liver"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs), [self.record_chr1, self.record_chr2])

    def test_metadata_filters_for_variant_calling(self):
        qs_tool = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_variant_calling_tool": "CallerB"},
            queryset=AlleleFrequency.objects.all(),
        ).qs
        qs_version = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_variant_calling_version": "2.1"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs_tool), [self.record_chr3])
        self.assertEqual(list(qs_version), [self.record_chr3])

    def test_metadata_filters_for_alignment(self):
        qs_tool = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_alignment_tool": "AlignerA"},
            queryset=AlleleFrequency.objects.all(),
        ).qs
        qs_ref = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_alignment_ref_genome": "GRCh37"},
            queryset=AlleleFrequency.objects.all(),
        ).qs
        qs_recal = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_alignment_recalibration": "BQSR"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs_tool), [self.record_chr1, self.record_chr2])
        self.assertEqual(list(qs_ref), [self.record_chr3])
        self.assertEqual(list(qs_recal), [self.record_chr1, self.record_chr2])

    def test_metadata_and_numeric_filters_combined(self):
        qs = filters.AlleleFrequencyFilter(
            {"sample_group_source_lab": "South", "af_min": 0.4},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        self.assertEqual(list(qs), [self.record_chr3])


class SampleGroupFilterTests(TestCase):
    """Validate :class:`app.filters.SampleGroupFilter` lookups."""

    def setUp(self):
        self.user = User.objects.create(username="sample-user")
        self.sample_origin = SampleOrigin.objects.create(
            tissue="Kidney",
            collection_method="Biopsy",
            storage_conditions="Frozen",
            time_stored="2 years",
        )
        self.other_origin = SampleOrigin.objects.create(
            tissue="Blood",
            collection_method="Draw",
            storage_conditions="Refrigerated",
        )
        self.group = SampleGroup.objects.create(
            name="Kidney Cohort",
            created_by=self.user.organization_profile,
            sample_origin=self.sample_origin,
        )
        self.other_group = SampleGroup.objects.create(
            name="Control Cohort",
            created_by=self.user.organization_profile,
            sample_origin=self.other_origin,
        )

    def test_search_matches_sample_origin_tissue(self):
        qs = filters.SampleGroupFilter(
            {"query": "Kidney"}, queryset=SampleGroup.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.group])

    def test_search_matches_sample_origin_collection_method(self):
        qs = filters.SampleGroupFilter(
            {"query": "Biopsy"}, queryset=SampleGroup.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.group])

    def test_search_matches_nested_sample_origin_time_stored(self):
        qs = filters.SampleGroupFilter(
            {"query": "2 years"}, queryset=SampleGroup.objects.all()
        ).qs

        self.assertEqual(list(qs), [self.group])
