"""Tests for the filter classes defined in :mod:`app.filters`."""

from types import SimpleNamespace

import pytest
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


@pytest.fixture
def allele_frequency_data(db):
    """Create allele frequency records shared across tests."""

    user = User.objects.create(username="filter-user")
    sample_origin = SampleOrigin.objects.create(tissue="Liver")
    other_origin = SampleOrigin.objects.create(tissue="Blood")
    alignment_primary = BioinfoAlignment.objects.create(
        tool="AlignerA",
        params="--fast",
        ref_genome_version="GRCh38",
        recalibration_settings="BQSR",
    )
    alignment_secondary = BioinfoAlignment.objects.create(
        tool="AlignerB",
        params="--sensitive",
        ref_genome_version="GRCh37",
        recalibration_settings="None",
    )
    variant_primary = BioinfoVariantCalling.objects.create(
        tool="CallerA",
        version="1.0",
        filtering_thresholds="QD>10",
        duplicate_handling="Remove",
        mq="60",
    )
    variant_secondary = BioinfoVariantCalling.objects.create(
        tool="CallerB",
        version="2.1",
        filtering_thresholds="QD>15",
        duplicate_handling="Mark",
        mq="50",
    )

    sample_group = SampleGroup.objects.create(
        name="Population A",
        created_by=user.organization_profile,
        source_lab="North Lab",
        sample_origin=sample_origin,
        bioinfo_alignment=alignment_primary,
        bioinfo_variant_calling=variant_primary,
    )
    other_group = SampleGroup.objects.create(
        name="Population B",
        created_by=user.organization_profile,
        source_lab="South Lab",
        sample_origin=other_origin,
        bioinfo_alignment=alignment_secondary,
        bioinfo_variant_calling=variant_secondary,
    )

    info_chr1 = Info.objects.create(
        af="0.12",
        dp="25",
        mq="55",
        additional={"QD": "15.0", "FS": "0.4", "SOR": "1.5"},
    )
    info_chr2 = Info.objects.create(
        af=None,
        dp=None,
        mq=None,
        additional=None,
    )
    info_chr3 = Info.objects.create(
        af="0.45",
        dp="100",
        mq="70",
        additional={"QD": "20.5", "FS": "3.1", "SOR": "0.9"},
    )

    record_chr1 = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom="chr1",
        pos=12345,
        variant_id="rs123",
        ref="A",
        alt="T",
        qual=75.0,
        filter="PASS",
        info=info_chr1,
    )
    record_chr2 = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom="chr2",
        pos=54321,
        variant_id="rs543",
        ref="G",
        alt="C",
        qual=15.0,
        filter="LowQual",
        info=info_chr2,
    )
    record_chr3 = AlleleFrequency.objects.create(
        sample_group=other_group,
        chrom="chr3",
        pos=88888,
        variant_id="rs789",
        ref="C",
        alt="G",
        qual=120.0,
        filter="PASS",
        info=info_chr3,
    )

    return SimpleNamespace(
        user=user,
        sample_origin=sample_origin,
        other_origin=other_origin,
        alignment_primary=alignment_primary,
        alignment_secondary=alignment_secondary,
        variant_primary=variant_primary,
        variant_secondary=variant_secondary,
        sample_group=sample_group,
        other_group=other_group,
        info_chr1=info_chr1,
        info_chr2=info_chr2,
        info_chr3=info_chr3,
        record_chr1=record_chr1,
        record_chr2=record_chr2,
        record_chr3=record_chr3,
    )


@pytest.mark.django_db
class TestAlleleFrequencyFilter:
    """Validate the behaviour of :class:`app.filters.AlleleFrequencyFilter`."""

    def test_search_matches_chrom_field(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"query": "chr1"}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr1]

    def test_search_matches_position_field(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"query": "12345"}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr1]

    def test_search_matches_variant_id_field(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"query": "rs543"}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr2]

    def test_chrom_filter_matches_exact_value(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"chrom": "chr3"}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr3]

    def test_chrom_filter_does_not_match_partial_values(self, allele_frequency_data):
        data = allele_frequency_data
        info_chr10 = Info.objects.create(
            af="0.05",
            dp="30",
            mq="45",
            additional={"QD": "12.0", "FS": "1.0", "SOR": "0.8"},
        )
        info_chr1 = Info.objects.create(
            af="0.08",
            dp="40",
            mq="50",
            additional={"QD": "13.0", "FS": "1.2", "SOR": "0.7"},
        )
        record_chr1 = AlleleFrequency.objects.create(
            sample_group=data.sample_group,
            chrom="1",
            pos=11112,
            variant_id="rs1",
            ref="C",
            alt="T",
            qual=55.0,
            filter="PASS",
            info=info_chr1,
        )
        record_chr10 = AlleleFrequency.objects.create(
            sample_group=data.sample_group,
            chrom="10",
            pos=22222,
            variant_id="rs10",
            ref="T",
            alt="C",
            qual=60.0,
            filter="PASS",
            info=info_chr10,
        )

        exact_match_qs = filters.AlleleFrequencyFilter(
            {"chrom": "1"}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert record_chr1 in exact_match_qs
        assert record_chr10 not in exact_match_qs

    def test_chrom_filter_normalizes_input(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"chrom": "  ChR1  "}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr1]

    def test_position_exact_filter(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"pos": data.record_chr2.pos}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr2]

    def test_reference_and_alternate_filters(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"ref": "G", "alt": "C"}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr2]

    def test_filter_pass_true_includes_only_pass_variants(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"filter_pass": True}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr1, data.record_chr3]

    def test_filter_pass_false_excludes_pass_variants(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"filter_pass": False}, queryset=AlleleFrequency.objects.all()
        ).qs

        assert list(qs) == [data.record_chr2]

    def test_af_range_filters_handle_placeholder_values(self, allele_frequency_data):
        data = allele_frequency_data
        placeholder_info = Info.objects.create(
            af=".",
            dp=".",
            mq=".",
            additional={"QD": ".", "FS": ".", "SOR": "."},
        )
        placeholder_record = AlleleFrequency.objects.create(
            sample_group=data.sample_group,
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

        assert list(qs) == [data.record_chr1]
        assert data.record_chr2 not in qs
        assert placeholder_record not in qs

    @pytest.mark.parametrize(
        ("field", "params", "expected_records"),
        [
            pytest.param(
                "position",
                {"pos_min": 20000, "pos_max": 60000},
                ("record_chr2",),
                id="position-range",
            ),
            pytest.param(
                "qual-min",
                {"qual_min": 100},
                ("record_chr3",),
                id="qual-min",
            ),
            pytest.param(
                "qual-max",
                {"qual_max": 50},
                ("record_chr2",),
                id="qual-max",
            ),
            pytest.param(
                "dp-min",
                {"dp_min": 50},
                ("record_chr3",),
                id="dp-min",
            ),
            pytest.param(
                "mq-max",
                {"mq_max": 60},
                ("record_chr1",),
                id="mq-max",
            ),
            pytest.param(
                "mq-min",
                {"mq_min": 60},
                ("record_chr3",),
                id="mq-min",
            ),
            pytest.param(
                "qd-range",
                {"qd_min": 18, "qd_max": 19},
                (),
                id="qd-range-no-match",
            ),
            pytest.param(
                "qd-max",
                {"qd_max": 16},
                ("record_chr1",),
                id="qd-max-only",
            ),
            pytest.param(
                "qd-min",
                {"qd_min": 18},
                ("record_chr3",),
                id="qd-min-only",
            ),
            pytest.param(
                "fs-max",
                {"fs_max": 1.0},
                ("record_chr1",),
                id="fs-max",
            ),
            pytest.param(
                "fs-min",
                {"fs_min": 2.0},
                ("record_chr3",),
                id="fs-min",
            ),
            pytest.param(
                "sor-min",
                {"sor_min": 1.0},
                ("record_chr1",),
                id="sor-min",
            ),
            pytest.param(
                "sor-max",
                {"sor_max": 1.0},
                ("record_chr3",),
                id="sor-max",
            ),
        ],
    )
    def test_numeric_range_filters(
        self, allele_frequency_data, field, params, expected_records
    ):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            params, queryset=AlleleFrequency.objects.all()
        ).qs

        expected_objects = [getattr(data, name) for name in expected_records]
        assert [obj.pk for obj in qs] == [obj.pk for obj in expected_objects], (
            f"Unexpected records for filter '{field}'"
        )

    def test_combined_numeric_and_string_filters(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"chrom": "chr1", "af_min": 0.1, "qd_max": 16},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        assert list(qs) == [data.record_chr1]

    def test_metadata_filters_for_source_lab(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"sample_group_source_lab": "South"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        assert list(qs) == [data.record_chr3]

    def test_metadata_filters_for_sample_origin(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"sample_group_sample_origin_tissue": "Liver"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        assert list(qs) == [data.record_chr1, data.record_chr2]

    def test_metadata_filters_for_variant_calling(self, allele_frequency_data):
        data = allele_frequency_data
        qs_tool = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_variant_calling_tool": "CallerB"},
            queryset=AlleleFrequency.objects.all(),
        ).qs
        qs_version = filters.AlleleFrequencyFilter(
            {"sample_group_bioinfo_variant_calling_version": "2.1"},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        assert list(qs_tool) == [data.record_chr3]
        assert list(qs_version) == [data.record_chr3]

    def test_metadata_filters_for_alignment(self, allele_frequency_data):
        data = allele_frequency_data
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

        assert list(qs_tool) == [data.record_chr1, data.record_chr2]
        assert list(qs_ref) == [data.record_chr3]
        assert list(qs_recal) == [data.record_chr1, data.record_chr2]

    def test_metadata_and_numeric_filters_combined(self, allele_frequency_data):
        data = allele_frequency_data
        qs = filters.AlleleFrequencyFilter(
            {"sample_group_source_lab": "South", "af_min": 0.4},
            queryset=AlleleFrequency.objects.all(),
        ).qs

        assert list(qs) == [data.record_chr3]


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
