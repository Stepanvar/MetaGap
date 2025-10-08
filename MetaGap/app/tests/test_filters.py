"""Tests for the filter classes defined in :mod:`app.filters`."""

from django.contrib.auth.models import User
from django.test import TestCase

from app import filters
from app.models import AlleleFrequency, SampleGroup, SampleOrigin


class AlleleFrequencyFilterTests(TestCase):
    """Validate the behaviour of :class:`app.filters.AlleleFrequencyFilter`."""

    def setUp(self):
        self.user = User.objects.create(username="filter-user")
        self.sample_group = SampleGroup.objects.create(
            name="Population A", created_by=self.user.organization_profile
        )

        self.record_chr1 = AlleleFrequency.objects.create(
            sample_group=self.sample_group,
            chrom="chr1",
            pos=12345,
            variant_id="rs123",
            ref="A",
            alt="T",
        )
        self.record_chr2 = AlleleFrequency.objects.create(
            sample_group=self.sample_group,
            chrom="chr2",
            pos=54321,
            variant_id="rs543",
            ref="G",
            alt="C",
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
