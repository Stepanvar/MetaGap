"""Additional tests for the filter classes defined in :mod:`app.filters` to increase coverage."""

from __future__ import annotations

from types import SimpleNamespace

import pytest
from django.contrib.auth.models import User

from app import filters
from app.models import (
    AlleleFrequency,
    BioinfoAlignment,
    BioinfoVariantCalling,
    Info,
    SampleGroup,
    SampleOrigin,
)

pytestmark = pytest.mark.django_db


# ----------------------------- Fixtures --------------------------------- #
@pytest.fixture
def additional_filter_data(db):
    """Create additional data for filter tests."""
    user = User.objects.create_user(username="filter-user-additional", password="x")

    # Sample origins
    sample_origin = SampleOrigin.objects.create(tissue="Brain")
    other_origin = SampleOrigin.objects.create(tissue="Muscle")

    # Alignment & variant-calling metadata for groups
    alignment_primary = BioinfoAlignment.objects.create(
        tool="AlignerC",
        params="--accurate",
        ref_genome_version="GRCh39",
        recalibration_settings="BQSR",
    )
    alignment_secondary = BioinfoAlignment.objects.create(
        tool="AlignerD",
        params="--fast",
        ref_genome_version="GRCh38",
        recalibration_settings="None",
    )
    variant_primary = BioinfoVariantCalling.objects.create(
        tool="CallerC",
        version="3.0",
        filtering_thresholds="QD>20",
        duplicate_handling="Remove",
        mq="70",
    )
    variant_secondary = BioinfoVariantCalling.objects.create(
        tool="CallerD",
        version="4.0",
        filtering_thresholds="QD>25",
        duplicate_handling="Mark",
        mq="80",
    )

    # Sample groups
    sample_group = SampleGroup.objects.create(
        name="Population C",
        created_by=user.organization_profile,
        source_lab="East Lab",
        sample_origin=sample_origin,
        bioinfo_alignment=alignment_primary,
        bioinfo_variant_calling=variant_primary,
    )
    other_group = SampleGroup.objects.create(
        name="Population D",
        created_by=user.organization_profile,
        source_lab="West Lab",
        sample_origin=other_origin,
        bioinfo_alignment=alignment_secondary,
        bioinfo_variant_calling=variant_secondary,
    )

    # INFO payloads
    info_chr1 = Info.objects.create(
        af="0.15",
        dp="30",
        mq="60",
        additional={"QD": "18.0", "FS": "0.5", "SOR": "1.8"},
    )
    info_chr2 = Info.objects.create(af=None, dp=None, mq=None, additional=None)
    info_chr3 = Info.objects.create(
        af="0.50",
        dp="110",
        mq="75",
        additional={"QD": "22.0", "FS": "4.0", "SOR": "1.0"},
    )

    # Variant records
    record_chr1 = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom="chr1",
        pos=15000,
        variant_id="rs150",
        ref="A",
        alt="G",
        qual=80.0,
        filter="PASS",
        info=info_chr1,
    )
    record_chr2 = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom="chr2",
        pos=60000,
        variant_id="rs600",
        ref="T",
        alt="C",
        qual=20.0,
        filter="LowQual",
        info=info_chr2,
    )
    record_chr3 = AlleleFrequency.objects.create(
        sample_group=other_group,
        chrom="chr4",
        pos=90000,
        variant_id="rs900",
        ref="C",
        alt="A",
        qual=125.0,
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


# Additional tests for AlleleFrequencyFilter
def test_universal_search_with_digit_position(additional_filter_data):
    """Test universal search with digit position."""
    data = additional_filter_data
    qs = filters.AlleleFrequencyFilter(
        {"query": "15000"}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert list(qs) == [data.record_chr1]


def test_universal_search_empty_value():
    """Test universal search with empty value."""
    qs = filters.AlleleFrequencyFilter(
        {"query": ""}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_normalize_chrom_none():
    """Test _normalize_chrom with None value."""
    filter_instance = filters.AlleleFrequencyFilter()
    result = filter_instance._normalize_chrom(None)
    assert result is None


def test_normalize_chrom_empty():
    """Test _normalize_chrom with empty value."""
    filter_instance = filters.AlleleFrequencyFilter()
    result = filter_instance._normalize_chrom("")
    assert result is None


def test_normalize_chrom_whitespace():
    """Test _normalize_chrom with whitespace."""
    filter_instance = filters.AlleleFrequencyFilter()
    result = filter_instance._normalize_chrom("  ")
    assert result is None


def test_filter_chrom_none_value(additional_filter_data):
    """Test filter_chrom with None value."""
    data = additional_filter_data
    qs = filters.AlleleFrequencyFilter(
        {"chrom": None}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_chrom_empty_value(additional_filter_data):
    """Test filter_chrom with empty value."""
    data = additional_filter_data
    qs = filters.AlleleFrequencyFilter(
        {"chrom": ""}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_pass_status_none_value(additional_filter_data):
    """Test filter_pass_status with None value."""
    data = additional_filter_data
    qs = filters.AlleleFrequencyFilter(
        {"filter_pass": None}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_pass_status_false(additional_filter_data):
    """Test filter_pass_status with False value."""
    data = additional_filter_data
    qs = filters.AlleleFrequencyFilter(
        {"filter_pass": False}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert data.record_chr2 in list(qs)
    assert data.record_chr1 not in list(qs)


def test_filter_info_numeric_invalid_parts(additional_filter_data):
    """Test filter_info_numeric with invalid parts."""
    data = additional_filter_data
    # Testing with an invalid field name that doesn't match the expected pattern
    filter_instance = filters.AlleleFrequencyFilter(data={'invalid_field': 10})
    result = filter_instance.qs
    # Should not cause errors


def test_filter_info_numeric_none_value(additional_filter_data):
    """Test filter_info_numeric with None value."""
    data = additional_filter_data
    # Test af_min filter with None value
    qs = filters.AlleleFrequencyFilter(
        {"af_min": None}, queryset=AlleleFrequency.objects.all()
    ).qs
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_info_numeric_invalid_field(additional_filter_data):
    """Test filter_info_numeric with invalid field."""
    data = additional_filter_data
    # Creating a custom filter instance to test invalid field
    filter_instance = filters.AlleleFrequencyFilter()
    qs = filter_instance.filter_info_numeric(
        queryset=AlleleFrequency.objects.all(),
        name="invalid_min",  # Invalid field name
        value=0.1
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


# Tests for AlleleFrequencySearchFilter
def test_allele_frequency_search_filter_init():
    """Test initialization of AlleleFrequencySearchFilter."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    assert filter_instance is not None


def test_filter_query_empty_value():
    """Test filter_query with empty value."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    qs = filter_instance.filter_query(
        queryset=AlleleFrequency.objects.all(),
        name="query",
        value=""
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_query_whitespace_value():
    """Test filter_query with whitespace value."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    qs = filter_instance.filter_query(
        queryset=AlleleFrequency.objects.all(),
        name="query",
        value="   "
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_pass_only_true_values(additional_filter_data):
    """Test filter_pass_only with various true values."""
    data = additional_filter_data
    
    # Test with True
    filter_instance = filters.AlleleFrequencySearchFilter({"pass_only": True})
    result = filter_instance.filter_pass_only(
        queryset=AlleleFrequency.objects.all(),
        name="pass_only",
        value=True
    )
    assert data.record_chr1 in list(result)
    assert data.record_chr2 not in list(result)


def test_filter_pass_only_false_values(additional_filter_data):
    """Test filter_pass_only with various false values."""
    data = additional_filter_data
    
    # Test with False
    filter_instance = filters.AlleleFrequencySearchFilter({"pass_only": False})
    result = filter_instance.filter_pass_only(
        queryset=AlleleFrequency.objects.all(),
        name="pass_only",
        value=False
    )
    assert list(result) == list(AlleleFrequency.objects.all())


def test_filter_info_float_empty_value(additional_filter_data):
    """Test filter_info_float with empty value."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    qs = filter_instance.filter_info_float(
        queryset=AlleleFrequency.objects.all(),
        name="af_min",
        value=""
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_info_float_none_value(additional_filter_data):
    """Test filter_info_float with None value."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    qs = filter_instance.filter_info_float(
        queryset=AlleleFrequency.objects.all(),
        name="af_min",
        value=None
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_info_int_empty_value(additional_filter_data):
    """Test filter_info_int with empty value."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    qs = filter_instance.filter_info_int(
        queryset=AlleleFrequency.objects.all(),
        name="dp_min",
        value=""
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_filter_info_int_none_value(additional_filter_data):
    """Test filter_info_int with None value."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    qs = filter_instance.filter_info_int(
        queryset=AlleleFrequency.objects.all(),
        name="dp_min",
        value=None
    )
    assert list(qs) == list(AlleleFrequency.objects.all())


def test_apply_bootstrap_widgets_checkbox():
    """Test _apply_bootstrap_widgets with checkbox input."""
    filter_instance = filters.AlleleFrequencySearchFilter()
    # This should apply without errors
    filter_instance._apply_bootstrap_widgets()


def test_merge_css_classes():
    """Test _merge_css_classes static method."""
    result = filters.AlleleFrequencySearchFilter._merge_css_classes("existing-class", "new-class")
    assert "existing-class" in result
    assert "new-class" in result


def test_merge_css_classes_duplicate():
    """Test _merge_css_classes with duplicate class."""
    result = filters.AlleleFrequencySearchFilter._merge_css_classes("existing-class", "existing-class")
    assert result.count("existing-class") == 1


def test_merge_css_classes_empty_existing():
    """Test _merge_css_classes with empty existing classes."""
    result = filters.AlleleFrequencySearchFilter._merge_css_classes("", "new-class")
    assert result == "new-class"


# Tests for SampleGroupFilter
def test_sample_group_universal_search_empty():
    """Test SampleGroupFilter universal search with empty value."""
    qs = filters.SampleGroupFilter(
        {"query": ""}, queryset=SampleGroup.objects.all()
    ).qs
    assert list(qs) == list(SampleGroup.objects.all())


def test_sample_group_universal_search_whitespace():
    """Test SampleGroupFilter universal search with whitespace value."""
    qs = filters.SampleGroupFilter(
        {"query": "   "}, queryset=SampleGroup.objects.all()
    ).qs
    assert list(qs) == list(SampleGroup.objects.all())