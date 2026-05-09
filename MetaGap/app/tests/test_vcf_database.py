"""Tests for the VCF database helpers defined in :mod:`app.services.vcf_database`."""

from __future__ import annotations

import json
from decimal import Decimal

import pytest
from django.core.exceptions import ValidationError
from django.db import IntegrityError

from app.models import AlleleFrequency, Format, Info, SampleGroup
from app.services.vcf_database import (
    VCFDatabaseWriter,
    convert_to_python_type,
    extract_info_field_metadata,
    extract_format_field_metadata,
    extract_contig_field_metadata,
    INFO_FIELD_FLOAT,
    INFO_FIELD_INT,
    INFO_FIELD_STRING,
    INFO_PLACEHOLDER_VALUES,
    INFO_FIELD_MAP,
    FORMAT_FIELD_MAP
)

pytestmark = pytest.mark.django_db


def test_convert_to_python_type_string():
    """Test convert_to_python_type with string type."""
    result = convert_to_python_type("hello", INFO_FIELD_STRING)
    assert result == "hello"


def test_convert_to_python_type_int():
    """Test convert_to_python_type with int type."""
    result = convert_to_python_type("123", INFO_FIELD_INT)
    assert result == 123


def test_convert_to_python_type_float():
    """Test convert_to_python_type with float type."""
    result = convert_to_python_type("12.34", INFO_FIELD_FLOAT)
    assert result == 12.34


def test_convert_to_python_type_empty_string():
    """Test convert_to_python_type with empty string."""
    result = convert_to_python_type("", INFO_FIELD_STRING)
    assert result == ""


def test_convert_to_python_type_placeholder():
    """Test convert_to_python_type with placeholder value."""
    result = convert_to_python_type(".", INFO_FIELD_STRING)
    assert result is None


def test_convert_to_python_type_invalid_int():
    """Test convert_to_python_type with invalid int."""
    result = convert_to_python_type("abc", INFO_FIELD_INT)
    assert result is None


def test_convert_to_python_type_invalid_float():
    """Test convert_to_python_type with invalid float."""
    result = convert_to_python_type("xyz", INFO_FIELD_FLOAT)
    assert result is None


def test_extract_info_field_metadata_basic():
    """Test extract_info_field_metadata with basic data."""
    header_info = {'AF': {'Number': 'A', 'Type': 'Float', 'Description': 'Allele Frequency'}}
    result = extract_info_field_metadata(header_info)
    assert 'af' in result


def test_extract_info_field_metadata_empty():
    """Test extract_info_field_metadata with empty dict."""
    result = extract_info_field_metadata({})
    assert result == {}


def test_extract_format_field_metadata_basic():
    """Test extract_format_field_metadata with basic data."""
    header_info = {'AD': {'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths'}}
    result = extract_format_field_metadata(header_info)
    assert 'ad' in result


def test_extract_format_field_metadata_empty():
    """Test extract_format_field_metadata with empty dict."""
    result = extract_format_field_metadata({})
    assert result == {}


def test_extract_contig_field_metadata_basic():
    """Test extract_contig_field_metadata with basic data."""
    header_info = {'1': {'length': 249250621, 'url': 'http://example.com/contig1'}}
    result = extract_contig_field_metadata(header_info)
    assert '1' in result


def test_extract_contig_field_metadata_empty():
    """Test extract_contig_field_metadata with empty dict."""
    result = extract_contig_field_metadata({})
    assert result == {}



def test_vcf_database_writer_init():
    """Test VCFDatabaseWriter initialization."""
    writer = VCFDatabaseWriter()
    assert writer is not None
    assert hasattr(writer, "metadata_field_aliases")
    assert hasattr(writer, "metadata_model_map")


def test_vcf_database_writer_serialize_alt():
    """VCFDatabaseWriter.serialize_alt joins alt alleles."""
    assert VCFDatabaseWriter.serialize_alt(["A", "T"]) == "A,T"
    assert VCFDatabaseWriter.serialize_alt([]) == ""
    assert VCFDatabaseWriter.serialize_alt(None) == ""


def test_vcf_database_writer_serialize_filter():
    """VCFDatabaseWriter.serialize_filter handles PASS and multi-value."""
    assert VCFDatabaseWriter.serialize_filter("PASS") == "PASS"
    assert VCFDatabaseWriter.serialize_filter(None) is None
    assert VCFDatabaseWriter.serialize_filter([]) is None


def test_vcf_database_writer_stringify():
    """VCFDatabaseWriter.stringify converts various types."""
    assert VCFDatabaseWriter.stringify(None) is None
    assert VCFDatabaseWriter.stringify("hello") == "hello"
    assert VCFDatabaseWriter.stringify([1, 2]) == "1,2"
    assert VCFDatabaseWriter.stringify(42) == "42"


def test_vcf_database_writer_create_info_empty():
    """create_info_instance returns None for empty/None info."""
    writer = VCFDatabaseWriter()
    assert writer.create_info_instance(None) is None
    assert writer.create_info_instance({}) is None


def test_vcf_database_writer_create_info_with_fields(db):
    """create_info_instance persists known INFO fields."""
    writer = VCFDatabaseWriter()
    info = writer.create_info_instance({"AF": "0.5", "DP": "100"})
    assert info is not None
    assert info.af == "0.5"
    assert info.dp == "100"


def test_vcf_database_writer_create_format_empty():
    """create_format_instance returns (None, None) for empty samples."""
    writer = VCFDatabaseWriter()
    result, name = writer.create_format_instance(None)
    assert result is None
    assert name is None


def test_info_field_maps():
    """Test the predefined INFO field maps."""
    assert 'af' in INFO_FIELD_MAP
    assert 'dp' in INFO_FIELD_MAP
    assert 'mq' in INFO_FIELD_MAP
    assert 'qd' in INFO_FIELD_MAP
    assert 'fs' in INFO_FIELD_MAP
    assert 'sor' in INFO_FIELD_MAP


def test_format_field_map():
    """Test the predefined FORMAT field map."""
    assert 'ad' in FORMAT_FIELD_MAP
    assert 'dp' in FORMAT_FIELD_MAP
    assert 'ft' in FORMAT_FIELD_MAP


def test_info_placeholder_values():
    """Test the INFO placeholder values."""
    assert '.' in INFO_PLACEHOLDER_VALUES
    assert '' in INFO_PLACEHOLDER_VALUES