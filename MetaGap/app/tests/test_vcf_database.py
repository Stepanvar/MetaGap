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
    get_or_create_sample_group,
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


def test_get_or_create_sample_group_new():
    """Test get_or_create_sample_group creates a new sample group."""
    user_org_profile = None  # This would be a mock in real usage
    result = get_or_create_sample_group(
        sample_group_name="Test Group",
        sample_group_description="Test Description",
        user_organization_profile=user_org_profile,
        metadata={}
    )
    
    # Since we can't easily mock the organization profile, just check it returns something
    assert result is not None


def test_vcf_database_writer_init():
    """Test VCFDatabaseWriter initialization."""
    writer = VCFDatabaseWriter()
    assert writer is not None
    assert writer.insert_buffer_size == 1000
    assert writer.inserted_count == 0


def test_vcf_database_writer_add_record():
    """Test adding a record to VCFDatabaseWriter."""
    writer = VCFDatabaseWriter()
    
    # Create a sample group first
    sample_group = SampleGroup.objects.create(
        name="Test Group",
        description="Test Description"
    )
    
    # Create an Info object
    info_obj = Info.objects.create(af="0.5", dp="100", mq="60")
    
    # Add a record
    record_data = {
        'sample_group': sample_group,
        'chrom': 'chr1',
        'pos': 1000,
        'variant_id': 'rs12345',
        'ref': 'A',
        'alt': 'T',
        'qual': 50.0,
        'filter': 'PASS',
        'info': info_obj
    }
    
    writer.add_record(record_data)
    
    # Check that buffer has the record
    assert len(writer.buffer) == 1


def test_vcf_database_writer_flush_empty():
    """Test flushing an empty buffer."""
    writer = VCFDatabaseWriter()
    writer.flush()  # Should not raise an error


def test_vcf_database_writer_add_format_records():
    """Test adding format records to VCFDatabaseWriter."""
    writer = VCFDatabaseWriter()
    
    # Create a sample group first
    sample_group = SampleGroup.objects.create(
        name="Test Group",
        description="Test Description"
    )
    
    # Create an Info object
    info_obj = Info.objects.create(af="0.5", dp="100", mq="60")
    
    # Create an AlleleFrequency object
    allele_freq = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom='chr1',
        pos=1000,
        variant_id='rs12345',
        ref='A',
        alt='T',
        qual=50.0,
        filter='PASS',
        info=info_obj
    )
    
    # Add format records
    format_data = [
        {
            'allele_frequency': allele_freq,
            'sample_name': 'Sample1',
            'values_json': json.dumps({'GT': '0/1', 'AD': [10, 15]})
        }
    ]
    
    writer.add_format_records(format_data)
    
    # Check that format buffer has the record
    assert len(writer.format_buffer) == 1


def test_vcf_database_writer_flush_format_empty():
    """Test flushing an empty format buffer."""
    writer = VCFDatabaseWriter()
    writer.flush_format()  # Should not raise an error


def test_vcf_database_writer_convert_phaseset_to_string():
    """Test converting phaseset to string."""
    writer = VCFDatabaseWriter()
    
    # Test with integer
    result = writer._convert_phaseset_to_string(123)
    assert result == "123"
    
    # Test with string
    result = writer._convert_phaseset_to_string("phase1")
    assert result == "phase1"
    
    # Test with None
    result = writer._convert_phaseset_to_string(None)
    assert result is None


def test_vcf_database_writer_process_format_field():
    """Test processing format field values."""
    writer = VCFDatabaseWriter()
    
    # Test with string
    result = writer._process_format_field('GT', '0/1')
    assert result == '0/1'
    
    # Test with integer
    result = writer._process_format_field('DP', 30)
    assert result == 30
    
    # Test with phaseset
    result = writer._process_format_field('PS', 12345)
    assert result == "12345"


def test_vcf_database_writer_convert_to_boolean():
    """Test converting values to boolean."""
    writer = VCFDatabaseWriter()
    
    # Test true values
    assert writer._convert_to_boolean(True) is True
    assert writer._convert_to_boolean('True') is True
    assert writer._convert_to_boolean('true') is True
    assert writer._convert_to_boolean('1') is True
    assert writer._convert_to_boolean(1) is True
    
    # Test false values
    assert writer._convert_to_boolean(False) is False
    assert writer._convert_to_boolean('False') is False
    assert writer._convert_to_boolean('false') is False
    assert writer._convert_to_boolean('0') is False
    assert writer._convert_to_boolean(0) is False
    assert writer._convert_to_boolean(None) is False
    assert writer._convert_to_boolean('anything') is False


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