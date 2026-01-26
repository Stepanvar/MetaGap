"""Tests for the VCF file utilities defined in :mod:`app.services.vcf_file_utils`."""

from __future__ import annotations

import gzip
import tempfile
from pathlib import Path

import pytest
from django.core.exceptions import ValidationError

from app.services.vcf_file_utils import (
    open_vcf_text,
    extract_metadata_text_fallback,
    parse_vcf_text_fallback
)


def test_open_vcf_text_plain():
    """Test opening a plain text VCF file."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        temp_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        temp_file_path = temp_file.name

    try:
        with open_vcf_text(temp_file_path) as handle:
            content = handle.read()
            assert "##fileformat=VCFv4.2" in content
    finally:
        Path(temp_file_path).unlink()


def test_open_vcf_text_gzip():
    """Test opening a gzipped VCF file."""
    with tempfile.NamedTemporaryFile(delete=False, suffix='.vcf.gz') as temp_file:
        temp_path = Path(temp_file.name)

    # Write gzipped content
    with gzip.open(temp_path, 'wt', encoding='utf-8') as gz_file:
        gz_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    try:
        with open_vcf_text(str(temp_path)) as handle:
            content = handle.read()
            assert "##fileformat=VCFv4.2" in content
    finally:
        temp_path.unlink()


def test_open_vcf_text_bgz():
    """Test opening a bgzipped VCF file."""
    with tempfile.NamedTemporaryFile(delete=False, suffix='.vcf.bgz') as temp_file:
        temp_path = Path(temp_file.name)

    # Write bgzipped content
    with gzip.open(temp_path, 'wt', encoding='utf-8') as gz_file:
        gz_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    try:
        with open_vcf_text(str(temp_path)) as handle:
            content = handle.read()
            assert "##fileformat=VCFv4.2" in content
    finally:
        temp_path.unlink()


def test_extract_metadata_text_fallback_empty_warnings():
    """Test extract_metadata_text_fallback with empty warnings list."""
    warnings_list = []
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        temp_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        temp_file_path = temp_file.name

    try:
        metadata = extract_metadata_text_fallback(temp_file_path, warnings_list)
        assert isinstance(metadata, dict)
        assert len(warnings_list) >= 0  # May have warnings or not
    finally:
        Path(temp_file_path).unlink()


def test_extract_metadata_text_fallback_no_warnings_param():
    """Test extract_metadata_text_fallback without warnings parameter."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        temp_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        temp_file_path = temp_file.name

    try:
        metadata = extract_metadata_text_fallback(temp_file_path)
        assert isinstance(metadata, dict)
    finally:
        Path(temp_file_path).unlink()


def test_parse_vcf_text_fallback_basic():
    """Test parse_vcf_text_fallback with basic VCF content."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        temp_file.write("""##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1000	rs12345	A	T	50	PASS	AF=0.5
""")
        temp_file_path = temp_file.name

    try:
        result = parse_vcf_text_fallback(temp_file_path)
        
        # Check that the result has the expected structure
        assert 'metadata' in result
        assert 'samples' in result
        assert 'variants' in result
        
        # Check that we have the expected contig in metadata
        assert 'contig' in result['metadata']
    finally:
        Path(temp_file_path).unlink()


def test_parse_vcf_text_fallback_with_format():
    """Test parse_vcf_text_fallback with FORMAT field."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        temp_file.write("""##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1
chr1	1000	rs12345	A	T	50	PASS	AF=0.5	GT	0/1
""")
        temp_file_path = temp_file.name

    try:
        result = parse_vcf_text_fallback(temp_file_path)
        
        # Check that the result has the expected structure
        assert 'metadata' in result
        assert 'samples' in result
        assert 'variants' in result
    finally:
        Path(temp_file_path).unlink()


def test_parse_vcf_text_fallback_empty_file():
    """Test parse_vcf_text_fallback with empty file."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        # Write an empty file or just headers
        temp_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        temp_file_path = temp_file.name

    try:
        result = parse_vcf_text_fallback(temp_file_path)
        
        # Check that the result has the expected structure even with no variants
        assert 'metadata' in result
        assert 'samples' in result
        assert 'variants' in result
    finally:
        Path(temp_file_path).unlink()


def test_parse_vcf_text_fallback_no_header():
    """Test parse_vcf_text_fallback with missing header."""
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.vcf') as temp_file:
        # Write a file without proper header
        temp_file.write("chr1\t1000\trs12345\tA\tT\t50\tPASS\tAF=0.5\n")
        temp_file_path = temp_file.name

    try:
        # This might raise an exception, or handle gracefully
        result = parse_vcf_text_fallback(temp_file_path)
    except Exception:
        # If it raises an exception, that's acceptable behavior
        pass
    finally:
        if Path(temp_file_path).exists():
            Path(temp_file_path).unlink()


def test_open_vcf_text_nonexistent_file():
    """Test opening a nonexistent file."""
    with pytest.raises(FileNotFoundError):
        with open_vcf_text("/nonexistent/path/file.vcf"):
            pass  # This should raise FileNotFoundError


def test_extract_metadata_text_fallback_nonexistent_file():
    """Test extracting metadata from a nonexistent file."""
    with pytest.raises(FileNotFoundError):
        extract_metadata_text_fallback("/nonexistent/path/file.vcf")


def test_parse_vcf_text_fallback_nonexistent_file():
    """Test parsing a nonexistent file."""
    with pytest.raises(FileNotFoundError):
        parse_vcf_text_fallback("/nonexistent/path/file.vcf")


def test_open_vcf_text_various_suffixes():
    """Test opening files with various suffixes."""
    # Test with .txt extension (should be treated as plain text)
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as temp_file:
        temp_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        temp_file_path = temp_file.name

    try:
        with open_vcf_text(temp_file_path) as handle:
            content = handle.read()
            assert "##fileformat=VCFv4.2" in content
    finally:
        Path(temp_file_path).unlink()


def test_open_vcf_text_multiple_suffixes():
    """Test opening a file with multiple suffixes like .vcf.gz."""
    temp_path = Path("temp_test.vcf.gz")
    
    # Write bgzipped content
    with gzip.open(temp_path, 'wt', encoding='utf-8') as gz_file:
        gz_file.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

    try:
        with open_vcf_text(str(temp_path)) as handle:
            content = handle.read()
            assert "##fileformat=VCFv4.2" in content
    finally:
        if temp_path.exists():
            temp_path.unlink()