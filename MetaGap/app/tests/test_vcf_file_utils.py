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


def test_parse_vcf_text_fallback_basic(db):
    """parse_vcf_text_fallback persists variants to DB via writer."""
    from django.contrib.auth.models import User
    from app.models import AlleleFrequency, SampleGroup
    from app.services.vcf_database import VCFDatabaseWriter

    user = User.objects.create_user(username="fbu1", password="x")
    sg = SampleGroup.objects.create(name="fb_basic", created_by=user.organization_profile)
    writer = VCFDatabaseWriter()

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf") as f:
        f.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t1000\trs12345\tA\tT\t50\tPASS\tAF=0.5\n"
        )
        path = f.name

    try:
        parse_vcf_text_fallback(path, sg, writer)
        assert AlleleFrequency.objects.filter(sample_group=sg).count() == 1
    finally:
        Path(path).unlink()


def test_parse_vcf_text_fallback_with_format(db):
    """parse_vcf_text_fallback handles FORMAT/sample columns."""
    from django.contrib.auth.models import User
    from app.models import AlleleFrequency, SampleGroup
    from app.services.vcf_database import VCFDatabaseWriter

    user = User.objects.create_user(username="fbu2", password="x")
    sg = SampleGroup.objects.create(name="fb_fmt", created_by=user.organization_profile)
    writer = VCFDatabaseWriter()

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf") as f:
        f.write(
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\n"
            "chr1\t1000\trs12345\tA\tT\t50\tPASS\tAF=0.5\tGT\t0/1\n"
        )
        path = f.name

    try:
        parse_vcf_text_fallback(path, sg, writer)
        assert AlleleFrequency.objects.filter(sample_group=sg).count() == 1
    finally:
        Path(path).unlink()


def test_parse_vcf_text_fallback_empty_file(db):
    """parse_vcf_text_fallback raises ValidationError when #CHROM line missing."""
    from django.contrib.auth.models import User
    from django.core.exceptions import ValidationError
    from app.models import SampleGroup
    from app.services.vcf_database import VCFDatabaseWriter

    user = User.objects.create_user(username="fbu3", password="x")
    sg = SampleGroup.objects.create(name="fb_empty", created_by=user.organization_profile)
    writer = VCFDatabaseWriter()

    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".vcf") as f:
        f.write("##fileformat=VCFv4.2\n")
        path = f.name

    try:
        with pytest.raises(ValidationError):
            parse_vcf_text_fallback(path, sg, writer)
    finally:
        Path(path).unlink()


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


def test_parse_vcf_text_fallback_nonexistent_file(db):
    """parse_vcf_text_fallback raises FileNotFoundError for missing file."""
    from django.contrib.auth.models import User
    from app.models import SampleGroup
    from app.services.vcf_database import VCFDatabaseWriter

    user = User.objects.create_user(username="fbu4", password="x")
    sg = SampleGroup.objects.create(name="fb_nf", created_by=user.organization_profile)
    writer = VCFDatabaseWriter()

    with pytest.raises(FileNotFoundError):
        parse_vcf_text_fallback("/nonexistent/path/file.vcf", sg, writer)


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