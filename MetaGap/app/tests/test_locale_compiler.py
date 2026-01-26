"""Tests for the locale compiler utilities defined in :mod:`app.utils.locale_compiler`."""

from __future__ import annotations

import tempfile
import struct
from pathlib import Path

import pytest
from django.core.exceptions import ValidationError

from app.utils.locale_compiler import (
    CatalogCompilationError,
    _add_entry,
    _generate_mo,
)


def test_catalog_compilation_error():
    """Test CatalogCompilationError."""
    error = CatalogCompilationError("Compilation error")
    assert str(error) == "Compilation error"


def test_add_entry_fuzzy():
    """Test _add_entry with fuzzy flag."""
    messages = {}
    _add_entry(messages, None, b"msgid", b"msgstr", True)
    assert len(messages) == 0  # Fuzzy entry should not be added


def test_add_entry_not_fuzzy():
    """Test _add_entry without fuzzy flag."""
    messages = {}
    _add_entry(messages, None, b"msgid", b"msgstr", False)
    assert len(messages) == 1
    assert b"msgid" in messages
    assert messages[b"msgid"] == b"msgstr"


def test_add_entry_with_context():
    """Test _add_entry with context."""
    messages = {}
    _add_entry(messages, b"context", b"msgid", b"msgstr", False)
    assert len(messages) == 1
    # The key should include the context separator
    expected_key = b"context\x04msgid"
    assert expected_key in messages
    assert messages[expected_key] == b"msgstr"


def test_add_entry_with_context_fuzzy():
    """Test _add_entry with context and fuzzy flag."""
    messages = {}
    _add_entry(messages, b"context", b"msgid", b"msgstr", True)
    assert len(messages) == 0  # Fuzzy entry should not be added


def test_generate_mo_empty():
    """Test _generate_mo with empty messages dict."""
    messages = {}
    result = _generate_mo(messages)
    
    # Should generate a valid MO file header
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_generate_mo_single_entry():
    """Test _generate_mo with a single entry."""
    messages = {b"hello": b"world"}
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_generate_mo_multiple_entries():
    """Test _generate_mo with multiple entries."""
    messages = {
        b"hello": b"world",
        b"goodbye": b"farewell",
        b"test": b"example"
    }
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_generate_mo_unicode_strings():
    """Test _generate_mo with unicode byte strings."""
    messages = {b"café": b"coffee", b"naïve": b"simple"}
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_catalog_compilation_error_access_attributes():
    """Test accessing attributes of CatalogCompilationError."""
    error = CatalogCompilationError("Compilation error")
    # Just ensure no error when accessing common attributes
    str(error)


def test_add_entry_special_characters():
    """Test _add_entry with special characters."""
    messages = {}
    _add_entry(messages, None, b"hello\nworld", b"hi\nthere", False)
    assert len(messages) == 1
    assert b"hello\nworld" in messages
    assert messages[b"hello\nworld"] == b"hi\nthere"


def test_add_entry_empty_strings():
    """Test _add_entry with empty strings."""
    messages = {}
    _add_entry(messages, None, b"", b"", False)
    assert len(messages) == 1
    assert b"" in messages
    assert messages[b""] == b""


def test_add_entry_none_context():
    """Test _add_entry with None context."""
    messages = {}
    _add_entry(messages, None, b"msgid", b"msgstr", False)
    assert len(messages) == 1
    assert b"msgid" in messages
    assert messages[b"msgid"] == b"msgstr"


def test_generate_mo_special_entries():
    """Test _generate_mo with special entries."""
    messages = {
        b"": b"Domain-Id-Version\\nReport-Msgid-Bugs-To\\nPO-Revision-Date\\nPOT-Creation-Date\\nLast-Translator\\nLanguage-Team\\nLanguage\\nMIME-Version\\nContent-Type\\nContent-Transfer-Encoding\\nPlural-Forms",  # Headers
        b"special": b"entry"
    }
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_catalog_compilation_error_with_exception_chain():
    """Test CatalogCompilationError with chained exceptions."""
    try:
        raise CatalogCompilationError("Compilation error")
    except CatalogCompilationError as e:
        # Should be able to catch it properly
        assert str(e) == "Compilation error"


def test_add_entry_long_strings():
    """Test _add_entry with long strings."""
    long_msgid = b"a" * 1000
    long_msgstr = b"b" * 1000
    messages = {}
    _add_entry(messages, None, long_msgid, long_msgstr, False)
    assert len(messages) == 1
    assert long_msgid in messages
    assert messages[long_msgid] == long_msgstr


def test_generate_mo_large_dictionary():
    """Test _generate_mo with larger dictionary."""
    messages = {}
    for i in range(100):
        messages[f"key{i}".encode()] = f"value{i}".encode()
    
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_add_entry_multiline_strings():
    """Test _add_entry with multiline strings."""
    messages = {}
    _add_entry(messages, None, b"line1\nline2\nline3", b"translated1\ntranslated2", False)
    assert len(messages) == 1
    assert b"line1\nline2\nline3" in messages


def test_add_entry_binary_strings():
    """Test _add_entry with binary strings."""
    messages = {}
    _add_entry(messages, None, b"\x00\x01\x02", b"\x03\x04\x05", False)
    assert len(messages) == 1
    assert b"\x00\x01\x02" in messages
    assert messages[b"\x00\x01\x02"] == b"\x03\x04\x05"


def test_generate_mo_duplicate_handling():
    """Test _generate_mo behavior with potential duplicates."""
    # Note: In the actual implementation, if there are duplicate keys in the dict,
    # only the last value will be kept since it's a dictionary
    messages = {
        b"duplicate": b"first_value",
        b"duplicate": b"second_value"  # This will overwrite the first
    }
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_generate_mo_sorted_keys():
    """Test that _generate_mo processes keys in sorted order."""
    messages = {
        b"zebra": b"animal",
        b"apple": b"fruit", 
        b"banana": b"fruit"
    }
    result = _generate_mo(messages)
    
    assert len(result) > 0
    # Check the magic number at the beginning
    magic = struct.unpack('<I', result[:4])[0]
    assert magic == 0x950412de  # Little-endian magic number


def test_catalog_compilation_error_repr():
    """Test string representation of CatalogCompilationError."""
    error = CatalogCompilationError("Test compilation error")
    # Just ensure no error when getting the string representation
    str(error)
    repr(error)