"""Unit tests for the ``VCFImporter`` helper methods."""

import os
import tempfile
from types import SimpleNamespace

from django.contrib.auth import get_user_model
from django.test import TestCase

from app.models import (
    BioinfoAlignment,
    BioinfoVariantCalling,
    Format,
    GenomeComplexity,
    Info,
    SampleGroup,
)
from app.services.vcf_importer import VCFImporter


class ImportHelpersTests(TestCase):
    """Validate mapping logic in the importer helper methods."""

    def setUp(self):
        user_model = get_user_model()
        self.user = user_model.objects.create_user(
            "import-helper", password="test-pass"
        )
        self.importer = VCFImporter(self.user)

    def test_create_info_instance_uses_field_map(self):
        """Mapped INFO keys populate model fields and extras go to additional."""

        info_payload = {
            "AC": 5,
            "DP": [10, 11],
            "QD": 12.34,
            "FS": "5.6",
            "SOR": 1.2,
            "Custom": "value",
        }

        info_instance = self.importer._create_info_instance(info_payload)

        self.assertIsNotNone(info_instance)
        self.assertEqual(Info.objects.count(), 1)
        self.assertEqual(info_instance.ac, "5")
        self.assertEqual(info_instance.dp, "10,11")
        self.assertEqual(info_instance.qd, "12.34")
        self.assertEqual(info_instance.fs, "5.6")
        self.assertEqual(info_instance.sor, "1.2")
        self.assertEqual(info_instance.additional, {"custom": "value"})

    def test_create_info_instance_discards_placeholder_values(self):
        """Placeholder INFO values are stored as NULL instead of literal dots."""

        info_payload = {"AF": ".", "DP": ".", "QD": ".", "FS": "."}

        info_instance = self.importer._create_info_instance(info_payload)

        self.assertIsNotNone(info_instance)
        self.assertEqual(Info.objects.count(), 1)
        self.assertIsNone(info_instance.af)
        self.assertIsNone(info_instance.dp)
        self.assertIsNone(info_instance.additional)

    def test_create_format_instance_uses_field_map(self):
        """Mapped FORMAT keys populate structured fields and extras remain JSON."""

        class MockSampleData(dict):
            phased = False

        sample_data = MockSampleData({"GT": (0, 1), "GQ": 99, "Extra": ["foo", "bar"]})
        samples = {"SAMPLE1": sample_data}

        format_instance, sample_name = self.importer._create_format_instance(samples)

        self.assertEqual(sample_name, "SAMPLE1")
        self.assertIsNotNone(format_instance)
        self.assertEqual(Format.objects.count(), 1)
        self.assertEqual(format_instance.genotype, "0/1")
        self.assertEqual(format_instance.payload["fields"], {"gt": "0/1", "gq": "99"})
        self.assertEqual(format_instance.payload["additional"], {"extra": "foo,bar"})

    def test_extract_metadata_text_fallback_reads_sample_header(self):
        """Metadata extraction lowers keys, derives name, and stops at #CHROM."""

        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".vcf", encoding="utf-8"
        ) as handle:
            handle.write("##fileformat=VCFv4.3\n")
            handle.write(
                "##SAMPLE=<ID=Sample01,Description=First sample,Project=Test,"
                "FlagWithoutValue=,MalformedEntry>\n"
            )
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            handle.write("1\t1000\t.\tA\tT\t.\tPASS\t.\n")
            handle.write("##SAMPLE=<ID=Ignored,Name=ShouldNotAppear>\n")

        try:
            metadata = self.importer._extract_metadata_text_fallback(handle.name)
        finally:
            os.remove(handle.name)

        self.assertEqual(metadata.get("id"), "Sample01")
        self.assertEqual(metadata.get("name"), "Sample01")
        self.assertEqual(metadata.get("description"), "First sample")
        self.assertEqual(metadata.get("project"), "Test")
        self.assertIn("flagwithoutvalue", metadata)
        self.assertNotIn("malformedentry", metadata)
        self.assertNotIn("shouldnotappear", metadata.values())

    def test_extract_section_data_supports_hyphenated_keys(self):
        """Hyphenated section prefixes map fields and feed remaining into additional."""

        raw_key = "Genome_Complexity-GC_Content"
        metadata = {
            raw_key.lower(): "41%",
            "genome_complexity-extra_metric": "observed",
            "unrelated": "ignore",
        }

        section_data, consumed, additional = self.importer._extract_section_data(
            metadata, "genome_complexity", GenomeComplexity
        )

        self.assertEqual(section_data["gc_content"], "41%")
        self.assertIn(raw_key.lower(), consumed)
        self.assertIn("genome_complexity-extra_metric", consumed)
        self.assertEqual(additional, {"extra_metric": "observed"})

    def test_extract_section_data_prefers_section_specific_tool_keys(self):
        """Section-specific aliases win over bare tool keys for overlapping sections."""

        metadata = {
            "bioinfoalignment_software": "BWA",
            "tool": "GATK",
            "bioinfoalignment_params": "-M",
            "bioinfovariantcalling_tool": "GATK",
        }

        alignment_data, alignment_consumed, _ = self.importer._extract_section_data(
            metadata, "bioinfo_alignment", BioinfoAlignment
        )
        variant_data, variant_consumed, _ = self.importer._extract_section_data(
            metadata, "bioinfo_variant_calling", BioinfoVariantCalling
        )

        self.assertEqual(alignment_data["tool"], "BWA")
        self.assertEqual(variant_data["tool"], "GATK")
        self.assertIn("bioinfoalignment_software", alignment_consumed)
        self.assertIn("bioinfovariantcalling_tool", variant_consumed)

    def test_extract_section_data_falls_back_to_section_value(self):
        """Section-level metadata populates the configured primary field when needed."""

        metadata = {"genome_complexity": "3.1Gb", "other": "value"}

        section_data, consumed, additional = self.importer._extract_section_data(
            metadata, "genome_complexity", GenomeComplexity
        )

        self.assertEqual(section_data["size"], "3.1Gb")
        self.assertIn("genome_complexity", consumed)
        self.assertIsNone(additional)

    def test_create_sample_group_recognizes_normalized_sample_headers(self):
        """Sample group creation handles punctuation and camel-case header keys."""

        class MockSampleRecord:
            def __init__(self, data):
                self.key = "SAMPLE"
                self._data = data

            def get(self, key):
                return self._data.get(key)

            def items(self):
                return self._data.items()

        mock_record = MockSampleRecord(
            {
                "ID": "GenomeCohort",
                "Center": "Genome Center",
                "Email?": "contact@example.com",
                "Phone?": "555-0100",
                "N": "8",
                "Inclusion?": "Adults",
            }
        )
        mock_vcf = SimpleNamespace(header=SimpleNamespace(records=[mock_record]))

        metadata = self.importer.extract_sample_group_metadata(mock_vcf)

        self.assertIn("source_lab", metadata)
        self.assertIn("contact_email", metadata)
        self.assertIn("contact_phone", metadata)
        self.assertIn("total_samples", metadata)
        self.assertIn("inclusion_criteria", metadata)
        self.assertNotIn("email?", metadata)

        user_model = get_user_model()
        user = user_model.objects.create_user("importer", password="test-pass")
        profile = user.organization_profile

        sample_group = self.importer._create_sample_group(
            metadata, "mock_file.vcf", profile
        )

        self.assertEqual(SampleGroup.objects.count(), 1)
        self.assertEqual(sample_group.name, "GenomeCohort")
        self.assertEqual(sample_group.source_lab, "Genome Center")
        self.assertEqual(sample_group.contact_email, "contact@example.com")
        self.assertEqual(sample_group.contact_phone, "555-0100")
        self.assertEqual(sample_group.total_samples, 8)
        self.assertEqual(sample_group.inclusion_criteria, "Adults")
