"""Unit tests for the ImportDataView helper methods."""

import os
import tempfile

from django.test import RequestFactory, TestCase

from app.models import Format, GenomeComplexity, Info
from app.views import ImportDataView


class ImportHelpersTests(TestCase):
    """Validate mapping logic in ImportDataView helper methods."""

    def setUp(self):
        self.factory = RequestFactory()
        self.request = self.factory.get("/import/")
        self.view = ImportDataView()
        self.view.request = self.request

    def test_create_info_instance_uses_field_map(self):
        """Mapped INFO keys populate model fields and extras go to additional."""

        info_payload = {
            "AC": 5,
            "DP": [10, 11],
            "Custom": "value",
        }

        info_instance = self.view._create_info_instance(info_payload)

        self.assertIsNotNone(info_instance)
        self.assertEqual(Info.objects.count(), 1)
        self.assertEqual(info_instance.ac, "5")
        self.assertEqual(info_instance.dp, "10,11")
        self.assertEqual(info_instance.additional, {"custom": "value"})

    def test_create_format_instance_uses_field_map(self):
        """Mapped FORMAT keys populate structured fields and extras remain JSON."""

        class MockSampleData(dict):
            phased = False

        sample_data = MockSampleData({"GT": (0, 1), "GQ": 99, "Extra": ["foo", "bar"]})
        samples = {"SAMPLE1": sample_data}

        format_instance, sample_name = self.view._create_format_instance(samples)

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
            metadata = self.view._extract_metadata_text_fallback(handle.name)
        finally:
            os.remove(handle.name)

        self.assertEqual(metadata.get("id"), "Sample01")
        self.assertEqual(metadata.get("name"), "Sample01")
        self.assertEqual(metadata.get("description"), "First sample")
        self.assertEqual(metadata.get("project"), "Test")
        self.assertIn("flagwithoutvalue", metadata)
        self.assertNotIn("malformedentry", metadata)
        self.assertNotIn("shouldnotappear", metadata.values())

    def test_extract_section_data_handles_hyphenated_keys(self):
        """Hyphenated section keys map to fields and extras become additional data."""

        hyphen_key = "Genome_Complexity-GC_Content".lower()
        metadata = {
            hyphen_key: "45%",
            "genome_complexity": "3.2 Gbp",
            "genome_complexity-extra_metric": "value",
        }

        section_data, consumed, additional = self.view._extract_section_data(
            metadata, "genome_complexity", GenomeComplexity
        )

        self.assertEqual(section_data.get("gc_content"), "45%")
        self.assertEqual(section_data.get("size"), "3.2 Gbp")
        self.assertIn(hyphen_key, consumed)
        self.assertIn("genome_complexity", consumed)
        self.assertNotIn("genome_complexity-extra_metric", consumed)
        self.assertEqual(additional, {"extra_metric": "value"})
