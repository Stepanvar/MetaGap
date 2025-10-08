"""Unit tests for the ImportDataView helper methods."""

from django.test import RequestFactory, TestCase

from app.models import Format, Info
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

    def test_build_additional_payload_coerces_and_filters_values(self):
        """Additional payload converts values and excludes consumed keys."""

        metadata = {
            "sample_json": '{"coverage": 99}',
            "sample_number": " 42 ",
            "sample_blank": "   ",
            "sample_raw": "hello",
            "sample_float": "3.14",
            "sample_consumed": "ignored",
            "sample.additional": '{"quality": "high"}',
            "other_section": "should be skipped",
        }

        additional = self.view._build_additional_payload(
            metadata,
            section="sample",
            consumed={"sample_consumed"},
        )

        self.assertIsNotNone(additional)
        self.assertEqual(additional["json"], {"coverage": 99})
        self.assertEqual(additional["number"], 42)
        self.assertIsNone(additional["blank"])
        self.assertEqual(additional["raw"], "hello")
        self.assertAlmostEqual(additional["float"], 3.14)
        self.assertEqual(additional["additional"], {"quality": "high"})
        self.assertNotIn("consumed", additional)
        self.assertNotIn("other_section", additional)

    def test_coerce_additional_value_behavior_matrix(self):
        """Direct coercion helper tests guard against regression."""

        coerce = ImportDataView._coerce_additional_value

        self.assertIsNone(coerce("   "))
        self.assertEqual(coerce('{"a": 1}'), {"a": 1})
        self.assertEqual(coerce("42"), 42)
        self.assertAlmostEqual(coerce("3.14"), 3.14)
        self.assertTrue(coerce("true"))
        self.assertEqual(coerce("not-json-or-number"), "not-json-or-number")
