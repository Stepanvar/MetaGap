"""Tests for search results and dashboard views."""

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from ..filters import AlleleFrequencySearchFilter
from ..models import (
    AlleleFrequency,
    BioinfoVariantCalling,
    Format,
    Info,
    SampleGroup,
    SampleOrigin,
)


class SearchResultsViewTests(TestCase):
    """Validate search behaviour, filter context, and table population."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="search_user",
            password="search-pass",
            email="search@example.com",
        )

        self.kidney_origin = SampleOrigin.objects.create(
            tissue="Kidney",
            collection_method="Biopsy",
            storage_conditions="Cryogenic",
        )
        self.liver_origin = SampleOrigin.objects.create(
            tissue="Liver",
            collection_method="Surgical",
            storage_conditions="Room Temperature",
        )

        self.variant_caller_a = BioinfoVariantCalling.objects.create(tool="CallerA")
        self.variant_caller_b = BioinfoVariantCalling.objects.create(tool="CallerB")

        self.kidney_group = SampleGroup.objects.create(
            name="Kidney Cohort",
            source_lab="Lab One",
            sample_origin=self.kidney_origin,
            bioinfo_variant_calling=self.variant_caller_a,
            created_by=self.user.organization_profile,
        )
        self.liver_group = SampleGroup.objects.create(
            name="Liver Cohort",
            source_lab="Lab Two",
            sample_origin=self.liver_origin,
            bioinfo_variant_calling=self.variant_caller_b,
            created_by=self.user.organization_profile,
        )

        self.format = Format.objects.create(genotype="0/1")

        self.kidney_info_high = Info.objects.create(
            af="0.15",
            ac="12",
            an="80",
            dp="55",
            mq="60",
        )
        self.kidney_info_low = Info.objects.create(
            af="0.05",
            ac="4",
            an="60",
            dp="35",
            mq="45",
        )
        self.liver_info = Info.objects.create(
            af="0.22",
            ac="18",
            an="90",
            dp="70",
            mq="58",
        )

        self.kidney_variant_pass = AlleleFrequency.objects.create(
            sample_group=self.kidney_group,
            chrom="1",
            pos=150,
            ref="A",
            alt="T",
            qual=180.5,
            filter="PASS",
            info=self.kidney_info_high,
            format=self.format,
        )
        self.kidney_variant_filtered = AlleleFrequency.objects.create(
            sample_group=self.kidney_group,
            chrom="2",
            pos=300,
            ref="C",
            alt="G",
            qual=75.0,
            filter="q10",
            info=self.kidney_info_low,
            format=self.format,
        )
        self.liver_variant = AlleleFrequency.objects.create(
            sample_group=self.liver_group,
            chrom="3",
            pos=400,
            ref="G",
            alt="A",
            qual=95.2,
            filter="PASS",
            info=self.liver_info,
            format=self.format,
        )

    def test_query_filters_variants_and_prioritises_columns(self) -> None:
        response = self.client.get(reverse("search_results"), {"query": "Kidney"})

        self.assertEqual(response.status_code, 200)

        table = response.context["table"]
        records = [row.record for row in table.rows]
        self.assertEqual(records, [self.kidney_variant_pass, self.kidney_variant_filtered])

        filterset = response.context["filter"]
        self.assertIsInstance(filterset, AlleleFrequencySearchFilter)
        self.assertEqual(filterset.form["query"].value(), "Kidney")

        column_names = [column.name for column in table.columns]
        self.assertEqual(
            column_names[:11],
            [
                "chrom",
                "pos",
                "ref",
                "alt",
                "qual",
                "filter",
                "info__af",
                "info__ac",
                "info__an",
                "info__dp",
                "info__mq",
            ],
        )
        self.assertNotIn("info__additional", column_names)
        self.assertNotIn("format__payload", column_names)

    def test_combined_filters_reduce_results(self) -> None:
        response = self.client.get(
            reverse("search_results"),
            {
                "chrom": "1",
                "pos_min": 100,
                "pos_max": 400,
                "pass_only": "True",
                "af_min": 0.1,
                "mq_min": 55,
                "sample_group_source_lab": "Lab One",
                "variant_calling_tool": "CallerA",
            },
        )

        self.assertEqual(response.status_code, 200)
        table = response.context["table"]
        records = [row.record for row in table.rows]
        self.assertEqual(records, [self.kidney_variant_pass])

    def test_pass_only_filter_excludes_non_pass_variants(self) -> None:
        response = self.client.get(
            reverse("search_results"),
            {"pass_only": "True"},
        )

        self.assertEqual(response.status_code, 200)
        table = response.context["table"]
        records = [row.record for row in table.rows]
        self.assertEqual(records, [self.kidney_variant_pass, self.liver_variant])

    def test_search_results_render_full_queryset_for_datatables(self) -> None:
        # Create more variants than the default django-tables2 page size to ensure
        # that pagination is disabled and the full queryset is rendered.
        for index in range(12):
            info = Info.objects.create(af=str(0.2 + index / 100))
            AlleleFrequency.objects.create(
                sample_group=self.kidney_group,
                chrom="5",
                pos=1000 + index,
                ref="A",
                alt="T",
                qual=50 + index,
                filter="PASS",
                info=info,
                format=self.format,
            )

        response = self.client.get(reverse("search_results"))

        self.assertEqual(response.status_code, 200)
        table = response.context["table"]
        self.assertEqual(len(table.rows), AlleleFrequency.objects.count())

    def test_search_results_use_advanced_filters_toggle(self) -> None:
        response = self.client.get(reverse("search_results"))

        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "data-bs-target=\"#advancedFiltersCollapse\"")
        self.assertContains(response, "Advanced filters")
        self.assertNotContains(response, "filter-card-header")


class DashboardViewTests(TestCase):
    """Exercise the dashboard context aggregation for organizations."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user_one = User.objects.create_user(
            username="dashboard_user_one",
            password="dashboard-pass",
            email="one@example.com",
        )
        self.user_two = User.objects.create_user(
            username="dashboard_user_two",
            password="dashboard-pass",
            email="two@example.com",
        )

        self.user_one_groups = [
            SampleGroup.objects.create(
                name=f"User One Group {index}",
                created_by=self.user_one.organization_profile,
            )
            for index in range(8)
        ]

        self.user_two_groups = [
            SampleGroup.objects.create(
                name=f"User Two Group {index}",
                created_by=self.user_two.organization_profile,
            )
            for index in range(2)
        ]

        # Populate allele frequencies for each organization.
        self.user_one_actions = [
            AlleleFrequency.objects.create(
                sample_group=group,
                chrom="1",
                pos=index + 1,
                ref="A",
                alt="T",
            )
            for index, group in enumerate(self.user_one_groups)
        ]

        self.user_two_actions = [
            AlleleFrequency.objects.create(
                sample_group=group,
                chrom="2",
                pos=index + 10,
                ref="C",
                alt="G",
            )
            for index, group in enumerate(self.user_two_groups)
        ]

    def test_dashboard_limits_and_filters_datasets_and_actions(self) -> None:
        self.client.force_login(self.user_one)

        response = self.client.get(reverse("dashboard"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "dashboard.html")

        recent_datasets = response.context["recent_datasets"]
        expected_datasets = list(
            SampleGroup.objects.filter(
                created_by=self.user_one.organization_profile
            ).order_by("-pk")[:6]
        )

        self.assertEqual(recent_datasets, expected_datasets)
        self.assertTrue(all(dataset.created_by == self.user_one.organization_profile for dataset in recent_datasets))
        self.assertLessEqual(len(recent_datasets), 6)

        recent_actions = response.context["recent_actions"]
        expected_actions = list(
            AlleleFrequency.objects.filter(
                sample_group__created_by=self.user_one.organization_profile
            ).order_by("-pk")[:6]
        )

        self.assertEqual(recent_actions, expected_actions)
        self.assertTrue(
            all(
                action.sample_group.created_by == self.user_one.organization_profile
                for action in recent_actions
            )
        )
        self.assertLessEqual(len(recent_actions), 6)

    def test_dashboard_respects_organization_for_different_user(self) -> None:
        self.client.force_login(self.user_two)

        response = self.client.get(reverse("dashboard"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "dashboard.html")

        self.assertEqual(
            response.context["recent_datasets"],
            list(
                SampleGroup.objects.filter(
                    created_by=self.user_two.organization_profile
                ).order_by("-pk")[:6]
            ),
        )
        self.assertEqual(
            response.context["recent_actions"],
            list(
                AlleleFrequency.objects.filter(
                    sample_group__created_by=self.user_two.organization_profile
                ).order_by("-pk")[:6]
            ),
        )

    def test_dashboard_redirects_anonymous_users_to_login(self) -> None:
        response = self.client.get(reverse("dashboard"))

        self.assertEqual(response.status_code, 302)
        login_url = reverse("login")
        self.assertTrue(response.url.startswith(login_url))
