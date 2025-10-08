"""Focused regression tests for the primary application views."""

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import RequestFactory, TestCase
from django.urls import reverse

from ..filters import SampleGroupFilter
from ..forms import ImportDataForm, SearchForm
from ..models import AlleleFrequency, SampleGroup, SampleOrigin
from ..views import EditProfileView


class HomePageViewTests(TestCase):
    """Ensure the landing page renders the expected search form."""

    def test_search_form_present_in_context(self) -> None:
        response = self.client.get(reverse("home"))

        self.assertEqual(response.status_code, 200)
        self.assertIn("form", response.context)
        self.assertIsInstance(response.context["form"], SearchForm)


class ProfileViewTests(TestCase):
    """Exercise the profile dashboard and its supporting context."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="profile_user",
            password="password123",
            email="user@example.com",
        )
        self.other_user = User.objects.create_user(
            username="other_user",
            password="password123",
            email="other@example.com",
        )

    def test_profile_lists_user_groups_and_import_form(self) -> None:
        SampleGroup.objects.create(
            name="VCF Cohort",
            created_by=self.user.organization_profile,
        )

        self.client.login(username="profile_user", password="password123")
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.context["organization_profile"], self.user.organization_profile)
        self.assertIsInstance(response.context["import_form"], ImportDataForm)
        self.assertEqual(response.context["import_form_action"], reverse("import_data"))
        self.assertEqual(response.context["import_form_enctype"], "multipart/form-data")

    def test_profile_context_includes_only_owned_sample_groups(self) -> None:
        SampleGroup.objects.create(
            name="Alpha",
            created_by=self.user.organization_profile,
        )
        SampleGroup.objects.create(
            name="Beta",
            created_by=self.user.organization_profile,
        )
        SampleGroup.objects.create(
            name="Gamma",
            created_by=self.other_user.organization_profile,
        )

        self.client.force_login(self.user)
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)

        expected_queryset = SampleGroup.objects.filter(
            created_by=self.user.organization_profile
        ).order_by("name")
        self.assertQuerySetEqual(
            response.context["sample_groups"],
            SampleGroup.objects.filter(
                created_by=self.user.organization_profile
            ).order_by("name"),
            transform=lambda group: group,
        )
        self.assertIsInstance(response.context["import_form"], ImportDataForm)
        self.assertEqual(response.context["import_form_action"], reverse("import_data"))
        self.assertEqual(
            response.context["import_form_enctype"],
            "multipart/form-data",
        )


class EditProfileViewTests(TestCase):
    """Validate the edit profile workflow and form configuration."""

    def setUp(self) -> None:
        super().setUp()
        self.factory = RequestFactory()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="edit_user",
            password="test-pass-123",
            email="initial@example.com",
            first_name="Initial",
            last_name="User",
        )

    def test_form_kwargs_uses_authenticated_user_instance(self) -> None:
        request = self.factory.get(reverse("edit_profile"))
        request.user = self.user

        view = EditProfileView()
        view.setup(request)

        form_kwargs = view.get_form_kwargs()

        self.assertIs(form_kwargs["instance"], self.user)

    def test_successful_profile_update_redirects_and_persists(self) -> None:
        self.client.force_login(self.user)

        response = self.client.post(
            reverse("edit_profile"),
            {
                "username": "edit_user",
                "email": "updated@example.com",
                "first_name": "Updated",
                "last_name": "Name",
                "organization_name": "Updated Org",
            },
        )

        self.assertRedirects(response, reverse("profile"))

        self.user.refresh_from_db()
        self.assertEqual(self.user.email, "updated@example.com")
        self.assertEqual(self.user.first_name, "Updated")
        self.assertEqual(self.user.last_name, "Name")
        self.assertEqual(
            self.user.organization_profile.organization_name,
            "Updated Org",
        )

        messages = list(get_messages(response.wsgi_request))
        self.assertIn(
            "Your profile has been updated.",
            [message.message for message in messages],
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

        self.kidney_group = SampleGroup.objects.create(
            name="Kidney Cohort",
            sample_origin=self.kidney_origin,
            created_by=self.user.organization_profile,
        )
        self.liver_group = SampleGroup.objects.create(
            name="Liver Cohort",
            sample_origin=self.liver_origin,
            created_by=self.user.organization_profile,
        )

    def test_search_query_filters_expected_sample_group(self) -> None:
        response = self.client.get(reverse("search_results"), {"query": "Kidney"})

        self.assertEqual(response.status_code, 200)

        table = response.context["table"]
        self.assertEqual([row.record for row in table.rows], [self.kidney_group])

        sample_filter = response.context["filter"]
        self.assertIsInstance(sample_filter, SampleGroupFilter)
        self.assertEqual(sample_filter.data.get("query"), "Kidney")
        self.assertQuerySetEqual(
            sample_filter.qs,
            [self.kidney_group],
            transform=lambda group: group,
        )

        form = response.context["form"]
        self.assertIsInstance(form, SearchForm)
        self.assertEqual(form.data.get("query"), "Kidney")

    def test_empty_query_returns_all_records(self) -> None:
        response = self.client.get(reverse("search_results"), {"query": ""})

        self.assertEqual(response.status_code, 200)

        table = response.context["table"]
        self.assertCountEqual(
            [row.record for row in table.rows],
            [self.kidney_group, self.liver_group],
        )

        sample_filter = response.context["filter"]
        self.assertIsInstance(sample_filter, SampleGroupFilter)
        self.assertEqual(sample_filter.data.get("query"), "")
        self.assertQuerySetEqual(
            sample_filter.qs.order_by("pk"),
            SampleGroup.objects.order_by("pk"),
            transform=lambda group: group,
        )

        form = response.context["form"]
        self.assertIsInstance(form, SearchForm)
        self.assertEqual(form.data.get("query"), "")


class DashboardViewTests(TestCase):
    """Confirm the dashboard only surfaces recent data from the logged-in user."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()

        self.user_one = User.objects.create_user(
            username="dashboard_user_one",
            password="dashboard-pass-one",
            email="one@example.com",
        )
        self.user_two = User.objects.create_user(
            username="dashboard_user_two",
            password="dashboard-pass-two",
            email="two@example.com",
        )

        # Create more than six datasets for the first user to exercise the cap.
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

class ImportDataViewTests(TestCase):
    """Validate the VCF import workflow end to end."""

    VCF_CONTENT = """##fileformat=VCFv4.2
##contig=<ID=1>
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=CLNSIG,Number=1,Type=String,Description="Clinical significance">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##SAMPLE=<ID=GroupA,Description=Imported group>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample001
1\t1234\trsTest\tA\tT\t99\tPASS\tAF=0.5;CLNSIG=Pathogenic\tGT:GQ\t0/1:99
"""

    def setUp(self) -> None:
        super().setUp()
        self.user = get_user_model().objects.create_user(
            username="vcf_user",
            email="vcf@example.com",
            password="import-pass",
        )

    def test_import_creates_sample_group_and_variant(self) -> None:
        self.client.login(username="vcf_user", password="import-pass")

        upload = SimpleUploadedFile(
            "import.vcf",
            self.VCF_CONTENT.encode("utf-8"),
            content_type="text/vcf",
        )

        response = self.client.post(reverse("import_data"), {"data_file": upload})

        self.assertRedirects(response, reverse("profile"))
        sample_group = SampleGroup.objects.get(name="GroupA")
        self.assertEqual(sample_group.allele_frequencies.count(), 1)
        allele = sample_group.allele_frequencies.get()
        self.assertEqual(allele.chrom, "1")
        self.assertEqual(allele.pos, 1234)
        self.assertEqual(allele.variant_id, "rsTest")
        self.assertEqual(allele.info.af, "0.5")
        self.assertEqual(allele.info.additional["clnsig"], "Pathogenic")
        self.assertEqual(allele.format.genotype, "0/1")
        self.assertEqual(allele.format.fields["gq"], "99")
        self.assertEqual(allele.format.additional["sample_id"], "Sample001")

        self.assertGreaterEqual(AlleleFrequency.objects.count(), 1)


class DeleteAccountViewTests(TestCase):
    """Exercise the account deletion endpoint."""

    def test_confirmed_post_deletes_user_and_logs_out(self) -> None:
        User = get_user_model()
        user = User.objects.create_user(
            username="delete_me",
            password="secure-pass",
            email="deleteme@example.com",
        )

        self.client.force_login(user)

        response = self.client.post(
            reverse("delete_account"),
            {"confirm": True},
            follow=True,
        )

        self.assertRedirects(response, reverse("home"))
        self.assertFalse(User.objects.filter(pk=user.pk).exists())
        self.assertNotIn("_auth_user_id", self.client.session)

        messages = list(get_messages(response.wsgi_request))
        self.assertIn(
            "Your account has been deleted.",
            [message.message for message in messages],
        )


class SampleGroupUpdateViewTests(TestCase):
    """Ensure sample group metadata can be edited securely."""

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="editor",
            email="editor@example.com",
            password="pass12345",
        )
        self.other_user = User.objects.create_user(
            username="intruder",
            email="intruder@example.com",
            password="pass12345",
        )
        self.sample_group = SampleGroup.objects.create(
            name="Editable Group",
            created_by=self.user.organization_profile,
            contact_email="old@example.com",
        )

    def test_login_required(self) -> None:
        response = self.client.get(
            reverse("sample_group_edit", args=[self.sample_group.pk])
        )

        self.assertEqual(response.status_code, 302)
        self.assertIn(reverse("login"), response.url)

    def test_user_can_update_owned_group(self) -> None:
        self.client.force_login(self.user)

        response = self.client.post(
            reverse("sample_group_edit", args=[self.sample_group.pk]),
            {
                "name": "Updated Group",
                "contact_email": "new@example.com",
                "total_samples": 42,
            },
            follow=True,
        )

        self.assertRedirects(response, reverse("profile"))
        self.sample_group.refresh_from_db()
        self.assertEqual(self.sample_group.name, "Updated Group")
        self.assertEqual(self.sample_group.contact_email, "new@example.com")
        self.assertEqual(self.sample_group.total_samples, 42)

    def test_user_cannot_edit_other_organisations_group(self) -> None:
        other_group = SampleGroup.objects.create(
            name="Locked Group",
            created_by=self.other_user.organization_profile,
        )

        self.client.force_login(self.user)
        response = self.client.get(reverse("sample_group_edit", args=[other_group.pk]))

        self.assertEqual(response.status_code, 404)
