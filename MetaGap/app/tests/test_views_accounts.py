"""Tests covering user-facing account and profile views."""

from __future__ import annotations

from unittest.mock import patch
from uuid import uuid4

from django.contrib.auth import get_user_model
from django.contrib.messages import get_messages
from django.test import RequestFactory, TestCase
from django.urls import reverse

from ..forms import CustomUserCreationForm, ImportDataForm, SearchForm
from ..models import SampleGroup
from ..views import EditProfileView, ProfileView


class HomePageViewTests(TestCase):
    """Ensure the landing page renders the expected search form."""

    def test_search_form_present_in_context(self) -> None:
        response = self.client.get(reverse("home"))

        self.assertEqual(response.status_code, 200)
        self.assertIn("form", response.context)
        self.assertIsInstance(response.context["form"], SearchForm)

    def test_search_form_preserves_submitted_filters(self) -> None:
        response = self.client.get(
            reverse("home"),
            {
                "chrom": "7",
                "pos_min": "1000",
                "ref": "C",
                "alt": "G",
                "pass_only": "on",
                "af_min": "0.25",
            },
        )

        self.assertEqual(response.status_code, 200)

        form = response.context["form"]
        self.assertTrue(form.is_bound)
        self.assertEqual(form["chrom"].value(), "7")
        self.assertEqual(form["pos_min"].value(), "1000")
        self.assertEqual(form["ref"].value(), "C")
        self.assertEqual(form["alt"].value(), "G")
        self.assertEqual(form["af_min"].value(), "0.25")
        self.assertEqual(form["pass_only"].value(), True)

        self.assertTrue(form.is_valid())
        self.assertEqual(form.cleaned_data["pos_min"], 1000)
        self.assertAlmostEqual(form.cleaned_data["af_min"], 0.25)
        self.assertTrue(form.cleaned_data["pass_only"])


class UserRegistrationViewTests(TestCase):
    """Validate the sign-up view renders the registration form template."""

    def test_signup_view_renders_template_with_form(self) -> None:
        response = self.client.get(reverse("signup"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "signup.html")
        self.assertIn("form", response.context)
        self.assertIsInstance(response.context["form"], CustomUserCreationForm)
        self.assertContains(response, "<h2>Sign Up</h2>", html=True)

    def test_signup_flow_creates_user_and_profile(self) -> None:
        User = get_user_model()
        unique_suffix = uuid4().hex[:8]
        form_payload = {
            "username": f"testuser_{unique_suffix}",
            "email": f"test{unique_suffix}@example.com",
            "organization_name": "Test Organization",
            "password1": "StrongPass123!",
            "password2": "StrongPass123!",
        }

        response = self.client.post(
            reverse("signup"),
            form_payload,
            follow=True,
        )

        self.assertEqual(response.status_code, 200)
        self.assertEqual(response.redirect_chain, [(reverse("login"), 302)])
        self.assertTemplateUsed(response, "login.html")

        try:
            created_user = User.objects.get(username=form_payload["username"])
            self.assertEqual(created_user.email, form_payload["email"])
            self.assertTrue(hasattr(created_user, "organization_profile"))
            self.assertEqual(
                created_user.organization_profile.organization_name,
                form_payload["organization_name"],
            )
        finally:
            User.objects.filter(username=form_payload["username"]).delete()


class StaticPageViewTests(TestCase):
    """Verify the informational static pages render their expected content."""

    def test_contact_page_renders_expected_information(self) -> None:
        response = self.client.get(reverse("contact"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "contact.html")
        self.assertContains(response, "support@metagap.org")
        self.assertContains(response, "GitHub Discussions")

    def test_about_page_renders_expected_information(self) -> None:
        response = self.client.get(reverse("about"))

        self.assertEqual(response.status_code, 200)
        self.assertTemplateUsed(response, "about.html")
        self.assertContains(response, "About")
        self.assertContains(response, "MetaGaP is an open platform")


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

    def test_profile_redirects_anonymous_users_to_login(self) -> None:
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 302)
        login_url = reverse("login")
        self.assertTrue(response.url.startswith(login_url))

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

    def test_profile_context_uses_mixin_helpers(self) -> None:
        self.client.force_login(self.user)
        sentinel_profile = object()
        queryset = SampleGroup.objects.filter(
            created_by=self.user.organization_profile
        )

        with patch.object(
            ProfileView, "get_organization_profile", return_value=sentinel_profile
        ) as mock_profile, patch.object(
            ProfileView, "get_owned_sample_groups", return_value=queryset
        ) as mock_groups:
            response = self.client.get(reverse("profile"))

        mock_profile.assert_called_once_with()
        mock_groups.assert_called_once_with()
        self.assertIs(response.context["organization_profile"], sentinel_profile)
        self.assertEqual(
            list(response.context["sample_groups"]),
            list(queryset.order_by("name")),
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
