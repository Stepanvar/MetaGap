"""Basic view smoke tests."""

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from app.models import SampleGroup


class ViewSmokeTests(TestCase):
    def test_home_page_renders(self):
        response = self.client.get(reverse("home"))
        self.assertEqual(response.status_code, 200)

    def test_about_page_renders(self):
        response = self.client.get(reverse("about"))
        self.assertEqual(response.status_code, 200)

    def test_contact_page_renders(self):
        response = self.client.get(reverse("contact"))
        self.assertEqual(response.status_code, 200)


class ProfileViewTests(TestCase):
    def setUp(self):
        super().setUp()
        User = get_user_model()
        self.user = User.objects.create_user(
            username="profile_user", password="password123", email="user@example.com"
        )
        self.other_user = User.objects.create_user(
            username="other_user", password="password123", email="other@example.com"
        )

    def test_profile_context_includes_owned_sample_groups(self):
        SampleGroup.objects.create(name="Alpha", created_by=self.user)
        SampleGroup.objects.create(name="Beta", created_by=self.user)
        SampleGroup.objects.create(name="Gamma", created_by=self.other_user)

        self.client.force_login(self.user)
        response = self.client.get(reverse("profile"))

        self.assertEqual(response.status_code, 200)

        sample_groups = list(response.context["sample_groups"])
        self.assertEqual([group.name for group in sample_groups], ["Alpha", "Beta"])

        self.assertEqual(
            response.context["organization_profile"],
            self.user.organization_profile,
        )
