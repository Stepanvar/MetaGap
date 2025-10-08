from django.contrib.auth import get_user_model
from django.test import TestCase

from app.models import OrganizationProfile


class UserProfileSignalTests(TestCase):
    def test_organization_profile_created_for_new_user(self):
        user = get_user_model().objects.create_user(
            username="signaltestuser",
            password="testpassword123",
        )

        self.assertTrue(
            OrganizationProfile.objects.filter(user=user).exists(),
            "OrganizationProfile should be created automatically for new users.",
        )

        user.refresh_from_db()
        self.assertTrue(
            hasattr(user, "organization_profile"),
            "User instance should have an organization_profile attribute.",
        )
        self.assertEqual(user.organization_profile.user, user)
