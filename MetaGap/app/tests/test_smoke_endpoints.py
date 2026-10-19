from django.test import Client, TestCase
from django.urls import reverse
from django.contrib.auth import get_user_model


class SmokeEndpointsTest(TestCase):
    """Endpoints the deploy smoke must hit. If these 404 the smoke is wrong."""

    @classmethod
    def setUpTestData(cls):
        User = get_user_model()
        cls.user = User.objects.create_user(username="smoke", password="x")

    def test_homepage_resolves(self):
        resp = self.client.get("/")
        self.assertIn(resp.status_code, (200, 302))

    def test_login_page_resolves(self):
        resp = self.client.get(reverse("login"))
        self.assertEqual(resp.status_code, 200)
