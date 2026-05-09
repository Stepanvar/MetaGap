"""Tests for the VCF import form view."""

from __future__ import annotations

import pytest
from django.contrib.auth.models import User
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import Client


pytestmark = pytest.mark.django_db


def _make_user(username):
    return User.objects.create_user(username=username, password="testpass123")


def test_import_view_get_authenticated():
    user = _make_user("import_get_user")
    client = Client()
    client.force_login(user)
    response = client.get("/en/import/")
    assert response.status_code == 200


def test_import_view_get_unauthenticated():
    client = Client()
    response = client.get("/en/import/")
    assert response.status_code == 302
    assert "login" in response["Location"]


def test_import_view_post_no_file():
    user = _make_user("import_nofile_user")
    client = Client()
    client.force_login(user)
    response = client.post("/en/import/", {})
    assert response.status_code == 200
    form = response.context["form"]
    assert form.errors


def test_import_view_post_invalid_file():
    user = _make_user("import_bad_user")
    client = Client()
    client.force_login(user)
    bad_file = SimpleUploadedFile("bad.vcf", b"not a vcf", content_type="text/plain")
    response = client.post("/en/import/", {"data_file": bad_file}, follow=True)
    assert response.status_code == 200
    messages = list(response.context.get("messages", []))
    assert any("error" in str(m.tags) for m in messages)
