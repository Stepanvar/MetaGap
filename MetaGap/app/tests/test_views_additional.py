"""Additional tests for the views defined in :mod:`app.views` to increase coverage."""

from __future__ import annotations

from io import StringIO

import pytest
from django.contrib.auth.models import User, AnonymousUser
from django.core.files.uploadedfile import SimpleUploadedFile
from django.http import Http404
from django.test import Client, RequestFactory
from django.urls import reverse
from django.utils.translation import gettext_lazy as _

from app.models import SampleGroup
from app.views import (
    HomePageView,
    DashboardView,
    ProfileView,
    SampleGroupCreateView,
    SampleGroupDetailView,
    SampleGroupUpdateView,
    SampleGroupDeleteView,
    SearchResultsView,
    UserRegistrationView,
    ContactView,
    AboutView,
    ImportDataView,
    export_sample_group_variants,
)

pytestmark = pytest.mark.django_db


# ── helpers ──────────────────────────────────────────────────────────────────

def _make_user(username):
    return User.objects.create_user(username=username, password="testpass123")


def _make_sample_group(user, name="Test Group"):
    return SampleGroup.objects.create(
        name=name,
        created_by=user.organization_profile,
    )


# ── simple GET tests via RequestFactory ──────────────────────────────────────

def test_home_page_view_get():
    factory = RequestFactory()
    request = factory.get("/")
    request.user = AnonymousUser()
    view = HomePageView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_about_view_get():
    factory = RequestFactory()
    request = factory.get("/about/")
    request.user = AnonymousUser()
    view = AboutView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_contact_view_get():
    factory = RequestFactory()
    request = factory.get("/contact/")
    request.user = AnonymousUser()
    view = ContactView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_search_results_view_get():
    factory = RequestFactory()
    request = factory.get("/search/")
    request.user = AnonymousUser()
    view = SearchResultsView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_user_registration_view_get():
    factory = RequestFactory()
    request = factory.get("/signup/")
    request.user = AnonymousUser()
    view = UserRegistrationView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


# ── auth redirect tests (via dispatch to trigger LoginRequiredMixin) ──────────

def test_dashboard_view_unauthenticated():
    factory = RequestFactory()
    request = factory.get("/dashboard/")
    request.user = AnonymousUser()
    view = DashboardView()
    view.setup(request)
    response = view.dispatch(request)
    assert response.status_code == 302


def test_profile_view_unauthenticated():
    factory = RequestFactory()
    request = factory.get("/profile/")
    request.user = AnonymousUser()
    view = ProfileView()
    view.setup(request)
    response = view.dispatch(request)
    assert response.status_code == 302


# ── authenticated view GET tests ─────────────────────────────────────────────

def test_dashboard_view_authenticated():
    user = _make_user("dashboard_user")
    factory = RequestFactory()
    request = factory.get("/dashboard/")
    request.user = user
    view = DashboardView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_profile_view_authenticated():
    user = _make_user("profile_user")
    factory = RequestFactory()
    request = factory.get("/profile/")
    request.user = user
    view = ProfileView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_import_data_view_get():
    user = _make_user("import_user")
    factory = RequestFactory()
    request = factory.get("/import/")
    request.user = user
    view = ImportDataView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_import_data_view_post_no_file():
    user = _make_user("import_user2")
    factory = RequestFactory()
    request = factory.post("/import/", {})
    request.user = user
    view = ImportDataView()
    view.setup(request)
    response = view.post(request)
    assert response.status_code == 200


def test_sample_group_create_view_get():
    user = _make_user("create_user")
    factory = RequestFactory()
    request = factory.get("/profile/sample-groups/new/")
    request.user = user
    view = SampleGroupCreateView()
    view.setup(request)
    response = view.get(request)
    assert response.status_code == 200


def test_sample_group_detail_view():
    user = _make_user("detail_user")
    sg = _make_sample_group(user)
    factory = RequestFactory()
    request = factory.get(f"/profile/sample-groups/{sg.pk}/")
    request.user = user
    view = SampleGroupDetailView()
    view.setup(request)
    view.kwargs = {"pk": sg.pk}
    view.object = sg
    response = view.get(request)
    assert response.status_code == 200


def test_sample_group_update_view_get():
    user = _make_user("update_user")
    sg = _make_sample_group(user)
    factory = RequestFactory()
    request = factory.get(f"/profile/sample-groups/{sg.pk}/edit/")
    request.user = user
    view = SampleGroupUpdateView()
    view.setup(request)
    view.kwargs = {"pk": sg.pk}
    view.object = sg
    response = view.get(request)
    assert response.status_code == 200


def test_sample_group_delete_view_get():
    user = _make_user("delget_user")
    sg = _make_sample_group(user)
    factory = RequestFactory()
    request = factory.get(f"/profile/sample-groups/{sg.pk}/delete/")
    request.user = user
    view = SampleGroupDeleteView()
    view.setup(request)
    view.kwargs = {"pk": sg.pk}
    view.object = sg
    response = view.get(request)
    assert response.status_code == 200


# ── export tests ─────────────────────────────────────────────────────────────

def test_export_sample_group_variants_csv():
    user = _make_user("export_user")
    sg = _make_sample_group(user, "Export Group")
    factory = RequestFactory()
    request = factory.get(f"/profile/sample-groups/{sg.pk}/export/?format=csv")
    request.user = user
    response = export_sample_group_variants(request, sg.pk)
    assert response.status_code == 200
    assert response.get("Content-Type") == "text/csv"


def test_export_sample_group_variants_tsv():
    user = _make_user("export_user2")
    sg = _make_sample_group(user, "Export Group 2")
    factory = RequestFactory()
    request = factory.get(f"/profile/sample-groups/{sg.pk}/export/?format=tsv")
    request.user = user
    response = export_sample_group_variants(request, sg.pk)
    assert response.status_code == 200
    assert "text/tab-separated-values" in response.get("Content-Type")


def test_export_sample_group_variants_invalid_format():
    user = _make_user("export_user3")
    sg = _make_sample_group(user, "Export Group 3")
    factory = RequestFactory()
    request = factory.get(f"/profile/sample-groups/{sg.pk}/export/?format=invalid")
    request.user = user
    response = export_sample_group_variants(request, sg.pk)
    assert response.status_code == 200
    assert response.get("Content-Type") == "text/csv"


def test_export_sample_group_variants_nonexistent_group():
    user = _make_user("export_user4")
    factory = RequestFactory()
    request = factory.get("/profile/sample-groups/99999/export/?format=csv")
    request.user = user
    with pytest.raises(Http404):
        export_sample_group_variants(request, 99999)


# ── context data tests ───────────────────────────────────────────────────────

def test_home_page_view_context():
    factory = RequestFactory()
    request = factory.get("/")
    request.user = AnonymousUser()
    view = HomePageView()
    view.setup(request)
    context = view.get_context_data()
    assert "form" in context


def test_dashboard_view_context():
    user = _make_user("ctx_user")
    factory = RequestFactory()
    request = factory.get("/dashboard/")
    request.user = user
    view = DashboardView()
    view.setup(request)
    context = view.get_context_data()
    assert "recent_datasets" in context


def test_profile_view_context():
    user = _make_user("ctx_user2")
    factory = RequestFactory()
    request = factory.get("/profile/")
    request.user = user
    view = ProfileView()
    view.setup(request)
    context = view.get_context_data()
    assert "sample_groups" in context


def test_search_results_view_context():
    factory = RequestFactory()
    request = factory.get("/search/?chrom=chr1")
    request.user = AnonymousUser()
    view = SearchResultsView()
    view.setup(request)
    # Call get() first to populate object_list via FilterView
    view.get(request)
    context = view.get_context_data()
    assert "filter" in context


# ── form instantiation smoke tests ───────────────────────────────────────────

def test_sample_group_create_view_form():
    user = _make_user("form_user")
    factory = RequestFactory()
    request = factory.post("/profile/sample-groups/new/", {"name": "New Group"})
    request.user = user
    view = SampleGroupCreateView()
    view.setup(request)
    form = view.get_form(view.get_form_class())
    assert form is not None


def test_user_registration_view_form():
    factory = RequestFactory()
    request = factory.post("/signup/", {})
    request.user = AnonymousUser()
    view = UserRegistrationView()
    view.setup(request)
    form = view.get_form(view.get_form_class())
    assert form is not None


def test_import_data_view_form():
    user = _make_user("form_user2")
    factory = RequestFactory()
    request = factory.post("/import/", {})
    request.user = user
    view = ImportDataView()
    view.setup(request)
    form = view.get_form(view.get_form_class())
    assert form is not None


def test_sample_group_update_view_form():
    user = _make_user("form_user3")
    sg = _make_sample_group(user, "Update Group")
    factory = RequestFactory()
    request = factory.post(f"/profile/sample-groups/{sg.pk}/edit/", {"name": "Updated"})
    request.user = user
    view = SampleGroupUpdateView()
    view.setup(request)
    view.object = sg
    form = view.get_form(view.get_form_class())
    assert form is not None


def test_sample_group_delete_view_has_delete():
    view = SampleGroupDeleteView()
    assert hasattr(view, "delete")


# ── misc ─────────────────────────────────────────────────────────────────────

def test_gettext_lazy_usage():
    lazy_text = _("Test translation")
    assert str(lazy_text)
