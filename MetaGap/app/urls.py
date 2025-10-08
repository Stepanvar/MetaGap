from django.urls import path
from django.conf import settings
from django.conf.urls.static import static
from django.contrib.auth import views as auth_views
from . import views

urlpatterns = [
    path('', views.HomePageView.as_view(), name='home'),
    path('dashboard/', views.DashboardView.as_view(), name='dashboard'),
    path('search/', views.SearchResultsView.as_view(), name='search_results'),
    path("contact/", views.ContactView.as_view(), name="contact"),
    path("about/", views.AboutView.as_view(), name="about"),
    path("signup/", views.UserRegistrationView.as_view(), name="signup"),
    path("login/", auth_views.LoginView.as_view(template_name="login.html"), name="login"),
    path("logout/", auth_views.LogoutView.as_view(next_page="/"), name="logout"),
    path("profile/", views.ProfileView.as_view(), name="profile"),
    path("profile/edit/", views.EditProfileView.as_view(), name="edit_profile"),
    path("profile/delete/", views.DeleteAccountView.as_view(), name="delete_account"),
    path(
        "profile/sample-groups/<int:pk>/",
        views.SampleGroupDetailView.as_view(),
        name="sample_group_detail",
    ),
    path(
        "profile/sample-groups/<int:pk>/edit/",
        views.SampleGroupUpdateView.as_view(),
        name="sample_group_edit",
    ),
    path("import/", views.ImportDataView.as_view(), name="import_data"),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT) + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)
