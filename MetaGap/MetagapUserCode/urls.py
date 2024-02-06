from django.urls import path
from django.contrib.auth import views as auth_views
from app import views
from django.contrib import admin

urlpatterns = [
    path("", views.HomePageView.as_view(), name="home"),
    path("results/", views.SearchResultView.as_view(), name="results"),
    path("contact/", views.contact, name="contact"),
    path("about/", views.about, name="about"),
    path("signup/", views.signup, name="signup"),
    path(
        "login/", auth_views.LoginView.as_view(template_name="login.html"), name="login"
    ),
    path("logout/", auth_views.LogoutView.as_view(next_page="home"), name="logout"),
    path("admin/", admin.site.urls),
    path("profile/", views.ProfileView.as_view(), name="profile"),
    path("adddata/", views.AddDataView.as_view(), name="adddata"),
]

# Ensure to update 'views.signup', 'views.contact', and 'views.about' to class-based views or function-based views accordingly.
# Also, make sure 'login.html' is the correct template path for your login view.
