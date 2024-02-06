from django.contrib import admin
from django.urls import path
from django.contrib.auth import views as auth_views
from app import views

urlpatterns = [
	path("", views.HomePageView.as_view(), name="home"),
	path("results/", views.SearchResultView.as_view(), name="results"),
	path(
		"contact/", views.ContactView.as_view(), name="contact"
	),
	path(
		"about/", views.AboutView.as_view(), name="about"
	),
	path(
		"signup/", views.UserRegistrationView.as_view(), name="signup"
	),
	path(
		"login/", auth_views.LoginView.as_view(template_name="login.html"), name="login"
	),
	path(
		"logout/", auth_views.LogoutView.as_view(next_page="/"), name="logout"
	),
	path("admin/", admin.site.urls),
	path("profile/", views.ProfileView.as_view(), name="profile"),
	path("adddata/", views.AddDataView.as_view(), name="adddata"),
]
