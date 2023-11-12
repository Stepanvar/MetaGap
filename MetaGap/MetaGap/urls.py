"""
Definition of urls for MetaGap.
"""

from django.urls import path
from django.contrib import admin
from app import views


urlpatterns = [
    path("", views.home, name="home"),
    path("results/", views.results, name="results"),
    path("contact/", views.contact, name="contact"),
    path("about/", views.about, name="about"),
    path("signup/", views.signup, name="signup"),
    path("login/", views.login, name="login"),
    path("logout/", views.logout, name="logout"),
    path("admin/", admin.site.urls),
    path("profile/", views.profile, name="profile"),
    path("adddata/", views.adddata, name="adddata"),
]
