from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path("", include("app.urls")),  # Include app's URL patterns
    path("admin/", admin.site.urls),
    path("admin/doc/", include("django.contrib.admindocs.urls")),  # If using admindocs
]
