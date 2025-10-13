from django.conf.urls.i18n import i18n_patterns
from django.contrib import admin
from django.urls import path, include

urlpatterns = [
    path("i18n/", include("django.conf.urls.i18n")),
]

urlpatterns += i18n_patterns(
    path("", include("app.urls")),  # Include app's URL patterns
    path("admin/", admin.site.urls),
    path("admin/doc/", include("django.contrib.admindocs.urls")),  # If using admindocs
)
