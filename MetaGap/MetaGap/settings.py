"""
Django settings for MetaGap project.

Based on 'django-admin startproject' using Django 2.1.2.

For more information on this file, see
https://docs.djangoproject.com/en/2.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/2.1/ref/settings/
"""

import os, dj_database_url


def _split_env_list(name: str, default: str) -> list[str]:
    """Return a clean list from a comma-separated environment variable."""

    return [value.strip() for value in os.getenv(name, default).split(",") if value.strip()]

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/2.1/howto/deployment/checklist/

DEBUG = os.getenv("DEBUG", "1") == "1"
SECRET_KEY = os.getenv("SECRET_KEY", "dev-secret")
ALLOWED_HOSTS = _split_env_list(
    "ALLOWED_HOSTS", "localhost,127.0.0.1,0.0.0.0,*.github.dev"
)
CSRF_TRUSTED_ORIGINS = _split_env_list(
    "CSRF_TRUSTED_ORIGINS",
    "https://*.github.dev,http://0.0.0.0:8000,http://localhost:8000",
)

SITE_NAME = os.getenv("SITE_NAME", "MetaGaP")
SITE_COLORS = {
    "primary": os.getenv("SITE_COLOR_PRIMARY", "#1f6feb"),
    "secondary": os.getenv("SITE_COLOR_SECONDARY", "#2ea043"),
}
DATABASES = {
    "default": dj_database_url.config(
        default=os.getenv("DATABASE_URL","sqlite:///db.sqlite3"),
        conn_max_age=600,
    )
}
# Application references
# https://docs.djangoproject.com/en/2.1/ref/settings/#std:setting-INSTALLED_APPS
INSTALLED_APPS = [
    # Django's built-in apps
    "django.contrib.admin",
    "django.contrib.admindocs",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",

    # Third-party apps
    "django_bootstrap5",
    "django_filters",
    "django_tables2",

    # Your apps
    "app",
]

# Middleware framework
# https://docs.djangoproject.com/en/2.1/topics/http/middleware/
MIDDLEWARE = [
	"django.middleware.security.SecurityMiddleware",
	"django.contrib.sessions.middleware.SessionMiddleware",
    "whitenoise.middleware.WhiteNoiseMiddleware",
	"django.middleware.common.CommonMiddleware",
	"django.middleware.csrf.CsrfViewMiddleware",
	"django.contrib.auth.middleware.AuthenticationMiddleware",
	"django.contrib.messages.middleware.MessageMiddleware",
	"django.middleware.clickjacking.XFrameOptionsMiddleware",

]

ROOT_URLCONF = "MetaGap.urls"

# Template configuration
# https://docs.djangoproject.com/en/2.1/topics/templates/
TEMPLATES = [
	{
		"BACKEND": "django.template.backends.django.DjangoTemplates",
		"DIRS": [],
		"APP_DIRS": True,
		"OPTIONS": {
                        "context_processors": [
                                "django.template.context_processors.debug",
                                "django.template.context_processors.request",
                                "django.contrib.auth.context_processors.auth",
                                "django.contrib.messages.context_processors.messages",
                                "MetaGap.context_processors.branding",
                        ],
		},
	},
]

WSGI_APPLICATION = "MetaGap.wsgi.application"

DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

# Password validation
# https://docs.djangoproject.com/en/2.1/ref/settings/#auth-password-validators
AUTH_PASSWORD_VALIDATORS = [
	{
		"NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
	},
	{
		"NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
	},
	{
		"NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
	},
	{
		"NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
	},
]

# Internationalization
# https://docs.djangoproject.com/en/2.1/topics/i18n/
LANGUAGE_CODE = "en-us"
TIME_ZONE = "UTC"
USE_I18N = True
USE_L10N = True
USE_TZ = True

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.1/howto/static-files/
STATIC_URL = "/static/"
STATIC_ROOT = os.path.join(BASE_DIR, "staticfiles")
STATICFILES_DIRS = [os.path.join(BASE_DIR, "app", "static")]
# Media files
MEDIA_URL = '/media/'
MEDIA_ROOT = os.path.join(BASE_DIR, 'media')
STORAGES = {
    "staticfiles": {
        "BACKEND": "whitenoise.storage.CompressedManifestStaticFilesStorage"
    }
}
LOGIN_URL = "login"
LOGIN_REDIRECT_URL = 'profile'
LOGOUT_REDIRECT_URL = 'home'
