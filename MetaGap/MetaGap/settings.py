"""
Django settings for MetaGap project.

Based on 'django-admin startproject' using Django 2.1.2.

For more information on this file, see
https://docs.djangoproject.com/en/2.1/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/2.1/ref/settings/
"""

import os, dj_database_url
from pathlib import Path

from django.utils.translation import gettext_lazy as _


def _split_env_list(name: str, default: str) -> list[str]:
    """Return a clean list from a comma-separated environment variable."""

    return [value.strip() for value in os.getenv(name, default).split(",") if value.strip()]


def _env_bool(name: str, default: str = "0") -> bool:
    """Parse a boolean environment flag."""

    return os.getenv(name, default).lower() in {"1", "true", "yes", "on"}

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
BASE_DIR = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

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
    "whitenoise.middleware.WhiteNoiseMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.locale.LocaleMiddleware",
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
LANGUAGE_CODE = "en"
LANGUAGES = [
    ("en", _("English")),
    ("ru", _("Russian")),
]
TIME_ZONE = "UTC"
USE_I18N = True
USE_L10N = True
USE_TZ = True
LOCALE_PATHS = [BASE_DIR / "locale"]
LANGUAGE_COOKIE_NAME = "metagap_language"

# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/2.1/howto/static-files/
STATIC_URL = "/static/"
STATIC_ROOT = os.path.join(BASE_DIR, "staticfiles")
# Media files
MEDIA_URL = '/media/'
MEDIA_ROOT = os.path.join(BASE_DIR, 'media')
STORAGES = {
    "default": {
        "BACKEND": "django.core.files.storage.FileSystemStorage",
        "OPTIONS": {"location": MEDIA_ROOT},
    },
    "staticfiles": {
        "BACKEND": "whitenoise.storage.CompressedManifestStaticFilesStorage"
    },
}
LOGIN_URL = "login"
LOGIN_REDIRECT_URL = 'profile'
LOGOUT_REDIRECT_URL = 'home'
WHITENOISE_MANIFEST_STRICT = False

# Logging
# ------
# When ``DEBUG`` is ``False`` Django will stop rendering its detailed error
# pages and users will only see the generic 500 page.  We still want the full
# traceback to appear in the application logs so developers can investigate
# issues without enabling debug mode in production.  The configuration below
# sends structured log messages to stdout and makes sure request exceptions are
# always recorded together with their traceback information.
LOG_LEVEL = os.getenv("LOG_LEVEL", "INFO").upper()
REQUEST_LOG_LEVEL = os.getenv("DJANGO_REQUEST_LOG_LEVEL", "ERROR").upper()
LOG_DIR = Path(os.getenv("LOG_DIR", BASE_DIR / "logs"))
LOG_DIR.mkdir(parents=True, exist_ok=True)
LOG_MAX_BYTES = int(os.getenv("LOG_MAX_BYTES", 5 * 1024 * 1024))
LOG_BACKUP_COUNT = int(os.getenv("LOG_BACKUP_COUNT", 5))
LOG_TO_CONSOLE = _env_bool("LOG_TO_CONSOLE", "1")
LOG_TO_FILE = _env_bool("LOG_TO_FILE", "0" if DEBUG else "1")
APP_LOG_FILE = LOG_DIR / os.getenv("APP_LOG_FILE", "application.log")
DJANGO_REQUEST_LOG_FILE = LOG_DIR / os.getenv(
    "DJANGO_REQUEST_LOG_FILE", "django_requests.log"
)

ENABLE_SMTP_LOGS = _env_bool("LOG_ENABLE_SMTP", "0")
SMTP_HOST = os.getenv("LOG_SMTP_HOST", "")
SMTP_PORT = int(os.getenv("LOG_SMTP_PORT", 587))
SMTP_USERNAME = os.getenv("LOG_SMTP_USERNAME")
SMTP_PASSWORD = os.getenv("LOG_SMTP_PASSWORD")
SMTP_FROM = os.getenv("LOG_SMTP_FROM", SMTP_USERNAME or "")
SMTP_TO = _split_env_list("LOG_SMTP_TO", "")
SMTP_SUBJECT = os.getenv("LOG_SMTP_SUBJECT", "MetaGap error notification")
SMTP_USE_TLS = _env_bool("LOG_SMTP_USE_TLS", "1")

APM_HOST = os.getenv("LOG_APM_HOST", "")
APM_URL = os.getenv("LOG_APM_URL", "/")
APM_METHOD = os.getenv("LOG_APM_METHOD", "POST")
APM_SECURE = _env_bool("LOG_APM_USE_TLS", "1")
APM_LEVEL = os.getenv("LOG_APM_LEVEL", "ERROR").upper()

handlers: dict[str, dict] = {}
if LOG_TO_CONSOLE:
    handlers["console"] = {
        "class": "logging.StreamHandler",
        "formatter": "verbose",
    }

if LOG_TO_FILE:
    handlers.update(
        {
            "app_file": {
                "class": "logging.handlers.RotatingFileHandler",
                "formatter": "verbose",
                "filename": APP_LOG_FILE,
                "maxBytes": LOG_MAX_BYTES,
                "backupCount": LOG_BACKUP_COUNT,
            },
            "django_request_file": {
                "class": "logging.handlers.RotatingFileHandler",
                "formatter": "verbose",
                "filename": DJANGO_REQUEST_LOG_FILE,
                "maxBytes": LOG_MAX_BYTES,
                "backupCount": LOG_BACKUP_COUNT,
            },
        }
    )

error_handlers: list[str] = []
if ENABLE_SMTP_LOGS and SMTP_HOST and SMTP_TO:
    handlers["smtp"] = {
        "class": "logging.handlers.SMTPHandler",
        "level": os.getenv("LOG_SMTP_LEVEL", "ERROR").upper(),
        "formatter": "verbose",
        "mailhost": (SMTP_HOST, SMTP_PORT),
        "fromaddr": SMTP_FROM,
        "toaddrs": SMTP_TO,
        "subject": SMTP_SUBJECT,
        "credentials": (SMTP_USERNAME, SMTP_PASSWORD)
        if SMTP_USERNAME and SMTP_PASSWORD
        else None,
        "secure": () if SMTP_USE_TLS else None,
    }
    error_handlers.append("smtp")

if APM_HOST:
    handlers["apm_http"] = {
        "class": "logging.handlers.HTTPHandler",
        "level": APM_LEVEL,
        "formatter": "verbose",
        "host": APM_HOST,
        "url": APM_URL,
        "method": APM_METHOD,
        "secure": APM_SECURE,
    }
    error_handlers.append("apm_http")

root_handlers: list[str] = []
if LOG_TO_CONSOLE:
    root_handlers.append("console")
if LOG_TO_FILE:
    root_handlers.append("app_file")
root_handlers.extend(error_handlers)

request_handlers: list[str] = []
if LOG_TO_CONSOLE:
    request_handlers.append("console")
if LOG_TO_FILE:
    request_handlers.append("django_request_file")
request_handlers.extend(error_handlers)

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "verbose": {
            "format": "%(asctime)s [%(levelname)s] %(name)s:%(lineno)d %(message)s",
        }
    },
    "handlers": handlers,
    "root": {"handlers": root_handlers, "level": LOG_LEVEL},
    "loggers": {
        "django.request": {
            "handlers": request_handlers,
            "level": REQUEST_LOG_LEVEL,
            "propagate": False,
        }
    },
}
