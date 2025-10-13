# Defines the app's AppConfig and imports app.signals in ready() to connect signal handlers.
from pathlib import Path
import logging
from django.apps import AppConfig

logger = logging.getLogger(__name__)


class AppConfig(AppConfig):
    name = 'app'

    def ready(self):
        super().ready()
        locale_root = Path(__file__).resolve().parent / 'locale'
        try:
            from .utils.locale_compiler import ensure_compiled_catalogs
        except ImportError:
            logger.warning('Unable to import locale compiler utility; skipping auto-compilation.')
        else:
            ensure_compiled_catalogs(locale_root)
        import app.signals
