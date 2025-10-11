# Defines the app's AppConfig and imports app.signals in ready() to connect signal handlers.
from django.apps import AppConfig

class AppConfig(AppConfig):
    name = 'app'

    def ready(self):
        import app.signals
