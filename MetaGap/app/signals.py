"""Connect post-save signal handlers that maintain user organization profiles."""

from django.apps import apps
from django.contrib.auth.models import User
from django.db.models.signals import post_save
from django.dispatch import receiver

@receiver(post_save, sender=User)
def create_user_profile(sender, instance, created, **kwargs):
    if kwargs.get("raw", False):
        return

    if created:
        OrganizationProfile = apps.get_model("app", "OrganizationProfile")
        OrganizationProfile.objects.create(user=instance)

@receiver(post_save, sender=User)
def save_user_profile(sender, instance, **kwargs):
    if kwargs.get("raw", False):
        return

    if hasattr(instance, "organization_profile"):
        instance.organization_profile.save()
