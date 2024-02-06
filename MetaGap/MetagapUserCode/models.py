"""
Definition of models.
"""

from django.db import models
from django.contrib.auth.models import User

# Create your models here.

from django.db import models
from django.contrib.auth.models import User


class PopulationFrequency(models.Model):
    name = models.CharField(max_length=255, blank=True)
    locus = models.CharField(max_length=255, blank=True)
    chromosome = models.IntegerField(default=0)
    position = models.IntegerField(default=0)
    minor_allele_frequency = models.FloatField(default=0.0)
    minor_allele = models.CharField(max_length=255, blank=True)
    major_allele = models.CharField(max_length=255, blank=True)
    major_allele_frequency = models.FloatField(default=0.0)
    description = models.TextField(blank=True)

    def __str__(self):
        return self.name


class Genotype(models.Model):
    name = models.CharField(max_length=255, blank=True)
    description = models.TextField(blank=True)

    def __str__(self):
        return self.name


class Common(models.Model):
    name = models.CharField(max_length=255, blank=True)
    description = models.TextField()
    contact_info = models.CharField(max_length=255)
    phen = models.ForeignKey(
        PopulationFrequency, on_delete=models.CASCADE, related_name="commons_phen"
    )
    gen = models.ForeignKey(
        Genotype, on_delete=models.CASCADE, related_name="commons_gen"
    )

    def __str__(self):
        return self.name


class UserData(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    commons = models.ManyToManyField(Common, related_name="user_data")

    def __str__(self):
        return self.user.username
