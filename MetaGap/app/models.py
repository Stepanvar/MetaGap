"""
Definition of models.
"""

from email.policy import default
from tkinter import TRUE
from django.db import models
from django.contrib.auth.models import User

# Create your models here.

class Phenotype(models.Model):
    name = models.CharField(max_length=255, blank=True)
    description = models.TextField(default='', blank=True)

class Genotype(models.Model):
    name = models.CharField(max_length=255, blank=True)
    description = models.TextField(default='', blank=True)

class Common(models.Model):
    name = models.CharField(max_length=255, blank=True)
    description = models.TextField()
    contact_info = models.CharField(max_length=255)
    phen = models.ForeignKey(Phenotype, on_delete=models.CASCADE)
    gen = models.ForeignKey(Genotype, on_delete=models.CASCADE)

class UserData(models.Model):
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    commons = models.ManyToManyField(Common)