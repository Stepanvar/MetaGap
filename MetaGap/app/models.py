"""
Definition of models.
"""

from django.db import models
from django.contrib.auth.models import User

# Create your models here.


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

	def get_fields(self):
		"""Returns a list of all field names on the instance."""
		return [
			field
			for field in self._meta.get_fields()
			if field.concrete and not field.many_to_many
		]

	def get_field_values(self):
		"""Returns a list of values for all fields on the instance, including detailed values for ForeignKey fields."""
		field_values = []
		for field in self.get_fields():
			value = getattr(self, field.name, "")
			if isinstance(value, models.Model):  # It's a ForeignKey
				# Fetch detailed values from the related model instance
				related_field_values = value.get_detailed_field_values()
				field_values.append(f"{field.name}: {related_field_values}")
			else:
				field_values.append(value)
		return field_values

	def get_detailed_field_values(self):
		"""Returns a string combining all field values of the model instance, for detailed representation."""
		details = []
		for field in self._meta.get_fields():
			# Skip reverse relations or many-to-many fields to avoid complex recursion
			if field.auto_created or field.many_to_many:
				continue
			value = getattr(self, field.name, "")
			details.append(f"{field.verbose_name.title()}: {value}")
		return '; '.join(details)

	def __str__(self):
		return self.name


class Genotype(models.Model):
	name = models.CharField(max_length=255, blank=True)
	description = models.TextField(blank=True)

	def get_fields(self):
		"""Returns a list of all field names on the instance."""
		return [
			field
			for field in self._meta.get_fields()
			if field.concrete and not field.many_to_many
		]

	def get_field_values(self):
		"""Returns a list of values for all fields on the instance."""
		field_values = []
		for field in self.get_fields():
			value = getattr(self, field.name, "")
			if isinstance(value, models.Model):  # It's a ForeignKey
				value = str(
					value
				)  # Use str to get a simple representation of the related model
			field_values.append(value)
		return field_values

	def __str__(self):
		return self.name


class Common(models.Model):
	name = models.CharField(max_length=255, blank=True)
	description = models.TextField()
	contact_info = models.CharField(max_length=255)
	population_frequency = models.ForeignKey(
		PopulationFrequency, on_delete=models.CASCADE, related_name="commons_population_frequency", blank=True, default=""
	)

	def get_fields(self):
		"""Returns a list of all field names on the instance."""
		return [
			field
			for field in self._meta.get_fields()
			if field.concrete and not field.many_to_many
		]

	def get_field_values(self):
		"""Returns a list of values for all fields on the instance, including detailed values for ForeignKey fields."""
		field_values = []
		for field in self.get_fields():
			value = getattr(self, field.name, "")
			if isinstance(value, models.Model):  # It's a ForeignKey
				# Fetch detailed values from the related model instance
				related_field_values = value.get_detailed_field_values()
				field_values.append(f"{field.name}: {related_field_values}")
			else:
				field_values.append(value)
		return field_values

	def get_detailed_field_values(self):
		"""Returns a string combining all field values of the model instance, for detailed representation."""
		details = []
		for field in self._meta.get_fields():
			# Skip reverse relations or many-to-many fields to avoid complex recursion
			if field.auto_created or field.many_to_many:
				continue
			value = getattr(self, field.name, "")
			details.append(f"{field.verbose_name.title()}: {value}")
		return '; '.join(details)

	def __str__(self):
		return self.name


class UserData(models.Model):
	user = models.ForeignKey(User, on_delete=models.CASCADE)
	commons = models.ManyToManyField(Common, related_name="user_data")

	def __str__(self):
		return self.user.username


# class DataSubmission(models.Model):
#	 user = models.ForeignKey(User, on_delete=models.CASCADE, related_name='submissions')
#	 title = models.CharField(max_length=255)
#	 description = models.TextField()
#	 file_upload = models.FileField(upload_to='submissions/', blank=True, null=True)
#	 created_at = models.DateTimeField(auto_now_add=True)
#	 updated_at = models.DateTimeField(auto_now=True)

#	 def __str__(self):
#		 return f'{self.title} by {self.user.username}'
