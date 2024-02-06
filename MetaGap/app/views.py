from tkinter import SEL
from django.db import models
from django.views.generic import ListView, CreateView, TemplateView
from django.contrib.auth.mixins import LoginRequiredMixin
from django.db.models.functions import Cast
from .models import PopulationFrequency, Common, UserData
from .forms import (
	CustomUserCreationForm,
	PopulationFrequencyForm,
	CommonForm,
	SearchForm,
)
from django.urls import reverse_lazy
from django.db.models import Q, CharField


class CommonListView(ListView):
	model = Common
	template_name = "common_list.html"
	context_object_name = "commons"


class CommonCreateView(CreateView):
	model = Common
	form_class = CommonForm
	template_name = "form_template.html"
	success_url = reverse_lazy("common_list")


# Additional views for search and profile management
class HomePageView(TemplateView):
	template_name = "index.html"

	def get_context_data(self, **kwargs):
		context = super().get_context_data(**kwargs)
		context["form"] = SearchForm()
		return context


class SearchResultView(ListView):
	template_name = "results.html"
	model = Common
	context_object_name = "results"
	def get_queryset(self):
		form = SearchForm(self.request.GET)
		if form.is_valid():
			query = form.cleaned_data['query']
			query_list = query.split()
			c_queries = Q()
			pf_queries = Q()
			for word in query_list:
				c_queries |= self.model_query(Common, word)
				pf_queries |= self.model_query(PopulationFrequency, word)
			com_results = Common.objects.filter(c_queries)
			pf_results = PopulationFrequency.objects.filter(pf_queries)
			return list(com_results) + list(pf_results)
		else:
			return Common.objects.none()
	def model_query(self, model, query_word):
		queries = Q()
		for field in model._meta.fields:
			if isinstance(field, (models.CharField, models.TextField)):
				queries |= Q(**{f"{field.name}__icontains": query_word})
			# else:
			# 	# Cast the field to a string and search
			# 	queries |= Q(**{f'{field.name}__icontains': Cast(field.name, output_field=CharField())})
		return queries
	def get_context_data(self, **kwargs):
		context = super().get_context_data(**kwargs)
		context['search_form'] = SearchForm(self.request.GET or None)
		return context

class ProfileView(LoginRequiredMixin, TemplateView):
	template_name = "profile.html"

	def get_context_data(self, **kwargs):
		context = super().get_context_data(**kwargs)
		user_data = UserData.objects.filter(user=self.request.user).first()
		if user_data:
			context["common_list"] = user_data.commons.all()
		return context


class AddDataView(LoginRequiredMixin, CreateView):
	model = Common
	template_name = "adddata.html"
	form_class = CommonForm
	second_form_class = PopulationFrequencyForm
	success_url = reverse_lazy("profile")

	def get_context_data(self, **kwargs):
		context = super().get_context_data(**kwargs)
		if "form2" not in context:
			context["form2"] = self.second_form_class()
		return context

	def post(self, request, *args, **kwargs):
		self.object = None
		form = self.get_form()
		form2 = self.second_form_class(request.POST)
		if form.is_valid() and form2.is_valid():
			return self.form_valid(form, form2)
		else:
			return self.form_invalid(form, form2)

	def form_valid(self, form, form2):
		pf = form2.save()
		com = form.save(commit=False)
		com.population_frequency = pf
		com.save()
		response = super().form_valid(form)
		# Handle UserData association as before
		user_data, created = UserData.objects.get_or_create(user=self.request.user)
		user_data.commons.add(self.object)
		return response

	def form_invalid(self, form, form2):
		return self.render_to_response(self.get_context_data(form=form, form2=form2))


class UserRegistrationView(CreateView):
	form_class = CustomUserCreationForm
	template_name = "signup.html"
	success_url = reverse_lazy("login")


class ContactView(TemplateView):
	template_name = "contact.html"


class AboutView(TemplateView):
	template_name = "about.html"


# Ensure to update URL configurations to include paths for these new views.
