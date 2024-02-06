from django.views.generic import ListView, CreateView, TemplateView, UpdateView
from django.contrib.auth.mixins import LoginRequiredMixin
from .models import PopulationFrequency, Genotype, Common
from .forms import PopulationFrequencyForm, GenotypeForm, CommonForm, SearchForm
from django.urls import reverse_lazy
from django.db.models import Q
from django.shortcuts import redirect


# Existing views
class PopulationFrequencyListView(ListView):
    model = PopulationFrequency
    template_name = "index.html"
    context_object_name = "population_frequencies"


class PopulationFrequencyCreateView(CreateView):
    model = PopulationFrequency
    form_class = PopulationFrequencyForm
    template_name = (
        "form_template.html"  # Ensure this template name matches your project's naming
    )
    success_url = reverse_lazy("population_frequency_list")  # Adjust as necessary


# Additional views for search and profile management
class HomePageView(TemplateView):
    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["form"] = SearchForm()
        return context


class SearchResultView(ListView):
    model = Common  # Adjust the model as necessary for search
    template_name = "results.html"

    def get_queryset(self):
        query = self.request.GET.get("q")
        return Common.objects.filter(
            Q(name__icontains=query) | Q(description__icontains=query)
        )


class ProfileView(LoginRequiredMixin, TemplateView):
    template_name = "profile.html"


class AddDataView(LoginRequiredMixin, CreateView):
    model = PopulationFrequency  # Adjust if you're using a different model for data addition
    form_class = PopulationFrequencyForm  # Adjust the form class as needed
    template_name = "adddata.html"
    success_url = reverse_lazy(
        "/profile/"
    )  # Adjust this to your specific success URL


# Ensure to update URL configurations to include paths for these new views.
