# views.py

from django.views.generic import FormView, ListView, CreateView, TemplateView
from django.contrib.auth.mixins import LoginRequiredMixin
from django.contrib.auth.decorators import login_required
from django.utils.decorators import method_decorator
from django.db.models import Q
from django.urls import reverse_lazy
from django.shortcuts import render
from django_tables2 import SingleTableView, RequestConfig
from .models import AlleleFrequency, SampleGroup
from .forms import CustomUserCreationForm, SearchForm, ImportDataForm
from django.contrib.auth.models import User
from django.contrib import messages
from django.views.generic.edit import FormView
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.files.storage import default_storage
from django.conf import settings
import os
from .forms import (
    CustomUserCreationForm,
    SearchForm,
    ImportDataForm,
    EditProfileForm,
    DeleteAccountForm,
    SampleGroupForm
)
# views.py
from django.shortcuts import render
from django_tables2 import RequestConfig
from .models import AlleleFrequency, SampleGroup
from django.views.generic import ListView
from django_tables2 import RequestConfig
from .tables import create_dynamic_table

class SampleGroupTableView(ListView):
    model = SampleGroup
    template_name = "sample_group_table.html"
    context_object_name = "sample_groups"
    paginate_by = 10

    def get_queryset(self):
        # Optimize query by selecting related AlleleFrequency objects
        return super().get_queryset().select_related("allele_frequency")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        # Get the dynamically created table class
        SampleGroupTable = create_combined_table(SampleGroup)
        table = SampleGroupTable(self.get_queryset())
        RequestConfig(self.request, paginate={"per_page": self.paginate_by}).configure(table)
        context["table"] = table
        return context

from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin

class SearchResultsView(SingleTableMixin, FilterView):
    template_name = 'results.html'
    model = SampleGroup
    context_object_name = 'sample_groups'
    paginate_by = 10

    # Dynamically create a table class that includes all related fields.
    table_class = create_dynamic_table(SampleGroup, table_name="CombinedTable", include_related=True)

    def get_queryset(self):
        # Adjust the queryset if needed (for example, using select_related/prefetch_related)
        return SampleGroup.objects.all()


class HomePageView(TemplateView):
    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["form"] = SearchForm()
        return context
    def get(self, request, *args, **kwargs):
        form = SearchForm()
        return self.render_to_response({'form': form})

class UserRegistrationView(CreateView):
    form_class = CustomUserCreationForm
    template_name = "signup.html"
    success_url = reverse_lazy("login")


class ContactView(TemplateView):
    template_name = "contact.html"


class AboutView(TemplateView):
    template_name = "about.html"

class ProfileView(LoginRequiredMixin, TemplateView):
    template_name = "profile.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["import_form"] = ImportDataForm()
        context["sample_groups"] = SampleGroup.objects.select_related(
            "created_by", "reference_genome_build", "allele_frequency"
        ).filter(created_by__user=self.request.user)
        return context


class ImportDataView(LoginRequiredMixin, FormView):
    template_name = "import_data.html"
    form_class = ImportDataForm
    success_url = reverse_lazy("profile")

    def form_valid(self, form):
        data_file = form.cleaned_data["data_file"]
        # Save the uploaded file temporarily
        file_path = default_storage.save("tmp/" + data_file.name, data_file)
        full_path = os.path.join(settings.MEDIA_ROOT, file_path)
        try:
            self.parse_vcf_file(full_path, self.request.user)
            messages.success(self.request, "Data imported successfully.")
        except Exception as e:
            messages.error(self.request, f"An error occurred: {e}")
        finally:
            # Delete the temporary file
            default_storage.delete(file_path)
        return super().form_valid(form)


def parse_vcf_file(self, file_path, user):
    """
    Parses the VCF file at file_path and saves data to the database.

    Args:
        file_path (str): The path to the VCF file.
        user (User): The user who uploaded the file.

    Raises:
        Exception: Any exception during parsing or saving data.
    """
    # Open the VCF file using pysam
    vcf_in = pysam.VariantFile(file_path)  # auto-detect input format (VCF or BCF)

    # Extract sample group metadata from VCF header
    sample_group_metadata = self.extract_sample_group_metadata(vcf_in)
    # Create a SampleGroup object
    sample_group = SampleGroup.objects.create(
        name=sample_group_metadata.get("name", "Sample Group"),
        created_by=user,
        doi=sample_group_metadata.get("doi"),
        source_lab=sample_group_metadata.get("source_lab"),
        contact_email=sample_group_metadata.get("contact_email"),
        # Include other metadata fields as needed
    )

    # Iterate over records and save AlleleFrequency objects
    allele_frequency_objects = []
    for record in vcf_in.fetch():
        info_dict = dict(record.info)
        allele_frequency = AlleleFrequency(
            sample_group=sample_group,
            chrom=record.chrom,
            pos=record.pos,
            variant_id=record.id,
            ref=record.ref,
            alt=",".join(record.alts),
            qual=record.qual,
            filter=",".join(record.filter.keys()) if record.filter.keys() else None,
            info=info_dict,
        )
        allele_frequency_objects.append(allele_frequency)

    # Bulk create allele frequencies for efficiency
    AlleleFrequency.objects.bulk_create(allele_frequency_objects)


def extract_sample_group_metadata(self, vcf_in):
    """
    Extracts sample group metadata from the VCF header.

    Args:
        vcf_in (pysam.VariantFile): The opened VCF file.

    Returns:
        dict: A dictionary containing metadata.
    """
    metadata = {}
    header = vcf_in.header
    # Custom metadata is often stored in the header lines starting with '##'
    # For example: ##SAMPLE=<ID=Sample1,Description="Sample description">
    for record in header.records:
        if record.key == "META":  # Replace 'META' with your actual metadata key
            for key, value in record.items():
                metadata[key.lower()] = value
    return metadata