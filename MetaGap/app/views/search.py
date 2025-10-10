"""Views that power searching and high-level navigation."""

from __future__ import annotations

from django.contrib.auth.mixins import LoginRequiredMixin
from django.views.generic import TemplateView
from django_filters.views import FilterView
from django_tables2.views import SingleTableMixin

from ..forms import SearchForm
from ..models import AlleleFrequency, SampleGroup
from ..tables import build_allele_frequency_table
from ..filters import AlleleFrequencySearchFilter


class SearchResultsView(SingleTableMixin, FilterView):
    template_name = "results.html"
    model = AlleleFrequency
    context_object_name = "allele_frequencies"
    paginate_by = 10
    filterset_class = AlleleFrequencySearchFilter
    ordering = ("chrom", "pos", "ref", "alt", "pk")

    # Dynamically create a table class prioritising variant descriptors first.
    table_class = build_allele_frequency_table()

    def get_queryset(self):
        base_queryset = (
            AlleleFrequency.objects.select_related(
                "info",
                "format",
                "sample_group",
                "sample_group__reference_genome_build",
                "sample_group__genome_complexity",
                "sample_group__sample_origin",
                "sample_group__material_type",
                "sample_group__library_construction",
                "sample_group__illumina_seq",
                "sample_group__ont_seq",
                "sample_group__pacbio_seq",
                "sample_group__iontorrent_seq",
                "sample_group__bioinfo_alignment",
                "sample_group__bioinfo_variant_calling",
                "sample_group__bioinfo_post_proc",
                "sample_group__input_quality",
                "sample_group__created_by",
            )
        )

        # Instantiate the filter set manually so that we can reuse it in the template context.
        self.filterset = self.filterset_class(
            self.request.GET or None, queryset=base_queryset.order_by(*self.ordering)
        )
        return self.filterset.qs.order_by(*self.ordering).distinct()

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault("filter", getattr(self, "filterset", None))
        if hasattr(self, "filterset"):
            context["filter_form"] = self.filterset.form
        return context


class HomePageView(TemplateView):
    """Landing page that exposes the primary search form."""

    template_name = "index.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["form"] = SearchForm()
        return context


class DashboardView(LoginRequiredMixin, TemplateView):
    """Authenticated dashboard summarising recent datasets and activity."""

    template_name = "dashboard.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        organization_profile = getattr(self.request.user, "organization_profile", None)

        dataset_queryset = (
            SampleGroup.objects.select_related("created_by", "created_by__user")
            .order_by("-pk")
        )
        action_queryset = (
            AlleleFrequency.objects.select_related(
                "sample_group",
                "sample_group__created_by",
                "sample_group__created_by__user",
            )
            .order_by("-pk")
        )

        if organization_profile is not None:
            dataset_queryset = dataset_queryset.filter(created_by=organization_profile)
            action_queryset = action_queryset.filter(
                sample_group__created_by=organization_profile
            )

        context.update(
            {
                "recent_datasets": list(dataset_queryset[:6]),
                "recent_actions": list(action_queryset[:6]),
            }
        )
        return context


__all__ = ["SearchResultsView", "HomePageView", "DashboardView"]
