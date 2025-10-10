"""Views related to user accounts and organisation profiles."""

from __future__ import annotations

from django.contrib import messages
from django.contrib.auth import logout
from django.contrib.auth.mixins import LoginRequiredMixin
from django.urls import reverse_lazy
from django.views.generic import CreateView, FormView, TemplateView

from ..forms import (
    CustomUserCreationForm,
    DeleteAccountForm,
    EditProfileForm,
    ImportDataForm,
)
from ..mixins import OrganizationSampleGroupMixin


class ProfileView(LoginRequiredMixin, OrganizationSampleGroupMixin, TemplateView):
    """Display the user's organisation profile and related import tools."""

    template_name = "profile.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        organization_profile = self.get_organization_profile()
        sample_groups = self.get_owned_sample_groups().order_by("name")

        context.update(
            {
                "organization_profile": organization_profile,
                "sample_groups": sample_groups,
                "import_form": ImportDataForm(),
                "import_form_action": reverse_lazy("import_data"),
                "import_form_enctype": "multipart/form-data",
            }
        )
        return context


class EditProfileView(LoginRequiredMixin, FormView):
    form_class = EditProfileForm
    template_name = "edit_profile.html"
    success_url = reverse_lazy("profile")

    def get_form_kwargs(self):
        kwargs = super().get_form_kwargs()
        kwargs["instance"] = self.request.user
        return kwargs

    def form_valid(self, form):
        form.save()
        messages.success(self.request, "Your profile has been updated.")
        return super().form_valid(form)


class DeleteAccountView(LoginRequiredMixin, FormView):
    form_class = DeleteAccountForm
    template_name = "confirm_delete.html"
    success_url = reverse_lazy("home")

    def form_valid(self, form):
        user = self.request.user
        logout(self.request)
        user.delete()
        messages.success(self.request, "Your account has been deleted.")
        return super().form_valid(form)

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context.setdefault(
            "form_extra_actions",
            [
                {
                    "url": reverse_lazy("profile"),
                    "class": "btn btn-secondary",
                    "label": "Cancel",
                    "text": "Cancel",
                }
            ],
        )
        return context


class UserRegistrationView(CreateView):
    form_class = CustomUserCreationForm
    template_name = "signup.html"
    success_url = reverse_lazy("login")


__all__ = [
    "ProfileView",
    "EditProfileView",
    "DeleteAccountView",
    "UserRegistrationView",
]
