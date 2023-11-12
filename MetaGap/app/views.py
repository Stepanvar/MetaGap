"""
Definition of views.
"""

from datetime import datetime
from django.contrib.auth.views import LoginView, LogoutView
from django.http import HttpRequest
from django.shortcuts import render, redirect   
from django.contrib.auth.decorators import login_required
from django.contrib.auth import authenticate, login, logout
from django.db.models import Q
from django import template
from app import forms
from .models import Phenotype, Genotype, Common, UserData
register = template.Library()
# try this
@register.filter(name='custom_date')
def custom_date(value):
    return value.strftime("%d %b, %Y")

def home(request):
    """Renders the home page."""
    assert isinstance(request, HttpRequest)
    if request.method == 'POST':
        query = request.POST.get('search', '')
        common_results = Common.objects.filter(
        Q(name__icontains=query) | Q(description__icontains=query) | Q(contact_info__icontains=query) |
        Q(phen__name__icontains=query) | Q(phen__description__icontains=query) |
        Q(gen__name__icontains=query) | Q(gen__description__icontains=query)
    )
        combined_results = []

        for common in common_results:
            result = {
                'name': common.name,
                'phen': common.phen.name,
                'gen': common.gen.name,
                'contact_info': common.contact_info,
                'desc': common.description
            }
            combined_results.append(result)
        return render(request, 'app/results.html', {'results': combined_results})
    else:
        return render(
        request,
        'app/index.html',
        {
            'title':'Home Page',
            'year': datetime.now().year,
        }
    )

def contact(request):
    """Renders the contact page."""
    assert isinstance(request, HttpRequest)
    return render(
        request,
        'app/contact.html',
        {
            'title':'Contact',
            'message':'Your contact page.',
            'year':datetime.now().year,
        }
    )

def about(request):
    """Renders the about page."""
    assert isinstance(request, HttpRequest)
    return render(
        request,
        'app/about.html',
        {
            'title':'About',
            'message':'Your application description page.',
            'year':datetime.now().year,
        }
    )

def results(request):
    results = request.GET.get('results')
    return render(request, 'app/results.html', {'title':'Results',
                                                'results': results})
 
def signup(request):
    if request.method == 'POST':
        form = forms.BootstrapUserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password1')
            user = authenticate(username=username, password=password)
            login(request)
            return redirect('app/home')
    else:
        form = forms.BootstrapUserCreationForm()
    return render(request, 'app/signup.html', {'form': form})

def login(request):
    return LoginView.as_view(
            template_name="app/login.html",
            authentication_form=forms.BootstrapAuthenticationForm,
            extra_context={
                "title": "Log in",
                "year": datetime.now().year,
            },
        ) (request)

def logout(request):
    return LogoutView.as_view(next_page="/") (request)

@login_required
def delete_user(request):
    user = request.user
    user.delete()
    return redirect("app/login")

@login_required(login_url='/login/')
def profile(request):
    if request.method == 'POST':
       return delete_user(request)
    user_data = UserData.objects.filter(user=request.user).first()
    if user_data is not None:
        common_list = user_data.commons.all()
    else:
        common_list = []
    return render(request, 'app/profile.html', {'title':'Profile',
                                                'common_list': common_list})

@login_required
def adddata(request):
    if request.method == 'POST':
        # Create and save the phenotype and genotype instances first
        phenotype = Phenotype(
            name=request.POST['phenotype_name'], 
            description=request.POST['phenotype_description']
        )
        phenotype.save()

        genotype = Genotype(
            name=request.POST['genotype_name'], 
            description=request.POST['genotype_description']
        )
        genotype.save()

        # Now create the common instance with references to the saved phenotype and genotype
        common = Common(
            name=request.POST['common_name'], 
            description=request.POST['common_description'], 
            contact_info=request.POST['common_contact_info'],
            phen=phenotype,
            gen=genotype
        )
        common.save()

        # Associate the common instance with the user
        user_data = UserData.objects.filter(user=request.user).first()
        user_data.commons.add(common)
        return redirect('profile')
    return render(request, 'app/adddata.html')
