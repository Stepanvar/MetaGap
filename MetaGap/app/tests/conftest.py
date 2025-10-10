"""Shared pytest fixtures for the app test suite."""

from __future__ import annotations

import os
from types import SimpleNamespace
from uuid import uuid4

import django
import pytest
from django.apps import apps
from django.conf import settings
from django.contrib.auth import get_user_model
from django.test.utils import setup_test_environment, teardown_test_environment

if not settings.configured:
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "MetaGap.settings")

if not apps.ready:
    django.setup()

if "testserver" not in settings.ALLOWED_HOSTS:
    settings.ALLOWED_HOSTS.append("testserver")


@pytest.fixture(scope="session", autouse=True)
def django_test_environment():
    """Ensure Django's test environment helpers are initialised."""

    setup_test_environment()
    yield
    teardown_test_environment()


@pytest.fixture
def allele_frequency_filter_data() -> SimpleNamespace:
    """Construct a set of allele frequency records used across filter tests."""

    User = get_user_model()
    user = User.objects.create_user(
        username=f"filter-user-{uuid4().hex}", password="pass12345"
    )

    from app.models import (
        AlleleFrequency,
        BioinfoAlignment,
        BioinfoVariantCalling,
        Info,
        SampleGroup,
        SampleOrigin,
    )

    AlleleFrequency.objects.all().delete()
    SampleGroup.objects.all().delete()
    Info.objects.all().delete()
    BioinfoAlignment.objects.all().delete()
    BioinfoVariantCalling.objects.all().delete()
    SampleOrigin.objects.all().delete()

    sample_origin = SampleOrigin.objects.create(tissue="Liver")
    other_origin = SampleOrigin.objects.create(tissue="Blood")

    alignment_primary = BioinfoAlignment.objects.create(
        tool="AlignerA",
        params="--fast",
        ref_genome_version="GRCh38",
        recalibration_settings="BQSR",
    )
    alignment_secondary = BioinfoAlignment.objects.create(
        tool="AlignerB",
        params="--sensitive",
        ref_genome_version="GRCh37",
        recalibration_settings="None",
    )

    variant_primary = BioinfoVariantCalling.objects.create(
        tool="CallerA",
        version="1.0",
        filtering_thresholds="QD>10",
        duplicate_handling="Remove",
        mq="60",
    )
    variant_secondary = BioinfoVariantCalling.objects.create(
        tool="CallerB",
        version="2.1",
        filtering_thresholds="QD>15",
        duplicate_handling="Mark",
        mq="50",
    )

    sample_group = SampleGroup.objects.create(
        name=f"Population A {uuid4().hex}",
        created_by=user.organization_profile,
        source_lab="North Lab",
        sample_origin=sample_origin,
        bioinfo_alignment=alignment_primary,
        bioinfo_variant_calling=variant_primary,
    )
    other_group = SampleGroup.objects.create(
        name=f"Population B {uuid4().hex}",
        created_by=user.organization_profile,
        source_lab="South Lab",
        sample_origin=other_origin,
        bioinfo_alignment=alignment_secondary,
        bioinfo_variant_calling=variant_secondary,
    )

    info_chr1 = Info.objects.create(
        af="0.12",
        dp="25",
        mq="55",
        additional={"QD": "15.0", "FS": "0.4", "SOR": "1.5"},
    )
    info_chr2 = Info.objects.create(
        af=None,
        dp=None,
        mq=None,
        additional=None,
    )
    info_chr3 = Info.objects.create(
        af="0.45",
        dp="100",
        mq="70",
        additional={"QD": "20.5", "FS": "3.1", "SOR": "0.9"},
    )

    record_chr1 = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom="chr1",
        pos=12345,
        variant_id="rs123",
        ref="A",
        alt="T",
        qual=75.0,
        filter="PASS",
        info=info_chr1,
    )
    record_chr2 = AlleleFrequency.objects.create(
        sample_group=sample_group,
        chrom="chr2",
        pos=54321,
        variant_id="rs543",
        ref="G",
        alt="C",
        qual=15.0,
        filter="LowQual",
        info=info_chr2,
    )
    record_chr3 = AlleleFrequency.objects.create(
        sample_group=other_group,
        chrom="chr3",
        pos=88888,
        variant_id="rs789",
        ref="C",
        alt="G",
        qual=120.0,
        filter="PASS",
        info=info_chr3,
    )

    return SimpleNamespace(
        user=user,
        sample_group=sample_group,
        other_group=other_group,
        record_chr1=record_chr1,
        record_chr2=record_chr2,
        record_chr3=record_chr3,
    )


@pytest.fixture
def sample_group_filter_data() -> SimpleNamespace:
    """Create sample groups with different origins for filter coverage."""

    from app.models import SampleGroup, SampleOrigin

    SampleGroup.objects.all().delete()
    SampleOrigin.objects.all().delete()

    User = get_user_model()
    user = User.objects.create_user(
        username=f"sample-user-{uuid4().hex}", password="pass12345"
    )

    sample_origin = SampleOrigin.objects.create(
        tissue="Kidney",
        collection_method="Biopsy",
        storage_conditions="Frozen",
        time_stored="2 years",
    )
    other_origin = SampleOrigin.objects.create(
        tissue="Blood",
        collection_method="Draw",
        storage_conditions="Refrigerated",
    )

    group = SampleGroup.objects.create(
        name=f"Kidney Cohort {uuid4().hex}",
        created_by=user.organization_profile,
        sample_origin=sample_origin,
    )
    other_group = SampleGroup.objects.create(
        name=f"Control Cohort {uuid4().hex}",
        created_by=user.organization_profile,
        sample_origin=other_origin,
    )

    return SimpleNamespace(
        group=group,
        other_group=other_group,
        sample_origin=sample_origin,
        other_origin=other_origin,
    )
