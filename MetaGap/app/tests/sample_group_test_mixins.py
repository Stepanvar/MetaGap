"""Shared mixins and helpers for sample group view tests."""

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.urls import NoReverseMatch, reverse

from ..models import (
    AlleleFrequency,
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    Format,
    GenomeComplexity,
    IlluminaSeq,
    Info,
    InputQuality,
    LibraryConstruction,
    MaterialType,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)


class SampleGroupTestDataMixin:
    """Utility helpers for constructing richly populated sample groups."""

    def create_sample_group_with_variant(self, owner, *, name="Detail Cohort"):
        reference = ReferenceGenomeBuild.objects.create(
            build_name="GRCh38",
            build_version="v1",
        )
        genome_complexity = GenomeComplexity.objects.create(
            size="3.2Gb",
            ploidy="Diploid",
            gc_content="41%",
        )
        sample_origin = SampleOrigin.objects.create(
            tissue="Lung",
            collection_method="Biopsy",
            storage_conditions="-80C",
        )
        material_type = MaterialType.objects.create(
            material_type="DNA",
            integrity_number="9.8",
        )
        library_construction = LibraryConstruction.objects.create(
            kit="MetaPrep",
            fragmentation="Acoustic",
            adapter_ligation_efficiency="95%",
        )
        illumina_seq = IlluminaSeq.objects.create(
            instrument="NovaSeq 6000",
            flow_cell="S4",
        )
        bioinfo_alignment = BioinfoAlignment.objects.create(
            tool="BWA",
            ref_genome_version="GRCh38",
        )
        bioinfo_variant_calling = BioinfoVariantCalling.objects.create(
            tool="GATK",
            version="4.2",
        )
        bioinfo_post_proc = BioinfoPostProc.objects.create(
            normalization="bcftools",
        )
        input_quality = InputQuality.objects.create(
            a260_a280=1.8,
            dna_concentration=15.2,
        )

        sample_group = SampleGroup.objects.create(
            name=name,
            created_by=owner.organization_profile,
            doi="10.1000/detail-cohort",
            source_lab="MetaLab",
            contact_email="lab@example.com",
            contact_phone="123-456-7890",
            total_samples=5,
            inclusion_criteria="Adults",
            exclusion_criteria="Under 18",
            comments="Rich metadata snapshot",
            reference_genome_build=reference,
            genome_complexity=genome_complexity,
            sample_origin=sample_origin,
            material_type=material_type,
            library_construction=library_construction,
            sequencing_platform=SampleGroup.SequencingPlatform.ILLUMINA,
            illumina_seq=illumina_seq,
            bioinfo_alignment=bioinfo_alignment,
            bioinfo_variant_calling=bioinfo_variant_calling,
            bioinfo_post_proc=bioinfo_post_proc,
            input_quality=input_quality,
        )

        info = Info.objects.create(
            af="0.5",
            ac="5",
            an="10",
            dp="150",
            mq="60",
            additional={"clinvar_significance": "Pathogenic"},
        )
        fmt = Format.objects.create(
            genotype="0/1",
            payload={
                "fields": {"gq": "99"},
                "additional": {"sample_id": "Sample001"},
            },
        )

        allele = AlleleFrequency.objects.create(
            sample_group=sample_group,
            chrom="1",
            pos=123456,
            variant_id="rsDetail1",
            ref="A",
            alt="T",
            qual=42.0,
            filter="PASS",
            info=info,
            format=fmt,
        )

        return sample_group, allele


class SampleGroupViewMatrixMixin(SampleGroupTestDataMixin):
    """Shared fixtures and helpers for exercising ownership-aware views."""

    DETAIL_ROUTE_CANDIDATES = (
        "sample_group_detail",
        "profile_sample_group_detail",
        "sample-group-detail",
    )
    UPDATE_ROUTE_CANDIDATES = (
        "sample_group_edit",
        "sample_group_update",
        "sample-group-update",
    )
    DELETE_ROUTE_CANDIDATES = (
        "sample_group_delete",
        "sample-group-delete",
    )

    def setUp(self) -> None:
        super().setUp()
        User = get_user_model()
        self.owner = User.objects.create_user(
            username="matrix_owner",
            password="matrix-pass",
            email="matrix-owner@example.com",
        )
        self.intruder = User.objects.create_user(
            username="matrix_intruder",
            password="matrix-pass",
            email="matrix-intruder@example.com",
        )
        self.owned_group, self.owned_allele = self.create_sample_group_with_variant(
            self.owner, name="Matrix Owner Cohort"
        )
        (
            self.intruder_group,
            self.intruder_allele,
        ) = self.create_sample_group_with_variant(
            self.intruder, name="Matrix Intruder Cohort"
        )
        # Provide commonly used aliases for compatibility with existing tests.
        self.sample_group = self.owned_group
        self.allele = self.owned_allele
        self.other_user = self.intruder

    def _resolve_url(self, candidates, *args) -> str:
        for name in candidates:
            try:
                return reverse(name, args=args)
            except NoReverseMatch:
                continue
        self.fail(f"None of the route names {candidates} resolved successfully")

    def get_detail_url(self) -> str:
        return self._resolve_url(self.DETAIL_ROUTE_CANDIDATES, self.owned_group.pk)

    def get_update_url(self) -> str:
        return self._resolve_url(self.UPDATE_ROUTE_CANDIDATES, self.owned_group.pk)

    def get_delete_url(self) -> str:
        return self._resolve_url(self.DELETE_ROUTE_CANDIDATES, self.owned_group.pk)

    def get_import_url(self) -> str:
        return reverse("import_data")
