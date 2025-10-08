from django.contrib import admin

from .models import (
    AlleleFrequency,
    BioinfoAlignment,
    BioinfoPostProc,
    BioinfoVariantCalling,
    Format,
    GenomeComplexity,
    IlluminaSeq,
    Info,
    InputQuality,
    IonTorrentSeq,
    LibraryConstruction,
    MaterialType,
    OntSeq,
    OrganizationProfile,
    PacBioSeq,
    PlatformIndependent,
    ReferenceGenomeBuild,
    SampleGroup,
    SampleOrigin,
)

admin.site.register(OrganizationProfile)
admin.site.register(ReferenceGenomeBuild)
admin.site.register(GenomeComplexity)
admin.site.register(SampleOrigin)
admin.site.register(MaterialType)
admin.site.register(LibraryConstruction)
admin.site.register(IlluminaSeq)
admin.site.register(OntSeq)
admin.site.register(PacBioSeq)
admin.site.register(IonTorrentSeq)
admin.site.register(PlatformIndependent)
admin.site.register(BioinfoAlignment)
admin.site.register(BioinfoVariantCalling)
admin.site.register(BioinfoPostProc)
admin.site.register(InputQuality)
admin.site.register(SampleGroup)
admin.site.register(AlleleFrequency)
admin.site.register(Info)
admin.site.register(Format)
