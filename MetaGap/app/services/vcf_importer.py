"""Utilities for importing VCF files into the relational schema."""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional

import pysam

from ..models import SampleGroup
from .vcf_database import VCFDatabaseWriter
from .vcf_file_utils import extract_metadata_text_fallback, parse_vcf_text_fallback
from .vcf_metadata import VCFMetadataParser

logger = logging.getLogger(__name__)


class VCFImporter:
    """Encapsulate the VCF parsing workflow."""

    def __init__(self, user: Any) -> None:
        self.user = user
        self.warnings: list[str] = []
        self.metadata_parser = VCFMetadataParser(self.warnings)
        self.database_writer = VCFDatabaseWriter()

    def import_file(self, file_path: str) -> SampleGroup:
        """Import the provided VCF file and return the created sample group."""

        organization_profile = getattr(self.user, "organization_profile", None)
        if organization_profile is None:
            raise ValueError("The current user does not have an organisation profile.")

        metadata: Dict[str, Any] = {}
        sample_group: Optional[SampleGroup] = None
        try:
            with pysam.VariantFile(file_path) as vcf_in:
                metadata = self.metadata_parser.extract_sample_group_metadata(vcf_in)
                sample_group = self.database_writer.create_sample_group(
                    metadata, file_path, organization_profile
                )
                self._populate_sample_group_from_pysam(vcf_in, sample_group)
        except (OSError, ValueError) as exc:
            warning = (
                f"Could not parse VCF metadata with pysam: {exc}. "
                "Falling back to a text parser."
            )
            logger.warning("%s", warning)
            self.warnings.append(warning)
            if sample_group is not None:
                sample_group.delete()
            metadata = extract_metadata_text_fallback(file_path)
            sample_group = self.database_writer.create_sample_group(
                metadata, file_path, organization_profile
            )
            parse_vcf_text_fallback(file_path, sample_group, self.database_writer)

        assert sample_group is not None
        return sample_group

    def _populate_sample_group_from_pysam(
        self, vcf_in: pysam.VariantFile, sample_group: SampleGroup
    ) -> None:
        for record in vcf_in.fetch():
            info_instance = self.database_writer.create_info_instance(record.info)
            format_instance, format_sample = self.database_writer.create_format_instance(
                record.samples
            )

            self.database_writer.create_allele_frequency(
                sample_group,
                chrom=record.chrom,
                pos=record.pos,
                variant_id=record.id,
                ref=record.ref,
                alt=self.database_writer.serialize_alt(record.alts),
                qual=record.qual,
                filter_value=self.database_writer.serialize_filter(record.filter),
                info=info_instance,
                format_instance=format_instance,
                format_sample=format_sample,
            )
