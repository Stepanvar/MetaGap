"""Utilities for importing VCF files into the relational schema."""

from __future__ import annotations

import logging
from typing import Any, Dict, Optional, Tuple

import pysam

from django.core.exceptions import ObjectDoesNotExist, ValidationError
from django.db import transaction

from ..models import Format, Info, OrganizationProfile, SampleGroup
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

        try:
            organization_profile = self.user.organization_profile
        except AttributeError as exc:  # pragma: no cover - defensive fallback
            raise ValidationError(
                "Please complete your organization profile before importing data."
            ) from exc
        except (OrganizationProfile.DoesNotExist, ObjectDoesNotExist) as exc:
            raise ValidationError(
                "Please complete your organization profile before importing data."
            ) from exc

        if organization_profile is None:
            raise ValidationError(
                "Please complete your organization profile before importing data."
            )

        metadata: Dict[str, Any] = {}
        sample_group: Optional[SampleGroup] = None
        with transaction.atomic():
            try:
                with pysam.VariantFile(file_path) as vcf_in:
                    metadata = self.extract_sample_group_metadata(vcf_in)
                    sample_group = self._create_sample_group(
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
                    sample_group = None
                metadata = self._extract_metadata_text_fallback(file_path)
                sample_group = self._create_sample_group(
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
            try:
                metadata = self._extract_metadata_text_fallback(file_path)
                sample_group = self._create_sample_group(
                    metadata, file_path, organization_profile
                )
                parse_vcf_text_fallback(
                    file_path, sample_group, self.database_writer
                )
            except (UnicodeDecodeError, ValueError, TypeError) as fallback_exc:
                if sample_group is not None:
                    sample_group.delete()
                raise ValidationError(
                    "The uploaded VCF file appears to be invalid or corrupted. "
                    "Please verify the file contents and try again."
                ) from fallback_exc

        assert sample_group is not None
        return sample_group

    def _populate_sample_group_from_pysam(
        self, vcf_in: pysam.VariantFile, sample_group: SampleGroup
    ) -> None:
        for record in vcf_in.fetch():
            info_instance = self._create_info_instance(record.info)
            format_instance, format_sample = self._create_format_instance(
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

    # ------------------------------------------------------------------
    # Wrapper helpers used by the unit tests
    # ------------------------------------------------------------------
    def _create_info_instance(self, info: Any) -> Optional[Info]:
        return self.database_writer.create_info_instance(info)

    def _create_format_instance(
        self, samples: Any
    ) -> Tuple[Optional[Format], Optional[str]]:
        return self.database_writer.create_format_instance(samples)

    def _create_sample_group(
        self, metadata: Dict[str, Any], file_path: str, organization_profile: Any
    ) -> SampleGroup:
        return self.database_writer.create_sample_group(
            metadata, file_path, organization_profile
        )

    def _extract_metadata_text_fallback(self, file_path: str) -> Dict[str, Any]:
        return extract_metadata_text_fallback(file_path)

    def _extract_section_data(
        self, metadata: Dict[str, Any], section: str, model_cls: Any
    ) -> Tuple[Dict[str, Any], set[str], Optional[Dict[str, Any]]]:
        return self.database_writer._extract_section_data(metadata, section, model_cls)

    def extract_sample_group_metadata(
        self, vcf_in: pysam.VariantFile
    ) -> Dict[str, Any]:
        return self.metadata_parser.extract_sample_group_metadata(vcf_in)
