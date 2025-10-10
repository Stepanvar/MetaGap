"""Compatibility facade for the MetaGap VCF merging utilities."""
from __future__ import annotations

from MetaGap.MetagapUserCode import merge_vcf as _impl

from MetaGap.MetagapUserCode.merge_vcf import *  # noqa: F401,F403

__all__ = getattr(_impl, "__all__", [])
