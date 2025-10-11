"""Filtering utilities for interpreting allowed FILTER values in the VCF merge workflow.

This module centralizes logic that handles the ``FILTER`` field of VCF records.
``DEFAULT_ALLOWED_FILTER_VALUES`` defines the canonical whitelist used when no
overrides are provided. Helper functions such as
:func:`prepare_allowed_filter_values` and :func:`prepare_allowed_filter_set`
normalize user-supplied configurations into ready-to-use tuples or sets for
membership checks during validation and merging.
"""


from __future__ import annotations

from typing import Optional, Sequence, Set, Tuple

DEFAULT_ALLOWED_FILTER_VALUES: Tuple[str, ...] = ("PASS",)


def prepare_allowed_filter_values(values: Optional[Sequence[str]]) -> Tuple[str, ...]:
    """Return a tuple of allowed FILTER codes, falling back to defaults when empty."""
    if values:
        return tuple(values)
    return DEFAULT_ALLOWED_FILTER_VALUES


def prepare_allowed_filter_set(values: Optional[Sequence[str]]) -> Set[str]:
    """Normalize *values* and drop sentinel placeholders like '.', '', or None."""
    normalized = prepare_allowed_filter_values(values)
    return {v for v in normalized if v not in {None, "", "."}}
