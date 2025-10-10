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
