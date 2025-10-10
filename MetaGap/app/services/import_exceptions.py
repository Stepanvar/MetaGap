"""Custom exceptions raised by the VCF importer services."""

from __future__ import annotations

from typing import Iterable, List


class ImporterError(Exception):
    """Base exception for importer failures."""

    def __init__(self, message: str, *, warnings: Iterable[str] | None = None) -> None:
        normalized_message = message or "An unexpected importer error occurred."
        super().__init__(normalized_message)
        self.user_message: str = normalized_message
        self.warnings: List[str] = list(warnings or [])


class ImporterValidationError(ImporterError):
    """Raised when uploaded data fails validation checks."""


class ImporterConfigurationError(ImporterValidationError):
    """Raised when the user's account or environment is not ready for imports."""


__all__ = [
    "ImporterConfigurationError",
    "ImporterError",
    "ImporterValidationError",
]
