"""Custom exceptions raised by the VCF importer services."""

from __future__ import annotations

from typing import Iterable, List

from django.utils.translation import gettext_lazy as _


GENERIC_FALLBACK_VALIDATION_MESSAGE_RAW = (
    "The uploaded VCF file appears to be invalid or corrupted. "
    "Please verify the file contents and try again."
)
GENERIC_FALLBACK_VALIDATION_MESSAGE = _(GENERIC_FALLBACK_VALIDATION_MESSAGE_RAW)


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
    "GENERIC_FALLBACK_VALIDATION_MESSAGE",
    "GENERIC_FALLBACK_VALIDATION_MESSAGE_RAW",
    "ImporterConfigurationError",
    "ImporterError",
    "ImporterValidationError",
]
