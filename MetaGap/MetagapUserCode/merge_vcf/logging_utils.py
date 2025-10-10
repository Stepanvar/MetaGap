"""Shared logging helpers for the VCF merging workflow."""

from __future__ import annotations

import logging
import os
from typing import Optional

LOG_FILE = "script_execution.log"
logger = logging.getLogger("vcf_merger")
formatter = logging.Formatter(
    "%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)


def configure_logging(verbose: bool, log_file: Optional[str] = None) -> None:
    """(Re)configure the module logger.

    Existing handlers are removed to avoid duplicate log entries. A file handler is
    always attached, while a ``StreamHandler`` is only added when ``verbose`` is
    ``True`` to allow optional console echoing.
    """

    destination = log_file or LOG_FILE

    # Ensure the destination directory exists when a non-empty path is provided.
    directory = os.path.dirname(destination)
    if directory:
        os.makedirs(directory, exist_ok=True)

    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        try:
            handler.close()
        except Exception:  # pragma: no cover - defensive cleanup
            pass

    file_handler = logging.FileHandler(destination)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    if verbose:
        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)


class MergeVCFError(RuntimeError):
    """Base exception for unrecoverable errors in the merge workflow."""


class ValidationError(MergeVCFError):
    """Raised when validation detects an unrecoverable problem."""


class MergeConflictError(MergeVCFError):
    """Raised when merging fails due to conflicting inputs or tooling errors."""


def log_message(message: str, verbose: bool = False) -> None:
    """Log *message* using the configured handlers."""

    logger.info(message)


def handle_critical_error(message: str, exc_cls=None) -> None:
    """Log and raise a fatal error before exiting."""

    log_message("CRITICAL ERROR: " + message)
    exception_class = exc_cls or MergeVCFError
    raise exception_class(message)


def handle_non_critical_error(message: str) -> None:
    """Log a recoverable error."""

    log_message("WARNING: " + message)


__all__ = [
    "LOG_FILE",
    "configure_logging",
    "logger",
    "log_message",
    "handle_critical_error",
    "MergeVCFError",
    "ValidationError",
    "MergeConflictError",
    "handle_non_critical_error",
]


configure_logging(verbose=False)
