"""Shared logging helpers for the VCF merging workflow."""

from __future__ import annotations

import logging

LOG_FILE = "script_execution.log"
logger = logging.getLogger("vcf_merger")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
fh = logging.FileHandler(LOG_FILE)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)


class MergeVCFError(RuntimeError):
    """Base exception for unrecoverable errors in the merge workflow."""


class ValidationError(MergeVCFError):
    """Raised when validation detects an unrecoverable problem."""


class MergeConflictError(MergeVCFError):
    """Raised when merging fails due to conflicting inputs or tooling errors."""


def log_message(message: str, verbose: bool = False) -> None:
    """Log *message* and optionally echo it to stdout."""

    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message: str, exc_cls=None) -> None:
    """Log and raise a fatal error before exiting."""

    log_message("CRITICAL ERROR: " + message)
    exception_class = exc_cls or MergeVCFError
    raise exception_class(message)


def handle_non_critical_error(message: str) -> None:
    """Log and print a recoverable error."""

    log_message("WARNING: " + message)
    print("Warning: " + message)


__all__ = [
    "LOG_FILE",
    "logger",
    "log_message",
    "handle_critical_error",
    "MergeVCFError",
    "ValidationError",
    "MergeConflictError",
    "handle_non_critical_error",
]
