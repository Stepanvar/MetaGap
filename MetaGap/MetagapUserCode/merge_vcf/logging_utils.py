"""Shared logging helpers for the VCF merging workflow.

The :func:`configure_logging` helper exposes a small, well-documented surface for
customising how the merger reports progress. Downstream tools can reuse the
default configuration or redirect logs to a specific location::

    from MetaGap.MetagapUserCode.merge_vcf.logging_utils import configure_logging

    configure_logging(log_level="WARNING", log_file="/tmp/merge.log")

The helper is idempotent and clears previously registered handlers so repeated
configuration does not accumulate duplicate outputs.
"""

from __future__ import annotations

import logging
import os
from typing import Iterable

LOG_FILE = "script_execution.log"
"""Default filename for log output when file logging is enabled."""

LOG_FORMAT = "%(asctime)s : %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

logger = logging.getLogger("vcf_merger")
logger.propagate = False


def _normalize_level(level: int | str) -> int:
    """Return a numeric logging level for *level*.

    Accepts either an integer logging constant or a string such as ``"INFO"``.
    """

    if isinstance(level, str):
        upper = level.upper()
        if upper not in logging._nameToLevel:  # type: ignore[attr-defined]
            raise ValueError(f"Unknown log level: {level}")
        return logging._nameToLevel[upper]  # type: ignore[attr-defined]
    return level


def _clear_handlers(existing: Iterable[logging.Handler]) -> None:
    for handler in list(existing):
        try:
            handler.close()
        finally:
            logger.removeHandler(handler)


def configure_logging(
    *,
    log_level: int | str = logging.DEBUG,
    log_file: str | os.PathLike[str] | None = LOG_FILE,
    enable_file_logging: bool = True,
) -> None:
    """Configure logging for the VCF merger.

    Parameters
    ----------
    log_level:
        Logging level to apply to the ``vcf_merger`` logger. Accepts both the
        numeric constants from :mod:`logging` and their string counterparts
        (e.g., ``"INFO"`` or ``"WARNING"``).
    log_file:
        Destination path for file-based logging. Ignored when
        ``enable_file_logging`` is ``False``. Paths provided as
        :class:`~os.PathLike` objects are resolved with :func:`os.fspath`.
    enable_file_logging:
        When ``True`` (the default) a :class:`logging.FileHandler` is attached
        in addition to the console stream handler. Set to ``False`` to disable
        file logging entirely.
    """

    resolved_level = _normalize_level(log_level)

    _clear_handlers(logger.handlers)

    logger.setLevel(resolved_level)

    formatter = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    if enable_file_logging and log_file is not None:
        file_handler = logging.FileHandler(os.fspath(log_file))
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

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
    "configure_logging",
    "handle_critical_error",
    "MergeVCFError",
    "ValidationError",
    "MergeConflictError",
    "handle_non_critical_error",
]


configure_logging()
