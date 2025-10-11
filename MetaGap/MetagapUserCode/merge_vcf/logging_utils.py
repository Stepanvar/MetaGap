"""Shared logging helpers for the VCF merging workflow.

By default, the module wires the ``vcf_merger`` logger to emit informative
messages both to stdout and to ``script_execution.log`` using a consistent
timestamped format. Importing the module triggers :func:`configure_logging`,
ensuring that merge scripts maintain a persistent file log while surfacing
progress and warnings in the console.

The :func:`configure_logging` helper provides an idempotent, well-documented
entry point for customising logging behaviour. Call it with ``log_level`` to
adjust verbosity, ``log_file`` to redirect output, disable either handler, or
set ``create_dirs=True`` to create missing directories automatically. Repeated
calls clear previous handlers, avoiding duplicate outputs. Example::

    from MetaGap.MetagapUserCode.merge_vcf.logging_utils import configure_logging
    configure_logging(log_level="WARNING", log_file="/tmp/merge.log")

For error handling, the module defines :class:`MergeVCFError` and specialised
subclasses for validation and merge conflicts. Helpers such as
:func:`handle_critical_error` and :func:`handle_non_critical_error` standardise
logging of fatal versus recoverable issues—critical failures are logged at the
``CRITICAL`` level and raised as exceptions, while non-critical conditions are
logged as warnings.
"""
from __future__ import annotations

import logging
import os
from typing import Iterable

LOG_FILE = "script_execution.log"
LOG_FORMAT = "%(asctime)s : %(message)s"
LOG_DATE_FORMAT = "%Y-%m-%d %H:%M:%S"

logger = logging.getLogger("vcf_merger")
logger.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s : %(levelname)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
fh = logging.FileHandler(LOG_FILE)
fh.setFormatter(formatter)
logger.addHandler(fh)
ch = logging.StreamHandler()
ch.setFormatter(formatter)
logger.addHandler(ch)
logger.propagate = False


def _normalize_level(level: int | str) -> int:
    """Return a numeric logging level for *level*."""
    if isinstance(level, str):
        name = level.upper()
        try:
            return logging._nameToLevel[name]  # type: ignore[attr-defined]
        except KeyError as exc:
            raise ValueError(f"Unknown log level: {level}") from exc
    return int(level)


def _clear_handlers(existing: Iterable[logging.Handler]) -> None:
    for h in list(existing):
        try:
            h.close()
        finally:
            logger.removeHandler(h)


def configure_logging(
    *,
    log_level: int | str = logging.INFO,
    log_file: str | os.PathLike[str] | None = LOG_FILE,
    enable_file_logging: bool = True,
    enable_console: bool = True,
    create_dirs: bool = True,
    **legacy: object,  # accepts legacy 'verbose='
) -> None:
    """Idempotent logger setup for the VCF merger.

    Supports legacy ``verbose=bool`` (maps to enable_console).
    """
    if "verbose" in legacy:
        enable_console = bool(legacy["verbose"])

    level = _normalize_level(log_level)
    _clear_handlers(logger.handlers)
    logger.setLevel(level)

    fmt = logging.Formatter(LOG_FORMAT, datefmt=LOG_DATE_FORMAT)

    if enable_file_logging and log_file:
        path = os.fspath(log_file)
        if create_dirs:
            d = os.path.dirname(path)
            if d:
                os.makedirs(d, exist_ok=True)
        fh = logging.FileHandler(path)
        fh.setFormatter(fmt)
        logger.addHandler(fh)

    if enable_console:
        sh = logging.StreamHandler()
        sh.setFormatter(fmt)
        logger.addHandler(sh)


class MergeVCFError(RuntimeError):
    """Base exception for unrecoverable errors in the merge workflow."""


class ValidationError(MergeVCFError):
    """Raised when validation detects an unrecoverable problem."""


class MergeConflictError(MergeVCFError):
    """Raised when merging fails due to conflicting inputs or tooling errors."""


def log_message(
    message: str,
    verbose: bool = False,
    level: int = logging.INFO,
    *,
    exc_info: BaseException | bool | None = None,
) -> None:
    """Log *message* at the requested level and optionally echo it to stdout."""

    logger.log(level, message, exc_info=exc_info)
    if verbose:
        print(message)


def handle_critical_error(
    message: str,
    exc_cls=None,
    *,
    exc_info: BaseException | bool | None = None,
) -> None:
    """Log and raise a fatal error before exiting."""

    log_message(message, level=logging.ERROR)
    logger.critical(message, exc_info=exc_info)
    exception_class = exc_cls or MergeVCFError
    if isinstance(exc_info, BaseException):
        raise exception_class(message) from exc_info
    raise exception_class(message)


def handle_non_critical_error(message: str) -> None:
    """Log and print a recoverable error."""

    log_message(message, level=logging.WARNING)
    print(message)


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

# Default configuration: file + console at INFO level.
configure_logging()
