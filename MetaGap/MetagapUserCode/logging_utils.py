"""Shared logging utilities for MetaGap scripts."""

from __future__ import annotations

import logging
import sys

LOG_FILE = "script_execution.log"


def _configure_logger() -> logging.Logger:
    """Configure and return the shared logger instance."""

    logger = logging.getLogger("vcf_merger")
    if not logger.handlers:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s : %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )

        file_handler = logging.FileHandler(LOG_FILE)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        stream_handler = logging.StreamHandler()
        stream_handler.setFormatter(formatter)
        logger.addHandler(stream_handler)

    return logger


logger = _configure_logger()


def log_message(message: str, verbose: bool = False) -> None:
    """Log a message and optionally echo it to stdout."""

    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message: str) -> None:
    """Log a critical error and terminate the program."""

    log_message("CRITICAL ERROR: " + message)
    print(f"A critical error occurred. Check {LOG_FILE} for details.")
    sys.exit(1)


def handle_non_critical_error(message: str) -> None:
    """Log a non-critical error and emit a warning to stdout."""

    log_message("WARNING: " + message)
    print("Warning: " + message)


__all__ = [
    "LOG_FILE",
    "logger",
    "log_message",
    "handle_critical_error",
    "handle_non_critical_error",
]
