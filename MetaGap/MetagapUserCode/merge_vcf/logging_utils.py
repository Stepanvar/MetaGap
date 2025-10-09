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


def log_message(message: str, verbose: bool = False) -> None:
    """Log *message* and optionally echo it to stdout."""

    logger.info(message)
    if verbose:
        print(message)


def handle_critical_error(message: str) -> None:
    """Log and raise a fatal error before exiting."""

    log_message("CRITICAL ERROR: " + message)
    print(f"A critical error occurred. Check {LOG_FILE} for details.")
    raise SystemExit(1)


def handle_non_critical_error(message: str) -> None:
    """Log and print a recoverable error."""

    log_message("WARNING: " + message)
    print("Warning: " + message)


__all__ = [
    "LOG_FILE",
    "logger",
    "log_message",
    "handle_critical_error",
    "handle_non_critical_error",
]
