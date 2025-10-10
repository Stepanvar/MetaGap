"""Tests for logging configuration helpers."""

from __future__ import annotations

import importlib.util
import logging
import sys
from functools import lru_cache
from pathlib import Path


@lru_cache
def _load_logging_utils():
    base_dir = Path(__file__).resolve().parents[1]
    module_path = base_dir / "merge_vcf" / "logging_utils.py"
    spec = importlib.util.spec_from_file_location("test_logging_utils", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    sys.modules.setdefault(spec.name, module)
    spec.loader.exec_module(module)
    return module


def test_configure_logging_controls_console_output(tmp_path, capsys):
    logging_utils = _load_logging_utils()
    log_file = tmp_path / "run.log"

    try:
        logging_utils.configure_logging(verbose=False, log_file=str(log_file))
        logging_utils.log_message("quiet message")
        captured = capsys.readouterr()
        assert captured.out == ""
        assert captured.err == ""

        logging_utils.configure_logging(verbose=True, log_file=str(log_file))
        logging_utils.log_message("loud message")
        captured = capsys.readouterr()
        assert captured.out == ""
        assert "loud message" in captured.err

        logging_utils.configure_logging(verbose=False, log_file=str(log_file))
        logging_utils.log_message("quiet again")
        captured = capsys.readouterr()
        assert captured.out == ""
        assert captured.err == ""
    finally:
        logging_utils.configure_logging(verbose=False)


def test_configure_logging_replaces_handlers(tmp_path):
    logging_utils = _load_logging_utils()
    log_file = tmp_path / "handlers.log"

    try:
        logging_utils.configure_logging(verbose=True, log_file=str(log_file))
        file_handlers = [
            handler
            for handler in logging_utils.logger.handlers
            if isinstance(handler, logging.FileHandler)
        ]
        stream_handlers = [
            handler
            for handler in logging_utils.logger.handlers
            if isinstance(handler, logging.StreamHandler)
            and not isinstance(handler, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert len(stream_handlers) == 1

        logging_utils.configure_logging(verbose=True, log_file=str(log_file))
        file_handlers = [
            handler
            for handler in logging_utils.logger.handlers
            if isinstance(handler, logging.FileHandler)
        ]
        stream_handlers = [
            handler
            for handler in logging_utils.logger.handlers
            if isinstance(handler, logging.StreamHandler)
            and not isinstance(handler, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert len(stream_handlers) == 1

        logging_utils.configure_logging(verbose=False, log_file=str(log_file))
        file_handlers = [
            handler
            for handler in logging_utils.logger.handlers
            if isinstance(handler, logging.FileHandler)
        ]
        stream_handlers = [
            handler
            for handler in logging_utils.logger.handlers
            if isinstance(handler, logging.StreamHandler)
            and not isinstance(handler, logging.FileHandler)
        ]
        assert len(file_handlers) == 1
        assert len(stream_handlers) == 0
    finally:
        logging_utils.configure_logging(verbose=False)
