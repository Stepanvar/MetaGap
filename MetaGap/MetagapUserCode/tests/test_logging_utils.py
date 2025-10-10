"""Tests for the reusable logging helpers."""

from __future__ import annotations

import logging
import importlib.util
import sys
from pathlib import Path


def _reload_logging_utils():
    # Reload the module to reset globals and default configuration between tests.
    module_path = Path(__file__).resolve().parents[1] / "merge_vcf" / "logging_utils.py"
    module_name = "metagap_logging_utils_test"
    sys.modules.pop(module_name, None)
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    default_log = Path(module.LOG_FILE)
    if default_log.exists():
        default_log.unlink()
    return module


def test_configure_logging_custom_file(tmp_path):
    module = _reload_logging_utils()
    log_path = tmp_path / "custom.log"

    module.configure_logging(log_level="WARNING", log_file=log_path)
    module.logger.warning("custom destination works")

    for handler in module.logger.handlers:
        handler.flush()

    contents = log_path.read_text(encoding="utf-8")
    assert "custom destination works" in contents


def test_configure_logging_disable_file_then_enable(tmp_path):
    module = _reload_logging_utils()

    module.configure_logging(log_level=logging.ERROR, enable_file_logging=False)
    module.logger.error("only the stream handler should process this")

    assert not any(
        isinstance(handler, logging.FileHandler) for handler in module.logger.handlers
    )

    custom_path = tmp_path / "final.log"
    module.configure_logging(log_level="INFO", log_file=custom_path)
    module.logger.info("reconfigured logging writes to file")

    file_handlers = [
        handler for handler in module.logger.handlers if isinstance(handler, logging.FileHandler)
    ]
    assert len(file_handlers) == 1

    for handler in file_handlers:
        handler.flush()

    assert "reconfigured logging writes to file" in custom_path.read_text(encoding="utf-8")
