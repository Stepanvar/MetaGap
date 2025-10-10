"""Tests for the reusable logging helpers."""

from __future__ import annotations

import logging
import importlib.util
import sys
from pathlib import Path

import pytest


def _reload_logging_utils():
    """
    Reload the logging_utils module in isolation to reset globals and configuration
    between tests. Also deletes the default log file if the module defines LOG_FILE.
    """
    module_path = Path(__file__).resolve().parents[1] / "merge_vcf" / "logging_utils.py"
    module_name = "metagap_logging_utils_test"
    sys.modules.pop(module_name, None)
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)  # type: ignore[assignment]

    # Remove default log file if present to avoid cross-test interference
    log_file_path = getattr(module, "LOG_FILE", None)
    if log_file_path:
        default_log = Path(log_file_path)
        if default_log.exists():
            default_log.unlink()

    return module


# ------------------------
# Tests for error handlers
# ------------------------

def test_handle_non_critical_error_logs_warning(caplog):
    mod = _reload_logging_utils()
    message = "recoverable condition"

    with caplog.at_level(logging.DEBUG, logger=mod.logger.name):
        mod.handle_non_critical_error(message)

    warnings = [rec for rec in caplog.records if rec.levelno == logging.WARNING]
    assert warnings, "Expected a warning log entry"
    assert any(rec.message == message for rec in warnings)


def test_handle_critical_error_logs_error_and_critical(caplog):
    mod = _reload_logging_utils()
    message = "fatal condition detected"
    root_exc = ValueError("boom")

    with caplog.at_level(logging.DEBUG, logger=mod.logger.name):
        with pytest.raises(mod.MergeVCFError) as excinfo:
            mod.handle_critical_error(message, exc_info=root_exc)

    # Raised error should chain the original exception
    assert excinfo.value.__cause__ is root_exc

    # The message should be logged exactly once at ERROR and once at CRITICAL
    error_levels = [
        rec.levelno
        for rec in caplog.records
        if rec.message == message and rec.levelno in {logging.ERROR, logging.CRITICAL}
    ]
    assert error_levels.count(logging.ERROR) == 1
    assert error_levels.count(logging.CRITICAL) == 1


# ------------------------
# Tests for configuration
# ------------------------

def test_configure_logging_custom_file(tmp_path):
    mod = _reload_logging_utils()
    log_path = tmp_path / "custom.log"

    mod.configure_logging(log_level="WARNING", log_file=log_path)
    mod.logger.warning("custom destination works")

    # Flush any handlers so content is written
    for handler in mod.logger.handlers:
        try:
            handler.flush()
        except Exception:
            pass

    contents = log_path.read_text(encoding="utf-8")
    assert "custom destination works" in contents


def test_configure_logging_disable_file_then_enable(tmp_path):
    mod = _reload_logging_utils()

    # Disable file logging and emit an error (should not create FileHandler)
    mod.configure_logging(log_level=logging.ERROR, enable_file_logging=False)
    mod.logger.error("only the stream handler should process this")

    assert not any(isinstance(h, logging.FileHandler) for h in mod.logger.handlers)

    # Re-enable with a custom file and lower level
    custom_path = tmp_path / "final.log"
    mod.configure_logging(log_level="INFO", log_file=custom_path)
    mod.logger.info("reconfigured logging writes to file")

    file_handlers = [h for h in mod.logger.handlers if isinstance(h, logging.FileHandler)]
    assert len(file_handlers) == 1

    for handler in file_handlers:
        try:
            handler.flush()
        except Exception:
            pass

    assert "reconfigured logging writes to file" in custom_path.read_text(encoding="utf-8")
