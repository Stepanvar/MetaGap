import importlib.util
import logging
from pathlib import Path

import pytest


def _load_logging_utils():
    module_path = Path(__file__).resolve().parents[1] / "merge_vcf" / "logging_utils.py"
    spec = importlib.util.spec_from_file_location("merge_vcf.logging_utils", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


logging_utils = _load_logging_utils()
MergeVCFError = logging_utils.MergeVCFError
handle_critical_error = logging_utils.handle_critical_error
handle_non_critical_error = logging_utils.handle_non_critical_error
logger = logging_utils.logger


def test_handle_non_critical_error_logs_warning(caplog):
    message = "recoverable condition"

    with caplog.at_level(logging.DEBUG, logger=logger.name):
        handle_non_critical_error(message)

    warnings = [
        record for record in caplog.records if record.levelno == logging.WARNING
    ]
    assert warnings, "Expected a warning log entry"
    assert any(record.message == message for record in warnings)


def test_handle_critical_error_logs_error_and_critical(caplog):
    message = "fatal condition detected"
    root_exc = ValueError("boom")

    with caplog.at_level(logging.DEBUG, logger=logger.name):
        with pytest.raises(MergeVCFError) as excinfo:
            handle_critical_error(message, exc_info=root_exc)

    assert excinfo.value.__cause__ is root_exc

    error_levels = [
        record.levelno
        for record in caplog.records
        if record.message == message and record.levelno in {logging.ERROR, logging.CRITICAL}
    ]
    assert error_levels.count(logging.ERROR) == 1
    assert error_levels.count(logging.CRITICAL) == 1
