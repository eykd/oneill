"""Tests for the command-line adapter."""

import runpy
from io import StringIO

import pytest

from oneill import __version__
from oneill.adapters.cli import main, run


def test_main_reports_version_to_the_injected_stream() -> None:
    out = StringIO()
    assert main([], out) == 0
    assert out.getvalue() == f"oneill {__version__}\n"


def test_main_defaults_to_stdout(capsys: pytest.CaptureFixture[str]) -> None:
    assert main([]) == 0
    assert capsys.readouterr().out == f"oneill {__version__}\n"


def test_unknown_option_exits_with_usage_error() -> None:
    with pytest.raises(SystemExit) as exc:
        main(["--bogus"], StringIO())
    assert exc.value.code == 2


def test_run_exits_with_main_status(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.argv", ["oneill"])
    with pytest.raises(SystemExit) as exc:
        run()
    assert exc.value.code == 0


def test_python_dash_m_entry_point(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.argv", ["oneill"])
    with pytest.raises(SystemExit) as exc:
        runpy.run_module("oneill", run_name="__main__")
    assert exc.value.code == 0
