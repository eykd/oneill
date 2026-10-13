import sys
from io import StringIO
from typing import TYPE_CHECKING

import pytest

from oneill.adapters import atmo_cli
from oneill.adapters.atmo_text import BLOCK_NAMES
from oneill.adapters.cli import main
from oneill.adapters.plots import PlottingUnavailableError

if TYPE_CHECKING:
    from pathlib import Path

SKELETON = (
    "intro prose\n\n"
    + "\n\n".join(f"<!-- gen:{name} -->\n<!-- /gen -->" for name in BLOCK_NAMES)
    + "\n\nclosing prose\n"
)


def run_cli(*argv: str) -> tuple[int, str]:
    out = StringIO()
    return main(["atmosphere", *argv], out), out.getvalue()


def test_prints_every_block_with_markers() -> None:
    status, out = run_cli()
    assert status == 0
    for name in BLOCK_NAMES:
        assert f"<!-- gen:{name} -->" in out
    assert out.count("<!-- /gen -->") == len(BLOCK_NAMES)


def test_write_doc_is_idempotent_and_keeps_prose(tmp_path: Path) -> None:
    doc = tmp_path / "doc.md"
    doc.write_text(SKELETON, encoding="utf-8")
    status, out = run_cli("--write-doc", str(doc))
    assert status == 0
    assert "updated" in out
    first = doc.read_bytes()
    assert first.startswith(b"intro prose\n\n<!-- gen:gas -->\n| Preset")
    assert first.endswith(b"<!-- /gen -->\n\nclosing prose\n")
    status, out = run_cli("--write-doc", str(doc))
    assert status == 0
    assert "already up to date" in out
    assert doc.read_bytes() == first


def test_unknown_marker_exits_with_status_two(tmp_path: Path) -> None:
    doc = tmp_path / "doc.md"
    doc.write_text("<!-- gen:bogus -->\n<!-- /gen -->\n", encoding="utf-8")
    with pytest.raises(SystemExit) as exc:
        run_cli("--write-doc", str(doc))
    assert exc.value.code == 2
    assert doc.read_text(encoding="utf-8") == "<!-- gen:bogus -->\n<!-- /gen -->\n"


def test_missing_doc_exits_with_status_two(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as exc:
        run_cli("--write-doc", str(tmp_path / "missing.md"))
    assert exc.value.code == 2


def test_figures_are_saved(tmp_path: Path) -> None:
    status, out = run_cli("--figures", str(tmp_path / "figs"))
    assert status == 0
    assert "4 files" in out
    names = sorted(p.name for p in (tmp_path / "figs").iterdir())
    assert names == [
        "atmo_1_profiles.png",
        "atmo_2_wind.png",
        "atmo_3_parcel.png",
        "atmo_4_spinup.png",
    ]
    assert all((tmp_path / "figs" / n).stat().st_size > 5000 for n in names)


def test_figures_without_matplotlib_exit_with_status_two(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    def unavailable() -> None:
        raise PlottingUnavailableError

    monkeypatch.setattr(atmo_cli, "check_available", unavailable)
    with pytest.raises(SystemExit) as exc:
        run_cli("--figures", str(tmp_path))
    assert exc.value.code == 2


def test_subcommand_reads_sys_argv_when_not_given_argv(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(sys, "argv", ["oneill", "atmosphere"])
    out = StringIO()
    assert main(None, out) == 0
    assert "<!-- gen:gas -->" in out.getvalue()


def test_original_cli_is_unchanged() -> None:
    out = StringIO()
    assert main(["--sections", "lateral"], out) == 0
    assert "O'NEILL CYLINDER PARAMETERS" in out.getvalue()
