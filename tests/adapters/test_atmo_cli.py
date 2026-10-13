import json
import shutil
import sys
from io import StringIO
from pathlib import Path

import pytest

from oneill.adapters import atmo_cli
from oneill.adapters.atmo_text import BLOCK_NAMES
from oneill.adapters.cli import main
from oneill.adapters.plots import PlottingUnavailableError

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


class TestRunL3:
    def tiny(self, tmp_path: Path, *extra: str) -> tuple[int, str, Path]:
        summary = tmp_path / "data" / "l3.json"
        status, out = run_cli(
            "--run-l3",
            str(tmp_path / "figs"),
            "--l3-summary",
            str(summary),
            "--l3-days",
            "0.01",
            "--l3-n-r",
            "16",
            "--l3-n-theta",
            "32",
            "--l3-sample-every",
            "200",
            *extra,
        )
        return status, out, summary

    def test_runs_a_case_writes_figures_and_merges_the_summary(
        self, tmp_path: Path
    ) -> None:
        status, out, summary = self.tiny(tmp_path, "--l3-case", "primary")
        assert status == 0
        assert "primary: 0.01 days" in out
        assert "min wall" in out
        names = sorted(p.name for p in (tmp_path / "figs").iterdir())
        assert names == [
            "l3_primary_meanwind.png",
            "l3_primary_snapshot.png",
            "l3_primary_spectrum.png",
        ]
        assert all((tmp_path / "figs" / n).stat().st_size > 3000 for n in names)
        data = json.loads(summary.read_text(encoding="utf-8"))
        assert data["schema"] == 1
        assert set(data["cases"]) == {"primary"}
        assert data["cases"]["primary"]["wall_seconds"] >= 0.0
        # a second case merges and keeps the first; re-running replaces
        before = data["cases"]["primary"]
        _, _, _ = self.tiny(tmp_path, "--l3-case", "cb_off")
        merged = json.loads(summary.read_text(encoding="utf-8"))
        assert set(merged["cases"]) == {"primary", "cb_off"}
        assert merged["cases"]["primary"] == before
        assert merged["cases"]["cb_off"]["coriolis_buoyancy"] is False

    def test_the_doc_block_renders_from_the_json(self, tmp_path: Path) -> None:
        _, _, summary = self.tiny(tmp_path, "--l3-case", "steady")
        doc = tmp_path / "doc.md"
        doc.write_text(SKELETON, encoding="utf-8")
        shutil.copy(summary, tmp_path / "data" / "l3_summary.json")
        status, _ = run_cli("--write-doc", str(doc))
        assert status == 0
        text = doc.read_text(encoding="utf-8")
        assert "| steady |" in text
        assert "No L3 runs are committed yet" not in text
        # explicit path wins over the default next to the doc
        empty = tmp_path / "empty.json"
        empty.write_text('{"schema": 1, "cases": {}}\n', encoding="utf-8")
        run_cli("--write-doc", str(doc), "--l3-summary", str(empty))
        assert "No L3 runs are committed yet" in doc.read_text(encoding="utf-8")

    def test_invalid_summary_file_exits_with_status_two(self, tmp_path: Path) -> None:
        bad = tmp_path / "bad.json"
        bad.write_text('{"nope": 1}\n', encoding="utf-8")
        doc = tmp_path / "doc.md"
        doc.write_text(SKELETON, encoding="utf-8")
        with pytest.raises(SystemExit) as exc:
            run_cli("--write-doc", str(doc), "--l3-summary", str(bad))
        assert exc.value.code == 2

    def test_run_l3_without_matplotlib_exits_with_status_two(
        self, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
    ) -> None:
        def unavailable() -> None:
            raise PlottingUnavailableError

        monkeypatch.setattr(atmo_cli, "check_available", unavailable)
        with pytest.raises(SystemExit) as exc:
            run_cli("--run-l3", str(tmp_path))
        assert exc.value.code == 2

    def test_default_summary_path_is_the_committed_one(self) -> None:
        assert atmo_cli.DEFAULT_L3_SUMMARY.as_posix() == "docs/data/l3_summary.json"
        assert atmo_cli.load_l3_summary(Path("/nonexistent/l3.json")) is None
        merged = atmo_cli.merge_l3_summary(None, [{"name": "a"}])
        assert merged == {"schema": 1, "cases": {"a": {"name": "a"}}}
