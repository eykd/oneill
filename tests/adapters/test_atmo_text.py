from dataclasses import replace
from pathlib import Path

import pytest

from oneill.adapters.atmo_cli import load_l3_summary
from oneill.adapters.atmo_text import (
    BLOCK_NAMES,
    _band,
    render_annulus,
    render_blocks,
    render_equilibrium,
    render_parcel,
    render_spinup,
    render_surface,
    replace_blocks,
)
from oneill.atmo.report import AtmosphereReport, build_report

DOC = Path(__file__).resolve().parents[2] / "docs" / "atmosphere.md"
L3_JSON = DOC.parent / "data" / "l3_summary.json"
BLOCKS = {name: f"content of {name}" for name in BLOCK_NAMES}


@pytest.fixture(scope="module")
def report() -> AtmosphereReport:
    return build_report()


class TestRendering:
    def test_every_block_is_markdown_with_a_caveat(
        self, report: AtmosphereReport
    ) -> None:
        blocks = render_blocks(report)
        assert tuple(blocks) == BLOCK_NAMES
        for name, text in blocks.items():
            if name != "annulus":  # no L3 summary passed: just the pointer
                assert text.count("| --- |") >= 1
            assert report.caveats[name] in text

    def test_key_numbers_appear(self, report: AtmosphereReport) -> None:
        blocks = render_blocks(report)
        assert "O'Neill 1974 at 290.0 K" in blocks["gas"]
        assert "SP-413 at 290.0 K" in blocks["gas"]
        assert "1828" in blocks["equilibrium"]
        assert "759x" in blocks["numbers"]
        assert "4.6 to 13.0" in blocks["surface"]
        assert "4.9 h" in blocks["parcel"]
        assert "canopy drag, C_D = 0.0178" in blocks["spinup"]
        assert "peak wind at 10 m" in blocks["spinup"]

    def test_missing_values_render_as_dashes_and_words(
        self, report: AtmosphereReport
    ) -> None:
        parcel = replace(report.parcels[0], deflection=None, drift_speed=None)
        text = render_parcel(replace(report, parcels=(parcel,)))
        assert "| - | - |" in text
        spin = replace(report.spin[0], time_to_90=None)
        assert "not reached" in render_spinup(replace(report, spin=(spin,)))
        surface = replace(report.surface[0], layer=None)
        assert "| none |" in render_surface(replace(report, surface=(surface,)))
        eq = report.equilibrium[0]
        cloud = replace(eq.clouds[0], height=None)
        shifted = replace(eq, clouds=(cloud, eq.clouds[1]))
        assert "above the axis" in render_equilibrium(
            replace(report, equilibrium=(shifted,))
        )

    def test_durations_use_sensible_units(self, report: AtmosphereReport) -> None:
        text = render_spinup(report)
        assert "60 min" in text
        assert "23.0 h" in text
        assert " d |" in text
        assert "n/a (no end caps)" in text
        assert "band over blend width 0.15-0.6" in text
        assert "Ekman Reynolds" in text


SUMMARY = {
    "name": "primary",
    "description": "day/night, neutral windows",
    "days": 3.0,
    "nu": 1.0,
    "r_in_fraction": 0.1,
    "coriolis_buoyancy": True,
    "stationary": False,
    "rms_wind_10m": 2.345,
    "gust_10m": 7.1,
    "zonal_wind_10m": -0.512,
    "m_peak": 4,
    "length_peak_m": 4500.0,
    "arrest_low_m": 700.0,
    "arrest_high_m": 1500.0,
    "lock_fraction": 0.4,
    "large_scale_fraction": 0.2,
    "nusselt": 12.34,
    "wall_seconds": 5400.0,
}


class TestAnnulusBlock:
    def test_without_runs_it_points_at_the_command(self) -> None:
        text = render_annulus(None, "caveat text")
        assert "--run-l3" in text
        assert "*caveat text*" in text
        assert render_annulus({"cases": {}}, "c").startswith("No L3 runs")

    def test_table_of_cases(self) -> None:
        extra = {
            **SUMMARY,
            "name": "custom",
            "description": "extra",
            "wall_seconds": None,
        }
        text = render_annulus(
            {"cases": {"primary": SUMMARY, "custom": extra}}, "caveat text"
        )
        assert "| primary | 3 | 1 | on | 0.1 | no | 2.35 | 7.10 | -0.51 | 4 |" in text
        assert "4.50" in text
        assert "0.7-1.5" in text
        assert "90 min" in text
        assert "| custom |" in text
        assert text.index("| primary |") < text.index("| custom |")
        assert "- **primary**: day/night, neutral windows" in text
        assert "*caveat text*" in text

    def test_render_blocks_includes_the_summary(self, report: AtmosphereReport) -> None:
        blocks = render_blocks(report, {"cases": {"primary": SUMMARY}})
        assert "| primary |" in blocks["annulus"]
        assert report.caveats["annulus"] in blocks["annulus"]


class TestBand:
    def test_range_in_one_unit(self) -> None:
        assert _band([22.8 * 3600.0, 24.0 * 3600.0, 23.7 * 3600.0]) == "22.8-24.0 h"

    def test_mixed_units(self) -> None:
        assert _band([30.0 * 3600.0, 3.0 * 86400.0 * 1.2]) == "30.0 h - 4 d"

    def test_identical_values_and_missing_values(self) -> None:
        assert _band([3600.0 * 5, 3600.0 * 5.01]) == "5.0 h"
        assert _band([None, 7200.0 * 6, None]) == "12.0 h"
        assert _band([None, None]) == "not reached"


class TestReplaceBlocks:
    doc = (
        "# Title\n\nprose before\n\n<!-- gen:gas -->\nold\n<!-- /gen -->\n\n"
        "middle prose with `code`\n\n<!-- gen:parcel -->\n<!-- /gen -->\n\nend\n"
    )

    def test_replaces_only_the_marked_content(self) -> None:
        out = replace_blocks(self.doc, BLOCKS)
        assert "content of gas" in out
        assert "content of parcel" in out
        assert "old" not in out
        assert out.startswith("# Title\n\nprose before\n\n<!-- gen:gas -->\n")
        assert (
            "\n<!-- /gen -->\n\nmiddle prose with `code`\n\n<!-- gen:parcel -->" in out
        )
        assert out.endswith("<!-- /gen -->\n\nend\n")

    def test_is_idempotent(self) -> None:
        once = replace_blocks(self.doc, BLOCKS)
        assert replace_blocks(once, BLOCKS) == once
        assert replace_blocks(once, BLOCKS).encode() == once.encode()

    def test_document_without_markers_is_untouched(self) -> None:
        assert replace_blocks("just prose\n", BLOCKS) == "just prose\n"

    def test_trailing_newlines_in_blocks_are_normalised(self) -> None:
        out = replace_blocks(self.doc, {**BLOCKS, "gas": "x\n\n"})
        assert "<!-- gen:gas -->\nx\n<!-- /gen -->" in out

    @pytest.mark.parametrize(
        ("doc", "message"),
        [
            ("<!-- gen:nope -->\n<!-- /gen -->", "unknown block 'nope'"),
            ("<!-- gen:gas -->\nno end", "never closed"),
            ("<!-- /gen -->", "without an opening"),
            (
                "<!-- gen:gas -->\n<!-- gen:parcel -->\n<!-- /gen -->",
                "not closed before",
            ),
        ],
    )
    def test_malformed_markers_are_errors(self, doc: str, message: str) -> None:
        with pytest.raises(ValueError, match=message):
            replace_blocks(doc, BLOCKS)


def test_committed_doc_is_up_to_date(report: AtmosphereReport) -> None:
    committed = DOC.read_text(encoding="utf-8")
    l3 = load_l3_summary(L3_JSON)
    assert replace_blocks(committed, render_blocks(report, l3)) == committed, (
        "docs/atmosphere.md is stale: run `uv run oneill atmosphere "
        "--write-doc docs/atmosphere.md`"
    )
