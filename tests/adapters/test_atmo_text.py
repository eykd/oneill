from dataclasses import replace
from pathlib import Path

import pytest

from oneill.adapters.atmo_text import (
    BLOCK_NAMES,
    render_blocks,
    render_equilibrium,
    render_parcel,
    render_spinup,
    render_surface,
    replace_blocks,
)
from oneill.atmo.report import AtmosphereReport, build_report

DOC = Path(__file__).resolve().parents[2] / "docs" / "atmosphere.md"
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
        assert "canopy drag, C_D = 2e-3" in blocks["spinup"]

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
        assert "22.7 h" in text
        assert " d |" in text


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
    assert replace_blocks(committed, render_blocks(report)) == committed, (
        "docs/atmosphere.md is stale: run `uv run oneill atmosphere "
        "--write-doc docs/atmosphere.md`"
    )
