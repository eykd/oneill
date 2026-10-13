"""``oneill atmosphere``: regenerate the atmosphere tables and figures."""

import argparse
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from oneill.adapters.atmo_text import BLOCK_NAMES, render_blocks, replace_blocks
from oneill.adapters.plots import (
    PlottingUnavailableError,
    check_available,
    draw_figures,
)
from oneill.atmo.figures import atmosphere_figures
from oneill.atmo.gas import ONEILL_1974
from oneill.atmo.report import FLOOR_TEMPERATURE, build_report

if TYPE_CHECKING:
    from collections.abc import Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the ``atmosphere`` argument parser."""
    ap = argparse.ArgumentParser(
        prog="oneill atmosphere",
        description=(
            "Compute the atmosphere tables for the default habitat (3.2 km, 20 km, "
            "1 g, both gas presets, 290 K). Without --write-doc the generated "
            "blocks are printed."
        ),
    )
    ap.add_argument(
        "--write-doc",
        metavar="PATH",
        help=(
            "rewrite the text between <!-- gen:NAME --> and <!-- /gen --> markers "
            f"in PATH (names: {', '.join(BLOCK_NAMES)}); other text is untouched"
        ),
    )
    ap.add_argument(
        "--figures", metavar="DIR", help="save the atmosphere figures as PNG into DIR"
    )
    return ap


def main(argv: Sequence[str], stdout: TextIO) -> int:
    """Run ``oneill atmosphere``; invalid input exits with status 2."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.figures:
            check_available()
        report = build_report()
        blocks = render_blocks(report)
        doc_path = Path(args.write_doc) if args.write_doc else None
        updated = None
        if doc_path is not None:
            current = doc_path.read_text(encoding="utf-8")
            updated = replace_blocks(current, blocks)
    except (ValueError, OSError, PlottingUnavailableError) as exc:
        parser.error(str(exc))

    if doc_path is None:
        stdout.write(
            "\n\n".join(
                f"<!-- gen:{name} -->\n{text}\n<!-- /gen -->"
                for name, text in blocks.items()
            )
            + "\n"
        )
    else:
        changed = updated != doc_path.read_text(encoding="utf-8")
        if changed and updated is not None:
            doc_path.write_text(updated, encoding="utf-8")
        stdout.write(f"{doc_path}: {'updated' if changed else 'already up to date'}\n")
    if args.figures:
        gas = ONEILL_1974.at_temperature(FLOOR_TEMPERATURE)
        saved = draw_figures(
            atmosphere_figures(report, gas), save_dir=Path(args.figures)
        )
        stdout.write(f"figures saved to {args.figures}/ ({len(saved)} files)\n")
    return 0
