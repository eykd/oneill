"""``oneill atmosphere``: regenerate the atmosphere tables and figures, run L3."""

import argparse
import json
import time
from pathlib import Path
from typing import TYPE_CHECKING, Any, TextIO

from oneill.adapters.atmo_text import BLOCK_NAMES, render_blocks, replace_blocks
from oneill.adapters.plots import (
    PlottingUnavailableError,
    check_available,
    draw_figures,
)
from oneill.atmo.figures import atmosphere_figures
from oneill.atmo.gas import ONEILL_1974
from oneill.atmo.l3_runs import (
    CASES,
    PRODUCTION_N_R,
    PRODUCTION_STRETCH,
    case_figures,
    run_case,
    summarize,
)
from oneill.atmo.report import FLOOR_TEMPERATURE, build_report

if TYPE_CHECKING:
    from collections.abc import Sequence

DEFAULT_L3_SUMMARY = Path("docs/data/l3_summary.json")
_SCHEMA = 1


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
    ap.add_argument(
        "--run-l3",
        metavar="OUT_DIR",
        help="run the L3 annulus cases, save their figures in OUT_DIR and merge "
        "their summaries into the JSON file (long: use --l3-case and --l3-days)",
    )
    ap.add_argument(
        "--l3-case",
        action="append",
        choices=sorted(CASES),
        help="run only this case (repeatable); default: every case",
    )
    ap.add_argument(
        "--l3-days", type=float, default=3.0, help="simulated days per case"
    )
    ap.add_argument(
        "--l3-summary",
        metavar="PATH",
        help=f"L3 summary JSON (default {DEFAULT_L3_SUMMARY})",
    )
    ap.add_argument("--l3-n-r", type=int, default=PRODUCTION_N_R, help="radial cells")
    ap.add_argument(
        "--l3-stretch",
        type=float,
        default=PRODUCTION_STRETCH,
        help="radial clustering toward the floor",
    )
    ap.add_argument("--l3-n-theta", type=int, default=256, help="azimuthal points")
    ap.add_argument(
        "--l3-sample-every", type=float, default=1800.0, help="diagnostic interval, s"
    )
    return ap


def load_l3_summary(path: Path) -> dict[str, Any] | None:
    """The committed L3 summary, or ``None`` if the file does not exist.

    Raises:
        ValueError: if the file is not a valid summary.
    """
    if not path.exists():
        return None
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict) or data.get("schema") != _SCHEMA:
        msg = f"{path} is not an L3 summary (schema {_SCHEMA})"
        raise ValueError(msg)
    return data


def merge_l3_summary(
    existing: dict[str, Any] | None, summaries: list[dict[str, Any]]
) -> dict[str, Any]:
    """``existing`` with ``summaries`` added or replaced by case name."""
    cases = dict(existing["cases"]) if existing else {}
    for summary in summaries:
        cases[summary["name"]] = summary
    return {"schema": _SCHEMA, "cases": cases}


def _write_json(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def _run_l3(args: argparse.Namespace, stdout: TextIO) -> None:
    """Run the requested cases, save figures and merge their summaries."""
    out = Path(args.run_l3)
    names = args.l3_case or list(CASES)
    summary_path = Path(args.l3_summary) if args.l3_summary else DEFAULT_L3_SUMMARY
    summaries = []
    for name in names:
        start = time.perf_counter()
        result = run_case(
            CASES[name],
            args.l3_days,
            n_r=args.l3_n_r,
            n_theta=args.l3_n_theta,
            stretch=args.l3_stretch,
            sample_every=args.l3_sample_every,
        )
        wall = time.perf_counter() - start
        summary = summarize(result)
        summary["wall_seconds"] = round(wall, 1)
        summaries.append(summary)
        saved = draw_figures(case_figures(result), save_dir=out)
        stdout.write(
            f"{name}: {summary['days']:g} days, {summary['steps']} steps, "
            f"{wall / 60.0:.1f} min wall, stationary={summary['stationary']}, "
            f"{len(saved)} figures in {out}/\n"
        )
    merged = merge_l3_summary(load_l3_summary(summary_path), summaries)
    _write_json(summary_path, merged)
    stdout.write(f"summary written to {summary_path}\n")


def main(argv: Sequence[str], stdout: TextIO) -> int:
    """Run ``oneill atmosphere``; invalid input exits with status 2."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.figures or args.run_l3:
            check_available()
        if args.run_l3:
            _run_l3(args, stdout)
            return 0
        report = build_report()
        doc_path = Path(args.write_doc) if args.write_doc else None
        summary_path = (
            Path(args.l3_summary)
            if args.l3_summary
            else (doc_path.parent / "data" / "l3_summary.json" if doc_path else None)
        )
        if summary_path is None or not summary_path.exists():
            summary_path = DEFAULT_L3_SUMMARY
        blocks = render_blocks(report, load_l3_summary(summary_path))
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
