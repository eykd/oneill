"""Command-line adapter: parses arguments and writes to an injected stream."""

import argparse
import sys
from pathlib import Path
from typing import TYPE_CHECKING, TextIO

from oneill.adapters.plots import (
    PlottingUnavailableError,
    check_available,
    draw_figures,
    show_all,
)
from oneill.adapters.text import (
    render_ballistics,
    render_custom,
    render_elevator,
    render_lateral,
    render_selftest,
    render_stairs,
    render_summary,
    render_vertical,
)
from oneill.core.cylinder import G_STANDARD, Cylinder
from oneill.core.vectors import LocalVector
from oneill.motion.figures import build_figures
from oneill.motion.scenarios import (
    ballistics_report,
    custom_report,
    elevator_report,
    lateral_report,
    stairs_report,
    vertical_report,
)
from oneill.motion.verification import all_passed, run_selftest

if TYPE_CHECKING:
    from collections.abc import Sequence

    from oneill.motion.figures import FigureSpec

SECTIONS = ("lateral", "vertical", "stairs", "elevator", "ballistics")
_DEFAULT_HEIGHT_FRACTIONS = (0.0, 5 / 32, 0.5, 0.875)
_DEFAULT_ELEVATOR_TOP_FRACTION = 0.9375


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    ap = argparse.ArgumentParser(
        prog="oneill",
        description="Fictitious forces inside a rotating O'Neill cylinder",
    )
    ap.add_argument("--radius", type=float, default=3200.0, help="m (default 3200)")
    ap.add_argument("--length", type=float, default=20000.0, help="m (default 20000)")
    spin = ap.add_mutually_exclusive_group()
    spin.add_argument(
        "--rpm", type=float, help="spin rate; default = whatever gives --gravity"
    )
    spin.add_argument(
        "--gravity",
        type=float,
        default=G_STANDARD,
        help="floor gravity m/s^2 (default 1 g)",
    )
    ap.add_argument(
        "--spin-accel",
        type=float,
        default=0.0,
        help="rad/s^2, for spin-up/down (Euler force)",
    )
    ap.add_argument(
        "--mass", type=float, default=70.0, help="person mass in kg (default 70)"
    )
    ap.add_argument(
        "--speeds", type=float, nargs="+", default=[1.4, 3.0, 6.0, 12.0, 30.0]
    )
    ap.add_argument(
        "--headings",
        type=float,
        nargs="+",
        default=[0, 45, 90, 135, 180, 225, 270, 315],
    )
    ap.add_argument(
        "--climb-rates", type=float, nargs="+", default=[0.3, 1.0, 3.0, 10.0]
    )
    ap.add_argument(
        "--heights",
        type=float,
        nargs="+",
        default=None,
        help="heights for the vertical table "
        "(default: 0, 5/32, 1/2, 7/8 of the radius)",
    )
    ap.add_argument(
        "--stair-speed", type=float, default=0.7, help="m/s along the stair incline"
    )
    ap.add_argument("--stair-slope", type=float, default=35.0, help="deg")
    ap.add_argument(
        "--elevator-top", type=float, default=None, help="default: 15/16 of the radius"
    )
    ap.add_argument("--elevator-vmax", type=float, default=10.0)
    ap.add_argument("--elevator-amax", type=float, default=1.0)
    ap.add_argument(
        "--velocity",
        type=float,
        nargs=3,
        action="append",
        metavar=("SPIN", "AXIAL", "UP"),
        help="custom constant velocity to analyse (repeatable)",
    )
    ap.add_argument(
        "--at-height", type=float, default=0.0, help="height for --velocity cases"
    )
    ap.add_argument(
        "--sections",
        nargs="+",
        default=["all"],
        choices=["all", *SECTIONS],
    )
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--plot", action="store_true", help="show matplotlib figures")
    ap.add_argument("--save-plots", metavar="DIR", help="save figures as PNG into DIR")
    return ap


def _habitat(args: argparse.Namespace) -> Cylinder:
    if args.rpm is not None:
        return Cylinder.from_rpm(args.radius, args.length, args.rpm, args.spin_accel)
    return Cylinder.from_gravity(
        args.radius, args.length, args.gravity, args.spin_accel
    )


def _heights(args: argparse.Namespace, cyl: Cylinder) -> list[float]:
    if args.heights is not None:
        return list(args.heights)
    return [f * cyl.radius for f in _DEFAULT_HEIGHT_FRACTIONS]


def _elevator_top(args: argparse.Namespace, cyl: Cylinder) -> float:
    if args.elevator_top is not None:
        return float(args.elevator_top)
    return _DEFAULT_ELEVATOR_TOP_FRACTION * cyl.radius


def _report_text(args: argparse.Namespace, cyl: Cylinder) -> tuple[list[str], bool]:
    """Render everything requested and say whether all checks passed.

    Raises:
        ValueError: for impossible inputs.
    """
    chunks = [render_summary(cyl.summary())]
    if args.selftest:
        results = run_selftest(cyl)
        chunks.append(render_selftest(results))
        return chunks, all_passed(results)
    heights = _heights(args, cyl)
    top = _elevator_top(args, cyl)

    want = set(args.sections)
    everything = "all" in want
    if everything or "lateral" in want:
        chunks.append(
            render_lateral(lateral_report(cyl, args.speeds, args.headings, args.mass))
        )
    if everything or "vertical" in want:
        chunks.append(
            render_vertical(vertical_report(cyl, args.climb_rates, heights, args.mass))
        )
    if everything or "stairs" in want:
        stairs = stairs_report(
            cyl, args.stair_speed, args.stair_slope, args.headings, args.mass
        )
        chunks.append(render_stairs(stairs))
    if everything or "elevator" in want:
        ride = elevator_report(
            cyl, top, args.elevator_vmax, args.elevator_amax, args.mass
        )
        chunks.append(render_elevator(ride))
    if everything or "ballistics" in want:
        chunks.append(render_ballistics(ballistics_report(cyl)))
    if args.velocity:
        velocities = [LocalVector(*v) for v in args.velocity]
        custom = custom_report(cyl, velocities, args.at_height, args.mass)
        chunks.append(render_custom(custom))
    return chunks, True


def main(argv: Sequence[str] | None = None, stdout: TextIO | None = None) -> int:
    """Run the CLI and return the process exit status.

    Invalid inputs exit with status 2 via ``parser.error``; nothing is written to
    ``stdout`` in that case. A failing ``--selftest`` returns 1.
    """
    out = sys.stdout if stdout is None else stdout
    parser = build_parser()
    args = parser.parse_args(argv)
    plotting = (args.plot or bool(args.save_plots)) and not args.selftest
    try:
        if plotting:
            check_available()
        cyl = _habitat(args)
        chunks, ok = _report_text(args, cyl)
        figures: tuple[FigureSpec, ...] = ()
        if plotting:
            figures = build_figures(
                cyl,
                args.speeds,
                args.climb_rates,
                _elevator_top(args, cyl),
                args.elevator_vmax,
                args.elevator_amax,
            )
    except (ValueError, PlottingUnavailableError) as exc:
        parser.error(str(exc))
    out.write("".join(chunks))
    if plotting:
        save_dir = Path(args.save_plots) if args.save_plots else None
        draw_figures(figures, save_dir=save_dir, show=show_all if args.plot else None)
        if save_dir is not None:
            out.write(f"\nplots saved to {args.save_plots}/\n")
    return 0 if ok else 1


def run() -> None:
    """Console-script entry point."""
    sys.exit(main())
