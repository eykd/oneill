"""Figure specifications: what to plot, as plain data.

Nothing here draws. ``adapters.plots`` turns a ``FigureSpec`` into pixels.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.motion.ballistics import ballistic
from oneill.motion.forces import analyze_path
from oneill.motion.paths import FloorWalk, StraightPath
from oneill.motion.profiles import ConstantSpeed, Trapezoid
from oneill.motion.scenarios import ballistic_cases

if TYPE_CHECKING:
    from collections.abc import Sequence

    from oneill.core.cylinder import Cylinder
    from oneill.core.vectors import Vec3

_HEADING_STEP = 5.0
_CLIMB_SAMPLES = 250
_CLIMB_MAX_FRACTION = 0.99
_ELEVATOR_SAMPLES = 600
_HIGH_DROP = 100.0


@dataclass(frozen=True, eq=False)
class Series:
    """One curve; ``fmt`` is a matplotlib format string ("" for the default)."""

    label: str
    x: Vec3
    y: Vec3
    fmt: str = ""


@dataclass(frozen=True, eq=False)
class Panel:
    """One set of axes."""

    xlabel: str
    ylabel: str
    series: tuple[Series, ...]
    title: str = ""
    hlines: tuple[float, ...] = ()
    xticks: tuple[float, ...] = ()
    xticklabels: tuple[str, ...] = ()
    legend_size: float | None = None


@dataclass(frozen=True, eq=False)
class FigureSpec:
    """A figure: a grid of panels, a file name and a size in inches."""

    filename: str
    panels: tuple[Panel, ...]
    layout: tuple[int, int]
    figsize: tuple[float, float]
    suptitle: str = ""
    sharex: bool = False


def _floor_weight(cyl: Cylinder, speeds: Sequence[float]) -> FigureSpec:
    headings = np.arange(0.0, 360.0 + _HEADING_STEP, _HEADING_STEP)
    series = []
    for v in speeds:
        weights = [
            analyze_path(cyl, FloorWalk(cyl, hd, v), [0.0]).felt_weight[0]
            for hd in headings
        ]
        series.append(
            Series(
                f"{v:g} m/s", headings, 100.0 * np.array(weights) / cyl.surface_gravity
            )
        )
    panel = Panel(
        xlabel="",
        ylabel="felt weight [% floor g]",
        series=tuple(series),
        title="Level motion: Coriolis only changes your weight",
        hlines=(100.0,),
        xticks=(0.0, 90.0, 180.0, 270.0, 360.0),
        xticklabels=("+axial", "spinward", "-axial", "antispinward", "+axial"),
    )
    return FigureSpec("1_floor_weight_vs_heading.png", (panel,), (1, 1), (8.0, 5.0))


def _climbing(cyl: Cylinder, rates: Sequence[float]) -> FigureSpec:
    heights = np.linspace(0.0, _CLIMB_MAX_FRACTION * cyl.radius, _CLIMB_SAMPLES)
    gravity = 100.0 * cyl.omega**2 * (cyl.radius - heights) / cyl.surface_gravity
    push = [Series("local gravity (centrifugal)", heights, gravity, "k--")]
    tilt = []
    for rate in rates:
        hists = [
            analyze_path(
                cyl,
                StraightPath(cyl, 0.0, 90.0, ConstantSpeed(rate), 1.0, start_height=h),
                [0.0],
            )
            for h in heights
        ]
        push.append(
            Series(
                f"spinward push @ {rate:g} m/s up",
                heights,
                np.array(
                    [100.0 * x.lateral_mag[0] / cyl.surface_gravity for x in hists]
                ),
            )
        )
        tilt.append(
            Series(
                f"climbing {rate:g} m/s",
                heights,
                np.array([x.tilt_deg[0] for x in hists]),
            )
        )
    return FigureSpec(
        "2_climbing_vs_height.png",
        (
            Panel("height above floor [m]", "% floor g", tuple(push), legend_size=8.0),
            Panel(
                "height above floor [m]",
                "tilt of apparent vertical [deg]",
                tuple(tilt),
                legend_size=8.0,
            ),
        ),
        (1, 2),
        (12.0, 4.5),
        suptitle="Climbing: sideways Coriolis push vs. fading gravity",
    )


def _elevator(cyl: Cylinder, top: float, v_max: float, a_max: float) -> FigureSpec:
    profile = Trapezoid(top, v_max, a_max)
    times = np.linspace(0.0, profile.duration, _ELEVATOR_SAMPLES)
    heights: list[Series] = []
    forces: list[Series] = []
    tilts: list[Series] = []
    for label, elevation, h0 in (("up", 90.0, 0.0), ("down", -90.0, top)):
        path = StraightPath(
            cyl, 0.0, elevation, profile, profile.duration, start_height=h0
        )
        hist = analyze_path(cyl, path, times)
        heights.append(Series(label, times, hist.heights))
        forces.append(
            Series(f"{label}: felt weight", times, hist.pct(hist.felt_weight))
        )
        forces.append(
            Series(
                f"{label}: spinward push", times, hist.pct(hist.apparent[:, 0]), "--"
            )
        )
        tilts.append(Series(label, times, hist.tilt_deg))
    return FigureSpec(
        "3_elevator.png",
        (
            Panel("", "height [m]", tuple(heights), legend_size=8.0),
            Panel("", "% floor g", tuple(forces), legend_size=8.0),
            Panel("time [s]", "tilt [deg]", tuple(tilts), legend_size=8.0),
        ),
        (3, 1),
        (9.0, 8.0),
        suptitle="Elevator ride",
        sharex=True,
    )


def _ballistics(cyl: Cylinder) -> FigureSpec:
    low: list[Series] = []
    high: list[Series] = []
    for case in ballistic_cases(cyl):
        res = ballistic(cyl, case.start_height, case.v_local, case.label)
        target = high if case.start_height >= _HIGH_DROP else low
        target.append(Series(case.label, res.traj_spin, res.traj_height))
    return FigureSpec(
        "4_ballistics.png",
        (
            Panel("spinward distance [m]", "height [m]", tuple(low), legend_size=7.0),
            Panel("spinward distance [m]", "height [m]", tuple(high), legend_size=7.0),
        ),
        (1, 2),
        (12.0, 4.5),
        suptitle="Free-flight trajectories seen by habitat residents",
    )


def build_figures(
    cyl: Cylinder,
    speeds: Sequence[float],
    climb_rates: Sequence[float],
    elevator_top: float,
    elevator_vmax: float,
    elevator_amax: float,
) -> tuple[FigureSpec, ...]:
    """The four standard figures for ``cyl``."""
    return (
        _floor_weight(cyl, speeds),
        _climbing(cyl, climb_rates),
        _elevator(cyl, elevator_top, elevator_vmax, elevator_amax),
        _ballistics(cyl),
    )
