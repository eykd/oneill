"""Scenario builders: the questions the simulator answers, as plain data.

Each ``*_report`` function evaluates one family of situations and returns a
frozen value object. Nothing here formats or prints; see ``adapters.text``.
Vectors are ``LocalVector`` (spinward, axial, up); ``*_pct`` fields are percent
of the floor gravity, and ``*_n`` fields are newtons.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import finite
from oneill.core.vectors import LocalVector
from oneill.motion.ballistics import ballistic
from oneill.motion.forces import analyze_path
from oneill.motion.headings import heading_name
from oneill.motion.paths import FloorWalk, StraightPath, velocity_path
from oneill.motion.profiles import ConstantSpeed, Trapezoid

if TYPE_CHECKING:
    from collections.abc import Sequence

    from oneill.core.cylinder import Cylinder
    from oneill.motion.forces import ForceHistory

_ELEVATOR_NUDGE = 0.05
_ELEVATOR_SAMPLES = 2000


def _local(row: object) -> LocalVector:
    return LocalVector.from_array(row)  # type: ignore[arg-type]


def _coriolis_pct_per_mps(cyl: Cylinder) -> float:
    return 100.0 * 2.0 * cyl.omega / cyl.surface_gravity


# =============================================================================
# 1. Lateral motion on a level floor
# =============================================================================
@dataclass(frozen=True)
class LateralRow:
    """One walk at one speed and heading."""

    heading_deg: float
    direction: str
    speed: float
    coriolis_pct: LocalVector
    coriolis_force_n: LocalVector
    curvature_pct: float
    felt_weight_pct: float
    tilt_deg: float


@dataclass(frozen=True)
class LateralReport:
    """Level-floor motion, grouped by speed."""

    height: float
    mass: float
    surface_gravity: float
    coriolis_pct_per_mps: float
    groups: tuple[tuple[LateralRow, ...], ...]


def lateral_report(
    cyl: Cylinder,
    speeds: Sequence[float],
    headings: Sequence[float],
    mass: float,
    height: float = 0.0,
) -> LateralReport:
    """Weight and Coriolis force while walking on a level floor."""
    finite("mass", mass)
    groups = []
    for v in speeds:
        rows = []
        for hd in headings:
            hist = analyze_path(cyl, FloorWalk(cyl, hd, v, height=height), [0.0])
            rows.append(
                LateralRow(
                    heading_deg=hd,
                    direction=heading_name(hd),
                    speed=v,
                    coriolis_pct=_local(hist.pct(hist.coriolis[0])),
                    coriolis_force_n=_local(mass * hist.coriolis[0]),
                    curvature_pct=float(hist.pct(hist.acc[0, 2])),
                    felt_weight_pct=float(hist.pct(hist.felt_weight[0])),
                    tilt_deg=float(hist.tilt_deg[0]),
                )
            )
        groups.append(tuple(rows))
    return LateralReport(
        height=height,
        mass=mass,
        surface_gravity=cyl.surface_gravity,
        coriolis_pct_per_mps=_coriolis_pct_per_mps(cyl),
        groups=tuple(groups),
    )


# =============================================================================
# 2. Vertical motion
# =============================================================================
@dataclass(frozen=True)
class VerticalRow:
    """One constant-rate climb or descent."""

    height: float
    rate: float
    local_g_pct: float
    coriolis_pct: LocalVector
    coriolis_force_n: LocalVector
    felt_weight_pct: float
    tilt_deg: float
    push: str
    push_force_n: float


@dataclass(frozen=True)
class VerticalGroup:
    """All rates at one starting height; ``rows`` is empty if the height is invalid."""

    height: float
    rows: tuple[VerticalRow, ...]

    @property
    def skipped(self) -> bool:
        """Whether the height was outside ``0 <= height < radius``."""
        return not self.rows


@dataclass(frozen=True)
class VerticalReport:
    """Ladders, ropes and elevators at constant speed, grouped by height."""

    mass: float
    surface_gravity: float
    coriolis_pct_per_mps: float
    groups: tuple[VerticalGroup, ...]


def vertical_report(
    cyl: Cylinder,
    rates: Sequence[float],
    heights: Sequence[float],
    mass: float,
) -> VerticalReport:
    """Coriolis push and apparent-vertical tilt while climbing or descending."""
    finite("mass", mass)
    groups = []
    for hgt in heights:
        if not 0 <= hgt < cyl.radius:
            groups.append(VerticalGroup(hgt, ()))
            continue
        rows = []
        for rate in [*rates, *(-r for r in rates)]:
            path = StraightPath(
                cyl,
                0.0,
                90.0 if rate > 0 else -90.0,
                ConstantSpeed(abs(rate)),
                1.0,
                start_height=hgt,
            )
            hist = analyze_path(cyl, path, [0.0])
            rows.append(
                VerticalRow(
                    height=hgt,
                    rate=rate,
                    local_g_pct=float(hist.pct(-hist.centrifugal[0, 2])),
                    coriolis_pct=_local(hist.pct(hist.coriolis[0])),
                    coriolis_force_n=_local(mass * hist.coriolis[0]),
                    felt_weight_pct=float(hist.pct(hist.felt_weight[0])),
                    tilt_deg=float(hist.tilt_deg[0]),
                    push="spinward" if hist.apparent[0, 0] > 0 else "antispinward",
                    push_force_n=float(mass * hist.lateral_mag[0]),
                )
            )
        groups.append(VerticalGroup(hgt, tuple(rows)))
    return VerticalReport(
        mass=mass,
        surface_gravity=cyl.surface_gravity,
        coriolis_pct_per_mps=_coriolis_pct_per_mps(cyl),
        groups=tuple(groups),
    )


# =============================================================================
# 3. Stairs
# =============================================================================
@dataclass(frozen=True)
class StairsRow:
    """One heading on a straight staircase."""

    heading_deg: float
    direction: str
    sense: str
    velocity: LocalVector
    coriolis_pct: LocalVector
    felt_weight_pct: float
    lateral_pct: tuple[float, float]
    tilt_deg: float


@dataclass(frozen=True)
class StairsReport:
    """A straight staircase climbed (up) and descended (down) in each heading."""

    speed: float
    slope_deg: float
    height: float
    mass: float
    groups: tuple[tuple[StairsRow, ...], ...]


def stairs_report(
    cyl: Cylinder,
    speed: float,
    slope: float,
    headings: Sequence[float],
    mass: float,
    height: float = 3.0,
) -> StairsReport:
    """What you feel on stairs: both a weight change and a sideways push."""
    finite("mass", mass)
    groups = []
    for sense, elevation in (("up", slope), ("down", -slope)):
        rows = []
        for hd in headings:
            path = StraightPath(
                cyl, hd, elevation, ConstantSpeed(speed), 1.0, start_height=height
            )
            hist = analyze_path(cyl, path, [0.0])
            lateral = hist.pct(hist.apparent[0, :2])
            rows.append(
                StairsRow(
                    heading_deg=hd,
                    direction=heading_name(hd),
                    sense=sense,
                    velocity=_local(hist.vel[0]),
                    coriolis_pct=_local(hist.pct(hist.coriolis[0])),
                    felt_weight_pct=float(hist.pct(hist.felt_weight[0])),
                    lateral_pct=(float(lateral[0]), float(lateral[1])),
                    tilt_deg=float(hist.tilt_deg[0]),
                )
            )
        groups.append(tuple(rows))
    return StairsReport(
        speed=speed, slope_deg=slope, height=height, mass=mass, groups=tuple(groups)
    )


# =============================================================================
# 4. Elevator rides
# =============================================================================
@dataclass(frozen=True)
class ElevatorRow:
    """The rider's state at one instant of a ride."""

    t: float
    height: float
    v_up: float
    a_up: float
    local_g_pct: float
    felt_weight_pct: float
    spin_push_pct: float
    side_force_n: float
    tilt_deg: float


@dataclass(frozen=True)
class ElevatorRide:
    """A full ride in one direction, sampled and summarised."""

    label: str
    rows: tuple[ElevatorRow, ...]
    max_tilt_deg: float
    max_tilt_height: float
    max_side_force_n: float
    weight_min_pct: float
    weight_max_pct: float


@dataclass(frozen=True)
class ElevatorReport:
    """Up and down rides between the floor and ``top``."""

    top: float
    v_max: float
    a_max: float
    mass: float
    duration: float
    rides: tuple[ElevatorRide, ...]


def _elevator_ride(
    cyl: Cylinder, label: str, path: StraightPath, duration: float, mass: float, n: int
) -> ElevatorRide:
    # Acceleration jumps at the ends of a ride, so sample just inside them.
    nudge = min(_ELEVATOR_NUDGE, duration / 100)
    times = np.linspace(0.0, duration, n)
    times[0] += nudge
    times[-1] -= nudge
    hist = analyze_path(cyl, path, times)
    rows = tuple(
        ElevatorRow(
            t=float(times[i]),
            height=float(hist.heights[i]),
            v_up=float(hist.vel[i, 2]),
            a_up=float(hist.acc[i, 2]),
            local_g_pct=float(hist.pct(-hist.centrifugal[i, 2])),
            felt_weight_pct=float(hist.pct(hist.felt_weight[i])),
            spin_push_pct=float(hist.pct(hist.apparent[i, 0])),
            side_force_n=float(mass * hist.apparent[i, 0]),
            tilt_deg=float(hist.tilt_deg[i]),
        )
        for i in range(len(times))
    )
    fine = analyze_path(
        cyl, path, np.linspace(nudge, duration - nudge, _ELEVATOR_SAMPLES)
    )
    worst = int(np.argmax(fine.tilt_deg))
    weights = fine.pct(fine.felt_weight)
    return ElevatorRide(
        label=label,
        rows=rows,
        max_tilt_deg=float(fine.tilt_deg[worst]),
        max_tilt_height=float(fine.heights[worst]),
        max_side_force_n=float(mass * fine.lateral_mag.max()),
        weight_min_pct=float(weights.min()),
        weight_max_pct=float(weights.max()),
    )


def elevator_report(
    cyl: Cylinder,
    top: float,
    v_max: float,
    a_max: float,
    mass: float,
    rows: int = 15,
) -> ElevatorReport:
    """Ride between the floor and ``top`` metres up, accelerating and braking.

    Raises:
        ValueError: if ``top`` is not strictly between the floor and the axis.
    """
    finite("mass", mass)
    if not 0 < top < cyl.radius:
        msg = "elevator top must be between 0 and the radius"
        raise ValueError(msg)
    profile = Trapezoid(top, v_max, a_max)
    duration = profile.duration
    rides = tuple(
        _elevator_ride(
            cyl,
            label,
            StraightPath(cyl, 0.0, elevation, profile, duration, start_height=h0),
            duration,
            mass,
            rows,
        )
        for label, elevation, h0 in (("UP", 90.0, 0.0), ("DOWN", -90.0, top))
    )
    return ElevatorReport(top, v_max, a_max, mass, duration, rides)


# =============================================================================
# 5. Free flight
# =============================================================================
@dataclass(frozen=True)
class BallisticCase:
    """A named launch: height above the floor and velocity in local components."""

    label: str
    start_height: float
    v_local: LocalVector


_SQRT_HALF = math.sqrt(0.5)
BALLISTIC_CASES = (
    BallisticCase("drop from hand (1.5 m)", 1.5, LocalVector(0.0, 0.0, 0.0)),
    BallisticCase("drop from 10 m", 10.0, LocalVector(0.0, 0.0, 0.0)),
    BallisticCase("drop from 100 m tower", 100.0, LocalVector(0.0, 0.0, 0.0)),
    BallisticCase("drop from 1000 m", 1000.0, LocalVector(0.0, 0.0, 0.0)),
    BallisticCase("jump straight up 2.5 m/s", 0.0, LocalVector(0.0, 0.0, 2.5)),
    BallisticCase(
        "throw 15 m/s 45deg spinward",
        1.5,
        LocalVector(15 * _SQRT_HALF, 0.0, 15 * _SQRT_HALF),
    ),
    BallisticCase(
        "throw 15 m/s 45deg antispinward",
        1.5,
        LocalVector(-15 * _SQRT_HALF, 0.0, 15 * _SQRT_HALF),
    ),
    BallisticCase(
        "throw 15 m/s 45deg +axial",
        1.5,
        LocalVector(0.0, 15 * _SQRT_HALF, 15 * _SQRT_HALF),
    ),
)


def ballistic_cases(cyl: Cylinder) -> tuple[BallisticCase, ...]:
    """The default free-flight cases whose launch points fit inside ``cyl``."""
    return tuple(c for c in BALLISTIC_CASES if c.start_height < cyl.radius)


@dataclass(frozen=True)
class BallisticRow:
    """One flight: exact versus the flat-ground naive prediction."""

    label: str
    start_height: float
    flight_time: float
    naive_time: float
    landing: tuple[float, float]
    naive_landing: tuple[float, float]
    deviation: tuple[float, float]
    hits_endcap: bool


@dataclass(frozen=True)
class BallisticsReport:
    """Drops, jumps and throws."""

    rows: tuple[BallisticRow, ...]


def ballistics_report(
    cyl: Cylinder, cases: Sequence[BallisticCase] | None = None
) -> BallisticsReport:
    """Exact landing points for ``cases`` (default: the standard set)."""
    chosen = ballistic_cases(cyl) if cases is None else cases
    rows = []
    for case in chosen:
        res = ballistic(cyl, case.start_height, case.v_local, case.label)
        rows.append(
            BallisticRow(
                label=case.label,
                start_height=case.start_height,
                flight_time=res.flight_time,
                naive_time=res.naive_time,
                landing=(res.landing_spin, res.landing_axial),
                naive_landing=(res.naive_spin, res.naive_axial),
                deviation=(
                    res.landing_spin - res.naive_spin,
                    res.landing_axial - res.naive_axial,
                ),
                hits_endcap=res.hits_endcap,
            )
        )
    return BallisticsReport(tuple(rows))


# =============================================================================
# 6. Custom velocities
# =============================================================================
@dataclass(frozen=True)
class CustomEntry:
    """Full force breakdown for one constant velocity."""

    velocity: LocalVector
    centrifugal_pct: LocalVector
    coriolis_pct: LocalVector
    coriolis_force_n: LocalVector
    euler_pct: LocalVector
    neg_path_accel_pct: LocalVector
    apparent_pct: LocalVector
    felt_weight_pct: float
    tilt_deg: float
    support_force_n: LocalVector


@dataclass(frozen=True)
class CustomReport:
    """Force breakdowns for user-chosen velocities at one height."""

    height: float
    mass: float
    entries: tuple[CustomEntry, ...]


def custom_report(
    cyl: Cylinder,
    velocities: Sequence[LocalVector],
    height: float,
    mass: float,
) -> CustomReport:
    """Break down every term for each constant ``velocity`` at ``height``."""
    finite("mass", mass)
    entries = []
    for v in velocities:
        hist: ForceHistory = analyze_path(cyl, velocity_path(cyl, v, height), [0.0])
        entries.append(
            CustomEntry(
                velocity=v,
                centrifugal_pct=_local(hist.pct(hist.centrifugal[0])),
                coriolis_pct=_local(hist.pct(hist.coriolis[0])),
                coriolis_force_n=_local(mass * hist.coriolis[0]),
                euler_pct=_local(hist.pct(hist.euler[0])),
                neg_path_accel_pct=_local(hist.pct(-hist.acc[0])),
                apparent_pct=_local(hist.pct(hist.apparent[0])),
                felt_weight_pct=float(hist.pct(hist.felt_weight[0])),
                tilt_deg=float(hist.tilt_deg[0]),
                support_force_n=_local(hist.support_force(mass)[0]),
            )
        )
    return CustomReport(height, mass, tuple(entries))
