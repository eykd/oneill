"""Built-in physics self-test: independent closed forms against the force model."""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.vectors import Z_HAT, LocalVector
from oneill.motion.ballistics import ballistic, integrate_rotating
from oneill.motion.forces import analyze_path
from oneill.motion.paths import FloorWalk, StraightPath
from oneill.motion.profiles import ConstantSpeed

if TYPE_CHECKING:
    from collections.abc import Iterable

    from oneill.core.cylinder import Cylinder

_WEIGHT_TOL = 1e-9
_CORIOLIS_TOL = 1e-9
_RK4_POSITION_TOL = 1e-6
_RK4_MAX_STEP = 1e-2
_RK4_FLIGHT_FRACTION = 0.8


@dataclass(frozen=True)
class CheckResult:
    """One self-test comparison."""

    name: str
    got: float
    want: float
    tol: float

    @property
    def passed(self) -> bool:
        """Whether ``got`` is within ``tol`` of ``want``."""
        return abs(self.got - self.want) <= self.tol


def all_passed(results: Iterable[CheckResult]) -> bool:
    """Whether every check passed (vacuously true for none)."""
    return all(r.passed for r in results)


def run_selftest(cyl: Cylinder) -> tuple[CheckResult, ...]:
    """Check the force model and the ballistics solver against closed forms."""
    results: list[CheckResult] = []
    radius, omega = cyl.radius, cyl.omega

    for v in (-50.0, -1.4, 0.0, 1.4, 50.0):
        hist = analyze_path(cyl, FloorWalk(cyl, 90.0, v), [0.0])
        results.append(
            CheckResult(
                f"floor weight at {v:+g} m/s spinward == (wR+v)^2/R",
                float(hist.felt_weight[0]),
                (omega * radius + v) ** 2 / radius,
                _WEIGHT_TOL,
            )
        )
    for v in (0.5, 10.0):
        path = StraightPath(cyl, 0.0, 90.0, ConstantSpeed(v), 1.0)
        hist = analyze_path(cyl, path, [0.0])
        results.append(
            CheckResult(
                f"Coriolis climbing {v:g} m/s == 2wv spinward",
                float(hist.coriolis[0, 0]),
                2.0 * omega * v,
                _CORIOLIS_TOL,
            )
        )
    for label, h0, v_launch in (
        (f"drop {radius / 32:g} m", radius / 32, LocalVector(0.0, 0.0, 0.0)),
        ("throw", 1.5, LocalVector(10.0, 3.0, 10.0)),
        (f"drop {5 * radius / 16:g} m", 5 * radius / 16, LocalVector(0.0, 0.0, 0.0)),
    ):
        res = ballistic(cyl, h0, v_launch)
        t_check = _RK4_FLIGHT_FRACTION * res.flight_time
        p0 = cyl.floor_point(0.0, 0.0, h0)
        v0 = cyl.local_basis(p0).T @ v_launch.to_array()
        p_rk, _ = integrate_rotating(cyl, p0, v0, t_check, max_step=_RK4_MAX_STEP)
        p_exact = cyl.inertial_to_habitat(
            p0 + (v0 + np.cross(omega * Z_HAT, p0)) * t_check, t_check
        )
        results.append(
            CheckResult(
                f"rotating-frame RK4 vs exact inertial ({label}) position error [m]",
                float(np.linalg.norm(p_rk - p_exact)),
                0.0,
                _RK4_POSITION_TOL,
            )
        )
    return tuple(results)
