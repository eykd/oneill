"""Free flight: drops, jumps and throws.

A released body feels no force, so in the inertial frame it moves in a straight
line. ``ballistic`` solves that exactly; ``integrate_rotating`` integrates the
same motion in the habitat frame with fictitious forces as an independent check.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

import numpy as np

from oneill.core.checks import finite
from oneill.core.integrate import rk4
from oneill.core.rotating_frame import fictitious_acc
from oneill.core.vectors import Z_HAT, LocalVector, Vec3, as_vec3

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from oneill.core.cylinder import Cylinder

_MIN_SAMPLES = 2
_REST_TOLERANCE = 1e-18


class ExtraAcceleration(Protocol):
    """An additional acceleration (drag, buoyancy, ...) acting on the body."""

    def __call__(self, t: float, p: Vec3, v: Vec3, /) -> Vec3:
        """Acceleration at time ``t`` for habitat position ``p`` and velocity ``v``."""
        ...


@dataclass(frozen=True, eq=False)
class BallisticResult:
    """Outcome of a free-flight calculation, as seen by habitat residents."""

    label: str
    start_height: float
    v_launch: LocalVector
    flight_time: float
    landing_spin: float  # along-floor distance from below launch (+ = spinward)
    landing_axial: float
    naive_time: float  # flat-ground, uniform-gravity prediction
    naive_spin: float
    naive_axial: float
    traj_t: Vec3
    traj_spin: Vec3
    traj_axial: Vec3
    traj_height: Vec3
    hits_endcap: bool


def _flight_time(a: float, b: float, c: float) -> float:
    """Positive root of ``a T^2 + b T + c = 0`` for ``a > 0`` and ``c <= 0``.

    Uses ``q = -(b + sign(b) sqrt(disc)) / 2`` and the roots ``q/a`` and ``c/q`` so
    that neither cancels catastrophically; ``c = 0, b >= 0`` (launched outward
    from the floor) gives ``T = 0``.
    """
    root = math.sqrt(max(b * b - 4.0 * a * c, 0.0))
    q = -0.5 * (b + (root if b >= 0.0 else -root))
    if q == 0.0:
        return 0.0
    return max(q / a, c / q)


def ballistic(
    cyl: Cylinder,
    start_height: float,
    v_local: LocalVector,
    label: str = "",
    *,
    theta0: float = 0.0,
    z0: float = 0.0,
    n: int = 400,
) -> BallisticResult:
    """Exact free flight from ``start_height`` with launch velocity ``v_local``.

    Raises:
        ValueError: if the launch point is outside ``0 <= start_height < radius``,
            ``n`` is below 2, or the body is at rest in the inertial frame.
    """
    if not (math.isfinite(start_height) and 0.0 <= start_height < cyl.radius):
        msg = "start_height must satisfy 0 <= start_height < radius"
        raise ValueError(msg)
    if n < _MIN_SAMPLES:
        msg = f"n must be at least {_MIN_SAMPLES}, got {n}"
        raise ValueError(msg)
    finite("theta0", theta0)
    finite("z0", z0)
    p0 = cyl.floor_point(theta0, z0, start_height)
    basis = cyl.local_basis(p0)
    v_rel = basis.T @ v_local.to_array()
    # Inertial velocity (the frames coincide at t = 0).
    big_v = v_rel + np.cross(cyl.omega_at(0.0) * Z_HAT, p0)
    a = float(big_v[0] ** 2 + big_v[1] ** 2)
    b = 2.0 * float(p0[0] * big_v[0] + p0[1] * big_v[1])
    c = float(p0[0] ** 2 + p0[1] ** 2 - cyl.radius**2)
    if a < _REST_TOLERANCE:
        msg = "object is at rest in the inertial frame and never reaches the floor"
        raise ValueError(msg)
    flight = _flight_time(a, b, c)

    ts = np.linspace(0.0, flight, n)
    pts = p0[None, :] + ts[:, None] * big_v[None, :]
    # An inertial straight chord sweeps at most pi radians. Subtract the habitat's
    # full rotation explicitly so long flights do not alias between samples.
    swept = np.arctan2(pts[:, 1], pts[:, 0]) - math.atan2(p0[1], p0[0])
    swept = (swept + math.pi) % (2.0 * math.pi) - math.pi
    spin = cyl.radius * (swept - (cyl.omega * ts + 0.5 * cyl.spin_accel * ts**2))
    height = cyl.radius - np.hypot(pts[:, 0], pts[:, 1])
    axial = pts[:, 2] - z0

    g = cyl.surface_gravity
    naive_t = (v_local.up + math.sqrt(v_local.up**2 + 2.0 * g * start_height)) / g
    return BallisticResult(
        label=label,
        start_height=start_height,
        v_launch=v_local,
        flight_time=flight,
        landing_spin=float(spin[-1]),
        landing_axial=float(axial[-1]),
        naive_time=naive_t,
        naive_spin=v_local.spinward * naive_t,
        naive_axial=v_local.axial * naive_t,
        traj_t=ts,
        traj_spin=spin,
        traj_axial=axial,
        traj_height=height,
        hits_endcap=bool(abs(pts[-1, 2]) > cyl.length / 2.0),
    )


def integrate_rotating(
    cyl: Cylinder,
    p0: ArrayLike,
    v0: ArrayLike,
    t_end: float,
    *,
    max_step: float = 1e-3,
    extra_acc: ExtraAcceleration | None = None,
) -> tuple[Vec3, Vec3]:
    """RK4 integration of a body directly in the habitat frame.

    The body feels centrifugal, Coriolis and Euler terms plus, optionally,
    ``extra_acc(t, p, v)``. Returns the final ``(position, velocity)``; the inputs
    are not modified.
    """
    pos0, vel0 = as_vec3(p0), as_vec3(v0)

    def deriv(t: float, y: Vec3) -> Vec3:
        p, v = y[:3], y[3:]
        acc = fictitious_acc(cyl, t, p, v)
        if extra_acc is not None:
            acc = acc + extra_acc(t, p, v)
        return np.concatenate([v, acc])

    y = rk4(deriv, 0.0, np.concatenate([pos0, vel0]), t_end, max_step)
    return y[:3], y[3:]
