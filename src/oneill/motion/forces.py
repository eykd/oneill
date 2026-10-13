"""Force analysis: what a person following a prescribed path feels."""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import positive
from oneill.core.rotating_frame import (
    centrifugal_acc,
    coriolis_acc,
    euler_acc,
    spin_rate,
)

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from oneill.core.cylinder import Cylinder
    from oneill.core.vectors import Vec3
    from oneill.motion.paths import Path


@dataclass(frozen=True, eq=False)
class ForceHistory:
    """Accelerations along a path, in local (spinward, axial, up) components.

    Vector fields have shape ``(n, 3)``; scalar fields have shape ``(n,)``.
    All accelerations are in m/s^2 and ``g_ref`` is the floor gravity used to
    express them as percentages.
    """

    times: Vec3
    heights: Vec3
    vel: Vec3  # m/s
    acc: Vec3  # kinematic acceleration of the path (habitat frame)
    centrifugal: Vec3
    coriolis: Vec3
    euler: Vec3
    apparent: Vec3  # what you feel as "gravity"
    g_ref: float

    def __post_init__(self) -> None:
        """Check that every field is consistent with ``times``."""
        n = len(self.times)
        if self.heights.shape != (n,):
            msg = f"heights must have shape ({n},), got {self.heights.shape}"
            raise ValueError(msg)
        for name in ("vel", "acc", "centrifugal", "coriolis", "euler", "apparent"):
            arr = getattr(self, name)
            if arr.shape != (n, 3):
                msg = f"{name} must have shape ({n}, 3), got {arr.shape}"
                raise ValueError(msg)
        positive("g_ref", self.g_ref)

    @property
    def felt_weight(self) -> Vec3:
        """Pressure toward the floor (positive = pressed down)."""
        return -self.apparent[:, 2]

    @property
    def lateral_mag(self) -> Vec3:
        """Magnitude of the sideways (spinward, axial) part of apparent gravity."""
        return np.hypot(self.apparent[:, 0], self.apparent[:, 1])

    @property
    def tilt_deg(self) -> Vec3:
        """Lean of the apparent vertical away from the local vertical."""
        return np.degrees(np.arctan2(self.lateral_mag, self.felt_weight))

    def support_force(self, mass: float) -> Vec3:
        """Force the floor/ladder/stair must supply for a body of ``mass`` kg."""
        return -mass * self.apparent

    def pct(self, x: ArrayLike) -> Vec3:
        """Express accelerations as a percentage of floor gravity."""
        return 100.0 * np.asarray(x) / self.g_ref


def analyze_path(cyl: Cylinder, path: Path, times: ArrayLike) -> ForceHistory:
    """Evaluate the fictitious and apparent accelerations along ``path``.

    Apparent gravity is ``centrifugal + coriolis + euler - a_path``; the support
    must supply ``-m`` times it.

    Raises:
        AxisSingularityError: if the path passes through the spin axis.
    """
    ts = np.atleast_1d(np.asarray(times, dtype=np.float64))
    pos = path.position(ts)
    vel = path.velocity(ts)
    acc = path.acceleration(ts)
    omega = spin_rate(cyl, ts)
    cen = centrifugal_acc(omega, pos)
    cor = coriolis_acc(omega, vel)
    eul = euler_acc(cyl.spin_accel, pos)
    app = cen + cor + eul - acc
    basis = cyl.local_basis(pos)

    def local(x: Vec3) -> Vec3:
        return np.asarray(np.einsum("nij,nj->ni", basis, x), dtype=np.float64)

    return ForceHistory(
        ts,
        np.asarray(cyl.height_of(pos), dtype=np.float64),
        local(vel),
        local(acc),
        local(cen),
        local(cor),
        local(eul),
        local(app),
        cyl.surface_gravity,
    )
