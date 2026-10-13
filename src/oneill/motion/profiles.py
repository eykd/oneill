"""Distance-versus-time profiles with analytic first and second derivatives.

Each profile gives ``s(t)`` (distance travelled), ``ds(t)`` (speed) and
``dds(t)`` (acceleration). All three accept a scalar or an array and return an
array of the same shape.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

import numpy as np

from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

    type Floats = NDArray[np.float64]
    type Mask = NDArray[np.bool_]


class Profile(Protocol):
    """Anything that can say how far, how fast and how hard along a line."""

    def s(self, t: ArrayLike) -> Floats:
        """Distance travelled at time ``t``."""
        ...

    def ds(self, t: ArrayLike) -> Floats:
        """Speed at time ``t``."""
        ...

    def dds(self, t: ArrayLike) -> Floats:
        """Acceleration at time ``t``."""
        ...


@dataclass(frozen=True, slots=True)
class ConstantSpeed:
    """Uniform motion, linear for all ``t`` (including before zero)."""

    speed: float

    def __post_init__(self) -> None:
        """Reject a non-finite speed."""
        finite("speed", self.speed)

    def s(self, t: ArrayLike) -> Floats:
        """Distance travelled at time ``t``."""
        return self.speed * np.asarray(t, dtype=np.float64)

    def ds(self, t: ArrayLike) -> Floats:
        """Speed at time ``t`` (constant)."""
        return np.full_like(np.asarray(t, dtype=np.float64), self.speed)

    def dds(self, t: ArrayLike) -> Floats:
        """Acceleration at time ``t`` (zero)."""
        return np.zeros_like(np.asarray(t, dtype=np.float64))


@dataclass(frozen=True, slots=True)
class Trapezoid:
    """Accelerate, cruise, decelerate over ``distance``.

    Phases are half-open except that the ride includes both of its ends:
    acceleration on ``[0, t_acc)``, cruise on ``[t_acc, t_acc + t_cruise)``,
    braking on ``[t_acc + t_cruise, duration]``. Outside the ride the body is at
    rest: at 0 before it, at ``distance`` after.
    When ``v_max`` is never reached the profile is a triangle (no cruise).
    """

    distance: float
    v_max: float
    a_max: float

    def __post_init__(self) -> None:
        """Require a positive distance, speed limit and acceleration limit."""
        positive("distance", self.distance)
        positive("v_max", self.v_max)
        positive("a_max", self.a_max)

    @property
    def t_acc(self) -> float:
        """Duration of the acceleration (and of the braking) phase."""
        t_full = self.v_max / self.a_max
        if self.a_max * t_full**2 >= self.distance:
            return math.sqrt(self.distance / self.a_max)
        return t_full

    @property
    def t_cruise(self) -> float:
        """Duration of the constant-speed phase (zero for a triangle)."""
        return max(0.0, (self.distance - self.a_max * self.t_acc**2) / self.v_max)

    @property
    def v_peak(self) -> float:
        """Highest speed reached."""
        return self.a_max * self.t_acc

    @property
    def duration(self) -> float:
        """Total ride time."""
        return 2.0 * self.t_acc + self.t_cruise

    def _phases(self, t: Floats) -> tuple[Mask, Mask, Mask, Mask]:
        t1 = self.t_acc
        t2 = t1 + self.t_cruise
        return (
            (t >= 0.0) & (t < t1),
            (t >= t1) & (t < t2),
            (t >= t2) & (t <= self.duration),
            t > self.duration,
        )

    def s(self, t: ArrayLike) -> Floats:
        """Distance travelled at time ``t``."""
        tt = np.asarray(t, dtype=np.float64)
        accel, cruise, brake, done = self._phases(tt)
        a = self.a_max
        d_acc = 0.5 * a * self.t_acc**2
        return np.select(
            [accel, cruise, brake, done],
            [
                0.5 * a * tt * tt,
                d_acc + self.v_peak * (tt - self.t_acc),
                self.distance - 0.5 * a * (self.duration - tt) ** 2,
                np.full_like(tt, self.distance),
            ],
            default=0.0,
        )

    def ds(self, t: ArrayLike) -> Floats:
        """Speed at time ``t``."""
        tt = np.asarray(t, dtype=np.float64)
        accel, cruise, brake, _ = self._phases(tt)
        return np.select(
            [accel, cruise, brake],
            [
                self.a_max * tt,
                np.full_like(tt, self.v_peak),
                self.a_max * (self.duration - tt),
            ],
            default=0.0,
        )

    def dds(self, t: ArrayLike) -> Floats:
        """Acceleration at time ``t``."""
        tt = np.asarray(t, dtype=np.float64)
        accel, _, brake, _ = self._phases(tt)
        return np.select(
            [accel, brake],
            [np.full_like(tt, self.a_max), np.full_like(tt, -self.a_max)],
            default=0.0,
        )
