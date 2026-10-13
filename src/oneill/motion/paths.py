"""Prescribed motion in the habitat frame.

Every path gives habitat-frame ``position``, ``velocity`` and ``acceleration``
as functions of time. Each accepts a scalar time (result shape ``(3,)``) or an
array of times (result shape ``(n, 3)``).
"""

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import numpy as np

from oneill.core.checks import finite, positive
from oneill.core.vectors import Vec3, as_vec3
from oneill.motion.headings import heading_components, heading_name
from oneill.motion.profiles import ConstantSpeed

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import ArrayLike

    from oneill.core.cylinder import Cylinder
    from oneill.core.vectors import LocalVector
    from oneill.motion.profiles import Profile

_DEFAULT_WALK_DURATION = 10.0
_DEFAULT_STEP = 1.0e-2
_LEVEL_TOLERANCE = 1e-12


class Path(Protocol):
    """A motion with a name, a natural duration and analytic kinematics."""

    @property
    def name(self) -> str:
        """Human-readable label."""
        ...

    @property
    def duration(self) -> float:
        """Natural length of the motion in seconds."""
        ...

    def position(self, t: ArrayLike) -> Vec3:
        """Habitat-frame position at time ``t``."""
        ...

    def velocity(self, t: ArrayLike) -> Vec3:
        """Habitat-frame velocity at time ``t``."""
        ...

    def acceleration(self, t: ArrayLike) -> Vec3:
        """Habitat-frame (kinematic) acceleration at time ``t``."""
        ...


@dataclass(frozen=True)
class FloorWalk:
    """Constant-speed walk on a level (curved) deck, following its curvature."""

    cylinder: Cylinder
    heading_deg: float
    speed: float
    duration: float = _DEFAULT_WALK_DURATION
    height: float = 0.0
    theta0: float = 0.0
    z0: float = 0.0
    name: str = ""

    def __post_init__(self) -> None:
        """Validate the parameters and fill in the default name."""
        finite("heading_deg", self.heading_deg)
        finite("speed", self.speed)
        finite("theta0", self.theta0)
        finite("z0", self.z0)
        positive("duration", self.duration)
        if not (
            math.isfinite(self.height) and 0.0 <= self.height < self.cylinder.radius
        ):
            msg = f"height must satisfy 0 <= height < radius, got {self.height!r}"
            raise ValueError(msg)
        if not self.name:
            label = f"walk {heading_name(self.heading_deg)} {self.speed:g} m/s"
            object.__setattr__(self, "name", label)

    @property
    def radius(self) -> float:
        """Radius of the deck being walked."""
        return self.cylinder.radius - self.height

    def _components(self) -> tuple[float, float]:
        comps = heading_components(self.heading_deg)
        return self.speed * comps.spinward, self.speed * comps.axial

    def position(self, t: ArrayLike) -> Vec3:
        """Habitat-frame position at time ``t``."""
        tt = np.asarray(t, dtype=np.float64)
        vs, va = self._components()
        r = self.radius
        theta = self.theta0 + vs * tt / r
        return np.stack(
            [r * np.cos(theta), r * np.sin(theta), self.z0 + va * tt], axis=-1
        )

    def velocity(self, t: ArrayLike) -> Vec3:
        """Habitat-frame velocity at time ``t``."""
        tt = np.asarray(t, dtype=np.float64)
        vs, va = self._components()
        theta = self.theta0 + vs * tt / self.radius
        return np.stack(
            [-vs * np.sin(theta), vs * np.cos(theta), np.full_like(tt, va)], axis=-1
        )

    def acceleration(self, t: ArrayLike) -> Vec3:
        """Habitat-frame acceleration at time ``t`` (centripetal)."""
        tt = np.asarray(t, dtype=np.float64)
        vs, _ = self._components()
        r = self.radius
        theta = self.theta0 + vs * tt / r
        centripetal = vs * vs / r
        return np.stack(
            [
                -centripetal * np.cos(theta),
                -centripetal * np.sin(theta),
                np.zeros_like(tt),
            ],
            axis=-1,
        )


@dataclass(frozen=True)
class StraightPath:
    """A straight line fixed in the habitat (ladder, elevator shaft, staircase).

    The direction is given in the local frame at the start point; ``profile``
    says how far along the line the body is at each time.
    """

    cylinder: Cylinder
    heading_deg: float
    elevation_deg: float
    profile: Profile
    duration: float
    start_height: float = 0.0
    theta0: float = 0.0
    z0: float = 0.0
    name: str = ""
    _origin: Vec3 = field(init=False, repr=False, compare=False)
    _direction: Vec3 = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """Validate, resolve the line in the habitat frame, and name it."""
        finite("heading_deg", self.heading_deg)
        finite("elevation_deg", self.elevation_deg)
        finite("theta0", self.theta0)
        finite("z0", self.z0)
        positive("duration", self.duration)
        if not (
            math.isfinite(self.start_height)
            and 0.0 <= self.start_height <= self.cylinder.radius
        ):
            msg = f"start_height must lie within the habitat, got {self.start_height!r}"
            raise ValueError(msg)
        origin = self.cylinder.floor_point(self.theta0, self.z0, self.start_height)
        basis = self.cylinder.local_basis(origin)
        local = heading_components(self.heading_deg, self.elevation_deg).to_array()
        object.__setattr__(self, "_origin", origin)
        object.__setattr__(self, "_direction", basis.T @ local)
        if not self.name:
            label = (
                f"straight {heading_name(self.heading_deg)} "
                f"elev {self.elevation_deg:g} deg"
            )
            object.__setattr__(self, "name", label)

    def position(self, t: ArrayLike) -> Vec3:
        """Habitat-frame position at time ``t``."""
        s = self.profile.s(t)
        return self._origin + s[..., None] * self._direction

    def velocity(self, t: ArrayLike) -> Vec3:
        """Habitat-frame velocity at time ``t``."""
        return self.profile.ds(t)[..., None] * self._direction

    def acceleration(self, t: ArrayLike) -> Vec3:
        """Habitat-frame acceleration at time ``t``."""
        return self.profile.dds(t)[..., None] * self._direction


@dataclass(frozen=True)
class NumericPath:
    """An arbitrary position callable, differentiated with a 5-point stencil.

    Truncation error is O(step^4), so the default step is ample for smooth
    motion. The callable takes a scalar time and returns a ``(3,)`` position.
    """

    func: Callable[[float], ArrayLike]
    duration: float
    name: str = "numeric path"
    step: float = _DEFAULT_STEP

    def __post_init__(self) -> None:
        """Validate the parameters."""
        positive("duration", self.duration)
        positive("step", self.step)

    def _sample(self, t: float) -> Vec3:
        return as_vec3(self.func(t))

    def _map(self, t: ArrayLike, fn: Callable[[float], Vec3]) -> Vec3:
        tt = np.asarray(t, dtype=np.float64)
        flat = np.stack([fn(float(x)) for x in tt.ravel()], axis=0)
        return flat.reshape((*tt.shape, 3))

    def position(self, t: ArrayLike) -> Vec3:
        """Habitat-frame position at time ``t``."""
        return self._map(t, self._sample)

    def velocity(self, t: ArrayLike) -> Vec3:
        """Habitat-frame velocity at time ``t`` (4th-order central difference)."""
        h = self.step

        def deriv(x: float) -> Vec3:
            f = self._sample
            return (-f(x + 2 * h) + 8 * f(x + h) - 8 * f(x - h) + f(x - 2 * h)) / (
                12 * h
            )

        return self._map(t, deriv)

    def acceleration(self, t: ArrayLike) -> Vec3:
        """Habitat-frame acceleration at time ``t`` (4th-order central difference)."""
        h = self.step

        def deriv(x: float) -> Vec3:
            f = self._sample
            return (
                -f(x + 2 * h) + 16 * f(x + h) - 30 * f(x) + 16 * f(x - h) - f(x - 2 * h)
            ) / (12 * h * h)

        return self._map(t, deriv)


def velocity_path(
    cylinder: Cylinder, v_local: LocalVector, height: float = 0.0, name: str = ""
) -> FloorWalk | StraightPath:
    """Constant velocity given as (spinward, axial, up).

    Level motion follows the floor; anything with a vertical part is a straight
    line fixed in the habitat.
    """
    speed = v_local.speed
    heading = math.degrees(math.atan2(v_local.spinward, v_local.axial))
    if abs(v_local.up) < _LEVEL_TOLERANCE:
        return FloorWalk(cylinder, heading, speed, height=height, name=name)
    elevation = math.degrees(
        math.atan2(v_local.up, math.hypot(v_local.spinward, v_local.axial))
    )
    return StraightPath(
        cylinder,
        heading,
        elevation,
        ConstantSpeed(speed),
        1.0,
        start_height=height,
        name=name,
    )
