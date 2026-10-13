"""The spinning habitat: parameters, geometry and the person-centred local frame."""

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import finite, positive
from oneill.core.vectors import Vec3, as_vec3

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

G_STANDARD = 9.80665
"""Standard gravity in m/s^2."""

_AXIS_TOLERANCE = 1e-9
_UNSET = math.nan


class AxisSingularityError(ValueError):
    """Raised when a local frame is requested on the spin axis (undefined there)."""


@dataclass(frozen=True, slots=True)
class HabitatSummary:
    """Headline numbers describing a habitat (all derived from a ``Cylinder``)."""

    radius: float
    length: float
    omega: float
    rpm: float
    period: float
    rim_speed: float
    surface_gravity: float
    gravity_in_g: float
    gradient_pct_per_m: float
    coriolis_scale: float
    coriolis_pct_g_per_mps: float
    plus_one_pct_speed: float
    minus_one_pct_speed: float
    weightless_speed: float
    spin_accel: float


@dataclass(frozen=True, slots=True)
class Cylinder:
    """An O'Neill cylinder spinning about its z axis.

    Attributes:
        radius: Axis to walking surface, in metres.
        length: End cap to end cap, in metres.
        omega: Spin rate in rad/s. Left unset (NaN), it defaults to whatever gives
            1 g at the floor.
        spin_accel: Angular acceleration in rad/s^2; nonzero only while spinning
            up or down, and the source of the Euler force.
    """

    radius: float = 3200.0
    length: float = 20000.0
    omega: float = field(default=_UNSET)
    spin_accel: float = 0.0

    def __post_init__(self) -> None:
        """Validate the parameters and fill in the 1 g default spin rate."""
        positive("radius", self.radius)
        positive("length", self.length)
        finite("spin_accel", self.spin_accel)
        if math.isnan(self.omega):
            object.__setattr__(self, "omega", math.sqrt(G_STANDARD / self.radius))
        else:
            positive("omega", self.omega)

    @classmethod
    def from_gravity(
        cls,
        radius: float = 3200.0,
        length: float = 20000.0,
        surface_g: float = G_STANDARD,
        spin_accel: float = 0.0,
    ) -> Cylinder:
        """Build a habitat spun to give ``surface_g`` m/s^2 at the floor."""
        positive("surface_g", surface_g)
        return cls(radius, length, math.sqrt(surface_g / radius), spin_accel)

    @classmethod
    def from_rpm(
        cls,
        radius: float = 3200.0,
        length: float = 20000.0,
        rpm: float = 0.5,
        spin_accel: float = 0.0,
    ) -> Cylinder:
        """Build a habitat spinning at ``rpm`` revolutions per minute."""
        return cls(radius, length, rpm * 2.0 * math.pi / 60.0, spin_accel)

    # --- derived quantities -------------------------------------------------
    @property
    def rpm(self) -> float:
        """Spin rate in revolutions per minute."""
        return self.omega * 60.0 / (2.0 * math.pi)

    @property
    def period(self) -> float:
        """Time for one revolution, in seconds."""
        return 2.0 * math.pi / self.omega

    @property
    def rim_speed(self) -> float:
        """Speed of the floor in the inertial frame, in m/s."""
        return self.omega * self.radius

    @property
    def surface_gravity(self) -> float:
        """Centrifugal acceleration at the floor, in m/s^2."""
        return self.omega**2 * self.radius

    def gravity_at_height(self, height: float) -> float:
        """Centrifugal acceleration ``height`` metres above the floor."""
        return self.omega**2 * (self.radius - height)

    def omega_at(self, t: float) -> float:
        """Spin rate at time ``t`` (rad/s)."""
        return self.omega + self.spin_accel * t

    def frame_angle(self, t: float) -> float:
        """Angle through which the habitat has turned by time ``t`` (rad)."""
        return self.omega * t + 0.5 * self.spin_accel * t * t

    def summary(self) -> HabitatSummary:
        """Headline numbers for display."""
        g = self.surface_gravity
        rim = self.rim_speed
        return HabitatSummary(
            radius=self.radius,
            length=self.length,
            omega=self.omega,
            rpm=self.rpm,
            period=self.period,
            rim_speed=rim,
            surface_gravity=g,
            gravity_in_g=g / G_STANDARD,
            gradient_pct_per_m=100.0 / self.radius,
            coriolis_scale=2.0 * self.omega,
            coriolis_pct_g_per_mps=100.0 * 2.0 * self.omega / g,
            plus_one_pct_speed=rim * (math.sqrt(1.01) - 1.0),
            minus_one_pct_speed=rim * (1.0 - math.sqrt(0.99)),
            weightless_speed=rim,
            spin_accel=self.spin_accel,
        )

    # --- geometry -------------------------------------------------------------
    def floor_point(
        self, theta: float = 0.0, z: float = 0.0, height: float = 0.0
    ) -> Vec3:
        """Habitat-frame position at angle ``theta``, axial ``z`` and given height."""
        r = self.radius - height
        return np.array([r * math.cos(theta), r * math.sin(theta), z])

    def local_basis(self, p: ArrayLike) -> Vec3:
        """Rows are the (spinward, axial, up) unit vectors at habitat position ``p``.

        ``p`` may be a single ``(3,)`` point or an ``(n, 3)`` stack, giving a
        ``(3, 3)`` or ``(n, 3, 3)`` result.

        Raises:
            AxisSingularityError: if any point lies on the spin axis.
        """
        pts = as_vec3(p)
        rxy = np.hypot(pts[..., 0], pts[..., 1])
        if np.any(rxy < _AXIS_TOLERANCE):
            msg = "local frame undefined exactly on the axis"
            raise AxisSingularityError(msg)
        out = np.stack([pts[..., 0] / rxy, pts[..., 1] / rxy, np.zeros_like(rxy)], -1)
        zeros = np.zeros_like(out)
        zeros[..., 2] = 1.0
        spin = np.stack([-out[..., 1], out[..., 0], np.zeros_like(rxy)], -1)
        return np.stack([spin, zeros, -out], axis=-2)

    def height_of(self, p: ArrayLike) -> float | Vec3:
        """Height above the floor of habitat position(s) ``p``."""
        pts = as_vec3(p)
        return self.radius - np.hypot(pts[..., 0], pts[..., 1])

    def inertial_to_habitat(self, point: ArrayLike, t: float) -> Vec3:
        """Express inertial-frame position(s) in the habitat frame at time ``t``.

        The two frames coincide at ``t = 0``.
        """
        pts = as_vec3(point)
        phi = self.frame_angle(t)
        c, s = math.cos(phi), math.sin(phi)
        x = pts[..., 0] * c + pts[..., 1] * s
        y = -pts[..., 0] * s + pts[..., 1] * c
        return np.stack([x, y, pts[..., 2]], axis=-1)
