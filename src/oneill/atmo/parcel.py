"""L2: a single air parcel moving through the still habitat atmosphere.

The parcel feels, per unit mass, in the habitat frame:

* centrifugal and Coriolis accelerations (:mod:`oneill.core.rotating_frame`);
* the environment pressure-gradient force ``-grad(p_env) / rho_p``, with the
  environment in hydrostatic balance (``dp/dr = rho_env omega^2 r``, from
  :mod:`oneill.atmo.equilibrium`), so the net radial force is the buoyancy
  ``omega^2 r (1 - rho_env / rho_p)`` directed toward the axis when the parcel
  is light;
* linear drag ``-k v`` (the air is at rest in the habitat frame).

Dry parcels conserve potential temperature ``theta``: at pressure ``p`` the parcel
temperature is ``theta (p / p_rim)^kappa`` and its density ``p / (R_s T_p)``.

Angular momentum about the axis is conserved when ``k = 0`` (all other forces are
radial), so a parcel rising by ``dr`` acquires a relative spinward velocity
``omega (R0^2 / r - r)``, about ``2 omega dr`` (11 m/s per 100 m). That is an
**upper bound for compact 3-D thermals**. For plumes uniform along the axis the
pressure field absorbs the Coriolis force and the drift is much smaller
(``docs/research/atmosphere-literature.md``, Q7-Q8).

The model uses :func:`oneill.core.integrate.rk4` directly (atmo may not import the
sibling ``motion`` package, where ``integrate_rotating`` lives). Moist latent
heating is not modelled: the cloud-base time is for a dry parcel.
"""

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from oneill.atmo.equilibrium import adiabatic_profile, lifting_condensation_level
from oneill.core.checks import finite, positive
from oneill.core.integrate import rk4
from oneill.core.rotating_frame import centrifugal_acc, coriolis_acc
from oneill.core.vectors import LocalVector

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import ArrayLike, NDArray

    from oneill.atmo.equilibrium import AtmosphereProfile
    from oneill.atmo.gas import Gas
    from oneill.core.cylinder import Cylinder
    from oneill.core.vectors import Vec3

    type Environment = Callable[[Cylinder, Gas, ArrayLike], AtmosphereProfile]

RISE = 100.0
"""Default rise in metres at which the deflection is reported."""


@dataclass(frozen=True, slots=True)
class Release:
    """How and where the parcel starts.

    Attributes:
        height: Starting height above the floor in m.
        velocity: Initial velocity relative to still air, as (spinward, axial, up).
        temperature_excess: Parcel minus environment temperature, in K, at
            ``reference_height``; sets the conserved potential temperature.
        reference_height: Height where the excess applies; defaults to ``height``.
            Releasing at another height than the reference displaces a neutral
            parcel from its level of neutral buoyancy.
        drag: Linear drag coefficient ``k`` in 1/s.
    """

    height: float = 0.0
    velocity: LocalVector = field(default_factory=lambda: LocalVector(0.0, 0.0, 0.0))
    temperature_excess: float = 0.0
    reference_height: float | None = None
    drag: float = 0.0

    def __post_init__(self) -> None:
        """Validate the drag coefficient and the finite excess."""
        finite("height", self.height)
        finite("temperature_excess", self.temperature_excess)
        if not (math.isfinite(self.drag) and self.drag >= 0.0):
            msg = f"drag must be non-negative and finite, got {self.drag!r}"
            raise ValueError(msg)


@dataclass(frozen=True, slots=True)
class RunSettings:
    """Time stepping and stopping rules.

    Attributes:
        t_end: Maximum simulated time in s.
        sample_dt: Output interval in s.
        max_step: RK4 step bound in s.
        stop_height: Stop once the parcel reaches this height (m), if given.
        rise: Rise in m (above the release height) at which the deflection and
            drift speed are read off.
    """

    t_end: float = 300.0
    sample_dt: float = 0.5
    max_step: float = 0.25
    stop_height: float | None = None
    rise: float = RISE


@dataclass(frozen=True, slots=True)
class ParcelResult:
    """A parcel trajectory and its headline numbers.

    Attributes:
        t: Sample times in s.
        position: Habitat-frame positions ``(n, 3)`` in m.
        velocity: Habitat-frame velocities relative to the air ``(n, 3)`` in m/s.
        height: Height above the floor in m.
        spinward: Along-floor spinward displacement ``R dtheta`` from the start, m.
        axial: Axial displacement in m.
        speed: Speed relative to the air in m/s (what a resident feels).
        max_rise: Largest height gain above the release point, in m.
        deflection: Spinward displacement (m) when the parcel first has risen
            ``settings.rise`` metres (100 m by default), or ``None`` if it never does.
            A buoyant parcel with little drag is trapped by the Coriolis force
            and oscillates, so for small ``k`` this is usually ``None``.
        drift_speed: Spinward velocity (m/s) at that point.
        peak_speed: Largest relative speed on the trajectory, m/s.
        lcl_height: True height of the cloud base (m), if a humidity was given and
            the cloud base is below the axis.
        time_to_lcl: Time to first reach ``lcl_height`` in s, or ``None``.
        ended_at_floor: ``True`` if the run stopped because the parcel hit the floor.
    """

    t: NDArray[np.float64]
    position: NDArray[np.float64]
    velocity: NDArray[np.float64]
    height: NDArray[np.float64]
    spinward: NDArray[np.float64]
    axial: NDArray[np.float64]
    speed: NDArray[np.float64]
    max_rise: float
    deflection: float | None
    drift_speed: float | None
    peak_speed: float
    lcl_height: float | None
    time_to_lcl: float | None
    ended_at_floor: bool


def _first_crossing(
    t: NDArray[np.float64], values: NDArray[np.float64], level: float
) -> tuple[float, int] | None:
    """Interpolated time and first-sample index at which ``values`` reach ``level``."""
    above = np.flatnonzero(values >= level)
    if above.size == 0:
        return None
    i = int(above[0])
    if i == 0:
        return float(t[0]), 0
    frac = (level - values[i - 1]) / (values[i] - values[i - 1])
    return float(t[i - 1] + frac * (t[i] - t[i - 1])), i


def _interp_at(
    t: NDArray[np.float64], series: NDArray[np.float64], time: float, i: int
) -> NDArray[np.float64]:
    frac = (time - t[i - 1]) / (t[i] - t[i - 1])
    return np.asarray(series[i - 1] + frac * (series[i] - series[i - 1]))


def simulate_parcel(
    cyl: Cylinder,
    gas: Gas,
    release: Release | None = None,
    *,
    environment: Environment = adiabatic_profile,
    settings: RunSettings | None = None,
    relative_humidity: float | None = None,
) -> ParcelResult:
    """Integrate a dry parcel through a still environment.

    Args:
        cyl: The habitat (constant spin).
        gas: Floor gas state.
        release: Initial conditions; defaults to a neutral parcel at rest on the floor.
        environment: Profile builder from :mod:`oneill.atmo.equilibrium`;
            :func:`~oneill.atmo.equilibrium.adiabatic_profile` (neutral, the default)
            or :func:`~oneill.atmo.equilibrium.isothermal_profile` (stable).
        settings: Time stepping and stopping rules.
        relative_humidity: If given, the floor humidity (fraction) used to find
            the cloud base and the dry-parcel time to reach it.

    Raises:
        ValueError: for a spinning-up habitat, or a release or step outside bounds.
    """
    rel = Release() if release is None else release
    cfg = RunSettings() if settings is None else settings
    positive("t_end", cfg.t_end)
    positive("sample_dt", cfg.sample_dt)
    positive("max_step", cfg.max_step)
    positive("rise", cfg.rise)
    if cyl.spin_accel != 0.0:
        msg = "the parcel model assumes a constant spin rate"
        raise ValueError(msg)

    radius = cyl.radius
    omega = cyl.omega
    ref_height = rel.height if rel.reference_height is None else rel.reference_height
    ref = environment(cyl, gas, ref_height)
    kappa = gas.kappa
    theta = (
        float(ref.temperature + rel.temperature_excess)
        * (gas.p_rim / float(ref.pressure)) ** kappa
    )

    def deriv(_t: float, y: Vec3) -> Vec3:
        p, v = y[:3], y[3:]
        r = float(np.hypot(p[0], p[1]))
        env = environment(cyl, gas, min(max(radius - r, 0.0), radius))
        press, rho_env = float(env.pressure), float(env.density)
        rho_parcel = press / (gas.r_specific * theta * (press / gas.p_rim) ** kappa)
        radial = np.array([p[0] / r, p[1] / r, 0.0])
        pressure_force = -(rho_env / rho_parcel) * omega**2 * r * radial
        acc = (
            centrifugal_acc(omega, p)
            + coriolis_acc(omega, v)
            + pressure_force
            - rel.drag * v
        )
        return np.concatenate([v, acc])

    p0 = cyl.floor_point(0.0, 0.0, rel.height)
    v0 = cyl.local_basis(p0).T @ rel.velocity.to_array()
    y = np.concatenate([p0, v0])
    times = [0.0]
    states = [y]
    ended_at_floor = False
    n_out = math.ceil(cfg.t_end / cfg.sample_dt)
    for k in range(1, n_out + 1):
        t_prev, t_next = times[-1], min(k * cfg.sample_dt, cfg.t_end)
        y = rk4(deriv, t_prev, y, t_next, cfg.max_step)
        if radius - float(np.hypot(y[0], y[1])) < 0.0:
            ended_at_floor = True
            break
        times.append(t_next)
        states.append(y)
        if (
            cfg.stop_height is not None
            and radius - float(np.hypot(y[0], y[1])) >= cfg.stop_height
        ):
            break

    return _summarise(
        cyl,
        gas,
        np.array(times),
        np.array(states),
        (rel.height, cfg.rise),
        (ended_at_floor, relative_humidity),
    )


def _summarise(
    cyl: Cylinder,
    gas: Gas,
    t: NDArray[np.float64],
    states: NDArray[np.float64],
    start: tuple[float, float],
    outcome: tuple[bool, float | None],
) -> ParcelResult:
    start_height, rise_wanted = start
    ended_at_floor, relative_humidity = outcome
    pos, vel = states[:, :3], states[:, 3:]
    radial = np.hypot(pos[:, 0], pos[:, 1])
    height = cyl.radius - radial
    angle = np.unwrap(np.arctan2(pos[:, 1], pos[:, 0]))
    spinward = cyl.radius * (angle - angle[0])
    tangential = (pos[:, 0] * vel[:, 1] - pos[:, 1] * vel[:, 0]) / radial
    speed = np.linalg.norm(vel, axis=1)

    deflection = drift = None
    rise = _first_crossing(t, height, start_height + rise_wanted)
    if rise is not None:
        time, i = rise
        deflection = float(_interp_at(t, spinward, time, i))
        drift = float(_interp_at(t, tangential, time, i))

    lcl_height = time_to_lcl = None
    if relative_humidity is not None:
        lcl_height = lifting_condensation_level(cyl, gas, relative_humidity).height
        if lcl_height is not None:
            hit = _first_crossing(t, height, lcl_height)
            time_to_lcl = None if hit is None else hit[0]

    return ParcelResult(
        t=t,
        position=pos,
        velocity=vel,
        height=height,
        spinward=spinward,
        axial=pos[:, 2] - pos[0, 2],
        speed=speed,
        max_rise=float(height.max() - start_height),
        deflection=deflection,
        drift_speed=drift,
        peak_speed=float(speed.max()),
        lcl_height=lcl_height,
        time_to_lcl=time_to_lcl,
        ended_at_floor=ended_at_floor,
    )
