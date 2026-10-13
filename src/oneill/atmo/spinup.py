"""L1: axisymmetric spin-up and spin-down of the air, one radial dimension.

``u(r, t)`` is the azimuthal wind relative to the hull. Angular momentum of the
air obeys, with the hull spinning at ``Omega(t)`` and accelerating at ``alpha``::

    du/dt = (1/r^2) d/dr( r^3 nu d(u/r)/dr ) - u/tau_E - drag - r alpha

The viscous term is written in stress form, so a solid-body ``u ~ r`` feels no
viscous force (for constant ``nu`` it equals ``nu (u'' + u'/r - u/r^2)``).
``u/tau_E`` is the Ekman pumping of the end caps (Greenspan & Howard 1963;
``tau_E = (L/2) / sqrt(nu Omega)``, ``docs/research/atmosphere-literature.md``
Q2) and can be switched off; it is written as the rate ``sqrt(nu Omega) / (L/2)``,
which stays finite as ``Omega -> 0`` (spin-up from rest). ``r alpha`` is the
apparent force from the accelerating hull.

**Floor.** By default the floor is a canopy-drag layer: a body force
``-C_D |u_r| u / z_r`` in the bottom ``z_r`` metres (default 10 m), where ``u_r`` is
the wind at ``z_r`` and ``C_D = (kappa / ln(z_r / z0))^2`` follows from a terrain
roughness (:data:`DEFAULT_FLOOR_DRAG`). The stress on the air is then
``C_D |u_r| u_r`` and the wind at ``z_r`` is the one the drag law is written for,
so quote winds there. A no-slip floor (``floor_drag=None``) is kept for analytic
tests only: with ``nu_t`` of order 1 m^2/s it leaves a sub-metre boundary layer
that no real surface has. ``floor_drag=0`` is free slip.

**Viscosity regimes** (Q4). Eddy viscosity is the larger ``nu_unstable`` where the
absolute angular momentum ``Gamma = r (Omega r + u)`` falls outward and the smaller
``nu_stable`` elsewhere, blended smoothly with the Rayleigh discriminant
``q = -(2 Omega_in + du/dr) / (2 |Omega_in|)`` (``Omega_in = Omega + u/r`` is the
air's inertial rotation rate; ``q = -1`` for any solid-body state, ``q > 0`` is
unstable) through ``(1 + tanh(q / 0.3)) / 2``. A hard switch on the previous
step's wind, with Crank-Nicolson, made the scheme ring at the sharp unstable layer
(odd-even errors of order 10 m/s). The blend, Picard iteration of the viscosity
within each step and an L-stable time scheme remove that; a small residual
2-cell stair (about 1 % of the peak wind) remains where the shear sits at the
marginal-stability threshold, as it does in rotating-channel DNS. Strong
(order 50 %) spin changes depend on the blend width (a closure choice); changes of
10 % or less do not.

**Validity.** The end-cap sink is the *linear* Greenspan-Howard result, valid for
Rossby number ``Ro = |u| / (2 Omega r) <~ 0.3``. :attr:`SpinupResult.max_rossby`
reports the largest ``Ro`` reached; fast spin changes exceed it (and turbulent
end caps are not linear), so those results are order-of-magnitude only.

The scheme is a finite-volume theta method (default backward Euler, which is
L-stable; ``theta=0.5`` Crank-Nicolson with two backward-Euler starting steps is
available for smooth problems) on a grid clustered toward the floor. It conserves
angular momentum exactly: :class:`Budget` holds the floor torque and end-cap sink
that close the books.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.atmo.surface import SUBURBAN, drag_coefficient
from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from oneill.core.cylinder import Cylinder

    type Stop = Callable[[float, NDArray[np.float64]], bool]

DRAG_HEIGHT = 10.0
"""Reference height ``z_r`` of the canopy-drag floor, in m."""

DEFAULT_FLOOR_DRAG = drag_coefficient(DRAG_HEIGHT, SUBURBAN.z0)
"""``C_D`` for the suburban roughness (z0 = 0.5 m) at 10 m: about 0.018."""

LINEAR_SINK_ROSSBY = 0.3
"""Largest Rossby number for which the linear end-cap sink is trusted."""

_MAX_STEPS = 200_000
_LAG_ITERATIONS = 80
_PICARD_ITERATIONS = 3
_STARTUP_STEPS = 2
_MIN_CELLS = 8
_BLEND_WIDTH = 0.3
_BLEND_FLOOR = 0.2
_TINY = 1e-12


@dataclass(frozen=True, slots=True)
class SpinupParams:
    """Physical and numerical parameters of the 1-D model.

    Attributes:
        nu_stable: Viscosity (molecular or eddy) where angular momentum rises
            outward, m^2/s.
        nu_unstable: Viscosity where angular momentum falls outward, m^2/s;
            defaults to ``nu_stable``.
        floor_drag: Canopy-drag coefficient ``C_D`` (default
            :data:`DEFAULT_FLOOR_DRAG`); ``None`` for a no-slip floor (analytic
            tests only), ``0`` for free slip.
        drag_height: Thickness ``z_r`` of the drag layer and height of the wind
            the drag law uses, m.
        endcap: Include the end-cap Ekman pumping sink.
        nu_endcap: Viscosity in ``tau_E``; defaults to ``nu_stable`` (the end-cap
            layer is taken to be in the stable regime).
        n_cells: Radial cells.
        stretch: Grid clustering toward the floor (``0`` is uniform; the smallest
            cell is about ``R/n`` times ``beta / cosh^2(beta)``).
        theta: Time-scheme weight, ``1`` backward Euler (default), ``0.5``
            Crank-Nicolson.
    """

    nu_stable: float
    nu_unstable: float | None = None
    floor_drag: float | None = DEFAULT_FLOOR_DRAG
    drag_height: float = DRAG_HEIGHT
    endcap: bool = True
    nu_endcap: float | None = None
    n_cells: int = 120
    stretch: float = 2.5
    theta: float = 1.0

    def __post_init__(self) -> None:
        """Validate the viscosities, drag, grid size and scheme weight."""
        positive("nu_stable", self.nu_stable)
        if self.nu_unstable is not None:
            positive("nu_unstable", self.nu_unstable)
        if self.nu_endcap is not None:
            positive("nu_endcap", self.nu_endcap)
        positive("drag_height", self.drag_height)
        if self.floor_drag is not None and not (
            math.isfinite(self.floor_drag) and self.floor_drag >= 0.0
        ):
            msg = f"floor_drag must be non-negative and finite, got {self.floor_drag!r}"
            raise ValueError(msg)
        if self.n_cells < _MIN_CELLS:
            msg = f"n_cells must be at least {_MIN_CELLS}, got {self.n_cells}"
            raise ValueError(msg)
        if not (math.isfinite(self.stretch) and self.stretch >= 0.0):
            msg = f"stretch must be non-negative and finite, got {self.stretch!r}"
            raise ValueError(msg)
        if not 0.5 <= self.theta <= 1.0:  # noqa: PLR2004
            msg = f"theta must lie in [0.5, 1], got {self.theta!r}"
            raise ValueError(msg)

    @property
    def nu_high(self) -> float:
        """Viscosity used in the unstable regime."""
        return self.nu_stable if self.nu_unstable is None else self.nu_unstable

    @property
    def nu_cap(self) -> float:
        """Viscosity used in the end-cap rate."""
        return self.nu_stable if self.nu_endcap is None else self.nu_endcap


@dataclass(frozen=True, slots=True)
class HullSpin:
    """Hull spin schedule: ``Omega`` ramps at ``alpha`` for ``ramp_time`` seconds.

    Attributes:
        omega_start: Spin rate at ``t = 0`` in rad/s.
        alpha: Angular acceleration during the ramp in rad/s^2.
        ramp_time: Duration of the ramp in s; ``None`` keeps accelerating forever.
    """

    omega_start: float = 0.0
    alpha: float = 0.0
    ramp_time: float | None = None

    def omega(self, t: float) -> float:
        """Hull spin rate at time ``t`` in rad/s."""
        end = t if self.ramp_time is None else min(t, self.ramp_time)
        return self.omega_start + self.alpha * end


@dataclass(frozen=True, slots=True)
class Budget:
    """Rates (per step) that close the angular-momentum books exactly.

    Each is a rate of change of ``sum(M_i u_i)`` where ``M_i = integral r^2 dr``
    over the cell, so ``sum(M u)`` is proportional to the relative angular
    momentum. ``dt * (wall + endcap + hull)`` equals the change over each step.

    Attributes:
        dt: Step sizes in s.
        wall: Torque exchanged with the floor (floor drag, or the no-slip wall).
        endcap: Torque removed by end-cap Ekman pumping (non-positive for ``u > 0``).
        hull: Apparent torque from the hull's angular acceleration.
    """

    dt: NDArray[np.float64]
    wall: NDArray[np.float64]
    endcap: NDArray[np.float64]
    hull: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class SpinupResult:
    """A spin-up or spin-down run.

    Attributes:
        t: Times in s, length ``n_steps + 1``.
        radius: Cell centres in m.
        faces: Cell faces in m.
        moment: Cell moments ``M_i = (r_{i+1}^3 - r_i^3) / 3``.
        u: Relative wind ``(n_steps + 1, n_cells)`` in m/s.
        omega: Hull spin rate at each time, rad/s.
        budget: Closing angular-momentum rates.
        wall_wind_is_zero: ``True`` for a no-slip floor.
        ramp_end: Time the hull stops accelerating, s (0 if it never did).
        drag_height: Height ``z_r`` of the floor-drag reference wind, m.
        endcap: Whether the end-cap sink was on.
    """

    t: NDArray[np.float64]
    radius: NDArray[np.float64]
    faces: NDArray[np.float64]
    moment: NDArray[np.float64]
    u: NDArray[np.float64]
    omega: NDArray[np.float64]
    budget: Budget
    wall_wind_is_zero: bool
    ramp_end: float
    drag_height: float = DRAG_HEIGHT
    endcap: bool = True

    @property
    def corotation(self) -> NDArray[np.float64]:
        """Fraction of solid-body angular momentum, ``1 + sum(M u) / (Omega R^4/4)``.

        Zero while the hull is not spinning.
        """
        quarter = self.faces[-1] ** 4 / 4.0
        rel = self.u @ self.moment
        out = np.zeros_like(self.omega)
        np.divide(rel, self.omega * quarter, out=out, where=self.omega > 0.0)
        return np.where(self.omega > 0.0, 1.0 + out, 0.0)

    @property
    def relative_angular_momentum(self) -> NDArray[np.float64]:
        """``sum(M u)`` at each time, proportional to the air's lag."""
        return np.asarray(self.u @ self.moment)

    @property
    def max_rossby(self) -> float:
        """Largest ``|u| / (2 Omega r)`` reached, with ``Omega`` the final hull rate."""
        omega = float(self.omega.max())
        return float(np.max(np.abs(self.u) / (2.0 * omega * self.radius)))

    @property
    def sink_valid(self) -> bool:
        """``True`` if the linear end-cap sink is trusted (Ro <= 0.3 throughout)."""
        return self.max_rossby <= LINEAR_SINK_ROSSBY


@dataclass(frozen=True, slots=True)
class LagWind:
    """Steady wind lag during a slow, steady hull spin change.

    Attributes:
        radius: Cell centres in m.
        u: Relative wind in m/s (negative while the hull spins up).
        faces: Cell faces in m.
        wall_wind_is_zero: ``True`` for a no-slip floor.
    """

    radius: NDArray[np.float64]
    u: NDArray[np.float64]
    faces: NDArray[np.float64]
    wall_wind_is_zero: bool

    @property
    def peak(self) -> float:
        """Largest ``|u|`` in m/s."""
        return float(np.max(np.abs(self.u)))

    def wind_at_height(self, height: float) -> float:
        """Wind ``height`` metres above the floor in m/s."""
        return _wind_at_height(
            self.radius, self.u, self.faces[-1], height, self.wall_wind_is_zero
        )


def _wind_at_height(
    radius: NDArray[np.float64],
    u: NDArray[np.float64],
    floor: float,
    height: float,
    no_slip: bool,  # noqa: FBT001
) -> float:
    positive("height", height)
    wall_value = 0.0 if no_slip else float(u[-1])
    return float(
        np.interp(
            floor - height,
            np.append(radius, floor),
            np.append(u, wall_value),
        )
    )


def peak_wind_at_height(result: SpinupResult, height: float) -> float:
    """Largest ``|u|`` at ``height`` metres above the floor after ``t = 0``, in m/s.

    The initial state is excluded: for a sudden change it is just the imposed
    jump. Quote winds at ``result.drag_height`` unless the floor is no-slip.
    """
    return max(
        abs(
            _wind_at_height(
                result.radius,
                row,
                result.faces[-1],
                height,
                result.wall_wind_is_zero,
            )
        )
        for row in result.u[1:]
    )


def _first_crossing(
    t: NDArray[np.float64],
    metric: NDArray[np.float64],
    level: float,
    earliest: float = 0.0,
) -> float | None:
    """Linearly interpolated time at which ``metric`` first reaches ``level``.

    Only samples at or after ``earliest`` count; the result is not quantised to
    the time step.
    """
    idx = np.flatnonzero((t >= earliest) & (metric >= level))
    if idx.size == 0:
        return None
    i = int(idx[0])
    if i == 0 or t[i - 1] < earliest:
        return float(t[i])
    frac = (level - metric[i - 1]) / (metric[i] - metric[i - 1])
    return float(t[i - 1] + frac * (t[i] - t[i - 1]))


def time_to_corotation(result: SpinupResult, fraction: float = 0.9) -> float | None:
    """First time at or after the hull ramp that co-rotation reaches ``fraction``.

    Linearly interpolated between steps; ``None`` if it never does.
    """
    return _first_crossing(result.t, result.corotation, fraction, result.ramp_end)


def time_to_adjust(result: SpinupResult, fraction: float = 0.9) -> float | None:
    """First time the initial lag has decayed by ``fraction`` (sudden changes).

    Uses ``1 - sum(M u) / sum(M u_0)``, interpolated between steps; ``None`` if it
    never does.
    """
    lag = result.relative_angular_momentum
    return _first_crossing(result.t, 1.0 - lag / lag[0], fraction)


# --- numerics ----------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class _Grid:
    faces: NDArray[np.float64]
    centres: NDArray[np.float64]
    moment: NDArray[np.float64]
    quartic: NDArray[np.float64]  # (r_{i+1}^4 - r_i^4) / 4


def _grid(radius: float, n: int, stretch: float) -> _Grid:
    xi = np.linspace(0.0, 1.0, n + 1)
    faces = (
        radius * xi
        if stretch == 0.0
        else radius * np.tanh(stretch * xi) / math.tanh(stretch)
    )
    return _Grid(
        faces=faces,
        centres=0.5 * (faces[:-1] + faces[1:]),
        moment=np.diff(faces**3) / 3.0,
        quartic=np.diff(faces**4) / 4.0,
    )


@dataclass(frozen=True, slots=True)
class _Operator:
    lower: NDArray[np.float64]
    diag: NDArray[np.float64]
    upper: NDArray[np.float64]
    wall: NDArray[np.float64]  # floor torque coefficients: wall rate = -wall . u


@dataclass(frozen=True, slots=True)
class _Model:
    grid: _Grid
    params: SpinupParams
    half_length: float
    hull: HullSpin

    def sink_rate(self, t: float) -> float:
        """End-cap Ekman rate ``sqrt(nu Omega) / (L/2)`` in 1/s (0 at Omega = 0)."""
        if not self.params.endcap:
            return 0.0
        omega = abs(self.hull.omega(t))
        return math.sqrt(self.params.nu_cap * omega) / self.half_length

    def viscosity_faces(
        self, u: NDArray[np.float64], omega: float
    ) -> NDArray[np.float64]:
        """Face viscosities blended smoothly by the Rayleigh discriminant.

        ``q = -(2 Omega + u/r + du/dr)`` is positive where ``Gamma = r (Omega r + u)``
        falls outward; the weight of ``nu_unstable`` is ``(1 + tanh(q / q0)) / 2``.
        """
        g = self.grid
        p = self.params
        if p.nu_high == p.nu_stable:
            return np.full(len(g.faces), p.nu_stable)
        q = np.zeros(len(g.faces))
        du = np.diff(u) / np.diff(g.centres)
        u_face = 0.5 * (u[1:] + u[:-1])
        spin_in = omega + u_face / g.faces[1:-1]  # inertial rotation rate of the air
        q[1:-1] = -(2.0 * spin_in + du) / (
            2.0 * abs(spin_in) + _BLEND_FLOOR * abs(omega)
        )
        if p.floor_drag is None:  # no-slip wall: u = 0 at r = R
            q[-1] = -(2.0 * omega - u[-1] / (g.faces[-1] - g.centres[-1])) / (
                2.0 * abs(omega) + _BLEND_FLOOR * abs(omega) + _TINY
            )
        weight = 0.5 * (1.0 + np.tanh(q / _BLEND_WIDTH))
        return p.nu_stable + (p.nu_high - p.nu_stable) * weight

    def operator(self, nu: NDArray[np.float64], u: NDArray[np.float64]) -> _Operator:
        """Viscous plus floor tridiagonal operator from the wind ``u``.

        Row ``i`` acting on a wind is the net stress flux into cell ``i`` minus the
        floor drag. The floor torque is ``-wall . u``.
        """
        g = self.grid
        rc, rf = g.centres, g.faces
        k = np.zeros(len(rf))
        k[1:-1] = rf[1:-1] ** 3 * nu[1:-1] / np.diff(rc)
        lower = k[:-1] / np.roll(rc, 1)
        upper = k[1:] / np.roll(rc, -1)
        diag = -(k[:-1] + k[1:]) / rc
        wall = np.zeros(len(rc))
        floor = self.params.floor_drag
        if floor is None:  # no-slip: stress flux into the wall
            wall[-1] = rf[-1] ** 3 * nu[-1] / (rf[-1] - rc[-1]) / rc[-1]
        else:  # drag layer: body force C_D |u_r| / z_r in the bottom z_r metres
            z_r = self.params.drag_height
            u_ref = float(np.interp(rf[-1] - z_r, rc, u))
            overlap = np.clip((rf[1:] - (rf[-1] - z_r)) / np.diff(rf), 0.0, 1.0)
            wall = floor * abs(u_ref) / z_r * overlap * g.moment
        diag = diag - wall
        lower[0] = 0.0
        upper[-1] = 0.0
        return _Operator(lower, diag, upper, wall)


def _thomas(
    lower: NDArray[np.float64],
    diag: NDArray[np.float64],
    upper: NDArray[np.float64],
    rhs: NDArray[np.float64],
) -> NDArray[np.float64]:
    """Solve a tridiagonal system (``lower[0]`` and ``upper[-1]`` unused)."""
    n = len(diag)
    a, b, c, d = lower.tolist(), diag.tolist(), upper.tolist(), rhs.tolist()
    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] -= w * c[i - 1]
        d[i] -= w * d[i - 1]
    x = [0.0] * n
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        x[i] = (d[i] - c[i] * x[i + 1]) / b[i]
    return np.array(x)


def _apply(op: _Operator, u: NDArray[np.float64]) -> NDArray[np.float64]:
    out = op.diag * u
    out[1:] += op.lower[1:] * u[:-1]
    out[:-1] += op.upper[:-1] * u[1:]
    return out


def simulate(  # noqa: PLR0913
    cyl: Cylinder,
    params: SpinupParams,
    hull: HullSpin,
    u0: NDArray[np.float64],
    t_end: float,
    *,
    dt0: float = 1.0,
    growth: float = 1.05,
    dt_max: float = 600.0,
    stop: Stop | None = None,
) -> SpinupResult:
    """March the model from ``u0`` to ``t_end`` (or until ``stop(t, u)``).

    Steps start at ``dt0`` and grow geometrically to ``dt_max``. The viscosity
    regime and floor drag are Picard-iterated (three passes) within each step.
    With ``theta < 1`` the first two steps are backward Euler.

    Raises:
        ValueError: for invalid time stepping or if the step budget is exceeded.
    """
    positive("t_end", t_end)
    positive("dt0", dt0)
    positive("dt_max", dt_max)
    if growth < 1.0 or not math.isfinite(growth):
        msg = f"growth must be at least 1, got {growth!r}"
        raise ValueError(msg)
    if t_end / dt_max > _MAX_STEPS:
        msg = "step budget exceeded; increase dt_max or shorten t_end"
        raise ValueError(msg)
    grid = _grid(cyl.radius, params.n_cells, params.stretch)
    model = _Model(grid, params, cyl.length / 2.0, hull)
    u = np.array(u0, dtype=np.float64)
    if u.shape != grid.centres.shape:
        msg = f"u0 must have shape {grid.centres.shape}, got {u.shape}"
        raise ValueError(msg)

    ts, us, oms = [0.0], [u], [hull.omega(0.0)]
    dts: list[float] = []
    wall_r: list[float] = []
    cap_r: list[float] = []
    hull_r: list[float] = []
    t, dt = 0.0, dt0
    mass = grid.moment
    while t < t_end and not (stop is not None and stop(t, u)):
        step = min(dt, t_end - t)
        theta = 1.0 if len(dts) < _STARTUP_STEPS else params.theta
        om0, om1 = hull.omega(t), hull.omega(t + step)
        c0, c1 = model.sink_rate(t), model.sink_rate(t + step)
        alpha_bar = (om1 - om0) / step
        source = -alpha_bar * grid.quartic
        u_iter = u
        for _ in range(_PICARD_ITERATIONS):
            op = model.operator(model.viscosity_faces(u_iter, om1), u_iter)
            explicit = _apply(op, u) - c0 * mass * u
            rhs = mass * u / step + (1.0 - theta) * explicit + source
            u_iter = _thomas(
                -theta * op.lower,
                mass / step - theta * (op.diag - c1 * mass),
                -theta * op.upper,
                rhs,
            )
        u_new = u_iter
        mixed = theta * u_new + (1.0 - theta) * u
        wall_r.append(-float(op.wall @ mixed))
        cap_r.append(
            -float(np.sum(mass * (theta * c1 * u_new + (1.0 - theta) * c0 * u)))
        )
        hull_r.append(-alpha_bar * float(grid.quartic.sum()))
        u, t = u_new, t + step
        dts.append(step)
        ts.append(t)
        us.append(u)
        oms.append(om1)
        dt = min(dt * growth, dt_max)

    return SpinupResult(
        t=np.array(ts),
        radius=grid.centres,
        faces=grid.faces,
        moment=grid.moment,
        u=np.array(us),
        omega=np.array(oms),
        budget=Budget(
            np.array(dts), np.array(wall_r), np.array(cap_r), np.array(hull_r)
        ),
        wall_wind_is_zero=params.floor_drag is None,
        ramp_end=0.0 if hull.ramp_time is None else hull.ramp_time,
        drag_height=params.drag_height,
        endcap=params.endcap,
    )


# --- scenarios ---------------------------------------------------------------------


def hull_spin_up(
    cyl: Cylinder,
    params: SpinupParams,
    alpha: float,
    *,
    fraction: float = 0.9,
    t_max: float = 1.0e7,
    dt_max: float = 600.0,
) -> SpinupResult:
    """Hull spun from rest to ``cyl.omega`` at steady ``alpha``; air starts at rest.

    The run stops once the hull has finished its ramp and the air has reached
    ``fraction`` of co-rotation, or at ``t_max``.
    """
    positive("alpha", alpha)
    ramp = cyl.omega / alpha
    hull = HullSpin(0.0, alpha, ramp)
    quarter = cyl.radius**4 / 4.0
    grid = _grid(cyl.radius, params.n_cells, params.stretch)

    def done(t: float, u: NDArray[np.float64]) -> bool:
        corot = 1.0 + float(u @ grid.moment) / (cyl.omega * quarter)
        return t >= ramp and corot >= fraction

    zeros = np.zeros(params.n_cells)
    return simulate(cyl, params, hull, zeros, t_max, dt_max=dt_max, stop=done)


def sudden_change(
    cyl: Cylinder,
    params: SpinupParams,
    delta_omega: float,
    *,
    fraction: float = 0.9,
    t_max: float = 1.0e7,
    dt_max: float = 600.0,
) -> SpinupResult:
    """The hull's spin rate jumps by ``delta_omega`` while the air keeps its motion.

    Positive is spin-up (air lags, ``u < 0``); negative is spin-down. The run
    stops once ``fraction`` of the initial lag has gone, or at ``t_max``.

    With end caps the sink rate scales as ``sqrt(Omega)`` of the *new* hull rate,
    so with equal viscosities a spin-down (``0.9 Omega``) adjusts more slowly than
    a spin-up (``1.1 Omega``): the asymmetry is the Ekman time, not instability.
    """
    finite("delta_omega", delta_omega)
    if delta_omega == 0.0 or cyl.omega + delta_omega <= 0.0:
        msg = "delta_omega must be nonzero and leave the hull spinning forward"
        raise ValueError(msg)
    hull = HullSpin(cyl.omega + delta_omega)
    grid = _grid(cyl.radius, params.n_cells, params.stretch)
    u0 = -delta_omega * grid.centres
    lag0 = float(u0 @ grid.moment)

    def done(_t: float, u: NDArray[np.float64]) -> bool:
        return 1.0 - float(u @ grid.moment) / lag0 >= fraction

    return simulate(cyl, params, hull, u0, t_max, dt_max=dt_max, stop=done)


def steady_lag(cyl: Cylinder, params: SpinupParams, alpha: float) -> LagWind:
    """Steady wind lag while the hull accelerates slowly at ``alpha``.

    Solves ``0 = viscous + sink - drag - r alpha`` by a fixed number of damped
    fixed-point iterations on the viscosity regime and the floor drag. With
    constant ``nu``, no end caps and a no-slip floor this is
    ``u = alpha r (r^2 - R^2) / (8 nu)``.
    """
    finite("alpha", alpha)
    hull = HullSpin(cyl.omega)
    grid = _grid(cyl.radius, params.n_cells, params.stretch)
    model = _Model(grid, params, cyl.length / 2.0, hull)
    sink = model.sink_rate(0.0)
    source = alpha * grid.quartic
    u = np.zeros(params.n_cells)
    for _ in range(_LAG_ITERATIONS):
        op = model.operator(model.viscosity_faces(u, cyl.omega), u)
        u_new = _thomas(op.lower, op.diag - sink * grid.moment, op.upper, source)
        u = 0.5 * (u + u_new)
    return LagWind(grid.centres, u, grid.faces, params.floor_drag is None)
