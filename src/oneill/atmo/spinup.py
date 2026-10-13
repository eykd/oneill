"""L1: axisymmetric spin-up and spin-down of the air, one radial dimension.

``u(r, t)`` is the azimuthal wind relative to the hull. Angular momentum of the
air obeys, with the hull spinning at ``Omega(t)`` and accelerating at ``alpha``::

    du/dt = (1/r^2) d/dr( r^3 nu d(u/r)/dr ) - u/tau_E - r alpha

The viscous term is written in stress form, so a solid-body ``u ~ r`` feels no
viscous force (for constant ``nu`` it equals ``nu (u'' + u'/r - u/r^2)``).
``u/tau_E`` is the Ekman pumping of the end caps (Greenspan & Howard 1963;
``tau_E = (L/2) / sqrt(nu Omega)``, ``docs/research/atmosphere-literature.md``
Q2) and can be switched off. ``r alpha`` is the apparent force from the
accelerating hull.

Eddy viscosity has two regimes (Q4): the larger ``nu_unstable`` wherever the
absolute angular momentum ``Gamma = r (Omega r + u)`` decreases outward
(``d Gamma^2 / dr < 0``, Rayleigh), the smaller ``nu_stable`` elsewhere. Spin-up
from the floor is stable and slow; spin-down is unstable and fast (Kaiser 2020).

The floor is either no-slip, or a quadratic stress ``nu r d(u/r)/dr = -C_D |u| u``
on the first cell (``C_D`` dimensionless; ``0`` is free slip). The scheme is a
finite-volume theta method (Crank-Nicolson, with two backward-Euler starting steps)
on a grid clustered toward the floor. It conserves angular momentum exactly:
:class:`Budget` holds the wall torque and end-cap sink that close the books.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    from oneill.core.cylinder import Cylinder

    type Stop = Callable[[float, NDArray[np.float64]], bool]

_MAX_STEPS = 200_000
_LAG_ITERATIONS = 80
_STARTUP_STEPS = 2
_MIN_CELLS = 8


@dataclass(frozen=True, slots=True)
class SpinupParams:
    """Physical and numerical parameters of the 1-D model.

    Attributes:
        nu_stable: Viscosity (molecular or eddy) where angular momentum rises
            outward, m^2/s. It also sets ``tau_E``.
        nu_unstable: Viscosity where angular momentum falls outward, m^2/s;
            defaults to ``nu_stable``.
        floor_drag: ``None`` for a no-slip floor, else the quadratic drag
            coefficient ``C_D`` (``0`` for free slip).
        endcap: Include the end-cap Ekman pumping sink.
        n_cells: Radial cells.
        stretch: Grid clustering toward the floor (``0`` is uniform; the smallest
            cell is about ``R/n`` times ``beta / cosh^2(beta)``).
    """

    nu_stable: float
    nu_unstable: float | None = None
    floor_drag: float | None = None
    endcap: bool = True
    n_cells: int = 120
    stretch: float = 2.5

    def __post_init__(self) -> None:
        """Validate the viscosities, drag and grid size."""
        positive("nu_stable", self.nu_stable)
        if self.nu_unstable is not None:
            positive("nu_unstable", self.nu_unstable)
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

    @property
    def nu_high(self) -> float:
        """Viscosity used in the unstable regime."""
        return self.nu_stable if self.nu_unstable is None else self.nu_unstable


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
        wall: Torque exchanged with the floor.
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
    """Largest ``|u|`` over the run at ``height`` metres above the floor, in m/s."""
    return max(
        abs(
            _wind_at_height(
                result.radius, row, result.faces[-1], height, result.wall_wind_is_zero
            )
        )
        for row in result.u
    )


def time_to_corotation(result: SpinupResult, fraction: float = 0.9) -> float | None:
    """First time at or after the hull ramp that co-rotation reaches ``fraction``.

    Returns ``None`` if it never does.
    """
    ok = (result.t >= result.ramp_end) & (result.corotation >= fraction)
    idx = np.flatnonzero(ok)
    return None if idx.size == 0 else float(result.t[idx[0]])


def time_to_adjust(result: SpinupResult, fraction: float = 0.9) -> float | None:
    """First time the initial lag has decayed by ``fraction`` (sudden changes).

    Uses ``1 - sum(M u) / sum(M u_0)``; returns ``None`` if it never does.
    """
    lag = result.relative_angular_momentum
    done = 1.0 - lag / lag[0]
    idx = np.flatnonzero(done >= fraction)
    return None if idx.size == 0 else float(result.t[idx[0]])


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
class _Model:
    grid: _Grid
    params: SpinupParams
    half_length: float
    hull: HullSpin

    def sink_rate(self, t: float) -> float:
        """End-cap Ekman rate ``sqrt(nu Omega) / (L/2)`` in 1/s."""
        if not self.params.endcap:
            return 0.0
        omega = abs(self.hull.omega(t))
        return math.sqrt(self.params.nu_stable * omega) / self.half_length

    def viscosity_faces(
        self, u: NDArray[np.float64], omega: float
    ) -> NDArray[np.float64]:
        """Face viscosities chosen by the Rayleigh criterion on ``Gamma^2``."""
        g = self.grid
        wall_u = 0.0 if self.params.floor_drag is None else float(u[-1])
        gamma = g.centres * (omega * g.centres + u)
        r_wall = g.faces[-1]
        gamma_all = np.append(gamma, r_wall * (omega * r_wall + wall_u))
        unstable = np.zeros(len(g.faces), dtype=bool)
        unstable[1:] = gamma_all[1:] ** 2 < gamma_all[:-1] ** 2
        return np.where(unstable, self.params.nu_high, self.params.nu_stable)

    def operator(
        self, nu: NDArray[np.float64], wall_u: float
    ) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], float]:
        """Viscous tridiagonal operator and the wall-flux coefficient.

        Returns ``(lower, diag, upper, wall_coeff)`` such that row ``i`` of the
        operator acting on ``u`` is the net stress flux into cell ``i``, and the
        wall flux is ``-wall_coeff * u[-1]``.
        """
        g = self.grid
        rc, rf = g.centres, g.faces
        k = np.zeros(len(rf))
        k[1:-1] = rf[1:-1] ** 3 * nu[1:-1] / np.diff(rc)
        if self.params.floor_drag is None:
            k_wall = rf[-1] ** 3 * nu[-1] / (rf[-1] - rc[-1])
            wall_coeff = k_wall / rc[-1]
        else:
            wall_coeff = rf[-1] ** 2 * self.params.floor_drag * abs(wall_u)
        lower = k[:-1] / np.roll(rc, 1)
        upper = k[1:] / np.roll(rc, -1)
        diag = -(k[:-1] + k[1:]) / rc
        diag[-1] -= wall_coeff
        lower[0] = 0.0
        upper[-1] = 0.0
        return lower, diag, upper, wall_coeff


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


def _apply(
    lower: NDArray[np.float64],
    diag: NDArray[np.float64],
    upper: NDArray[np.float64],
    u: NDArray[np.float64],
) -> NDArray[np.float64]:
    out = diag * u
    out[1:] += lower[1:] * u[:-1]
    out[:-1] += upper[:-1] * u[1:]
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

    Steps start at ``dt0`` and grow geometrically to ``dt_max``; the first two
    are backward Euler to damp the start-up transient.

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
    model = _Model(
        _grid(cyl.radius, params.n_cells, params.stretch),
        params,
        cyl.length / 2.0,
        hull,
    )
    g = model.grid
    u = np.array(u0, dtype=np.float64)
    if u.shape != g.centres.shape:
        msg = f"u0 must have shape {g.centres.shape}, got {u.shape}"
        raise ValueError(msg)

    ts, us, oms = [0.0], [u], [hull.omega(0.0)]
    dts: list[float] = []
    wall_r: list[float] = []
    cap_r: list[float] = []
    hull_r: list[float] = []
    t, dt = 0.0, dt0
    while t < t_end and not (stop is not None and stop(t, u)):
        step = min(dt, t_end - t)
        theta = 1.0 if len(dts) < _STARTUP_STEPS else 0.5
        om0, om1 = hull.omega(t), hull.omega(t + step)
        nu = model.viscosity_faces(u, om0)
        lo, di, up, wall_coeff = model.operator(nu, float(u[-1]))
        c0, c1 = model.sink_rate(t), model.sink_rate(t + step)
        alpha_bar = (om1 - om0) / step
        source = -alpha_bar * g.quartic
        mass = g.moment
        explicit = _apply(lo, di, up, u) - c0 * mass * u
        rhs = mass * u / step + (1.0 - theta) * explicit + source
        u_new = _thomas(
            -theta * lo,
            mass / step - theta * (di - c1 * mass),
            -theta * up,
            rhs,
        )
        wall_r.append(-wall_coeff * (theta * u_new[-1] + (1.0 - theta) * u[-1]))
        cap_r.append(
            -float(np.sum(mass * (theta * c1 * u_new + (1.0 - theta) * c0 * u)))
        )
        hull_r.append(-alpha_bar * float(g.quartic.sum()))
        u, t = u_new, t + step
        dts.append(step)
        ts.append(t)
        us.append(u)
        oms.append(om1)
        dt = min(dt * growth, dt_max)

    return SpinupResult(
        t=np.array(ts),
        radius=g.centres,
        faces=g.faces,
        moment=g.moment,
        u=np.array(us),
        omega=np.array(oms),
        budget=Budget(
            np.array(dts), np.array(wall_r), np.array(cap_r), np.array(hull_r)
        ),
        wall_wind_is_zero=params.floor_drag is None,
        ramp_end=0.0 if hull.ramp_time is None else hull.ramp_time,
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

    Solves ``0 = viscous + sink - r alpha`` by a fixed number of fixed-point
    iterations on the viscosity regime and the quadratic floor drag. With
    constant ``nu``, no end caps and no-slip this is
    ``u = alpha r (r^2 - R^2) / (8 nu)``.
    """
    finite("alpha", alpha)
    hull = HullSpin(cyl.omega)
    model = _Model(
        _grid(cyl.radius, params.n_cells, params.stretch),
        params,
        cyl.length / 2.0,
        hull,
    )
    g = model.grid
    sink = model.sink_rate(0.0)
    source = alpha * g.quartic
    u = np.zeros(params.n_cells)
    for _ in range(_LAG_ITERATIONS):
        nu = model.viscosity_faces(u, cyl.omega)
        lo, di, up, _ = model.operator(nu, float(u[-1]))
        u_new = _thomas(lo, di - sink * g.moment, up, source)
        u = 0.5 * (u + u_new) if params.floor_drag else u_new
    return LagWind(g.centres, u, g.faces, params.floor_drag is None)
