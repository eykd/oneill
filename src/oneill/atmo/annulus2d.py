"""L3: 2-D Boussinesq convection in the (r, theta) annulus of a spinning habitat.

Design and derivations: ``docs/design/annulus2d.md``. The flow is invariant along
the axis. With ``u_r = psi_theta / r``, ``u_theta = -psi_r`` and ``zeta = -lap psi``::

    d zeta/dt + u.grad zeta = nu lap zeta + beta Omega^2 dT/dtheta - zeta/tau_E
                              - 2 Omega beta dT/dt           (Coriolis-buoyancy)
                              + curl(drag)
    d T/dt + u.grad T = kappa lap T                            (+ wall flux)
    d ubar/dt = -(1/r^2) d(r^2 <u_r u_theta>)/dr + F_cb + visc - ubar/tau_E + drag

``T`` is the potential-temperature anomaly. Gravity is centrifugal, ``g = Omega^2 r``,
pointing outward, so warm air feels a force toward the axis and the baroclinic
torque is ``+ beta Omega^2 dT/dtheta``. Coriolis cancels in 2-D incompressible flow.
The mean (m = 0) wind is advanced as ``ubar`` because ``psi_0`` has no boundary
value. Coriolis acts on the *mass flux*, whose divergence is ``beta dT/dt`` (a
closed fluid conserves mass), so the extended-Boussinesq Coriolis-buoyancy source is
``-2 Omega beta dT/dt`` in the vorticity equation, built from the model's own
temperature tendency (``F_cb`` is its mean part, anomalies about the area mean). It
vanishes in a steady state; there is no steady retrograde mean wind from it. It
moves angular momentum between relative motion and the moment of inertia, so
:meth:`Annulus.absolute_angular_momentum` changes only through wall and end-cap
torques.

Numerics: rfft in theta (2/3 dealiasing, an optional polar filter and an optional
exponential spectral filter), finite volume in r on a grid clustered toward the floor,
variable-step IMEX SBDF2 (diffusion and end-cap damping implicit, everything else
explicit, backward Euler for the first step), tridiagonal solves batched over Fourier
modes. Advection is in flux form: conservative in r (third-order upwind-biased face
values, zero flux through both walls) and spectral in theta, so the theta-mean heat
content changes only by the wall flux, exactly. The mean-wind Reynolds stress is also
a flux, so angular momentum telescopes. Energy is *not* conserved exactly. The
floor is the L1 canopy-drag layer; the vorticity wall condition is free slip
(no-slip, a three-point Taylor formula, for validation). Dirichlet walls use the
second-order three-point wall derivative, so the onset thresholds converge at
second order.
"""

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, NamedTuple

import numpy as np

from oneill.atmo.spinup import DEFAULT_FLOOR_DRAG, DRAG_HEIGHT
from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import ArrayLike, NDArray

    from oneill.atmo.gas import Gas
    from oneill.core.cylinder import Cylinder

    type Array = NDArray[np.float64]
    type CArray = NDArray[np.complex128]

DEFAULT_CFL = 0.15
"""Default CFL number: 0.4 left a ~4 % energy gain when an inviscid jet saturated."""

_STEP_RATIO_LIMIT = 1.0 + math.sqrt(2.0)
_MIN_CELLS = 8
_MIN_THETA = 8
_FILTER_ORDER = 8


@dataclass(frozen=True, slots=True)
class Forcing:
    """Surface heating of the outer wall (the floor).

    Attributes:
        land_flux: Peak land heat flux in W/m^2 (daytime, or constant if steady).
        night_flux: Land flux at the night peak in W/m^2 (negative); used with a
            ``period``.
        period: Day length in s for a sine day/night cycle; ``None`` is steady.
        strips: Number of land strips (equal number of window strips).
        window: ``"neutral"`` (window flux balances the land flux, zero net flux),
            ``"zero"`` (windows do not exchange heat), or ``"uniform"`` (no
            strips: the land flux everywhere, an axisymmetric case).
        phase: Azimuth in rad of the start of the first land strip.
        edge_width: Width ``w`` in rad of a ``tanh(delta_theta / w)`` smoothing of the
            strip edges (``0``, the default, keeps the step function). A width of
            about two grid spacings avoids exciting the grid scale in theta.
    """

    land_flux: float = 100.0
    night_flux: float = -30.0
    period: float | None = None
    strips: int = 3
    window: str = "neutral"
    phase: float = 0.0
    edge_width: float = 0.0

    def __post_init__(self) -> None:
        """Validate the strip count, period, window mode and edge width."""
        if not (math.isfinite(self.edge_width) and self.edge_width >= 0.0):
            msg = f"edge_width must be non-negative and finite, got {self.edge_width!r}"
            raise ValueError(msg)
        finite("land_flux", self.land_flux)
        finite("night_flux", self.night_flux)
        if self.strips < 1:
            msg = f"strips must be at least 1, got {self.strips}"
            raise ValueError(msg)
        if self.period is not None:
            positive("period", self.period)
        if self.window not in ("neutral", "zero", "uniform"):
            msg = f"window must be 'neutral', 'zero' or 'uniform', got {self.window!r}"
            raise ValueError(msg)

    def land_level(self, t: float) -> float:
        """Land flux in W/m^2 at time ``t``."""
        if self.period is None:
            return self.land_flux
        s = math.sin(2.0 * math.pi * t / self.period)
        return self.land_flux * s if s > 0.0 else self.night_flux * -s

    def flux(self, theta: ArrayLike, t: float) -> Array:
        """Wall heat flux (W/m^2, into the air) at azimuths ``theta``."""
        theta = np.asarray(theta, dtype=np.float64)
        level = self.land_level(t)
        if self.window == "uniform":
            return np.full_like(theta, level)
        if self.edge_width > 0.0:
            wave = np.sin(self.strips * (theta - self.phase))
            land = 0.5 * (1.0 + np.tanh(wave / (self.strips * self.edge_width)))
        else:
            sector = np.floor((theta - self.phase) / (math.pi / self.strips))
            land = (sector.astype(int) % 2 == 0).astype(np.float64)
        window = 0.0
        if self.window == "neutral" and np.sum(1.0 - land) > 0.0:
            # exactly zero net flux on the grid
            window = -level * float(np.sum(land) / np.sum(1.0 - land))
        return np.asarray(land * level + (1.0 - land) * window)


@dataclass(frozen=True, slots=True)
class AnnulusConfig:
    """Parameters of the annulus model.

    Attributes:
        radius: Outer radius ``R`` in m.
        length: Habitat length in m (sets the end-cap rate).
        omega: Spin rate in rad/s.
        r_in_fraction: Inner wall radius over ``R``.
        n_r: Radial cells.
        n_theta: Azimuthal points (modes up to ``n_theta / 3`` after dealiasing).
        stretch: Radial clustering toward the floor (0 uniform).
        nu: Eddy viscosity in m^2/s.
        prandtl: Eddy Prandtl number (``kappa = nu / prandtl``).
        beta: Thermal expansion coefficient in 1/K.
        rho_cp: ``rho c_p`` in J/(m^3 K), converts heat flux to temperature flux.
        floor_drag: Canopy-drag coefficient ``C_D`` (``0`` disables the drag).
        drag_height: Drag layer thickness and reference height in m.
        endcap: Include end-cap Ekman damping ``1/tau_E``.
        wall: Vorticity wall condition, ``"free"`` or ``"noslip"``.
        advection: Include the nonlinear advection terms.
        coriolis_buoyancy: Include the mass-consistent Coriolis-buoyancy force.
        polar_filter: Taper the retained azimuthal modes toward the inner wall.
        tau_relax: Relaxation time of the area-mean temperature in s, or ``None``.
        outer_t: Outer thermal wall, ``"flux"`` (forced by ``forcing``) or ``"fixed"``.
        inner_t: Inner thermal wall, ``"flux"`` (insulating unless ``flux_inner``) or
            ``"fixed"``.
        t_outer: Outer wall temperature in K for ``"fixed"``.
        t_inner: Inner wall temperature in K for ``"fixed"``.
        flux_inner: ``rho c_p kappa dT/dr`` at the inner wall in W/m^2 for ``"flux"``:
            positive when ``T`` increases outward, i.e. heat flows *out of* the fluid
            into the inner wall (the outer-wall forcing is the opposite: positive
            heats the fluid).
        forcing: Surface heating pattern.
        modes: If given, only these azimuthal wavenumbers (and 0) are retained
            (linear stability tests).
        filter_strength: Strength ``alpha`` of an exponential spectral filter
            ``exp(-alpha (m / m_top)^8)`` applied to ``zeta`` and ``T`` every
            step (0 disables it); a scale-selective stand-in for the unresolved
            small scales at the small eddy viscosities of the production runs.
    """

    radius: float = 3200.0
    length: float = 20000.0
    omega: float = 0.05535863189241584
    r_in_fraction: float = 0.1
    n_r: int = 64
    n_theta: int = 256
    stretch: float = 1.8
    nu: float = 1.0
    prandtl: float = 1.0
    beta: float = 1.0 / 290.0
    rho_cp: float = 603.0
    floor_drag: float = DEFAULT_FLOOR_DRAG
    drag_height: float = DRAG_HEIGHT
    endcap: bool = True
    wall: str = "free"
    advection: bool = True
    coriolis_buoyancy: bool = True
    polar_filter: bool = True
    tau_relax: float | None = None
    outer_t: str = "flux"
    inner_t: str = "flux"
    t_outer: float = 0.0
    t_inner: float = 0.0
    flux_inner: float = 0.0
    forcing: Forcing = field(default_factory=Forcing)
    modes: tuple[int, ...] | None = None
    filter_strength: float = 0.0

    def __post_init__(self) -> None:
        """Validate the parameters."""
        for name in ("radius", "length", "omega", "nu", "prandtl", "beta", "rho_cp"):
            positive(name, getattr(self, name))
        positive("drag_height", self.drag_height)
        finite("floor_drag", self.floor_drag)
        if self.floor_drag < 0.0:
            msg = f"floor_drag must be non-negative, got {self.floor_drag!r}"
            raise ValueError(msg)
        if not 0.0 < self.r_in_fraction < 1.0:
            msg = f"r_in_fraction must lie in (0, 1), got {self.r_in_fraction!r}"
            raise ValueError(msg)
        if self.n_r < _MIN_CELLS or self.n_theta < _MIN_THETA:
            msg = f"need n_r >= {_MIN_CELLS} and n_theta >= {_MIN_THETA}"
            raise ValueError(msg)
        if self.wall not in ("free", "noslip"):
            msg = f"wall must be 'free' or 'noslip', got {self.wall!r}"
            raise ValueError(msg)
        for name in (self.outer_t, self.inner_t):
            if name not in ("flux", "fixed"):
                msg = f"thermal wall must be 'flux' or 'fixed', got {name!r}"
                raise ValueError(msg)
        if not (math.isfinite(self.stretch) and self.stretch >= 0.0):
            msg = f"stretch must be non-negative and finite, got {self.stretch!r}"
            raise ValueError(msg)
        if self.tau_relax is not None:
            positive("tau_relax", self.tau_relax)

    @classmethod
    def from_habitat(cls, cyl: Cylinder, gas: Gas, **kwargs: object) -> AnnulusConfig:
        """Config for ``cyl`` and ``gas`` (``rho c_p`` and ``beta`` from the gas)."""
        return cls(
            radius=cyl.radius,
            length=cyl.length,
            omega=cyl.omega,
            beta=1.0 / gas.temperature,
            rho_cp=gas.rho_rim * gas.cp,
            **kwargs,  # type: ignore[arg-type]
        )

    @property
    def kappa(self) -> float:
        """Eddy thermal diffusivity in m^2/s."""
        return self.nu / self.prandtl

    @property
    def endcap_rate(self) -> float:
        """End-cap damping rate ``1/tau_E = sqrt(nu Omega) / (L/2)`` in 1/s."""
        if not self.endcap:
            return 0.0
        return math.sqrt(self.nu * self.omega) / (self.length / 2.0)


@dataclass(eq=False)
class AnnulusState:
    """Evolving fields (spectral in theta, cell-centred in r).

    ``zeta`` and ``temp`` have shape ``(n_r, n_theta // 2 + 1)`` (``norm="forward"``
    rfft: column 0 is the theta-mean); ``ubar`` is the theta-mean azimuthal wind.
    """

    t: float
    zeta: CArray
    temp: CArray
    ubar: Array
    previous: tuple[CArray, CArray, Array] | None = None
    last_terms: tuple[CArray, CArray, Array] | None = None
    last_dt: float | None = None
    steps: int = 0


@dataclass(frozen=True, slots=True)
class Terms:
    """Pieces of the mean-wind forcing (rates in m/s^2) for budget checks."""

    reynolds: Array
    coriolis_buoyancy: Array
    drag: Array


def _d_weights(
    centres: Array, r_in: float, radius: float
) -> tuple[Array, Array, Array]:
    """Second-order nonuniform first-derivative matrix on cell centres.

    Returns ``(D, d_in, d_out)`` so the derivative is ``D f + d_in f_in + d_out f_out``
    where ``f_in`` and ``f_out`` are the wall values.
    """
    n = len(centres)
    x = np.concatenate([[r_in], centres, [radius]])
    d = np.zeros((n, n))
    d_in = np.zeros(n)
    d_out = np.zeros(n)
    for i in range(n):
        h1 = x[i + 1] - x[i]
        h2 = x[i + 2] - x[i + 1]
        w_minus = -h2 / (h1 * (h1 + h2))
        w_zero = (h2 - h1) / (h1 * h2)
        w_plus = h1 / (h2 * (h1 + h2))
        if i == 0:
            d_in[0] = w_minus
        else:
            d[i, i - 1] = w_minus
        d[i, i] = w_zero
        if i == n - 1:
            d_out[i] = w_plus
        else:
            d[i, i + 1] = w_plus
    return d, d_in, d_out


def _lagrange_weights(xs: Array, centre: int) -> Array:
    """Weights ``w`` with ``f'(xs[centre]) ~ sum w_k f(xs[k])`` through ``xs``."""
    n = len(xs)
    x0 = xs[centre]
    w = np.zeros(n)
    for k in range(n):
        if k == centre:
            w[k] = sum(1.0 / (x0 - xs[j]) for j in range(n) if j != centre)
            continue
        num = 1.0
        den = 1.0
        for j in range(n):
            if j != k:
                den *= xs[k] - xs[j]
                if j != centre:
                    num *= x0 - xs[j]
        w[k] = num / den
    return w


def _interp_weights(xs: Array, x0: float) -> Array:
    """Lagrange weights ``w`` with ``f(x0) ~ sum w_k f(xs[k])``."""
    n = len(xs)
    w = np.ones(n)
    for k in range(n):
        for j in range(n):
            if j != k:
                w[k] *= (x0 - xs[j]) / (xs[k] - xs[j])
    return w


def _face_weights(
    centres: Array, faces: Array, r_in: float, radius: float
) -> tuple[Array, Array, Array]:
    """Face interpolation matrices on ``[f_in, f, f_out]`` for the interior faces.

    Returns ``(upwind_plus, upwind_minus, linear)``, each ``(n-1, n+2)``: the first
    two are third-order upwind-biased (flow toward larger / smaller ``r``), the last
    is the linear interpolation of the two adjacent cells.
    """
    n = len(centres)
    x = np.concatenate([[r_in], centres, [radius]])
    plus = np.zeros((n - 1, n + 2))
    minus = np.zeros((n - 1, n + 2))
    linear = np.zeros((n - 1, n + 2))
    for j in range(1, n):
        xf = faces[j]
        ip = [j - 1, j, j + 1]
        im = [j, j + 1, j + 2]
        plus[j - 1, ip] = _interp_weights(x[ip], xf)
        minus[j - 1, im] = _interp_weights(x[im], xf)
        il = [j, j + 1]
        linear[j - 1, il] = _interp_weights(x[il], xf)
    return plus, minus, linear


def _wall_weights(x0: float, x1: float, x2: float) -> tuple[float, float, float]:
    """Weights of ``f'(x0) ~ c0 f0 + c1 f1 + c2 f2`` (quadratic through three nodes)."""
    c0 = (2.0 * x0 - x1 - x2) / ((x0 - x1) * (x0 - x2))
    c1 = (x0 - x2) / ((x1 - x0) * (x1 - x2))
    c2 = (x0 - x1) / ((x2 - x0) * (x2 - x1))
    return c0, c1, c2


def _wall_derivative(
    x0: float, x1: float, x2: float, f0: CArray, f1: CArray, f2: CArray
) -> CArray:
    """Derivative at ``x0`` of the quadratic through three nodes."""
    c0 = (2.0 * x0 - x1 - x2) / ((x0 - x1) * (x0 - x2))
    c1 = (x0 - x2) / ((x1 - x0) * (x1 - x2))
    c2 = (x0 - x1) / ((x2 - x0) * (x2 - x1))
    return np.asarray(c0 * f0 + c1 * f1 + c2 * f2)


type Factors = tuple[CArray, CArray, CArray]


class _Factors(NamedTuple):
    zeta: Factors
    temp: Factors
    mean: Factors


def _factor(
    lower: CArray, diag: CArray, upper: CArray
) -> tuple[CArray, CArray, CArray]:
    """Thomas factorisation along axis 0 for many columns at once."""
    n = diag.shape[0]
    denom = np.empty_like(diag)
    c_prime = np.empty_like(upper)
    denom[0] = diag[0]
    c_prime[0] = upper[0] / denom[0]
    for i in range(1, n):
        denom[i] = diag[i] - lower[i] * c_prime[i - 1]
        c_prime[i] = upper[i] / denom[i]
    return lower, denom, c_prime


def _solve(factors: tuple[CArray, CArray, CArray], rhs: CArray) -> CArray:
    lower, denom, c_prime = factors
    n = denom.shape[0]
    y = np.array(rhs, dtype=np.complex128)
    y[0] = y[0] / denom[0]
    for i in range(1, n):
        y[i] = (y[i] - lower[i] * y[i - 1]) / denom[i]
    for i in range(n - 2, -1, -1):
        y[i] = y[i] - c_prime[i] * y[i + 1]
    return y


class Annulus:
    """The annulus model: grid, operators, time stepping and diagnostics."""

    def __init__(self, config: AnnulusConfig) -> None:
        """Build the grid, derivative matrices, masks and radial operators."""
        self.config = config
        c = config
        self.r_in = c.r_in_fraction * c.radius
        n = c.n_r
        xi = np.linspace(0.0, 1.0, n + 1)
        s = xi if c.stretch == 0.0 else np.tanh(c.stretch * xi) / math.tanh(c.stretch)
        self.faces: Array = self.r_in + (c.radius - self.r_in) * s
        self.r: Array = 0.5 * (self.faces[:-1] + self.faces[1:])
        self.width: Array = np.diff(self.faces)
        self.moment: Array = np.diff(self.faces**3) / 3.0
        self.area: Array = np.diff(self.faces**2) / 2.0
        self.n_theta = c.n_theta
        self.k = c.n_theta // 2 + 1
        self.theta: Array = 2.0 * np.pi * np.arange(c.n_theta) / c.n_theta
        self.m: Array = np.arange(self.k, dtype=np.float64)
        self.d_mat, self.d_in, self.d_out = _d_weights(self.r, self.r_in, c.radius)
        self.f_plus, self.f_minus, self.f_linear = _face_weights(
            self.r, self.faces, self.r_in, c.radius
        )
        self._lap_coefficients()
        self._build_mask()
        self.drag_weight = np.clip(
            (self.faces[1:] - (c.radius - c.drag_height)) / self.width, 0.0, 1.0
        )
        near = int(np.searchsorted(self.r, c.radius - c.drag_height) - 1)
        self._ref_index = int(np.clip(near, 0, n - 2))
        r0, r1 = self.r[self._ref_index], self.r[self._ref_index + 1]
        self._ref_weight = float(
            np.clip((c.radius - c.drag_height - r0) / (r1 - r0), 0.0, 1.0)
        )
        self._poisson = self._poisson_factors()
        self._cache: dict[tuple[float, float], _Factors] = {}

    # --- setup ----------------------------------------------------------------
    def _lap_coefficients(self) -> None:
        r, rf = self.r, self.faces
        face_a = rf[1:-1] / np.diff(r)
        self._lo0 = np.zeros(self.config.n_r)
        self._up0 = np.zeros(self.config.n_r)
        self._lo0[1:] = face_a / self.area[1:]
        self._up0[:-1] = face_a / self.area[:-1]
        # second-order Dirichlet closure: wall derivative from (wall, c0, c1)
        self._d_in = _wall_weights(self.r_in, r[0], r[1])
        self._d_out = _wall_weights(self.config.radius, r[-1], r[-2])
        self.w_in = -self.r_in * self._d_in[0] / self.area[0]
        self.w_out = self.config.radius * self._d_out[0] / self.area[-1]

    def bands(
        self, *, dirichlet_in: bool, dirichlet_out: bool
    ) -> tuple[Array, Array, Array]:
        """Radial Laplacian ``(1/r) d(r f')/dr`` as ``(lower, diag, upper)`` bands.

        A wall is Neumann (zero flux) unless Dirichlet, in which case the wall
        derivative is the second-order three-point one-sided formula and the wall
        value enters through ``w_in`` / ``w_out`` (rate per unit wall value).
        """
        lo, up = self._lo0.copy(), self._up0.copy()
        di = -(lo + up)
        if dirichlet_in:
            di[0] += -self.r_in * self._d_in[1] / self.area[0]
            up[0] += -self.r_in * self._d_in[2] / self.area[0]
        if dirichlet_out:
            di[-1] += self.config.radius * self._d_out[1] / self.area[-1]
            lo[-1] += self.config.radius * self._d_out[2] / self.area[-1]
        return lo, di, up

    def _t_bands(self) -> tuple[Array, Array, Array]:
        c = self.config
        return self.bands(
            dirichlet_in=c.inner_t == "fixed", dirichlet_out=c.outer_t == "fixed"
        )

    def apply_laplacian(
        self,
        hat: CArray,
        bands: tuple[Array, Array, Array],
        f_in: CArray,
        f_out: CArray,
    ) -> CArray:
        """Apply ``bands`` (and the azimuthal ``-m^2/r^2`` term) to ``hat``.

        ``f_in`` / ``f_out`` are the Dirichlet wall values (zero for a Neumann wall,
        where the coefficient is not used).
        """
        lo, di, up = bands
        out = di[:, None] * hat - (self.m[None, :] ** 2 / self.r[:, None] ** 2) * hat
        out[1:] += lo[1:, None] * hat[:-1]
        out[:-1] += up[:-1, None] * hat[1:]
        out[0] += self.w_in * f_in
        out[-1] += self.w_out * f_out
        return np.asarray(out)

    def _build_mask(self) -> None:
        c = self.config
        m_trunc = c.n_theta // 3
        cut = np.full(c.n_r, float(m_trunc))
        if c.polar_filter:
            cut = np.maximum(4.0, np.floor(m_trunc * self.r / c.radius))
        mask = self.m[None, :] <= cut[:, None]
        if c.modes is not None:
            keep = np.isin(self.m, (0, *c.modes))
            mask = mask & keep[None, :]
        self.mask: NDArray[np.bool_] = mask
        top = np.maximum(
            1.0, np.minimum(float(m_trunc), cut) if c.polar_filter else float(m_trunc)
        )
        ramp = (
            self.m[None, :] / np.broadcast_to(top, (c.n_r,))[:, None]
        ) ** _FILTER_ORDER
        self.damping: Array = np.exp(-c.filter_strength * ramp)

    def _poisson_factors(self) -> tuple[CArray, CArray, CArray]:
        lo, di, up = self.bands(dirichlet_in=True, dirichlet_out=True)
        diag = di[:, None] - self.m[None, :] ** 2 / self.r[:, None] ** 2
        return _factor(
            np.repeat(lo[:, None], self.k, axis=1).astype(np.complex128),
            diag.astype(np.complex128),
            np.repeat(up[:, None], self.k, axis=1).astype(np.complex128),
        )

    # --- transforms -----------------------------------------------------------
    def to_phys(self, hat: CArray) -> Array:
        """Physical values on the ``(n_r, n_theta)`` grid."""
        return np.asarray(np.fft.irfft(hat, n=self.n_theta, axis=1, norm="forward"))

    def to_hat(self, phys: Array) -> CArray:
        """Spectral coefficients (column 0 is the theta-mean)."""
        return np.asarray(np.fft.rfft(phys, axis=1, norm="forward"))

    def _radial(self, hat: CArray, f_in: CArray, f_out: CArray) -> CArray:
        return np.asarray(
            self.d_mat @ hat
            + self.d_in[:, None] * f_in[None, :]
            + self.d_out[:, None] * f_out[None, :]
        )

    def _divergence(
        self, hat: CArray, f_in: CArray, f_out: CArray, u_r_face: Array, u_th: Array
    ) -> CArray:
        """Flux-form ``div(u f)``: conservative in ``r``, spectral in ``theta``.

        The radial flux uses third-order upwind-biased face values and vanishes on
        both walls, so the theta-mean of ``f`` integrates exactly to the wall flux.
        """
        ext = np.concatenate([f_in[None, :], hat, f_out[None, :]], axis=0)
        f_face = np.where(
            u_r_face > 0.0,
            self.to_phys(self.f_plus @ ext),
            self.to_phys(self.f_minus @ ext),
        )
        flux = np.zeros((self.config.n_r + 1, self.n_theta))
        flux[1:-1] = self.faces[1:-1, None] * u_r_face * f_face
        flux_hat = self.to_hat(flux)
        radial = (flux_hat[1:] - flux_hat[:-1]) / self.area[:, None]
        f_phys = self.to_phys(hat)
        azimuthal = 1j * self.m[None, :] * self.to_hat(u_th * f_phys) / self.r[:, None]
        return np.asarray(radial + azimuthal)

    def _with_mean_vorticity(
        self, state: AnnulusState, z_in: CArray, z_out: CArray
    ) -> tuple[CArray, CArray, CArray]:
        """``zeta`` with the vorticity of the mean wind in column 0 (and walls).

        The transport of the mean vorticity by ``u_r`` is the Rayleigh term of
        shear instability; it must be advected although ``ubar`` is evolved apart.
        """
        c = self.config
        zero = np.zeros(self.k, dtype=np.complex128)
        ru = np.zeros((c.n_r, self.k), dtype=np.complex128)
        ru[:, 0] = self.r * state.ubar
        ru_in, ru_out = zero.copy(), zero.copy()
        ru_in[0] = self.r_in * state.ubar[0]
        ru_out[0] = c.radius * state.ubar[-1]
        full = np.array(state.zeta)
        full[:, 0] = (self._radial(ru, ru_in, ru_out)[:, 0] / self.r).real
        w_in, w_out = np.array(z_in), np.array(z_out)
        w_in[0] = 2.0 * state.ubar[0] / self.r_in
        w_out[0] = 2.0 * state.ubar[-1] / c.radius
        return full, w_in, w_out

    def _u_r_face(self, psi: CArray) -> Array:
        """Radial wind at the interior faces from the streamfunction."""
        psi_face = self.f_linear @ np.concatenate(
            [np.zeros((1, self.k)), psi, np.zeros((1, self.k))], axis=0
        )
        return self.to_phys(1j * self.m[None, :] * psi_face / self.faces[1:-1, None])

    # --- fields ---------------------------------------------------------------
    def psi(self, zeta: CArray) -> CArray:
        """Streamfunction from ``zeta = -lap psi`` with psi = 0 on both walls."""
        rhs = -np.array(zeta, dtype=np.complex128)
        rhs[:, 0] = 0.0
        out = _solve(self._poisson, rhs)
        out[:, 0] = 0.0
        return out

    def velocities(self, state: AnnulusState) -> tuple[Array, Array, CArray]:
        """Physical ``(u_r, u_theta)`` and the streamfunction."""
        psi = self.psi(state.zeta)
        u_r_hat = 1j * self.m[None, :] * psi / self.r[:, None]
        zero = np.zeros(self.k, dtype=np.complex128)
        u_th_hat = -self._radial(psi, zero, zero)
        u_th_hat[:, 0] = state.ubar
        return self.to_phys(u_r_hat), self.to_phys(u_th_hat), psi

    def _zeta_walls(self, psi: CArray) -> tuple[CArray, CArray]:
        c = self.config
        r = self.r
        if c.wall == "noslip":
            # psi = psi_rr x^2 / 2 + psi_rrr x^3 / 6 near a rigid wall (psi = psi_r = 0)
            d1, d2 = r[0] - self.r_in, r[1] - self.r_in
            z_in = (
                -2.0 * (psi[0] * d2**3 - psi[1] * d1**3) / (d1**2 * d2**2 * (d2 - d1))
            )
            d1, d2 = c.radius - r[-1], c.radius - r[-2]
            z_out = (
                -2.0 * (psi[-1] * d2**3 - psi[-2] * d1**3) / (d1**2 * d2**2 * (d2 - d1))
            )
        else:
            zero = np.zeros(self.k, dtype=np.complex128)
            dpsi_in = _wall_derivative(self.r_in, r[0], r[1], zero, psi[0], psi[1])
            dpsi_out = _wall_derivative(c.radius, r[-1], r[-2], zero, psi[-1], psi[-2])
            z_in = -2.0 * dpsi_in / self.r_in
            z_out = -2.0 * dpsi_out / c.radius
        z_in[0] = 0.0
        z_out[0] = 0.0
        return z_in, z_out

    def _flux_hat(self, t: float) -> CArray:
        """Outer-wall temperature flux (K m/s) as theta modes."""
        c = self.config
        q = c.forcing.flux(self.theta, t) / c.rho_cp
        return np.asarray(np.fft.rfft(q, norm="forward"))

    def _t_walls(self, temp: CArray, t: float) -> tuple[CArray, CArray]:
        c = self.config
        kappa = c.kappa
        zero = np.zeros(self.k, dtype=np.complex128)
        if c.outer_t == "fixed":
            t_out = zero.copy()
            t_out[0] = c.t_outer
        else:
            t_out = temp[-1] + self._flux_hat(t) / kappa * (c.radius - self.r[-1])
        if c.inner_t == "fixed":
            t_in = zero.copy()
            t_in[0] = c.t_inner
        else:
            q_in = zero.copy()
            q_in[0] = c.flux_inner / c.rho_cp
            t_in = temp[0] - q_in / kappa * (self.r[0] - self.r_in)
        return t_in, t_out

    # --- explicit terms ------------------------------------------------------
    def _u_ref(self, u_th: Array) -> Array:
        j, w = self._ref_index, self._ref_weight
        return np.asarray((1.0 - w) * u_th[j] + w * u_th[j + 1])

    def _drag_force(self, u_th: Array) -> Array:
        c = self.config
        u_ref = self._u_ref(u_th)
        stress = c.floor_drag * np.abs(u_ref) * u_ref / c.drag_height
        return np.asarray(-self.drag_weight[:, None] * stress[None, :])

    def explicit_terms(
        self, state: AnnulusState
    ) -> tuple[tuple[CArray, CArray, Array], Terms, tuple[CArray, CArray]]:
        """Explicit terms at ``state``, the mean-wind pieces and the wall vorticity."""
        c = self.config
        u_r, u_th, psi = self.velocities(state)
        z_in, z_out = self._zeta_walls(psi)
        t_in, t_out = self._t_walls(state.temp, state.t)
        i_m = 1j * self.m[None, :]
        adv_t_hat = np.zeros_like(state.temp)
        adv_z_hat = np.zeros_like(state.zeta)
        if c.advection:
            u_face = self._u_r_face(psi)
            adv_t_hat = self._divergence(state.temp, t_in, t_out, u_face, u_th)
            zeta_full, zf_in, zf_out = self._with_mean_vorticity(state, z_in, z_out)
            adv_z_hat = self._divergence(zeta_full, zf_in, zf_out, u_face, u_th)
        r_col = self.r[:, None]
        n_zeta = -adv_z_hat + c.beta * c.omega**2 * i_m * state.temp
        n_temp = np.where(self.mask, -adv_t_hat, 0.0)
        if c.tau_relax is not None:
            mean_t = float(
                np.sum(self.area * state.temp[:, 0].real) / np.sum(self.area)
            )
            n_temp[:, 0] -= mean_t / c.tau_relax
        cb = np.zeros(c.n_r)
        if c.coriolis_buoyancy:
            d_t = self._diffusive_tendency(state, n_temp)
            n_zeta = n_zeta - 2.0 * c.omega * c.beta * d_t
            cb = self._mean_coriolis_force(d_t[:, 0].real)
        force = self._drag_force(u_th)
        if c.floor_drag > 0.0:
            force_hat = self.to_hat(force)
            zero = np.zeros(self.k, dtype=np.complex128)
            f_out = force_hat[-1]
            curl = self._radial(r_col * force_hat, zero, c.radius * f_out) / r_col
            n_zeta = n_zeta + curl
        n_zeta[:, 0] = 0.0
        n_zeta = np.where(self.mask, n_zeta, 0.0)

        flux_uv = np.mean(u_r * u_th, axis=1)
        phi = np.zeros(c.n_r + 1)
        phi[1:-1] = self.faces[1:-1] ** 2 * 0.5 * (flux_uv[1:] + flux_uv[:-1])
        reynolds = -np.diff(phi) / self.moment
        drag = np.mean(force, axis=1) if c.floor_drag > 0.0 else np.zeros(c.n_r)
        terms = Terms(np.asarray(reynolds), cb, np.asarray(drag))
        n_ubar = terms.reynolds + terms.coriolis_buoyancy + terms.drag
        return (n_zeta, n_temp, n_ubar), terms, (z_in, z_out)

    def _diffusive_tendency(self, state: AnnulusState, n_temp: CArray) -> CArray:
        """Full ``dT/dt``: the explicit terms plus diffusion and the wall flux."""
        t_in, t_out = self._t_fixed_values()
        lap = self.apply_laplacian(state.temp, self._t_bands(), t_in, t_out)
        total = self.config.kappa * lap + self._flux_bc(state.t) + n_temp
        return np.asarray(np.where(self.mask, total, 0.0))

    def temperature_tendency(self, state: AnnulusState) -> CArray:
        """The model's own ``dT/dt`` (spectral) at ``state``."""
        (_, n_t, _), _, _ = self.explicit_terms(state)
        return self._diffusive_tendency(state, n_t)

    def _mean_coriolis_force(self, d_t0: Array) -> Array:
        """Mean azimuthal Coriolis-buoyancy force from the mean temperature tendency.

        ``(1/r) d(r F)/dr = -2 Omega beta dT'/dt`` with ``T'`` the anomaly about the
        area mean (the closed fluid's mass is fixed), integrated from the inner wall.
        """
        c = self.config
        g = d_t0 - float(np.sum(self.area * d_t0) / np.sum(self.area))
        cum = np.concatenate([[0.0], np.cumsum(self.area * g)])
        return np.asarray(-c.omega * c.beta * (cum[:-1] + cum[1:]) / self.r)

    @property
    def inertia_weights(self) -> Array:
        """Weights ``P`` with ``sum(P dT) = -dI'/dt`` for the moment of inertia."""
        c_i = self.moment / self.r
        tail = np.concatenate([np.cumsum(c_i[::-1])[::-1][1:], [0.0]])
        return np.asarray(self.area * (c_i + 2.0 * tail))

    def absolute_angular_momentum(self, state: AnnulusState) -> float:
        """Relative angular momentum plus the moment-of-inertia correction.

        ``2 pi [sum(M ubar) + Omega beta sum(P (T0 - <T0>))]``: with the consistent
        Coriolis-buoyancy force this changes only through wall and end-cap torques.
        """
        c = self.config
        t0 = state.temp[:, 0].real
        t_prime = t0 - float(np.sum(self.area * t0) / np.sum(self.area))
        return float(
            2.0
            * np.pi
            * (
                np.sum(self.moment * state.ubar)
                + c.omega * c.beta * np.sum(self.inertia_weights * t_prime)
            )
        )

    # --- implicit pieces ------------------------------------------------------
    def _factors(self, dt: float, omega: float | None) -> _Factors:
        """Implicit-operator factors for step ``dt`` and step ratio ``omega``.

        ``omega = dt / previous dt`` (``None`` for the first, backward Euler, step);
        the variable-step BDF2 weight on the new level is ``(1 + 2 w)/((1 + w) dt)``.
        """
        key = (round(dt, 9), 0.0 if omega is None else round(omega, 9))
        if key in self._cache:
            return self._cache[key]
        c = self.config
        a0 = 1.0 / dt if omega is None else (1.0 + 2.0 * omega) / ((1.0 + omega) * dt)
        m2r2 = self.m[None, :] ** 2 / self.r[:, None] ** 2

        def build(
            bands: tuple[Array, Array, Array], coeff: float, sink: float
        ) -> tuple[CArray, CArray, CArray]:
            lo, di, up = bands
            return _factor(
                np.repeat(-coeff * lo[:, None], self.k, axis=1).astype(np.complex128),
                (a0 + sink - coeff * (di[:, None] - m2r2)).astype(np.complex128),
                np.repeat(-coeff * up[:, None], self.k, axis=1).astype(np.complex128),
            )

        fz = build(
            self.bands(dirichlet_in=True, dirichlet_out=True), c.nu, c.endcap_rate
        )
        ft = build(self._t_bands(), c.kappa, 0.0)
        # mean wind: viscous flux operator in stress form, free slip at both walls
        k_face = np.zeros(c.n_r + 1)
        k_face[1:-1] = self.faces[1:-1] ** 3 * c.nu / np.diff(self.r)
        lo_u = k_face[:-1] / np.roll(self.r, 1)
        up_u = k_face[1:] / np.roll(self.r, -1)
        di_u = -(k_face[:-1] + k_face[1:]) / self.r
        lo_u[0] = 0.0
        up_u[-1] = 0.0
        mass = self.moment
        fu = _factor(
            (-lo_u)[:, None].astype(np.complex128),
            (a0 * mass - di_u + c.endcap_rate * mass)[:, None].astype(np.complex128),
            (-up_u)[:, None].astype(np.complex128),
        )
        out = _Factors(fz, ft, fu)
        self._cache[key] = out
        return out

    def _bc_zeta(self, z_in: CArray, z_out: CArray) -> CArray:
        bc = np.zeros((self.config.n_r, self.k), dtype=np.complex128)
        bc[0] += self.config.nu * self.w_in * z_in
        bc[-1] += self.config.nu * self.w_out * z_out
        return bc

    def _flux_bc(self, t: float) -> CArray:
        """Heat-flux (Neumann) wall terms of ``dT/dt``; none on a fixed-T wall."""
        c = self.config
        bc = np.zeros((c.n_r, self.k), dtype=np.complex128)
        if c.outer_t == "flux":
            bc[-1] += c.radius * self._flux_hat(t) / self.area[-1]
        if c.inner_t == "flux":
            q_in = np.zeros(self.k, dtype=np.complex128)
            q_in[0] = c.flux_inner / c.rho_cp
            bc[0] -= self.r_in * q_in / self.area[0]
        return bc

    def _bc_temp(self, t: float) -> CArray:
        """Right-hand side of the implicit temperature solve from the walls."""
        c = self.config
        t_in, t_out = self._t_fixed_values()
        bc = self._flux_bc(t)
        if c.outer_t == "fixed":
            bc[-1] += c.kappa * self.w_out * t_out
        if c.inner_t == "fixed":
            bc[0] += c.kappa * self.w_in * t_in
        return bc

    # --- stepping -------------------------------------------------------------
    def initial_state(self, amplitude: float = 1.0e-3, seed: int = 0) -> AnnulusState:
        """Quiescent air with a small random temperature perturbation."""
        rng = np.random.default_rng(seed)
        noise = rng.standard_normal((self.config.n_r, self.config.n_theta))
        temp = self.to_hat(amplitude * noise)
        temp = np.where(self.mask, temp, 0.0)
        temp[:, 0] = 0.0
        zeta = np.zeros_like(temp)
        return AnnulusState(0.0, zeta, temp, np.zeros(self.config.n_r))

    def step(self, state: AnnulusState, dt: float) -> None:
        """Advance ``state`` by ``dt`` seconds (variable-step SBDF2).

        The first step is backward Euler. ``dt`` may change between steps; keep
        increases to a factor of about 1.25 per step (see :meth:`suggest_dt`).
        """
        positive("dt", dt)
        first = state.previous is None or state.last_dt is None
        omega = None if first else dt / float(state.last_dt)  # type: ignore[arg-type]
        if omega is not None and omega >= _STEP_RATIO_LIMIT:
            msg = (
                f"step ratio {omega:.2f} must stay below 1 + sqrt(2) "
                "(variable-step BDF2 zero-stability)"
            )
            raise ValueError(msg)
        factors = self._factors(dt, omega)
        (n_z, n_t, n_u), _, (z_in, z_out) = self.explicit_terms(state)
        mass = self.moment
        if omega is None or state.previous is None or state.last_terms is None:
            rhs_z = state.zeta / dt + n_z
            rhs_t = state.temp / dt + n_t
            rhs_u = mass * (state.ubar / dt + n_u)
        else:
            p_z, p_t, p_u = state.previous
            l_z, l_t, l_u = state.last_terms
            w1, w2 = 1.0 + omega, omega**2 / (1.0 + omega)
            rhs_z = (w1 * state.zeta - w2 * p_z) / dt + w1 * n_z - omega * l_z
            rhs_t = (w1 * state.temp - w2 * p_t) / dt + w1 * n_t - omega * l_t
            rhs_u = mass * ((w1 * state.ubar - w2 * p_u) / dt + w1 * n_u - omega * l_u)
        zeta = _solve(factors.zeta, rhs_z + self._bc_zeta(z_in, z_out))
        zeta[:, 0] = 0.0
        temp = _solve(factors.temp, rhs_t + self._bc_temp(state.t + dt))
        ubar = _solve(factors.mean, rhs_u[:, None].astype(np.complex128))[:, 0].real
        zeta = np.asarray(np.where(self.mask, zeta, 0.0) * self.damping)
        temp = np.asarray(np.where(self.mask, temp, 0.0) * self.damping)
        state.previous = (state.zeta, state.temp, state.ubar)
        state.last_terms = (n_z, n_t, n_u)
        state.zeta, state.temp, state.ubar = zeta, temp, ubar
        state.t += dt
        state.last_dt = dt
        state.steps += 1

    def suggest_dt(
        self,
        state: AnnulusState,
        cfl: float = DEFAULT_CFL,
        dt_max: float = 10.0,
        ratio: float = 1.25,
    ) -> float:
        """A CFL time step on a geometric ladder, rising at most one rung per step.

        Keeping ``dt`` on the ladder ``dt_max / ratio**k`` lets the implicit
        operator factors be cached.
        """
        target = self.stable_dt(state, cfl, dt_max)
        k = max(0, math.ceil(math.log(dt_max / target) / math.log(ratio)))
        dt = dt_max / ratio**k
        if state.last_dt is not None and dt > state.last_dt * ratio * (1.0 + 1e-9):
            dt = state.last_dt * ratio
        return dt

    def advance(
        self,
        state: AnnulusState,
        t_end: float,
        *,
        cfl: float = DEFAULT_CFL,
        dt_max: float = 10.0,
        every: float | None = None,
        callback: Callable[[AnnulusState], None] | None = None,
    ) -> None:
        """Integrate to ``t_end`` on the CFL ladder; optional periodic callback."""
        next_call = state.t + every if every else math.inf
        while state.t < t_end - 1e-9:
            dt = min(self.suggest_dt(state, cfl, dt_max), t_end - state.t)
            self.step(state, dt)
            if callback is not None and state.t >= next_call - 1e-9:
                callback(state)
                next_call += every or 0.0

    def stable_dt(
        self, state: AnnulusState, cfl: float = DEFAULT_CFL, dt_max: float = 10.0
    ) -> float:
        """CFL time step from the current winds (radial and azimuthal)."""
        u_r, u_th, _ = self.velocities(state)
        c = self.config
        dr = np.minimum(self.width, np.r_[np.diff(self.r), self.width[-1]])
        r_eff = np.maximum(self.r, 1.0) * 2.0 * np.pi / c.n_theta
        if c.polar_filter:
            r_eff = np.full_like(r_eff, c.radius * 2.0 * np.pi / c.n_theta)
        rate = np.max(np.abs(u_r) / dr[:, None]) + np.max(np.abs(u_th) / r_eff[:, None])
        return float(min(dt_max, cfl / max(rate, 1e-12)))

    # --- linear stability -----------------------------------------------------
    def _mode_tendency(self, state: AnnulusState, m: int) -> tuple[CArray, CArray]:
        """Semi-discrete ``d zeta_m/dt`` and ``d T_m/dt`` (no time stepping)."""
        c = self.config
        (n_z, n_t, _), _, (z_in, z_out) = self.explicit_terms(state)
        t_in, t_out = self._t_fixed_values()
        lap_z = self.apply_laplacian(
            state.zeta, self.bands(dirichlet_in=True, dirichlet_out=True), z_in, z_out
        )
        lap_t = self.apply_laplacian(state.temp, self._t_bands(), t_in, t_out)
        return (
            c.nu * lap_z[:, m] + n_z[:, m] - c.endcap_rate * state.zeta[:, m],
            c.kappa * lap_t[:, m] + n_t[:, m],
        )

    def _t_fixed_values(self) -> tuple[CArray, CArray]:
        c = self.config
        zero = np.zeros(self.k, dtype=np.complex128)
        t_in, t_out = zero.copy(), zero.copy()
        if c.inner_t == "fixed":
            t_in[0] = c.t_inner
        if c.outer_t == "fixed":
            t_out[0] = c.t_outer
        return t_in, t_out

    def max_growth_rate(self, base: AnnulusState, m: int) -> float:
        """Largest real eigenvalue of the semi-discrete linearisation about ``base``.

        Probes the model's own operators with small perturbations of azimuthal mode
        ``m`` (vorticity and temperature), so it tests the spatial discretisation
        without time stepping. ``base`` supplies the mean temperature and wind.
        """
        n = self.config.n_r
        eps = 1.0e-6
        f0 = np.concatenate(self._mode_tendency(base, m))
        jac = np.zeros((2 * n, 2 * n), dtype=np.complex128)
        for j in range(2 * n):
            probe = AnnulusState(
                base.t,
                np.zeros_like(base.zeta),
                np.array(base.temp),
                np.array(base.ubar),
            )
            if j < n:
                probe.zeta[j, m] = eps
            else:
                probe.temp[j - n, m] += eps
            jac[:, j] = (np.concatenate(self._mode_tendency(probe, m)) - f0) / eps
        return float(np.max(np.linalg.eigvals(jac).real))

    # --- diagnostics ----------------------------------------------------------
    def energy(self, state: AnnulusState) -> float:
        """Kinetic energy ``int (u_r^2 + u_theta^2) / 2 dA`` per unit length."""
        u_r, u_th, _ = self.velocities(state)
        per_r = 0.5 * np.mean(u_r**2 + u_th**2, axis=1)
        return float(2.0 * np.pi * np.sum(self.area * per_r))

    def mean_wind_energy(self, state: AnnulusState) -> float:
        """Energy of the theta-mean flow alone."""
        return float(np.pi * np.sum(self.area * state.ubar**2))

    def enstrophy(self, state: AnnulusState) -> float:
        """``int zeta^2 / 2 dA`` including the mean vorticity of ``ubar``."""
        zero = np.zeros(self.k, dtype=np.complex128)
        ru = np.zeros((self.config.n_r, self.k), dtype=np.complex128)
        ru[:, 0] = self.r * state.ubar
        wall_out = zero.copy()
        wall_out[0] = self.config.radius * state.ubar[-1]
        wall_in = zero.copy()
        wall_in[0] = self.r_in * state.ubar[0]
        z_hat = np.array(state.zeta)
        z_hat[:, 0] = (self._radial(ru, wall_in, wall_out)[:, 0] / self.r).real
        zeta = self.to_phys(z_hat)
        return float(np.pi * np.sum(self.area * np.mean(zeta**2, axis=1)))

    def angular_momentum(self, state: AnnulusState) -> float:
        """Relative angular momentum ``2 pi int r^2 ubar dr`` (per unit density)."""
        return float(2.0 * np.pi * np.sum(self.moment * state.ubar))

    def heat_content(self, state: AnnulusState) -> float:
        """Integral of the temperature anomaly ``int T dA``."""
        return float(2.0 * np.pi * np.sum(self.area * state.temp[:, 0].real))

    def mean_wind_terms(self, state: AnnulusState) -> Terms:
        """The mean-wind forcing pieces at ``state`` (for budget checks)."""
        return self.explicit_terms(state)[1]

    def wind_at_height(self, state: AnnulusState, height: float) -> Array:
        """Azimuthal wind profile ``u_theta(theta)`` at ``height`` above the floor."""
        _, u_th, _ = self.velocities(state)
        r0 = self.config.radius - positive("height", height)
        j = int(np.clip(np.searchsorted(self.r, r0) - 1, 0, self.config.n_r - 2))
        w = float(np.clip((r0 - self.r[j]) / (self.r[j + 1] - self.r[j]), 0.0, 1.0))
        return np.asarray((1.0 - w) * u_th[j] + w * u_th[j + 1])

    def spectrum(self, state: AnnulusState, r_min: float, r_max: float) -> Array:
        """Kinetic energy by azimuthal wavenumber for ``r_min <= r <= r_max``."""
        psi = self.psi(state.zeta)
        zero = np.zeros(self.k, dtype=np.complex128)
        u_r = 1j * self.m[None, :] * psi / self.r[:, None]
        u_th = -self._radial(psi, zero, zero)
        sel = (self.r >= r_min) & (self.r <= r_max)
        e = 0.5 * (np.abs(u_r[sel]) ** 2 + np.abs(u_th[sel]) ** 2)
        weight = self.area[sel][:, None] * 2.0 * np.pi
        spec = 2.0 * np.sum(weight * e, axis=0)
        spec[0] = float(np.pi * np.sum(self.area[sel] * state.ubar[sel] ** 2))
        return np.asarray(spec)
