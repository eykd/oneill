"""The L3 production cases: configuration, a run, its summary and its figures.

Everything here is pure (no files, no clock): the adapter times the run, writes the
JSON summary and saves the figures. ``summarize`` reduces a run to the numbers in the
``gen:annulus`` doc block, taken from the second half of the record.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import numpy as np

from oneill.atmo.annulus2d import Annulus, AnnulusConfig, Forcing
from oneill.atmo.annulus_diagnostics import (
    Sample,
    large_scale_fraction,
    lock_fraction,
    spectrum_scale,
    stationarity,
    take_sample,
    wind_statistics,
)
from oneill.atmo.gas import SP413
from oneill.core.cylinder import Cylinder
from oneill.core.figspec import FigureSpec, Mesh, Panel, Series

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    type Array = NDArray[np.float64]

DAY = 86400.0
PRODUCTION_N_R = 96
PRODUCTION_STRETCH = 2.4
"""96 cells clustered to a 2.4 m floor cell: the 8.9 m grid (64 cells) put only about
two cells across the 17-20 m heated layer; the winds changed by 3-4 % between them."""
FLOOR_TEMPERATURE = 290.0
SAMPLE_EVERY = 1800.0
ENERGY_POINTS = 60
PROFILE_POINTS = 24
ARREST_LOW, ARREST_HIGH = 0.08, 0.18
DAY_FLUX, NIGHT_FLUX, STEADY_FLUX = 200.0, -30.0, 100.0


@dataclass(frozen=True, slots=True)
class Case:
    """One production run.

    Attributes:
        name: Short identifier (CLI ``--l3-case``).
        description: One line for the doc.
        nu: Eddy viscosity in m^2/s.
        coriolis_buoyancy: Include the mass-consistent Coriolis-buoyancy force.
        r_in_fraction: Inner wall radius over ``R``.
        heating: ``"cycle"`` (day/night, flux-neutral windows), ``"steady"`` or
            ``"uniform_sink"`` (day/night, windows inert, relax the area mean).
    """

    name: str
    description: str
    nu: float = 1.0
    coriolis_buoyancy: bool = True
    r_in_fraction: float = 0.1
    heating: str = "cycle"


CASES: dict[str, Case] = {
    c.name: c
    for c in (
        Case(
            "primary", "day/night, flux-neutral windows, Coriolis-buoyancy on, nu_t = 1"
        ),
        Case(
            "steady",
            "steady 100 W/m^2 land flux, otherwise as primary",
            heating="steady",
        ),
        Case(
            "cb_off", "as primary with Coriolis-buoyancy off", coriolis_buoyancy=False
        ),
        Case("nu10", "as primary with nu_t = 10", nu=10.0),
        Case("rin005", "as primary with the inner wall at 0.05 R", r_in_fraction=0.05),
        Case("rin02", "as primary with the inner wall at 0.2 R", r_in_fraction=0.2),
        Case(
            "uniform_sink",
            "inert windows with a uniform sink on the area mean (tau_T = 1 day)",
            heating="uniform_sink",
        ),
    )
}


def case_config(
    case: Case,
    cyl: Cylinder,
    n_r: int = PRODUCTION_N_R,
    n_theta: int = 256,
    stretch: float = PRODUCTION_STRETCH,
    edge_cells: float = 0.0,
) -> AnnulusConfig:
    """The model configuration of ``case`` for habitat ``cyl``.

    ``edge_cells`` smooths the strip edges over that many azimuthal grid spacings
    (``0``: step edges).
    """
    width = edge_cells * 2.0 * math.pi / n_theta
    gas = SP413.at_temperature(FLOOR_TEMPERATURE)  # the half-atmosphere habitat
    if case.heating == "steady":
        forcing = Forcing(land_flux=STEADY_FLUX, edge_width=width)
    elif case.heating == "uniform_sink":
        forcing = Forcing(
            DAY_FLUX, NIGHT_FLUX, period=DAY, window="zero", edge_width=width
        )
    else:
        forcing = Forcing(DAY_FLUX, NIGHT_FLUX, period=DAY, edge_width=width)
    return AnnulusConfig.from_habitat(
        cyl,
        gas,
        nu=case.nu,
        r_in_fraction=case.r_in_fraction,
        coriolis_buoyancy=case.coriolis_buoyancy,
        forcing=forcing,
        tau_relax=DAY if case.heating == "uniform_sink" else None,
        n_r=n_r,
        n_theta=n_theta,
        stretch=stretch,
    )


@dataclass(frozen=True, eq=False)
class CaseResult:
    """A finished run: sample series and the final fields."""

    case: Case
    config: AnnulusConfig
    samples: tuple[Sample, ...]
    steps: int
    radius: Array
    theta: Array
    zeta: Array
    temp: Array


def run_case(  # noqa: PLR0913
    case: Case,
    days: float,
    *,
    cyl: Cylinder | None = None,
    n_r: int = PRODUCTION_N_R,
    n_theta: int = 256,
    stretch: float = PRODUCTION_STRETCH,
    edge_cells: float = 0.0,
    sample_every: float = SAMPLE_EVERY,
    seed: int = 0,
    dt_max: float = 10.0,
    progress: Callable[[float, float], None] | None = None,
) -> CaseResult:
    """Integrate ``case`` for ``days`` simulated days, sampling the diagnostics.

    ``progress(t, energy)`` is called at each sample.
    """
    habitat = Cylinder() if cyl is None else cyl
    config = case_config(case, habitat, n_r, n_theta, stretch, edge_cells)
    model = Annulus(config)
    state = model.initial_state(seed=seed)
    samples: list[Sample] = []

    def record(s: Any) -> None:  # noqa: ANN401
        sample = take_sample(model, s)
        samples.append(sample)
        if progress is not None:
            progress(sample.t, sample.energy)

    record(state)
    model.advance(state, days * DAY, dt_max=dt_max, every=sample_every, callback=record)
    if samples[-1].t < state.t - 1e-6:
        record(state)
    return CaseResult(
        case=case,
        config=config,
        samples=tuple(samples),
        steps=state.steps,
        radius=np.array(model.r),
        theta=np.array(model.theta),
        zeta=model.to_phys(state.zeta),
        temp=model.to_phys(state.temp),
    )


def _first_cell(cfg: AnnulusConfig) -> float:
    """Width in m of the cell next to the floor."""
    gap = cfg.radius * (1.0 - cfg.r_in_fraction)
    if cfg.stretch == 0.0:
        return gap / cfg.n_r
    return gap * (
        1.0 - math.tanh(cfg.stretch * (1.0 - 1.0 / cfg.n_r)) / math.tanh(cfg.stretch)
    )


def _round(x: float) -> float:
    return float(f"{x:.6g}")


def _downsample(t: Array, y: Array, n: int) -> list[list[float]]:
    idx = np.unique(np.linspace(0, len(t) - 1, min(n, len(t))).astype(int))
    return [[_round(float(t[i] / DAY)), _round(float(y[i]))] for i in idx]


def summarize(result: CaseResult) -> dict[str, Any]:
    """JSON-ready summary of a run (statistics over the second half of the record)."""
    cfg = result.config
    samples = result.samples
    t = np.array([s.t for s in samples])
    energy = np.array([s.energy for s in samples])
    tail = t >= 0.5 * t[-1]
    tail_samples = [s for s, keep in zip(samples, tail, strict=True) if keep]
    spectrum = np.mean([s.spectrum for s in tail_samples], axis=0)
    r_ref = 0.75 * cfg.radius
    scale = spectrum_scale(spectrum, r_ref)
    winds = [wind_statistics(s.wind) for s in tail_samples]
    rms10 = float(np.mean([w[0] for w in winds]))
    gust10 = float(np.mean([w[1] for w in winds]))
    zonal10 = float(np.mean([np.mean(s.wind) for s in tail_samples]))
    eps = float(np.mean([s.dissipation for s in tail_samples]))
    mu = cfg.endcap_rate + cfg.floor_drag * rms10 / cfg.radius
    arrest = math.sqrt(eps / mu**3) if eps > 0.0 and mu > 0.0 else 0.0
    ubar = np.mean([s.ubar for s in tail_samples], axis=0)
    pick = np.unique(np.linspace(0, len(result.radius) - 1, PROFILE_POINTS).astype(int))
    stat = stationarity(t, energy)
    mean_energy = float(np.mean([s.energy for s in tail_samples]))
    return {
        "name": result.case.name,
        "description": result.case.description,
        "days": _round(float(t[-1] / DAY)),
        "steps": result.steps,
        "n_r": cfg.n_r,
        "first_cell_m": _round(_first_cell(cfg)),
        "n_theta": cfg.n_theta,
        "nu": cfg.nu,
        "r_in_fraction": cfg.r_in_fraction,
        "coriolis_buoyancy": cfg.coriolis_buoyancy,
        "stationary": stat.stationary,
        "energy_drift": _round(stat.drift),
        "energy": _round(mean_energy),
        "mean_wind_energy_fraction": _round(
            float(np.mean([s.mean_energy for s in tail_samples])) / mean_energy
            if mean_energy > 0.0
            else 0.0
        ),
        "m_peak": scale.m_peak,
        "length_peak_m": _round(scale.length_peak),
        "length_centroid_m": _round(scale.length_centroid),
        "lock_fraction": _round(lock_fraction(spectrum, cfg.forcing.strips)),
        "large_scale_fraction": _round(large_scale_fraction(spectrum)),
        "rms_wind_10m": _round(rms10),
        "gust_10m": _round(gust10),
        "zonal_wind_10m": _round(zonal10),
        "nusselt": _round(float(np.mean([s.nusselt for s in tail_samples]))),
        "dissipation": _round(eps),
        "mu": _round(mu),
        "arrest_low_m": _round(ARREST_LOW * arrest),
        "arrest_high_m": _round(ARREST_HIGH * arrest),
        "mean_wind_profile": [
            [_round(float(result.radius[i])), _round(float(ubar[i]))] for i in pick
        ],
        "energy_series": _downsample(t, energy, ENERGY_POINTS),
    }


# --- figures ----------------------------------------------------------------------


def _polar(result: CaseResult, field: Array) -> Mesh:
    theta = np.append(result.theta, 2.0 * math.pi)
    wrapped = np.concatenate([field, field[:, :1]], axis=1)
    r = result.radius[:, None]
    return Mesh(
        "",
        r * np.cos(theta)[None, :],
        r * np.sin(theta)[None, :],
        wrapped,
        cmap="RdBu_r",
        symmetric=True,
    )


def case_figures(result: CaseResult) -> tuple[FigureSpec, ...]:
    """Vorticity and temperature snapshot, spectrum against time, mean-wind profile."""
    name = result.case.name
    zeta = _polar(result, result.zeta)
    temp = _polar(result, result.temp)
    snapshot = FigureSpec(
        f"l3_{name}_snapshot.png",
        (
            Panel(
                "x [m]",
                "y [m]",
                (),
                title="vorticity [1/s]",
                equal_aspect=True,
                meshes=(
                    Mesh(
                        "zeta [1/s]",
                        zeta.x,
                        zeta.y,
                        zeta.values,
                        "RdBu_r",
                        symmetric=True,
                    ),
                ),
            ),
            Panel(
                "x [m]",
                "y [m]",
                (),
                title="temperature anomaly [K]",
                equal_aspect=True,
                meshes=(
                    Mesh(
                        "T [K]", temp.x, temp.y, temp.values, "RdBu_r", symmetric=True
                    ),
                ),
            ),
        ),
        (1, 2),
        (12.0, 6.0),
        suptitle=f"{name}: end of run (floor at the rim)",
    )
    t_days = np.array([s.t for s in result.samples]) / DAY
    spec = np.array([s.spectrum[1:] for s in result.samples]).T
    m = np.arange(1, spec.shape[0] + 1, dtype=np.float64)
    top = min(spec.shape[0], 40)
    grid_t, grid_m = np.meshgrid(t_days, m[:top])
    log_spec = np.log10(np.maximum(spec[:top], 1e-30))
    log_spec = np.maximum(log_spec, log_spec.max() - 8.0)  # show eight decades
    spectrum = FigureSpec(
        f"l3_{name}_spectrum.png",
        (
            Panel(
                "time [days]",
                "azimuthal wavenumber m",
                (),
                title="log10 kinetic energy by wavenumber",
                meshes=(Mesh("log10 E", grid_t, grid_m, log_spec, "magma"),),
            ),
        ),
        (1, 1),
        (8.0, 5.0),
        suptitle=f"{name}: azimuthal spectrum",
    )
    n = len(result.samples)
    tail = result.samples[n // 2 :]
    ubar = np.mean([s.ubar for s in tail], axis=0)
    last = result.samples[-1].ubar
    meanwind = FigureSpec(
        f"l3_{name}_meanwind.png",
        (
            Panel(
                "mean azimuthal wind [m/s]",
                "radius [m]",
                (
                    Series("second-half mean", ubar, result.radius),
                    Series("final", last, result.radius, "--"),
                ),
                title="theta-mean wind (positive is spinward)",
                legend_size=8.0,
            ),
        ),
        (1, 1),
        (6.0, 6.0),
        suptitle=f"{name}: zonal wind",
    )
    return snapshot, spectrum, meanwind
