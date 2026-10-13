"""The atmosphere report: every number in the write-up, computed from the models.

:func:`build_report` runs the gas, L0, L2b, L2 and L1 models for a habitat and both
gas presets and returns frozen result objects. Nothing here prints, formats or
draws; the adapters turn these into markdown and figures. Every row carries the
label of the assumption it depends on (gas preset and temperature, floor
condition), and :attr:`AtmosphereReport.caveats` holds the standing caveats of
each section.
"""

import math
from dataclasses import dataclass
from functools import cache

import numpy as np

from oneill.atmo.equilibrium import (
    CloudBase,
    adiabatic_axis_ratio,
    adiabatic_temperature_drop,
    axis_geopotential,
    dry_lapse_rate,
    isothermal_axis_ratio,
    lifting_condensation_level,
    total_mass,
)
from oneill.atmo.gas import ONEILL_1974, SP413, Gas
from oneill.atmo.numbers import DimensionlessTable, dimensionless_table
from oneill.atmo.parcel import (
    Release,
    RunSettings,
    simulate_parcel,
    steady_rise_speed,
)
from oneill.atmo.spinup import (
    DEFAULT_FLOOR_DRAG,
    DRAG_HEIGHT,
    LINEAR_SINK_EPSILON,
    SpinupParams,
    hull_spin_up,
    peak_wind_at_height,
    steady_lag,
    sudden_change,
    time_to_adjust,
    time_to_corotation,
)
from oneill.atmo.surface import (
    CANOPY_NOTE,
    FOREST,
    ORCHARD,
    SUBURBAN,
    Stability,
    Terrain,
    dynamic_pressure,
    friction_velocity,
    gust_speed,
    rotational_stability_profile,
    spinward_fraction,
    unstable_layer,
)
from oneill.core.cylinder import Cylinder

EARTH_OMEGA = 7.2921159e-5
"""Earth's sidereal spin rate in rad/s."""

FLOOR_TEMPERATURE = 290.0
"""Floor temperature in K used for every number in the report."""

PRESETS = (ONEILL_1974, SP413)
TERRAINS = (SUBURBAN, ORCHARD, FOREST)
HUMIDITIES = (0.5, 0.7)
WIND_SPEED = 3.0
WIND_HEIGHT = 20.0
STABILITY_HEIGHTS = (1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 30.0)
HEADINGS = (("spinward", 90.0), ("antispinward", 270.0), ("axial", 0.0))
BUOYANCY_EXCESSES = (1.0, 3.0)
PARCEL_DRAGS = (0.0, 0.01, 0.05)
RAMPS = (3600.0, 86400.0)
EDDY_VISCOSITIES = (1.0, 10.0)
FLOOR_LABEL = f"canopy drag, C_D = {DEFAULT_FLOOR_DRAG:.4f} (suburb, z0 = 0.5 m)"
BLEND_WIDTHS = (0.15, 0.3, 0.6)
BLEND_DT_MAX = 150.0
"""Time step for the blend rows (50 % times drift ~3 % at the default 600 s)."""
"""Closure blend widths tried in the sensitivity rows."""
STEP = 0.1
"""Fractional hull spin change used in the sudden-change rows."""
LAG_FRACTION_PER_DAY = 0.01
SNAPSHOT_HOURS = (1.0, 3.0, 6.0, 12.0, 24.0)
PARCEL_SECONDS = 300.0
DEFLECTION_RISE = 2.0

CAVEATS = {
    "annulus": (
        "Two-dimensional (axis-invariant) Boussinesq convection with constant eddy "
        "viscosity, a free-slip inner wall at the stated radius (an artifact: the "
        "habitat has no inner wall) and a canopy-drag floor. It can support or "
        "undermine the 2-D hypothesis but cannot show where axis-invariance breaks "
        "down. The arrest scale uses the run's measured dissipation. Energy is "
        "conserved only to a few percent when eddies saturate."
    ),
    "gas": (
        "Dry mixture ideal gas; water vapour is a passive trace carried in the "
        "SP-413 total pressure but not in R_s, c_p or the density (about 2 % bias). "
        "Viscosity is Wilke mixing of per-species Sutherland laws."
    ),
    "equilibrium": (
        "Still, horizontally uniform air in hydrostatic balance; the isothermal and "
        "dry-adiabatic profiles bracket reality. Cloud base uses Bolton's formula "
        "for air with a dry-adiabatic ascent; true height converts the geopotential, "
        "Earth-equivalent height is phi / g."
    ),
    "numbers": (
        "Viscous scales use the half-length L/2 = 10 km. Eddy viscosities are "
        "assumed constants, not derived. N is for an isothermal atmosphere at the "
        "floor."
    ),
    "surface": (
        f"{CANOPY_NOTE.capitalize()}. The canopy attenuation follows from shear "
        "continuity and is low for the suburb (0.87, against 1-3 in urban fits). The "
        "20 m reference height lies inside the forest's roughness sublayer, so the "
        "forest's u* and gusts are overestimates."
    ),
    "parcel": (
        "A point-mass parcel with no perturbation pressure or added mass: a ceiling "
        "for compact thermals, not a prediction. Real thermals are elongated along "
        "the axis by inertial waves and pressure absorbs the Coriolis force."
    ),
    "spinup": (
        "One-dimensional axisymmetric model with constant eddy viscosities and a "
        "canopy-drag floor (C_D from the suburban log law at 10 m); winds are quoted "
        "at that 10 m reference height, the one the drag law is written for. The "
        "end-cap sink is the linear Greenspan-Howard result, valid for "
        "epsilon = |u| / (Omega r) <= 0.3; rows with a larger peak epsilon (fast "
        "ramps, the 50 % step) are order-of-magnitude only, and turbulent end caps "
        "are not linear. Adjustment times depend on the blend width of the "
        "viscosity closure; the sensitivity band is shown, and 50 % times drift "
        "about 3 % with the production time step."
    ),
}


@dataclass(frozen=True, slots=True)
class GasRow:
    """Gas properties for one preset at the floor temperature."""

    label: str
    temperature: float
    p_dry: float
    p_vapour: float
    molar_mass: float
    r_specific: float
    cp: float
    gamma: float
    viscosity: float
    nu: float
    density: float


@dataclass(frozen=True, slots=True)
class EquilibriumRow:
    """L0 equilibrium numbers for one preset."""

    label: str
    phi_axis: float
    axis_ratio_isothermal: float
    axis_ratio_adiabatic: float
    dp_isothermal: float
    dp_adiabatic: float
    delta_t_adiabatic: float
    floor_lapse_rate: float
    mass_isothermal: float
    clouds: tuple[CloudBase, ...]
    humidities: tuple[float, ...]


@dataclass(frozen=True, slots=True)
class NumbersRow:
    """L0 dimensionless numbers for one preset."""

    label: str
    table: DimensionlessTable
    coriolis_vs_earth: float


@dataclass(frozen=True, slots=True)
class SurfaceRow:
    """L2b wind numbers for one terrain and preset."""

    terrain: Terrain
    gas_label: str
    ustar: float
    gust_speed: float
    gust_factor: float
    q_mean: float
    q_gust: float
    layer: tuple[float, float] | None


@dataclass(frozen=True, slots=True)
class StabilityRow:
    """Rotational stability class by height for one terrain and wind direction."""

    terrain: str
    heading: str
    heights: tuple[float, ...]
    classes: tuple[Stability, ...]


@dataclass(frozen=True, eq=False)
class ParcelRow:
    """L2 parcel numbers for one release.

    Attributes:
        temperature_excess: Parcel minus environment temperature in K.
        drag: Drag coefficient in 1/s.
        max_rise: Largest height gain within the short run, m.
        peak_speed: Largest speed relative to the air in the short run, m/s.
        deflection: Spinward displacement at ``DEFLECTION_RISE`` m of rise, m.
        drift_speed: Spinward velocity at that point, m/s.
        steady_rise: Long-time rise ``k b / (k^2 + 4 omega^2)`` in m/s (0 for k=0).
        steady_drift: Long-time spinward drift ``2 omega w / k`` in m/s.
        time_to_lcl: Dry-parcel time to the cloud base in s, if computed.
        t: Sample times of the short run, s.
        height: Heights of the short run, m.
        spinward: Spinward displacements of the short run, m.
    """

    temperature_excess: float
    drag: float
    max_rise: float
    peak_speed: float
    deflection: float | None
    drift_speed: float | None
    steady_rise: float
    steady_drift: float
    time_to_lcl: float | None
    t: np.ndarray
    height: np.ndarray
    spinward: np.ndarray


@dataclass(frozen=True, slots=True)
class SpinRow:
    """Hull spin-up from rest: one ramp, viscosity, end-cap and floor choice."""

    ramp: float
    nu: float
    endcap: bool
    floor: str
    time_to_90: float | None
    peak_wind_ref: float
    peak_wind_interior: float
    max_epsilon: float
    sink_valid: bool


@dataclass(frozen=True, slots=True)
class SuddenRow:
    """Sudden hull spin change with a Rayleigh-switched viscosity."""

    delta_fraction: float
    nu_stable: float
    nu_unstable: float
    endcap: bool
    time_to_adjust: float | None
    max_epsilon: float
    sink_valid: bool


@dataclass(frozen=True, slots=True)
class BlendRow:
    """Sensitivity of a sudden-change adjustment time to the closure blend width."""

    delta_fraction: float
    widths: tuple[float, ...]
    times: tuple[float | None, ...]
    max_epsilon: float
    sink_valid: bool


@dataclass(frozen=True, slots=True)
class LagRow:
    """Steady wind lag for a hull spin change of 1 % per day."""

    nu: float
    endcap: bool
    floor: str
    peak: float
    wind_ref: float


@dataclass(frozen=True, slots=True)
class MolecularRow:
    """Spin-up with molecular viscosity only (end caps on, canopy-drag floor).

    Illustrative only: the laminar end-cap sink is invalid at the Ekman-layer
    Reynolds number ``Re = U / sqrt(nu Omega)`` (in the thousands), where the
    layers are turbulent, and the 1 h ramp also exceeds ``epsilon = 0.3``.
    """

    label: str
    nu: float
    time_to_90: float | None
    peak_wind_ref: float
    ekman_reynolds: float


@dataclass(frozen=True, eq=False)
class Snapshot:
    """Wind profile ``u(r)`` at one time of the reference spin-up."""

    hours: float
    radius: np.ndarray
    u: np.ndarray


@dataclass(frozen=True, eq=False)
class AtmosphereReport:
    """Everything the write-up needs, for one habitat."""

    cylinder: Cylinder
    gas: tuple[GasRow, ...]
    equilibrium: tuple[EquilibriumRow, ...]
    numbers: tuple[NumbersRow, ...]
    surface: tuple[SurfaceRow, ...]
    stability: tuple[StabilityRow, ...]
    parcels: tuple[ParcelRow, ...]
    spin: tuple[SpinRow, ...]
    sudden: tuple[SuddenRow, ...]
    blend: tuple[BlendRow, ...]
    lag: tuple[LagRow, ...]
    molecular: tuple[MolecularRow, ...]
    snapshots: tuple[Snapshot, ...]
    corotation_hours: np.ndarray
    corotation: np.ndarray
    caveats: dict[str, str]


def _gas_rows(gases: tuple[Gas, ...]) -> tuple[GasRow, ...]:
    return tuple(
        GasRow(
            label=g.label,
            temperature=g.temperature,
            p_dry=g.p_rim,
            p_vapour=g.p_vapour,
            molar_mass=g.molar_mass,
            r_specific=g.r_specific,
            cp=g.cp,
            gamma=g.gamma,
            viscosity=g.viscosity(),
            nu=g.nu,
            density=g.rho_rim,
        )
        for g in gases
    )


def _equilibrium_rows(
    cyl: Cylinder, gases: tuple[Gas, ...]
) -> tuple[EquilibriumRow, ...]:
    rows = []
    for g in gases:
        iso, adi = isothermal_axis_ratio(cyl, g), adiabatic_axis_ratio(cyl, g)
        rows.append(
            EquilibriumRow(
                label=g.label,
                phi_axis=axis_geopotential(cyl),
                axis_ratio_isothermal=iso,
                axis_ratio_adiabatic=adi,
                dp_isothermal=g.p_rim * (1.0 - iso),
                dp_adiabatic=g.p_rim * (1.0 - adi),
                delta_t_adiabatic=adiabatic_temperature_drop(cyl, g),
                floor_lapse_rate=dry_lapse_rate(cyl, g, 0.0),
                mass_isothermal=total_mass(cyl, g.p_rim, iso),
                clouds=tuple(
                    lifting_condensation_level(cyl, g, rh) for rh in HUMIDITIES
                ),
                humidities=HUMIDITIES,
            )
        )
    return tuple(rows)


def _surface_rows(
    cyl: Cylinder, gases: tuple[Gas, ...]
) -> tuple[tuple[SurfaceRow, ...], tuple[StabilityRow, ...]]:
    surface: list[SurfaceRow] = []
    stability: list[StabilityRow] = []
    for terrain in TERRAINS:
        ustar = friction_velocity(WIND_SPEED, WIND_HEIGHT, terrain)
        gust = gust_speed(WIND_SPEED, ustar)
        layer = unstable_layer(ustar, terrain, cyl.omega)
        surface.extend(
            SurfaceRow(
                terrain=terrain,
                gas_label=g.label,
                ustar=ustar,
                gust_speed=gust,
                gust_factor=gust / WIND_SPEED,
                q_mean=dynamic_pressure(g.rho_rim, WIND_SPEED),
                q_gust=dynamic_pressure(g.rho_rim, gust),
                layer=layer,
            )
            for g in gases
        )
        for name, heading in HEADINGS:
            profile = rotational_stability_profile(
                STABILITY_HEIGHTS,
                ustar,
                terrain,
                cyl.omega,
                spinward_fraction(heading),
            )
            stability.append(
                StabilityRow(
                    terrain.name,
                    name,
                    STABILITY_HEIGHTS,
                    tuple(Stability(int(c)) for c in profile.classification),
                )
            )
    return tuple(surface), tuple(stability)


def _parcel_rows(cyl: Cylinder, gas: Gas) -> tuple[ParcelRow, ...]:
    drags = (*PARCEL_DRAGS, 2.0 * cyl.omega)
    lcl = lifting_condensation_level(cyl, gas, HUMIDITIES[0]).height
    short = RunSettings(
        t_end=PARCEL_SECONDS, sample_dt=1.0, max_step=1.0, rise=DEFLECTION_RISE
    )
    rows = []
    for excess in BUOYANCY_EXCESSES:
        buoyancy = cyl.surface_gravity * excess / gas.temperature
        for k in drags:
            release = Release(temperature_excess=excess, drag=k)
            run = simulate_parcel(cyl, gas, release, settings=short)
            t_lcl = None
            if k == drags[-1] and lcl is not None:
                long = RunSettings(
                    t_end=60000.0, sample_dt=50.0, max_step=5.0, stop_height=lcl
                )
                t_lcl = simulate_parcel(
                    cyl,
                    gas,
                    release,
                    settings=long,
                    relative_humidity=HUMIDITIES[0],
                ).time_to_lcl
            w = steady_rise_speed(buoyancy, k, cyl.omega)
            rows.append(
                ParcelRow(
                    temperature_excess=excess,
                    drag=k,
                    max_rise=run.max_rise,
                    peak_speed=run.peak_speed,
                    deflection=run.deflection,
                    drift_speed=run.drift_speed,
                    steady_rise=w,
                    steady_drift=0.0 if k == 0.0 else 2.0 * cyl.omega * w / k,
                    time_to_lcl=t_lcl,
                    t=run.t,
                    height=run.height,
                    spinward=run.spinward,
                )
            )
    return tuple(rows)


def _spin_rows(
    cyl: Cylinder,
) -> tuple[tuple[SpinRow, ...], tuple[Snapshot, ...], np.ndarray, np.ndarray]:
    rows = []
    snapshots: tuple[Snapshot, ...] = ()
    hours = corot = np.zeros(0)
    for ramp in RAMPS:
        for nu in EDDY_VISCOSITIES:
            for endcap in (True, False):
                params = SpinupParams(nu, endcap=endcap, n_cells=100)
                run = hull_spin_up(cyl, params, cyl.omega / ramp, t_max=3.0e7)
                rows.append(
                    SpinRow(
                        ramp=ramp,
                        nu=nu,
                        endcap=endcap,
                        floor=FLOOR_LABEL,
                        time_to_90=time_to_corotation(run),
                        peak_wind_ref=peak_wind_at_height(run, DRAG_HEIGHT),
                        peak_wind_interior=float(np.abs(run.u).max()),
                        max_epsilon=run.max_epsilon,
                        sink_valid=run.sink_valid,
                    )
                )
                if (ramp, nu, endcap) == (RAMPS[0], 1.0, True):
                    snapshots = tuple(
                        Snapshot(
                            h,
                            run.radius,
                            run.u[int(np.argmin(np.abs(run.t - 3600.0 * h)))],
                        )
                        for h in SNAPSHOT_HOURS
                    )
                    hours, corot = run.t / 3600.0, run.corotation
    return tuple(rows), snapshots, hours, corot


def _sudden_rows(cyl: Cylinder) -> tuple[SuddenRow, ...]:
    rows = []
    for endcap in (True, False):
        for sign in (-1.0, 1.0):
            params = SpinupParams(1.0, nu_unstable=10.0, endcap=endcap, n_cells=100)
            run = sudden_change(cyl, params, sign * STEP * cyl.omega, t_max=3.0e7)
            rows.append(
                SuddenRow(
                    sign * STEP,
                    1.0,
                    10.0,
                    endcap,
                    time_to_adjust(run),
                    run.max_epsilon,
                    run.sink_valid,
                )
            )
    return tuple(rows)


def _blend_rows(cyl: Cylinder) -> tuple[BlendRow, ...]:
    rows = []
    for delta in (-STEP, STEP, -0.5):
        times, epsilons = [], []
        for width in BLEND_WIDTHS:
            params = SpinupParams(1.0, nu_unstable=10.0, n_cells=100, blend_width=width)
            run = sudden_change(
                cyl, params, delta * cyl.omega, t_max=3.0e7, dt_max=BLEND_DT_MAX
            )
            times.append(time_to_adjust(run))
            epsilons.append(run.max_epsilon)
        eps = max(epsilons)
        rows.append(
            BlendRow(delta, BLEND_WIDTHS, tuple(times), eps, eps <= LINEAR_SINK_EPSILON)
        )
    return tuple(rows)


def _lag_rows(cyl: Cylinder) -> tuple[LagRow, ...]:
    alpha = LAG_FRACTION_PER_DAY * cyl.omega / 86400.0
    rows = []
    for nu in EDDY_VISCOSITIES:
        for endcap in (True, False):
            lag = steady_lag(cyl, SpinupParams(nu, endcap=endcap, n_cells=200), alpha)
            rows.append(
                LagRow(
                    nu,
                    endcap,
                    FLOOR_LABEL,
                    lag.peak,
                    lag.wind_at_height(DRAG_HEIGHT),
                )
            )
    return tuple(rows)


def _molecular_rows(cyl: Cylinder, gases: tuple[Gas, ...]) -> tuple[MolecularRow, ...]:
    rows = []
    for g in gases:
        run = hull_spin_up(
            cyl,
            SpinupParams(g.nu, endcap=True, n_cells=100),
            cyl.omega / RAMPS[0],
            t_max=5.0e8,
            dt_max=86400.0,
        )
        rows.append(
            MolecularRow(
                g.label,
                g.nu,
                time_to_corotation(run),
                peak_wind_at_height(run, DRAG_HEIGHT),
                peak_wind_at_height(run, DRAG_HEIGHT) / math.sqrt(g.nu * cyl.omega),
            )
        )
    return tuple(rows)


def build_report(cyl: Cylinder | None = None) -> AtmosphereReport:
    """Compute every table of the atmosphere write-up.

    The result is cached per habitat; the default is 3.2 km, 20 km, 1 g.
    """
    return _build(Cylinder() if cyl is None else cyl)


@cache
def _build(habitat: Cylinder) -> AtmosphereReport:
    gases = tuple(g.at_temperature(FLOOR_TEMPERATURE) for g in PRESETS)
    surface, stability = _surface_rows(habitat, gases)
    spin, snapshots, hours, corot = _spin_rows(habitat)
    return AtmosphereReport(
        cylinder=habitat,
        gas=_gas_rows(gases),
        equilibrium=_equilibrium_rows(habitat, gases),
        numbers=tuple(
            NumbersRow(
                g.label,
                dimensionless_table(habitat, g),
                habitat.omega / EARTH_OMEGA,
            )
            for g in gases
        ),
        surface=surface,
        stability=stability,
        parcels=_parcel_rows(habitat, gases[0]),
        spin=spin,
        sudden=_sudden_rows(habitat),
        blend=_blend_rows(habitat),
        lag=_lag_rows(habitat),
        molecular=_molecular_rows(habitat, gases),
        snapshots=snapshots,
        corotation_hours=hours,
        corotation=corot,
        caveats=dict(CAVEATS),
    )
