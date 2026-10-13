"""Figure specifications for the atmosphere report (data only; nothing draws)."""

from typing import TYPE_CHECKING

import numpy as np

from oneill.atmo.equilibrium import (
    adiabatic_profile,
    isothermal_profile,
    lifting_condensation_level,
)
from oneill.atmo.report import (
    HUMIDITIES,
    TERRAINS,
    WIND_HEIGHT,
    WIND_SPEED,
)
from oneill.atmo.surface import (
    SUBURBAN,
    friction_velocity,
    unstable_layer,
    wind_shear,
    wind_speed,
)
from oneill.core.figspec import Band, FigureSpec, Panel, Series

if TYPE_CHECKING:
    from oneill.atmo.gas import Gas
    from oneill.atmo.report import AtmosphereReport
    from oneill.core.cylinder import Cylinder

_SAMPLES = 200
_WIND_TOP = 60.0


def profiles_figure(cyl: Cylinder, gas: Gas) -> FigureSpec:
    """Pressure and temperature against height, with the cloud-base levels."""
    h = np.linspace(0.0, cyl.radius, _SAMPLES)
    iso = isothermal_profile(cyl, gas, h)
    adi = adiabatic_profile(cyl, gas, h)
    clouds = tuple(
        c.height
        for c in (lifting_condensation_level(cyl, gas, rh) for rh in HUMIDITIES)
        if c.height is not None
    )
    pressure = Panel(
        "p / p_rim",
        "height above floor [m]",
        (
            Series("isothermal", iso.pressure / gas.p_rim, h),
            Series("dry adiabatic", adi.pressure / gas.p_rim, h, "--"),
        ),
        title="pressure",
        hlines=clouds,
        legend_size=8.0,
    )
    temperature = Panel(
        "temperature [K]",
        "height above floor [m]",
        (
            Series("isothermal", iso.temperature, h),
            Series("dry adiabatic", adi.temperature, h, "--"),
        ),
        title="temperature (cloud base for RH 50 %, 70 % as lines)",
        hlines=clouds,
        legend_size=8.0,
    )
    return FigureSpec(
        "atmo_1_profiles.png",
        (pressure, temperature),
        (1, 2),
        (10.0, 5.0),
        suptitle=f"Hydrostatic profiles: {gas.label}",
    )


def wind_figure(cyl: Cylinder) -> FigureSpec:
    """Wind profile over three terrains, with the spinward unstable layer."""
    z = np.linspace(0.5, _WIND_TOP, _SAMPLES)
    speeds, shears = [], []
    for terrain in TERRAINS:
        ustar = friction_velocity(WIND_SPEED, WIND_HEIGHT, terrain)
        speeds.append(Series(terrain.name, wind_speed(z, ustar, terrain), z))
        shears.append(Series(terrain.name, z, wind_shear(z, ustar, terrain)))
    ustar = friction_velocity(WIND_SPEED, WIND_HEIGHT, SUBURBAN)
    layer = unstable_layer(ustar, SUBURBAN, cyl.omega)
    bands = () if layer is None else (Band(*layer, "suburb: spinward wind unstable"),)
    speed_panel = Panel(
        "wind speed [m/s]",
        "height above floor [m]",
        tuple(speeds),
        title=f"{WIND_SPEED:g} m/s at {WIND_HEIGHT:g} m",
        bands=bands,
        legend_size=8.0,
    )
    shear_panel = Panel(
        "height above floor [m]",
        "dU/dz [1/s]",
        tuple(shears),
        title="shear against 2 omega (line)",
        hlines=(2.0 * cyl.omega,),
        legend_size=8.0,
    )
    return FigureSpec(
        "atmo_2_wind.png",
        (speed_panel, shear_panel),
        (1, 2),
        (10.0, 5.0),
        suptitle="Surface-layer wind and rotational stability (area-averaged canopy)",
    )


def parcel_figure(report: AtmosphereReport) -> FigureSpec:
    """Warm-parcel trajectories (a ceiling), for the +3 K release."""
    rows = [p for p in report.parcels if p.temperature_excess == 3.0]  # noqa: PLR2004
    rise = Panel(
        "time [s]",
        "height gained [m]",
        tuple(
            Series(f"k = {p.drag:.3f} /s", p.t, p.height - p.height[0]) for p in rows
        ),
        title="rise",
        legend_size=8.0,
    )
    drift = Panel(
        "time [s]",
        "spinward displacement [m]",
        tuple(Series(f"k = {p.drag:.3f} /s", p.t, p.spinward) for p in rows),
        title="spinward drift",
        legend_size=8.0,
    )
    return FigureSpec(
        "atmo_3_parcel.png",
        (rise, drift),
        (1, 2),
        (10.0, 5.0),
        suptitle="+3 K parcel, no perturbation pressure (upper bound)",
    )


def spinup_figure(report: AtmosphereReport) -> FigureSpec:
    """Wind snapshots and co-rotation for the 1 h spin-up, nu_t = 1, end caps on."""
    snapshots = Panel(
        "radius [m]",
        "relative wind u [m/s]",
        tuple(Series(f"t = {s.hours:g} h", s.radius, s.u) for s in report.snapshots),
        title="u(r) snapshots (no-slip floor)",
        legend_size=8.0,
    )
    corotation = Panel(
        "time [h]",
        "co-rotation fraction",
        (Series("co-rotation", report.corotation_hours, report.corotation),),
        title="fraction of solid-body angular momentum",
        hlines=(0.9,),
        legend_size=8.0,
    )
    return FigureSpec(
        "atmo_4_spinup.png",
        (snapshots, corotation),
        (1, 2),
        (10.0, 5.0),
        suptitle="Hull spin-up over 1 h, nu_t = 1 m2/s, end caps on",
    )


def atmosphere_figures(report: AtmosphereReport, gas: Gas) -> tuple[FigureSpec, ...]:
    """All four atmosphere figures, in file order."""
    return (
        profiles_figure(report.cylinder, gas),
        wind_figure(report.cylinder),
        parcel_figure(report),
        spinup_figure(report),
    )
