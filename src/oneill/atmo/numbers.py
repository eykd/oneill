"""L0: the dimensionless numbers that say which physics matters in the habitat.

Rotation enters every estimate through the Coriolis parameter ``f = 2 omega``
(``docs/research/atmosphere-literature.md``; ``model-parameters.md`` L0 tables).
Greenspan-Howard (1963) and Stewartson (1957) lengths use the *half-length*
``H = L / 2`` of the cylinder, because spin-up proceeds from both end caps.

All functions are pure; :func:`dimensionless_table` bundles the headline numbers
into frozen dataclasses for the reporting layers to render.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from collections.abc import Sequence

    from oneill.atmo.gas import Gas
    from oneill.core.cylinder import Cylinder

STROUHAL_BLUFF_BODY = 0.2
"""Strouhal number for the shedding estimate (round bodies; buildings 0.12-0.14)."""


@dataclass(frozen=True, slots=True)
class ViscousScales:
    """Everything that depends on the (molecular or eddy) viscosity.

    Attributes:
        nu: Kinematic viscosity in m^2/s.
        ekman_number: ``nu / (omega H^2)`` with ``H`` the half-length.
        ekman_depth: ``sqrt(nu / omega)`` in m.
        stewartson_quarter: ``H E^(1/4)`` in m.
        stewartson_third: ``H E^(1/3)`` in m.
        spin_up_time: Greenspan-Howard time ``H / sqrt(nu omega)`` in s.
    """

    nu: float
    ekman_number: float
    ekman_depth: float
    stewartson_quarter: float
    stewartson_third: float
    spin_up_time: float


@dataclass(frozen=True, slots=True)
class DimensionlessTable:
    """Headline dimensionless numbers for one habitat, gas and flow.

    Attributes:
        omega: Spin rate in rad/s.
        coriolis_parameter: ``2 omega`` in 1/s.
        inertial_period: ``pi / omega`` in s.
        rossby: ``U / (2 omega L)`` for the flow scale.
        obstacle_rossby: ``U / (2 omega D)`` for the obstacle diameter.
        brunt_vaisala: ``N`` in 1/s at the chosen height, isothermal atmosphere.
        rossby_radius: ``N H / (2 omega)`` in m.
        shedding_frequency: Von Karman shedding frequency ``St U / D`` in Hz.
        shedding_to_coriolis: ``2 pi f_shed / (2 omega)``.
        molecular: Viscous scales for the gas's molecular viscosity.
        eddy: Viscous scales for each requested eddy viscosity.
    """

    omega: float
    coriolis_parameter: float
    inertial_period: float
    rossby: float
    obstacle_rossby: float
    brunt_vaisala: float
    rossby_radius: float
    shedding_frequency: float
    shedding_to_coriolis: float
    molecular: ViscousScales
    eddy: tuple[ViscousScales, ...]


def rossby_number(speed: float, length: float, omega: float) -> float:
    """Rossby number ``U / (2 omega L)``: inertia over Coriolis."""
    positive("length", length)
    positive("omega", omega)
    return finite("speed", speed) / (2.0 * omega * length)


def ekman_number(nu: float, omega: float, half_length: float) -> float:
    """Ekman number ``nu / (omega H^2)`` (Greenspan 1968; ``H`` the half-length)."""
    return positive("nu", nu) / (
        positive("omega", omega) * positive("half_length", half_length) ** 2
    )


def ekman_depth(nu: float, omega: float) -> float:
    """Ekman layer depth ``sqrt(nu / omega)`` in m.

    This is the textbook ``sqrt(2 nu / f)`` with the Coriolis parameter
    ``f = 2 omega``; it reproduces the research table (0.023 m molecular at half
    an atmosphere, 4.2 m for ``nu_t = 1``, 13 m for ``nu_t = 10``). Writing
    ``sqrt(2 nu / omega)`` would be the same formula with ``f`` taken as
    ``omega``, a factor sqrt(2) too deep.
    """
    return math.sqrt(positive("nu", nu) / positive("omega", omega))


def stewartson_quarter_thickness(nu: float, omega: float, half_length: float) -> float:
    """Stewartson ``E^(1/4)`` layer thickness ``H E^(1/4)`` (Stewartson 1957)."""
    return float(half_length * ekman_number(nu, omega, half_length) ** 0.25)


def stewartson_third_thickness(nu: float, omega: float, half_length: float) -> float:
    """Stewartson ``E^(1/3)`` layer thickness ``H E^(1/3)`` (Stewartson 1957)."""
    return float(half_length * ekman_number(nu, omega, half_length) ** (1.0 / 3.0))


def spin_up_time(nu: float, omega: float, half_length: float) -> float:
    """Greenspan-Howard (1963) spin-up time ``H / sqrt(nu omega)`` in s.

    ``H`` is the half-length ``L / 2`` (10 km for the default habitat).
    """
    positive("half_length", half_length)
    return half_length / math.sqrt(positive("nu", nu) * positive("omega", omega))


def inertial_period(omega: float) -> float:
    """Period ``pi / omega`` of inertial oscillations (frequency ``2 omega``), in s."""
    return math.pi / positive("omega", omega)


def brunt_vaisala_squared(
    g_local: float, temperature: float, dt_dh: float, cp: float
) -> float:
    """Squared buoyancy frequency ``N^2 = (g/T)(g/c_p + dT/dh)`` in 1/s^2.

    ``g`` is the local centrifugal gravity ``omega^2 r`` and ``h`` the height
    above the floor; ``N^2`` is zero on the dry adiabat and negative where the
    lapse rate is superadiabatic.
    """
    positive("g_local", g_local)
    positive("temperature", temperature)
    positive("cp", cp)
    return g_local / temperature * (g_local / cp + finite("dt_dh", dt_dh))


def brunt_vaisala(g_local: float, temperature: float, dt_dh: float, cp: float) -> float:
    """Buoyancy frequency ``N`` in 1/s (see :func:`brunt_vaisala_squared`).

    Raises:
        ValueError: if the layer is convectively unstable (``N^2 < 0``).
    """
    n2 = brunt_vaisala_squared(g_local, temperature, dt_dh, cp)
    if n2 < 0.0:
        msg = f"layer is convectively unstable (N^2 = {n2!r})"
        raise ValueError(msg)
    return math.sqrt(n2)


def rossby_deformation_radius(n: float, depth: float, omega: float) -> float:
    """Rossby deformation radius ``N H / (2 omega)`` in m."""
    return finite("n", n) * positive("depth", depth) / (2.0 * positive("omega", omega))


def shedding_frequency(
    speed: float, diameter: float, strouhal: float = STROUHAL_BLUFF_BODY
) -> float:
    """Von Karman shedding frequency ``St U / D`` in Hz (Okajima 1982 for ``St``)."""
    return (
        positive("strouhal", strouhal)
        * positive("speed", speed)
        / positive("diameter", diameter)
    )


def viscous_scales(nu: float, omega: float, half_length: float) -> ViscousScales:
    """Bundle the viscosity-dependent scales for ``nu``."""
    return ViscousScales(
        nu=nu,
        ekman_number=ekman_number(nu, omega, half_length),
        ekman_depth=ekman_depth(nu, omega),
        stewartson_quarter=stewartson_quarter_thickness(nu, omega, half_length),
        stewartson_third=stewartson_third_thickness(nu, omega, half_length),
        spin_up_time=spin_up_time(nu, omega, half_length),
    )


def dimensionless_table(
    cyl: Cylinder,
    gas: Gas,
    *,
    wind_speed: float = 3.0,
    length_scale: float = 1000.0,
    obstacle_diameter: float = 10.0,
    nu_eddy: Sequence[float] = (1.0, 10.0),
    height: float = 0.0,
    depth: float | None = None,
) -> DimensionlessTable:
    """Compute the headline dimensionless numbers.

    The defaults are illustrative (a 3 m/s wind, a 1 km flow scale, a 10 m
    building). The buoyancy frequency is for an isothermal atmosphere at the gas
    temperature, ``N^2 = g^2 / (c_p T)``, at ``height`` above the floor.

    Args:
        cyl: The habitat.
        gas: Floor gas state, supplying ``nu``, ``c_p`` and ``T``.
        wind_speed: Flow speed ``U`` in m/s.
        length_scale: Horizontal flow scale ``L`` in m for the Rossby number.
        obstacle_diameter: Obstacle size ``D`` in m for shedding and obstacle Ro.
        nu_eddy: Eddy viscosities in m^2/s to tabulate.
        height: Height above the floor in m at which to evaluate ``N``.
        depth: Vertical scale ``H`` for the Rossby radius; defaults to the
            floor-to-axis depth ``R``.
    """
    omega = cyl.omega
    half_length = cyl.length / 2.0
    f_shed = shedding_frequency(wind_speed, obstacle_diameter)
    n = brunt_vaisala(cyl.gravity_at_height(height), gas.temperature, 0.0, gas.cp)
    return DimensionlessTable(
        omega=omega,
        coriolis_parameter=2.0 * omega,
        inertial_period=inertial_period(omega),
        rossby=rossby_number(wind_speed, length_scale, omega),
        obstacle_rossby=rossby_number(wind_speed, obstacle_diameter, omega),
        brunt_vaisala=n,
        rossby_radius=rossby_deformation_radius(
            n, cyl.radius if depth is None else depth, omega
        ),
        shedding_frequency=f_shed,
        shedding_to_coriolis=2.0 * math.pi * f_shed / (2.0 * omega),
        molecular=viscous_scales(gas.nu, omega, half_length),
        eddy=tuple(viscous_scales(nu, omega, half_length) for nu in nu_eddy),
    )
