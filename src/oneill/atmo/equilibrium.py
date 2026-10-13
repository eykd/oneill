"""L0: hydrostatic equilibrium of the air in a spinning cylinder.

The working coordinate is the centrifugal geopotential per unit mass,
``phi = omega^2 (R^2 - r^2) / 2``, which plays the role of ``g z`` on Earth:
hydrostatic balance is ``dp = -rho dphi`` and the dry adiabat is
``T = T_rim - phi / c_p``. Every Earth formula therefore carries over exactly when
written in ``phi`` (``docs/research/atmosphere-literature.md`` Q1; plan Part 2). A
height ``h`` above the floor has ``r = R - h``.

Habitat height and *Earth-equivalent* height ``phi / g`` differ because the
centrifugal gravity falls to zero at the axis, so a cloud base quoted with Earth
lapse rates (as O'Neill's was) must be converted through ``phi``.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import positive
from oneill.core.cylinder import G_STANDARD

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

    from oneill.atmo.gas import Gas
    from oneill.core.cylinder import Cylinder

_AXIS_SLACK = 1e-9


@dataclass(frozen=True, slots=True)
class AtmosphereProfile:
    """Thermodynamic state sampled at heights above the floor.

    Attributes:
        height: Height above the floor in m.
        phi: Centrifugal geopotential in J/kg.
        pressure: Pressure in Pa.
        temperature: Temperature in K.
        density: Density in kg/m^3.
    """

    height: NDArray[np.float64]
    phi: NDArray[np.float64]
    pressure: NDArray[np.float64]
    temperature: NDArray[np.float64]
    density: NDArray[np.float64]


@dataclass(frozen=True, slots=True)
class CloudBase:
    """Lifting condensation level of a rising dry-adiabatic parcel.

    Attributes:
        label: Gas preset and floor temperature the number belongs to.
        temperature: Parcel temperature at the LCL in K.
        phi: Geopotential of the LCL in J/kg.
        earth_equivalent_height: ``phi / g`` in m: the height at which Earth's
            gravity would give the same ``phi`` (what an Earth-lapse-rate
            estimate such as O'Neill's returns).
        height: True height above the floor in m, or ``None`` when the LCL lies
            above the axis (the parcel never condenses before reaching it).
    """

    label: str
    temperature: float
    phi: float
    earth_equivalent_height: float
    height: float | None


def _phi(cyl: Cylinder, height: ArrayLike) -> NDArray[np.float64]:
    """Geopotential ``omega^2 (2 R h - h^2) / 2`` for heights in ``[0, R]``."""
    h = np.asarray(height, dtype=np.float64)
    if np.any(h < 0.0) or np.any(h > cyl.radius * (1.0 + _AXIS_SLACK)):
        msg = f"height must lie between the floor and the axis (0..{cyl.radius} m)"
        raise ValueError(msg)
    return np.asarray(0.5 * cyl.omega**2 * (2.0 * cyl.radius * h - h * h))


def geopotential(cyl: Cylinder, height: float) -> float:
    """Geopotential ``omega^2 (R^2 - r^2) / 2`` at ``height`` metres above the floor.

    Raises:
        ValueError: if ``height`` is outside ``0..R``.
    """
    return float(_phi(cyl, height))


def axis_geopotential(cyl: Cylinder) -> float:
    """Geopotential of the axis, ``omega^2 R^2 / 2`` (15.7 kJ/kg at 3.2 km, 1 g)."""
    return 0.5 * cyl.omega**2 * cyl.radius**2


def height_from_geopotential(cyl: Cylinder, phi: float) -> float:
    """Invert :func:`geopotential`: ``h = R - sqrt(R^2 - 2 phi / omega^2)``.

    Raises:
        ValueError: if ``phi`` is negative or beyond the axis geopotential.
    """
    top = axis_geopotential(cyl)
    if not 0.0 <= phi <= top * (1.0 + _AXIS_SLACK):
        msg = f"phi must lie in 0..{top} J/kg (floor to axis), got {phi!r}"
        raise ValueError(msg)
    radicand = max(cyl.radius**2 - 2.0 * phi / cyl.omega**2, 0.0)
    return cyl.radius - math.sqrt(radicand)


def earth_equivalent_height(phi: float) -> float:
    """Height on Earth (standard gravity) with geopotential ``phi``: ``phi / g``."""
    return phi / G_STANDARD


def dry_lapse_rate(cyl: Cylinder, gas: Gas, height: float) -> float:
    """Dry-adiabatic lapse rate ``-dT/dh = omega^2 r / c_p`` in K/m.

    It is ``g / c_p`` (9.76 K/km for air at 1 g) at the floor and falls linearly
    to zero at the axis.
    """
    geopotential(cyl, height)  # validates the range
    return cyl.omega**2 * (cyl.radius - height) / gas.cp


def isothermal_profile(
    cyl: Cylinder, gas: Gas, heights: ArrayLike
) -> AtmosphereProfile:
    """Isothermal atmosphere: ``p = p_rim exp(-phi / (R_s T))``."""
    phi = _phi(cyl, heights)
    t = gas.temperature
    p = gas.p_rim * np.exp(-phi / (gas.r_specific * t))
    return AtmosphereProfile(
        height=np.asarray(heights, dtype=np.float64),
        phi=phi,
        pressure=p,
        temperature=np.full_like(phi, t),
        density=p / (gas.r_specific * t),
    )


def adiabatic_profile(cyl: Cylinder, gas: Gas, heights: ArrayLike) -> AtmosphereProfile:
    """Dry-adiabatic atmosphere: ``T = T_rim - phi / c_p``, Poisson for ``p``.

    Raises:
        ValueError: if the profile would reach absolute zero before the axis.
    """
    phi = _phi(cyl, heights)
    t = gas.temperature - phi / gas.cp
    if np.any(t <= 0.0):
        msg = "dry adiabat reaches absolute zero: habitat too deep for this gas"
        raise ValueError(msg)
    p = gas.p_rim * (t / gas.temperature) ** (gas.cp / gas.r_specific)
    return AtmosphereProfile(
        height=np.asarray(heights, dtype=np.float64),
        phi=phi,
        pressure=p,
        temperature=t,
        density=p / (gas.r_specific * t),
    )


def isothermal_axis_ratio(cyl: Cylinder, gas: Gas) -> float:
    """``p_axis / p_rim`` for the isothermal profile (0.828 for air at 290 K)."""
    return math.exp(-axis_geopotential(cyl) / (gas.r_specific * gas.temperature))


def adiabatic_axis_ratio(cyl: Cylinder, gas: Gas) -> float:
    """``p_axis / p_rim`` for the dry-adiabatic profile (0.823 for air at 290 K)."""
    theta = 1.0 - axis_geopotential(cyl) / (gas.cp * gas.temperature)
    if theta <= 0.0:
        msg = "dry adiabat reaches absolute zero: habitat too deep for this gas"
        raise ValueError(msg)
    return float(theta ** (gas.cp / gas.r_specific))


def adiabatic_temperature_drop(cyl: Cylinder, gas: Gas) -> float:
    """Dry-adiabatic temperature drop from floor to axis, ``phi_axis / c_p`` in K."""
    return axis_geopotential(cyl) / gas.cp


def total_mass(cyl: Cylinder, p_rim: float, axis_ratio: float) -> float:
    """Total air mass in kg for any hydrostatic profile with the given axis ratio.

    From ``dp = -rho dphi`` and the shell volume ``dV = 2 pi L dphi / omega^2``,
    the mass is ``2 pi L (p_rim - p_axis) / omega^2``, independent of the
    profile's shape. End-cap curvature is ignored.

    Args:
        cyl: The habitat.
        p_rim: Floor pressure in Pa.
        axis_ratio: ``p_axis / p_rim`` in ``(0, 1]``.
    """
    positive("p_rim", p_rim)
    _check_ratio(axis_ratio)
    return 2.0 * math.pi * cyl.length * p_rim * (1.0 - axis_ratio) / cyl.omega**2


def rim_pressure_from_mass(cyl: Cylinder, mass: float, axis_ratio: float) -> float:
    """Invert :func:`total_mass`: the floor pressure that holds ``mass`` kg of air."""
    positive("mass", mass)
    _check_ratio(axis_ratio)
    if axis_ratio == 1.0:
        msg = "axis_ratio of 1 carries no mass information (cylinder with no spin)"
        raise ValueError(msg)
    return mass * cyl.omega**2 / (2.0 * math.pi * cyl.length * (1.0 - axis_ratio))


def _check_ratio(axis_ratio: float) -> None:
    if not 0.0 < axis_ratio <= 1.0:
        msg = f"axis_ratio must lie in (0, 1], got {axis_ratio!r}"
        raise ValueError(msg)


def lifting_condensation_level(
    cyl: Cylinder, gas: Gas, relative_humidity: float
) -> CloudBase:
    """Cloud base of a floor parcel lifted dry-adiabatically, via Bolton's formula.

    Bolton (1980, Mon. Weather Rev. 108:1046, eq. 22) gives the LCL temperature
    ``T_L = 1 / (1 / (T - 55) - ln(RH) / 2840) + 55`` (K). The parcel cools
    dry-adiabatically, so the LCL geopotential is ``c_p (T - T_L)``, then
    converted to true height and to Earth-equivalent height.

    Bolton's constants embed air's ``R_s / c_p = 0.2854`` (and the Magnus
    saturation curve), so for the SP-413 mixture (``R_s / c_p`` about 0.284) the
    LCL temperature is used unchanged; the effect on the result is well under
    1 %. The result carries ``gas.label`` so the preset and temperature travel
    with the number.

    Args:
        cyl: The habitat.
        gas: Gas at the floor temperature (the parcel's starting point).
        relative_humidity: Floor relative humidity as a fraction in ``(0, 1]``.

    Raises:
        ValueError: if ``relative_humidity`` is outside ``(0, 1]``.
    """
    if not 0.0 < relative_humidity <= 1.0:
        msg = f"relative_humidity must lie in (0, 1], got {relative_humidity!r}"
        raise ValueError(msg)
    t = gas.temperature
    t_lcl = 1.0 / (1.0 / (t - 55.0) - math.log(relative_humidity) / 2840.0) + 55.0
    phi = gas.cp * (t - t_lcl)
    height = (
        height_from_geopotential(cyl, phi) if phi <= axis_geopotential(cyl) else None
    )
    return CloudBase(
        label=gas.label,
        temperature=t_lcl,
        phi=phi,
        earth_equivalent_height=earth_equivalent_height(phi),
        height=height,
    )
