"""Ideal-gas mixture properties of the habitat air.

Everything here is a pure function of the composition (mole fractions), the rim
pressure and the temperature. Two presets cover the designs the research compares:
O'Neill's 1974 sea-level air and the NASA SP-413 Stanford-torus atmosphere
(``docs/research/atmosphere-literature.md`` Q1, ``model-parameters.md`` gas table).
"""

import math
from dataclasses import dataclass, replace

from oneill.core.checks import finite, positive

R_UNIVERSAL = 8.314462618
"""Molar gas constant in J/(mol K) (CODATA 2018 exact value)."""

_SUM_TOLERANCE = 1e-9


@dataclass(frozen=True, slots=True)
class _Species:
    """Pure-species data: molar mass, molar c_p and Sutherland viscosity constants."""

    molar_mass: float  # kg/mol
    cp_molar: float  # J/(mol K), ideal gas near 298 K (NIST/JANAF)
    mu_ref: float  # Pa s at t_ref
    t_ref: float  # K
    sutherland: float  # K

    def viscosity(self, t: float) -> float:
        """Sutherland's law for this species, in Pa s."""
        ratio = (self.t_ref + self.sutherland) / (t + self.sutherland)
        return float(self.mu_ref * (t / self.t_ref) ** 1.5 * ratio)


# Sutherland constants: White, Viscous Fluid Flow (N2, O2, CO2 at 300.55 K); argon
# from the usual tabulation (2.125e-5 Pa s at 273.15 K, S = 144 K).
_O2 = _Species(31.998e-3, 29.376, 2.018e-5, 300.55, 127.0)
_N2 = _Species(28.0134e-3, 29.124, 1.781e-5, 300.55, 111.0)
_AR = _Species(39.948e-3, 20.786, 2.125e-5, 273.15, 144.0)
_CO2 = _Species(44.0095e-3, 37.129, 1.480e-5, 300.55, 240.0)
_SPECIES = (_O2, _N2, _AR, _CO2)


def _wilke(x: tuple[float, ...], mu: tuple[float, ...], m: tuple[float, ...]) -> float:
    """Wilke (1950) mixing rule for the viscosity of a gas mixture."""
    total = 0.0
    for i, xi in enumerate(x):
        if xi == 0.0:
            continue
        denom = 0.0
        for j, xj in enumerate(x):
            phi = (1.0 + math.sqrt(mu[i] / mu[j]) * (m[j] / m[i]) ** 0.25) ** 2 / (
                math.sqrt(8.0 * (1.0 + m[i] / m[j]))
            )
            denom += xj * phi
        total += xi * mu[i] / denom
    return total


@dataclass(frozen=True, slots=True)
class Gas:
    """A dry ideal-gas mixture at the rim pressure and temperature.

    Water vapour is not a component of the mixture; it is a passive trace carried
    in ``p_vapour`` and enters the dynamics through the humidity in
    :mod:`oneill.atmo.equilibrium`.

    Attributes:
        x_o2: Mole fraction of oxygen.
        x_n2: Mole fraction of nitrogen.
        x_ar: Mole fraction of argon.
        x_co2: Mole fraction of carbon dioxide.
        p_rim: Pressure of the dry mixture at the floor (rim), in Pa.
        temperature: Air temperature at the floor, in K.
        p_vapour: Water-vapour partial pressure at the floor, in Pa. Carried for
            reporting (:attr:`p_total`) but not in ``R_s``, ``c_p`` or the
            density, a bias of about 2 % for SP-413.
        name: Label carried into reports (``"custom"`` unless a preset).
    """

    x_o2: float
    x_n2: float
    x_ar: float = 0.0
    x_co2: float = 0.0
    p_rim: float = 101325.0
    temperature: float = 288.15
    p_vapour: float = 0.0
    name: str = "custom"

    def __post_init__(self) -> None:
        """Validate fractions (non-negative, summing to one), pressures and T."""
        for name, value in (
            ("x_o2", self.x_o2),
            ("x_n2", self.x_n2),
            ("x_ar", self.x_ar),
            ("x_co2", self.x_co2),
        ):
            if not (math.isfinite(value) and value >= 0.0):
                msg = f"{name} must be non-negative and finite, got {value!r}"
                raise ValueError(msg)
        total = self.x_o2 + self.x_n2 + self.x_ar + self.x_co2
        if abs(total - 1.0) > _SUM_TOLERANCE:
            msg = f"mole fractions must sum to 1, got {total!r}"
            raise ValueError(msg)
        positive("p_rim", self.p_rim)
        positive("temperature", self.temperature)
        finite("p_vapour", self.p_vapour)
        if self.p_vapour < 0.0:
            msg = f"p_vapour must be non-negative, got {self.p_vapour!r}"
            raise ValueError(msg)

    def at_temperature(self, temperature: float) -> Gas:
        """Return the same mixture and rim pressure at another temperature (K)."""
        return replace(self, temperature=temperature)

    @property
    def label(self) -> str:
        """Preset name and floor temperature, for labelling reported numbers."""
        return f"{self.name} at {self.temperature:.1f} K"

    @property
    def p_total(self) -> float:
        """Dry-mixture plus water-vapour pressure at the floor, in Pa."""
        return self.p_rim + self.p_vapour

    @property
    def _fractions(self) -> tuple[float, float, float, float]:
        return (self.x_o2, self.x_n2, self.x_ar, self.x_co2)

    @property
    def molar_mass(self) -> float:
        """Mean molar mass in kg/mol."""
        return sum(
            x * sp.molar_mass for x, sp in zip(self._fractions, _SPECIES, strict=True)
        )

    @property
    def r_specific(self) -> float:
        """Specific gas constant R/M in J/(kg K)."""
        return R_UNIVERSAL / self.molar_mass

    @property
    def cp(self) -> float:
        """Specific heat at constant pressure in J/(kg K) (mole-weighted)."""
        molar = sum(
            x * sp.cp_molar for x, sp in zip(self._fractions, _SPECIES, strict=True)
        )
        return molar / self.molar_mass

    @property
    def cv(self) -> float:
        """Specific heat at constant volume in J/(kg K)."""
        return self.cp - self.r_specific

    @property
    def gamma(self) -> float:
        """Ratio of specific heats c_p / c_v."""
        return self.cp / self.cv

    @property
    def kappa(self) -> float:
        """Poisson exponent R_s / c_p."""
        return self.r_specific / self.cp

    def density(self, pressure: float, temperature: float) -> float:
        """Ideal-gas density p / (R_s T) in kg/m^3."""
        return pressure / (self.r_specific * temperature)

    @property
    def rho_rim(self) -> float:
        """Density at the floor in kg/m^3."""
        return self.density(self.p_rim, self.temperature)

    def viscosity(self, temperature: float | None = None) -> float:
        """Dynamic viscosity in Pa s: Wilke's rule over per-species Sutherland laws.

        Args:
            temperature: Temperature in K; defaults to the floor temperature.
        """
        t = self.temperature
        if temperature is not None:
            t = positive("temperature", temperature)
        return _wilke(
            self._fractions,
            tuple(sp.viscosity(t) for sp in _SPECIES),
            tuple(sp.molar_mass for sp in _SPECIES),
        )

    @property
    def nu(self) -> float:
        """Kinematic viscosity at the floor in m^2/s."""
        return self.viscosity() / self.rho_rim


def _sp413_fractions() -> tuple[float, float, float]:
    """Dry mole fractions (O2, N2, CO2) from SP-413 Table 4-3 partial pressures."""
    p_o2, p_n2, p_co2 = 22.7, 26.6, 0.4  # kPa; CO2 at its "< 0.4" bound
    total = p_o2 + p_n2 + p_co2
    return p_o2 / total, p_n2 / total, p_co2 / total


ONEILL_1974 = Gas(
    name="O'Neill 1974",
    x_o2=0.2095,
    x_n2=0.7808,
    x_ar=0.0093,
    x_co2=0.0004,
    p_rim=101325.0,
    temperature=288.15,
)
"""O'Neill 1974: sea-level Earth air at 1 atm (Physics Today 27(9):32; Q1)."""

_X_O2, _X_N2, _X_CO2 = _sp413_fractions()
SP413 = Gas(
    name="SP-413",
    x_o2=_X_O2,
    x_n2=_X_N2,
    x_ar=0.0,
    x_co2=_X_CO2,
    p_rim=49.7e3,
    temperature=293.15,
    p_vapour=1.0e3,
)
"""NASA SP-413 Table 4-3 Stanford-torus air at 20 degC (Q1).

The dry mixture is the table's partial pressures (O2 22.7, N2 26.6, CO2 0.4 kPa,
the last at its upper bound), 49.7 kPa; the table's 1.0 kPa of water vapour is
carried separately in ``p_vapour`` (total 50.7 kPa against the quoted 50.8).
"""
