"""Ideal-gas mixture properties of the habitat air.

Everything here is a pure function of the composition (mole fractions), the rim
pressure and the temperature. Two presets cover the designs the research compares:
O'Neill's 1974 sea-level air and the NASA SP-413 Stanford-torus atmosphere
(``docs/research/atmosphere-literature.md`` Q1, ``model-parameters.md`` gas table).
"""

import math
from dataclasses import dataclass, replace

from oneill.core.checks import positive

R_UNIVERSAL = 8.314462618
"""Molar gas constant in J/(mol K) (CODATA 2018 exact value)."""

_SUM_TOLERANCE = 1e-9

# Molar mass (kg/mol) and molar c_p (J/(mol K), ideal gas near 298 K); NIST
# Chemistry WebBook / JANAF tabulations.
_M_O2, _CP_O2 = 31.998e-3, 29.376
_M_N2, _CP_N2 = 28.0134e-3, 29.124
_M_AR, _CP_AR = 39.948e-3, 20.786
_M_CO2, _CP_CO2 = 44.0095e-3, 37.129

# Sutherland's law with the standard air calibration (reference viscosity at
# 273.15 K and Sutherland temperature).
_MU_REF = 1.716e-5
_T_REF = 273.15
_SUTHERLAND_S = 110.4


@dataclass(frozen=True, slots=True)
class Gas:
    """A dry ideal-gas mixture at the rim pressure and temperature.

    Water vapour is not a component here; humidity enters through the moist
    thermodynamics in :mod:`oneill.atmo.equilibrium`.

    Attributes:
        x_o2: Mole fraction of oxygen.
        x_n2: Mole fraction of nitrogen.
        x_ar: Mole fraction of argon.
        x_co2: Mole fraction of carbon dioxide.
        p_rim: Total pressure at the floor (rim), in Pa.
        temperature: Air temperature at the floor, in K.
    """

    x_o2: float
    x_n2: float
    x_ar: float = 0.0
    x_co2: float = 0.0
    p_rim: float = 101325.0
    temperature: float = 288.15

    def __post_init__(self) -> None:
        """Validate fractions (non-negative, summing to one), pressure and T."""
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

    def at_temperature(self, temperature: float) -> Gas:
        """Return the same mixture and rim pressure at another temperature (K)."""
        return replace(self, temperature=temperature)

    @property
    def molar_mass(self) -> float:
        """Mean molar mass in kg/mol."""
        return (
            self.x_o2 * _M_O2
            + self.x_n2 * _M_N2
            + self.x_ar * _M_AR
            + self.x_co2 * _M_CO2
        )

    @property
    def r_specific(self) -> float:
        """Specific gas constant R/M in J/(kg K)."""
        return R_UNIVERSAL / self.molar_mass

    @property
    def cp(self) -> float:
        """Specific heat at constant pressure in J/(kg K) (mole-weighted)."""
        molar = (
            self.x_o2 * _CP_O2
            + self.x_n2 * _CP_N2
            + self.x_ar * _CP_AR
            + self.x_co2 * _CP_CO2
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
        """Dynamic viscosity in Pa s from Sutherland's law.

        The constants are the standard *air* calibration; they are applied to any
        mixture here because O2, N2 and Ar viscosities differ by under 15 % and the
        SP-413 mixture is O2/N2.

        Args:
            temperature: Temperature in K; defaults to the floor temperature.
        """
        t = self.temperature
        if temperature is not None:
            t = positive("temperature", temperature)
        ratio = (_T_REF + _SUTHERLAND_S) / (t + _SUTHERLAND_S)
        return float(_MU_REF * (t / _T_REF) ** 1.5 * ratio)

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
    x_o2=_X_O2,
    x_n2=_X_N2,
    x_ar=0.0,
    x_co2=_X_CO2,
    p_rim=50.8e3,
    temperature=293.15,
)
"""NASA SP-413 Table 4-3 Stanford-torus air: 50.8 kPa total, 20 degC (Q1).

The dry mole fractions are the table's partial pressures (O2 22.7, N2 26.6, CO2
0.4 kPa, the last at its upper bound) normalised; the 1.0 kPa of water vapour is
carried by the relative humidity in the moist profile, not by this dry mixture.
"""
