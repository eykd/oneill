import math
from dataclasses import FrozenInstanceError, replace
from typing import Any

import pytest
from hypothesis import given
from hypothesis import strategies as st

from oneill.atmo.gas import ONEILL_1974, R_UNIVERSAL, SP413, Gas


class TestEarthAir:
    def test_specific_gas_constant(self) -> None:
        assert ONEILL_1974.r_specific == pytest.approx(287.05, rel=1e-4)

    def test_cp(self) -> None:
        assert ONEILL_1974.cp == pytest.approx(1005.0, rel=2e-3)

    def test_gamma_and_kappa(self) -> None:
        assert ONEILL_1974.gamma == pytest.approx(1.40, abs=2e-3)
        assert ONEILL_1974.kappa == pytest.approx(0.2857, abs=1e-3)
        assert ONEILL_1974.cv == pytest.approx(ONEILL_1974.cp - 287.05, rel=1e-3)

    def test_viscosity_at_288k(self) -> None:
        assert ONEILL_1974.viscosity(288.0) == pytest.approx(1.79e-5, rel=1e-2)

    def test_viscosity_defaults_to_floor_temperature(self) -> None:
        assert ONEILL_1974.viscosity() == ONEILL_1974.viscosity(288.15)

    def test_sea_level_density_and_nu(self) -> None:
        assert ONEILL_1974.rho_rim == pytest.approx(1.225, rel=1e-3)
        assert ONEILL_1974.nu == pytest.approx(1.46e-5, rel=1e-2)

    def test_molar_mass(self) -> None:
        assert ONEILL_1974.molar_mass == pytest.approx(28.965e-3, rel=1e-3)
        assert R_UNIVERSAL / ONEILL_1974.molar_mass == ONEILL_1974.r_specific


class TestSp413:
    def test_dry_partials_and_vapour(self) -> None:
        assert SP413.p_rim == 49.7e3
        assert SP413.p_vapour == 1.0e3
        assert SP413.p_total == pytest.approx(50.7e3)
        assert SP413.x_o2 * SP413.p_rim == pytest.approx(22.7e3, rel=1e-12)
        assert SP413.x_n2 * SP413.p_rim == pytest.approx(26.6e3, rel=1e-12)
        assert SP413.x_o2 == pytest.approx(0.45, abs=0.01)

    def test_is_half_an_atmosphere(self) -> None:
        assert SP413.p_total / ONEILL_1974.p_rim == pytest.approx(0.5, abs=0.01)

    def test_fractions_sum_to_one(self) -> None:
        assert SP413.x_o2 + SP413.x_n2 + SP413.x_co2 == pytest.approx(1.0, abs=1e-15)

    def test_half_atm_kinematic_viscosity_and_heat_capacity(self) -> None:
        # research table: nu ~ 3.0e-5 m^2/s, rho c_p ~ 603 J/(m^3 K) at 288 K
        gas = SP413.at_temperature(288.0)
        assert gas.nu == pytest.approx(3.0e-5, rel=0.04)
        assert gas.rho_rim * gas.cp == pytest.approx(603.0, rel=0.02)

    def test_o2_rich_mixture_is_not_less_viscous_than_air(self) -> None:
        # O2 is more viscous than N2, so the mixture mu exceeds air's; at equal
        # pressure the denser mixture's nu is about the same as air's.
        air = replace(ONEILL_1974, p_rim=SP413.p_rim, temperature=288.0)
        mix = SP413.at_temperature(288.0)
        assert mix.viscosity() > air.viscosity()
        assert mix.nu == pytest.approx(air.nu, rel=0.02)


class TestWilkeViscosity:
    def test_pure_species_limits(self) -> None:
        # reference viscosities of the Sutherland tabulation
        assert Gas(1.0, 0.0).viscosity(300.55) == pytest.approx(2.018e-5, rel=1e-12)
        assert Gas(0.0, 1.0).viscosity(300.55) == pytest.approx(1.781e-5, rel=1e-12)
        assert Gas(0.0, 0.0, 1.0).viscosity(273.15) == pytest.approx(
            2.125e-5, rel=1e-12
        )
        assert Gas(0.0, 0.0, 0.0, 1.0).viscosity(300.55) == pytest.approx(
            1.48e-5, rel=1e-12
        )

    def test_air_matches_published_values(self) -> None:
        # CRC / US Standard Atmosphere: 1.789e-5 at 288.15 K, 1.846e-5 at 300 K
        assert ONEILL_1974.viscosity(288.15) == pytest.approx(1.789e-5, rel=1e-2)
        assert ONEILL_1974.viscosity(300.0) == pytest.approx(1.846e-5, rel=1e-2)

    @given(st.floats(0.0, 1.0), st.floats(200.0, 500.0))
    def test_binary_mixture_lies_between_pure_species(self, x: float, t: float) -> None:
        mix = Gas(x, 1.0 - x).viscosity(t)
        pure = sorted((Gas(1.0, 0.0).viscosity(t), Gas(0.0, 1.0).viscosity(t)))
        assert pure[0] * (1 - 1e-12) <= mix <= pure[1] * (1 + 1e-12)

    def test_label(self) -> None:
        assert ONEILL_1974.label == "O'Neill 1974 at 288.1 K"
        assert SP413.at_temperature(290.0).label == "SP-413 at 290.0 K"
        assert Gas(0.2, 0.8).name == "custom"

    def test_negative_vapour_pressure_rejected(self) -> None:
        with pytest.raises(ValueError, match="p_vapour"):
            Gas(0.2, 0.8, p_vapour=-1.0)


class TestGasObject:
    def test_at_temperature_keeps_composition_and_pressure(self) -> None:
        hot = ONEILL_1974.at_temperature(300.0)
        assert hot.temperature == 300.0
        assert hot.p_rim == ONEILL_1974.p_rim
        assert hot.x_ar == ONEILL_1974.x_ar

    def test_frozen(self) -> None:
        with pytest.raises(FrozenInstanceError):
            ONEILL_1974.p_rim = 1.0  # type: ignore[misc]

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"x_o2": 0.5, "x_n2": 0.4},
            {"x_o2": -0.1, "x_n2": 1.1},
            {"x_o2": math.nan, "x_n2": 0.5},
            {"x_o2": 0.2, "x_n2": 0.8, "p_rim": 0.0},
            {"x_o2": 0.2, "x_n2": 0.8, "p_rim": math.inf},
            {"x_o2": 0.2, "x_n2": 0.8, "temperature": -1.0},
        ],
    )
    def test_rejects_invalid(self, kwargs: dict[str, Any]) -> None:
        with pytest.raises(ValueError, match="must"):
            Gas(**kwargs)

    def test_rejects_non_positive_viscosity_temperature(self) -> None:
        with pytest.raises(ValueError, match="temperature"):
            ONEILL_1974.viscosity(0.0)

    @given(st.floats(1e3, 1e6), st.floats(150.0, 400.0))
    def test_density_is_ideal_gas_law(self, p: float, t: float) -> None:
        rho = ONEILL_1974.density(p, t)
        assert p == pytest.approx(rho * ONEILL_1974.r_specific * t, rel=1e-12)

    @given(st.floats(150.0, 600.0))
    def test_viscosity_increases_with_temperature(self, t: float) -> None:
        assert ONEILL_1974.viscosity(t + 1.0) > ONEILL_1974.viscosity(t)
