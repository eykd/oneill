import math
from dataclasses import FrozenInstanceError, replace

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
        assert ONEILL_1974.viscosity(288.0) == pytest.approx(1.79e-5, rel=5e-3)

    def test_viscosity_defaults_to_floor_temperature(self) -> None:
        assert ONEILL_1974.viscosity() == ONEILL_1974.viscosity(288.15)

    def test_sea_level_density_and_nu(self) -> None:
        assert ONEILL_1974.rho_rim == pytest.approx(1.225, rel=1e-3)
        assert ONEILL_1974.nu == pytest.approx(1.46e-5, rel=1e-2)

    def test_molar_mass(self) -> None:
        assert ONEILL_1974.molar_mass == pytest.approx(28.965e-3, rel=1e-3)
        assert R_UNIVERSAL / ONEILL_1974.molar_mass == ONEILL_1974.r_specific


class TestSp413:
    def test_pressure_and_oxygen_partial_pressure(self) -> None:
        assert SP413.p_rim == 50.8e3
        assert SP413.x_o2 * 49.7 == pytest.approx(22.7, rel=1e-12)
        assert SP413.x_o2 == pytest.approx(0.45, abs=0.01)

    def test_is_half_an_atmosphere(self) -> None:
        assert SP413.p_rim / ONEILL_1974.p_rim == pytest.approx(0.5, abs=0.01)

    def test_half_atm_air_matches_research_viscosity(self) -> None:
        # The research table's 3.0e-5 m^2/s is *air* at half an atmosphere.
        half_atm_air = replace(ONEILL_1974, p_rim=50.8e3, temperature=288.0)
        assert half_atm_air.nu == pytest.approx(3.0e-5, rel=0.05)
        assert half_atm_air.rho_rim * half_atm_air.cp == pytest.approx(603.0, rel=0.03)

    def test_oxygen_rich_mixture_is_denser_so_less_viscous_than_air(self) -> None:
        # mu is the same Sutherland value; rho scales with molar mass.
        air = replace(ONEILL_1974, p_rim=50.8e3, temperature=288.0)
        mix = SP413.at_temperature(288.0)
        expected = air.nu * ONEILL_1974.molar_mass / SP413.molar_mass
        assert mix.nu == pytest.approx(expected, rel=1e-12)
        assert mix.nu == pytest.approx(2.81e-5, rel=0.01)

    def test_fractions_sum_to_one(self) -> None:
        assert SP413.x_o2 + SP413.x_n2 + SP413.x_co2 == pytest.approx(1.0, abs=1e-15)


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
    def test_rejects_invalid(self, kwargs: dict[str, float]) -> None:
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
