import math
from dataclasses import replace
from typing import TYPE_CHECKING

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.atmo.equilibrium import (
    AtmosphereProfile,
    adiabatic_axis_ratio,
    adiabatic_profile,
    adiabatic_temperature_drop,
    axis_geopotential,
    dry_lapse_rate,
    earth_equivalent_height,
    geopotential,
    height_from_geopotential,
    isothermal_axis_ratio,
    isothermal_profile,
    lifting_condensation_level,
    rim_pressure_from_mass,
    total_mass,
)
from oneill.atmo.gas import ONEILL_1974, SP413, Gas
from oneill.core.cylinder import G_STANDARD, Cylinder

if TYPE_CHECKING:
    from collections.abc import Callable


AIR_290 = ONEILL_1974.at_temperature(290.0)
SP_290 = SP413.at_temperature(290.0)
G = G_STANDARD


def earth_isothermal(gas: Gas, z: np.ndarray) -> np.ndarray:
    return np.asarray(gas.p_rim * np.exp(-G * z / (gas.r_specific * gas.temperature)))


def earth_adiabatic(gas: Gas, z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    t = gas.temperature - G * z / gas.cp
    return gas.p_rim * (t / gas.temperature) ** (gas.cp / gas.r_specific), t


def habitat_z(cyl: Cylinder, h: np.ndarray) -> np.ndarray:
    """Earth-equivalent height of habitat height ``h``, derived independently."""
    return np.asarray(cyl.omega**2 * (2.0 * cyl.radius * h - h**2) / 2.0 / G)


class TestGeopotential:
    def test_default_depth(self, default_cyl: Cylinder) -> None:
        assert axis_geopotential(default_cyl) == pytest.approx(15.7e3, rel=2e-3)
        assert earth_equivalent_height(axis_geopotential(default_cyl)) == (
            pytest.approx(1600.0, rel=0.01)
        )

    def test_floor_and_axis(self, cyl: Cylinder) -> None:
        assert geopotential(cyl, 0.0) == 0.0
        assert geopotential(cyl, cyl.radius) == pytest.approx(axis_geopotential(cyl))

    def test_small_height_is_g_times_height(self, cyl: Cylinder) -> None:
        h = 1e-6 * cyl.radius
        assert geopotential(cyl, h) == pytest.approx(cyl.surface_gravity * h, rel=1e-6)

    @given(st.floats(0.0, 1.0))
    @example(0.0)
    @example(1.0)
    def test_round_trip(self, fraction: float) -> None:
        cyl = Cylinder()
        h = fraction * cyl.radius
        assert height_from_geopotential(cyl, geopotential(cyl, h)) == pytest.approx(
            h, abs=1e-6
        )

    @pytest.mark.parametrize("height", [-1.0, 3201.0])
    def test_height_out_of_range(self, default_cyl: Cylinder, height: float) -> None:
        with pytest.raises(ValueError, match="between the floor and the axis"):
            geopotential(default_cyl, height)

    @pytest.mark.parametrize("phi", [-1.0, 16000.0])
    def test_phi_out_of_range(self, default_cyl: Cylinder, phi: float) -> None:
        with pytest.raises(ValueError, match="floor to axis"):
            height_from_geopotential(default_cyl, phi)

    def test_axis_phi_with_rounding_slack(self, default_cyl: Cylinder) -> None:
        top = axis_geopotential(default_cyl)
        assert height_from_geopotential(default_cyl, top * (1 + 1e-12)) == (
            pytest.approx(3200.0)
        )


class TestAxisRatios:
    def test_isothermal_default(self, default_cyl: Cylinder) -> None:
        assert isothermal_axis_ratio(default_cyl, AIR_290) == pytest.approx(
            0.83, abs=0.005
        )

    def test_adiabatic_default(self, default_cyl: Cylinder) -> None:
        assert adiabatic_axis_ratio(default_cyl, AIR_290) == pytest.approx(
            0.822, abs=0.003
        )

    def test_adiabatic_temperature_drop(self, default_cyl: Cylinder) -> None:
        assert adiabatic_temperature_drop(default_cyl, AIR_290) == pytest.approx(
            15.6, abs=0.1
        )

    def test_oneill_quoted_ratio(self, default_cyl: Cylinder) -> None:
        # O'Neill 1974 quotes p_axis/p_rim ~ 0.825 for sea-level air.
        ratio = isothermal_axis_ratio(default_cyl, ONEILL_1974)
        assert ratio == pytest.approx(0.825, abs=0.005)

    @given(st.floats(2e3, 2e5))
    def test_ratio_independent_of_rim_pressure(self, p_rim: float) -> None:
        cyl = Cylinder()
        gas = replace(AIR_290, p_rim=p_rim)
        assert isothermal_axis_ratio(cyl, gas) == isothermal_axis_ratio(cyl, AIR_290)
        assert adiabatic_axis_ratio(cyl, gas) == adiabatic_axis_ratio(cyl, AIR_290)
        iso = isothermal_profile(cyl, gas, [0.0, cyl.radius]).pressure
        adi = adiabatic_profile(cyl, gas, [0.0, cyl.radius]).pressure
        assert iso[1] / iso[0] == pytest.approx(isothermal_axis_ratio(cyl, gas))
        assert adi[1] / adi[0] == pytest.approx(adiabatic_axis_ratio(cyl, gas))

    def test_adiabatic_ratio_too_deep(self) -> None:
        deep = Cylinder.from_gravity(radius=1.0e5)
        with pytest.raises(ValueError, match="absolute zero"):
            adiabatic_axis_ratio(deep, AIR_290)
        with pytest.raises(ValueError, match="absolute zero"):
            adiabatic_profile(deep, AIR_290, [0.0, 1.0e5])


class TestProfiles:
    heights_fraction = np.linspace(0.0, 1.0, 41)

    def test_floor_state_is_the_gas(self, cyl: Cylinder) -> None:
        for build in (isothermal_profile, adiabatic_profile):
            prof = build(cyl, AIR_290, [0.0])
            assert prof.pressure[0] == AIR_290.p_rim
            assert prof.temperature[0] == AIR_290.temperature
            assert prof.density[0] == pytest.approx(AIR_290.rho_rim)

    @pytest.mark.parametrize("build", [isothermal_profile, adiabatic_profile])
    def test_hydrostatic_balance_dp_dr(
        self,
        cyl: Cylinder,
        build: Callable[[Cylinder, Gas, np.ndarray], AtmosphereProfile],
    ) -> None:
        r = cyl.radius * np.linspace(0.05, 0.95, 19)
        d = 1e-4 * cyl.radius
        outer = build(cyl, AIR_290, cyl.radius - (r + d)).pressure
        inner = build(cyl, AIR_290, cyl.radius - (r - d)).pressure
        dp_dr = (outer - inner) / (2.0 * d)
        mid = build(cyl, AIR_290, cyl.radius - r).density
        assert dp_dr == pytest.approx(mid * cyl.omega**2 * r, rel=1e-6)

    def test_isothermal_is_isothermal(self, cyl: Cylinder) -> None:
        prof = isothermal_profile(cyl, AIR_290, self.heights_fraction * cyl.radius)
        assert np.all(prof.temperature == 290.0)

    def test_adiabatic_cools_to_axis(self, cyl: Cylinder) -> None:
        prof = adiabatic_profile(cyl, AIR_290, self.heights_fraction * cyl.radius)
        assert np.all(np.diff(prof.temperature) < 0.0)
        assert np.all(np.diff(prof.pressure) < 0.0)
        drop = prof.temperature[0] - prof.temperature[-1]
        assert drop == pytest.approx(adiabatic_temperature_drop(cyl, AIR_290))

    @pytest.mark.parametrize("gas", [AIR_290, SP_290])
    def test_matches_earth_isothermal_at_equivalent_height(
        self, cyl: Cylinder, gas: Gas
    ) -> None:
        h = self.heights_fraction * cyl.radius
        prof = isothermal_profile(cyl, gas, h)
        assert prof.pressure == pytest.approx(
            earth_isothermal(gas, habitat_z(cyl, h)), rel=1e-12
        )

    @pytest.mark.parametrize("gas", [AIR_290, SP_290])
    def test_matches_earth_adiabatic_at_equivalent_height(
        self, cyl: Cylinder, gas: Gas
    ) -> None:
        h = self.heights_fraction * cyl.radius
        prof = adiabatic_profile(cyl, gas, h)
        p_earth, t_earth = earth_adiabatic(gas, habitat_z(cyl, h))
        assert prof.pressure == pytest.approx(p_earth, rel=1e-12)
        assert prof.temperature == pytest.approx(t_earth, rel=1e-12)

    @given(st.floats(50.0, 20000.0), st.floats(0.0, 1.0))
    def test_earth_match_property_over_radius(
        self, radius: float, fraction: float
    ) -> None:
        cyl = Cylinder.from_gravity(radius=radius)
        h = np.array([fraction * radius])
        prof = adiabatic_profile(cyl, AIR_290, h)
        p_earth, _ = earth_adiabatic(AIR_290, habitat_z(cyl, h))
        assert prof.pressure == pytest.approx(p_earth, rel=1e-9)

    def test_density_is_ideal_gas(self, default_cyl: Cylinder) -> None:
        prof = adiabatic_profile(default_cyl, SP_290, [0.0, 1000.0, 3200.0])
        assert prof.density == pytest.approx(
            prof.pressure / (SP_290.r_specific * prof.temperature)
        )

    def test_height_out_of_range(self, default_cyl: Cylinder) -> None:
        with pytest.raises(ValueError, match="between the floor and the axis"):
            isothermal_profile(default_cyl, AIR_290, [-5.0])


class TestLapseRate:
    def test_floor_value(self, default_cyl: Cylinder) -> None:
        rate = dry_lapse_rate(default_cyl, AIR_290, 0.0)
        assert rate * 1000.0 == pytest.approx(9.76, abs=0.02)

    def test_zero_at_axis(self, default_cyl: Cylinder) -> None:
        assert dry_lapse_rate(default_cyl, AIR_290, 3200.0) == 0.0

    def test_matches_finite_difference_of_profile(self, default_cyl: Cylinder) -> None:
        h, d = 1500.0, 1.0
        prof = adiabatic_profile(default_cyl, AIR_290, [h - d, h + d])
        slope = -(prof.temperature[1] - prof.temperature[0]) / (2.0 * d)
        assert slope == pytest.approx(dry_lapse_rate(default_cyl, AIR_290, h), rel=1e-6)

    def test_height_out_of_range(self, default_cyl: Cylinder) -> None:
        with pytest.raises(ValueError, match="between the floor and the axis"):
            dry_lapse_rate(default_cyl, AIR_290, -1.0)


class TestMass:
    @pytest.mark.parametrize("gas", [AIR_290, SP_290])
    def test_matches_numerical_integral_isothermal(
        self, cyl: Cylinder, gas: Gas
    ) -> None:
        h = np.linspace(0.0, cyl.radius, 20001)
        prof = isothermal_profile(cyl, gas, h)
        r = cyl.radius - h
        numeric = cyl.length * np.trapezoid(prof.density * 2.0 * math.pi * r, h)
        formula = total_mass(cyl, gas.p_rim, isothermal_axis_ratio(cyl, gas))
        assert formula == pytest.approx(numeric, rel=1e-8)

    def test_matches_numerical_integral_adiabatic(self, cyl: Cylinder) -> None:
        h = np.linspace(0.0, cyl.radius, 20001)
        prof = adiabatic_profile(cyl, AIR_290, h)
        r = cyl.radius - h
        numeric = cyl.length * np.trapezoid(prof.density * 2.0 * math.pi * r, h)
        formula = total_mass(cyl, AIR_290.p_rim, adiabatic_axis_ratio(cyl, AIR_290))
        assert formula == pytest.approx(numeric, rel=1e-8)

    def test_default_habitat_mass_order_of_magnitude(
        self, default_cyl: Cylinder
    ) -> None:
        mass = total_mass(
            default_cyl, 101325.0, isothermal_axis_ratio(default_cyl, AIR_290)
        )
        volume = math.pi * 3200.0**2 * 20000.0
        assert mass / volume == pytest.approx(1.1, rel=0.1)  # kg/m^3 mean density

    @given(st.floats(1e3, 2e5), st.floats(0.05, 0.99))
    def test_round_trip(self, p_rim: float, ratio: float) -> None:
        cyl = Cylinder()
        mass = total_mass(cyl, p_rim, ratio)
        assert rim_pressure_from_mass(cyl, mass, ratio) == pytest.approx(p_rim)

    @pytest.mark.parametrize("ratio", [0.0, -0.2, 1.1])
    def test_bad_ratio(self, default_cyl: Cylinder, ratio: float) -> None:
        with pytest.raises(ValueError, match="axis_ratio"):
            total_mass(default_cyl, 1e5, ratio)
        with pytest.raises(ValueError, match="axis_ratio"):
            rim_pressure_from_mass(default_cyl, 1e9, ratio)

    def test_no_spin_limit_is_rejected(self, default_cyl: Cylinder) -> None:
        with pytest.raises(ValueError, match="no mass information"):
            rim_pressure_from_mass(default_cyl, 1e9, 1.0)

    def test_non_positive_inputs(self, default_cyl: Cylinder) -> None:
        with pytest.raises(ValueError, match="p_rim"):
            total_mass(default_cyl, 0.0, 0.8)
        with pytest.raises(ValueError, match="mass"):
            rim_pressure_from_mass(default_cyl, -1.0, 0.8)


def magnus_dew_point(t_kelvin: float, rh: float) -> float:
    """Independent dew point (K) from the Magnus formula."""
    tc = t_kelvin - 273.15
    gamma = math.log(rh) + 17.67 * tc / (243.5 + tc)
    return 273.15 + 243.5 * gamma / (17.67 - gamma)


def bolton_lcl_temperature(t_kelvin: float, rh: float) -> float:
    """Independent LCL temperature: Magnus dew point + Bolton eq. 15."""
    t_dew = magnus_dew_point(t_kelvin, rh)
    return 1.0 / (1.0 / (t_dew - 56.0) + math.log(t_kelvin / t_dew) / 800.0) + 56.0


class TestCloudBase:
    @pytest.mark.parametrize("rh", [0.3, 0.5, 0.7, 0.9])
    def test_matches_independent_lcl_temperature(
        self, default_cyl: Cylinder, rh: float
    ) -> None:
        cloud = lifting_condensation_level(default_cyl, AIR_290, rh)
        assert cloud.temperature == pytest.approx(
            bolton_lcl_temperature(290.0, rh), abs=0.15
        )

    def test_fifty_percent_in_o_neill_band(self, default_cyl: Cylinder) -> None:
        cloud = lifting_condensation_level(default_cyl, AIR_290, 0.5)
        assert cloud.earth_equivalent_height == pytest.approx(1300.0, abs=100.0)
        assert cloud.height is not None
        assert 1400.0 < cloud.height < 2100.0

    @pytest.mark.parametrize("rh", [0.5, 0.7, 0.9])
    def test_rule_of_thumb_125_m_per_kelvin_of_dew_point_spread(
        self, default_cyl: Cylinder, rh: float
    ) -> None:
        cloud = lifting_condensation_level(default_cyl, AIR_290, rh)
        spread = 290.0 - magnus_dew_point(290.0, rh)
        assert cloud.earth_equivalent_height == pytest.approx(125.0 * spread, rel=0.1)

    def test_saturated_air_condenses_at_the_floor(self, default_cyl: Cylinder) -> None:
        cloud = lifting_condensation_level(default_cyl, AIR_290, 1.0)
        assert cloud.phi == 0.0
        assert cloud.height == 0.0
        assert cloud.temperature == pytest.approx(290.0)

    def test_dry_air_never_condenses_before_the_axis(
        self, default_cyl: Cylinder
    ) -> None:
        cloud = lifting_condensation_level(default_cyl, AIR_290, 0.05)
        assert cloud.height is None
        assert cloud.phi > axis_geopotential(default_cyl)
        assert cloud.earth_equivalent_height == pytest.approx(cloud.phi / G)

    @given(st.floats(0.05, 0.99))
    def test_wetter_air_has_a_lower_cloud_base(self, rh: float) -> None:
        cyl = Cylinder()
        wetter = lifting_condensation_level(cyl, AIR_290, min(rh + 0.01, 1.0))
        drier = lifting_condensation_level(cyl, AIR_290, rh)
        assert wetter.phi < drier.phi

    @given(st.floats(0.5, 1.0))
    def test_true_height_exceeds_earth_equivalent(self, rh: float) -> None:
        cloud = lifting_condensation_level(Cylinder(), AIR_290, rh)
        assert cloud.height is not None
        assert cloud.height >= cloud.earth_equivalent_height - 1e-9

    def test_height_is_phi_converted(self, default_cyl: Cylinder) -> None:
        cloud = lifting_condensation_level(default_cyl, SP_290, 0.5)
        assert cloud.height is not None
        assert geopotential(default_cyl, cloud.height) == pytest.approx(cloud.phi)

    @pytest.mark.parametrize("rh", [0.0, -0.1, 1.2])
    def test_bad_humidity(self, default_cyl: Cylinder, rh: float) -> None:
        with pytest.raises(ValueError, match="relative_humidity"):
            lifting_condensation_level(default_cyl, AIR_290, rh)
