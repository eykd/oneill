import math
from typing import TYPE_CHECKING

import pytest
from hypothesis import given
from hypothesis import strategies as st

from oneill.atmo.equilibrium import adiabatic_profile, dry_lapse_rate
from oneill.atmo.gas import ONEILL_1974, SP413
from oneill.atmo.numbers import (
    brunt_vaisala,
    brunt_vaisala_squared,
    dimensionless_table,
    ekman_depth,
    ekman_number,
    inertial_period,
    rossby_deformation_radius,
    rossby_number,
    shedding_frequency,
    spin_up_time,
    stewartson_quarter_thickness,
    stewartson_third_thickness,
    viscous_scales,
)

if TYPE_CHECKING:
    from collections.abc import Callable

    from oneill.core.cylinder import Cylinder

OMEGA = 0.0554
HALF = 10_000.0
HOUR = 3600.0


class TestElementary:
    def test_rossby(self) -> None:
        assert rossby_number(3.0, 1000.0, OMEGA) == pytest.approx(0.02708, rel=1e-3)

    @given(st.floats(0.01, 100.0), st.floats(1.0, 1e5), st.floats(1e-3, 1.0))
    def test_rossby_is_ratio_of_inertia_to_coriolis(
        self, u: float, length: float, omega: float
    ) -> None:
        assert rossby_number(u, length, omega) == pytest.approx(
            (u**2 / length) / (2.0 * omega * u)
        )

    def test_inertial_period(self) -> None:
        assert inertial_period(OMEGA) == pytest.approx(56.7, abs=0.05)

    @pytest.mark.parametrize(
        ("nu", "depth"), [(3.0e-5, 0.023), (1.0, 4.2), (10.0, 13.4)]
    )
    def test_ekman_depth_matches_research_table(self, nu: float, depth: float) -> None:
        assert ekman_depth(nu, OMEGA) == pytest.approx(depth, abs=0.05, rel=0.02)

    def test_ekman_depth_is_sqrt_two_nu_over_f(self) -> None:
        assert ekman_depth(1.0, OMEGA) == pytest.approx(math.sqrt(2.0 / (2.0 * OMEGA)))

    @pytest.mark.parametrize(("nu", "hours"), [(0.1, 37.0), (1.0, 11.8), (10.0, 3.7)])
    def test_spin_up_time_matches_research_table(self, nu: float, hours: float) -> None:
        assert spin_up_time(nu, OMEGA, HALF) / HOUR == pytest.approx(hours, rel=0.01)

    def test_molecular_spin_up_is_months(self) -> None:
        days = spin_up_time(3.0e-5, OMEGA, HALF) / 86400.0
        assert days == pytest.approx(90.0, rel=0.01)

    def test_stewartson_layers(self) -> None:
        assert stewartson_quarter_thickness(3.0e-5, OMEGA, HALF) == pytest.approx(
            15.0, rel=0.03
        )
        assert stewartson_quarter_thickness(1.0, OMEGA, HALF) == pytest.approx(
            200.0, rel=0.05
        )
        e = ekman_number(1.0, OMEGA, HALF)
        assert e == pytest.approx(1.0 / (OMEGA * HALF**2))
        assert stewartson_third_thickness(1.0, OMEGA, HALF) == pytest.approx(
            HALF * e ** (1.0 / 3.0)
        )
        assert stewartson_third_thickness(1.0, OMEGA, HALF) < (
            stewartson_quarter_thickness(1.0, OMEGA, HALF)
        )  # E < 1, so E^(1/3) < E^(1/4)

    @given(st.floats(1e-6, 100.0), st.floats(1e-3, 1.0), st.floats(10.0, 1e5))
    def test_spin_up_is_ekman_time_over_sqrt_e(
        self, nu: float, omega: float, half: float
    ) -> None:
        # tau = H / sqrt(nu omega) = (1 / omega) E^(-1/2)
        e = ekman_number(nu, omega, half)
        assert spin_up_time(nu, omega, half) == pytest.approx(e**-0.5 / omega)

    def test_shedding(self) -> None:
        f = shedding_frequency(3.0, 10.0)
        assert f == pytest.approx(0.06)
        assert 2.0 * math.pi * f / (2.0 * OMEGA) == pytest.approx(3.4, abs=0.05)
        assert shedding_frequency(3.0, 10.0, strouhal=0.13) == pytest.approx(0.039)

    def test_rossby_radius(self) -> None:
        assert rossby_deformation_radius(0.018, 3200.0, OMEGA) == pytest.approx(
            0.018 * 3200.0 / (2.0 * OMEGA)
        )

    @pytest.mark.parametrize(
        "call",
        [
            lambda: rossby_number(1.0, 0.0, OMEGA),
            lambda: ekman_number(0.0, OMEGA, HALF),
            lambda: ekman_depth(1.0, 0.0),
            lambda: spin_up_time(1.0, OMEGA, 0.0),
            lambda: shedding_frequency(1.0, -1.0),
            lambda: rossby_deformation_radius(0.01, 0.0, OMEGA),
            lambda: brunt_vaisala_squared(0.0, 290.0, 0.0, 1005.0),
        ],
    )
    def test_invalid_arguments(self, call: Callable[[], float]) -> None:
        with pytest.raises(ValueError, match="must be"):
            call()


class TestBruntVaisala:
    def test_zero_on_dry_adiabat(self, cyl: Cylinder) -> None:
        for h in (0.0, 0.3 * cyl.radius, 0.9 * cyl.radius):
            g = cyl.gravity_at_height(h)
            dt_dh = -dry_lapse_rate(cyl, ONEILL_1974, h)
            n2 = brunt_vaisala_squared(g, 290.0, dt_dh, ONEILL_1974.cp)
            assert n2 == pytest.approx(0.0, abs=1e-12)

    def test_isothermal_value(self, default_cyl: Cylinder) -> None:
        g = default_cyl.gravity_at_height(0.0)
        n = brunt_vaisala(g, 290.0, 0.0, ONEILL_1974.cp)
        assert n == pytest.approx(g / math.sqrt(ONEILL_1974.cp * 290.0))
        assert n == pytest.approx(0.0182, abs=2e-4)

    def test_profile_finite_difference_matches(self, default_cyl: Cylinder) -> None:
        h, d = 800.0, 1.0
        gas = ONEILL_1974.at_temperature(290.0)
        prof = adiabatic_profile(default_cyl, gas, [h - d, h + d])
        dt_dh = (prof.temperature[1] - prof.temperature[0]) / (2.0 * d)
        n2 = brunt_vaisala_squared(
            default_cyl.gravity_at_height(h), 290.0, dt_dh, gas.cp
        )
        assert abs(n2) < 1e-10

    def test_unstable_layer_raises(self) -> None:
        assert brunt_vaisala_squared(9.8, 290.0, -0.02, 1005.0) < 0.0
        with pytest.raises(ValueError, match="convectively unstable"):
            brunt_vaisala(9.8, 290.0, -0.02, 1005.0)

    def test_n_is_far_below_coriolis(self, default_cyl: Cylinder) -> None:
        n = brunt_vaisala(9.8, 290.0, 0.0, 1005.0)
        assert n < 2.0 * default_cyl.omega


class TestTable:
    def test_default_table(self, default_cyl: Cylinder) -> None:
        gas = SP413.at_temperature(290.0)
        table = dimensionless_table(default_cyl, gas)
        assert table.omega == default_cyl.omega
        assert table.coriolis_parameter == 2.0 * default_cyl.omega
        assert table.inertial_period == pytest.approx(56.7, abs=0.1)
        assert table.rossby == pytest.approx(3.0 / (2.0 * default_cyl.omega * 1000.0))
        assert table.obstacle_rossby == pytest.approx(
            3.0 / (2.0 * default_cyl.omega * 10.0)
        )
        assert table.shedding_frequency == pytest.approx(0.06)
        assert table.shedding_to_coriolis == pytest.approx(3.4, abs=0.05)
        assert table.molecular.nu == gas.nu
        assert [s.nu for s in table.eddy] == [1.0, 10.0]
        assert table.eddy[0].ekman_depth == pytest.approx(4.2, abs=0.1)
        assert table.eddy[1].spin_up_time / HOUR == pytest.approx(3.7, rel=0.02)
        assert table.rossby_radius == pytest.approx(
            table.brunt_vaisala * 3200.0 / (2.0 * default_cyl.omega)
        )

    def test_custom_arguments(self, default_cyl: Cylinder) -> None:
        table = dimensionless_table(
            default_cyl,
            ONEILL_1974,
            wind_speed=5.0,
            length_scale=200.0,
            obstacle_diameter=20.0,
            nu_eddy=[0.5],
            height=1000.0,
            depth=500.0,
        )
        assert len(table.eddy) == 1
        assert table.rossby_radius == pytest.approx(
            table.brunt_vaisala * 500.0 / (2.0 * default_cyl.omega)
        )
        floor = dimensionless_table(default_cyl, ONEILL_1974)
        assert table.brunt_vaisala < floor.brunt_vaisala  # weaker g higher up

    def test_viscous_scales_bundle(self) -> None:
        s = viscous_scales(1.0, OMEGA, HALF)
        assert s.ekman_depth == ekman_depth(1.0, OMEGA)
        assert s.spin_up_time == spin_up_time(1.0, OMEGA, HALF)
