import math
from dataclasses import FrozenInstanceError

import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from oneill.atmo.surface import (
    FOREST,
    KAPPA,
    ORCHARD,
    SUBURBAN,
    Stability,
    Terrain,
    bradshaw_number,
    classify_rotational_stability,
    dynamic_pressure,
    friction_velocity,
    gust_speed,
    log_law_wind,
    rotational_stability_profile,
    spinward_fraction,
    unstable_layer,
    wind_shear,
    wind_speed,
)
from oneill.core.cylinder import Cylinder

OMEGA = Cylinder().omega  # 1 g, R = 3.2 km: 0.05536 1/s
TERRAINS = [SUBURBAN, ORCHARD, FOREST]
terrains = st.sampled_from(TERRAINS)


class TestTerrain:
    @pytest.mark.parametrize("terrain", TERRAINS)
    def test_presets_sit_inside_research_ranges(self, terrain: Terrain) -> None:
        assert 0.5 <= terrain.d / terrain.h_c <= 0.8
        assert 0.05 <= terrain.z0 / terrain.h_c <= 0.13
        assert terrain.attenuation > 0.0

    def test_suburb_matches_research_table(self) -> None:
        assert 0.3 <= SUBURBAN.z0 <= 0.8
        assert 3.0 <= SUBURBAN.d <= 7.0
        assert SUBURBAN.h_c == 10.0

    def test_forest_and_orchard_roughness(self) -> None:
        assert 1.0 <= FOREST.z0 <= 2.0
        assert 0.5 <= ORCHARD.z0 <= 1.0

    def test_attenuation_values(self) -> None:
        # hand values: 10 / (5 ln 10), 6 / (2 ln(10/3)), 15 / (5 ln(10/3))
        assert SUBURBAN.attenuation == pytest.approx(0.8686, abs=1e-3)
        assert ORCHARD.attenuation == pytest.approx(2.4919, abs=1e-3)
        assert FOREST.attenuation == pytest.approx(2.4919, abs=1e-3)

    def test_frozen(self) -> None:
        with pytest.raises(FrozenInstanceError):
            SUBURBAN.z0 = 1.0  # type: ignore[misc]

    @pytest.mark.parametrize(
        ("z0", "d", "h_c"),
        [(0.0, 5.0, 10.0), (0.5, 0.0, 10.0), (0.5, 5.0, -1.0), (5.0, 5.0, 10.0)],
    )
    def test_invalid(self, z0: float, d: float, h_c: float) -> None:
        with pytest.raises(ValueError, match=r"must be|canopy top"):
            Terrain("bad", z0, d, h_c)


class TestLogLaw:
    def test_zero_at_d_plus_z0(self) -> None:
        u = log_law_wind(SUBURBAN.d + SUBURBAN.z0, 0.4, SUBURBAN.z0, SUBURBAN.d)
        assert float(u) == pytest.approx(0.0, abs=1e-12)

    def test_negative_between_d_and_d_plus_z0(self) -> None:
        assert float(log_law_wind(5.2, 0.4, 0.5, 5.0)) < 0.0

    def test_one_e_fold_per_u_star_over_kappa(self) -> None:
        z = 5.0 + 0.5 * math.e
        assert float(log_law_wind(z, 0.4, 0.5, 5.0)) == pytest.approx(0.4 / KAPPA)

    @pytest.mark.parametrize("z", [5.0, 3.0])
    def test_needs_height_above_d(self, z: float) -> None:
        with pytest.raises(ValueError, match="above the displacement"):
            log_law_wind(z, 0.4, 0.5, 5.0)

    def test_invalid_z0(self) -> None:
        with pytest.raises(ValueError, match="z0"):
            log_law_wind(10.0, 0.4, 0.0, 5.0)


class TestCompositeProfile:
    @given(terrains, st.floats(0.05, 1.5))
    def test_continuous_speed_and_shear_at_canopy_top(
        self, terrain: Terrain, ustar: float
    ) -> None:
        eps = 1e-9
        lo = terrain.h_c - eps
        hi = terrain.h_c + eps
        assert float(wind_speed(lo, ustar, terrain)) == pytest.approx(
            float(wind_speed(hi, ustar, terrain)), rel=1e-7
        )
        assert float(wind_shear(lo, ustar, terrain)) == pytest.approx(
            float(wind_shear(hi, ustar, terrain)), rel=1e-7
        )

    @given(terrains, st.floats(0.05, 1.5))
    def test_log_law_above_canopy(self, terrain: Terrain, ustar: float) -> None:
        z = np.linspace(terrain.h_c, 10.0 * terrain.h_c, 30)
        assert wind_speed(z, ustar, terrain) == pytest.approx(
            log_law_wind(z, ustar, terrain.z0, terrain.d)
        )

    @given(terrains)
    def test_exponential_decay_in_canopy(self, terrain: Terrain) -> None:
        z = np.linspace(0.1, terrain.h_c, 40)
        u = wind_speed(z, 0.4, terrain)
        assert np.all(np.diff(u) > 0.0)
        assert np.all(u > 0.0)
        ratio = u[-1] / u[0]
        expected = math.exp(terrain.attenuation * (1.0 - 0.1 / terrain.h_c))
        assert ratio == pytest.approx(expected)

    @given(terrains, st.floats(0.05, 1.5), st.floats(0.5, 60.0))
    def test_analytic_shear_matches_finite_difference(
        self, terrain: Terrain, ustar: float, z: float
    ) -> None:
        d = 1e-6
        if abs(z - terrain.h_c) < 10 * d:
            return
        hi = wind_speed(z + d, ustar, terrain)
        lo = wind_speed(z - d, ustar, terrain)
        assert float(wind_shear(z, ustar, terrain)) == pytest.approx(
            float((hi - lo) / (2 * d)), rel=1e-5
        )

    @pytest.mark.parametrize("z", [0.0, -1.0])
    def test_heights_must_be_positive(self, z: float) -> None:
        with pytest.raises(ValueError, match="positive"):
            wind_speed(z, 0.3, SUBURBAN)
        with pytest.raises(ValueError, match="positive"):
            wind_shear(z, 0.3, SUBURBAN)

    @given(terrains, st.floats(0.01, 2.0))
    def test_speed_scales_with_ustar(self, terrain: Terrain, ustar: float) -> None:
        z = np.array([2.0, terrain.h_c, 3.0 * terrain.h_c])
        assert wind_speed(z, ustar, terrain) == pytest.approx(
            ustar * wind_speed(z, 1.0, terrain)
        )


class TestFrictionVelocity:
    @given(terrains, st.floats(0.1, 40.0), st.floats(1.0, 20.0))
    def test_round_trip(self, terrain: Terrain, speed: float, above: float) -> None:
        height = terrain.h_c + above
        ustar = friction_velocity(speed, height, terrain)
        assert float(wind_speed(height, ustar, terrain)) == pytest.approx(speed)

    def test_round_trip_inside_canopy(self) -> None:
        ustar = friction_velocity(0.5, 4.0, SUBURBAN)
        assert float(wind_speed(4.0, ustar, SUBURBAN)) == pytest.approx(0.5)

    def test_suburb_ratio(self) -> None:
        ratio = friction_velocity(1.0, SUBURBAN.d + 12.5, SUBURBAN)
        assert ratio == pytest.approx(0.12, abs=0.01)

    def test_forest_ratio(self) -> None:
        assert friction_velocity(1.0, FOREST.d + 15.0, FOREST) == pytest.approx(
            0.17, abs=0.01
        )

    def test_open_grass_ratio(self) -> None:
        # 10 m over z0 = 0.03 m grass: u*/U = kappa / ln(10 / 0.03) ~ 0.07
        speed = float(log_law_wind(10.0, 1.0, 0.03, 1e-12))
        assert 1.0 / speed == pytest.approx(0.07, abs=0.005)

    @pytest.mark.parametrize(("speed", "height"), [(0.0, 20.0), (3.0, 0.0)])
    def test_invalid(self, speed: float, height: float) -> None:
        with pytest.raises(ValueError, match="must be"):
            friction_velocity(speed, height, SUBURBAN)


class TestGustsAndLoads:
    def test_suburban_gust_factor_in_research_range(self) -> None:
        mean = 3.0
        ustar = friction_velocity(mean, SUBURBAN.d + 12.5, SUBURBAN)
        factor = gust_speed(mean, ustar) / mean
        assert 1.7 <= factor <= 2.0

    def test_open_country_gust_factor(self) -> None:
        assert 1.4 <= gust_speed(1.0, 0.07) <= 1.5

    def test_custom_peak_factor(self) -> None:
        assert gust_speed(2.0, 0.1, peak_factor=0.0) == 2.0
        assert gust_speed(2.0, 0.1, peak_factor=2.0, sigma_ratio=2.0) == (
            pytest.approx(2.4)
        )

    def test_dynamic_pressure(self) -> None:
        assert dynamic_pressure(1.2, 10.0) == pytest.approx(60.0)

    def test_half_atmosphere_halves_the_load(self) -> None:
        assert dynamic_pressure(0.6, 10.0) == pytest.approx(
            0.5 * dynamic_pressure(1.2, 10.0)
        )

    def test_invalid(self) -> None:
        with pytest.raises(ValueError, match="density"):
            dynamic_pressure(0.0, 1.0)
        with pytest.raises(ValueError, match="ustar"):
            gust_speed(1.0, 0.0)


class TestRotationalStability:
    def test_spinward_fraction(self) -> None:
        assert spinward_fraction(90.0) == 1.0
        assert spinward_fraction(270.0) == -1.0
        assert spinward_fraction(0.0) == 0.0
        assert spinward_fraction(180.0) == 0.0
        assert spinward_fraction(30.0) == pytest.approx(0.5)

    def test_bradshaw_number_values(self) -> None:
        # S = -1/2 (U' = 4 omega) is the most destabilised case: B = -1/4.
        assert bradshaw_number(4.0 * OMEGA, OMEGA) == pytest.approx(-0.25)
        assert bradshaw_number(2.0 * OMEGA, OMEGA) == 0.0
        # antispinward U' = -omega: S = 2, B = 2 * 3
        assert bradshaw_number(-OMEGA, OMEGA) == pytest.approx(6.0)

    def test_bradshaw_undefined_without_shear(self) -> None:
        with pytest.raises(ValueError, match="undefined"):
            bradshaw_number(0.0, OMEGA)

    @pytest.mark.parametrize(
        ("factor", "expected"),
        [
            (0.5, Stability.STABLE),
            (1.0 - 1e-9, Stability.STABLE),
            (1.0, Stability.NEUTRAL),
            (1.0 + 1e-9, Stability.UNSTABLE),
            (2.0, Stability.UNSTABLE),
            (4.0, Stability.UNSTABLE),
            (100.0, Stability.UNSTABLE),
            (0.0, Stability.NEUTRAL),
            (-0.5, Stability.STABLE),
            (-3.0, Stability.STABLE),
        ],
    )
    def test_classify_by_shear_over_two_omega(
        self, factor: float, expected: Stability
    ) -> None:
        assert classify_rotational_stability(factor * 2.0 * OMEGA, OMEGA) is expected

    @given(
        st.floats(-5.0, 5.0).filter(lambda f: f == 0.0 or abs(f) > 1e-6),
        st.floats(1e-3, 1.0),
    )
    def test_classification_matches_rayleigh(self, factor: float, omega: float) -> None:
        # Rayleigh: angular momentum falls outward (unstable) iff U' > 2 omega.
        result = classify_rotational_stability(factor * 2.0 * omega, omega)
        if factor > 1.0 + 1e-9:
            assert result is Stability.UNSTABLE
        elif factor < 1.0 - 1e-9 and factor != 0.0:
            assert result is Stability.STABLE

    @pytest.mark.parametrize("terrain", TERRAINS)
    def test_spinward_unstable_only_in_the_bottom_layer(self, terrain: Terrain) -> None:
        z = np.linspace(0.05, 200.0, 4000)
        prof = rotational_stability_profile(z, 0.35, terrain, OMEGA, 1.0)
        unstable = prof.classification == Stability.UNSTABLE
        assert np.array_equal(unstable, prof.shear > 2.0 * OMEGA)
        assert np.array_equal(unstable, prof.bradshaw < 0.0)
        idx = np.flatnonzero(unstable)
        assert idx.size > 0
        assert np.all(np.diff(idx) == 1)  # one contiguous layer
        assert np.all(prof.classification[idx[-1] + 1 :] == Stability.STABLE)
        assert z[idx[-1]] < 15.0 + terrain.d

    @pytest.mark.parametrize("terrain", TERRAINS)
    def test_analytic_layer_matches_grid(self, terrain: Terrain) -> None:
        layer = unstable_layer(0.35, terrain, OMEGA)
        assert layer is not None
        low, high = layer
        z = np.linspace(0.01, 100.0, 100000)
        prof = rotational_stability_profile(z, 0.35, terrain, OMEGA, 1.0)
        flagged = z[prof.classification == Stability.UNSTABLE]
        assert flagged.min() == pytest.approx(max(low, 0.01), abs=2e-3)
        assert flagged.max() == pytest.approx(high, abs=2e-3)

    def test_antispinward_is_always_stable(self) -> None:
        z = np.linspace(0.05, 200.0, 500)
        prof = rotational_stability_profile(z, 0.5, SUBURBAN, OMEGA, -1.0)
        assert np.all(prof.classification == Stability.STABLE)
        assert np.all(prof.bradshaw > 0.0)

    def test_axial_is_neutral(self) -> None:
        z = np.linspace(0.05, 200.0, 500)
        axial = spinward_fraction(180.0)
        prof = rotational_stability_profile(z, 0.5, SUBURBAN, OMEGA, axial)
        assert np.all(prof.classification == Stability.NEUTRAL)
        assert np.all(np.isnan(prof.bradshaw))

    def test_oblique_wind_uses_only_the_spinward_component(self) -> None:
        z = np.array([8.0])
        full = rotational_stability_profile(z, 0.35, SUBURBAN, OMEGA, 1.0)
        half = rotational_stability_profile(z, 0.35, SUBURBAN, OMEGA, 0.5)
        assert half.shear == pytest.approx(0.5 * full.shear)

    @pytest.mark.parametrize("fraction", [-1.1, 1.5])
    def test_fraction_range(self, fraction: float) -> None:
        with pytest.raises(ValueError, match="spinward"):
            rotational_stability_profile([5.0], 0.3, SUBURBAN, OMEGA, fraction)


class TestRayleighCriterion:
    """Pin the classifier to the physical frame via angular momentum."""

    @staticmethod
    def angular_momentum(r: float, shear: float) -> float:
        # absolute angular momentum per unit mass of air at radius r = R - z, with
        # spinward relative wind U(z) = shear * z
        z = CYL_R - r
        return r * (OMEGA * r + shear * z)

    @given(
        st.floats(-6.0, 6.0).filter(lambda f: abs(f - 1.0) > 0.03),
        st.floats(1.0, 20.0),
    )
    def test_unstable_iff_angular_momentum_decreases_outward(
        self, factor: float, z: float
    ) -> None:
        shear = factor * 2.0 * OMEGA
        r, dr = CYL_R - z, 1e-3
        dl_dr = (
            self.angular_momentum(r + dr, shear) - self.angular_momentum(r - dr, shear)
        ) / (2.0 * dr)
        unstable = classify_rotational_stability(shear, OMEGA) is Stability.UNSTABLE
        assert unstable == (dl_dr < 0.0)

    def test_antispinward_and_axial_have_no_decrease(self) -> None:
        assert classify_rotational_stability(-0.3, OMEGA) is Stability.STABLE
        assert classify_rotational_stability(0.0, OMEGA) is Stability.NEUTRAL


CYL_R = Cylinder().radius


class TestUnstableLayer:
    def test_suburb_three_metres_per_second(self) -> None:
        ustar = friction_velocity(3.0, 20.0, SUBURBAN)
        layer = unstable_layer(ustar, SUBURBAN, OMEGA)
        assert layer is not None
        low, high = layer
        assert high == pytest.approx(SUBURBAN.d + ustar / (2.0 * OMEGA * KAPPA))
        assert high > SUBURBAN.h_c
        assert 0.0 < low < SUBURBAN.h_c

    def test_weak_wind_has_no_unstable_layer(self) -> None:
        assert unstable_layer(0.01, SUBURBAN, OMEGA) is None

    def test_strong_wind_is_unstable_to_the_floor(self) -> None:
        layer = unstable_layer(2.0, SUBURBAN, OMEGA)
        assert layer is not None
        assert layer[0] == 0.0

    @given(terrains, st.floats(0.02, 2.0))
    def test_layer_grows_with_wind(self, terrain: Terrain, ustar: float) -> None:
        weak = unstable_layer(ustar, terrain, OMEGA)
        strong = unstable_layer(1.5 * ustar, terrain, OMEGA)
        if weak is not None:
            assert strong is not None
            assert strong[1] >= weak[1]
            assert strong[0] <= weak[0]

    def test_invalid(self) -> None:
        with pytest.raises(ValueError, match="ustar"):
            unstable_layer(0.0, SUBURBAN, OMEGA)
        with pytest.raises(ValueError, match="omega"):
            unstable_layer(0.3, SUBURBAN, 0.0)
        with pytest.raises(ValueError, match="omega"):
            classify_rotational_stability(1.0, 0.0)
