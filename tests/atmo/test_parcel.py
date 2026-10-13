import math
from dataclasses import replace

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from oneill.atmo.equilibrium import (
    adiabatic_profile,
    isothermal_profile,
    lifting_condensation_level,
)
from oneill.atmo.gas import ONEILL_1974
from oneill.atmo.numbers import brunt_vaisala
from oneill.atmo.parcel import (
    ParcelResult,
    Release,
    RunSettings,
    simulate_parcel,
)
from oneill.core.cylinder import Cylinder
from oneill.core.vectors import LocalVector

CYL = Cylinder()
OMEGA = CYL.omega
GAS = ONEILL_1974.at_temperature(290.0)
FAST = RunSettings(t_end=120.0, sample_dt=0.5, max_step=0.5)


def run(release: Release | None = None, **kwargs: object) -> ParcelResult:
    return simulate_parcel(CYL, GAS, release, **kwargs)  # type: ignore[arg-type]


def crossings(t: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Linearly interpolated zero-crossing times of ``x``."""
    i = np.flatnonzero(x[:-1] * x[1:] < 0.0)
    return np.asarray(t[i] - x[i] * (t[i + 1] - t[i]) / (x[i + 1] - x[i]))


class TestNeutralParcel:
    @pytest.mark.parametrize("drag", [0.0, 0.05])
    @pytest.mark.parametrize("height", [0.0, 1500.0])
    def test_at_rest_stays_at_rest(self, height: float, drag: float) -> None:
        result = run(Release(height=height, drag=drag), settings=FAST)
        assert np.max(np.abs(result.velocity)) < 1e-8
        assert np.max(np.abs(result.height - height)) < 1e-6

    @settings(max_examples=8)
    @given(st.floats(0.0, 3000.0), st.floats(0.0, 0.2))
    def test_at_rest_property(self, height: float, drag: float) -> None:
        result = run(Release(height=height, drag=drag), settings=FAST)
        assert result.peak_speed < 1e-8

    def test_default_release_is_neutral_at_rest(self) -> None:
        result = simulate_parcel(CYL, GAS)
        assert result.peak_speed < 1e-8
        assert result.max_rise < 1e-6

    @pytest.mark.parametrize(
        "kick", [LocalVector(0.0, 0.0, 1.0), LocalVector(0.5, 0.3, -0.7)]
    )
    def test_kick_traces_inertial_circle(self, kick: LocalVector) -> None:
        period = math.pi / OMEGA
        cfg = RunSettings(t_end=period, sample_dt=0.5, max_step=0.05)
        result = run(Release(height=1000.0, velocity=kick), settings=cfg)
        p0 = result.position[0]
        v0 = CYL.local_basis(p0).T @ kick.to_array()
        phi = 2.0 * OMEGA * result.t
        sin, cos = np.sin(phi), np.cos(phi)
        expected = np.stack(
            [
                p0[0] + (v0[0] * sin + v0[1] * (1 - cos)) / (2 * OMEGA),
                p0[1] + (v0[1] * sin - v0[0] * (1 - cos)) / (2 * OMEGA),
                p0[2] + v0[2] * result.t,
            ],
            axis=1,
        )
        assert result.position == pytest.approx(expected, abs=1e-5)
        # circle radius v/(2 omega), period pi/omega: back at the start
        horizontal = math.hypot(v0[0], v0[1])
        centre = np.array([p0[0] + v0[1] / (2 * OMEGA), p0[1] - v0[0] / (2 * OMEGA)])
        dist = np.hypot(
            result.position[:, 0] - centre[0], result.position[:, 1] - centre[1]
        )
        assert dist == pytest.approx(horizontal / (2 * OMEGA), abs=1e-5)
        assert result.position[-1] == pytest.approx(
            p0 + np.array([0.0, 0.0, v0[2] * period]), abs=1e-5
        )
        assert result.speed == pytest.approx(np.linalg.norm(v0), abs=1e-6)

    def test_kicked_neutral_parcel_oscillates_at_two_omega(self) -> None:
        cfg = RunSettings(t_end=400.0, sample_dt=0.25, max_step=0.25)
        kick = LocalVector(0.0, 0.0, 0.1)
        result = run(Release(height=1000.0, velocity=kick), settings=cfg)
        cross = crossings(result.t, result.height - 1000.0)
        measured = math.pi * (len(cross) - 1) / (cross[-1] - cross[0])
        assert measured == pytest.approx(2.0 * OMEGA, rel=1e-4)


class TestStableStratification:
    def test_oscillates_at_sqrt_n2_plus_4_omega2(self) -> None:
        cfg = RunSettings(t_end=400.0, sample_dt=0.25, max_step=0.25)
        result = run(
            Release(height=1000.0, velocity=LocalVector(0.0, 0.0, 0.1)),
            environment=isothermal_profile,
            settings=cfg,
        )
        cross = crossings(result.t, result.height - 1000.0)
        measured = math.pi * (len(cross) - 1) / (cross[-1] - cross[0])
        n = brunt_vaisala(CYL.gravity_at_height(1000.0), 290.0, 0.0, GAS.cp)
        expected = math.sqrt(n**2 + 4.0 * OMEGA**2)
        assert measured == pytest.approx(expected, rel=1e-3)
        assert measured > 2.0 * OMEGA * 1.005  # stratification is resolved


class TestWarmParcel:
    def test_conserves_angular_momentum_and_drifts_spinward(self) -> None:
        cfg = RunSettings(t_end=60.0, sample_dt=0.25, max_step=0.25)
        result = run(Release(temperature_excess=1.0), settings=cfg)
        x, y = result.position[:, 0], result.position[:, 1]
        vx, vy = result.velocity[:, 0], result.velocity[:, 1]
        ang_mom = x * vy - y * vx + OMEGA * (x**2 + y**2)
        assert ang_mom == pytest.approx(OMEGA * CYL.radius**2, rel=1e-9)
        r = CYL.radius - result.height
        tangential = (x * vy - y * vx) / r
        assert tangential == pytest.approx(OMEGA * (CYL.radius**2 / r - r), abs=1e-6)
        assert result.spinward[-1] > 0.0

    def test_drift_speed_is_two_omega_per_metre_of_rise(self) -> None:
        cfg = RunSettings(t_end=100.0, sample_dt=0.05, max_step=0.25, rise=2.0)
        result = run(Release(temperature_excess=1.0), settings=cfg)
        assert result.drift_speed == pytest.approx(
            2.0 * OMEGA * 2.0 * (1 + 2.0 / (2 * CYL.radius)), rel=2e-3
        )
        assert result.deflection is not None
        assert result.deflection > 0.0

    def test_coriolis_traps_the_thermal(self) -> None:
        # constant buoyancy b in a rotating frame gives a cycloid: the parcel only
        # rises b / (2 omega^2) = R dT / (2 T), then drifts spinward at b / (2 omega).
        cfg = RunSettings(t_end=600.0, sample_dt=0.5, max_step=0.25)
        result = run(Release(temperature_excess=1.0), settings=cfg)
        b = CYL.surface_gravity * 1.0 / 290.0
        assert result.max_rise == pytest.approx(b / (2 * OMEGA**2), rel=0.01)
        assert result.spinward[-1] / 600.0 == pytest.approx(b / (2 * OMEGA), rel=0.05)
        assert result.deflection is None
        assert result.drift_speed is None

    def test_drag_shrinks_the_deflection_monotonically(self) -> None:
        results = [
            run(
                Release(temperature_excess=1.0, drag=k),
                settings=RunSettings(
                    t_end=300.0, sample_dt=0.1, max_step=0.25, rise=2.0
                ),
            )
            for k in (0.0, 0.01, 0.05)
        ]
        deflections = [r.deflection for r in results]
        drifts = [r.drift_speed for r in results]
        peaks = [r.peak_speed for r in results]
        assert None not in deflections
        assert None not in drifts
        assert sorted(deflections, reverse=True) == deflections  # type: ignore[type-var]
        assert len(set(deflections)) == 3
        assert sorted(drifts, reverse=True) == drifts  # type: ignore[type-var]
        assert len(set(drifts)) == 3
        assert sorted(peaks, reverse=True) == peaks
        assert len(set(peaks)) == 3

    def test_strong_drag_lets_the_parcel_rise(self) -> None:
        cfg = RunSettings(t_end=900.0, sample_dt=1.0, max_step=0.5)
        result = run(Release(temperature_excess=3.0, drag=0.05), settings=cfg)
        assert result.deflection is not None
        assert result.drift_speed is not None
        assert 0.0 < result.drift_speed < 2.0 * OMEGA * 100.0


class TestOutputs:
    def test_stop_height(self) -> None:
        cfg = RunSettings(t_end=900.0, sample_dt=1.0, max_step=0.5, stop_height=20.0)
        result = run(Release(temperature_excess=3.0, drag=0.05), settings=cfg)
        assert result.height[-1] >= 20.0 > result.height[-2]
        assert result.t[-1] < 900.0

    def test_hitting_the_floor_ends_the_run(self) -> None:
        result = run(Release(velocity=LocalVector(0.0, 0.0, -2.0)), settings=FAST)
        assert result.ended_at_floor
        assert len(result.t) == 1

    def test_does_not_flag_floor_for_normal_run(self) -> None:
        assert not run(settings=FAST).ended_at_floor

    def test_axial_displacement(self) -> None:
        result = run(
            Release(velocity=LocalVector(0.0, 0.4, 0.0), height=500.0), settings=FAST
        )
        assert result.axial[-1] == pytest.approx(0.4 * result.t[-1], rel=1e-9)

    def test_time_to_cloud_base(self) -> None:
        rh = 0.99
        lcl = lifting_condensation_level(CYL, GAS, rh).height
        assert lcl is not None
        cfg = RunSettings(t_end=900.0, sample_dt=1.0, max_step=0.5)
        result = run(
            Release(temperature_excess=3.0, drag=0.05),
            settings=cfg,
            relative_humidity=rh,
        )
        assert result.lcl_height == lcl
        assert result.time_to_lcl is not None
        j = int(np.searchsorted(result.t, result.time_to_lcl))
        assert result.height[j - 1] < lcl <= result.height[j]

    def test_saturated_air_reaches_cloud_base_at_once(self) -> None:
        result = run(settings=FAST, relative_humidity=1.0)
        assert result.lcl_height == 0.0
        assert result.time_to_lcl == 0.0

    def test_cloud_base_not_reached(self) -> None:
        result = run(
            Release(temperature_excess=1.0), settings=FAST, relative_humidity=0.5
        )
        assert result.lcl_height is not None
        assert result.time_to_lcl is None

    def test_dry_air_has_no_cloud_base(self) -> None:
        result = run(settings=FAST, relative_humidity=0.05)
        assert result.lcl_height is None
        assert result.time_to_lcl is None

    def test_adiabatic_default_matches_explicit_profile(self) -> None:
        explicit = run(
            Release(temperature_excess=1.0),
            environment=adiabatic_profile,
            settings=FAST,
        )
        default = run(Release(temperature_excess=1.0), settings=FAST)
        assert default.height == pytest.approx(explicit.height)


class TestValidation:
    def test_rejects_spin_up(self) -> None:
        with pytest.raises(ValueError, match="constant spin"):
            simulate_parcel(Cylinder.from_gravity(spin_accel=1e-4), GAS)

    @pytest.mark.parametrize("drag", [-0.1, math.nan, math.inf])
    def test_rejects_bad_drag(self, drag: float) -> None:
        with pytest.raises(ValueError, match="drag"):
            Release(drag=drag)

    def test_rejects_non_finite_excess(self) -> None:
        with pytest.raises(ValueError, match="temperature_excess"):
            Release(temperature_excess=math.nan)

    @pytest.mark.parametrize("field", ["t_end", "sample_dt", "max_step", "rise"])
    def test_rejects_bad_settings(self, field: str) -> None:
        with pytest.raises(ValueError, match=field):
            simulate_parcel(CYL, GAS, settings=replace(RunSettings(), **{field: 0.0}))
