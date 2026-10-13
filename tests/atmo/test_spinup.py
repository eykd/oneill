import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.special import j1, jn_zeros

from oneill.atmo.numbers import spin_up_time
from oneill.atmo.spinup import (
    HullSpin,
    SpinupParams,
    hull_spin_up,
    peak_wind_at_height,
    simulate,
    steady_lag,
    sudden_change,
    time_to_adjust,
    time_to_corotation,
)
from oneill.core.cylinder import Cylinder

CYL = Cylinder()
SMALL = Cylinder(radius=100.0, length=1000.0)  # fast diffusion: tests stay short


def analytic_lag(r: np.ndarray, radius: float, alpha: float, nu: float) -> np.ndarray:
    """Steady lag u = alpha r (r^2 - R^2) / (8 nu), derived from nu L[u] = alpha r."""
    return alpha * r * (r**2 - radius**2) / (8.0 * nu)


class TestSteadyLag:
    @pytest.mark.parametrize("stretch", [0.0, 2.5])
    @pytest.mark.parametrize("nu", [1.0, 10.0])
    def test_matches_closed_form(self, nu: float, stretch: float) -> None:
        params = SpinupParams(nu, endcap=False, n_cells=200, stretch=stretch)
        alpha = 1.0e-8
        lag = steady_lag(CYL, params, alpha)
        expected = analytic_lag(lag.radius, CYL.radius, alpha, nu)
        assert lag.u == pytest.approx(expected, abs=1e-4 * np.abs(expected).max())

    def test_lags_for_spin_up_and_leads_for_spin_down(self) -> None:
        params = SpinupParams(1.0, endcap=False)
        assert np.all(steady_lag(CYL, params, 1e-8).u < 0.0)
        assert np.all(steady_lag(CYL, params, -1e-8).u > 0.0)

    def test_scales_inversely_with_viscosity(self) -> None:
        a = steady_lag(CYL, SpinupParams(1.0, endcap=False), 1e-8)
        b = steady_lag(CYL, SpinupParams(4.0, endcap=False), 1e-8)
        assert b.u == pytest.approx(a.u / 4.0, rel=1e-9)

    def test_end_caps_limit_the_lag(self) -> None:
        with_caps = steady_lag(CYL, SpinupParams(1.0), 1e-8)
        without = steady_lag(CYL, SpinupParams(1.0, endcap=False), 1e-8)
        assert with_caps.peak < 0.2 * without.peak
        # interior balance alpha r = -u / tau_E far from the wall
        tau = spin_up_time(1.0, CYL.omega, CYL.length / 2.0)
        mid = int(np.argmin(np.abs(with_caps.radius - 0.3 * CYL.radius)))
        assert with_caps.u[mid] == pytest.approx(
            -1e-8 * with_caps.radius[mid] * tau, rel=0.02
        )

    def test_time_marching_converges_to_the_closed_form(self) -> None:
        params = SpinupParams(1.0, endcap=False, n_cells=160, stretch=0.0)
        alpha = 1.0e-6
        hull = HullSpin(SMALL.omega, alpha)
        run = simulate(
            SMALL, params, hull, np.zeros(160), 6000.0, dt0=5.0, growth=1.0, dt_max=5.0
        )
        expected = analytic_lag(run.radius, SMALL.radius, alpha, 1.0)
        assert run.u[-1] == pytest.approx(expected, abs=1e-3 * np.abs(expected).max())

    def test_quadratic_floor_drag_torque_balance(self) -> None:
        # wall torque R^2 C_D u^2 must carry the whole source alpha R^4 / 4
        c_d, alpha = 2.0e-3, 1.0e-8
        params = SpinupParams(1.0, floor_drag=c_d, endcap=False)
        lag = steady_lag(CYL, params, alpha)
        assert abs(lag.u[-1]) == pytest.approx(
            CYL.radius * math.sqrt(alpha / (4.0 * c_d)), rel=1e-6
        )
        assert lag.u[-1] < 0.0

    def test_wind_at_height(self) -> None:
        params = SpinupParams(1.0, endcap=False)
        lag = steady_lag(CYL, params, 1e-8)
        assert lag.wind_at_height(0.01) == pytest.approx(0.0, abs=0.05 * lag.peak)
        assert abs(lag.wind_at_height(500.0)) > abs(lag.wind_at_height(2.0))
        drag = steady_lag(CYL, SpinupParams(1.0, floor_drag=1e-3, endcap=False), 1e-8)
        assert drag.wind_at_height(0.01) == pytest.approx(drag.u[-1])
        with pytest.raises(ValueError, match="height"):
            lag.wind_at_height(0.0)


def mode_shape(r: np.ndarray, j: float) -> np.ndarray:
    """No-slip eigenmode J1(j r / R) of the viscous operator on SMALL."""
    return np.asarray(j1(j * r / SMALL.radius))


class TestModes:
    @pytest.mark.parametrize("mode", [1, 2, 3])
    def test_bessel_mode_decay_rates(self, mode: int) -> None:
        nu = 1.0
        j = float(jn_zeros(1, mode)[mode - 1])
        n = 240
        params = SpinupParams(nu, endcap=False, n_cells=n, stretch=0.0)
        rate = nu * (j / SMALL.radius) ** 2
        t_end = 1.0 / rate
        dt = t_end / 400.0
        centres = (np.arange(n) + 0.5) * SMALL.radius / n
        run = simulate(
            SMALL,
            params,
            HullSpin(SMALL.omega),
            mode_shape(centres, j),
            t_end,
            dt0=dt,
            growth=1.0,
            dt_max=dt,
        )
        phi = mode_shape(run.radius, j)
        amplitude = run.u @ (run.moment * phi) / (phi @ (run.moment * phi))
        measured = -math.log(amplitude[-1] / amplitude[0]) / run.t[-1]
        assert measured == pytest.approx(rate, rel=2e-3)
        # the mode keeps its shape
        assert run.u[-1] == pytest.approx(
            amplitude[-1] * phi, abs=1e-3 * np.abs(phi).max()
        )

    def test_solid_body_is_an_exact_steady_state(self) -> None:
        params = SpinupParams(5.0, floor_drag=0.0, endcap=False, n_cells=40)
        run = simulate(SMALL, params, HullSpin(SMALL.omega), np.zeros(40), 100.0)
        assert np.max(np.abs(run.u)) == 0.0
        rigid = 0.003 * run.radius
        free = simulate(SMALL, params, HullSpin(SMALL.omega), rigid, 100.0)
        assert free.u[-1] == pytest.approx(rigid, abs=1e-12)


class TestEndCaps:
    def test_interior_decays_at_one_over_tau_e_with_half_length(self) -> None:
        # free-slip floor + solid-body wind: only the end-cap sink acts
        params = SpinupParams(1.0, floor_drag=0.0, n_cells=40)
        tau = spin_up_time(1.0, CYL.omega, CYL.length / 2.0)
        grid = simulate(CYL, params, HullSpin(CYL.omega), np.zeros(40), 1.0)
        u0 = 0.01 * grid.radius
        run = simulate(
            CYL,
            params,
            HullSpin(CYL.omega),
            u0,
            tau,
            dt0=tau / 300,
            growth=1.0,
            dt_max=tau / 300,
        )
        assert run.u[-1] == pytest.approx(u0 * math.exp(-1.0), rel=1e-4)

    def test_longer_habitat_pumps_more_slowly(self) -> None:
        params = SpinupParams(1.0, floor_drag=0.0, n_cells=20)
        probe = Cylinder(radius=500.0, length=2000.0)
        long = Cylinder(radius=500.0, length=4000.0)
        rates = []
        for cyl in (probe, long):
            grid = simulate(cyl, params, HullSpin(cyl.omega), np.zeros(20), 1.0)
            u0 = 0.01 * grid.radius
            run = simulate(
                cyl,
                params,
                HullSpin(cyl.omega),
                u0,
                3000.0,
                dt0=10.0,
                growth=1.0,
                dt_max=10.0,
            )
            rates.append(-math.log(run.u[-1, 5] / u0[5]) / 3000.0)
        assert rates[0] / rates[1] == pytest.approx(2.0, rel=1e-3)

    def test_greenspan_howard_interior_spin_up(self) -> None:
        # small step in hull rate; mid-radius wind decays on tau_E (diffusion from
        # the wall has not reached it yet)
        tau = spin_up_time(1.0, CYL.omega, CYL.length / 2.0)
        run = sudden_change(CYL, SpinupParams(1.0), 1e-4, t_max=tau, dt_max=tau / 200)
        mid = int(np.argmin(np.abs(run.radius - 0.5 * CYL.radius)))
        ratio = run.u[-1, mid] / run.u[0, mid]
        assert ratio == pytest.approx(math.exp(-run.t[-1] / tau), rel=0.05)


class TestBudget:
    @pytest.mark.parametrize("floor", [None, 3.0e-3])
    def test_angular_momentum_budget_closes_exactly(self, floor: float | None) -> None:
        params = SpinupParams(1.0, nu_unstable=8.0, floor_drag=floor, n_cells=60)
        run = hull_spin_up(SMALL, params, SMALL.omega / 200.0, fraction=0.95)
        b = run.budget
        change = np.diff(run.u @ run.moment)
        flux = b.dt * (b.wall + b.endcap + b.hull)
        assert change == pytest.approx(flux, rel=1e-9, abs=1e-9 * np.abs(flux).max())

    def test_total_torque_matches_sudden_change(self) -> None:
        params = SpinupParams(1.0, endcap=False, n_cells=60)
        run = sudden_change(SMALL, params, -0.05 * SMALL.omega)
        b = run.budget
        lag0 = run.relative_angular_momentum[0]
        # the floor removes exactly the initial excess; there is no other sink
        assert np.sum(b.dt * b.wall) == pytest.approx(
            run.relative_angular_momentum[-1] - lag0, rel=1e-9
        )
        assert np.all(b.endcap == 0.0)
        assert np.all(b.hull == 0.0)

    def test_hull_term_equals_apparent_torque(self) -> None:
        params = SpinupParams(1.0, n_cells=30)
        alpha = 1e-4
        run = hull_spin_up(SMALL, params, alpha, fraction=0.5)
        ramp_steps = run.t[1:] <= run.ramp_end + 1e-9
        assert run.budget.hull[ramp_steps] == pytest.approx(
            -alpha * SMALL.radius**4 / 4.0
        )


class TestRegimes:
    def test_spin_down_is_faster_than_spin_up_when_regimes_differ(self) -> None:
        params = SpinupParams(1.0, nu_unstable=10.0, endcap=False, n_cells=100)
        down = time_to_adjust(sudden_change(SMALL, params, -0.5 * SMALL.omega))
        up = time_to_adjust(sudden_change(SMALL, params, +0.5 * SMALL.omega))
        assert down is not None
        assert up is not None
        assert up / down > 2.0

    def test_symmetric_when_regimes_agree(self) -> None:
        params = SpinupParams(3.0, endcap=False, n_cells=100)
        down = time_to_adjust(sudden_change(SMALL, params, -0.5 * SMALL.omega))
        up = time_to_adjust(sudden_change(SMALL, params, +0.5 * SMALL.omega))
        assert down == pytest.approx(up, rel=1e-6)

    @settings(max_examples=6, deadline=None)
    @given(st.floats(2.0, 20.0))
    def test_larger_unstable_viscosity_speeds_spin_down_only(
        self, ratio: float
    ) -> None:
        base = SpinupParams(1.0, endcap=False, n_cells=60)
        fast = SpinupParams(1.0, nu_unstable=ratio, endcap=False, n_cells=60)
        d0 = time_to_adjust(sudden_change(SMALL, base, -0.5 * SMALL.omega))
        d1 = time_to_adjust(sudden_change(SMALL, fast, -0.5 * SMALL.omega))
        u0 = time_to_adjust(sudden_change(SMALL, base, 0.5 * SMALL.omega))
        u1 = time_to_adjust(sudden_change(SMALL, fast, 0.5 * SMALL.omega))
        assert d0 is not None
        assert d1 is not None
        assert d1 < d0
        assert u1 == pytest.approx(u0, rel=1e-9)


class TestHullSpinUp:
    def test_reaches_ninety_percent_corotation(self) -> None:
        run = hull_spin_up(CYL, SpinupParams(1.0), CYL.omega / 3600.0)
        t90 = time_to_corotation(run)
        assert t90 is not None
        assert t90 >= 3600.0 - 1e-9
        assert run.ramp_end == pytest.approx(3600.0)
        assert run.corotation[0] == 0.0
        assert run.corotation[-1] >= 0.9
        assert run.omega[-1] == pytest.approx(CYL.omega)
        assert np.all(run.u[1:-1].min(axis=1) <= 0.0)  # the air only lags

    def test_end_caps_speed_up_corotation(self) -> None:
        alpha = CYL.omega / 3600.0
        with_caps = time_to_corotation(hull_spin_up(CYL, SpinupParams(1.0), alpha))
        without = time_to_corotation(
            hull_spin_up(CYL, SpinupParams(1.0, endcap=False), alpha)
        )
        assert with_caps is not None
        assert without is not None
        assert with_caps < 0.2 * without

    def test_more_viscosity_means_faster_spin_up_without_caps(self) -> None:
        alpha = CYL.omega / 3600.0
        slow = time_to_corotation(
            hull_spin_up(CYL, SpinupParams(1.0, endcap=False), alpha)
        )
        fast = time_to_corotation(
            hull_spin_up(CYL, SpinupParams(10.0, endcap=False), alpha)
        )
        assert slow is not None
        assert fast is not None
        assert slow / fast == pytest.approx(10.0, rel=0.35)

    def test_faster_spin_up_gives_a_stronger_wind(self) -> None:
        params = SpinupParams(1.0)
        gentle = hull_spin_up(CYL, params, CYL.omega / 86400.0)
        brisk = hull_spin_up(CYL, params, CYL.omega / 3600.0)
        assert peak_wind_at_height(brisk, 2.0) > peak_wind_at_height(gentle, 2.0)
        assert peak_wind_at_height(brisk, 2.0) < peak_wind_at_height(brisk, 200.0)

    def test_not_reached_returns_none(self) -> None:
        run = hull_spin_up(CYL, SpinupParams(1.0), CYL.omega / 3600.0, t_max=100.0)
        assert time_to_corotation(run) is None

    def test_sudden_change_not_reached_returns_none(self) -> None:
        run = sudden_change(CYL, SpinupParams(1.0), 1e-4, t_max=10.0)
        assert time_to_adjust(run) is None


class TestValidation:
    def test_params(self) -> None:
        with pytest.raises(ValueError, match="nu_stable"):
            SpinupParams(0.0)
        with pytest.raises(ValueError, match="nu_unstable"):
            SpinupParams(1.0, nu_unstable=-1.0)
        with pytest.raises(ValueError, match="floor_drag"):
            SpinupParams(1.0, floor_drag=-0.1)
        with pytest.raises(ValueError, match="n_cells"):
            SpinupParams(1.0, n_cells=3)
        with pytest.raises(ValueError, match="stretch"):
            SpinupParams(1.0, stretch=-1.0)

    def test_nu_high_defaults_to_stable(self) -> None:
        assert SpinupParams(2.0).nu_high == 2.0
        assert SpinupParams(2.0, nu_unstable=9.0).nu_high == 9.0

    def test_simulate_arguments(self) -> None:
        params = SpinupParams(1.0, n_cells=20)
        hull = HullSpin(CYL.omega)
        zeros = np.zeros(20)
        with pytest.raises(ValueError, match="growth"):
            simulate(CYL, params, hull, zeros, 10.0, growth=0.5)
        with pytest.raises(ValueError, match="u0"):
            simulate(CYL, params, hull, np.zeros(5), 10.0)
        with pytest.raises(ValueError, match="step budget"):
            simulate(CYL, params, hull, zeros, 1.0e9, dt_max=1.0)
        with pytest.raises(ValueError, match="t_end"):
            simulate(CYL, params, hull, zeros, 0.0)

    def test_scenarios(self) -> None:
        params = SpinupParams(1.0, n_cells=20)
        with pytest.raises(ValueError, match="alpha"):
            hull_spin_up(CYL, params, 0.0)
        with pytest.raises(ValueError, match="delta_omega"):
            sudden_change(CYL, params, 0.0)
        with pytest.raises(ValueError, match="delta_omega"):
            sudden_change(CYL, params, -2.0 * CYL.omega)

    def test_hull_schedule(self) -> None:
        hull = HullSpin(1.0, 0.5, 4.0)
        assert hull.omega(2.0) == 2.0
        assert hull.omega(10.0) == 3.0
        assert HullSpin(1.0, 0.5).omega(10.0) == 6.0
