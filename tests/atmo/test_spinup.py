import math

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy.special import j1, jn_zeros

from oneill.atmo.numbers import spin_up_time
from oneill.atmo.spinup import (
    DEFAULT_FLOOR_DRAG,
    DRAG_HEIGHT,
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
from oneill.atmo.surface import SUBURBAN, drag_coefficient
from oneill.core.cylinder import Cylinder

CYL = Cylinder()
SMALL = Cylinder(radius=100.0, length=1000.0)  # fast diffusion: tests stay short


def noslip(nu: float, **kwargs: object) -> SpinupParams:
    """Analytic-test parameters: a no-slip floor (not physical at ``nu_t`` ~ 1)."""
    return SpinupParams(nu, floor_drag=None, **kwargs)  # type: ignore[arg-type]


def analytic_lag(r: np.ndarray, radius: float, alpha: float, nu: float) -> np.ndarray:
    """Steady lag u = alpha r (r^2 - R^2) / (8 nu), derived from nu L[u] = alpha r."""
    return alpha * r * (r**2 - radius**2) / (8.0 * nu)


class TestDefaults:
    def test_default_floor_is_the_suburban_log_law_drag(self) -> None:
        assert DRAG_HEIGHT == 10.0
        assert pytest.approx(0.0178, abs=1e-4) == DEFAULT_FLOOR_DRAG
        assert SpinupParams(1.0).floor_drag == DEFAULT_FLOOR_DRAG
        assert drag_coefficient(10.0, SUBURBAN.z0) == DEFAULT_FLOOR_DRAG
        # u*^2 / U^2 from the log law
        ustar = 0.4 / math.log(10.0 / 0.5)
        assert ustar**2 == pytest.approx(DEFAULT_FLOOR_DRAG)


class TestSteadyLag:
    @pytest.mark.parametrize("stretch", [0.0, 2.5])
    @pytest.mark.parametrize("nu", [1.0, 10.0])
    def test_matches_closed_form(self, nu: float, stretch: float) -> None:
        params = noslip(nu, endcap=False, n_cells=200, stretch=stretch)
        alpha = 1.0e-8
        lag = steady_lag(CYL, params, alpha)
        expected = analytic_lag(lag.radius, CYL.radius, alpha, nu)
        assert lag.u == pytest.approx(expected, abs=1e-4 * np.abs(expected).max())

    def test_lags_for_spin_up_and_leads_for_spin_down(self) -> None:
        params = noslip(1.0, endcap=False)
        assert np.all(steady_lag(CYL, params, 1e-8).u < 0.0)
        assert np.all(steady_lag(CYL, params, -1e-8).u > 0.0)

    def test_scales_inversely_with_viscosity(self) -> None:
        a = steady_lag(CYL, noslip(1.0, endcap=False), 1e-8)
        b = steady_lag(CYL, noslip(4.0, endcap=False), 1e-8)
        assert b.u == pytest.approx(a.u / 4.0, rel=1e-9)

    def test_end_caps_limit_the_lag(self) -> None:
        with_caps = steady_lag(CYL, noslip(1.0), 1e-8)
        without = steady_lag(CYL, noslip(1.0, endcap=False), 1e-8)
        assert with_caps.peak < 0.2 * without.peak
        # interior balance alpha r = -u / tau_E far from the wall
        tau = spin_up_time(1.0, CYL.omega, CYL.length / 2.0)
        mid = int(np.argmin(np.abs(with_caps.radius - 0.3 * CYL.radius)))
        assert with_caps.u[mid] == pytest.approx(
            -1e-8 * with_caps.radius[mid] * tau, rel=0.02
        )

    def test_time_marching_converges_to_the_closed_form(self) -> None:
        params = noslip(1.0, endcap=False, n_cells=160, stretch=0.0)
        alpha = 1.0e-6
        hull = HullSpin(SMALL.omega, alpha)
        run = simulate(
            SMALL, params, hull, np.zeros(160), 6000.0, dt0=5.0, growth=1.0, dt_max=5.0
        )
        expected = analytic_lag(run.radius, SMALL.radius, alpha, 1.0)
        assert run.u[-1] == pytest.approx(expected, abs=3e-3 * np.abs(expected).max())

    def test_floor_drag_torque_balance(self) -> None:
        # the drag layer carries the whole source: R^2 C_D u_r^2 ~ alpha R^4 / 4
        c_d, alpha = 2.0e-3, 1.0e-8
        params = SpinupParams(1.0, floor_drag=c_d, endcap=False)
        lag = steady_lag(CYL, params, alpha)
        u_ref = lag.wind_at_height(DRAG_HEIGHT)
        assert abs(u_ref) == pytest.approx(
            CYL.radius * math.sqrt(alpha / (4.0 * c_d)), rel=0.05
        )
        assert u_ref < 0.0

    def test_wind_at_height(self) -> None:
        lag = steady_lag(CYL, noslip(1.0, endcap=False), 1e-8)
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
        params = noslip(nu, endcap=False, n_cells=n, stretch=0.0, theta=0.5)
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

    def test_backward_euler_decay_is_first_order_accurate(self) -> None:
        j = float(jn_zeros(1, 1)[0])
        n = 120
        params = noslip(1.0, endcap=False, n_cells=n, stretch=0.0)
        rate = (j / SMALL.radius) ** 2
        centres = (np.arange(n) + 0.5) * SMALL.radius / n
        errors = []
        for steps in (100, 200):
            dt = 1.0 / rate / steps
            run = simulate(
                SMALL,
                params,
                HullSpin(SMALL.omega),
                mode_shape(centres, j),
                1.0 / rate,
                dt0=dt,
                growth=1.0,
                dt_max=dt,
            )
            phi = mode_shape(run.radius, j)
            amp = run.u @ (run.moment * phi) / (phi @ (run.moment * phi))
            errors.append(abs(-math.log(amp[-1] / amp[0]) / run.t[-1] - rate) / rate)
        assert errors[0] / errors[1] == pytest.approx(2.0, rel=0.35)

    def test_solid_body_is_an_exact_steady_state(self) -> None:
        params = SpinupParams(5.0, floor_drag=0.0, endcap=False, n_cells=40)
        run = simulate(SMALL, params, HullSpin(SMALL.omega), np.zeros(40), 100.0)
        assert np.max(np.abs(run.u)) == 0.0
        rigid = 0.003 * run.radius
        free = simulate(SMALL, params, HullSpin(SMALL.omega), rigid, 100.0)
        assert free.u[-1] == pytest.approx(rigid, abs=1e-12)

    def test_hull_acceleration_drives_the_air_the_right_way(self) -> None:
        # free slip, no caps: nothing touches the air, so it stays inertially at
        # rest and the relative wind is exactly -(Omega - Omega_0) r. A flipped
        # sign of the hull term would give +(Omega - Omega_0) r.
        params = SpinupParams(3.0, floor_drag=0.0, endcap=False, n_cells=40)
        alpha = 1.0e-3
        run = simulate(SMALL, params, HullSpin(SMALL.omega, alpha), np.zeros(40), 200.0)
        quarter = SMALL.radius**4 / 4.0
        assert run.u[-1] @ run.moment == pytest.approx(
            -alpha * run.t[-1] * quarter, rel=1e-9
        )
        # cell winds are r^2-weighted means of -alpha t r, close to the centre value
        assert run.u[-1] == pytest.approx(-alpha * run.t[-1] * run.radius, rel=3e-2)
        assert np.all(run.u[-1] < 0.0)
        assert np.all(run.budget.hull < 0.0)
        assert np.all(run.budget.wall == 0.0)

    def test_spin_down_drives_the_air_forward(self) -> None:
        params = SpinupParams(3.0, floor_drag=0.0, endcap=False, n_cells=40)
        run = simulate(
            SMALL, params, HullSpin(SMALL.omega, -1.0e-3), np.zeros(40), 100.0
        )
        assert np.all(run.u[-1] > 0.0)
        assert np.all(run.budget.hull > 0.0)


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
        assert run.u[-1] == pytest.approx(u0 * math.exp(-1.0), rel=3e-3)

    def test_longer_habitat_pumps_more_slowly(self) -> None:
        params = SpinupParams(1.0, floor_drag=0.0, n_cells=20)
        rates = []
        for length in (2000.0, 4000.0):
            cyl = Cylinder(radius=500.0, length=length)
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
        assert rates[0] / rates[1] == pytest.approx(2.0, rel=1e-2)

    def test_cap_viscosity_can_be_set_separately(self) -> None:
        rates = []
        for nu_cap in (None, 4.0):
            params = SpinupParams(1.0, floor_drag=0.0, nu_endcap=nu_cap, n_cells=20)
            grid = simulate(CYL, params, HullSpin(CYL.omega), np.zeros(20), 1.0)
            u0 = 0.01 * grid.radius
            run = simulate(
                CYL,
                params,
                HullSpin(CYL.omega),
                u0,
                1000.0,
                dt0=10.0,
                growth=1.0,
                dt_max=10.0,
            )
            rates.append(-math.log(run.u[-1, 5] / u0[5]) / 1000.0)
        assert rates[1] / rates[0] == pytest.approx(2.0, rel=1e-2)

    def test_greenspan_howard_interior_spin_up(self) -> None:
        # small step in hull rate; mid-radius wind decays on tau_E (diffusion from
        # the floor has not reached it yet)
        tau = spin_up_time(1.0, CYL.omega, CYL.length / 2.0)
        run = sudden_change(CYL, SpinupParams(1.0), 1e-4, t_max=tau, dt_max=tau / 400)
        mid = int(np.argmin(np.abs(run.radius - 0.5 * CYL.radius)))
        ratio = run.u[-1, mid] / run.u[0, mid]
        assert ratio == pytest.approx(math.exp(-run.t[-1] / tau), rel=0.05)

    def test_sink_vanishes_without_rotation(self) -> None:
        # spin-up from rest: the rate sqrt(nu Omega)/(L/2) is finite at Omega = 0
        params = SpinupParams(1.0, floor_drag=0.0, n_cells=20)
        run = simulate(SMALL, params, HullSpin(0.0), np.zeros(20), 50.0)
        assert np.all(np.isfinite(run.u))
        assert np.all(run.budget.endcap == 0.0)

    def test_ekman_time_scales_with_the_new_hull_rate(self) -> None:
        # equal viscosities, caps on: spin-down (0.9 Omega) adjusts more slowly
        # than spin-up (1.1 Omega) because tau_E ~ 1 / sqrt(Omega)
        params = noslip(1.0, n_cells=100)
        down = time_to_adjust(sudden_change(CYL, params, -0.1 * CYL.omega))
        up = time_to_adjust(sudden_change(CYL, params, 0.1 * CYL.omega))
        assert down is not None
        assert up is not None
        assert down / up == pytest.approx(math.sqrt(1.1 / 0.9), rel=0.03)


class TestBudget:
    @pytest.mark.parametrize("theta", [1.0, 0.5])
    @pytest.mark.parametrize("floor", [None, 0.02])
    def test_angular_momentum_budget_closes_exactly(
        self, floor: float | None, theta: float
    ) -> None:
        params = SpinupParams(
            1.0, nu_unstable=8.0, floor_drag=floor, n_cells=60, theta=theta
        )
        run = hull_spin_up(SMALL, params, SMALL.omega / 200.0, fraction=0.95)
        b = run.budget
        change = np.diff(run.u @ run.moment)
        flux = b.dt * (b.wall + b.endcap + b.hull)
        assert change == pytest.approx(flux, rel=1e-9, abs=1e-9 * np.abs(flux).max())

    def test_total_torque_matches_sudden_change(self) -> None:
        params = noslip(1.0, endcap=False, n_cells=60)
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
    """Rayleigh-switched viscosity: converged, smooth, and only mildly asymmetric."""

    def adjust(self, delta: float, dt_max: float = 20.0, n_cells: int = 100) -> float:
        params = noslip(1.0, nu_unstable=10.0, endcap=False, n_cells=n_cells)
        t = time_to_adjust(
            sudden_change(SMALL, params, delta * SMALL.omega, dt_max=dt_max, t_max=5e4)
        )
        assert t is not None
        return t

    def test_time_to_adjust_converges_under_dt_refinement(self) -> None:
        # first-order (backward Euler) convergence: dt versus dt / 4
        for delta in (-0.1, 0.1):
            coarse, mid, fine = (self.adjust(delta, dt) for dt in (20.0, 5.0, 1.25))
            assert coarse == pytest.approx(mid, rel=0.01)
            assert mid == pytest.approx(fine, rel=0.01)
        strong = [self.adjust(-0.5, dt) for dt in (5.0, 1.25)]
        assert strong[0] == pytest.approx(strong[1], rel=0.02)

    def test_no_odd_even_oscillation_near_the_wall(self) -> None:
        # a hard regime switch with Crank-Nicolson rang at +-12 m/s with ~50 slope
        # sign changes; the profile of a strong spin-down must be smooth
        params = noslip(1.0, nu_unstable=10.0, endcap=False, n_cells=100)
        run = sudden_change(SMALL, params, -0.5 * SMALL.omega, dt_max=20.0, t_max=5e4)
        worst = 0
        for row in run.u[::5]:
            slope = np.diff(row)
            signs = np.sign(slope[np.abs(slope) > 1e-9 * np.abs(row).max()])
            worst = max(worst, int(np.count_nonzero(np.diff(signs))))
        assert worst <= 4
        assert np.max(np.abs(run.u)) <= 0.5 * SMALL.omega * SMALL.radius * (1 + 1e-6)

    def test_small_changes_are_nearly_symmetric(self) -> None:
        # the interior is Rayleigh-stable in a relative solid-body change; only a
        # thin wall layer sees the larger viscosity
        down, up = self.adjust(-0.1), self.adjust(0.1)
        assert down / up == pytest.approx(1.0, abs=0.03)

    def test_strong_spin_down_is_faster_when_the_regimes_differ(self) -> None:
        assert self.adjust(-0.5) < 0.95 * self.adjust(0.5)

    def test_symmetric_when_regimes_agree(self) -> None:
        params = noslip(3.0, endcap=False, n_cells=100)
        down = time_to_adjust(sudden_change(SMALL, params, -0.5 * SMALL.omega))
        up = time_to_adjust(sudden_change(SMALL, params, +0.5 * SMALL.omega))
        assert down == pytest.approx(up, rel=1e-6)

    @settings(max_examples=5, deadline=None)
    @given(st.floats(2.0, 20.0))
    def test_unstable_viscosity_matters_mainly_for_spin_down(
        self, ratio: float
    ) -> None:
        base = noslip(1.0, endcap=False, n_cells=60)
        fast = noslip(1.0, nu_unstable=ratio, endcap=False, n_cells=60)
        d0 = time_to_adjust(sudden_change(SMALL, base, -0.5 * SMALL.omega))
        d1 = time_to_adjust(sudden_change(SMALL, fast, -0.5 * SMALL.omega))
        u0 = time_to_adjust(sudden_change(SMALL, base, 0.5 * SMALL.omega))
        u1 = time_to_adjust(sudden_change(SMALL, fast, 0.5 * SMALL.omega))
        assert d0 is not None
        assert d1 is not None
        assert d1 < d0
        assert u1 is not None
        assert u0 is not None
        assert (u0 - u1) / u0 < 0.5 * (d0 - d1) / d0  # spin-up barely notices


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
        assert slow / fast == pytest.approx(10.0, rel=0.4)

    def test_faster_spin_up_gives_a_stronger_wind(self) -> None:
        params = SpinupParams(1.0)
        gentle = hull_spin_up(CYL, params, CYL.omega / 86400.0)
        brisk = hull_spin_up(CYL, params, CYL.omega / 3600.0)
        z_r = DRAG_HEIGHT
        assert peak_wind_at_height(brisk, z_r) > peak_wind_at_height(gentle, z_r)
        assert peak_wind_at_height(brisk, z_r) < peak_wind_at_height(brisk, 200.0)

    def test_floor_drag_caps_the_wind_at_the_reference_height(self) -> None:
        alpha = CYL.omega / 3600.0
        draggy = hull_spin_up(CYL, SpinupParams(1.0), alpha)
        slippery = hull_spin_up(CYL, SpinupParams(1.0, floor_drag=0.0), alpha)
        assert peak_wind_at_height(draggy, 10.0) < peak_wind_at_height(slippery, 10.0)

    def test_rossby_number_flags_the_linear_sink(self) -> None:
        fast = hull_spin_up(CYL, SpinupParams(1.0), CYL.omega / 3600.0)
        slow = hull_spin_up(CYL, SpinupParams(10.0), CYL.omega / 86400.0)
        assert fast.max_rossby == pytest.approx(0.55, abs=0.05)
        assert not fast.sink_valid
        assert slow.max_rossby < 0.3
        assert slow.sink_valid

    def test_corotation_time_never_precedes_the_end_of_the_ramp(self) -> None:
        run = hull_spin_up(CYL, SpinupParams(1.0), CYL.omega / 3600.0)
        t_any = time_to_corotation(run, fraction=0.0)
        assert t_any is not None
        assert t_any >= run.ramp_end
        assert t_any <= run.ramp_end + 700.0  # the first sample after the ramp

    def test_crossing_times_are_interpolated_between_steps(self) -> None:
        run = hull_spin_up(CYL, SpinupParams(1.0), CYL.omega / 3600.0)
        t90 = time_to_corotation(run)
        assert t90 is not None
        assert t90 not in run.t  # lies strictly between two samples

    def test_not_reached_returns_none(self) -> None:
        run = hull_spin_up(CYL, SpinupParams(1.0), CYL.omega / 3600.0, t_max=100.0)
        assert time_to_corotation(run) is None

    def test_sudden_change_not_reached_returns_none(self) -> None:
        run = sudden_change(CYL, SpinupParams(1.0), 1e-4, t_max=10.0)
        assert time_to_adjust(run) is None

    def test_peak_wind_excludes_the_initial_jump(self) -> None:
        run = sudden_change(SMALL, SpinupParams(1.0, n_cells=60), 0.1 * SMALL.omega)
        r_ref = SMALL.radius - DRAG_HEIGHT
        jump = abs(float(np.interp(r_ref, run.radius, run.u[0])))
        assert jump == pytest.approx(0.1 * SMALL.omega * r_ref, rel=1e-3)
        assert peak_wind_at_height(run, DRAG_HEIGHT) < jump


class TestValidation:
    def test_params(self) -> None:
        with pytest.raises(ValueError, match="nu_stable"):
            SpinupParams(0.0)
        with pytest.raises(ValueError, match="nu_unstable"):
            SpinupParams(1.0, nu_unstable=-1.0)
        with pytest.raises(ValueError, match="nu_endcap"):
            SpinupParams(1.0, nu_endcap=0.0)
        with pytest.raises(ValueError, match="drag_height"):
            SpinupParams(1.0, drag_height=0.0)
        with pytest.raises(ValueError, match="floor_drag"):
            SpinupParams(1.0, floor_drag=-0.1)
        with pytest.raises(ValueError, match="n_cells"):
            SpinupParams(1.0, n_cells=3)
        with pytest.raises(ValueError, match="stretch"):
            SpinupParams(1.0, stretch=-1.0)
        with pytest.raises(ValueError, match="theta"):
            SpinupParams(1.0, theta=0.4)

    def test_nu_high_and_nu_cap_default_to_stable(self) -> None:
        assert SpinupParams(2.0).nu_high == 2.0
        assert SpinupParams(2.0, nu_unstable=9.0).nu_high == 9.0
        assert SpinupParams(2.0).nu_cap == 2.0
        assert SpinupParams(2.0, nu_endcap=5.0).nu_cap == 5.0

    def test_drag_coefficient_needs_height_above_z0(self) -> None:
        with pytest.raises(ValueError, match="exceed"):
            drag_coefficient(0.5, 0.5)

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
