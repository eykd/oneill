import math

import numpy as np
import pytest

from oneill.atmo.annulus2d import (
    DEFAULT_CFL,
    Annulus,
    AnnulusConfig,
    AnnulusState,
    Forcing,
    _lagrange_weights,
)
from oneill.atmo.gas import ONEILL_1974
from oneill.atmo.numbers import spin_up_time
from oneill.atmo.spinup import SpinupParams, sudden_change
from oneill.core.cylinder import Cylinder

OMEGA = Cylinder().omega


def small(**kwargs: object) -> AnnulusConfig:
    """A fast configuration (24 x 32) for the default habitat."""
    base: dict[str, object] = {"n_r": 24, "n_theta": 32}
    base.update(kwargs)
    return AnnulusConfig(**base)  # type: ignore[arg-type]


def blob(model: Annulus, r0: float, theta0: float, sr: float, sth: float) -> np.ndarray:
    """A smooth Gaussian bump on the model grid (physical space)."""
    r = model.r[:, None]
    th = model.theta[None, :]
    dth = np.angle(np.exp(1j * (th - theta0)))
    return np.asarray(np.exp(-(((r - r0) / sr) ** 2) - ((dth * r0) / sth) ** 2))


class TestForcing:
    def test_strips_alternate_land_and_window(self) -> None:
        f = Forcing(land_flux=100.0)
        theta = (np.arange(6) + 0.5) * math.pi / 3.0
        assert f.flux(theta, 0.0) == pytest.approx([100, -100, 100, -100, 100, -100])

    def test_neutral_windows_balance_the_land_exactly(self) -> None:
        theta = np.linspace(0.0, 2.0 * math.pi, 600, endpoint=False)
        for t in (0.0, 1000.0, 30000.0):
            f = Forcing(land_flux=200.0, night_flux=-30.0, period=86400.0)
            assert np.mean(f.flux(theta, t)) == pytest.approx(0.0, abs=1e-9)

    def test_zero_windows_and_uniform(self) -> None:
        theta = (np.arange(6) + 0.5) * math.pi / 3.0
        assert Forcing(window="zero").flux(theta, 0.0) == pytest.approx(
            [100, 0, 100, 0, 100, 0]
        )
        assert Forcing(window="uniform").flux(theta, 0.0) == pytest.approx(100.0)

    def test_day_night_cycle(self) -> None:
        f = Forcing(land_flux=200.0, night_flux=-30.0, period=86400.0)
        assert f.land_level(0.0) == 0.0
        assert f.land_level(21600.0) == pytest.approx(200.0)
        assert f.land_level(64800.0) == pytest.approx(-30.0)
        assert Forcing(land_flux=70.0).land_level(12345.0) == 70.0

    def test_smoothed_edges_keep_the_pattern_and_the_zero_net_flux(self) -> None:
        theta = np.linspace(0.0, 2.0 * math.pi, 720, endpoint=False)
        hard = Forcing(land_flux=200.0).flux(theta, 0.0)
        smooth = Forcing(land_flux=200.0, edge_width=0.05).flux(theta, 0.0)
        assert np.mean(smooth) == pytest.approx(0.0, abs=1e-9)
        assert np.max(np.abs(np.diff(smooth))) < 0.5 * np.max(np.abs(np.diff(hard)))
        interior = np.abs(np.sin(3.0 * theta)) > 0.6
        assert smooth[interior] == pytest.approx(hard[interior], abs=2.0)
        # the edge value is the mean of land and window
        edge = Forcing(land_flux=200.0, edge_width=0.05).flux(
            np.array([math.pi / 3.0]), 0.0
        )
        assert edge[0] == pytest.approx(0.0, abs=1e-9)
        tiny = Forcing(land_flux=200.0, edge_width=1e-6).flux(theta, 0.0)
        off_edge = np.abs(np.sin(3.0 * theta)) > 0.01
        assert tiny[off_edge] == pytest.approx(hard[off_edge], abs=1e-3)

    def test_phase_shifts_the_pattern(self) -> None:
        theta = np.array([0.1, 0.1 + math.pi / 3.0])
        shifted = Forcing(phase=math.pi / 3.0)
        assert shifted.flux(theta, 0.0) == pytest.approx([-100.0, 100.0])

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"strips": 0},
            {"period": 0.0},
            {"window": "bogus"},
            {"land_flux": math.nan},
            {"edge_width": -0.1},
        ],
    )
    def test_invalid(self, kwargs: dict[str, object]) -> None:
        with pytest.raises(
            ValueError, match=r"must|strips|window|land_flux|edge_width"
        ):
            Forcing(**kwargs)  # type: ignore[arg-type]


class TestConfig:
    def test_defaults_and_derived_rates(self) -> None:
        cfg = AnnulusConfig()
        assert cfg.kappa == 1.0
        tau = spin_up_time(1.0, cfg.omega, cfg.length / 2.0)
        assert cfg.endcap_rate == pytest.approx(1.0 / tau)
        assert AnnulusConfig(endcap=False).endcap_rate == 0.0

    def test_from_habitat_takes_beta_and_rho_cp_from_the_gas(self) -> None:
        gas = ONEILL_1974.at_temperature(290.0)
        cfg = AnnulusConfig.from_habitat(Cylinder(), gas, nu=3.0)
        assert cfg.beta == pytest.approx(1.0 / 290.0)
        assert cfg.rho_cp == pytest.approx(gas.rho_rim * gas.cp)
        assert cfg.omega == Cylinder().omega
        assert cfg.nu == 3.0

    @pytest.mark.parametrize(
        "kwargs",
        [
            {"radius": -1.0},
            {"floor_drag": -0.1},
            {"r_in_fraction": 1.0},
            {"n_r": 2},
            {"n_theta": 4},
            {"wall": "bogus"},
            {"outer_t": "bogus"},
            {"inner_t": "bogus"},
            {"stretch": -1.0},
            {"tau_relax": 0.0},
        ],
    )
    def test_invalid(self, kwargs: dict[str, object]) -> None:
        with pytest.raises(ValueError, match=r"must|need"):
            AnnulusConfig(**kwargs)  # type: ignore[arg-type]


class TestGridAndOperators:
    def test_grid_clusters_toward_the_floor(self) -> None:
        model = Annulus(AnnulusConfig())
        assert model.faces[0] == pytest.approx(320.0)
        assert model.faces[-1] == pytest.approx(3200.0)
        assert model.width[-1] < 0.3 * model.width[0]
        assert 6.0 < model.width[-1] < 12.0  # the plan's ~8 m first cell
        assert np.all(model.width > 0.0)

    def test_uniform_grid_without_stretch(self) -> None:
        model = Annulus(small(stretch=0.0))
        assert model.width == pytest.approx(model.width[0])

    def test_radial_derivative_is_second_order_exact_on_quadratics(self) -> None:
        model = Annulus(small())
        f = (model.r**2)[:, None] * np.ones((1, model.k))
        f_in = np.full(model.k, model.r_in**2, dtype=complex)
        f_out = np.full(model.k, model.config.radius**2, dtype=complex)
        d = model._radial(f.astype(complex), f_in, f_out)  # noqa: SLF001
        assert d[:, 0].real == pytest.approx(2.0 * model.r)

    def test_face_interpolation_is_third_order_upwind_and_linear(self) -> None:
        model = Annulus(small(stretch=0.0))
        x = np.concatenate([[model.r_in], model.r, [model.config.radius]])
        faces = model.faces[1:-1]
        quad = x**2
        assert model.f_plus @ quad == pytest.approx(faces**2, rel=1e-9)
        assert model.f_minus @ quad == pytest.approx(faces**2, rel=1e-9)
        assert model.f_linear @ x == pytest.approx(faces)
        assert model.f_linear.sum(axis=1) == pytest.approx(1.0)

    def test_lagrange_weights_differentiate_polynomials(self) -> None:
        xs = np.array([0.0, 1.0, 2.5, 4.0])
        w = _lagrange_weights(xs, 1)
        assert w @ xs**2 == pytest.approx(2.0)
        assert w @ np.ones(4) == pytest.approx(0.0, abs=1e-12)

    def test_poisson_solver_matches_an_analytic_mode_to_second_order(self) -> None:
        errors = []
        for n in (24, 48):
            model = Annulus(small(n_r=n, n_theta=16, stretch=0.0, radius=1.0))
            r, r_in = model.r, model.r_in
            k = math.pi / (1.0 - r_in)
            psi = np.sin(k * (r - r_in))
            m = 2
            zeta = -(-(k**2) * psi + k * np.cos(k * (r - r_in)) / r - m**2 * psi / r**2)
            z = np.zeros((n, model.k), dtype=complex)
            z[:, m] = zeta
            errors.append(np.max(np.abs(model.psi(z)[:, m] - psi)))
        assert errors[0] / errors[1] == pytest.approx(4.0, rel=0.1)
        assert errors[1] < 1e-3

    def test_transforms_round_trip_and_mean_is_column_zero(self) -> None:
        model = Annulus(small())
        rng = np.random.default_rng(1)
        x = rng.standard_normal((model.config.n_r, model.n_theta))
        x_hat = model.to_hat(x)
        assert x_hat[:, 0].real == pytest.approx(x.mean(axis=1))
        assert model.to_phys(x_hat) == pytest.approx(x)

    def test_radial_cut_of_the_polar_filter(self) -> None:
        model = Annulus(small(n_theta=64))
        assert model.mask[0].sum() < model.mask[-1].sum()
        full = Annulus(small(n_theta=64, polar_filter=False))
        assert full.mask.all(axis=0).sum() == full.mask[0].sum()
        only = Annulus(small(n_theta=64, modes=(3,), polar_filter=False))
        assert np.flatnonzero(only.mask[0]).tolist() == [0, 3]


def run_steps(model: Annulus, state: AnnulusState, dt: float, n: int) -> None:
    for _ in range(n):
        model.step(state, dt)


class TestHeatBudget:
    def test_wall_flux_is_the_only_heat_source_exactly(self) -> None:
        cfg = small(forcing=Forcing(land_flux=100.0, window="zero"), nu=1.0)
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        run_steps(model, state, 5.0, 40)
        mean_flux = np.mean(cfg.forcing.flux(model.theta, 0.0)) / cfg.rho_cp
        expected = 2.0 * math.pi * cfg.radius * mean_flux * state.t
        assert model.heat_content(state) == pytest.approx(expected, rel=1e-9)

    def test_flux_neutral_windows_leave_the_heat_content_unchanged(self) -> None:
        model = Annulus(small(forcing=Forcing(land_flux=200.0)))
        state = model.initial_state(amplitude=0.0)
        run_steps(model, state, 5.0, 40)
        assert model.heat_content(state) == pytest.approx(0.0, abs=1e-6)

    def test_relaxation_removes_the_area_mean_at_the_given_rate(self) -> None:
        cfg = small(
            forcing=Forcing(land_flux=0.0, window="zero"),
            tau_relax=500.0,
            advection=False,
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        state.temp[:, 0] = 2.0
        run_steps(model, state, 1.0, 500)
        mean = np.sum(model.area * state.temp[:, 0].real) / np.sum(model.area)
        assert mean == pytest.approx(2.0 * math.exp(-1.0), rel=2e-3)

    def test_axisymmetric_heating_drives_no_flow_and_a_conduction_profile(
        self,
    ) -> None:
        errors = []
        for n_r in (24, 48):
            cfg = small(
                n_r=n_r,
                forcing=Forcing(land_flux=200.0, window="uniform"),
                nu=1.0,
                coriolis_buoyancy=False,
            )
            model = Annulus(cfg)
            state = model.initial_state(amplitude=0.0)
            run_steps(model, state, 10.0, 100)
            assert np.max(np.abs(state.zeta)) == 0.0
            assert np.max(np.abs(state.temp[:, 1:])) < 1e-12
            assert np.max(np.abs(state.ubar)) == 0.0
            # semi-infinite conduction: wall temperature 2 q sqrt(kappa t / pi)
            q = 200.0 / cfg.rho_cp
            t_wall = state.temp[-1, 0].real + q / cfg.kappa * (cfg.radius - model.r[-1])
            exact = 2.0 * q * math.sqrt(cfg.kappa * state.t / math.pi)
            errors.append(abs(t_wall / exact - 1.0))
        assert errors[0] < 0.08
        assert errors[1] < 0.5 * errors[0]  # converges under radial refinement

    def test_fixed_temperature_walls_and_an_inner_flux(self) -> None:
        cfg = small(
            outer_t="fixed",
            inner_t="fixed",
            t_outer=1.0,
            t_inner=-1.0,
            advection=False,
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        run_steps(model, state, 20.0, 400)
        profile = state.temp[:, 0].real
        assert profile[-1] > 0.0 > profile[0]
        assert np.all(np.diff(profile) > 0.0)
        flux_cfg = small(
            inner_t="flux", flux_inner=50.0, forcing=Forcing(land_flux=0.0)
        )
        fm = Annulus(flux_cfg)
        fs = fm.initial_state(amplitude=0.0)
        run_steps(fm, fs, 5.0, 40)
        expected = 2.0 * math.pi * fm.r_in * (50.0 / flux_cfg.rho_cp) * fs.t
        assert fm.heat_content(fs) == pytest.approx(-expected, rel=1e-9)


class TestDragAndDamping:
    def test_end_cap_damping_decays_the_mean_wind_at_one_over_tau_e(self) -> None:
        cfg = small(
            length=2000.0,
            nu=1.0,
            floor_drag=0.0,
            advection=False,
            n_theta=8,
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = 0.001 * model.r
        tau = 1.0 / cfg.endcap_rate
        run_steps(model, state, tau / 400.0, 400)
        assert state.ubar == pytest.approx(0.001 * model.r * math.exp(-1.0), rel=2e-3)

    def test_vorticity_modes_decay_at_the_same_rate(self) -> None:
        cfg = small(
            length=2000.0,
            nu=1e-3,
            floor_drag=0.0,
            advection=False,
            n_theta=16,
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        shape = np.sin(math.pi * (model.r - model.r_in) / (cfg.radius - model.r_in))
        state.zeta[:, 2] = 1e-5 * shape
        tau = 1.0 / cfg.endcap_rate
        run_steps(model, state, tau / 400.0, 400)
        ratio = np.max(np.abs(state.zeta[:, 2])) / 1e-5
        assert ratio == pytest.approx(math.exp(-1.0), rel=5e-3)

    def test_longer_habitat_damps_more_slowly(self) -> None:
        rates = []
        for length in (2000.0, 4000.0):
            cfg = small(
                length=length,
                floor_drag=0.0,
                advection=False,
                n_theta=8,
                forcing=Forcing(land_flux=0.0),
            )
            model = Annulus(cfg)
            state = model.initial_state(amplitude=0.0)
            state.ubar[:] = 0.001 * model.r
            run_steps(model, state, 10.0, 200)
            rates.append(-math.log(state.ubar[10] / (0.001 * model.r[10])) / state.t)
        assert rates[0] / rates[1] == pytest.approx(2.0, rel=1e-2)

    def test_floor_drag_acts_on_the_bottom_layer_only(self) -> None:
        cfg = small(
            endcap=False,
            advection=False,
            n_theta=8,
            nu=1e-3,
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = 0.5 * model.r / cfg.radius
        run_steps(model, state, 5.0, 20)
        slowed = state.ubar < 0.5 * model.r / cfg.radius - 1e-12
        assert slowed[-1]
        assert not slowed[:-3].any()

    def test_no_drag_means_no_drag_torque(self) -> None:
        model = Annulus(small(floor_drag=0.0))
        state = model.initial_state()
        assert np.all(model.mean_wind_terms(state).drag == 0.0)


def jet_state(model: Annulus, centre: float, width: float, peak: float) -> AnnulusState:
    state = model.initial_state(amplitude=0.0)
    state.ubar[:] = peak / np.cosh((model.r - centre) / width) ** 2
    return state


def inviscid(**kwargs: object) -> AnnulusConfig:
    """No viscosity worth speaking of, no drag, no caps, no buoyancy source."""
    base: dict[str, object] = {
        "nu": 1e-6,
        "floor_drag": 0.0,
        "endcap": False,
        "polar_filter": False,
        "forcing": Forcing(land_flux=0.0),
        "coriolis_buoyancy": False,
        "radius": 1000.0,
        "r_in_fraction": 0.2,
        "stretch": 0.0,
        "n_r": 48,
        "n_theta": 64,
    }
    base.update(kwargs)
    return AnnulusConfig(**base)  # type: ignore[arg-type]


class TestInviscidDynamics:
    def test_smooth_vortices_conserve_energy_and_enstrophy(self) -> None:
        model = Annulus(inviscid())
        state = model.initial_state(amplitude=0.0)
        shape = np.sin(math.pi * (model.r - model.r_in) / (1000.0 - model.r_in))
        for m, amp in ((2, 4e-3), (3, 3e-3)):
            state.zeta[:, m] = amp * shape
        e0, z0 = model.energy(state), model.enstrophy(state)
        for _ in range(60):
            model.step(state, model.suggest_dt(state, dt_max=20.0))
        assert model.energy(state) == pytest.approx(e0, rel=5e-3)
        assert model.enstrophy(state) == pytest.approx(z0, rel=2e-2)

    def test_shear_instability_moves_energy_from_the_mean_wind_into_eddies(
        self,
    ) -> None:
        model = Annulus(inviscid())
        state = jet_state(model, 600.0, 60.0, 1.0)
        shape = np.sin(math.pi * (model.r - model.r_in) / (1000.0 - model.r_in))
        state.zeta[:, 3] = 1e-7 * shape
        total0 = model.energy(state)
        mean0 = model.mean_wind_energy(state)
        assert total0 == pytest.approx(mean0, rel=1e-6)
        while state.t < 6000.0:
            model.step(state, model.suggest_dt(state, dt_max=20.0))
        mean = model.mean_wind_energy(state)
        eddy = model.energy(state) - mean
        assert eddy > 0.1 * total0  # grew from ~1e-14 of the total
        assert mean < mean0  # and took it from the jet
        # at the default CFL the saturated jet keeps its energy to well under 1 %
        assert model.energy(state) == pytest.approx(total0, rel=1e-2)

    def test_mean_wind_alone_is_a_steady_solution_without_viscosity(self) -> None:
        model = Annulus(inviscid())
        state = jet_state(model, 600.0, 60.0, 1.0)
        before = state.ubar.copy()
        for _ in range(20):
            model.step(state, 10.0)
        assert state.ubar == pytest.approx(before, abs=1e-6)
        assert np.max(np.abs(state.zeta)) == 0.0

    def test_polar_filter_is_part_of_the_energy_bookkeeping(self) -> None:
        # with the filter on, high azimuthal modes near the inner wall are removed
        filtered = Annulus(inviscid(polar_filter=True))
        plain = Annulus(inviscid(polar_filter=False))
        assert filtered.mask.sum() < plain.mask.sum()


class TestAngularMomentum:
    def test_reynolds_torque_vanishes_and_cb_torque_is_the_inertia_change(self) -> None:
        cfg = small(n_theta=48, forcing=Forcing(land_flux=200.0))
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.05, seed=3)
        for _ in range(30):
            model.step(state, 5.0)
            terms = model.mean_wind_terms(state)
            assert np.sum(model.moment * terms.reynolds) == pytest.approx(
                0.0, abs=1e-9 * (1.0 + np.sum(model.moment * np.abs(terms.reynolds)))
            )
            d_t0 = model.temperature_tendency(state)[:, 0].real
            g = d_t0 - np.sum(model.area * d_t0) / np.sum(model.area)
            expected = -cfg.omega * cfg.beta * np.sum(model.inertia_weights * g)
            assert np.sum(model.moment * terms.coriolis_buoyancy) == pytest.approx(
                expected, rel=1e-9, abs=1e-18
            )

    def test_vorticity_source_is_minus_two_omega_beta_dt_dt(self) -> None:
        on = Annulus(small(n_theta=48, forcing=Forcing(land_flux=200.0)))
        off = Annulus(
            small(n_theta=48, forcing=Forcing(land_flux=200.0), coriolis_buoyancy=False)
        )
        state = on.initial_state(amplitude=0.05, seed=4)
        for _ in range(20):
            on.step(state, 5.0)
        (n_on, _, _), _, _ = on.explicit_terms(state)
        (n_off, _, _), _, _ = off.explicit_terms(state)
        d_t = on.temperature_tendency(state)
        expected = -2.0 * on.config.omega * on.config.beta * d_t
        expected[:, 0] = 0.0
        assert n_on - n_off == pytest.approx(expected, rel=1e-9, abs=1e-18)
        assert np.max(np.abs(expected)) > 0.0

    def test_a_steady_state_feels_no_coriolis_buoyancy(self) -> None:
        cfg = small(
            nu=500.0,
            advection=False,
            outer_t="fixed",
            inner_t="fixed",
            t_outer=1.0,
            t_inner=0.0,
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        state.temp[:, 0] = np.linspace(2.0, -1.0, cfg.n_r)  # far from steady
        early = np.max(np.abs(model.mean_wind_terms(state).coriolis_buoyancy))
        run_steps(model, state, 300.0, 600)
        late = np.max(np.abs(model.mean_wind_terms(state).coriolis_buoyancy))
        assert early > 1e-9
        assert late < 1e-3 * early
        assert np.max(np.abs(model.temperature_tendency(state)[:, 0])) < 1e-4

    def test_the_mean_wind_with_coriolis_buoyancy_differs_only_while_heat_moves(
        self,
    ) -> None:
        cfg = small(
            n_theta=8,
            advection=False,
            floor_drag=0.0,
            endcap=False,
            forcing=Forcing(land_flux=100.0, window="uniform"),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        run_steps(model, state, 10.0, 60)
        # uniform heating of a closed fluid expands it: the mean wind responds ...
        assert np.max(np.abs(state.ubar)) > 0.0
        # ... but the absolute angular momentum is untouched (no wall or cap torque)
        relative = abs(model.angular_momentum(state))
        assert relative > 0.0
        assert abs(model.absolute_angular_momentum(state)) < 1e-4 * relative  # O(dt)

    def test_absolute_angular_momentum_changes_only_by_wall_and_cap_torques(
        self,
    ) -> None:
        cfg = small(n_theta=48, forcing=Forcing(land_flux=200.0), length=4000.0)
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.05, seed=5)
        dt = 4.0
        for _ in range(40):  # let a mean wind develop
            model.step(state, dt)
        ledger = []
        for _ in range(60):
            l_before = model.absolute_angular_momentum(state)
            rel_before = model.angular_momentum(state)
            terms_before = model.mean_wind_terms(state)
            model.step(state, dt)
            l_after = model.absolute_angular_momentum(state)
            rel_after = model.angular_momentum(state)
            terms_after = model.mean_wind_terms(state)
            drag = (
                2.0
                * math.pi
                * np.sum(model.moment * 0.5 * (terms_before.drag + terms_after.drag))
            )
            cap = -cfg.endcap_rate * 0.5 * (rel_before + rel_after)
            ledger.append((l_after - l_before, dt * (drag + cap), abs(rel_after)))
        change = np.array([x[0] for x in ledger])
        torque = np.array([x[1] for x in ledger])
        scale = max(x[2] for x in ledger)
        assert np.max(np.abs(change - torque)) < 2e-3 * scale

    def test_omega_only_enters_through_beta_omega_squared_and_tau_e(self) -> None:
        base = small(n_theta=32, coriolis_buoyancy=False, length=4000.0)
        faster = small(
            n_theta=32,
            coriolis_buoyancy=False,
            omega=2.0 * base.omega,
            beta=base.beta / 4.0,
            length=4000.0 * math.sqrt(2.0),
        )
        assert faster.endcap_rate == pytest.approx(base.endcap_rate)
        a, b = Annulus(base), Annulus(faster)
        sa, sb = (
            a.initial_state(amplitude=0.05, seed=2),
            b.initial_state(amplitude=0.05, seed=2),
        )
        for _ in range(40):
            a.step(sa, 6.0)
            b.step(sb, 6.0)
        assert sb.temp == pytest.approx(sa.temp, rel=1e-9, abs=1e-12)
        assert sb.zeta == pytest.approx(sa.zeta, rel=1e-9, abs=1e-14)
        assert sb.ubar == pytest.approx(sa.ubar, rel=1e-9, abs=1e-14)


class TestPhysicalSigns:
    def test_a_warm_patch_rises_toward_the_axis_and_a_cold_one_sinks(self) -> None:
        cfg = small(
            n_r=32,
            n_theta=48,
            nu=5.0,
            forcing=Forcing(land_flux=0.0),
            floor_drag=0.0,
            endcap=False,
        )
        results = []
        for sign in (1.0, -1.0):
            model = Annulus(cfg)
            state = model.initial_state(amplitude=0.0)
            patch = blob(model, 2000.0, math.pi, 150.0, 150.0)
            state.temp = np.where(model.mask, model.to_hat(sign * 2.0 * patch), 0.0)
            for _ in range(15):
                model.step(state, 4.0)
            u_r, _, _ = model.velocities(state)
            results.append(float(np.sum(u_r * patch) / np.sum(patch)))
        warm, cold = results
        assert warm < 0.0  # inward: toward the axis, "up"
        assert cold > 0.0
        assert warm == pytest.approx(-cold, rel=1e-2)  # antisymmetric to leading order


class TestAgainstTheSpinupModel:
    def test_axisymmetric_spin_change_matches_l1(self) -> None:
        cyl = Cylinder()
        delta = 0.1 * cyl.omega
        hull = cyl.omega + delta
        l1 = sudden_change(
            cyl,
            SpinupParams(1.0, n_cells=100),
            delta,
            t_max=72000.0,
            fraction=2.0,
        )
        lag = l1.relative_angular_momentum
        cfg = AnnulusConfig(
            omega=hull,
            n_theta=8,
            r_in_fraction=0.02,
            nu=1.0,
            forcing=Forcing(land_flux=0.0),
            n_r=100,
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = -delta * model.r
        l0 = model.angular_momentum(state)
        for hours in (5.0, 10.0, 20.0):
            model.advance(state, hours * 3600.0, dt_max=300.0)
            l1_value = float(np.interp(hours * 3600.0, l1.t, lag / lag[0]))
            assert model.angular_momentum(state) / l0 == pytest.approx(
                l1_value, rel=0.03
            )


class TestLambOseen:
    def test_mean_wind_follows_the_viscous_vortex_solution(self) -> None:
        nu, gamma = 1.0, 2.0e4
        cfg = AnnulusConfig(
            n_theta=8,
            n_r=400,
            stretch=0.0,
            r_in_fraction=0.005,
            nu=nu,
            endcap=False,
            floor_drag=0.0,
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)

        def vortex(t: float) -> np.ndarray:
            return np.asarray(
                gamma
                / (2.0 * math.pi * model.r)
                * (1.0 - np.exp(-(model.r**2) / (4.0 * nu * t)))
            )

        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = vortex(1000.0)
        for _ in range(50):
            model.step(state, 20.0)
        core = model.r < 400.0
        exact = vortex(2000.0)
        assert state.t == pytest.approx(1000.0)
        assert state.ubar[core] == pytest.approx(exact[core], rel=0.02)


class TestTransientHeating:
    def test_the_mean_wind_is_driven_only_while_the_heating_shape_changes(self) -> None:
        # Uniform wall heating of a closed fluid: the temperature first adjusts
        # (dT'/dt has a changing shape), then it just warms (dT/dt uniform, anomaly
        # zero). The Coriolis-buoyancy force follows dT'/dt, so the mean wind is
        # driven early and afterwards only decays at 1/tau_E. A wind still being
        # forced in the frozen state would be a bug.
        cfg = AnnulusConfig(
            n_r=16,
            n_theta=8,
            nu=5.0e4,
            advection=False,
            floor_drag=0.0,
            length=2.0e6,
            endcap=True,
            forcing=Forcing(land_flux=100.0, window="uniform"),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        early_force = []
        for _ in range(30):
            run_steps(model, state, 40.0, 1)
            early_force.append(
                np.max(np.abs(model.mean_wind_terms(state).coriolis_buoyancy))
            )
        assert max(early_force) > 0.0
        assert np.all(state.ubar > 0.0)  # spinward: the heated fluid expands inward
        run_steps(model, state, 40.0, 100)  # 4000 s: many diffusion times
        frozen = state.ubar.copy()
        assert np.max(
            np.abs(model.mean_wind_terms(state).coriolis_buoyancy)
        ) < 1e-3 * max(early_force)
        tau = 1.0 / cfg.endcap_rate
        run_steps(model, state, 40.0, round(0.5 * tau / 40.0))  # half a decay time
        assert state.ubar == pytest.approx(frozen * math.exp(-0.5), rel=1e-2)


class TestAdvanceAndDiagnostics:
    def test_advance_reaches_the_end_time_and_calls_back(self) -> None:
        model = Annulus(small(forcing=Forcing(land_flux=100.0)))
        state = model.initial_state(amplitude=0.01)
        calls = []
        model.advance(
            state, 200.0, dt_max=5.0, every=50.0, callback=lambda s: calls.append(s.t)
        )
        assert state.t == pytest.approx(200.0)
        assert len(calls) >= 3
        assert calls == sorted(calls)
        plain = model.initial_state(amplitude=0.01)
        model.advance(plain, 30.0, dt_max=5.0)
        assert plain.t == pytest.approx(30.0)

    def test_suggest_dt_sits_on_a_ladder_and_rises_slowly(self) -> None:
        model = Annulus(small(forcing=Forcing(land_flux=0.0)))
        state = model.initial_state(amplitude=0.0)
        assert model.suggest_dt(state, dt_max=8.0) == 8.0  # quiescent: the cap
        state.last_dt = 1.0
        assert model.suggest_dt(state, dt_max=8.0) == pytest.approx(1.25)
        state.last_dt = 8.0 / 1.25**6
        state.ubar[:] = 20.0
        small_dt = model.suggest_dt(state, dt_max=8.0)
        ladder = 8.0 / 1.25 ** np.arange(0, 40)
        assert np.min(np.abs(ladder - small_dt)) < 1e-9
        assert small_dt < 8.0

    def test_step_ratios_beyond_zero_stability_are_rejected(self) -> None:
        model = Annulus(small(forcing=Forcing(land_flux=0.0)))
        state = model.initial_state(amplitude=0.0)
        model.step(state, 1.0)
        model.step(state, 1.0)
        model.step(state, 2.0)  # ratio 2 is fine
        with pytest.raises(ValueError, match="step ratio"):
            model.step(state, 5.0)

    def test_default_cfl_is_conservative(self) -> None:
        assert DEFAULT_CFL == 0.15
        model = Annulus(small(forcing=Forcing(land_flux=0.0)))
        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = 5.0
        default = model.stable_dt(state, dt_max=1e9)
        assert default == pytest.approx(
            0.15 / 0.4 * model.stable_dt(state, cfl=0.4, dt_max=1e9)
        )

    def test_variable_steps_stay_accurate(self) -> None:
        # SBDF2 with a changing step reproduces exponential decay
        cfg = small(
            length=2000.0,
            floor_drag=0.0,
            advection=False,
            n_theta=8,
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = 0.001 * model.r
        tau = 1.0 / cfg.endcap_rate
        steps = [tau / 400.0] * 100 + [tau / 800.0] * 100 + [tau / 500.0] * 150
        for dt in steps:
            model.step(state, dt)
        expected = math.exp(-state.t / tau)
        assert state.ubar[3] / (0.001 * model.r[3]) == pytest.approx(expected, rel=2e-3)

    def test_solid_body_wind_has_uniform_vorticity(self) -> None:
        model = Annulus(small(forcing=Forcing(land_flux=0.0)))
        state = model.initial_state(amplitude=0.0)
        state.ubar[:] = 0.002 * model.r
        expected = (2.0 * 0.002) ** 2 * np.pi * np.sum(model.area)
        assert model.enstrophy(state) == pytest.approx(expected, rel=2e-2)
        assert model.angular_momentum(state) == pytest.approx(
            2.0 * math.pi * 0.002 * np.sum(model.moment * model.r)
        )

    def test_spectrum_and_wind_profile(self) -> None:
        model = Annulus(small(n_theta=48))
        state = model.initial_state(amplitude=0.0)
        shape = np.sin(math.pi * (model.r - model.r_in) / (3200.0 - model.r_in))
        state.zeta[:, 4] = 1e-3 * shape
        state.ubar[:] = 0.5
        spec = model.spectrum(state, 0.0, 3200.0)
        assert np.argmax(spec[1:]) + 1 == 4
        assert spec[0] == pytest.approx(math.pi * np.sum(model.area) * 0.25, rel=1e-9)
        wind = model.wind_at_height(state, 10.0)
        assert wind.shape == (model.n_theta,)
        assert np.all(np.isfinite(wind))
        with pytest.raises(ValueError, match="height"):
            model.wind_at_height(state, 0.0)

    def test_noslip_walls_and_polar_filter_off_run(self) -> None:
        cfg = small(wall="noslip", polar_filter=False, filter_strength=0.1)
        model = Annulus(cfg)
        state = model.initial_state(amplitude=0.01)
        for _ in range(20):
            model.step(state, 2.0)
        assert np.all(np.isfinite(state.zeta))
        assert model.energy(state) >= 0.0


GAP_RADIUS = 100.0
ETA = 0.99
GAP = GAP_RADIUS * (1.0 - ETA)
R_MEAN = GAP_RADIUS * (1.0 + ETA) / 2.0
G_EFF = 10.0


def secant_zero(f: object, lo: float, hi: float) -> float:
    """Zero of a nearly linear function from two evaluations."""
    f_lo, f_hi = f(lo), f(hi)  # type: ignore[operator]
    return float(lo - f_lo * (hi - lo) / (f_hi - f_lo))


def rigid_critical_ra(
    n_r: int, wall: str, k_crit: float, bracket: tuple[float, float]
) -> float:
    """Critical Ra of the thin annulus (eta = 0.99), fixed-T walls."""
    m = round(k_crit * R_MEAN / GAP)

    def growth(ra: float) -> float:
        cfg = AnnulusConfig(
            radius=GAP_RADIUS,
            r_in_fraction=ETA,
            omega=math.sqrt(G_EFF / R_MEAN),
            beta=1.0,
            nu=1.0,
            n_r=n_r,
            n_theta=2 * (3 * m // 2 + 2),
            stretch=0.0,
            wall=wall,
            outer_t="fixed",
            inner_t="fixed",
            t_outer=ra / (G_EFF * GAP**3),
            t_inner=0.0,
            floor_drag=0.0,
            endcap=False,
            coriolis_buoyancy=False,
            polar_filter=False,
            modes=(m,),
            forcing=Forcing(land_flux=0.0),
        )
        model = Annulus(cfg)
        base = model.initial_state(amplitude=0.0)
        base.temp[:, 0] = (
            cfg.t_outer
            * np.log(model.r / model.r_in)
            / math.log(GAP_RADIUS / model.r_in)
        )
        return model.max_growth_rate(base, m)

    return secant_zero(growth, *bracket)


def flux_critical_ra(n_r: int) -> float:
    """Critical flux Rayleigh number of the thin annulus, rigid walls."""

    def growth(ra: float) -> float:
        gradient = ra / (G_EFF * GAP**4)  # nu = kappa = beta = 1
        q = gradient * R_MEAN / GAP_RADIUS
        cfg = AnnulusConfig(
            radius=GAP_RADIUS,
            r_in_fraction=ETA,
            omega=math.sqrt(G_EFF / R_MEAN),
            beta=1.0,
            nu=1.0,
            rho_cp=1.0,
            n_r=n_r,
            n_theta=16,
            stretch=0.0,
            wall="noslip",
            outer_t="flux",
            inner_t="flux",
            flux_inner=q * GAP_RADIUS / (ETA * GAP_RADIUS),
            floor_drag=0.0,
            endcap=False,
            coriolis_buoyancy=False,
            polar_filter=False,
            modes=(1,),
            forcing=Forcing(land_flux=q, window="uniform"),
        )
        model = Annulus(cfg)
        base = model.initial_state(amplitude=0.0)
        base.temp[:, 0] = gradient * R_MEAN * np.log(model.r / model.r_in)
        return model.max_growth_rate(base, 1)

    return secant_zero(growth, 650.0, 800.0)


class TestConvectionOnset:
    """Raw thresholds converge at second order, with no extrapolation."""

    def test_free_slip_onset_is_657_5(self) -> None:
        coarse = rigid_critical_ra(24, "free", 2.221, (600.0, 720.0))
        fine = rigid_critical_ra(48, "free", 2.221, (600.0, 720.0))
        assert fine == pytest.approx(657.5, rel=0.003)  # + 0.2 % from curvature
        assert abs(fine - 658.9) < 0.4 * abs(coarse - 658.9)

    def test_rigid_fixed_temperature_onset_is_1708(self) -> None:
        coarse = rigid_critical_ra(24, "noslip", 3.117, (1650.0, 1770.0))
        fine = rigid_critical_ra(48, "noslip", 3.117, (1650.0, 1770.0))
        assert fine == pytest.approx(1708.0, rel=0.003)
        assert (1708.0 - fine) < 0.4 * (1708.0 - coarse)  # second order: ~4x

    def test_rigid_fixed_flux_onset_is_720(self) -> None:
        coarse, fine = flux_critical_ra(12), flux_critical_ra(24)
        assert fine == pytest.approx(720.0, rel=0.003)
        assert (720.0 - fine) < 0.4 * (720.0 - coarse)

    def test_the_stable_side_of_the_threshold_decays(self) -> None:
        # with Ra well below critical every mode is damped
        cfg = small(
            forcing=Forcing(land_flux=0.0),
            outer_t="fixed",
            inner_t="fixed",
            t_outer=1e-3,
        )
        model = Annulus(cfg)
        base = model.initial_state(amplitude=0.0)
        assert model.max_growth_rate(base, 3) < 0.0


class TestWallVorticityClosure:
    """Exact streamfunction -> discrete Poisson -> wall vorticity converges."""

    @staticmethod
    def wall_error(n_r: int, wall: str) -> float:
        cfg = AnnulusConfig(
            radius=1.0,
            r_in_fraction=0.3,
            n_r=n_r,
            n_theta=16,
            stretch=0.0,
            wall=wall,
            floor_drag=0.0,
        )
        model = Annulus(cfg)
        m = 2
        length = 1.0 - model.r_in
        s = (model.r - model.r_in) / length
        if wall == "noslip":  # psi = psi_r = 0 on both walls
            psi = (s * (1.0 - s)) ** 2
            dpsi = 2.0 * s * (1.0 - s) * (1.0 - 2.0 * s) / length
            d2psi = (2.0 - 12.0 * s * (1.0 - s)) / length**2
            exact_in = -2.0 / length**2
            exact_out = -2.0 / length**2
        else:  # free slip needs only psi = 0: exact zeta_w = -2 psi_r / r
            psi = np.sin(math.pi * s)
            dpsi = math.pi * np.cos(math.pi * s) / length
            d2psi = -((math.pi / length) ** 2) * psi
            exact_in = -2.0 * (math.pi / length) / model.r_in
            exact_out = -2.0 * (-math.pi / length) / 1.0
        zeta = -(d2psi + dpsi / model.r - m**2 * psi / model.r**2)
        z = np.zeros((n_r, model.k), dtype=complex)
        z[:, m] = zeta
        z_in, z_out = model._zeta_walls(model.psi(z))  # noqa: SLF001
        return float(
            max(
                abs(z_in[m] - exact_in) / abs(exact_in),
                abs(z_out[m] - exact_out) / abs(exact_out),
            )
        )

    def test_free_slip_wall_vorticity_converges_at_second_order(self) -> None:
        errors = [self.wall_error(n, "free") for n in (12, 24, 48)]
        assert errors[1] < 0.3 * errors[0]
        assert errors[2] < 0.3 * errors[1]
        assert errors[2] < 5e-4

    def test_no_slip_wall_vorticity_converges_at_first_order(self) -> None:
        # the cell-centred psi error near a rigid wall is not smooth enough for the
        # Taylor formula to be second order; the onset thresholds still are
        errors = [self.wall_error(n, "noslip") for n in (24, 48, 96)]
        assert errors[1] < 0.65 * errors[0]
        assert errors[2] < 0.65 * errors[1]
        assert errors[2] < 0.05
