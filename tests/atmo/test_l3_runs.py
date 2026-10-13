import json
import math

import numpy as np
import pytest

from oneill.atmo.annulus2d import Annulus, AnnulusConfig, Forcing
from oneill.atmo.annulus_diagnostics import (
    large_scale_fraction,
    lock_fraction,
    spectrum_scale,
    stationarity,
    take_sample,
    wind_statistics,
)
from oneill.atmo.l3_runs import (
    CASES,
    DAY,
    CaseResult,
    case_config,
    case_figures,
    run_case,
    summarize,
)
from oneill.core.cylinder import Cylinder


class TestDiagnostics:
    def model(self) -> Annulus:
        cfg = AnnulusConfig(
            n_r=16, n_theta=32, forcing=Forcing(land_flux=200.0), nu=3.0
        )
        return Annulus(cfg)

    def test_quiescent_sample_is_zero(self) -> None:
        model = self.model()
        state = model.initial_state(amplitude=0.0)
        s = take_sample(model, state)
        assert s.energy == 0.0
        assert s.mean_energy == 0.0
        assert s.nusselt == 0.0
        assert s.dissipation == 0.0
        assert s.spectrum.shape == (model.k,)
        assert np.all(s.wind == 0.0)
        assert s.ubar.shape == (model.config.n_r,)

    def test_a_convecting_sample_has_energy_heat_exchange_and_dissipation(self) -> None:
        model = self.model()
        state = model.initial_state(amplitude=0.05, seed=1)
        model.advance(state, 900.0, dt_max=10.0)
        s = take_sample(model, state)
        assert s.energy > 0.0
        assert s.enstrophy > 0.0
        assert s.nusselt > 0.0
        assert s.dissipation == pytest.approx(
            model.config.nu * 2.0 * s.enstrophy / (2.0 * math.pi * np.sum(model.area)),
            rel=1e-9,
        )
        assert np.sum(s.spectrum) == pytest.approx(s.energy, rel=0.5)

    def test_spectrum_scale(self) -> None:
        spec = np.zeros(10)
        spec[3] = 4.0
        spec[5] = 4.0
        scale = spectrum_scale(spec, 1000.0)
        assert scale.m_peak == 3
        assert scale.m_centroid == pytest.approx(4.0)
        assert scale.length_peak == pytest.approx(2.0 * math.pi * 1000.0 / 3.0)
        assert scale.length_centroid == pytest.approx(2.0 * math.pi * 1000.0 / 4.0)
        empty = spectrum_scale(np.zeros(5), 1000.0)
        assert empty.m_peak == 1
        assert empty.m_centroid == 1.0

    def test_lock_and_large_scale_fractions(self) -> None:
        spec = np.array([9.0, 1.0, 1.0, 6.0, 1.0, 1.0])
        assert lock_fraction(spec, 3) == pytest.approx(6.0 / 10.0)
        assert lock_fraction(spec, 1) == pytest.approx(0.1)
        assert large_scale_fraction(spec, 2) == pytest.approx(0.2)
        silent = np.array([1.0, 0.0, 0.0, 0.0])
        assert lock_fraction(silent) == 0.0
        assert large_scale_fraction(silent) == 0.0

    def test_wind_statistics(self) -> None:
        rms, gust = wind_statistics(np.array([3.0, -3.0, 3.0, -3.0]))
        assert rms == pytest.approx(3.0)
        assert gust == pytest.approx(3.0 + 3.0 * 3.0)
        rms0, gust0 = wind_statistics(np.full(8, 2.0))
        assert rms0 == pytest.approx(2.0)
        assert gust0 == pytest.approx(2.0)

    def test_stationarity_detector(self) -> None:
        t = np.arange(40.0)
        assert stationarity(t, np.full(40, 5.0)).stationary
        flat = 5.0 + 0.1 * np.sin(t)
        assert stationarity(t, flat).stationary
        growing = stationarity(t, 1.0 + t)
        assert not growing.stationary
        assert growing.drift > 0.1
        decaying = stationarity(t, 40.0 - t)
        assert decaying.drift < 0.0
        assert not stationarity(t[:3], np.ones(3)).stationary  # too short
        assert not stationarity(t[:1], np.ones(1)).stationary


class TestCases:
    def test_the_production_set(self) -> None:
        assert set(CASES) == {
            "primary",
            "steady",
            "cb_off",
            "nu10",
            "rin005",
            "rin02",
            "uniform_sink",
        }
        cyl = Cylinder()
        primary = case_config(CASES["primary"], cyl)
        assert primary.nu == 1.0
        assert primary.coriolis_buoyancy
        assert primary.forcing.period == DAY
        assert primary.forcing.window == "neutral"
        assert primary.tau_relax is None
        assert primary.rho_cp == pytest.approx(
            603.0, rel=0.03
        )  # the half-atmosphere air
        assert not case_config(CASES["cb_off"], cyl).coriolis_buoyancy
        assert case_config(CASES["nu10"], cyl).nu == 10.0
        assert case_config(CASES["rin005"], cyl).r_in_fraction == 0.05
        assert case_config(CASES["rin02"], cyl).r_in_fraction == 0.2
        steady = case_config(CASES["steady"], cyl)
        assert steady.forcing.period is None
        sink = case_config(CASES["uniform_sink"], cyl)
        assert sink.forcing.window == "zero"
        assert sink.tau_relax == DAY


@pytest.fixture(scope="module")
def tiny() -> CaseResult:
    return run_case(
        CASES["primary"], 0.04, n_r=16, n_theta=32, sample_every=300.0, dt_max=10.0
    )


class TestRun:
    def test_samples_cover_the_run(self, tiny: CaseResult) -> None:
        times = [s.t for s in tiny.samples]
        assert times[0] == 0.0
        assert times[-1] == pytest.approx(0.04 * DAY, abs=10.0)
        assert times == sorted(times)
        assert len(times) >= 10
        assert tiny.steps > 100
        assert tiny.zeta.shape == (16, 32)
        assert tiny.temp.shape == (16, 32)

    def test_progress_callback(self) -> None:
        seen: list[float] = []
        run_case(
            CASES["steady"],
            0.002,
            n_r=16,
            n_theta=32,
            sample_every=100.0,
            progress=lambda t, _e: seen.append(t),
        )
        assert seen[0] == 0.0
        assert len(seen) >= 2

    def test_a_run_ending_on_a_sample_is_not_sampled_twice(self) -> None:
        end = 0.002 * DAY
        result = run_case(
            CASES["steady"], 0.002, n_r=16, n_theta=32, sample_every=end / 2.0
        )
        times = [s.t for s in result.samples]
        assert times[-1] == pytest.approx(end)
        assert len(set(np.round(times, 6))) == len(times)

    def test_summary_is_json_ready_and_complete(self, tiny: CaseResult) -> None:
        summary = summarize(tiny)
        text = json.dumps(summary)  # no numpy types survive
        assert json.loads(text) == summary
        for key in (
            "days",
            "steps",
            "stationary",
            "energy_drift",
            "m_peak",
            "length_peak_m",
            "lock_fraction",
            "large_scale_fraction",
            "rms_wind_10m",
            "gust_10m",
            "zonal_wind_10m",
            "nusselt",
            "dissipation",
            "arrest_low_m",
            "arrest_high_m",
            "mean_wind_profile",
            "energy_series",
        ):
            assert key in summary
        assert summary["name"] == "primary"
        assert summary["arrest_low_m"] < summary["arrest_high_m"]
        assert summary["arrest_high_m"] == pytest.approx(
            summary["arrest_low_m"] * 0.18 / 0.08, rel=1e-4
        )
        assert 1 <= len(summary["energy_series"]) <= 60
        assert (
            len(summary["mean_wind_profile"]) == 16
            or len(summary["mean_wind_profile"]) <= 24
        )
        assert summary["rms_wind_10m"] >= 0.0
        assert 0.0 <= summary["lock_fraction"] <= 1.0

    def test_quiet_run_has_no_arrest_scale(self) -> None:
        quiet = run_case(
            CASES["steady"], 0.0005, n_r=16, n_theta=32, sample_every=20.0, seed=0
        )
        summary = summarize(quiet)
        assert summary["mean_wind_energy_fraction"] >= 0.0
        assert summary["arrest_low_m"] >= 0.0

    def test_figures_are_finite_data_specs(self, tiny: CaseResult) -> None:
        figures = case_figures(tiny)
        assert [f.filename for f in figures] == [
            "l3_primary_snapshot.png",
            "l3_primary_spectrum.png",
            "l3_primary_meanwind.png",
        ]
        snapshot, spectrum, meanwind = figures
        assert len(snapshot.panels) == 2
        for panel in snapshot.panels:
            mesh = panel.meshes[0]
            assert mesh.x.shape == mesh.values.shape == (16, 33)
            assert np.all(np.isfinite(mesh.values))
            assert panel.equal_aspect
        # the polar mesh closes: first and last columns coincide
        mesh = snapshot.panels[0].meshes[0]
        assert mesh.x[:, 0] == pytest.approx(mesh.x[:, -1])
        spec_mesh = spectrum.panels[0].meshes[0]
        assert spec_mesh.values.shape == spec_mesh.x.shape
        assert np.all(np.isfinite(spec_mesh.values))
        assert len(meanwind.panels[0].series) == 2
