import math

import numpy as np
import pytest

from oneill.atmo.equilibrium import isothermal_axis_ratio, lifting_condensation_level
from oneill.atmo.figures import atmosphere_figures
from oneill.atmo.gas import ONEILL_1974
from oneill.atmo.report import (
    CAVEATS,
    FLOOR_LABEL,
    FLOOR_TEMPERATURE,
    AtmosphereReport,
    build_report,
)
from oneill.core.cylinder import Cylinder


@pytest.fixture(scope="module")
def report() -> AtmosphereReport:
    return build_report()


def test_build_report_is_cached_per_habitat(report: AtmosphereReport) -> None:
    assert build_report() is report
    assert build_report(Cylinder()) is report


def test_caveats_cover_every_section(report: AtmosphereReport) -> None:
    assert set(report.caveats) == {
        "gas",
        "equilibrium",
        "numbers",
        "surface",
        "parcel",
        "spinup",
    }
    assert report.caveats == CAVEATS
    assert report.caveats is not CAVEATS  # a copy: callers cannot corrupt it


def test_gas_rows_are_labelled_with_preset_and_temperature(
    report: AtmosphereReport,
) -> None:
    assert [g.label for g in report.gas] == [
        "O'Neill 1974 at 290.0 K",
        "SP-413 at 290.0 K",
    ]
    assert report.gas[0].r_specific == pytest.approx(287.05, rel=1e-4)
    assert report.gas[1].p_dry == 49.7e3
    assert report.gas[1].p_vapour == 1.0e3
    assert all(g.temperature == FLOOR_TEMPERATURE for g in report.gas)


def test_equilibrium_matches_the_research_numbers(report: AtmosphereReport) -> None:
    oneill = report.equilibrium[0]
    assert oneill.axis_ratio_isothermal == isothermal_axis_ratio(
        Cylinder(), ONEILL_1974.at_temperature(290.0)
    )
    assert oneill.axis_ratio_isothermal == pytest.approx(0.83, abs=0.005)
    assert oneill.delta_t_adiabatic == pytest.approx(15.6, abs=0.1)
    assert oneill.floor_lapse_rate * 1e3 == pytest.approx(9.76, abs=0.02)
    assert [c.label for c in oneill.clouds] == [oneill.label, oneill.label]
    expected = lifting_condensation_level(
        Cylinder(), ONEILL_1974.at_temperature(290.0), 0.5
    )
    assert oneill.clouds[0].height == expected.height


def test_numbers_rows(report: AtmosphereReport) -> None:
    assert len(report.numbers) == 2
    assert report.numbers[0].coriolis_vs_earth == pytest.approx(759.0, abs=1.0)
    assert [s.nu for s in report.numbers[0].table.eddy] == [1.0, 10.0]


def test_surface_rows_cover_each_terrain_and_preset(report: AtmosphereReport) -> None:
    assert len(report.surface) == 6
    suburb = report.surface[0]
    assert suburb.layer == pytest.approx((4.64, 12.97), abs=0.02)
    assert suburb.q_mean == pytest.approx(2.0 * report.surface[1].q_mean, rel=0.03)
    assert len(report.stability) == 9
    spinward = report.stability[0]
    assert spinward.heading == "spinward"
    assert [c.name for c in spinward.classes][:3] == ["STABLE", "STABLE", "UNSTABLE"]
    assert {c.name for c in report.stability[2].classes} == {"NEUTRAL"}


def test_parcel_rows(report: AtmosphereReport) -> None:
    assert len(report.parcels) == 8
    free = report.parcels[0]
    assert free.drag == 0.0
    assert free.steady_rise == 0.0
    assert free.time_to_lcl is None
    assert free.max_rise == pytest.approx(5.5, abs=0.1)
    damped = report.parcels[3]
    assert damped.time_to_lcl == pytest.approx(17758.0, rel=0.01)
    assert report.parcels[7].time_to_lcl == pytest.approx(5925.0, rel=0.01)
    assert len(damped.t) == len(damped.height) == len(damped.spinward)


def test_spin_rows_state_the_floor_condition_and_validity(
    report: AtmosphereReport,
) -> None:
    assert len(report.spin) == 8
    assert {r.floor for r in report.spin} == {FLOOR_LABEL}
    assert {r.floor for r in report.lag} == {FLOOR_LABEL}
    assert "C_D = 0.0178" in FLOOR_LABEL
    row = next(
        r for r in report.spin if (r.ramp, r.nu, r.endcap) == (3600.0, 1.0, True)
    )
    assert row.time_to_90 is not None
    assert row.time_to_90 / 3600.0 == pytest.approx(23.0, abs=0.5)
    assert row.peak_wind_ref == pytest.approx(18.3, abs=1.0)
    assert row.max_rossby > 0.3
    assert not row.sink_valid
    gentle = next(
        r for r in report.spin if (r.ramp, r.nu, r.endcap) == (86400.0, 10.0, True)
    )
    assert gentle.sink_valid


def test_sudden_changes_are_nearly_symmetric_and_spin_down_is_slower(
    report: AtmosphereReport,
) -> None:
    down, up = report.sudden[0], report.sudden[1]
    assert down.delta_fraction < 0.0 < up.delta_fraction
    assert down.time_to_adjust is not None
    assert up.time_to_adjust is not None
    # with end caps tau_E ~ 1 / sqrt(new rate): spin-down (0.9) is the slower
    assert down.time_to_adjust / up.time_to_adjust == pytest.approx(1.08, abs=0.04)


def test_blend_width_sensitivity_is_reported(report: AtmosphereReport) -> None:
    assert [b.delta_fraction for b in report.blend] == [-0.1, 0.1, -0.5]
    for row in report.blend:
        assert row.widths == (0.15, 0.3, 0.6)
        assert all(t is not None for t in row.times)
    mild, strong = report.blend[0], report.blend[2]
    assert mild.times[0] is not None
    assert mild.times[2] is not None
    assert strong.times[0] is not None
    assert strong.times[2] is not None
    assert abs(strong.times[2] / strong.times[0] - 1.0) > abs(
        mild.times[2] / mild.times[0] - 1.0
    )


def test_reference_run_data(report: AtmosphereReport) -> None:
    assert len(report.snapshots) == 5
    assert report.corotation[0] == 0.0
    assert report.corotation[-1] >= 0.9
    assert len(report.corotation) == len(report.corotation_hours)
    assert [m.label for m in report.molecular] == [g.label for g in report.gas]


def test_every_number_is_finite(report: AtmosphereReport) -> None:
    for row in report.spin:
        assert math.isfinite(row.peak_wind_ref)
    for lag in report.lag:
        assert lag.peak > 0.0


def test_figures(report: AtmosphereReport) -> None:
    gas = ONEILL_1974.at_temperature(290.0)
    figures = atmosphere_figures(report, gas)
    assert [f.filename for f in figures] == [
        "atmo_1_profiles.png",
        "atmo_2_wind.png",
        "atmo_3_parcel.png",
        "atmo_4_spinup.png",
    ]
    for fig in figures:
        assert len(fig.panels) == 2
        for panel in fig.panels:
            assert panel.series
            for s in panel.series:
                assert np.all(np.isfinite(s.x))
                assert np.all(np.isfinite(s.y))
    wind = figures[1].panels[0]
    assert len(wind.bands) == 1
    assert wind.bands[0].lo < wind.bands[0].hi
