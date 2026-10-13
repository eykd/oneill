import numpy as np
import pytest

from oneill.core.cylinder import G_STANDARD, Cylinder
from oneill.motion.figures import FigureSpec, Panel, Series, build_figures

SPEEDS = (1.4, 6.0)
RATES = (0.3, 3.0)


@pytest.fixture
def figures(default_cyl: Cylinder) -> dict[str, FigureSpec]:
    specs = build_figures(default_cyl, SPEEDS, RATES, 3000.0, 10.0, 1.0)
    return {s.filename: s for s in specs}


def test_four_named_figures_in_order(default_cyl: Cylinder) -> None:
    specs = build_figures(default_cyl, SPEEDS, RATES, 3000.0, 10.0, 1.0)
    assert [s.filename for s in specs] == [
        "1_floor_weight_vs_heading.png",
        "2_climbing_vs_height.png",
        "3_elevator.png",
        "4_ballistics.png",
    ]
    assert [len(s.panels) for s in specs] == [1, 2, 3, 2]
    assert [s.layout for s in specs] == [(1, 1), (1, 2), (3, 1), (1, 2)]


def test_every_series_is_finite_and_well_shaped(figures: dict[str, FigureSpec]) -> None:
    for spec in figures.values():
        for panel in spec.panels:
            assert panel.series
            for series in panel.series:
                assert series.x.shape == series.y.shape
                assert series.x.ndim == 1
                assert np.all(np.isfinite(series.x))
                assert np.all(np.isfinite(series.y))


class TestFloorWeight:
    def test_one_curve_per_speed_over_a_full_circle_of_headings(
        self, figures: dict[str, FigureSpec]
    ) -> None:
        panel = figures["1_floor_weight_vs_heading.png"].panels[0]
        assert [s.label for s in panel.series] == ["1.4 m/s", "6 m/s"]
        for series in panel.series:
            assert series.x[0] == 0.0
            assert series.x[-1] == 360.0
            assert len(series.x) == 73
        assert panel.hlines == (100.0,)
        assert panel.xticklabels == (
            "+axial",
            "spinward",
            "-axial",
            "antispinward",
            "+axial",
        )
        assert panel.xticks == (0.0, 90.0, 180.0, 270.0, 360.0)

    def test_curve_follows_the_inertial_frame_weight_formula(
        self, default_cyl: Cylinder, figures: dict[str, FigureSpec]
    ) -> None:
        series = figures["1_floor_weight_vs_heading.png"].panels[0].series[1]
        v = 6.0
        expected = (
            100.0
            * (default_cyl.rim_speed + v * np.sin(np.radians(series.x))) ** 2
            / default_cyl.radius
            / G_STANDARD
        )
        np.testing.assert_allclose(series.y, expected, rtol=1e-9)
        # heavier spinward, lighter antispinward, unchanged along the axis
        assert series.y[18] > 100.0 > series.y[54]
        assert series.y[0] == pytest.approx(100.0, rel=1e-9)


class TestClimbing:
    def test_local_gravity_curve_and_per_rate_curves(
        self, default_cyl: Cylinder, figures: dict[str, FigureSpec]
    ) -> None:
        spec = figures["2_climbing_vs_height.png"]
        push, tilt = spec.panels
        assert [s.label for s in push.series] == [
            "local gravity (centrifugal)",
            "spinward push @ 0.3 m/s up",
            "spinward push @ 3 m/s up",
        ]
        assert [s.label for s in tilt.series] == [
            "climbing 0.3 m/s",
            "climbing 3 m/s",
        ]
        gravity = push.series[0]
        assert gravity.fmt == "k--"
        assert gravity.y[0] == pytest.approx(100.0)
        assert gravity.x[-1] == pytest.approx(0.99 * default_cyl.radius)
        assert gravity.y[-1] == pytest.approx(1.0, rel=1e-9)

    def test_push_is_constant_and_tilt_grows_with_height(
        self, default_cyl: Cylinder, figures: dict[str, FigureSpec]
    ) -> None:
        push, tilt = figures["2_climbing_vs_height.png"].panels
        expected_push = 100.0 * 2.0 * default_cyl.omega * 3.0 / G_STANDARD
        np.testing.assert_allclose(push.series[2].y, expected_push, rtol=1e-6)
        assert np.all(np.diff(tilt.series[1].y) > 0.0)
        assert tilt.series[1].y[-1] > 45.0


class TestElevator:
    def test_up_and_down_curves_in_three_stacked_panels(
        self, figures: dict[str, FigureSpec]
    ) -> None:
        spec = figures["3_elevator.png"]
        heights, forces, tilt = spec.panels
        assert spec.sharex is True
        assert [s.label for s in heights.series] == ["up", "down"]
        assert [s.label for s in forces.series] == [
            "up: felt weight",
            "up: spinward push",
            "down: felt weight",
            "down: spinward push",
        ]
        assert [s.label for s in tilt.series] == ["up", "down"]
        up, down = heights.series
        assert up.y[0] == pytest.approx(0.0, abs=1e-9)
        assert up.y[-1] == pytest.approx(3000.0, abs=1e-6)
        assert down.y[0] == pytest.approx(3000.0, abs=1e-6)
        assert up.x[-1] == pytest.approx(310.0)
        assert forces.series[1].fmt == "--"
        assert tilt.xlabel == "time [s]"

    def test_push_is_positive_going_up_and_negative_going_down(
        self, figures: dict[str, FigureSpec]
    ) -> None:
        forces = figures["3_elevator.png"].panels[1]
        assert forces.series[1].y.max() > 0.0
        assert forces.series[3].y.min() < 0.0


class TestBallistics:
    def test_cases_split_between_panels_by_launch_height(
        self, figures: dict[str, FigureSpec]
    ) -> None:
        low, high = figures["4_ballistics.png"].panels
        assert "drop from hand (1.5 m)" in [s.label for s in low.series]
        assert [s.label for s in high.series] == [
            "drop from 100 m tower",
            "drop from 1000 m",
        ]
        for series in high.series:
            assert series.y[0] == pytest.approx(series.y.max())
            assert series.y[-1] == pytest.approx(0.0, abs=1e-6)

    def test_cases_that_do_not_fit_are_omitted(self) -> None:
        small = Cylinder.from_gravity(radius=500.0, length=2000.0)
        spec = build_figures(small, SPEEDS, RATES, 400.0, 10.0, 1.0)[3]
        labels = [s.label for p in spec.panels for s in p.series]
        assert "drop from 1000 m" not in labels
        assert "drop from 100 m tower" in labels


def test_elevator_shaft_must_fit_inside_the_habitat(default_cyl: Cylinder) -> None:
    with pytest.raises(ValueError, match="height"):
        build_figures(default_cyl, SPEEDS, RATES, 4000.0, 10.0, 1.0)


def test_value_objects_default_to_plain_styles() -> None:
    series = Series("a", np.zeros(2), np.zeros(2))
    assert series.fmt == ""
    panel = Panel("x", "y", (series,))
    assert (panel.hlines, panel.xticks, panel.xticklabels) == ((), (), ())
