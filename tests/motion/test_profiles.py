import math

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.motion.profiles import ConstantSpeed, Trapezoid

RIDES = [
    (3000.0, 10.0, 1.0),
    (50.0, 10.0, 1.0),
    (100.0, 10.0, 1.0),
    (1.0, 100.0, 0.5),
    (93.75, 10.0, 1.0),
    (0.001, 5.0, 2.0),
]


class TestConstantSpeed:
    def test_is_linear_for_all_times(self) -> None:
        prof = ConstantSpeed(2.5)
        assert prof.s(0.0) == 0.0
        assert prof.s(4.0) == 10.0
        assert prof.s(-0.01) == pytest.approx(-0.025)

    def test_derivatives(self) -> None:
        prof = ConstantSpeed(-2.5)
        np.testing.assert_array_equal(prof.ds(np.array([-1.0, 0.0, 9.0])), -2.5)
        np.testing.assert_array_equal(prof.dds(np.array([-1.0, 0.0, 9.0])), 0.0)
        assert prof.ds(3.0).shape == ()

    def test_array_input_keeps_shape(self) -> None:
        ts = np.linspace(0.0, 1.0, 5)
        assert ConstantSpeed(1.0).s(ts).shape == (5,)

    @pytest.mark.parametrize("bad", [math.nan, math.inf])
    def test_rejects_non_finite_speed(self, bad: float) -> None:
        with pytest.raises(ValueError, match="speed"):
            ConstantSpeed(bad)


class TestTrapezoid:
    @pytest.mark.parametrize(
        ("distance", "v_max", "a_max", "expected_duration"),
        [
            (3000.0, 10.0, 1.0, 310.0),
            (50.0, 10.0, 1.0, 2.0 * math.sqrt(50.0)),
            (100.0, 10.0, 1.0, 20.0),
            (1.0, 100.0, 0.5, 2.0 * math.sqrt(2.0)),
        ],
    )
    def test_duration(
        self, distance: float, v_max: float, a_max: float, expected_duration: float
    ) -> None:
        assert Trapezoid(distance, v_max, a_max).duration == pytest.approx(
            expected_duration, rel=1e-12
        )

    def test_known_waypoints_of_the_reference_ride(self) -> None:
        prof = Trapezoid(3000.0, 10.0, 1.0)
        assert prof.s(10.0) == pytest.approx(50.0)
        assert prof.s(155.0) == pytest.approx(1500.0)
        assert prof.s(300.0) == pytest.approx(2950.0)
        assert prof.s(prof.duration) == pytest.approx(3000.0)

    def test_phase_derivatives_of_the_reference_ride(self) -> None:
        prof = Trapezoid(3000.0, 10.0, 1.0)
        t = np.array([-1.0, 0.0, 5.0, 10.0, 150.0, 300.0, 305.0, 310.0, 320.0])
        np.testing.assert_allclose(prof.ds(t), [0, 0, 5, 10, 10, 10, 5, 0, 0])
        np.testing.assert_allclose(prof.dds(t), [0, 1, 1, 0, 0, -1, -1, 0, 0])

    @pytest.mark.parametrize(("distance", "v_max", "a_max"), RIDES)
    def test_profile_respects_its_limits(
        self, distance: float, v_max: float, a_max: float
    ) -> None:
        prof = Trapezoid(distance, v_max, a_max)
        assert prof.s(0.0) == 0.0
        assert prof.s(-5.0) == 0.0
        assert prof.s(prof.duration) == pytest.approx(distance, rel=1e-12)
        assert prof.s(prof.duration + 5.0) == distance

        ts = np.linspace(0.0, prof.duration, 2001)
        speeds = prof.ds(ts)
        assert speeds.min() >= -1e-12
        assert speeds.max() <= min(v_max, math.sqrt(distance * a_max)) * (1.0 + 1e-9)
        assert np.abs(prof.dds(ts)).max() <= a_max * (1.0 + 1e-12)
        np.testing.assert_allclose(
            prof.s(ts) + prof.s(ts[::-1]), distance, rtol=1e-9, atol=1e-12
        )

    @given(st.floats(0.01, 1e4), st.floats(0.1, 100.0), st.floats(0.05, 20.0))
    @example(100.0, 10.0, 1.0)
    @example(50.0, 10.0, 1.0)
    @example(3000.0, 10.0, 1.0)
    def test_derivatives_are_the_derivatives_of_the_position(
        self, distance: float, v_max: float, a_max: float
    ) -> None:
        prof = Trapezoid(distance, v_max, a_max)
        t_acc = min(v_max / a_max, math.sqrt(distance / a_max))
        h = 1e-3 * t_acc
        # Probe the middle of each phase so no kink falls inside the stencil.
        probes = [0.5 * t_acc, prof.duration - 0.5 * t_acc]
        if prof.t_cruise > 4.0 * h:
            probes.append(t_acc + 0.5 * prof.t_cruise)
        # Finite differences amplify float round-off in ds ~ a*duration by 1/h.
        slack = 1e-5 + 1e-15 * a_max * prof.duration / h
        for t in probes:
            num_v = (prof.s(t + h) - prof.s(t - h)) / (2.0 * h)
            num_a = (prof.ds(t + h) - prof.ds(t - h)) / (2.0 * h)
            assert num_v == pytest.approx(prof.ds(t), rel=1e-6, abs=1e-6)
            assert num_a == pytest.approx(prof.dds(t), abs=slack)

    @pytest.mark.parametrize(
        ("distance", "v_max", "a_max"),
        [
            (0.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0),
            (1.0, 0.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, 0.0),
            (1.0, 1.0, -1.0),
            (math.nan, 1.0, 1.0),
        ],
    )
    def test_rejects_non_positive_parameters(
        self, distance: float, v_max: float, a_max: float
    ) -> None:
        with pytest.raises(ValueError, match="positive"):
            Trapezoid(distance, v_max, a_max)
