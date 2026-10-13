import math
from typing import TYPE_CHECKING

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from oneill.core.integrate import rk4

if TYPE_CHECKING:
    from collections.abc import Callable

    from oneill.core.vectors import Vec3


def decay(rate: float) -> Callable[[float, Vec3], Vec3]:
    def deriv(_t: float, y: Vec3) -> Vec3:
        return -rate * y

    return deriv


def oscillator(_t: float, y: Vec3) -> Vec3:
    return np.array([y[1], -y[0]])


@settings(max_examples=15)
@given(st.floats(0.1, 3.0), st.floats(0.5, 3.0), st.floats(-5.0, 5.0))
def test_exponential_decay_matches_the_closed_form(
    rate: float, t_end: float, y0: float
) -> None:
    y = rk4(decay(rate), 0.0, np.array([y0]), t_end, max_step=1e-2)
    assert y[0] == pytest.approx(y0 * math.exp(-rate * t_end), rel=1e-7, abs=1e-9)


def test_harmonic_oscillator_returns_after_a_period() -> None:
    y = rk4(oscillator, 0.0, np.array([1.0, 0.0]), 2.0 * math.pi, max_step=1e-2)
    np.testing.assert_allclose(y, [1.0, 0.0], atol=1e-8)


def test_error_shrinks_with_the_fourth_power_of_the_step() -> None:
    exact = math.exp(-2.0)
    errors = [
        abs(rk4(decay(1.0), 0.0, np.array([1.0]), 2.0, max_step=h)[0] - exact)
        for h in (0.2, 0.1)
    ]
    assert errors[0] / errors[1] == pytest.approx(16.0, rel=0.2)


def test_time_dependent_derivative_sees_the_current_time() -> None:
    # y' = 2t  ->  y = t^2, which RK4 integrates exactly.
    y = rk4(lambda t, _y: np.array([2.0 * t]), 1.0, np.array([1.0]), 4.0, max_step=0.7)
    assert y[0] == pytest.approx(16.0, rel=1e-12)


def test_backward_integration_undoes_forward_integration() -> None:
    y0 = np.array([1.0, 0.5])
    forward = rk4(oscillator, 0.0, y0, 3.0, max_step=1e-2)
    back = rk4(oscillator, 3.0, forward, 0.0, max_step=1e-2)
    np.testing.assert_allclose(back, y0, atol=1e-9)


def test_zero_length_interval_returns_an_untouched_copy() -> None:
    y0 = np.array([1, 2, 3])
    y = rk4(oscillator, 5.0, y0, 5.0, max_step=0.1)
    np.testing.assert_array_equal(y, [1.0, 2.0, 3.0])
    assert y is not y0
    assert y0.dtype.kind == "i"


def test_step_count_never_exceeds_max_step() -> None:
    calls: list[float] = []

    def deriv(t: float, y: Vec3) -> Vec3:
        calls.append(t)
        return np.zeros_like(y)

    rk4(deriv, 0.0, np.array([0.0]), 1.0, max_step=0.3)
    assert len(calls) == 4 * 4  # ceil(1 / 0.3) = 4 steps, 4 stages each
    assert max(calls) == pytest.approx(1.0)


@pytest.mark.parametrize("bad", [0.0, -1.0, math.nan, math.inf])
def test_rejects_bad_max_step(bad: float) -> None:
    with pytest.raises(ValueError, match="max_step"):
        rk4(oscillator, 0.0, np.array([1.0, 0.0]), 1.0, max_step=bad)


@pytest.mark.parametrize("bad", [math.nan, math.inf])
def test_rejects_non_finite_times(bad: float) -> None:
    with pytest.raises(ValueError, match="t_end"):
        rk4(oscillator, 0.0, np.array([1.0, 0.0]), bad, max_step=0.1)
    with pytest.raises(ValueError, match="t0"):
        rk4(oscillator, bad, np.array([1.0, 0.0]), 1.0, max_step=0.1)
