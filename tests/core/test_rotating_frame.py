import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.core.cylinder import Cylinder
from oneill.core.rotating_frame import (
    centrifugal_acc,
    coriolis_acc,
    euler_acc,
    fictitious_acc,
    spin_rate,
)
from oneill.core.vectors import Z_HAT

vec = st.tuples(st.floats(-5e3, 5e3), st.floats(-5e3, 5e3), st.floats(-5e3, 5e3)).map(
    np.array
)
rates = st.floats(-2.0, 2.0)


@given(rates, vec)
@example(0.05, np.array([3200.0, 0.0, 0.0]))
def test_centrifugal_matches_the_double_cross_product(
    omega: float, p: np.ndarray
) -> None:
    w = omega * Z_HAT
    np.testing.assert_allclose(
        centrifugal_acc(omega, p), -np.cross(w, np.cross(w, p)), atol=1e-9, rtol=1e-12
    )


@given(rates, vec)
def test_coriolis_matches_the_cross_product(omega: float, v: np.ndarray) -> None:
    np.testing.assert_allclose(
        coriolis_acc(omega, v), -2.0 * np.cross(omega * Z_HAT, v), atol=1e-9
    )


@given(rates, vec)
def test_euler_matches_the_cross_product(alpha: float, p: np.ndarray) -> None:
    np.testing.assert_allclose(
        euler_acc(alpha, p), -np.cross(alpha * Z_HAT, p), atol=1e-9
    )


@given(rates, vec)
def test_coriolis_does_no_work_and_has_no_axial_part(
    omega: float, v: np.ndarray
) -> None:
    cor = coriolis_acc(omega, v)
    assert float(np.dot(cor, v)) == pytest.approx(0.0, abs=1e-6 * (1 + float(v @ v)))
    assert cor[2] == 0.0


def test_stacked_inputs_broadcast_against_per_sample_rates() -> None:
    p = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    omega = np.array([0.5, 1.0])
    out = centrifugal_acc(omega, p)
    assert out.shape == (2, 3)
    np.testing.assert_allclose(out[1], centrifugal_acc(1.0, p[1]))


def test_spin_rate_is_linear_in_time() -> None:
    cyl = Cylinder.from_gravity(spin_accel=1e-4)
    ts = np.array([0.0, 10.0, 100.0])
    np.testing.assert_allclose(spin_rate(cyl, ts), cyl.omega + 1e-4 * ts, rtol=1e-14)
    assert spin_rate(cyl, 5.0) == pytest.approx(cyl.omega_at(5.0))


def test_fictitious_acc_is_the_sum_of_the_three_terms() -> None:
    cyl = Cylinder.from_gravity(spin_accel=1e-4)
    p = cyl.floor_point(0.3, 5.0, 200.0)
    v = np.array([1.0, -2.0, 3.0])
    w = cyl.omega_at(7.0)
    np.testing.assert_allclose(
        fictitious_acc(cyl, 7.0, p, v),
        centrifugal_acc(w, p) + coriolis_acc(w, v) + euler_acc(1e-4, p),
    )
