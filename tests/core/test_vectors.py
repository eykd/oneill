import numpy as np
import pytest
from hypothesis import given
from hypothesis import strategies as st

from oneill.core.checks import finite, positive
from oneill.core.vectors import Z_HAT, LocalVector, as_vec3

components = st.floats(-1e6, 1e6)


def test_z_hat_is_read_only_unit_axis() -> None:
    np.testing.assert_array_equal(Z_HAT, [0.0, 0.0, 1.0])
    with pytest.raises(ValueError, match="read-only"):
        Z_HAT[0] = 1.0


def test_as_vec3_accepts_single_and_stacked_vectors() -> None:
    assert as_vec3([1, 2, 3]).dtype == np.float64
    assert as_vec3([[1, 2, 3], [4, 5, 6]]).shape == (2, 3)


@pytest.mark.parametrize("bad", [5.0, [1.0, 2.0], [[1.0, 2.0, 3.0, 4.0]]])
def test_as_vec3_rejects_wrong_shapes(bad: object) -> None:
    with pytest.raises(ValueError, match="length 3"):
        as_vec3(bad)  # type: ignore[arg-type]


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_as_vec3_rejects_non_finite(bad: float) -> None:
    with pytest.raises(ValueError, match="finite"):
        as_vec3([0.0, bad, 1.0])


@given(components, components, components)
def test_local_vector_array_round_trip(a: float, b: float, c: float) -> None:
    v = LocalVector(a, b, c)
    np.testing.assert_array_equal(v.to_array(), [a, b, c])
    assert LocalVector.from_array(v.to_array()) == v
    assert v.speed == pytest.approx(float(np.linalg.norm([a, b, c])), rel=1e-12)


def test_local_vector_rejects_non_finite_and_stacks() -> None:
    with pytest.raises(ValueError, match="axial must be finite"):
        LocalVector(0.0, float("nan"), 0.0)
    with pytest.raises(ValueError, match="single vector"):
        LocalVector.from_array([[1.0, 2.0, 3.0]])


def test_finite_and_positive() -> None:
    assert finite("x", 2) == 2.0
    assert positive("x", 2) == 2.0
    for bad in (float("nan"), float("inf")):
        with pytest.raises(ValueError, match="must be finite"):
            finite("x", bad)
        with pytest.raises(ValueError, match="positive and finite"):
            positive("x", bad)
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError, match="positive and finite"):
            positive("x", bad)
