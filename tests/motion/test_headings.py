import math

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.motion.headings import heading_components, heading_name


@pytest.mark.parametrize(
    ("heading", "elevation", "expected"),
    [
        (0.0, 0.0, (0.0, 1.0, 0.0)),
        (90.0, 0.0, (1.0, 0.0, 0.0)),
        (180.0, 0.0, (0.0, -1.0, 0.0)),
        (270.0, 0.0, (-1.0, 0.0, 0.0)),
        (360.0, 0.0, (0.0, 1.0, 0.0)),
        (-90.0, 0.0, (-1.0, 0.0, 0.0)),
        (45.0, 0.0, (math.sqrt(0.5), math.sqrt(0.5), 0.0)),
        (0.0, 90.0, (0.0, 0.0, 1.0)),
        (123.0, -90.0, (0.0, 0.0, -1.0)),
        (
            90.0,
            35.0,
            (math.cos(math.radians(35.0)), 0.0, math.sin(math.radians(35.0))),
        ),
    ],
)
def test_heading_components(
    heading: float, elevation: float, expected: tuple[float, float, float]
) -> None:
    vec = heading_components(heading, elevation)
    np.testing.assert_allclose(vec.to_array(), expected, atol=1e-12)


@given(st.floats(-720.0, 720.0), st.floats(-90.0, 90.0))
@example(0.0, 90.0)
@example(359.0, -90.0)
def test_heading_components_are_unit_vectors(heading: float, elevation: float) -> None:
    assert heading_components(heading, elevation).speed == pytest.approx(1.0, abs=1e-12)


@given(st.floats(-720.0, 720.0), st.floats(-90.0, 90.0))
def test_spinward_part_of_a_level_heading_is_sin_heading(
    heading: float, elevation: float
) -> None:
    vec = heading_components(heading, elevation)
    assert vec.spinward == pytest.approx(
        math.cos(math.radians(elevation)) * math.sin(math.radians(heading)),
        abs=1e-12,
    )
    assert vec.up == pytest.approx(math.sin(math.radians(elevation)), abs=1e-12)


@pytest.mark.parametrize(
    ("heading", "expected"),
    [
        (0, "+axial"),
        (45.0, "+axial&spinward"),
        (90, "spinward"),
        (135, "-axial&spinward"),
        (180, "-axial"),
        (225, "-axial&antispin"),
        (270, "antispinward"),
        (315, "+axial&antispin"),
        (360, "+axial"),
        (450, "spinward"),
        (-90, "antispinward"),
        (-1e-12, "+axial"),
        (30, "30 deg"),
        (90.5, "90.5 deg"),
        (390, "30 deg"),
        (-30, "330 deg"),
    ],
)
def test_heading_name(heading: float, expected: str) -> None:
    assert heading_name(heading) == expected


@given(st.integers(-1000, 1000), st.sampled_from([0.0, 0.5]), st.integers(-5, 5))
def test_heading_name_is_periodic(whole: int, frac: float, turns: int) -> None:
    heading = whole + frac
    assert heading_name(heading) == heading_name(heading + 360.0 * turns)
