"""Shared pytest configuration: Hypothesis profiles."""

import os
from typing import TYPE_CHECKING

import pytest
from hypothesis import HealthCheck, settings

from oneill.core.cylinder import Cylinder

if TYPE_CHECKING:
    from collections.abc import Callable

settings.register_profile(
    "dev",
    max_examples=25,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
settings.register_profile(
    "ci",
    max_examples=300,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow],
)
settings.load_profile(os.environ.get("HYPOTHESIS_PROFILE", "dev"))


CYLINDER_FACTORIES: dict[str, Callable[[], Cylinder]] = {
    "oneill-3.2km-1g": Cylinder,
    "oneill-3.2km-0.5rpm": lambda: Cylinder.from_rpm(3200.0, 20000.0, 0.5),
    "torus-900m-1rpm": lambda: Cylinder.from_rpm(900.0, 2000.0, 1.0),
    "station-250m-2rpm": lambda: Cylinder.from_rpm(250.0, 500.0, 2.0),
    "mega-10km-lunar-g": lambda: Cylinder.from_gravity(10000.0, 50000.0, 1.62),
}


@pytest.fixture(params=list(CYLINDER_FACTORIES))
def cyl(request: pytest.FixtureRequest) -> Cylinder:
    """A constant-spin habitat spanning very different sizes, spins and gravities.

    Any test requesting this fixture runs once per habitat, guarding against
    formulas that only work for the default 3.2 km / 1 g design.
    """
    return CYLINDER_FACTORIES[request.param]()


@pytest.fixture
def default_cyl() -> Cylinder:
    """The design-brief habitat: R = 3.2 km, L = 20 km, exactly 1 g at the floor."""
    return Cylinder()


@pytest.fixture
def spinup_cyl() -> Cylinder:
    """The default habitat while spinning up (alpha = +1e-4 rad/s^2)."""
    return Cylinder.from_gravity(spin_accel=1.0e-4)


@pytest.fixture(params=["constant-spin", "spinning-up"])
def flight_cyl(request: pytest.FixtureRequest) -> Cylinder:
    """Default-sized habitat, once with constant spin and once spinning up."""
    alpha = 0.0 if request.param == "constant-spin" else 1.0e-4
    return Cylinder.from_gravity(spin_accel=alpha)
