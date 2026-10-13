"""Fictitious forces and free flight inside a spinning O'Neill cylinder.

The public API re-exported here is everything most callers need; see the
``core``, ``motion`` and ``adapters`` sub-packages for the full layers.
"""

from importlib.metadata import version

from oneill.core.cylinder import (
    G_STANDARD,
    AxisSingularityError,
    Cylinder,
    HabitatSummary,
)
from oneill.core.integrate import rk4
from oneill.core.vectors import LocalVector
from oneill.motion.ballistics import (
    BallisticResult,
    ExtraAcceleration,
    ballistic,
    integrate_rotating,
)
from oneill.motion.forces import ForceHistory, analyze_path
from oneill.motion.headings import heading_components, heading_name
from oneill.motion.paths import (
    FloorWalk,
    NumericPath,
    Path,
    StraightPath,
    velocity_path,
)
from oneill.motion.profiles import ConstantSpeed, Profile, Trapezoid
from oneill.motion.scenarios import (
    ballistics_report,
    custom_report,
    elevator_report,
    lateral_report,
    stairs_report,
    vertical_report,
)
from oneill.motion.verification import CheckResult, all_passed, run_selftest

__version__ = version("oneill")

__all__ = [
    "G_STANDARD",
    "AxisSingularityError",
    "BallisticResult",
    "CheckResult",
    "ConstantSpeed",
    "Cylinder",
    "ExtraAcceleration",
    "FloorWalk",
    "ForceHistory",
    "HabitatSummary",
    "LocalVector",
    "NumericPath",
    "Path",
    "Profile",
    "StraightPath",
    "Trapezoid",
    "__version__",
    "all_passed",
    "analyze_path",
    "ballistic",
    "ballistics_report",
    "custom_report",
    "elevator_report",
    "heading_components",
    "heading_name",
    "integrate_rotating",
    "lateral_report",
    "rk4",
    "run_selftest",
    "stairs_report",
    "velocity_path",
    "vertical_report",
]
