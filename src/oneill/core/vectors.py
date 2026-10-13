"""Vector types and conventions.

Habitat-frame vectors are plain ``(3,)`` or ``(n, 3)`` float arrays ``(x, y, z)``
with ``z`` along the cylinder axis. Person-centred vectors are written
``(spinward, axial, up)``: spinward is the way the floor moves, axial is along
the spin vector, and up points at the axis. That triad is left-handed.
"""

import math
from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray

type Vec3 = NDArray[np.float64]

Z_HAT: Vec3 = np.array([0.0, 0.0, 1.0])
Z_HAT.flags.writeable = False


def as_vec3(value: ArrayLike) -> Vec3:
    """Convert ``value`` to a finite float array whose last axis has length 3."""
    arr = np.asarray(value, dtype=np.float64)
    if arr.ndim == 0 or arr.shape[-1] != 3:  # noqa: PLR2004
        msg = f"expected a vector (or stack of vectors) of length 3, got {arr.shape}"
        raise ValueError(msg)
    if not np.all(np.isfinite(arr)):
        msg = "vector components must be finite"
        raise ValueError(msg)
    return arr


@dataclass(frozen=True, slots=True)
class LocalVector:
    """A vector in a person's own (spinward, axial, up) components."""

    spinward: float
    axial: float
    up: float

    def __post_init__(self) -> None:
        """Reject non-finite components."""
        for name in ("spinward", "axial", "up"):
            if not math.isfinite(getattr(self, name)):
                msg = f"{name} must be finite, got {getattr(self, name)!r}"
                raise ValueError(msg)

    @classmethod
    def from_array(cls, value: ArrayLike) -> LocalVector:
        """Build from any length-3 sequence ordered (spinward, axial, up)."""
        arr = as_vec3(value)
        if arr.ndim != 1:
            msg = f"expected a single vector, got shape {arr.shape}"
            raise ValueError(msg)
        return cls(float(arr[0]), float(arr[1]), float(arr[2]))

    def to_array(self) -> Vec3:
        """Return the components as a ``(3,)`` float array."""
        return np.array([self.spinward, self.axial, self.up])

    @property
    def speed(self) -> float:
        """Euclidean length of the vector."""
        return math.sqrt(self.spinward**2 + self.axial**2 + self.up**2)
