"""Small argument validators shared by the value objects."""

import math


def finite(name: str, value: float) -> float:
    """Return ``value`` as a float, raising ``ValueError`` unless it is finite."""
    if not math.isfinite(value):
        msg = f"{name} must be finite, got {value!r}"
        raise ValueError(msg)
    return float(value)


def positive(name: str, value: float) -> float:
    """Return ``value`` as a float, raising ``ValueError`` unless finite and > 0."""
    if not (math.isfinite(value) and value > 0.0):
        msg = f"{name} must be positive and finite, got {value!r}"
        raise ValueError(msg)
    return float(value)
