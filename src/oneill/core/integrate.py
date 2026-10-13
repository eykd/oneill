"""A small generic fixed-step RK4 integrator.

Shared by the free-body ballistics cross-check and (later) the air-parcel model.
"""

import math
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from collections.abc import Callable

    from numpy.typing import NDArray

    type State = NDArray[np.float64]
    type Derivative = Callable[[float, State], State]


def rk4(
    deriv: Derivative, t0: float, y0: State, t_end: float, max_step: float
) -> State:
    """Integrate ``dy/dt = deriv(t, y)`` from ``t0`` to ``t_end``; return ``y(t_end)``.

    The interval is split into equal steps no longer than ``max_step``. Integrating
    backwards (``t_end < t0``) is allowed. ``y0`` is never modified.
    """
    finite("t0", t0)
    finite("t_end", t_end)
    positive("max_step", max_step)
    y = np.array(y0, dtype=np.float64)
    n = max(1, math.ceil(abs(t_end - t0) / max_step))
    if t_end == t0:
        return y
    h = (t_end - t0) / n
    t = t0
    for _ in range(n):
        k1 = deriv(t, y)
        k2 = deriv(t + 0.5 * h, y + 0.5 * h * k1)
        k3 = deriv(t + 0.5 * h, y + 0.5 * h * k2)
        k4 = deriv(t + h, y + h * k3)
        y = y + h / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t += h
    return y
