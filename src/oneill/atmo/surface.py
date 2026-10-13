"""L2b: surface-layer wind over buildings and trees, and rotational stability.

Frame (``docs/research/atmosphere-literature.md``, "Sign convention"): ``x``
spinward, ``y`` axial, ``z`` up (toward the axis), right-handed, so the spin
vector is ``Omega = -Omega y_hat``. Note that :class:`oneill.core.vectors.LocalVector`
stores ``(spinward, axial, up)`` with axial along ``+Z_HAT``, a left-handed triple;
nothing here depends on that sign because only the spinward component of the wind
shear matters to the Bradshaw criterion.

The wind profile is a log law above the canopy and an exponential (Cionco-type)
profile inside it, joined with continuous speed *and* shear. The exponential
attenuation coefficient is therefore not a free parameter: it is fixed by shear
continuity at canopy top (:attr:`Terrain.attenuation`).

Caveat: the canopy attenuation implied by shear continuity is small for the
suburb (0.87, against roughly 1-3 in Macdonald-type urban fits), and the whole
profile is an area-averaged canopy result. The unstable band it predicts is a
property of that horizontally averaged profile; near-wall shear around individual
buildings is not represented.

Rotational stability follows Brethouwer (2017, JFM 828:424): with
``S = -2 Omega / U'`` and ``U' = dU/dz`` of the spinward component, the Bradshaw
number ``B = S (S + 1)`` is negative (shear turbulence *augmented*) exactly for
``-1 < S < 0``, i.e. ``U' > 2 Omega``.
"""

import math
from dataclasses import dataclass
from enum import IntEnum
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import finite, positive

if TYPE_CHECKING:
    from numpy.typing import ArrayLike, NDArray

KAPPA = 0.40
"""von Karman constant (model-parameters.md, L2b)."""

SIGMA_U_OVER_USTAR = 2.35
"""Streamwise turbulence intensity ``sigma_u / u*`` (2.2-2.5; midpoint)."""

PEAK_FACTOR = 3.0
"""Gaussian peak factor for a ~3 s gust in a 10 min record (our estimate, MED)."""


CANOPY_NOTE = (
    "area-averaged canopy profile; near-wall shear around individual buildings "
    "is not represented"
)
"""Caveat attached to stability results."""


@dataclass(frozen=True, slots=True)
class Terrain:
    """Surface roughness description: ``(z0, d, h_c)`` plus a name.

    Attributes:
        name: Label for reports.
        z0: Roughness length in m.
        d: Zero-plane displacement height in m.
        h_c: Canopy (obstacle) height in m.
    """

    name: str
    z0: float
    d: float
    h_c: float

    def __post_init__(self) -> None:
        """Require positive lengths with the canopy top more than ``z0`` above ``d``."""
        positive("z0", self.z0)
        positive("d", self.d)
        positive("h_c", self.h_c)
        if self.h_c - self.d <= self.z0:
            msg = "canopy top must sit more than z0 above the displacement height"
            raise ValueError(msg)

    @property
    def attenuation(self) -> float:
        """Exponential canopy attenuation ``a = h_c / ((h_c - d) ln((h_c - d)/z0))``.

        Chosen so that ``dU/dz`` is continuous across the canopy top. For the
        suburb this is 0.87, low next to Macdonald-type urban fits of about 1-3,
        so in-canopy winds there are probably overestimated.
        """
        gap = self.h_c - self.d
        return self.h_c / (gap * math.log(gap / self.z0))


# Presets from model-parameters.md (L2b surface layer): suburb z0 0.3-0.8 m,
# d 3-7 m (h ~ 10 m); forest/orchard z0 ~ 0.1 h, d ~ 2/3 h (Raupach 1994).
SUBURBAN = Terrain("suburban 2-3 storey", z0=0.5, d=5.0, h_c=10.0)
"""Low-rise suburb: h = 10 m, d = 0.5 h, z0 = 0.5 m (Grimmond & Oke 1999; Q9)."""
ORCHARD = Terrain("orchard", z0=0.6, d=4.0, h_c=6.0)
"""Orchard: h = 6 m, d = 2/3 h, z0 = 0.1 h (Raupach 1994; Q9)."""
FOREST = Terrain("forest", z0=1.5, d=10.0, h_c=15.0)
"""Forest: h = 15 m, d = 2/3 h, z0 = 0.1 h (Raupach 1994; Q9)."""


class Stability(IntEnum):
    """Rotational (Bradshaw) stability class of shear turbulence."""

    UNSTABLE = -1
    NEUTRAL = 0
    STABLE = 1


@dataclass(frozen=True, slots=True)
class StabilityProfile:
    """Rotational stability along a height profile.

    Attributes:
        z: Heights above the floor in m.
        shear: Spinward-component shear ``U'`` in 1/s.
        bradshaw: ``B = S (S + 1)``; NaN where there is no spinward shear.
        classification: :class:`Stability` values as an integer array.
        note: Reminder that the profile is an area-averaged canopy result.
    """

    z: NDArray[np.float64]
    shear: NDArray[np.float64]
    bradshaw: NDArray[np.float64]
    classification: NDArray[np.int8]
    note: str = CANOPY_NOTE


def _heights(z: ArrayLike) -> NDArray[np.float64]:
    arr = np.asarray(z, dtype=np.float64)
    if not np.all(arr > 0.0):
        msg = "heights must be positive"
        raise ValueError(msg)
    return arr


def log_law_wind(
    z: ArrayLike, ustar: float, z0: float, d: float
) -> NDArray[np.float64]:
    """Neutral log-law speed ``(u*/kappa) ln((z - d)/z0)`` in m/s.

    It is zero at ``z = d + z0`` (and negative below, where the law does not apply).

    Raises:
        ValueError: if any ``z <= d``.
    """
    arr = np.asarray(z, dtype=np.float64)
    positive("z0", z0)
    if not np.all(arr > d):
        msg = "log law needs heights above the displacement height d"
        raise ValueError(msg)
    return np.asarray(finite("ustar", ustar) / KAPPA * np.log((arr - d) / z0))


def wind_speed(z: ArrayLike, ustar: float, terrain: Terrain) -> NDArray[np.float64]:
    """Mean wind speed at heights ``z``: log law above ``h_c``, exponential below.

    Inside the canopy ``U(z) = U(h_c) exp(a (z/h_c - 1))`` with
    :attr:`Terrain.attenuation` ``a``.
    """
    arr = _heights(z)
    top = np.maximum(arr, terrain.h_c)
    u_top = ustar / KAPPA * math.log((terrain.h_c - terrain.d) / terrain.z0)
    above = ustar / KAPPA * np.log((top - terrain.d) / terrain.z0)
    below = u_top * np.exp(
        terrain.attenuation * (np.minimum(arr, terrain.h_c) / terrain.h_c - 1.0)
    )
    return np.asarray(np.where(arr >= terrain.h_c, above, below))


def wind_shear(z: ArrayLike, ustar: float, terrain: Terrain) -> NDArray[np.float64]:
    """Analytic ``dU/dz`` (1/s) of :func:`wind_speed`: ``u*/(kappa (z-d))`` above."""
    arr = _heights(z)
    above = ustar / (KAPPA * (np.maximum(arr, terrain.h_c) - terrain.d))
    below = terrain.attenuation / terrain.h_c * wind_speed(arr, ustar, terrain)
    return np.asarray(np.where(arr >= terrain.h_c, above, below))


def friction_velocity(speed: float, height: float, terrain: Terrain) -> float:
    """Friction velocity ``u*`` giving ``speed`` m/s at ``height`` (inverse profile).

    Over the suburb ``u*/U`` is about 0.12 at 10-15 m above ``d`` (Q9).
    """
    positive("speed", speed)
    return speed / float(wind_speed(positive("height", height), 1.0, terrain))


def gust_speed(
    mean_speed: float,
    ustar: float,
    peak_factor: float = PEAK_FACTOR,
    sigma_ratio: float = SIGMA_U_OVER_USTAR,
) -> float:
    """Peak gust ``U + k sigma_u`` with ``sigma_u = 2.35 u*`` (Q9, MED).

    Over a 10 m suburb this gives a gust factor of about 1.9, inside the research
    range of 1.7-2.0.
    """
    return mean_speed + peak_factor * sigma_ratio * positive("ustar", ustar)


def dynamic_pressure(density: float, speed: float) -> float:
    """Dynamic pressure ``rho U^2 / 2`` in Pa, the basis of wind loads."""
    return 0.5 * positive("density", density) * finite("speed", speed) ** 2


def spinward_fraction(heading_deg: float) -> float:
    """Spinward component of a unit horizontal wind blowing toward ``heading_deg``.

    Headings follow the package convention (0 = +axial, 90 = spinward,
    180 = -axial, 270 = antispinward); the result is ``sin(heading)`` snapped
    so exactly axial headings give an exact zero.
    """
    return round(math.sin(math.radians(finite("heading_deg", heading_deg))), 12)


def bradshaw_number(shear: float, omega: float) -> float:
    """Bradshaw number ``B = S (S + 1)`` with ``S = -2 omega / U'``.

    Raises:
        ValueError: if ``shear`` is zero (``S`` undefined).
    """
    if shear == 0.0:
        msg = "Bradshaw number is undefined without shear"
        raise ValueError(msg)
    s = -2.0 * positive("omega", omega) / shear
    return s * (s + 1.0)


def _bradshaw(shear: NDArray[np.float64], omega: float) -> NDArray[np.float64]:
    """Vectorised ``B = S (S + 1)``, NaN where there is no shear."""
    s = np.full_like(shear, np.nan)
    np.divide(-2.0 * omega, shear, out=s, where=shear != 0.0)
    with np.errstate(over="ignore"):  # denormal shear: B overflows to +inf (stable)
        return np.asarray(s * (s + 1.0))


def _classify(bradshaw: NDArray[np.float64]) -> NDArray[np.int8]:
    out = np.zeros(bradshaw.shape, dtype=np.int8)
    out[bradshaw > 0.0] = Stability.STABLE
    out[bradshaw < 0.0] = Stability.UNSTABLE
    return out


def classify_rotational_stability(shear: float, omega: float) -> Stability:
    """Classify spinward-component shear ``U'`` against the spin rate.

    ``U' > 2 omega`` is UNSTABLE (anticyclonic, ``-1 < S < 0``); ``0 < U' < 2 omega``
    and every ``U' < 0`` (antispinward, cyclonic) are STABLE; zero shear (axial
    wind) and the marginal ``U' = 2 omega`` are NEUTRAL.
    """
    positive("omega", omega)
    shear_arr = np.asarray(finite("shear", shear))
    return Stability(int(_classify(_bradshaw(shear_arr, omega))))


def rotational_stability_profile(
    z: ArrayLike,
    ustar: float,
    terrain: Terrain,
    omega: float,
    spinward: float,
) -> StabilityProfile:
    """Stability class at heights ``z`` for a wind with a given spinward fraction.

    Args:
        z: Heights above the floor in m.
        ustar: Friction velocity in m/s.
        terrain: Roughness description.
        omega: Spin rate in rad/s.
        spinward: Spinward component of the wind direction in ``[-1, 1]``:
            +1 spinward, -1 antispinward, 0 axial (see :func:`spinward_fraction`).
    """
    positive("omega", omega)
    if not -1.0 <= spinward <= 1.0:
        msg = f"spinward must lie in [-1, 1], got {spinward!r}"
        raise ValueError(msg)
    heights = _heights(z)
    shear = spinward * wind_shear(heights, ustar, terrain)
    bradshaw = _bradshaw(shear, omega)
    return StabilityProfile(
        z=heights,
        shear=shear,
        bradshaw=bradshaw,
        classification=_classify(bradshaw),
    )


def unstable_layer(
    ustar: float, terrain: Terrain, omega: float
) -> tuple[float, float] | None:
    """Heights ``(z_low, z_high)`` between which a spinward wind is unstable.

    The layer is where ``U' > 2 omega``. Above the canopy the shear is
    ``u*/(kappa (z - d))``, so the top is ``z_high = d + u*/(2 omega kappa)``
    (2-7 m above ``d`` for ``u*`` of 0.1-0.3 m/s, Q5). Inside the canopy the shear
    falls off exponentially, giving ``z_low = h_c (1 + ln(2 omega / U'_top)/a)``,
    floored at 0.

    This is an area-averaged canopy result (:data:`CANOPY_NOTE`).

    Returns:
        ``None`` if the shear at canopy top is already at most ``2 omega``.
    """
    positive("omega", omega)
    positive("ustar", ustar)
    top_shear = ustar / (KAPPA * (terrain.h_c - terrain.d))
    if top_shear <= 2.0 * omega:
        return None
    z_low = terrain.h_c * (
        1.0 + math.log(2.0 * omega / top_shear) / terrain.attenuation
    )
    z_high = terrain.d + ustar / (2.0 * omega * KAPPA)
    return max(z_low, 0.0), z_high
