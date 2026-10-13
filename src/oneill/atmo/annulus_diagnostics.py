"""Diagnostics of the 2-D annulus convection runs (pure functions of model state).

These answer hypothesis H1 (``docs/design/annulus2d.md`` section 6): the dominant
eddy scale, the wind at 10 m, gusts, heat exchange, lock-in to the three strips, the
mean (zonal) wind, and whether a run has reached a statistically stationary state.
"""

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from oneill.core.checks import positive

if TYPE_CHECKING:
    from numpy.typing import NDArray

    from oneill.atmo.annulus2d import Annulus, AnnulusState

    type Array = NDArray[np.float64]

WIND_HEIGHT = 10.0
"""Height in m of the wind sample (the canopy-drag reference height)."""

PEAK_FACTOR = 3.0
"""Gaussian peak factor of the gust estimate (as in :mod:`oneill.atmo.surface`)."""

_TINY = 1.0e-30


@dataclass(frozen=True, eq=False)
class Sample:
    """One snapshot of the diagnostics.

    Attributes:
        t: Time in s.
        energy: Total kinetic energy ``int (u^2)/2 dA`` per unit length.
        mean_energy: Energy of the theta-mean wind alone.
        enstrophy: ``int zeta^2/2 dA``.
        spectrum: Kinetic energy by azimuthal wavenumber over the outer half of the
            annulus (entry 0 is the mean wind).
        wind: Azimuthal wind ``u_theta(theta)`` at :data:`WIND_HEIGHT`, m/s.
        ubar: Theta-mean azimuthal wind profile, m/s.
        nusselt: Axis-rim heat-exchange number (see :func:`take_sample`).
        dissipation: Viscous dissipation ``nu <zeta^2>`` per unit mass, W/kg.
    """

    t: float
    energy: float
    mean_energy: float
    enstrophy: float
    spectrum: Array
    wind: Array
    ubar: Array
    nusselt: float
    dissipation: float


@dataclass(frozen=True, slots=True)
class Scale:
    """Dominant azimuthal scale of a spectrum.

    Attributes:
        m_peak: Wavenumber of the spectral maximum (``m >= 1``).
        m_centroid: Energy-weighted mean wavenumber.
        length_peak: ``2 pi r / m_peak`` in m.
        length_centroid: ``2 pi r / m_centroid`` in m.
    """

    m_peak: int
    m_centroid: float
    length_peak: float
    length_centroid: float


@dataclass(frozen=True, slots=True)
class Stationarity:
    """Result of the stationarity test on a time series.

    Attributes:
        stationary: ``True`` if the mean of the last quarter differs from the mean
            of the preceding quarter by less than ``tol``.
        drift: That relative difference (positive if the series is growing).
    """

    stationary: bool
    drift: float


def take_sample(model: Annulus, state: AnnulusState) -> Sample:
    """Diagnostics of ``state``.

    The heat-exchange number is ``Nu = q / (kappa dT / (R - r_in))`` where ``q`` is
    the mean over theta of the positive part of ``-u_r T'`` at mid radius (``T'`` is
    the departure from the theta-mean; ``-u_r`` is inward, the "up" direction) and
    ``dT`` is the peak-to-peak temperature contrast around the outer wall cell. It is
    zero without convection and large when plumes carry heat across the annulus.
    """
    c = model.config
    u_r, _, _ = model.velocities(state)
    temp = model.to_phys(state.temp)
    j = int(np.argmin(np.abs(model.r - 0.5 * (model.r_in + c.radius))))
    t_prime = temp[j] - temp[j].mean()
    q = float(np.mean(np.maximum(-u_r[j] * t_prime, 0.0)))
    contrast = float(np.ptp(temp[-1]))
    nusselt = q * (c.radius - model.r_in) / (c.kappa * max(contrast, _TINY))
    zeta_rms2 = model.enstrophy(state) / (math.pi * np.sum(model.area))  # <zeta^2>
    return Sample(
        t=state.t,
        energy=model.energy(state),
        mean_energy=model.mean_wind_energy(state),
        enstrophy=model.enstrophy(state),
        spectrum=model.spectrum(state, 0.5 * c.radius, c.radius),
        wind=model.wind_at_height(state, WIND_HEIGHT),
        ubar=np.array(state.ubar),
        nusselt=nusselt,
        dissipation=c.nu * zeta_rms2,
    )


def spectrum_scale(spectrum: Array, r_ref: float) -> Scale:
    """Peak and centroid scale of the eddy part (``m >= 1``) of ``spectrum``."""
    eddy = np.asarray(spectrum[1:], dtype=np.float64)
    total = float(eddy.sum())
    m = np.arange(1, len(eddy) + 1, dtype=np.float64)
    peak = int(np.argmax(eddy)) + 1
    centroid = float(np.sum(m * eddy) / total) if total > 0.0 else float(peak)
    two_pi_r = 2.0 * math.pi * positive("r_ref", r_ref)
    return Scale(peak, centroid, two_pi_r / peak, two_pi_r / centroid)


def lock_fraction(spectrum: Array, strips: int = 3) -> float:
    """Fraction of the eddy energy in the strip wavenumber (``m = strips``).

    Near 1 means the circulation is locked to the heating pattern; a small value with
    energy at ``m = 1, 2`` is the signature of an inverse cascade to a large vortex.
    """
    eddy = float(np.sum(spectrum[1:]))
    return float(spectrum[strips] / eddy) if eddy > 0.0 else 0.0


def large_scale_fraction(spectrum: Array, m_max: int = 2) -> float:
    """Fraction of the eddy energy at ``1 <= m <= m_max`` (condensate indicator)."""
    eddy = float(np.sum(spectrum[1:]))
    return float(np.sum(spectrum[1 : m_max + 1]) / eddy) if eddy > 0.0 else 0.0


def wind_statistics(wind: Array) -> tuple[float, float]:
    """RMS and gust (mean ``|u|`` plus 3 sigma) of the azimuthal wind samples."""
    rms = float(np.sqrt(np.mean(wind**2)))
    gust = float(np.mean(np.abs(wind)) + PEAK_FACTOR * np.std(wind))
    return rms, gust


def stationarity(t: Array, series: Array, tol: float = 0.1) -> Stationarity:
    """Is ``series`` statistically stationary at its end?

    Compares the mean over the last quarter of the record with the mean over the
    quarter before it.
    """
    n = len(series)
    quarter = max(1, n // 4)
    last = float(np.mean(series[n - quarter :]))
    before = (
        float(np.mean(series[max(0, n - 2 * quarter) : n - quarter]))
        if n >= 2 * quarter
        else last
    )
    scale = max(abs(last), abs(before), _TINY)
    drift = (last - before) / scale
    return Stationarity(abs(drift) < tol and len(t) >= 4, drift)  # noqa: PLR2004
