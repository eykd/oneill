"""Fictitious accelerations in the co-rotating habitat frame.

For a frame spinning at ``omega`` about z (angular acceleration ``alpha``), a
free body of unit mass obeys ``a = centrifugal + coriolis + euler`` with::

    centrifugal = -omega x (omega x p)  =  omega^2 (x, y, 0)
    coriolis    = -2 omega x v          =  2 omega (v_y, -v_x, 0)
    euler       = -alpha x p            =  alpha (y, -x, 0)

All functions are vectorised: ``omega`` and ``alpha`` may be scalars or arrays
that broadcast against the leading axes of ``p`` / ``v`` (shape ``(..., 3)``).
"""

from typing import TYPE_CHECKING

import numpy as np

from oneill.core.vectors import Vec3, as_vec3

if TYPE_CHECKING:
    from numpy.typing import ArrayLike

    from oneill.core.cylinder import Cylinder


def _planar(x: ArrayLike, y: ArrayLike) -> Vec3:
    """Stack in-plane components into vectors with a zero axial part."""
    xx, yy = np.broadcast_arrays(np.asarray(x, np.float64), np.asarray(y, np.float64))
    return np.stack([xx, yy, np.zeros_like(xx)], axis=-1)


def spin_rate(cyl: Cylinder, t: ArrayLike) -> Vec3:
    """Spin rate ``omega(t)`` at scalar or array time ``t`` (rad/s)."""
    return cyl.omega + cyl.spin_accel * np.asarray(t, dtype=np.float64)


def centrifugal_acc(omega: ArrayLike, p: ArrayLike) -> Vec3:
    """Centrifugal acceleration at habitat position(s) ``p``."""
    pts = as_vec3(p)
    w2 = np.asarray(omega, dtype=np.float64) ** 2
    return _planar(w2 * pts[..., 0], w2 * pts[..., 1])


def coriolis_acc(omega: ArrayLike, v: ArrayLike) -> Vec3:
    """Coriolis acceleration for habitat-frame velocity(ies) ``v``."""
    vel = as_vec3(v)
    w = np.asarray(omega, dtype=np.float64)
    return _planar(2.0 * w * vel[..., 1], -2.0 * w * vel[..., 0])


def euler_acc(alpha: ArrayLike, p: ArrayLike) -> Vec3:
    """Euler acceleration at habitat position(s) ``p`` (zero at constant spin)."""
    pts = as_vec3(p)
    a = np.asarray(alpha, dtype=np.float64)
    return _planar(a * pts[..., 1], -a * pts[..., 0])


def fictitious_acc(cyl: Cylinder, t: ArrayLike, p: ArrayLike, v: ArrayLike) -> Vec3:
    """Total fictitious acceleration of a free body at time ``t``."""
    return (
        centrifugal_acc(spin_rate(cyl, t), p)
        + coriolis_acc(spin_rate(cyl, t), v)
        + euler_acc(cyl.spin_accel, p)
    )
