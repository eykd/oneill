"""Compass headings for lateral motion.

0 is +axial, 90 is spinward, 180 is -axial and 270 is antispinward.
"""

import math

from oneill.core.vectors import LocalVector

HEADING_NAMES = {
    0: "+axial",
    45: "+axial&spinward",
    90: "spinward",
    135: "-axial&spinward",
    180: "-axial",
    225: "-axial&antispin",
    270: "antispinward",
    315: "+axial&antispin",
}

_NAME_TOLERANCE = 1e-9


def heading_name(heading_deg: float) -> str:
    """Name a heading, wrapping modulo 360 degrees; unnamed ones print in degrees."""
    wrapped = heading_deg % 360.0
    nearest = round(wrapped)
    if abs(wrapped - nearest) < _NAME_TOLERANCE and nearest % 360 in HEADING_NAMES:
        return HEADING_NAMES[nearest % 360]
    return f"{wrapped:g} deg"


def heading_components(heading_deg: float, elevation_deg: float = 0.0) -> LocalVector:
    """Unit direction as (spinward, axial, up).

    Heading 0 is +axial and 90 is spinward; positive elevation points upward
    (toward the axis).
    """
    h, e = math.radians(heading_deg), math.radians(elevation_deg)
    return LocalVector(
        math.cos(e) * math.sin(h), math.cos(e) * math.cos(h), math.sin(e)
    )
