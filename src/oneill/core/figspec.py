"""Plain-data figure specifications shared by the motion and atmosphere layers.

Nothing here draws. ``adapters.plots`` turns a ``FigureSpec`` into pixels.
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from oneill.core.vectors import Vec3


@dataclass(frozen=True, eq=False)
class Series:
    """One curve; ``fmt`` is a matplotlib format string ("" for the default)."""

    label: str
    x: Vec3
    y: Vec3
    fmt: str = ""


@dataclass(frozen=True, eq=False)
class Band:
    """A shaded horizontal band between two y values, with a legend label."""

    lo: float
    hi: float
    label: str


@dataclass(frozen=True, eq=False)
class Panel:
    """One set of axes."""

    xlabel: str
    ylabel: str
    series: tuple[Series, ...]
    title: str = ""
    hlines: tuple[float, ...] = ()
    xticks: tuple[float, ...] = ()
    xticklabels: tuple[str, ...] = ()
    legend_size: float | None = None
    bands: tuple[Band, ...] = ()


@dataclass(frozen=True, eq=False)
class FigureSpec:
    """A figure: a grid of panels, a file name and a size in inches."""

    filename: str
    panels: tuple[Panel, ...]
    layout: tuple[int, int]
    figsize: tuple[float, float]
    suptitle: str = ""
    sharex: bool = False
