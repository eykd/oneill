"""Plot adapter: draws ``FigureSpec``s with matplotlib.

matplotlib is an optional dependency (``pip install 'oneill[plot]'``) and is
imported lazily, only here, so importing this module is always cheap and safe.
"""

from typing import TYPE_CHECKING, Any

import numpy as np

if TYPE_CHECKING:
    from collections.abc import Callable, Iterable
    from pathlib import Path
    from types import ModuleType

    from oneill.motion.figures import FigureSpec

_DPI = 120
_INSTALL_HINT = "plotting needs matplotlib: pip install 'oneill[plot]'"


class PlottingUnavailableError(Exception):
    """Raised when matplotlib is not installed."""


def _pyplot(*, interactive: bool) -> ModuleType:
    try:
        import matplotlib as mpl  # noqa: PLC0415

        if not interactive:
            mpl.use("Agg")
        import matplotlib.pyplot as plt  # noqa: PLC0415
    except ImportError as exc:
        raise PlottingUnavailableError(_INSTALL_HINT) from exc
    return plt


def check_available() -> None:
    """Raise ``PlottingUnavailableError`` unless matplotlib can be imported."""
    try:
        import matplotlib as mpl  # noqa: F401, PLC0415
    except ImportError as exc:
        raise PlottingUnavailableError(_INSTALL_HINT) from exc


def show_all() -> None:
    """Show every open figure in an interactive window and block until closed."""
    _pyplot(interactive=True).show()


def _draw(plt: ModuleType, spec: FigureSpec) -> Any:  # noqa: ANN401
    rows, cols = spec.layout
    fig, axes = plt.subplots(rows, cols, figsize=spec.figsize, sharex=spec.sharex)
    for ax, panel in zip(np.atleast_1d(axes).ravel(), spec.panels, strict=True):
        for s in panel.series:
            ax.plot(s.x, s.y, s.fmt, label=s.label)
        for level in panel.hlines:
            ax.axhline(level, color="k", lw=0.5)
        if panel.xticks:
            ax.set_xticks(panel.xticks)
            ax.set_xticklabels(panel.xticklabels)
        ax.set_xlabel(panel.xlabel)
        ax.set_ylabel(panel.ylabel)
        if panel.title:
            ax.set_title(panel.title)
        if panel.series:
            ax.legend(fontsize=panel.legend_size)
        ax.grid(alpha=0.3)
    if spec.suptitle:
        fig.suptitle(spec.suptitle)
    return fig


def draw_figures(
    specs: Iterable[FigureSpec],
    *,
    save_dir: Path | None = None,
    show: Callable[[], None] | None = None,
) -> tuple[Path, ...]:
    """Draw ``specs``; optionally save PNGs into ``save_dir`` and call ``show``.

    Without ``show`` the non-interactive Agg backend is forced. All figures are
    open when ``show`` is called and closed before this function returns.

    Returns:
        The saved file paths (empty if ``save_dir`` is None).

    Raises:
        PlottingUnavailableError: if matplotlib is not installed.
    """
    plt = _pyplot(interactive=show is not None)
    figures = []
    saved = []
    try:
        for spec in specs:
            fig = _draw(plt, spec)
            figures.append(fig)
            if save_dir is not None:
                save_dir.mkdir(parents=True, exist_ok=True)
                path = save_dir / spec.filename
                fig.savefig(path, dpi=_DPI, bbox_inches="tight")
                saved.append(path)
        if show is not None:
            show()
    finally:
        for fig in figures:
            plt.close(fig)
    return tuple(saved)
