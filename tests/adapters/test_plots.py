import sys
from typing import TYPE_CHECKING

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import numpy as np
import pytest

from oneill.adapters import plots
from oneill.adapters.plots import (
    PlottingUnavailableError,
    check_available,
    draw_figures,
    show_all,
)
from oneill.motion.figures import FigureSpec, Panel, Series, build_figures

if TYPE_CHECKING:
    from pathlib import Path

    from oneill.core.cylinder import Cylinder

NAMES = [
    "1_floor_weight_vs_heading.png",
    "2_climbing_vs_height.png",
    "3_elevator.png",
    "4_ballistics.png",
]


@pytest.fixture
def specs(default_cyl: Cylinder) -> tuple[FigureSpec, ...]:
    return build_figures(default_cyl, (1.4, 6.0), (0.3, 3.0), 3000.0, 10.0, 1.0)


@pytest.fixture(autouse=True)
def _close_figures() -> None:
    plt.close("all")


def ink_fraction(path: Path) -> float:
    """Fraction of pixels that are neither white nor transparent."""
    img = mpimg.imread(path)
    rgb = img[..., :3]
    return float(np.mean(np.any(rgb < 0.97, axis=-1)))


def test_saves_four_non_blank_pngs(
    tmp_path: Path, specs: tuple[FigureSpec, ...]
) -> None:
    saved = draw_figures(specs, save_dir=tmp_path / "out" / "nested")
    assert [p.name for p in saved] == NAMES
    for path in saved:
        assert path.read_bytes().startswith(b"\x89PNG")
        img = mpimg.imread(path)
        assert img.shape[0] > 200
        assert img.shape[1] > 400
        assert ink_fraction(path) > 0.01
        assert len(np.unique(img.reshape(-1, img.shape[-1]), axis=0)) > 20


def test_without_save_dir_or_show_nothing_is_written(
    tmp_path: Path, specs: tuple[FigureSpec, ...], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    assert draw_figures(specs) == ()
    assert list(tmp_path.iterdir()) == []


def test_non_interactive_use_forces_the_agg_backend(
    specs: tuple[FigureSpec, ...], monkeypatch: pytest.MonkeyPatch
) -> None:
    import matplotlib as mpl  # noqa: PLC0415

    used: list[str] = []
    monkeypatch.setattr(mpl, "use", used.append)
    draw_figures(specs)
    assert used == ["Agg"]


def test_show_callback_is_called_once_after_drawing_and_keeps_the_backend(
    specs: tuple[FigureSpec, ...], monkeypatch: pytest.MonkeyPatch
) -> None:
    import matplotlib as mpl  # noqa: PLC0415

    used: list[str] = []
    monkeypatch.setattr(mpl, "use", used.append)
    seen: list[int] = []
    draw_figures(specs, show=lambda: seen.append(len(plt.get_fignums())))
    assert used == []
    assert seen == [4]


def test_figures_are_closed_afterwards(specs: tuple[FigureSpec, ...]) -> None:
    draw_figures(specs, show=lambda: None)
    assert plt.get_fignums() == []


def test_show_all_delegates_to_pyplot(monkeypatch: pytest.MonkeyPatch) -> None:
    calls: list[str] = []
    monkeypatch.setattr(plt, "show", lambda: calls.append("shown"))
    show_all()
    assert calls == ["shown"]


def test_single_panel_and_custom_styles_render(tmp_path: Path) -> None:
    ts = np.linspace(0.0, 1.0, 20)
    spec = FigureSpec(
        "one.png",
        (
            Panel(
                "x",
                "y",
                (Series("a", ts, ts, "k--"), Series("b", ts, ts**2)),
                title="t",
                hlines=(0.5,),
                xticks=(0.0, 1.0),
                xticklabels=("lo", "hi"),
            ),
        ),
        (1, 1),
        (4.0, 3.0),
    )
    (path,) = draw_figures([spec], save_dir=tmp_path)
    assert ink_fraction(path) > 0.01


class TestAvailability:
    def test_check_passes_when_matplotlib_is_installed(self) -> None:
        check_available()

    def test_missing_matplotlib_raises_a_helpful_error(
        self, monkeypatch: pytest.MonkeyPatch, specs: tuple[FigureSpec, ...]
    ) -> None:
        monkeypatch.setitem(sys.modules, "matplotlib", None)
        with pytest.raises(PlottingUnavailableError, match=r"oneill\[plot\]"):
            check_available()
        with pytest.raises(PlottingUnavailableError, match="matplotlib"):
            draw_figures(specs)

    def test_importing_the_module_does_not_import_matplotlib(self) -> None:
        assert "matplotlib" not in vars(plots)


def test_a_panel_without_series_still_renders(tmp_path: Path) -> None:
    spec = FigureSpec("empty.png", (Panel("x", "y", ()),), (1, 1), (3.0, 2.0))
    (path,) = draw_figures([spec], save_dir=tmp_path)
    assert path.exists()
