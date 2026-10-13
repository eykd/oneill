import re
from pathlib import Path
from typing import TYPE_CHECKING

import oneill

if TYPE_CHECKING:
    import pytest

README = Path(__file__).parent.parent / "README.md"


def test_every_exported_name_exists_and_is_unique() -> None:
    assert len(set(oneill.__all__)) == len(oneill.__all__)
    for name in oneill.__all__:
        assert hasattr(oneill, name), name


def test_all_is_sorted_constants_first_then_alphabetical() -> None:
    names = oneill.__all__
    assert names == sorted(names, key=lambda n: (not n.isupper(), n))


def test_version_is_a_string() -> None:
    assert re.fullmatch(r"\d+\.\d+\.\d+.*", oneill.__version__)


def test_readme_library_example_runs_and_prints_what_it_claims(
    capsys: pytest.CaptureFixture[str],
) -> None:
    text = README.read_text()
    block = re.search(r"```python\n(.*?)```", text, re.DOTALL)
    assert block is not None
    code = block.group(1)
    exec(compile(code, "README.md", "exec"), {})  # noqa: S102
    out = capsys.readouterr().out
    assert "10.36" in out
    assert "4.62" in out
    assert "-17.23" in out
