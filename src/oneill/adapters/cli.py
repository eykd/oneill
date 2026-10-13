"""Command-line adapter: parses arguments and writes to an injected stream."""

import argparse
import sys
from typing import TYPE_CHECKING, TextIO

from oneill import __version__

if TYPE_CHECKING:
    from collections.abc import Sequence


def build_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    return argparse.ArgumentParser(
        prog="oneill",
        description="Fictitious forces inside a rotating O'Neill cylinder",
    )


def main(argv: Sequence[str] | None = None, stdout: TextIO | None = None) -> int:
    """Run the CLI and return the process exit status."""
    out = sys.stdout if stdout is None else stdout
    build_parser().parse_args(argv)
    out.write(f"oneill {__version__}\n")
    return 0


def run() -> None:
    """Console-script entry point."""
    sys.exit(main())
