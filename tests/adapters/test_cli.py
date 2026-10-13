"""Tests for the command-line adapter."""

import runpy
from io import StringIO

import pytest

from oneill.adapters import cli
from oneill.adapters.cli import main, run
from oneill.motion.verification import CheckResult


def run_cli(*argv: str) -> tuple[int, str]:
    out = StringIO()
    status = main(list(argv), out)
    return status, out.getvalue()


class TestReports:
    def test_default_run_prints_every_section(self) -> None:
        status, out = run_cli()
        assert status == 0
        for token in (
            "O'NEILL CYLINDER PARAMETERS",
            "1. MOVING LATERALLY",
            "2. ASCENDING",
            "3. STAIRS",
            "4. ELEVATOR RIDE",
            "5. FREE FLIGHT",
        ):
            assert token in out
        assert "6. CUSTOM" not in out

    def test_default_elevator_goes_to_15_16_of_the_radius(self) -> None:
        _, out = run_cli("--sections", "elevator")
        assert "4. ELEVATOR RIDE floor <-> 3000 m" in out

    def test_default_vertical_heights_are_radius_fractions(self) -> None:
        _, out = run_cli("--sections", "vertical")
        assert "(skipping height" not in out
        for height in ("0", "500", "1600", "2800"):
            assert f"\n{height:>7} " in out

    def test_explicit_heights_replace_the_defaults(self) -> None:
        _, out = run_cli("--sections", "vertical", "--heights", "0", "123")
        assert "\n    123 " in out
        assert "\n   2800 " not in out

    def test_section_selection(self) -> None:
        _, out = run_cli("--sections", "vertical")
        assert "2. ASCENDING" in out
        assert "1. MOVING LATERALLY" not in out
        assert "5. FREE FLIGHT" not in out

    def test_several_sections_keep_canonical_order(self) -> None:
        _, out = run_cli("--sections", "ballistics", "lateral")
        assert out.index("1. MOVING LATERALLY") < out.index("5. FREE FLIGHT")

    def test_rpm_option_overrides_gravity(self) -> None:
        _, out = run_cli("--rpm", "0.5", "--sections", "lateral")
        assert "0.5000 rpm" in out

    def test_rpm_and_gravity_are_mutually_exclusive(self) -> None:
        with pytest.raises(SystemExit) as exc:
            run_cli("--rpm", "0.5", "--gravity", "9.81")
        assert exc.value.code == 2

    def test_spin_up_shows_the_euler_line(self) -> None:
        _, out = run_cli("--spin-accel", "1e-4", "--sections", "lateral")
        assert "spin acceleration" in out

    def test_custom_velocities_accept_negative_components(self) -> None:
        _, out = run_cli(
            "--sections", "lateral",
            "--velocity", "2", "1", "0.5",
            "--velocity", "-3", "0", "0",
        )  # fmt: skip
        assert "6. CUSTOM VELOCITIES" in out
        assert out.count("APPARENT gravity") == 2

    def test_custom_velocities_honour_at_height(self) -> None:
        _, out = run_cli(
            "--sections", "lateral", "--velocity", "1", "0", "0", "--at-height", "50"
        )
        assert "6. CUSTOM VELOCITIES at height 50 m" in out

    @pytest.mark.parametrize(
        "argv",
        [
            ["--radius", "100", "--length", "300"],
            ["--radius", "1000"],
            ["--radius", "8000", "--length", "40000"],
            ["--radius", "1000", "--gravity", "3.71", "--mass", "90"],
        ],
    )
    def test_full_report_works_for_other_habitat_sizes(self, argv: list[str]) -> None:
        status, out = run_cli(*argv)
        assert status == 0
        assert "4. ELEVATOR RIDE" in out
        assert "5. FREE FLIGHT" in out

    def test_small_habitat_drops_ballistic_cases_that_do_not_fit(self) -> None:
        _, out = run_cli("--radius", "500", "--sections", "ballistics")
        assert "drop from 100 m tower" in out
        assert "drop from 1000 m" not in out

    def test_default_output_is_independent_of_how_it_is_captured(
        self, capsys: pytest.CaptureFixture[str]
    ) -> None:
        assert main([]) == 0
        captured = capsys.readouterr().out
        assert captured == run_cli()[1]


class TestSelfTest:
    def test_selftest_passes_and_skips_the_report_sections(self) -> None:
        status, out = run_cli("--selftest")
        assert status == 0
        assert "ALL PASSED" in out
        assert "O'NEILL CYLINDER PARAMETERS" in out
        assert "1. MOVING LATERALLY" not in out

    def test_failing_selftest_exits_non_zero(
        self, monkeypatch: pytest.MonkeyPatch
    ) -> None:
        monkeypatch.setattr(
            cli, "run_selftest", lambda _cyl: (CheckResult("bad", 2.0, 1.0, 0.1),)
        )
        status, out = run_cli("--selftest")
        assert status == 1
        assert "SOME CHECKS FAILED" in out


class TestErrors:
    @pytest.mark.parametrize(
        ("argv", "message"),
        [
            (["--radius", "-5"], "radius"),
            (["--length", "0"], "length"),
            (["--gravity", "-1"], "surface_g"),
            (["--rpm", "0"], "omega"),
            (["--elevator-top", "99999", "--sections", "elevator"], "elevator top"),
            (["--elevator-vmax", "0", "--sections", "elevator"], "v_max"),
            (["--mass", "nan"], "mass"),
            (["--at-height", "5000", "--velocity", "1", "0", "0"], "height"),
        ],
    )
    def test_domain_errors_are_usage_errors_and_leave_stdout_empty(
        self,
        argv: list[str],
        message: str,
        capsys: pytest.CaptureFixture[str],
    ) -> None:
        out = StringIO()
        with pytest.raises(SystemExit) as exc:
            main(argv, out)
        assert exc.value.code == 2
        assert out.getvalue() == ""
        assert message in capsys.readouterr().err

    def test_unknown_option_exits_with_usage_error(self) -> None:
        with pytest.raises(SystemExit) as exc:
            main(["--bogus"], StringIO())
        assert exc.value.code == 2


class TestEntryPoints:
    def test_run_exits_with_main_status(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr("sys.argv", ["oneill", "--sections", "lateral"])
        with pytest.raises(SystemExit) as exc:
            run()
        assert exc.value.code == 0

    def test_python_dash_m_entry_point(self, monkeypatch: pytest.MonkeyPatch) -> None:
        monkeypatch.setattr("sys.argv", ["oneill", "--sections", "lateral"])
        with pytest.raises(SystemExit) as exc:
            runpy.run_module("oneill", run_name="__main__")
        assert exc.value.code == 0
