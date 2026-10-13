import math
import re
from typing import TYPE_CHECKING

from hypothesis import given
from hypothesis import strategies as st

from oneill.adapters.text import (
    fmt3,
    render_ballistics,
    render_custom,
    render_elevator,
    render_lateral,
    render_selftest,
    render_stairs,
    render_summary,
    render_vertical,
    section,
)
from oneill.core.vectors import LocalVector
from oneill.motion.scenarios import (
    BallisticCase,
    ballistics_report,
    custom_report,
    elevator_report,
    lateral_report,
    stairs_report,
    vertical_report,
)
from oneill.motion.verification import CheckResult, run_selftest

if TYPE_CHECKING:
    from oneill.core.cylinder import Cylinder


class TestFmt3:
    def test_known_vector(self) -> None:
        assert fmt3([1.5, -2.25, 0.0]) == "(  +1.500,  -2.250,  +0.000)"

    def test_values_that_round_to_zero_never_print_as_negative_zero(self) -> None:
        assert fmt3([0.04, -0.04, 0.06], p=1, w=5) == "( +0.0, +0.0, +0.1)"
        assert "-0.000" not in fmt3([-1e-9, -0.0, 1.0])

    @given(st.lists(st.floats(-999.0, 999.0), min_size=1, max_size=4))
    def test_fixed_width_and_sign_property(self, values: list[float]) -> None:
        text = fmt3(values)
        fields = text[1:-1].split(",")
        assert len(fields) == len(values)
        assert all(len(f) == 8 and f.strip()[0] in "+-" for f in fields)
        assert "-0.000" not in text


def test_section_is_a_banner_with_a_leading_blank_line() -> None:
    assert section("TITLE") == "\n" + "=" * 110 + "\nTITLE\n" + "=" * 110 + "\n"


class TestSummary:
    def test_mentions_key_figures(self, default_cyl: Cylinder) -> None:
        text = render_summary(default_cyl.summary())
        for token in (
            "O'NEILL CYLINDER PARAMETERS",
            "3,200.0 m",
            "20,000.0 m",
            "0.5286 rpm",
            "period 113.5 s",
            "177.1 m/s (638 km/h)",
            "1.0000 g",
            "Coriolis scale",
            "weightless at        177.1 m/s antispinward",
        ):
            assert token in text
        assert "spin acceleration" not in text

    def test_spin_up_line_appears_only_when_relevant(
        self, spinup_cyl: Cylinder
    ) -> None:
        assert "spin acceleration    +1.000e-04 rad/s^2" in render_summary(
            spinup_cyl.summary()
        )


class TestLateral:
    def test_golden_values(self, default_cyl: Cylinder) -> None:
        out = render_lateral(lateral_report(default_cyl, [1.4], [90, 270], 70.0))
        for token in (
            "101.587",
            "98.426",
            "-10.85",
            "+10.85",
            "spinward",
            "antispinward",
            "TAKEAWAY",
            "1. MOVING LATERALLY ON A LEVEL FLOOR (height 0 m)",
            "Forces for a 70 kg person.",
        ):
            assert token in out

    def test_row_layout_is_exact(self, default_cyl: Cylinder) -> None:
        out = render_lateral(lateral_report(default_cyl, [1.4], [90], 70.0))
        row = next(line for line in out.splitlines() if "101.587" in line)
        assert row == (
            "  90  spinward           1.40  (  +0.000,  +0.000,  -1.581)  "
            "(   +0.00,   +0.00,  -10.85)   0.006   101.587   0.00"
        )

    def test_group_per_speed_ends_with_a_blank_line(
        self, default_cyl: Cylinder
    ) -> None:
        out = render_lateral(lateral_report(default_cyl, [1.0, 2.0], [0, 90], 70.0))
        assert out.count("\n\n") == 3  # banner, notes->table, two groups' separators


class TestVertical:
    def test_golden_values_and_height_filtering(self, default_cyl: Cylinder) -> None:
        heights = [0.0, -5.0, default_cyl.radius, default_cyl.radius + 10.0]
        out = render_vertical(vertical_report(default_cyl, [1.0], heights, 70.0))
        assert "+1.129" in out
        assert "-1.129" in out
        assert out.count("skipping height") == 3
        assert "(skipping height -5.0: outside 0..radius)" in out
        assert "spinward (" in out
        assert "antispinward (" in out
        assert "TAKEAWAY" in out


class TestStairs:
    def test_prints_one_row_per_heading_for_both_directions(
        self, default_cyl: Cylinder
    ) -> None:
        out = render_stairs(
            stairs_report(default_cyl, 0.7, 35.0, [0, 90, 180, 270], 70.0)
        )
        assert out.count(" up ") >= 4
        assert out.count(" down ") >= 4
        assert "TAKEAWAY" in out
        assert (
            "3. STAIRS: 0.7 m/s along a 35 deg straight staircase (evaluated at 3 m)"
            in out
        )

    def test_lateral_column_is_a_clean_two_vector(self, default_cyl: Cylinder) -> None:
        out = render_stairs(stairs_report(default_cyl, 0.7, 35.0, [90], 70.0))
        row = next(line for line in out.splitlines() if " up " in line)
        assert re.search(r"\(\s*\+\d\.\d{3},\s*[+-]\d\.\d{3}\)  \s*\d\.\d\d$", row)
        assert " )" not in out


class TestElevator:
    def test_both_rides_are_tabulated_with_the_duration(
        self, default_cyl: Cylinder
    ) -> None:
        out = render_elevator(elevator_report(default_cyl, 3000.0, 10.0, 1.0, 70.0))
        assert "UP ride, duration 310.0 s" in out
        assert "DOWN ride, duration 310.0 s" in out
        assert "max tilt" in out
        assert "(v_max 10 m/s, a_max 1 m/s^2, 70 kg rider)" in out
        assert "felt weight range" in out


class TestBallistics:
    def test_all_default_cases_appear(self, default_cyl: Cylinder) -> None:
        out = render_ballistics(ballistics_report(default_cyl))
        for label in (
            "drop from hand (1.5 m)",
            "drop from 1000 m",
            "throw 15 m/s 45deg +axial",
        ):
            assert label in out
        assert "TAKEAWAY" in out
        assert "[hits endcap]" not in out

    def test_endcap_flag_is_shown(self, default_cyl: Cylinder) -> None:
        case = BallisticCase("long throw", 1.5, LocalVector(0.0, 5000.0, 10.0))
        out = render_ballistics(ballistics_report(default_cyl, [case]))
        assert out.rstrip().splitlines()[-3].endswith("[hits endcap]")


class TestCustom:
    def test_one_block_per_velocity(self, default_cyl: Cylinder) -> None:
        velocities = [LocalVector(2.0, 1.0, 0.5), LocalVector(-3.0, 0.0, 0.0)]
        out = render_custom(custom_report(default_cyl, velocities, 0.0, 70.0))
        assert out.count("APPARENT gravity") == 2
        assert out.count("support force") == 2
        assert "6. CUSTOM VELOCITIES at height 0 m" in out
        assert "velocity (   +2.00,   +1.00,   +0.50) m/s" in out


class TestSelfTest:
    def test_passing_run(self, default_cyl: Cylinder) -> None:
        out = render_selftest(run_selftest(default_cyl))
        assert "SELF-TEST" in out
        assert out.count("[PASS]") == 10
        assert "[FAIL]" not in out
        assert "ALL PASSED" in out

    def test_failing_run(self) -> None:
        out = render_selftest([CheckResult("bad", 2.0, 1.0, 0.1)])
        assert "[FAIL] bad: got 2, expected 1" in out
        assert out.rstrip().endswith("SOME CHECKS FAILED")

    def test_numbers_use_seven_significant_digits(self) -> None:
        out = render_selftest([CheckResult("pi", math.pi, math.pi, 0.0)])
        assert "got 3.141593, expected 3.141593" in out
