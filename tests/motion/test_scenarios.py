import math

import pytest
from hypothesis import given
from hypothesis import strategies as st

from oneill.core.cylinder import G_STANDARD, Cylinder
from oneill.core.vectors import LocalVector
from oneill.motion.scenarios import (
    BALLISTIC_CASES,
    BallisticCase,
    ballistic_cases,
    ballistics_report,
    custom_report,
    elevator_report,
    lateral_report,
    stairs_report,
    vertical_report,
)


class TestLateral:
    def test_golden_values_for_walking_1_4_mps(self, default_cyl: Cylinder) -> None:
        """1.4 m/s: 101.587 % weight spinward, 98.426 % antispinward, 10.85 N."""
        report = lateral_report(default_cyl, [1.4], [90, 270], 70.0)
        spin, anti = report.groups[0]
        assert spin.felt_weight_pct == pytest.approx(101.587, abs=5e-4)
        assert anti.felt_weight_pct == pytest.approx(98.426, abs=5e-4)
        assert spin.coriolis_force_n.up == pytest.approx(-10.85, abs=5e-3)
        assert anti.coriolis_force_n.up == pytest.approx(+10.85, abs=5e-3)
        assert (spin.direction, anti.direction) == ("spinward", "antispinward")

    def test_structure_groups_by_speed(self, default_cyl: Cylinder) -> None:
        report = lateral_report(default_cyl, [1.0, 2.0, 3.0], [0, 90], 70.0, height=5.0)
        assert [len(g) for g in report.groups] == [2, 2, 2]
        assert [g[0].speed for g in report.groups] == [1.0, 2.0, 3.0]
        assert report.height == 5.0
        assert report.mass == 70.0
        assert report.surface_gravity == default_cyl.surface_gravity
        assert report.coriolis_pct_per_mps == pytest.approx(
            100.0 * 2.0 * default_cyl.omega / G_STANDARD
        )

    @given(st.floats(-30.0, 30.0), st.floats(0.0, 359.0))
    def test_row_matches_the_inertial_frame_weight(
        self, speed: float, heading: float
    ) -> None:
        cyl = Cylinder()
        row = lateral_report(cyl, [speed], [heading], 70.0).groups[0][0]
        v_spin = speed * math.sin(math.radians(heading))
        expected = 100.0 * (cyl.rim_speed + v_spin) ** 2 / cyl.radius / G_STANDARD
        assert row.felt_weight_pct == pytest.approx(expected, rel=1e-9)
        assert row.tilt_deg == pytest.approx(0.0, abs=1e-6)
        assert row.coriolis_pct.spinward == 0.0
        assert row.coriolis_pct.axial == 0.0
        assert row.curvature_pct == pytest.approx(
            100.0 * v_spin**2 / cyl.radius / G_STANDARD, abs=1e-9
        )


class TestVertical:
    def test_golden_sideways_push_and_skipped_heights(
        self, default_cyl: Cylinder
    ) -> None:
        heights = [0.0, -5.0, default_cyl.radius, default_cyl.radius + 10.0]
        report = vertical_report(default_cyl, [1.0], heights, 70.0)
        assert [g.skipped for g in report.groups] == [False, True, True, True]
        up, down = report.groups[0].rows
        assert up.coriolis_pct.spinward == pytest.approx(1.129, abs=5e-4)
        assert down.coriolis_pct.spinward == pytest.approx(-1.129, abs=5e-4)
        assert (up.push, down.push) == ("spinward", "antispinward")
        assert [g.height for g in report.groups] == heights

    def test_rows_cover_ascending_then_descending_rates(
        self, default_cyl: Cylinder
    ) -> None:
        group = vertical_report(default_cyl, [0.5, 2.0], [100.0], 80.0).groups[0]
        assert [r.rate for r in group.rows] == [0.5, 2.0, -0.5, -2.0]

    @given(st.floats(0.1, 20.0), st.floats(0.0, 0.95))
    def test_tilt_and_local_gravity_follow_the_closed_form(
        self, rate: float, frac: float
    ) -> None:
        cyl = Cylinder()
        h = frac * cyl.radius
        row = vertical_report(cyl, [rate], [h], 70.0).groups[0].rows[0]
        local_g = cyl.gravity_at_height(h)
        assert row.local_g_pct == pytest.approx(100.0 * local_g / G_STANDARD)
        assert row.felt_weight_pct == pytest.approx(row.local_g_pct, rel=1e-9)
        expected = math.degrees(math.atan2(2.0 * cyl.omega * rate, local_g))
        assert row.tilt_deg == pytest.approx(expected, abs=1e-7)
        assert row.push_force_n == pytest.approx(
            70.0 * math.hypot(2.0 * cyl.omega * rate, 0.0), rel=1e-9
        )


class TestStairs:
    def test_up_then_down_groups_with_one_row_per_heading(
        self, default_cyl: Cylinder
    ) -> None:
        report = stairs_report(default_cyl, 0.7, 35.0, [0, 90, 180, 270], 70.0)
        up, down = report.groups
        assert [r.sense for r in up] == ["up"] * 4
        assert [r.sense for r in down] == ["down"] * 4
        assert [r.direction for r in up] == [
            "+axial",
            "spinward",
            "-axial",
            "antispinward",
        ]
        assert (report.speed, report.slope_deg, report.height) == (0.7, 35.0, 3.0)

    def test_values_match_the_superposition_of_level_and_vertical_motion(
        self, default_cyl: Cylinder
    ) -> None:
        up_spinward = stairs_report(default_cyl, 0.7, 35.0, [90], 70.0).groups[0][0]
        e = math.radians(35.0)
        v_spin, v_up = 0.7 * math.cos(e), 0.7 * math.sin(e)
        w = default_cyl.omega
        assert up_spinward.velocity.spinward == pytest.approx(v_spin)
        assert up_spinward.velocity.up == pytest.approx(v_up)
        assert up_spinward.lateral_pct[0] == pytest.approx(
            100.0 * 2.0 * w * v_up / G_STANDARD, rel=1e-6
        )
        assert up_spinward.lateral_pct[1] == pytest.approx(0.0, abs=1e-9)
        local_g = default_cyl.gravity_at_height(3.0)
        assert up_spinward.felt_weight_pct == pytest.approx(
            100.0 * (local_g + 2.0 * w * v_spin) / G_STANDARD, rel=1e-9
        )


class TestElevator:
    def test_rides_have_the_expected_duration_and_directions(
        self, default_cyl: Cylinder
    ) -> None:
        report = elevator_report(default_cyl, 3000.0, 10.0, 1.0, 70.0)
        assert report.duration == pytest.approx(310.0)
        assert [r.label for r in report.rides] == ["UP", "DOWN"]
        up, down = report.rides
        assert len(up.rows) == len(down.rows) == 15
        assert up.rows[0].height < down.rows[0].height

    def test_rows_are_sampled_just_inside_the_ride(self, default_cyl: Cylinder) -> None:
        up = elevator_report(default_cyl, 3000.0, 10.0, 1.0, 70.0).rides[0]
        assert up.rows[0].t == pytest.approx(0.05)
        assert up.rows[-1].t == pytest.approx(310.0 - 0.05)
        assert up.rows[0].a_up == pytest.approx(1.0)
        assert up.rows[-1].a_up == pytest.approx(-1.0)

    def test_summary_statistics_bound_the_rows(self, default_cyl: Cylinder) -> None:
        up = elevator_report(default_cyl, 3000.0, 10.0, 1.0, 70.0).rides[0]
        assert up.max_tilt_deg >= max(r.tilt_deg for r in up.rows) - 1e-9
        assert up.max_tilt_height > 2900.0
        assert up.weight_min_pct <= min(r.felt_weight_pct for r in up.rows)
        assert up.weight_max_pct >= max(r.felt_weight_pct for r in up.rows)
        assert up.max_side_force_n >= max(abs(r.side_force_n) for r in up.rows) - 1e-6

    def test_short_ride_that_never_reaches_cruise_speed(
        self, default_cyl: Cylinder
    ) -> None:
        report = elevator_report(default_cyl, 50.0, 10.0, 1.0, 70.0, rows=5)
        assert report.duration == pytest.approx(2.0 * math.sqrt(50.0))
        assert len(report.rides[0].rows) == 5

    @pytest.mark.parametrize("top", [0.0, -1.0, 3200.0, 4000.0])
    def test_rejects_impossible_shafts(self, default_cyl: Cylinder, top: float) -> None:
        with pytest.raises(ValueError, match="elevator top"):
            elevator_report(default_cyl, top, 10.0, 1.0, 70.0)


class TestBallistics:
    def test_default_cases_all_fit_the_default_habitat(
        self, default_cyl: Cylinder
    ) -> None:
        assert ballistic_cases(default_cyl) == BALLISTIC_CASES
        report = ballistics_report(default_cyl)
        assert [r.label for r in report.rows] == [c.label for c in BALLISTIC_CASES]

    def test_cases_that_do_not_fit_are_left_out(self) -> None:
        small = Cylinder.from_gravity(radius=500.0, length=2000.0)
        labels = [c.label for c in ballistic_cases(small)]
        assert "drop from 100 m tower" in labels
        assert "drop from 1000 m" not in labels
        assert len(ballistics_report(small).rows) == len(labels)

    def test_rows_hold_exact_naive_and_deviation(self, default_cyl: Cylinder) -> None:
        row = ballistics_report(default_cyl).rows[0]
        assert row.landing[0] < 0.0
        assert row.deviation == (
            row.landing[0] - row.naive_landing[0],
            row.landing[1] - row.naive_landing[1],
        )
        assert row.flight_time > row.naive_time
        assert row.hits_endcap is False

    def test_custom_cases_and_endcap_flag(self, default_cyl: Cylinder) -> None:
        case = BallisticCase("long throw", 1.5, LocalVector(0.0, 5000.0, 10.0))
        row = ballistics_report(default_cyl, [case]).rows[0]
        assert row.label == "long throw"
        assert row.hits_endcap is True


class TestCustom:
    def test_breakdown_is_internally_consistent(self, default_cyl: Cylinder) -> None:
        report = custom_report(
            default_cyl,
            [LocalVector(2.0, 1.0, 0.5), LocalVector(-3.0, 0.0, 0.0)],
            0.0,
            70.0,
        )
        assert len(report.entries) == 2
        e = report.entries[0]
        parts = (
            e.centrifugal_pct.up
            + e.coriolis_pct.up
            + e.euler_pct.up
            + e.neg_path_accel_pct.up
        )
        assert e.apparent_pct.up == pytest.approx(parts, abs=1e-9)
        assert e.felt_weight_pct == pytest.approx(-e.apparent_pct.up)
        assert e.support_force_n.up == pytest.approx(
            -70.0 * e.apparent_pct.up / 100.0 * G_STANDARD, rel=1e-9
        )
        assert e.coriolis_force_n.spinward == pytest.approx(
            70.0 * 2.0 * default_cyl.omega * 0.5, rel=1e-9
        )
        assert report.height == 0.0
        assert report.mass == 70.0

    def test_euler_appears_only_while_spinning_up(self, spinup_cyl: Cylinder) -> None:
        entry = custom_report(spinup_cyl, [LocalVector(0.0, 0.0, 0.0)], 0.0, 70.0)
        assert entry.entries[0].euler_pct.spinward < 0.0


@pytest.mark.parametrize(
    "builder",
    [
        lambda c: lateral_report(c, [1.0], [0.0], math.nan),
        lambda c: vertical_report(c, [1.0], [0.0], math.inf),
        lambda c: stairs_report(c, 0.7, 35.0, [0.0], math.nan),
        lambda c: elevator_report(c, 100.0, 10.0, 1.0, math.nan),
        lambda c: custom_report(c, [], 0.0, math.nan),
    ],
)
def test_every_builder_rejects_a_non_finite_mass(
    builder: object, default_cyl: Cylinder
) -> None:
    with pytest.raises(ValueError, match="mass"):
        builder(default_cyl)  # type: ignore[operator]
