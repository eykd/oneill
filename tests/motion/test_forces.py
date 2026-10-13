import math
from types import SimpleNamespace

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.core.cylinder import G_STANDARD, AxisSingularityError, Cylinder
from oneill.core.vectors import LocalVector
from oneill.motion.forces import ForceHistory, analyze_path
from oneill.motion.paths import (
    FloorWalk,
    NumericPath,
    Path,
    StraightPath,
    velocity_path,
)
from oneill.motion.profiles import ConstantSpeed, Trapezoid

ACC_TOL = 1e-9  # m/s^2: paths are analytic, so only round-off remains
speeds = st.floats(-60.0, 60.0)
velocities = st.tuples(
    st.floats(-30.0, 30.0), st.floats(-30.0, 30.0), st.floats(-30.0, 30.0)
)
FREE_FALL_TOL = 1e-6  # numeric path (5-point stencil, step 1e-2)


def probe(cyl: Cylinder, path: Path, t: float = 0.0) -> ForceHistory:
    """Analyse ``path`` at the single instant ``t``."""
    return analyze_path(cyl, path, [t])


def straight(
    cyl: Cylinder,
    heading: float,
    elevation: float,
    speed: float,
    start_height: float = 0.0,
) -> StraightPath:
    return StraightPath(
        cyl, heading, elevation, ConstantSpeed(speed), 1.0, start_height=start_height
    )


@pytest.fixture
def elevator(default_cyl: Cylinder) -> SimpleNamespace:
    """A floor-to-3000 m ride (v_max 10 m/s, a_max 1 m/s^2) in the default habitat."""
    prof = Trapezoid(3000.0, 10.0, 1.0)
    return SimpleNamespace(
        profile=prof,
        duration=prof.duration,
        up=StraightPath(default_cyl, 0.0, 90.0, prof, prof.duration),
        down=StraightPath(
            default_cyl, 0.0, -90.0, prof, prof.duration, start_height=3000.0
        ),
    )


class TestStandingStill:
    @pytest.mark.parametrize("height_frac", [0.0, 0.1, 0.5, 0.9, 0.999])
    def test_only_centrifugal_force_acts(
        self, cyl: Cylinder, height_frac: float
    ) -> None:
        h = height_frac * cyl.radius
        hist = probe(cyl, FloorWalk(cyl, 0.0, 0.0, height=h))
        np.testing.assert_allclose(
            hist.centrifugal[0],
            [0.0, 0.0, -cyl.gravity_at_height(h)],
            rtol=1e-9,
            atol=1e-12,
        )
        np.testing.assert_allclose(hist.coriolis[0], 0.0, atol=1e-12)
        np.testing.assert_allclose(hist.euler[0], 0.0, atol=0.0)
        np.testing.assert_allclose(hist.acc[0], 0.0, atol=1e-12)
        assert hist.felt_weight[0] == pytest.approx(cyl.gravity_at_height(h), rel=1e-9)
        assert hist.tilt_deg[0] == pytest.approx(0.0, abs=1e-9)
        assert hist.heights[0] == pytest.approx(h, abs=1e-9)

    @given(
        st.floats(-20.0, 20.0),
        st.floats(-1e4, 1e4),
        st.floats(0.0, 360.0),
        st.floats(0.5, 20.0),
    )
    @example(math.pi, 9999.0, 60.0, 3.0)
    @example(-2.5, -9000.0, 60.0, 3.0)
    def test_forces_do_not_depend_on_where_you_are(
        self, theta0: float, z0: float, heading: float, speed: float
    ) -> None:
        cyl = Cylinder()
        reference = probe(cyl, FloorWalk(cyl, heading, speed))
        hist = probe(cyl, FloorWalk(cyl, heading, speed, theta0=theta0, z0=z0))
        np.testing.assert_allclose(hist.apparent, reference.apparent, atol=1e-9)
        np.testing.assert_allclose(hist.coriolis, reference.coriolis, atol=1e-9)


class TestLevelMotion:
    @pytest.mark.parametrize("v_spin", [-50.0, -10.0, -1.4, 0.0, 1.4, 10.0, 50.0])
    def test_weight_matches_inertial_frame_centripetal_formula(
        self, cyl: Cylinder, v_spin: float
    ) -> None:
        hist = probe(cyl, FloorWalk(cyl, 90.0, v_spin))
        assert hist.felt_weight[0] == pytest.approx(
            (cyl.rim_speed + v_spin) ** 2 / cyl.radius, abs=ACC_TOL
        )

    @given(
        st.floats(10.0, 1.0e5),
        st.floats(0.05, 3.0),
        speeds,
        st.floats(0.0, 0.9),
    )
    @example(3200.0, 0.5, 6.0, 0.5)
    def test_weight_formula_holds_for_any_habitat_and_deck(
        self, radius: float, rpm: float, v_spin: float, frac: float
    ) -> None:
        cyl = Cylinder.from_rpm(radius, 1000.0, rpm)
        r = radius * (1.0 - frac)
        hist = probe(cyl, FloorWalk(cyl, 90.0, v_spin, height=frac * radius))
        expected = (cyl.omega * r + v_spin) ** 2 / r
        assert hist.felt_weight[0] == pytest.approx(expected, rel=1e-9, abs=1e-9)

    @pytest.mark.parametrize("heading", [0, 45, 90, 135, 180, 225, 270, 315, 30, 200.5])
    @pytest.mark.parametrize("speed", [1.4, 6.0])
    def test_coriolis_is_purely_vertical_and_set_by_spinward_speed(
        self, cyl: Cylinder, heading: float, speed: float
    ) -> None:
        hist = probe(cyl, FloorWalk(cyl, heading, speed))
        v_spin = speed * math.sin(math.radians(heading))
        np.testing.assert_allclose(
            hist.coriolis[0], [0.0, 0.0, -2.0 * cyl.omega * v_spin], atol=ACC_TOL
        )
        assert hist.coriolis[0, 1] == 0.0

    @pytest.mark.parametrize("heading", [0, 45, 90, 135, 180, 225, 270, 315])
    def test_no_sideways_push_and_no_tilt_when_walking(
        self, cyl: Cylinder, heading: float
    ) -> None:
        hist = probe(cyl, FloorWalk(cyl, heading, 6.0))
        assert hist.lateral_mag[0] == pytest.approx(0.0, abs=1e-9)
        assert hist.tilt_deg[0] == pytest.approx(0.0, abs=1e-7)

    @pytest.mark.parametrize("heading", [0.0, 180.0])
    @pytest.mark.parametrize("speed", [1.4, 12.0, 100.0])
    def test_axial_motion_changes_nothing(
        self, cyl: Cylinder, heading: float, speed: float
    ) -> None:
        hist = probe(cyl, FloorWalk(cyl, heading, speed))
        assert hist.felt_weight[0] == pytest.approx(cyl.surface_gravity, abs=ACC_TOL)
        np.testing.assert_allclose(hist.coriolis[0], 0.0, atol=1e-9)

    def test_spinward_antispinward_asymmetry(self, cyl: Cylinder) -> None:
        v = 6.0
        heavy = probe(cyl, FloorWalk(cyl, 90.0, v)).felt_weight[0]
        light = probe(cyl, FloorWalk(cyl, 270.0, v)).felt_weight[0]
        assert heavy > cyl.surface_gravity > light
        assert heavy - light == pytest.approx(4.0 * cyl.omega * v, abs=2 * ACC_TOL)
        assert heavy + light - 2.0 * cyl.surface_gravity == pytest.approx(
            2.0 * v * v / cyl.radius, abs=2 * ACC_TOL
        )

    def test_diagonal_walk_equals_its_spinward_component(self, cyl: Cylinder) -> None:
        diagonal = probe(cyl, FloorWalk(cyl, 45.0, 6.0)).felt_weight[0]
        spinward = probe(cyl, FloorWalk(cyl, 90.0, 6.0 * math.sqrt(0.5))).felt_weight[0]
        assert diagonal == pytest.approx(spinward, abs=ACC_TOL)

    @pytest.mark.parametrize("height_frac", [0.25, 0.5])
    @pytest.mark.parametrize("v_spin", [-6.0, -1.4, 1.4, 6.0])
    def test_upper_decks_use_the_local_radius(
        self, cyl: Cylinder, height_frac: float, v_spin: float
    ) -> None:
        r = cyl.radius * (1.0 - height_frac)
        hist = probe(cyl, FloorWalk(cyl, 90.0, v_spin, height=height_frac * cyl.radius))
        assert hist.felt_weight[0] == pytest.approx(
            (cyl.omega * r + v_spin) ** 2 / r, abs=ACC_TOL
        )

    def test_weightless_when_cancelling_the_spin(self, cyl: Cylinder) -> None:
        hist = probe(cyl, FloorWalk(cyl, 270.0, cyl.rim_speed))
        assert hist.felt_weight[0] == pytest.approx(0.0, abs=ACC_TOL)

    def test_twice_rim_speed_antispinward_restores_full_weight(
        self, cyl: Cylinder
    ) -> None:
        hist = probe(cyl, FloorWalk(cyl, 270.0, 2.0 * cyl.rim_speed))
        assert hist.felt_weight[0] == pytest.approx(cyl.surface_gravity, abs=ACC_TOL)


class TestVerticalMotion:
    @pytest.mark.parametrize("rate", [0.3, 1.0, 10.0, -0.3, -1.0, -10.0])
    @pytest.mark.parametrize("height_frac", [0.0, 0.5, 0.9])
    def test_coriolis_is_sideways_and_independent_of_height(
        self, cyl: Cylinder, rate: float, height_frac: float
    ) -> None:
        path = straight(
            cyl,
            0.0,
            math.copysign(90.0, rate),
            abs(rate),
            start_height=height_frac * cyl.radius,
        )
        hist = probe(cyl, path)
        np.testing.assert_allclose(
            hist.coriolis[0], [2.0 * cyl.omega * rate, 0.0, 0.0], atol=ACC_TOL
        )
        assert (hist.apparent[0, 0] > 0) == (rate > 0)

    @pytest.mark.parametrize("rate", [0.3, 10.0, -10.0])
    @pytest.mark.parametrize("height_frac", [0.0, 0.5, 0.9])
    def test_constant_speed_climbing_does_not_change_weight(
        self, cyl: Cylinder, rate: float, height_frac: float
    ) -> None:
        h = height_frac * cyl.radius
        path = straight(cyl, 0.0, math.copysign(90.0, rate), abs(rate), start_height=h)
        assert probe(cyl, path).felt_weight[0] == pytest.approx(
            cyl.gravity_at_height(h), abs=ACC_TOL
        )

    @pytest.mark.parametrize("rate", [1.0, 10.0])
    @pytest.mark.parametrize("height_frac", [0.0, 0.5, 0.9])
    def test_tilt_of_apparent_vertical(
        self, cyl: Cylinder, rate: float, height_frac: float
    ) -> None:
        r = cyl.radius * (1.0 - height_frac)
        path = straight(cyl, 0.0, 90.0, rate, start_height=cyl.radius - r)
        expected = math.degrees(math.atan2(2.0 * cyl.omega * rate, cyl.omega**2 * r))
        assert probe(cyl, path).tilt_deg[0] == pytest.approx(expected, abs=1e-7)

    def test_near_the_axis_the_sideways_push_dominates(self, cyl: Cylinder) -> None:
        path = straight(cyl, 0.0, 90.0, 10.0, start_height=cyl.radius - 1.0)
        hist = probe(cyl, path)
        assert hist.felt_weight[0] == pytest.approx(cyl.omega**2 * 1.0, abs=ACC_TOL)
        assert hist.tilt_deg[0] > 89.0


class TestCombinedMotion:
    @pytest.mark.parametrize(
        "velocity",
        [
            (1.0, 0.0, 1.0),
            (-2.0, 3.0, 0.5),
            (0.0, 1.4, -0.4),
            (5.0, -5.0, 5.0),
            (-0.57, 0.0, 0.4),
            (0.0, 0.0, -3.0),
        ],
    )
    def test_general_coriolis_formula(
        self, cyl: Cylinder, velocity: tuple[float, float, float]
    ) -> None:
        vs, _, vu = velocity
        path = velocity_path(cyl, LocalVector(*velocity), height=0.1 * cyl.radius)
        hist = probe(cyl, path)
        np.testing.assert_allclose(
            hist.coriolis[0],
            [2.0 * cyl.omega * vu, 0.0, -2.0 * cyl.omega * vs],
            atol=ACC_TOL,
        )

    @given(velocities, st.floats(0.0, 0.9))
    @example((1.0, 2.0, 0.5), 0.1)
    @example((0.0, 0.0, 4.0), 0.0)
    @example((-7.0, 1.0, -2.0), 0.5)
    def test_coriolis_is_perpendicular_to_velocity(
        self, velocity: tuple[float, float, float], frac: float
    ) -> None:
        cyl = Cylinder()
        path = velocity_path(cyl, LocalVector(*velocity), height=frac * cyl.radius)
        hist = probe(cyl, path)
        scale = 1.0 + float(np.dot(hist.vel[0], hist.vel[0]))
        assert float(np.dot(hist.coriolis[0], hist.vel[0])) == pytest.approx(
            0.0, abs=1e-9 * scale
        )

    @pytest.mark.parametrize("heading", [0.0, 45.0, 90.0, 200.0, 270.0])
    @pytest.mark.parametrize("slope", [35.0, -35.0, 60.0])
    def test_stairs_are_the_superposition_of_level_and_vertical_motion(
        self, cyl: Cylinder, heading: float, slope: float
    ) -> None:
        v = 0.7
        e = math.radians(slope)
        stairs = probe(cyl, straight(cyl, heading, slope, v, start_height=3.0))
        level = probe(cyl, straight(cyl, heading, 0.0, v * math.cos(e), 3.0))
        vertical = probe(
            cyl,
            straight(cyl, 0.0, math.copysign(90.0, slope), v * abs(math.sin(e)), 3.0),
        )
        np.testing.assert_allclose(
            stairs.coriolis[0],
            level.coriolis[0] + vertical.coriolis[0],
            atol=ACC_TOL,
        )

    @pytest.mark.parametrize(
        ("heading", "slope"),
        [
            (90.0, 35.0),
            (90.0, -35.0),
            (270.0, 35.0),
            (270.0, -35.0),
            (0.0, 35.0),
            (180.0, -35.0),
            (45.0, 35.0),
        ],
    )
    def test_what_you_feel_on_a_straight_staircase(
        self, cyl: Cylinder, heading: float, slope: float
    ) -> None:
        v, h0 = 0.7, 3.0
        e, hd = math.radians(slope), math.radians(heading)
        hist = probe(cyl, straight(cyl, heading, slope, v, start_height=h0))
        v_spin, v_up = v * math.cos(e) * math.sin(hd), v * math.sin(e)
        assert hist.felt_weight[0] - cyl.gravity_at_height(h0) == pytest.approx(
            2.0 * cyl.omega * v_spin, abs=ACC_TOL
        )
        assert hist.apparent[0, 0] == pytest.approx(2.0 * cyl.omega * v_up, abs=ACC_TOL)
        assert hist.apparent[0, 1] == pytest.approx(0.0, abs=ACC_TOL)


class TestEulerForce:
    @pytest.mark.parametrize("alpha", [1.0e-4, -1.0e-4])
    def test_standing_person_during_spin_change(self, alpha: float) -> None:
        c = Cylinder.from_gravity(spin_accel=alpha)
        times = np.array([0.0, 10.0, 100.0])
        hist = analyze_path(c, FloorWalk(c, 0.0, 0.0), times)
        np.testing.assert_allclose(
            hist.euler, np.tile([-alpha * c.radius, 0.0, 0.0], (3, 1)), atol=1e-12
        )
        np.testing.assert_allclose(
            hist.felt_weight, (c.omega + alpha * times) ** 2 * c.radius, rtol=1e-12
        )
        np.testing.assert_allclose(hist.apparent[:, 0], -alpha * c.radius, atol=1e-12)

    def test_no_euler_force_at_constant_spin(self, cyl: Cylinder) -> None:
        hist = probe(cyl, velocity_path(cyl, LocalVector(3.0, 2.0, 1.0)))
        np.testing.assert_array_equal(hist.euler, np.zeros((1, 3)))


class TestAnalysisBookkeeping:
    def test_apparent_gravity_is_the_sum_of_its_parts(
        self, spinup_cyl: Cylinder
    ) -> None:
        prof = Trapezoid(500.0, 8.0, 1.5)
        path = StraightPath(spinup_cyl, 30.0, 40.0, prof, prof.duration)
        hist = analyze_path(spinup_cyl, path, np.linspace(0.5, prof.duration - 0.5, 25))
        np.testing.assert_allclose(
            hist.apparent,
            hist.centrifugal + hist.coriolis + hist.euler - hist.acc,
            atol=1e-12,
        )

    def test_support_force_balances_apparent_gravity(
        self, default_cyl: Cylinder
    ) -> None:
        hist = probe(default_cyl, FloorWalk(default_cyl, 0.0, 0.0))
        np.testing.assert_allclose(
            hist.support_force(70.0)[0],
            [0.0, 0.0, 70.0 * G_STANDARD],
            rtol=1e-9,
            atol=1e-9,
        )
        np.testing.assert_array_equal(hist.support_force(3.0), -3.0 * hist.apparent)

    def test_free_fall_setup_is_not_trivial(self, default_cyl: Cylinder) -> None:
        hist = probe(default_cyl, FloorWalk(default_cyl, 0.0, 0.0, height=200.0))
        assert np.linalg.norm(hist.centrifugal[0]) > 1.0

    @pytest.mark.parametrize(
        "v_rel",
        [(0.0, 0.0, 0.0), (10.0, 3.0, 10.0), (-40.0, -8.0, 2.0), (0.0, 25.0, -6.0)],
    )
    @pytest.mark.parametrize("t", [0.0, 2.0, 5.0])
    def test_free_fall_feels_weightless(
        self, flight_cyl: Cylinder, v_rel: tuple[float, float, float], t: float
    ) -> None:
        """A force-free body (straight line in the inertial frame) feels no gravity."""
        p0 = flight_cyl.floor_point(0.3, 5.0, 200.0)
        v0 = flight_cyl.local_basis(p0).T @ np.array(v_rel, float)
        big_v = v0 + np.cross(flight_cyl.omega * np.array([0.0, 0.0, 1.0]), p0)
        path = NumericPath(
            lambda s: flight_cyl.inertial_to_habitat(p0 + big_v * s, s), duration=5.0
        )
        hist = analyze_path(flight_cyl, path, [t])
        np.testing.assert_allclose(hist.apparent[0], 0.0, atol=FREE_FALL_TOL)

    def test_array_shapes_and_scalar_time(self, default_cyl: Cylinder) -> None:
        path = FloorWalk(default_cyl, 90.0, 2.0)
        times = np.linspace(0.0, 6.0, 7)
        hist = analyze_path(default_cyl, path, times)
        for name in ("vel", "acc", "centrifugal", "coriolis", "euler", "apparent"):
            assert getattr(hist, name).shape == (7, 3)
        for name in ("heights", "felt_weight", "lateral_mag", "tilt_deg"):
            assert getattr(hist, name).shape == (7,)
        np.testing.assert_array_equal(hist.times, times)
        assert analyze_path(default_cyl, path, 2.0).apparent.shape == (1, 3)

    def test_derived_properties_of_force_history(self) -> None:
        zeros = np.zeros((3, 3))
        apparent = np.array([[3.0, 4.0, -5.0], [0.0, 0.0, -9.0], [0.0, 2.0, 2.0]])
        hist = ForceHistory(
            np.arange(3.0),
            np.zeros(3),
            zeros,
            zeros,
            zeros,
            zeros,
            zeros,
            apparent,
            g_ref=10.0,
        )
        np.testing.assert_allclose(hist.felt_weight, [5.0, 9.0, -2.0])
        np.testing.assert_allclose(hist.lateral_mag, [5.0, 0.0, 2.0])
        np.testing.assert_allclose(hist.tilt_deg, [45.0, 0.0, 135.0])
        assert hist.pct(5.0) == pytest.approx(50.0)
        np.testing.assert_allclose(hist.support_force(2.0), -2.0 * apparent)

    def test_force_history_validates_shapes_and_reference(self) -> None:
        zeros = np.zeros((3, 3))
        with pytest.raises(ValueError, match="shape"):
            ForceHistory(
                np.arange(3.0),
                np.zeros(3),
                zeros,
                zeros,
                zeros,
                zeros,
                zeros,
                np.zeros((2, 3)),
                g_ref=10.0,
            )
        with pytest.raises(ValueError, match="shape"):
            ForceHistory(
                np.arange(3.0),
                np.zeros(2),
                zeros,
                zeros,
                zeros,
                zeros,
                zeros,
                zeros,
                g_ref=10.0,
            )
        with pytest.raises(ValueError, match="g_ref"):
            ForceHistory(
                np.arange(3.0),
                np.zeros(3),
                zeros,
                zeros,
                zeros,
                zeros,
                zeros,
                zeros,
                g_ref=0.0,
            )

    def test_pulled_down_a_shaft_faster_than_gravity_gives_negative_weight(
        self, default_cyl: Cylinder
    ) -> None:
        class Plunge:
            """s = 10 t^2: a constant 20 m/s^2 acceleration down the shaft."""

            def s(self, t: object) -> np.ndarray:
                return 10.0 * np.asarray(t, float) ** 2

            def ds(self, t: object) -> np.ndarray:
                return 20.0 * np.asarray(t, float)

            def dds(self, t: object) -> np.ndarray:
                return np.full_like(np.asarray(t, float), 20.0)

        path = StraightPath(default_cyl, 0.0, -90.0, Plunge(), 2.0, start_height=100.0)
        hist = probe(default_cyl, path, t=1.0)
        expected = default_cyl.gravity_at_height(90.0) - 20.0
        assert hist.felt_weight[0] == pytest.approx(expected, abs=ACC_TOL)
        assert hist.felt_weight[0] < 0.0
        assert hist.tilt_deg[0] > 90.0

    def test_numeric_and_analytic_paths_agree(self, default_cyl: Cylinder) -> None:
        walk = FloorWalk(default_cyl, 90.0, 6.0)
        analytic = probe(default_cyl, walk)
        numeric = probe(
            default_cyl, NumericPath(walk.position, duration=1.0, step=1e-2)
        )
        np.testing.assert_allclose(analytic.apparent, numeric.apparent, atol=1e-7)

    def test_path_through_the_axis_is_rejected(self, default_cyl: Cylinder) -> None:
        on_axis = NumericPath(lambda t: np.array([0.0, 0.0, t]), duration=1.0)
        with pytest.raises(AxisSingularityError):
            analyze_path(default_cyl, on_axis, [0.5])


class TestElevatorRide:
    @pytest.mark.parametrize(
        ("direction", "t", "height", "v_up", "a_up"),
        [
            ("up", 5.0, 12.5, 5.0, 1.0),
            ("up", 150.0, 1450.0, 10.0, 0.0),
            ("up", 305.0, 2987.5, 5.0, -1.0),
            ("down", 5.0, 2987.5, -5.0, -1.0),
            ("down", 150.0, 1550.0, -10.0, 0.0),
            ("down", 305.0, 12.5, -5.0, 1.0),
        ],
    )
    def test_each_phase_of_the_ride(
        self,
        default_cyl: Cylinder,
        elevator: SimpleNamespace,
        direction: str,
        t: float,
        height: float,
        v_up: float,
        a_up: float,
    ) -> None:
        hist = probe(default_cyl, getattr(elevator, direction), t=t)
        assert hist.heights[0] == pytest.approx(height, abs=1e-9)
        assert hist.vel[0, 2] == pytest.approx(v_up, abs=1e-9)
        assert hist.acc[0, 2] == pytest.approx(a_up, abs=ACC_TOL)
        assert hist.felt_weight[0] == pytest.approx(
            default_cyl.gravity_at_height(height) + a_up, abs=ACC_TOL
        )
        assert hist.apparent[0, 0] == pytest.approx(
            2.0 * default_cyl.omega * v_up, abs=ACC_TOL
        )
        assert hist.apparent[0, 1] == pytest.approx(0.0, abs=ACC_TOL)

    @pytest.mark.parametrize("direction", ["up", "down"])
    def test_sideways_push_tracks_vertical_speed_for_the_whole_ride(
        self, default_cyl: Cylinder, elevator: SimpleNamespace, direction: str
    ) -> None:
        times = np.linspace(0.0, elevator.duration, 311)
        hist = analyze_path(default_cyl, getattr(elevator, direction), times)
        np.testing.assert_allclose(
            hist.apparent[:, 0], 2.0 * default_cyl.omega * hist.vel[:, 2], atol=1e-9
        )
        np.testing.assert_allclose(hist.apparent[:, 1], 0.0, atol=1e-9)

    def test_tilt_grows_as_the_car_approaches_the_axis(
        self, default_cyl: Cylinder, elevator: SimpleNamespace
    ) -> None:
        hist = analyze_path(default_cyl, elevator.up, np.linspace(20.0, 290.0, 28))
        assert np.all(np.diff(hist.tilt_deg) > 0.0)
        assert hist.tilt_deg[0] < 10.0 < hist.tilt_deg[-1]
