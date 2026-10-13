import math

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.core.cylinder import G_STANDARD, AxisSingularityError, Cylinder
from oneill.core.vectors import LocalVector
from oneill.motion.paths import (
    FloorWalk,
    NumericPath,
    Path,
    StraightPath,
    velocity_path,
)
from oneill.motion.profiles import ConstantSpeed, Trapezoid

headings = st.floats(-400.0, 400.0)
speeds = st.floats(-40.0, 40.0)
times = st.floats(0.0, 500.0)


class TestFloorWalk:
    @pytest.mark.parametrize("heading", [0.0, 45.0, 90.0, 270.0])
    @pytest.mark.parametrize("height_frac", [0.0, 0.4])
    def test_follows_the_curved_floor(
        self, cyl: Cylinder, heading: float, height_frac: float
    ) -> None:
        h = height_frac * cyl.radius
        path = FloorWalk(cyl, heading, 5.0, height=h)
        for t in (0.0, 1.0, 60.0, 3600.0):
            assert cyl.height_of(path.position(t)) == pytest.approx(h, abs=1e-8)

    def test_axial_walk_is_a_straight_line_along_the_axis(
        self, default_cyl: Cylinder
    ) -> None:
        path = FloorWalk(default_cyl, 0.0, 1.5, z0=-100.0)
        p0, p1 = path.position(0.0), path.position(40.0)
        np.testing.assert_allclose(p1[:2], p0[:2], atol=1e-9)
        assert p1[2] == pytest.approx(-100.0 + 60.0)
        np.testing.assert_allclose(path.velocity(7.0), [0.0, 0.0, 1.5], atol=1e-12)
        np.testing.assert_allclose(path.acceleration(7.0), 0.0, atol=1e-12)

    @pytest.mark.parametrize("speed", [5.0, -5.0])
    def test_walking_once_around_returns_to_the_start(
        self, cyl: Cylinder, speed: float
    ) -> None:
        path = FloorWalk(cyl, 90.0, speed)
        lap_time = 2.0 * math.pi * cyl.radius / abs(speed)
        np.testing.assert_allclose(
            path.position(lap_time), path.position(0.0), atol=1e-6
        )

    def test_default_name_and_duration(self, default_cyl: Cylinder) -> None:
        path = FloorWalk(default_cyl, 90.0, 1.4)
        assert path.name == "walk spinward 1.4 m/s"
        assert path.duration == 10.0
        assert FloorWalk(default_cyl, 90.0, 1.4, name="x", duration=3.0).name == "x"

    @given(headings, speeds, times, st.floats(0.0, 0.9))
    @example(90.0, 30.0, 0.0, 0.0)
    @example(270.0, -12.0, 100.0, 0.5)
    def test_velocity_and_acceleration_match_numeric_differences(
        self, heading: float, speed: float, t: float, frac: float
    ) -> None:
        cyl = Cylinder()
        walk = FloorWalk(cyl, heading, speed, height=frac * cyl.radius)
        numeric = NumericPath(walk.position, duration=1.0, step=1e-2)
        np.testing.assert_allclose(
            numeric.velocity(t), walk.velocity(t), atol=1e-7, rtol=1e-9
        )
        np.testing.assert_allclose(
            numeric.acceleration(t), walk.acceleration(t), atol=1e-7, rtol=1e-9
        )

    @given(headings, speeds, times)
    def test_acceleration_is_centripetal(
        self, heading: float, speed: float, t: float
    ) -> None:
        cyl = Cylinder()
        walk = FloorWalk(cyl, heading, speed)
        vs = speed * math.sin(math.radians(heading))
        acc = walk.acceleration(t)
        pos = walk.position(t)
        assert acc[2] == 0.0
        assert float(np.hypot(acc[0], acc[1])) == pytest.approx(
            vs * vs / cyl.radius, rel=1e-12, abs=1e-15
        )
        assert float(acc[0] * pos[0] + acc[1] * pos[1]) <= 1e-9

    def test_array_times_give_stacked_vectors(self, default_cyl: Cylinder) -> None:
        walk = FloorWalk(default_cyl, 60.0, 3.0)
        ts = np.array([0.0, 1.0, 2.0])
        for fn in (walk.position, walk.velocity, walk.acceleration):
            out = fn(ts)
            assert out.shape == (3, 3)
            np.testing.assert_array_equal(out[1], fn(1.0))

    @pytest.mark.parametrize("bad_height", [-1.0, 3200.0, 4000.0, math.nan])
    def test_rejects_heights_outside_the_habitat(self, bad_height: float) -> None:
        with pytest.raises(ValueError, match="height"):
            FloorWalk(Cylinder(), 90.0, 1.0, height=bad_height)

    @pytest.mark.parametrize(
        ("kwargs", "match"),
        [
            ({"heading_deg": math.nan}, "heading_deg"),
            ({"speed": math.inf}, "speed"),
            ({"theta0": math.nan}, "theta0"),
            ({"z0": math.nan}, "z0"),
            ({"duration": 0.0}, "duration"),
        ],
    )
    def test_rejects_bad_scalars(self, kwargs: dict[str, float], match: str) -> None:
        args = {"heading_deg": 90.0, "speed": 1.0} | kwargs
        with pytest.raises(ValueError, match=match):
            FloorWalk(Cylinder(), **args)  # type: ignore[arg-type]


class TestStraightPath:
    @pytest.mark.parametrize(
        ("elevation", "start_height"), [(90.0, 0.0), (-90.0, 500.0)]
    )
    def test_vertical_path_changes_height_by_distance_travelled(
        self, cyl: Cylinder, elevation: float, start_height: float
    ) -> None:
        start_height = min(start_height, 0.5 * cyl.radius)
        path = StraightPath(
            cyl,
            0.0,
            elevation,
            ConstantSpeed(2.0),
            10.0,
            start_height=start_height,
        )
        sign = 1.0 if elevation > 0 else -1.0
        assert cyl.height_of(path.position(10.0)) == pytest.approx(
            start_height + sign * 20.0, abs=1e-9
        )

    def test_horizontal_straight_line_leaves_the_curved_floor(
        self, cyl: Cylinder
    ) -> None:
        d = 100.0
        path = StraightPath(cyl, 90.0, 0.0, ConstantSpeed(1.0), d)
        height = cyl.height_of(path.position(d))
        assert height < 0.0
        assert height == pytest.approx(cyl.radius - math.hypot(cyl.radius, d), rel=1e-9)

    def test_cannot_start_on_the_axis(self, default_cyl: Cylinder) -> None:
        with pytest.raises(AxisSingularityError):
            StraightPath(
                default_cyl,
                0.0,
                90.0,
                ConstantSpeed(1.0),
                1.0,
                start_height=default_cyl.radius,
            )

    @pytest.mark.parametrize("bad", [-1.0, 3201.0, math.nan])
    def test_rejects_start_heights_outside_the_habitat(self, bad: float) -> None:
        with pytest.raises(ValueError, match="start_height"):
            StraightPath(
                Cylinder(), 0.0, 90.0, ConstantSpeed(1.0), 1.0, start_height=bad
            )

    def test_derivatives_follow_the_profile(self, default_cyl: Cylinder) -> None:
        prof = Trapezoid(3000.0, 10.0, 1.0)
        path = StraightPath(
            default_cyl, 0.0, 90.0, prof, prof.duration, start_height=0.0
        )
        up = default_cyl.local_basis(path.position(0.0)).T @ np.array([0.0, 0.0, 1.0])
        np.testing.assert_allclose(path.velocity(5.0), 5.0 * up, atol=1e-12)
        np.testing.assert_allclose(path.acceleration(5.0), 1.0 * up, atol=1e-12)
        np.testing.assert_allclose(path.acceleration(150.0), 0.0, atol=1e-12)
        out = path.position(np.array([0.0, 5.0]))
        assert out.shape == (2, 3)

    @given(
        headings,
        st.floats(-90.0, 90.0),
        st.floats(0.0, 0.9),
        st.floats(0.0, 200.0),
    )
    @example(30.0, 40.0, 0.0, 12.0)
    def test_velocity_and_acceleration_match_numeric_differences(
        self, heading: float, elevation: float, frac: float, t: float
    ) -> None:
        cyl = Cylinder()
        prof = Trapezoid(100000.0, 10.0, 1.0)
        path = StraightPath(
            cyl, heading, elevation, prof, 10.0, start_height=frac * cyl.radius
        )
        # Stay well inside the cruise phase so no profile kink is in the stencil.
        t += 20.0
        numeric = NumericPath(path.position, duration=1.0, step=1e-2)
        np.testing.assert_allclose(numeric.velocity(t), path.velocity(t), atol=1e-6)
        np.testing.assert_allclose(
            numeric.acceleration(t), path.acceleration(t), atol=1e-6
        )

    def test_default_name(self, default_cyl: Cylinder) -> None:
        path = StraightPath(default_cyl, 0.0, 90.0, ConstantSpeed(1.0), 1.0)
        assert path.name == "straight +axial elev 90 deg"
        named = StraightPath(
            default_cyl, 0.0, 90.0, ConstantSpeed(1.0), 1.0, name="ladder"
        )
        assert named.name == "ladder"


class TestNumericPath:
    def test_exact_for_quadratic_motion(self) -> None:
        p0 = np.array([3000.0, 10.0, 5.0])
        v = np.array([1.0, -2.0, 3.0])
        a = np.array([0.1, 0.2, -0.3])

        def pos(t: float) -> np.ndarray:
            return p0 + v * t + 0.5 * a * t * t

        path = NumericPath(pos, duration=5.0)
        for t in (0.0, 2.0, 100.0):
            np.testing.assert_allclose(path.position(t), pos(t))
            np.testing.assert_allclose(path.velocity(t), v + a * t, atol=1e-6)
            np.testing.assert_allclose(path.acceleration(t), a, atol=1e-6)

    def test_array_times_and_metadata(self) -> None:
        path = NumericPath(lambda t: np.array([t, 2.0 * t, 0.0]), 4.0, name="line")
        ts = np.array([0.0, 1.0, 2.0])
        assert path.velocity(ts).shape == (3, 3)
        np.testing.assert_allclose(path.velocity(ts), [[1.0, 2.0, 0.0]] * 3, atol=1e-9)
        assert (path.name, path.duration) == ("line", 4.0)

    def test_default_name(self) -> None:
        assert (
            NumericPath(lambda t: np.array([t, 0.0, 0.0]), 1.0).name == "numeric path"
        )

    @pytest.mark.parametrize("step", [0.0, -1e-2, math.nan])
    def test_rejects_bad_step(self, step: float) -> None:
        with pytest.raises(ValueError, match="step"):
            NumericPath(lambda t: np.array([t, 0.0, 0.0]), 1.0, step=step)

    def test_step_size_does_not_matter_when_converged(
        self, default_cyl: Cylinder
    ) -> None:
        walk = FloorWalk(default_cyl, 90.0, 6.0)
        coarse = NumericPath(walk.position, 1.0, step=1e-2)
        fine = NumericPath(walk.position, 1.0, step=1e-3)
        np.testing.assert_allclose(
            coarse.acceleration(0.0), fine.acceleration(0.0), atol=1e-5
        )

    def test_conforms_to_the_path_protocol(self, default_cyl: Cylinder) -> None:
        paths: list[Path] = [
            FloorWalk(default_cyl, 90.0, 1.0),
            StraightPath(default_cyl, 0.0, 90.0, ConstantSpeed(1.0), 1.0),
            NumericPath(lambda t: np.array([t, 0.0, 0.0]), 1.0),
        ]
        assert all(p.duration > 0 for p in paths)


class TestVelocityPath:
    @pytest.mark.parametrize(
        "velocity",
        [
            (1.4, 0.0, 0.0),
            (0.0, -2.0, 0.0),
            (3.0, 4.0, 0.0),
            (0.0, 0.0, 1.0),
            (0.0, 0.0, -3.0),
            (1.0, 2.0, 0.5),
            (-0.57, 0.0, 0.4),
            (5.0, -5.0, -5.0),
        ],
    )
    def test_reproduces_the_requested_velocity(
        self, cyl: Cylinder, velocity: tuple[float, float, float]
    ) -> None:
        path = velocity_path(cyl, LocalVector(*velocity), height=0.1 * cyl.radius)
        pos = path.position(0.0)
        local = cyl.local_basis(pos) @ path.velocity(0.0)
        np.testing.assert_allclose(local, velocity, atol=1e-9)

    def test_level_motion_follows_floor_but_climbing_is_straight(
        self, default_cyl: Cylinder
    ) -> None:
        level = velocity_path(default_cyl, LocalVector(30.0, 0.0, 0.0))
        climb = velocity_path(default_cyl, LocalVector(30.0, 0.0, 1.0))
        assert isinstance(level, FloorWalk)
        assert isinstance(climb, StraightPath)
        assert default_cyl.height_of(level.position(60.0)) == pytest.approx(
            0.0, abs=1e-8
        )
        assert default_cyl.height_of(climb.position(60.0)) != pytest.approx(
            60.0, abs=1e-3
        )

    def test_zero_velocity_is_standing_still(self, default_cyl: Cylinder) -> None:
        path = velocity_path(default_cyl, LocalVector(0.0, 0.0, 0.0))
        np.testing.assert_array_equal(path.position(0.0), path.position(100.0))
        assert G_STANDARD > 0

    def test_name_is_forwarded(self, default_cyl: Cylinder) -> None:
        assert (
            velocity_path(default_cyl, LocalVector(1.0, 0.0, 0.0), name="a").name == "a"
        )
        climb = velocity_path(default_cyl, LocalVector(1.0, 0.0, 1.0), name="b")
        assert climb.name == "b"
