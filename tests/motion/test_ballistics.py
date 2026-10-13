import itertools
import math

import numpy as np
import pytest
from hypothesis import example, given, settings
from hypothesis import strategies as st

from oneill.core.cylinder import G_STANDARD, Cylinder
from oneill.core.vectors import LocalVector, Vec3
from oneill.motion.ballistics import (
    BallisticResult,
    ballistic,
    integrate_rotating,
)

POS_TOL = 1e-5  # m, agreement between independent trajectory solvers
REST = LocalVector(0.0, 0.0, 0.0)


def wrap_angle(angle: float) -> float:
    """Wrap an angle in radians into ``[-pi, pi)``."""
    return (angle + math.pi) % (2.0 * math.pi) - math.pi


def drop_reference(cyl: Cylinder, height: float) -> tuple[float, float]:
    """Exact (flight time, spinward landing offset) for release at rest.

    Inertial frame: the object leaves radius ``r0 = R - h`` tangentially at
    ``omega r0`` and travels the chord ``sqrt(R^2 - r0^2)`` in a straight line
    while the floor turns through ``omega T``.
    """
    r0 = cyl.radius - height
    chord = math.sqrt(cyl.radius**2 - r0**2)
    t = chord / (cyl.omega * r0)
    return t, cyl.radius * (math.atan2(chord, r0) - cyl.omega * t)


def jump_reference(cyl: Cylinder, v_up: float) -> tuple[float, float]:
    """Exact (flight time, spinward landing offset) for a vertical jump.

    A chord leaving a circle at angle ``beta = atan(v_up / (omega R))`` to the
    tangent subtends ``2 beta`` at the centre.
    """
    x = v_up / cyl.rim_speed
    t = 2.0 * cyl.radius * v_up / (v_up**2 + cyl.rim_speed**2)
    return t, cyl.radius * (2.0 * math.atan(x) - cyl.omega * t)


class TestDroppedObjects:
    @pytest.mark.parametrize(
        "height_frac", [1e-6, 5e-4, 0.03, 0.3, 0.9, 0.9997, 0.999997]
    )
    def test_matches_inertial_frame_closed_form(
        self, cyl: Cylinder, height_frac: float
    ) -> None:
        """From millimetres above the floor to millimetres from the axis.

        The last cases fly for hundreds to hundreds of thousands of revolutions.
        """
        h = height_frac * cyl.radius
        t_ref, landing_ref = drop_reference(cyl, h)
        res = ballistic(cyl, h, REST)
        assert res.flight_time == pytest.approx(t_ref, rel=1e-9)
        assert res.landing_spin == pytest.approx(landing_ref, rel=1e-9, abs=1e-9)
        assert res.landing_axial == pytest.approx(0.0, abs=1e-12)

    @settings(max_examples=60)
    @given(st.floats(10.0, 1.0e5), st.floats(0.05, 3.0), st.floats(1e-7, 0.999999))
    @example(3200.0, 0.5, 0.03)
    @example(3200.0, 0.5, 0.999997)
    def test_closed_form_property(self, radius: float, rpm: float, frac: float) -> None:
        cyl = Cylinder.from_rpm(radius, 100.0, rpm)
        t_ref, landing_ref = drop_reference(cyl, frac * radius)
        res = ballistic(cyl, frac * radius, REST)
        assert res.flight_time == pytest.approx(t_ref, rel=1e-8)
        assert res.landing_spin == pytest.approx(landing_ref, rel=1e-8, abs=1e-8)

    def test_always_lands_antispinward_and_further_from_higher_up(
        self, default_cyl: Cylinder
    ) -> None:
        heights = [0.5, 1.5, 10.0, 100.0, 1000.0, 3000.0]
        landings = [ballistic(default_cyl, h, REST).landing_spin for h in heights]
        assert all(x < 0.0 for x in landings)
        assert all(b < a for a, b in itertools.pairwise(landings))

    @pytest.mark.parametrize("height", [1.5, 100.0, 1000.0])
    def test_falls_more_slowly_than_uniform_gravity_predicts(
        self, default_cyl: Cylinder, height: float
    ) -> None:
        res = ballistic(default_cyl, height, REST)
        assert res.naive_time == pytest.approx(
            math.sqrt(2.0 * height / G_STANDARD), rel=1e-12
        )
        assert res.naive_spin == 0.0
        assert res.naive_axial == 0.0
        assert res.flight_time > res.naive_time

    def test_small_drop_agrees_with_first_order_coriolis_theory(
        self, default_cyl: Cylinder
    ) -> None:
        res = ballistic(default_cyl, 1.5, REST)
        first_order = -default_cyl.omega * G_STANDARD * res.naive_time**3 / 3.0
        assert res.landing_spin == pytest.approx(first_order, rel=5e-3)

    def test_trajectory_is_well_formed(self, default_cyl: Cylinder) -> None:
        res = ballistic(default_cyl, 100.0, REST, n=250)
        assert (
            len(res.traj_t)
            == len(res.traj_spin)
            == len(res.traj_axial)
            == len(res.traj_height)
            == 250
        )
        assert res.traj_t[0] == 0.0
        assert res.traj_t[-1] == pytest.approx(res.flight_time)
        assert (res.traj_spin[0], res.traj_axial[0]) == (0.0, 0.0)
        assert res.traj_height[0] == pytest.approx(100.0, abs=1e-9)
        assert res.traj_height[-1] == pytest.approx(0.0, abs=1e-6)
        assert np.all(np.diff(res.traj_height) < 0.0)
        assert res.traj_spin[-1] == res.landing_spin

    @pytest.mark.parametrize("alpha", [1.0e-4, -1.0e-4])
    def test_spin_up_shifts_the_landing_point_by_the_extra_floor_travel(
        self, alpha: float
    ) -> None:
        """The falling object is force-free; only the floor moves differently."""
        base = ballistic(Cylinder.from_gravity(), 100.0, REST)
        res = ballistic(Cylinder.from_gravity(spin_accel=alpha), 100.0, REST)
        assert res.flight_time == pytest.approx(base.flight_time, rel=1e-12)
        shift = 0.5 * alpha * base.flight_time**2 * 3200.0
        assert res.landing_spin == pytest.approx(base.landing_spin - shift, abs=1e-6)


class TestJumpsAndThrows:
    @pytest.mark.parametrize("v_up", [0.5, 2.5, 10.0, 100.0, 1000.0])
    def test_vertical_jump_matches_inertial_frame_closed_form(
        self, cyl: Cylinder, v_up: float
    ) -> None:
        t_ref, landing_ref = jump_reference(cyl, v_up)
        res = ballistic(cyl, 0.0, LocalVector(0.0, 0.0, v_up))
        assert res.flight_time == pytest.approx(t_ref, rel=1e-9)
        assert res.landing_spin == pytest.approx(landing_ref, rel=1e-9, abs=1e-8)
        assert res.naive_time == pytest.approx(
            2.0 * v_up / cyl.surface_gravity, rel=1e-12
        )

    @settings(max_examples=60)
    @given(st.floats(10.0, 1.0e5), st.floats(0.05, 3.0), st.floats(0.01, 500.0))
    def test_jump_closed_form_property(
        self, radius: float, rpm: float, v_up: float
    ) -> None:
        cyl = Cylinder.from_rpm(radius, 100.0, rpm)
        t_ref, landing_ref = jump_reference(cyl, v_up)
        res = ballistic(cyl, 0.0, LocalVector(0.0, 0.0, v_up))
        assert res.flight_time == pytest.approx(t_ref, rel=1e-8)
        assert res.landing_spin == pytest.approx(landing_ref, rel=1e-8, abs=1e-7)

    def test_jump_lands_slightly_spinward_as_first_order_theory_predicts(
        self, default_cyl: Cylinder
    ) -> None:
        v = 2.5
        res = ballistic(default_cyl, 0.0, LocalVector(0.0, 0.0, v))
        assert res.landing_spin > 0.0
        assert res.landing_spin == pytest.approx(
            4.0 / 3.0 * default_cyl.omega * v**3 / G_STANDARD**2, rel=1e-2
        )

    def test_jump_apex_height(self, default_cyl: Cylinder) -> None:
        v = 2.5
        res = ballistic(default_cyl, 0.0, LocalVector(0.0, 0.0, v), n=2001)
        x = v / default_cyl.rim_speed
        exact = default_cyl.radius * (1.0 - 1.0 / math.sqrt(1.0 + x * x))
        assert res.traj_height.max() == pytest.approx(exact, abs=1e-6)
        assert res.traj_height.max() == pytest.approx(
            v * v / (2.0 * G_STANDARD), rel=1e-3
        )

    @pytest.mark.parametrize(
        ("start_height", "v_local"),
        [
            pytest.param(1.5, (0.0, 0.0, 0.0), id="drop-from-hand"),
            pytest.param(100.0, (0.0, 0.0, 0.0), id="drop-100m"),
            pytest.param(1000.0, (0.0, 0.0, 0.0), id="drop-1000m"),
            pytest.param(0.0, (0.0, 0.0, 2.5), id="vertical-jump"),
            pytest.param(1.5, (10.6, 0.0, 10.6), id="throw-spinward"),
            pytest.param(1.5, (-10.6, 0.0, 10.6), id="throw-antispinward"),
            pytest.param(1.5, (0.0, 10.6, 10.6), id="throw-axial"),
            pytest.param(1.5, (7.0, -5.0, 12.0), id="throw-diagonal"),
            pytest.param(50.0, (-30.0, 2.0, -5.0), id="thrown-down-from-balcony"),
            pytest.param(2.0, (60.0, 0.0, 0.0), id="fast-level-spinward"),
        ],
    )
    def test_inertial_solution_agrees_with_rotating_frame_integration(
        self,
        flight_cyl: Cylinder,
        start_height: float,
        v_local: tuple[float, float, float],
    ) -> None:
        """Exact inertial straight line vs RK4 with all three fictitious forces."""
        res = ballistic(flight_cyl, start_height, LocalVector(*v_local))
        p0 = flight_cyl.floor_point(0.0, 0.0, start_height)
        v0 = flight_cyl.local_basis(p0).T @ np.array(v_local, float)
        p_end, _ = integrate_rotating(
            flight_cyl, p0, v0, res.flight_time, max_step=0.01
        )
        assert math.hypot(p_end[0], p_end[1]) == pytest.approx(
            flight_cyl.radius, abs=POS_TOL
        )
        assert flight_cyl.radius * math.atan2(p_end[1], p_end[0]) == pytest.approx(
            res.landing_spin, abs=POS_TOL
        )
        assert p_end[2] == pytest.approx(res.landing_axial, abs=POS_TOL)

    @settings(max_examples=5)
    @given(
        st.floats(-30.0, 30.0),
        st.floats(-30.0, 30.0),
        st.floats(1.0, 30.0),
        st.floats(0.0, 50.0),
    )
    def test_solvers_agree_for_random_throws(
        self, vs: float, va: float, vu: float, h0: float
    ) -> None:
        cyl = Cylinder.from_gravity(spin_accel=5e-5)
        res = ballistic(cyl, h0, LocalVector(vs, va, vu))
        p0 = cyl.floor_point(0.0, 0.0, h0)
        v0 = cyl.local_basis(p0).T @ np.array([vs, va, vu])
        p_end, _ = integrate_rotating(cyl, p0, v0, res.flight_time, max_step=0.01)
        assert cyl.radius * math.atan2(p_end[1], p_end[0]) == pytest.approx(
            res.landing_spin, abs=POS_TOL
        )
        assert p_end[2] == pytest.approx(res.landing_axial, abs=POS_TOL)

    @pytest.mark.parametrize("v_axial", [-20.0, 3.0, 15.0])
    def test_axial_motion_is_unaffected_by_rotation(
        self, default_cyl: Cylinder, v_axial: float
    ) -> None:
        res = ballistic(default_cyl, 1.5, LocalVector(4.0, v_axial, 6.0))
        assert res.landing_axial == pytest.approx(v_axial * res.flight_time, rel=1e-12)
        np.testing.assert_allclose(
            res.traj_axial, v_axial * res.traj_t, rtol=1e-12, atol=1e-12
        )

    @given(
        st.floats(-20.0, 20.0),
        st.floats(-1e4, 1e4),
        st.sampled_from([(100.0, (0.0, 0.0, 0.0)), (1.5, (7.0, -5.0, 12.0))]),
    )
    @example(math.pi, 1234.0, (100.0, (0.0, 0.0, 0.0)))
    @example(-2.0, 1234.0, (1.5, (7.0, -5.0, 12.0)))
    @example(1.5 * math.pi, 1234.0, (1.5, (7.0, -5.0, 12.0)))
    def test_result_does_not_depend_on_launch_position(
        self,
        theta0: float,
        z0: float,
        case: tuple[float, tuple[float, float, float]],
    ) -> None:
        start_height, v_local = case
        cyl = Cylinder()
        vec = LocalVector(*v_local)
        reference = ballistic(cyl, start_height, vec)
        res = ballistic(cyl, start_height, vec, theta0=theta0, z0=z0)
        assert res.flight_time == pytest.approx(reference.flight_time, rel=1e-9)
        assert res.landing_spin == pytest.approx(reference.landing_spin, abs=1e-6)
        assert res.landing_axial == pytest.approx(reference.landing_axial, abs=1e-6)

    @pytest.mark.parametrize(
        "v_local", [(3.0, 1.0, -2.0), (5.0, 0.0, 0.0), (0.0, 4.0, 0.0)]
    )
    def test_launch_from_the_floor_without_upward_speed_lands_immediately(
        self, default_cyl: Cylinder, v_local: tuple[float, float, float]
    ) -> None:
        res = ballistic(default_cyl, 0.0, LocalVector(*v_local))
        assert res.flight_time == pytest.approx(0.0, abs=1e-9)
        assert res.landing_spin == pytest.approx(0.0, abs=1e-9)
        assert res.landing_axial == pytest.approx(0.0, abs=1e-9)
        assert np.all(np.isfinite(res.traj_height))

    def test_object_at_rest_in_the_inertial_frame_never_lands(
        self, default_cyl: Cylinder
    ) -> None:
        r = default_cyl.radius - 1.5
        with pytest.raises(ValueError, match="rest in the inertial frame"):
            ballistic(default_cyl, 1.5, LocalVector(-default_cyl.omega * r, 0.0, 0.0))

    def test_almost_hovering_object_drifts_for_many_revolutions(
        self, default_cyl: Cylinder
    ) -> None:
        """1 mm/s short of hovering: ~27 hours, ~860 habitat revolutions."""
        h, slip = 1.5, 1.0e-3
        r = default_cyl.radius - h
        chord = math.sqrt(default_cyl.radius**2 - r**2)
        t_ref = chord / slip
        landing_ref = default_cyl.radius * (
            math.atan2(chord, r) - default_cyl.omega * t_ref
        )
        res = ballistic(
            default_cyl, h, LocalVector(-default_cyl.omega * r + slip, 0.0, 0.0)
        )
        assert res.flight_time == pytest.approx(t_ref, rel=1e-8)
        assert res.landing_spin == pytest.approx(landing_ref, rel=1e-8)
        assert res.traj_height.min() >= -1e-6
        assert res.traj_height.max() <= h + 1e-9

    def test_trajectory_straight_through_the_axis(self, default_cyl: Cylinder) -> None:
        """Cancel the spin and jump: a diameter through the axis, no NaNs."""
        v_up = 50.0
        res = ballistic(
            default_cyl, 0.0, LocalVector(-default_cyl.rim_speed, 0.0, v_up), n=401
        )
        t_ref = 2.0 * default_cyl.radius / v_up
        assert res.flight_time == pytest.approx(t_ref, rel=1e-9)
        landing_angle = res.landing_spin / default_cyl.radius
        assert wrap_angle(
            landing_angle - (math.pi - default_cyl.omega * t_ref)
        ) == pytest.approx(0.0, abs=1e-9)
        assert np.all(np.isfinite(res.traj_spin))
        assert res.traj_height.max() == pytest.approx(default_cyl.radius, abs=1e-6)

    @pytest.mark.parametrize(
        ("z0_offset", "v_axial", "expected"),
        [
            (None, 10.0, False),
            (-1.0, 10.0, True),
            (-1.0, -10.0, False),
            (+1.0, -10.0, True),
        ],
    )
    def test_end_cap_detection(
        self,
        default_cyl: Cylinder,
        z0_offset: float | None,
        v_axial: float,
        expected: bool,  # noqa: FBT001
    ) -> None:
        half = default_cyl.length / 2.0
        if z0_offset is None:
            z0 = 0.0
        elif z0_offset < 0:
            z0 = half + z0_offset
        else:
            z0 = -half + z0_offset
        res = ballistic(default_cyl, 1.5, LocalVector(0.0, v_axial, 5.0), z0=z0)
        assert res.hits_endcap is expected

    @pytest.mark.parametrize("bad_height", [-1.0, -1e-9, 3200.0, 5000.0, math.nan])
    def test_launch_point_must_be_inside_the_hull_and_off_the_axis(
        self, default_cyl: Cylinder, bad_height: float
    ) -> None:
        with pytest.raises(ValueError, match="start_height"):
            ballistic(default_cyl, bad_height, LocalVector(0.0, 0.0, 1.0))

    def test_rejects_too_few_trajectory_samples(self, default_cyl: Cylinder) -> None:
        with pytest.raises(ValueError, match="at least 2"):
            ballistic(default_cyl, 1.5, REST, n=1)

    def test_naive_flat_ground_prediction_and_metadata(
        self, default_cyl: Cylinder
    ) -> None:
        res = ballistic(
            default_cyl, 1.5, LocalVector(3.0, 4.0, 2.0), label="test throw", n=2
        )
        g = G_STANDARD
        t_naive = (2.0 + math.sqrt(4.0 + 2.0 * g * 1.5)) / g
        assert res.naive_time == pytest.approx(t_naive, rel=1e-12)
        assert (res.naive_spin, res.naive_axial) == pytest.approx(
            (3.0 * t_naive, 4.0 * t_naive), rel=1e-12
        )
        assert res.label == "test throw"
        assert res.start_height == 1.5
        assert res.v_launch == LocalVector(3.0, 4.0, 2.0)
        assert len(res.traj_t) == 2
        assert isinstance(res, BallisticResult)


class TestRotatingFrameIntegrator:
    def test_zero_duration_returns_untouched_copies(
        self, default_cyl: Cylinder
    ) -> None:
        p0, v0 = np.array([3000, 10, 5]), np.array([1, 2, 3])
        p, v = integrate_rotating(default_cyl, p0, v0, 0.0)
        np.testing.assert_array_equal(p, [3000.0, 10.0, 5.0])
        np.testing.assert_array_equal(v, [1.0, 2.0, 3.0])
        assert p is not p0
        assert v is not v0
        assert p0.dtype.kind == "i"

    def test_body_at_rest_in_inertial_frame_circles_antispinward(
        self, flight_cyl: Cylinder
    ) -> None:
        """Centrifugal + Coriolis combine into a net centripetal pull of omega^2 r."""
        p0 = flight_cyl.floor_point(0.4, 7.0, 200.0)
        v0 = -np.cross(flight_cyl.omega * np.array([0.0, 0.0, 1.0]), p0)
        t = flight_cyl.period / 4.0
        p, _ = integrate_rotating(flight_cyl, p0, v0, t, max_step=0.05)
        np.testing.assert_allclose(
            p, flight_cyl.inertial_to_habitat(p0, t), atol=POS_TOL
        )

    @pytest.mark.parametrize(
        "v0", [(0.0, 0.0, 0.0), (10.0, 3.0, 10.0), (-100.0, 0.0, 50.0)]
    )
    def test_jacobi_integral_is_conserved(
        self, default_cyl: Cylinder, v0: tuple[float, float, float]
    ) -> None:
        """At constant spin E = v^2/2 - omega^2 rho^2/2 is conserved."""

        def jacobi(p: Vec3, v: Vec3) -> float:
            kinetic = 0.5 * float(np.dot(v, v))
            return kinetic - 0.5 * default_cyl.omega**2 * float(p[0] ** 2 + p[1] ** 2)

        p0, vel0 = default_cyl.floor_point(0.0, 0.0, 500.0), np.array(v0)
        p1, v1 = integrate_rotating(default_cyl, p0, vel0, 10.0, max_step=0.02)
        assert jacobi(p1, v1) == pytest.approx(jacobi(p0, vel0), rel=1e-9)

    def test_extra_acceleration_adds_to_the_fictitious_forces(
        self, default_cyl: Cylinder
    ) -> None:
        """A constant axial push is not touched by rotation: dz = v t + a t^2 / 2."""
        push = 0.3

        def extra(_t: float, _p: Vec3, _v: Vec3) -> Vec3:
            return np.array([0.0, 0.0, push])

        p0 = default_cyl.floor_point(0.0, 0.0, 100.0)
        v0 = np.array([0.0, 0.0, 2.0])
        p, v = integrate_rotating(
            default_cyl, p0, v0, 10.0, max_step=0.05, extra_acc=extra
        )
        assert p[2] == pytest.approx(2.0 * 10.0 + 0.5 * push * 100.0, abs=1e-9)
        assert v[2] == pytest.approx(2.0 + push * 10.0, abs=1e-12)

    def test_extra_acceleration_receives_time_position_and_velocity(
        self, default_cyl: Cylinder
    ) -> None:
        """Linear drag on the axial velocity decays as exp(-k t)."""
        k = 0.2
        seen: list[float] = []

        def drag(t: float, _p: Vec3, v: Vec3) -> Vec3:
            seen.append(t)
            return np.array([0.0, 0.0, -k * v[2]])

        p0 = default_cyl.floor_point(0.0, 0.0, 100.0)
        _, v = integrate_rotating(
            default_cyl,
            p0,
            np.array([0.0, 0.0, 5.0]),
            4.0,
            max_step=0.02,
            extra_acc=drag,
        )
        assert v[2] == pytest.approx(5.0 * math.exp(-k * 4.0), rel=1e-9)
        assert min(seen) == 0.0
        assert max(seen) == pytest.approx(4.0)
