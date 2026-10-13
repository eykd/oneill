import math

import numpy as np
import pytest
from hypothesis import example, given
from hypothesis import strategies as st

from oneill.core.cylinder import (
    G_STANDARD,
    AxisSingularityError,
    Cylinder,
    HabitatSummary,
)
from oneill.core.vectors import Z_HAT

radii = st.floats(1.0, 1.0e5)
angles = st.floats(-20.0, 20.0)
fractions = st.floats(0.0, 0.999)


class TestParameters:
    def test_default_matches_design_brief(self, default_cyl: Cylinder) -> None:
        assert default_cyl.radius == 3200.0
        assert default_cyl.length == 20000.0
        assert default_cyl.surface_gravity == pytest.approx(G_STANDARD, rel=1e-12)
        assert default_cyl.rpm == pytest.approx(0.5286, abs=5e-4)
        assert default_cyl.period == pytest.approx(113.5, abs=0.05)
        assert default_cyl.rim_speed == pytest.approx(177.15, abs=0.01)
        assert default_cyl.spin_accel == 0.0

    @pytest.mark.parametrize(
        ("radius", "g"),
        [(3200.0, G_STANDARD), (100.0, 3.71), (8000.0, 1.62), (0.5, 9.81), (1e6, 9.81)],
    )
    def test_from_gravity_round_trip(self, radius: float, g: float) -> None:
        c = Cylinder.from_gravity(radius, 1000.0, g)
        assert c.surface_gravity == pytest.approx(g, rel=1e-12)
        assert c.gravity_at_height(0.0) == pytest.approx(g, rel=1e-12)

    @given(radii, st.floats(0.1, 50.0))
    def test_from_gravity_round_trip_property(self, radius: float, g: float) -> None:
        c = Cylinder.from_gravity(radius, 10.0, g)
        assert c.surface_gravity == pytest.approx(g, rel=1e-12)

    @pytest.mark.parametrize("rpm", [0.01, 0.5, 1.0, 2.0, 30.0])
    def test_from_rpm_round_trip(self, rpm: float) -> None:
        c = Cylinder.from_rpm(3200.0, 20000.0, rpm)
        assert c.rpm == pytest.approx(rpm, rel=1e-12)
        assert c.period == pytest.approx(60.0 / rpm, rel=1e-12)

    @given(radii, st.floats(0.01, 60.0))
    def test_from_rpm_round_trip_property(self, radius: float, rpm: float) -> None:
        c = Cylinder.from_rpm(radius, 10.0, rpm)
        assert c.rpm == pytest.approx(rpm, rel=1e-12)

    def test_from_helpers_forward_spin_accel(self) -> None:
        assert Cylinder.from_gravity(spin_accel=2e-4).spin_accel == 2e-4
        assert Cylinder.from_rpm(spin_accel=-2e-4).spin_accel == -2e-4

    def test_derived_quantities_are_mutually_consistent(self, cyl: Cylinder) -> None:
        assert cyl.rim_speed**2 / cyl.radius == pytest.approx(
            cyl.surface_gravity, rel=1e-12
        )
        assert cyl.rim_speed * cyl.period == pytest.approx(
            2.0 * math.pi * cyl.radius, rel=1e-12
        )
        assert cyl.rpm * cyl.period == pytest.approx(60.0, rel=1e-12)

    @pytest.mark.parametrize("radius", [10.0, 500.0, 3200.0, 50000.0])
    def test_default_spin_gives_one_g_for_any_radius(self, radius: float) -> None:
        assert Cylinder(radius=radius).surface_gravity == pytest.approx(
            G_STANDARD, rel=1e-12
        )

    def test_explicit_omega_is_respected(self) -> None:
        c = Cylinder(radius=100.0, length=50.0, omega=0.3)
        assert c.omega == 0.3
        assert c.surface_gravity == pytest.approx(9.0, rel=1e-12)

    @pytest.mark.parametrize(
        ("radius", "length"),
        [(0.0, 100.0), (-1.0, 100.0), (100.0, 0.0), (100.0, -5.0)],
    )
    def test_rejects_non_physical_dimensions(
        self, radius: float, length: float
    ) -> None:
        with pytest.raises(ValueError, match="positive"):
            Cylinder(radius=radius, length=length)

    @pytest.mark.parametrize("bad", [math.inf, -math.inf, 0.0, -0.1])
    def test_rejects_non_physical_spin(self, bad: float) -> None:
        with pytest.raises(ValueError, match="omega"):
            Cylinder(omega=bad)

    def test_rejects_non_finite_spin_accel(self) -> None:
        with pytest.raises(ValueError, match="spin_accel"):
            Cylinder(spin_accel=math.nan)

    @pytest.mark.parametrize("bad", [0.0, -1.0, math.nan])
    def test_from_gravity_rejects_non_positive_gravity(self, bad: float) -> None:
        with pytest.raises(ValueError, match="surface_g"):
            Cylinder.from_gravity(surface_g=bad)

    def test_is_immutable_and_hashable(self, default_cyl: Cylinder) -> None:
        with pytest.raises(AttributeError):
            default_cyl.radius = 1.0  # type: ignore[misc]
        assert default_cyl == Cylinder()
        assert hash(default_cyl) == hash(Cylinder())

    @pytest.mark.parametrize("frac", [0.0, 0.25, 0.5, 0.9, 1.0])
    def test_gravity_falls_linearly_to_zero_on_the_axis(
        self, cyl: Cylinder, frac: float
    ) -> None:
        expected = cyl.surface_gravity * (1.0 - frac)
        assert cyl.gravity_at_height(frac * cyl.radius) == pytest.approx(
            expected, rel=1e-12, abs=1e-12
        )

    def test_constant_spin_kinematics(self, cyl: Cylinder) -> None:
        assert cyl.omega_at(123.4) == cyl.omega
        assert cyl.frame_angle(0.0) == 0.0
        assert cyl.frame_angle(cyl.period) == pytest.approx(2.0 * math.pi, rel=1e-12)

    @pytest.mark.parametrize("alpha", [1.0e-4, -1.0e-4])
    def test_spin_up_kinematics(self, alpha: float) -> None:
        c = Cylinder.from_gravity(spin_accel=alpha)
        t, h = 37.0, 1.0e-3
        assert c.omega_at(t) == pytest.approx(c.omega + alpha * t, rel=1e-14)
        d_angle = (c.frame_angle(t + h) - c.frame_angle(t - h)) / (2.0 * h)
        assert d_angle == pytest.approx(c.omega_at(t), rel=1e-9)


class TestSummary:
    def test_default_summary_values(self, default_cyl: Cylinder) -> None:
        s = default_cyl.summary()
        assert isinstance(s, HabitatSummary)
        assert s.radius == 3200.0
        assert s.length == 20000.0
        assert s.omega == default_cyl.omega
        assert s.rpm == default_cyl.rpm
        assert s.period == default_cyl.period
        assert s.rim_speed == default_cyl.rim_speed
        assert s.surface_gravity == default_cyl.surface_gravity
        assert s.gravity_in_g == pytest.approx(1.0, rel=1e-12)
        assert s.gradient_pct_per_m == pytest.approx(100.0 / 3200.0)
        assert s.coriolis_scale == pytest.approx(2.0 * default_cyl.omega)
        assert s.coriolis_pct_g_per_mps == pytest.approx(
            100.0 * 2.0 * default_cyl.omega / G_STANDARD
        )
        assert s.spin_accel == 0.0

    def test_one_percent_weight_speeds_solve_the_weight_formula(
        self, cyl: Cylinder
    ) -> None:
        s = cyl.summary()
        w = cyl.rim_speed
        heavy = (w + s.plus_one_pct_speed) ** 2 / w**2
        light = (w - s.minus_one_pct_speed) ** 2 / w**2
        assert heavy == pytest.approx(1.01, rel=1e-12)
        assert light == pytest.approx(0.99, rel=1e-12)
        assert s.weightless_speed == cyl.rim_speed


class TestGeometry:
    @given(angles, st.floats(-1e4, 1e4), fractions)
    @example(0.0, 0.0, 0.0)
    @example(1.0, -500.0, 0.999)
    @example(-2.0, 0.0, 0.5)
    def test_floor_point_and_height_round_trip(
        self, theta: float, z: float, frac: float
    ) -> None:
        cyl = Cylinder()
        h = frac * cyl.radius
        p = cyl.floor_point(theta, z, h)
        assert math.hypot(p[0], p[1]) == pytest.approx(
            cyl.radius - h, rel=1e-12, abs=1e-9
        )
        assert p[2] == z
        assert cyl.height_of(p) == pytest.approx(h, abs=1e-9)

    @pytest.mark.parametrize("theta", [0.0, 0.7, math.pi / 2, math.pi, -2.0, 5.5])
    def test_local_basis_is_orthonormal(self, cyl: Cylinder, theta: float) -> None:
        basis = cyl.local_basis(cyl.floor_point(theta, 12.0, 0.0))
        np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=1e-12)

    @given(angles, st.floats(-1e4, 1e4), fractions)
    def test_local_basis_orthonormal_and_left_handed_property(
        self, theta: float, z: float, frac: float
    ) -> None:
        cyl = Cylinder()
        basis = cyl.local_basis(cyl.floor_point(theta, z, frac * cyl.radius))
        np.testing.assert_allclose(basis @ basis.T, np.eye(3), atol=1e-12)
        assert np.linalg.det(basis) == pytest.approx(-1.0, abs=1e-12)

    @pytest.mark.parametrize("theta", [0.0, 0.7, math.pi / 2, math.pi, -2.0, 5.5])
    def test_local_basis_directions_have_physical_meaning(
        self, cyl: Cylinder, theta: float
    ) -> None:
        p = cyl.floor_point(theta, -40.0, 0.3 * cyl.radius)
        spin, axial, up = cyl.local_basis(p)
        on_axis = p + up * math.hypot(p[0], p[1])
        assert math.hypot(on_axis[0], on_axis[1]) == pytest.approx(0.0, abs=1e-9)
        v_floor = np.cross(cyl.omega * Z_HAT, p)
        np.testing.assert_allclose(spin, v_floor / np.linalg.norm(v_floor), atol=1e-12)
        np.testing.assert_array_equal(axial, Z_HAT)

    def test_local_frame_is_left_handed(self, default_cyl: Cylinder) -> None:
        basis = default_cyl.local_basis(default_cyl.floor_point(0.3))
        assert np.linalg.det(basis) == pytest.approx(-1.0, abs=1e-12)

    def test_local_basis_does_not_depend_on_axial_position(
        self, default_cyl: Cylinder
    ) -> None:
        b1 = default_cyl.local_basis(default_cyl.floor_point(1.1, -9000.0))
        b2 = default_cyl.local_basis(default_cyl.floor_point(1.1, +9000.0))
        np.testing.assert_array_equal(b1, b2)

    def test_local_basis_accepts_a_stack_of_points(self, default_cyl: Cylinder) -> None:
        pts = np.array([default_cyl.floor_point(t) for t in (0.0, 1.0, 2.0)])
        stacked = default_cyl.local_basis(pts)
        assert stacked.shape == (3, 3, 3)
        for i in range(3):
            np.testing.assert_array_equal(stacked[i], default_cyl.local_basis(pts[i]))

    @pytest.mark.parametrize("point", [(0.0, 0.0, 0.0), (0.0, 0.0, 123.0)])
    def test_local_basis_is_undefined_on_the_axis(
        self, default_cyl: Cylinder, point: tuple[float, float, float]
    ) -> None:
        with pytest.raises(AxisSingularityError):
            default_cyl.local_basis(np.array(point))

    def test_axis_error_is_a_value_error(self) -> None:
        assert issubclass(AxisSingularityError, ValueError)

    def test_height_of_handles_stacks(self, default_cyl: Cylinder) -> None:
        pts = np.array([default_cyl.floor_point(0.0, 0.0, h) for h in (0, 10, 20)])
        np.testing.assert_allclose(default_cyl.height_of(pts), [0, 10, 20], atol=1e-9)


class TestInertialToHabitat:
    def test_identity_at_time_zero(self, flight_cyl: Cylinder) -> None:
        pt = np.array([123.0, -456.0, 78.0])
        np.testing.assert_array_equal(flight_cyl.inertial_to_habitat(pt, 0.0), pt)

    def test_quarter_turn(self, default_cyl: Cylinder) -> None:
        pt = np.array([default_cyl.radius, 0.0, 42.0])
        out = default_cyl.inertial_to_habitat(pt, default_cyl.period / 4.0)
        np.testing.assert_allclose(out, [0.0, -default_cyl.radius, 42.0], atol=1e-9)

    @given(st.floats(0.0, 2000.0), st.floats(-5e3, 5e3), st.floats(-5e3, 5e3))
    @example(17.0, 1200.0, -700.0)
    def test_is_a_rotation_about_the_axis(self, t: float, x: float, y: float) -> None:
        cyl = Cylinder.from_gravity(spin_accel=1.0e-4)
        pt = np.array([x, y, 55.0])
        out = cyl.inertial_to_habitat(pt, t)
        assert math.hypot(out[0], out[1]) == pytest.approx(
            math.hypot(x, y), rel=1e-12, abs=1e-12
        )
        assert out[2] == 55.0

    def test_accepts_a_stack_of_points(self, default_cyl: Cylinder) -> None:
        pts = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        out = default_cyl.inertial_to_habitat(pts, 7.0)
        for i in range(2):
            np.testing.assert_allclose(
                out[i], default_cyl.inertial_to_habitat(pts[i], 7.0)
            )
