import pytest

from oneill.core.cylinder import Cylinder
from oneill.motion.verification import CheckResult, all_passed, run_selftest


class TestCheckResult:
    @pytest.mark.parametrize(
        ("got", "want", "tol", "passed"),
        [
            (1.0, 1.0, 0.0, True),
            (1.0, 1.1, 0.2, True),
            (1.0, 1.1, 0.05, False),
            (-3.0, 3.0, 1.0, False),
        ],
    )
    def test_passed_means_within_tolerance(
        self,
        got: float,
        want: float,
        tol: float,
        passed: bool,  # noqa: FBT001
    ) -> None:
        assert CheckResult("x", got, want, tol).passed is passed

    def test_all_passed_needs_every_check(self) -> None:
        good = CheckResult("good", 1.0, 1.0, 0.1)
        bad = CheckResult("bad", 2.0, 1.0, 0.1)
        assert all_passed([good, good])
        assert not all_passed([good, bad])
        assert all_passed([])


class TestSelfTest:
    def test_default_habitat_passes_every_check(self, default_cyl: Cylinder) -> None:
        results = run_selftest(default_cyl)
        assert len(results) == 10
        assert all_passed(results), [r for r in results if not r.passed]

    def test_check_names_describe_the_physics(self, default_cyl: Cylinder) -> None:
        names = [r.name for r in run_selftest(default_cyl)]
        assert names[0] == "floor weight at -50 m/s spinward == (wR+v)^2/R"
        assert "Coriolis climbing 0.5 m/s == 2wv spinward" in names
        assert sum("rotating-frame RK4 vs exact inertial" in n for n in names) == 3

    def test_small_station_and_spin_up_also_pass(self) -> None:
        station = Cylinder.from_rpm(250.0, 500.0, 2.0)
        assert all_passed(run_selftest(station))
        spinning_up = Cylinder.from_gravity(spin_accel=1e-4)
        assert all_passed(run_selftest(spinning_up))
