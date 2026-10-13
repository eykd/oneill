"""Plain-text rendering of reports.

Every ``render_*`` function takes a value object from ``oneill.motion`` and
returns a string; nothing here prints.
"""

from typing import TYPE_CHECKING

from oneill.core.cylinder import G_STANDARD

if TYPE_CHECKING:
    from collections.abc import Iterable

    from oneill.core.cylinder import HabitatSummary
    from oneill.motion.scenarios import (
        BallisticsReport,
        CustomReport,
        ElevatorReport,
        LateralReport,
        StairsReport,
        VerticalReport,
    )
    from oneill.motion.verification import CheckResult

_RULE = "=" * 110
_KMH_PER_MPS = 3.6


def fmt3(values: Iterable[float], p: int = 3, w: int = 8) -> str:
    """Format a vector as ``(+a.aaa, +b.bbb, ...)`` with signs and fixed width.

    Values that round to zero print as ``+0`` rather than ``-0``.
    """
    parts = []
    for x in values:
        value = 0.0 if abs(x) < 0.5 * 10 ** (-p) else x
        parts.append(f"{value:+{w}.{p}f}")
    return "(" + ",".join(parts) + ")"


def section(title: str) -> str:
    """A banner heading, preceded by a blank line."""
    return f"\n{_RULE}\n{title}\n{_RULE}\n"


def _join(lines: Iterable[str]) -> str:
    """Join print-style lines (each may contain embedded newlines)."""
    return "".join(f"{line}\n" for line in lines)


def render_summary(s: HabitatSummary) -> str:
    """The habitat parameter block."""
    lines = [
        f"radius               {s.radius:,.1f} m",
        f"length               {s.length:,.1f} m",
        (
            f"spin rate            {s.omega:.6f} rad/s = {s.rpm:.4f} rpm "
            f"(period {s.period:.1f} s)"
        ),
        (
            f"rim speed            {s.rim_speed:.1f} m/s "
            f"({s.rim_speed * _KMH_PER_MPS:.0f} km/h)"
        ),
        (
            f"floor gravity        {s.surface_gravity:.4f} m/s^2 = "
            f"{s.surface_gravity / G_STANDARD:.4f} g"
        ),
        (
            f"gravity gradient     -{s.gradient_pct_per_m:.4f} % per metre of height "
            f"(-{100.0 * s.gradient_pct_per_m:.2f} % at 100 m, 0 on the axis)"
        ),
        (
            f"Coriolis scale       2*omega = {s.coriolis_scale:.5f} m/s^2 per m/s = "
            f"{s.coriolis_pct_g_per_mps:.3f} %g per m/s"
        ),
        (
            f"+1% weight at        {s.plus_one_pct_speed:.2f} m/s spinward;  "
            f"-1% at {s.minus_one_pct_speed:.2f} m/s antispinward"
        ),
        (
            f"weightless at        {s.weightless_speed:.1f} m/s antispinward "
            "(you cancel the spin)"
        ),
    ]
    if s.spin_accel:
        lines.append(
            f"spin acceleration    {s.spin_accel:+.3e} rad/s^2 (Euler force active)"
        )
    return section("O'NEILL CYLINDER PARAMETERS") + _join(lines)


def render_lateral(report: LateralReport) -> str:
    """Section 1: moving laterally on a level floor."""
    g = report.surface_gravity
    mass = report.mass
    lines = [
        (
            f"Vectors are (spinward, axial, up). %g = percent of floor gravity "
            f"({g:.3f} m/s^2). Forces for a {mass:g} kg person."
        ),
        (
            "'curv' = extra weight from following the curved floor (v_spin^2/r); "
            "'tilt' = lean of the apparent vertical.\n"
        ),
        (
            f"{'hdg':>4}  {'direction':<16} {'speed':>6}  "
            f"{'Coriolis accel [%g]':^28}  {'Coriolis force [N]':^28}  "
            f"{'curv':>6}  {'felt wt':>8}  {'tilt':>5}"
        ),
        (
            f"{'deg':>4}  {'':<16} {'m/s':>6}  {'(spin, axial, up)':^28}  "
            f"{'(spin, axial, up)':^28}  {'%g':>6}  {'%g':>8}  {'deg':>5}"
        ),
    ]
    for group in report.groups:
        lines.extend(
            f"{r.heading_deg:4.0f}  {r.direction:<16} {r.speed:6.2f}  "
            f"{fmt3(r.coriolis_pct.to_array())}  "
            f"{fmt3(r.coriolis_force_n.to_array(), p=2)}  "
            f"{r.curvature_pct:6.3f}  {r.felt_weight_pct:8.3f}  {r.tilt_deg:5.2f}"
            for r in group
        )
        lines.append("")
    lines.append(
        f"TAKEAWAY: level motion gives a purely VERTICAL Coriolis force: each m/s "
        f"spinward adds ~{report.coriolis_pct_per_mps:.2f}% to your weight, each m/s "
        "antispinward removes it, axial motion gives none. Exactly, weight = "
        "(omega*r + v_spin)^2 / r. There is no sideways deflection while walking."
    )
    title = f"1. MOVING LATERALLY ON A LEVEL FLOOR (height {report.height:g} m)"
    return section(title) + _join(lines)


def render_vertical(report: VerticalReport) -> str:
    """Section 2: ascending and descending vertically."""
    lines = [
        (
            "rate + = up (toward the axis). 'local g' = centrifugal gravity at that "
            "height. 'push' = direction you are shoved; the ladder/car wall must push "
            "back.\n"
        ),
        (
            f"{'height':>7} {'rate':>6}  {'local g':>7}  {'Coriolis accel [%g]':^28}  "
            f"{'Coriolis force [N]':^28}  {'felt wt':>7}  {'tilt':>6}  push"
        ),
        (
            f"{'m':>7} {'m/s':>6}  {'%g':>7}  {'(spin, axial, up)':^28}  "
            f"{'(spin, axial, up)':^28}  {'%g':>7}  {'deg':>6}"
        ),
    ]
    for group in report.groups:
        if group.skipped:
            lines.append(f"(skipping height {group.height}: outside 0..radius)")
            continue
        lines.extend(
            f"{r.height:7.0f} {r.rate:+6.2f}  {r.local_g_pct:7.2f}  "
            f"{fmt3(r.coriolis_pct.to_array())}  "
            f"{fmt3(r.coriolis_force_n.to_array(), p=2)}  "
            f"{r.felt_weight_pct:7.2f}  {r.tilt_deg:6.2f}  "
            f"{r.push} ({r.push_force_n:.1f} N)"
            for r in group.rows
        )
        lines.append("")
    lines.append(
        f"TAKEAWAY: each m/s of climbing pushes you SPINWARD by "
        f"~{report.coriolis_pct_per_mps:.2f}% g (descending: antispinward). "
        "The push does not weaken with height, but local gravity does, so near the "
        "axis the apparent vertical tilts strongly."
    )
    title = (
        "2. ASCENDING / DESCENDING VERTICALLY (ladder, elevator, rope) "
        "at constant speed"
    )
    return section(title) + _join(lines)


def render_stairs(report: StairsReport) -> str:
    """Section 3: stairs."""
    lines = [
        (
            f"{'hdg':>4}  {'direction':<16} {'sense':<5} {'velocity [m/s]':^28}  "
            f"{'Coriolis accel [%g]':^28}  {'felt wt':>7}  "
            f"{'lateral push [%g]':^19}  {'tilt':>5}"
        ),
        (
            f"{'':>4}  {'':<16} {'':<5} {'(spin, axial, up)':^28}  "
            f"{'(spin, axial, up)':^28}  {'%g':>7}  {'(spin, axial)':^19}  {'deg':>5}"
        ),
    ]
    for group in report.groups:
        lines.extend(
            f"{r.heading_deg:4.0f}  {r.direction:<16} {r.sense:<5} "
            f"{fmt3(r.velocity.to_array())}  {fmt3(r.coriolis_pct.to_array())}  "
            f"{r.felt_weight_pct:7.3f}  {fmt3(r.lateral_pct)}  {r.tilt_deg:5.2f}"
            for r in group
        )
        lines.append("")
    lines.append(
        "TAKEAWAY: stairs combine both effects. The horizontal part changes your "
        "weight (spinward heavier); the vertical part pushes you sideways "
        "(up -> spinward, down -> antispinward)."
    )
    title = (
        f"3. STAIRS: {report.speed:g} m/s along a {report.slope_deg:g} deg straight "
        f"staircase (evaluated at {report.height:g} m)"
    )
    return section(title) + _join(lines)


def render_elevator(report: ElevatorReport) -> str:
    """Section 4: elevator rides."""
    lines: list[str] = []
    for ride in report.rides:
        lines.append(f"\n{ride.label} ride, duration {report.duration:.1f} s")
        lines.append(
            f"{'t[s]':>7} {'height':>8} {'v_up':>7} {'a_up':>6} {'local g%':>9} "
            f"{'felt wt%':>9} {'spin push%':>11} {'side N':>8} {'tilt deg':>9}"
        )
        lines.extend(
            f"{r.t:7.1f} {r.height:8.1f} {r.v_up:+7.2f} {r.a_up:+6.2f} "
            f"{r.local_g_pct:9.2f} {r.felt_weight_pct:9.2f} "
            f"{r.spin_push_pct:+11.3f} {r.side_force_n:+8.1f} {r.tilt_deg:9.2f}"
            for r in ride.rows
        )
        lines.append(
            f"  max tilt {ride.max_tilt_deg:.1f} deg at {ride.max_tilt_height:.0f} m; "
            f"max side force {ride.max_side_force_n:.1f} N; "
            f"felt weight range {ride.weight_min_pct:.1f}..{ride.weight_max_pct:.1f} %g"
        )
    title = (
        f"4. ELEVATOR RIDE floor <-> {report.top:g} m  (v_max {report.v_max:g} m/s, "
        f"a_max {report.a_max:g} m/s^2, {report.mass:g} kg rider)"
    )
    return section(title) + _join(lines)


def render_ballistics(report: BallisticsReport) -> str:
    """Section 5: drops, jumps and throws."""
    lines = [
        (
            "Landing = along-floor distance from the point under the launch, "
            "(spinward, axial) in metres. 'naive' = flat ground with uniform floor "
            "gravity.\n"
        ),
        (
            f"{'case':<32} {'h0':>6} {'T exact':>8} {'T naive':>8}  "
            f"{'landing exact':^20}  {'landing naive':^20}  {'deviation':^20}"
        ),
    ]
    lines.extend(
        f"{r.label:<32} {r.start_height:6.1f} {r.flight_time:8.3f} "
        f"{r.naive_time:8.3f}  {fmt3(r.landing, p=3, w=9)}  "
        f"{fmt3(r.naive_landing, p=3, w=9)}  {fmt3(r.deviation, p=3, w=9)}"
        + ("  [hits endcap]" if r.hits_endcap else "")
        for r in report.rows
    )
    lines.append(
        "\nTAKEAWAY: falling objects land antispinward of the point below them; a "
        "vertical jump lands slightly spinward."
    )
    title = (
        "5. FREE FLIGHT: drops, jumps, throws "
        "(exact straight-line solution in the inertial frame)"
    )
    return section(title) + _join(lines)


def render_custom(report: CustomReport) -> str:
    """Section 6: custom velocities."""
    lines: list[str] = []
    for e in report.entries:
        lines.extend(
            [
                f"\nvelocity {fmt3(e.velocity.to_array(), p=2)} m/s",
                f"  centrifugal      {fmt3(e.centrifugal_pct.to_array())}",
                (
                    f"  Coriolis         {fmt3(e.coriolis_pct.to_array())}   = "
                    f"{fmt3(e.coriolis_force_n.to_array(), p=2)} N"
                ),
                f"  Euler            {fmt3(e.euler_pct.to_array())}",
                f"  -path accel      {fmt3(e.neg_path_accel_pct.to_array())}",
                (
                    f"  APPARENT gravity {fmt3(e.apparent_pct.to_array())}  "
                    f"felt weight {e.felt_weight_pct:.3f} %g, tilt {e.tilt_deg:.2f} deg"
                ),
                f"  support force    {fmt3(e.support_force_n.to_array(), p=2)} N",
            ]
        )
    title = (
        f"6. CUSTOM VELOCITIES at height {report.height:g} m "
        "(components spinward, axial, up; %g)"
    )
    return section(title) + _join(lines)


def render_selftest(results: Iterable[CheckResult]) -> str:
    """The self-test results, one line per check, with a verdict."""
    checks = list(results)
    lines = [
        f"  [{'PASS' if r.passed else 'FAIL'}] {r.name}: "
        f"got {r.got:.7g}, expected {r.want:.7g}"
        for r in checks
    ]
    verdict = "ALL PASSED" if all(r.passed for r in checks) else "SOME CHECKS FAILED"
    lines.append(f"\n{verdict}")
    return section("SELF-TEST") + _join(lines)
