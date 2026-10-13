"""Markdown rendering of the atmosphere report, and generated-block substitution.

``render_blocks`` turns an :class:`~oneill.atmo.report.AtmosphereReport` into the
named markdown blocks of ``docs/atmosphere.md``; ``replace_blocks`` rewrites only
the text between ``<!-- gen:NAME -->`` and ``<!-- /gen -->`` markers. Nothing here
touches the file system.
"""

import re
from typing import TYPE_CHECKING

from oneill.atmo.report import (
    DEFLECTION_RISE,
    HUMIDITIES,
    STABILITY_HEIGHTS,
    WIND_HEIGHT,
    WIND_SPEED,
)
from oneill.atmo.spinup import DRAG_HEIGHT
from oneill.atmo.surface import Stability

if TYPE_CHECKING:
    from collections.abc import Callable, Mapping, Sequence

    from oneill.atmo.report import AtmosphereReport

BLOCK_NAMES = ("gas", "equilibrium", "numbers", "surface", "parcel", "spinup")

_MARKER = re.compile(r"<!-- (/gen|gen:[A-Za-z0-9_-]+) -->")
_OPEN_PREFIX = "gen:"
_CLOSE = "/gen"
_CLASS_LETTER = {
    Stability.UNSTABLE: "U",
    Stability.NEUTRAL: "N",
    Stability.STABLE: "S",
}
_HOUR = 3600.0
_DAY = 86400.0
_TWO_HOURS = 2.0 * _HOUR
_TWO_DAYS = 2.0 * _DAY


def _table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    lines = ["| " + " | ".join(headers) + " |"]
    lines.append("|" + "|".join(" --- " for _ in headers) + "|")
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)


def _num(value: float | None, digits: int = 2) -> str:
    if value is None:
        return "-"
    return f"{value:.{digits}f}"


def _sci(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}e}".replace("e-0", "e-").replace("e+0", "e+")


def _duration(seconds: float | None) -> str:
    """Seconds as minutes, hours or days, whichever reads best."""
    if seconds is None:
        return "not reached"
    if seconds < _TWO_HOURS:
        return f"{seconds / 60.0:.0f} min"
    if seconds < _TWO_DAYS:
        return f"{seconds / _HOUR:.1f} h"
    return f"{seconds / _DAY:.0f} d"


def _caveat(report: AtmosphereReport, name: str) -> str:
    return f"*{report.caveats[name]}*"


def render_gas(report: AtmosphereReport) -> str:
    """Gas properties of both presets."""
    rows = [
        [
            g.label,
            _num(g.p_dry / 1e3, 1),
            _num(g.p_vapour / 1e3, 1),
            _num(g.molar_mass * 1e3, 2),
            _num(g.r_specific, 1),
            _num(g.cp, 0),
            _num(g.gamma, 3),
            _num(g.viscosity * 1e5, 3),
            _num(g.nu * 1e5, 3),
            _num(g.density, 3),
        ]
        for g in report.gas
    ]
    headers = [
        "Preset",
        "p dry (kPa)",
        "p vapour (kPa)",
        "M (g/mol)",
        "R_s (J/kg/K)",
        "c_p (J/kg/K)",
        "gamma",
        "mu (1e-5 Pa s)",
        "nu (1e-5 m2/s)",
        "rho (kg/m3)",
    ]
    return f"{_table(headers, rows)}\n\n{_caveat(report, 'gas')}"


def render_equilibrium(report: AtmosphereReport) -> str:
    """Pressure, temperature and cloud-base numbers."""
    rows = [
        [
            e.label,
            _num(e.phi_axis / 1e3, 2),
            _num(e.axis_ratio_isothermal, 3),
            _num(e.axis_ratio_adiabatic, 3),
            _num(e.dp_isothermal / 1e3, 1),
            _num(e.dp_adiabatic / 1e3, 1),
            _num(e.delta_t_adiabatic, 1),
            _num(e.floor_lapse_rate * 1e3, 2),
            _num(e.mass_isothermal / 1e11, 2),
        ]
        for e in report.equilibrium
    ]
    main = _table(
        [
            "Preset",
            "phi axis (kJ/kg)",
            "p_axis/p_rim isothermal",
            "p_axis/p_rim adiabatic",
            "dp isothermal (kPa)",
            "dp adiabatic (kPa)",
            "dT adiabatic (K)",
            "floor lapse (K/km)",
            "air mass (1e11 kg)",
        ],
        rows,
    )
    clouds = [
        [
            c.label,
            f"{rh:.0%}",
            _num(c.temperature - 273.15, 1),
            "above the axis" if c.height is None else _num(c.height, 0),
            _num(c.earth_equivalent_height, 0),
        ]
        for e in report.equilibrium
        for rh, c in zip(e.humidities, e.clouds, strict=True)
    ]
    cloud_table = _table(
        [
            "Gas at floor",
            "RH",
            "T at cloud base (degC)",
            "cloud base, true height (m)",
            "cloud base, Earth-equivalent (m)",
        ],
        clouds,
    )
    return f"{main}\n\n{cloud_table}\n\n{_caveat(report, 'equilibrium')}"


def render_numbers(report: AtmosphereReport) -> str:
    """Dimensionless numbers and viscous scales."""
    first = report.numbers[0].table
    head = _table(
        [
            "Gas at floor",
            "2 omega (1/s)",
            "vs Earth",
            "inertial period (s)",
            "Ro (U=3 m/s, L=1 km)",
            "Ro (D=10 m)",
            "N (1/s)",
            "Rossby radius (m)",
            "shedding f (Hz)",
            "shedding / 2 omega",
        ],
        [
            [
                n.label,
                _num(n.table.coriolis_parameter, 4),
                f"{n.coriolis_vs_earth:.0f}x",
                _num(n.table.inertial_period, 1),
                _num(n.table.rossby, 3),
                _num(n.table.obstacle_rossby, 2),
                _num(n.table.brunt_vaisala, 4),
                _num(n.table.rossby_radius, 0),
                _num(n.table.shedding_frequency, 3),
                _num(n.table.shedding_to_coriolis, 1),
            ]
            for n in report.numbers
        ],
    )
    viscous = [(f"molecular, {n.label}", n.table.molecular) for n in report.numbers] + [
        (f"eddy, nu_t = {s.nu:g}", s) for s in first.eddy
    ]
    scales = _table(
        [
            "Viscosity",
            "nu (m2/s)",
            "Ekman number",
            "Ekman depth (m)",
            "E^1/4 H (m)",
            "E^1/3 H (m)",
            "spin-up time",
        ],
        [
            [
                label,
                _sci(s.nu),
                _sci(s.ekman_number),
                _num(s.ekman_depth, 3),
                _num(s.stewartson_quarter, 1),
                _num(s.stewartson_third, 1),
                _duration(s.spin_up_time),
            ]
            for label, s in viscous
        ],
    )
    return f"{head}\n\n{scales}\n\n{_caveat(report, 'numbers')}"


def render_surface(report: AtmosphereReport) -> str:
    """Wind loads and the rotational stability of the surface layer."""
    caption = (
        f"Mean wind {WIND_SPEED:g} m/s at {WIND_HEIGHT:g} m. "
        "Unstable layer: spinward wind, heights where dU/dz > 2 omega."
    )
    loads = _table(
        [
            "Terrain",
            "Gas at floor",
            "z0 (m)",
            "d (m)",
            "h_c (m)",
            "u* (m/s)",
            "gust (m/s)",
            "gust factor",
            "q mean (Pa)",
            "q gust (Pa)",
            "unstable layer (m)",
        ],
        [
            [
                s.terrain.name,
                s.gas_label,
                _num(s.terrain.z0, 1),
                _num(s.terrain.d, 1),
                _num(s.terrain.h_c, 1),
                _num(s.ustar, 3),
                _num(s.gust_speed, 2),
                _num(s.gust_factor, 2),
                _num(s.q_mean, 1),
                _num(s.q_gust, 1),
                "none" if s.layer is None else f"{s.layer[0]:.1f} to {s.layer[1]:.1f}",
            ]
            for s in report.surface
        ],
    )
    heights = [f"{z:g} m" for z in STABILITY_HEIGHTS]
    stability = _table(
        ["Terrain", "Wind toward", *heights],
        [
            [
                row.terrain,
                row.heading,
                *(_CLASS_LETTER[c] for c in row.classes),
            ]
            for row in report.stability
        ],
    )
    legend = "U unstable (shear turbulence augmented), S stable, N neutral."
    return (
        f"{caption}\n\n{loads}\n\n{stability}\n\n{legend}\n\n"
        f"{_caveat(report, 'surface')}"
    )


def render_parcel(report: AtmosphereReport) -> str:
    """Buoyant parcel numbers: short-run and long-time behaviour."""
    caption = (
        f"Warm parcel released at the floor in dry-adiabatic air (O'Neill air, "
        f"{report.gas[0].temperature:.0f} K). Short run 300 s; deflection read at "
        f"{DEFLECTION_RISE:g} m of rise; cloud base at RH {HUMIDITIES[0]:.0%}, "
        "computed from a long run where shown."
    )
    rows = [
        [
            f"{p.temperature_excess:g}",
            _num(p.drag, 3),
            _num(p.max_rise, 1),
            _num(p.peak_speed, 2),
            _num(p.deflection, 2),
            _num(p.drift_speed, 3),
            _num(p.steady_rise, 3),
            _num(p.steady_drift, 2),
            _duration(p.time_to_lcl) if p.time_to_lcl is not None else "-",
        ]
        for p in report.parcels
    ]
    table = _table(
        [
            "dT (K)",
            "drag k (1/s)",
            "max rise in 300 s (m)",
            "peak speed (m/s)",
            "deflection at 2 m rise (m)",
            "drift at 2 m rise (m/s)",
            "steady rise (m/s)",
            "steady drift (m/s)",
            "time to cloud base",
        ],
        rows,
    )
    return f"{caption}\n\n{table}\n\n{_caveat(report, 'parcel')}"


def _band(times: Sequence[float | None]) -> str:
    """Range of durations over a closure sweep, e.g. ``22.8-24.0 h``."""
    values = [t for t in times if t is not None]
    if not values:
        return "not reached"
    lo, hi = _duration(min(values)), _duration(max(values))
    if lo == hi:
        return lo
    unit = lo.split(" ")[-1]
    if hi.endswith(unit):
        return f"{lo.rsplit(' ', 1)[0]}-{hi}"
    return f"{lo} - {hi}"


def _validity(endcap: bool, valid: bool) -> str:  # noqa: FBT001
    if not endcap:
        return "n/a (no end caps)"
    return "yes" if valid else "no (order of magnitude)"


def render_spinup(report: AtmosphereReport) -> str:
    """Spin-up, spin-down and lag-wind numbers."""
    floor = report.spin[0].floor
    spin = _table(
        [
            "Hull ramp",
            "nu_t (m2/s)",
            "End caps",
            "time to 90 % co-rotation",
            f"peak wind at {DRAG_HEIGHT:g} m (m/s)",
            "peak wind anywhere (m/s)",
            "peak epsilon",
            "linear cap sink valid",
        ],
        [
            [
                _duration(r.ramp),
                f"{r.nu:g}",
                "on" if r.endcap else "off",
                _duration(r.time_to_90),
                _num(r.peak_wind_ref, 1),
                _num(r.peak_wind_interior, 1),
                _num(r.max_epsilon, 2),
                _validity(r.endcap, r.sink_valid),
            ]
            for r in report.spin
        ],
    )
    sudden = _table(
        [
            "Hull spin change",
            "nu_t stable / unstable",
            "End caps",
            "time to 90 % adjustment",
            "epsilon",
            "linear cap sink valid",
        ],
        [
            [
                f"{s.delta_fraction:+.0%}",
                f"{s.nu_stable:g} / {s.nu_unstable:g}",
                "on" if s.endcap else "off",
                _duration(s.time_to_adjust),
                _num(s.max_epsilon, 2),
                _validity(s.endcap, s.sink_valid),
            ]
            for s in report.sudden
        ],
    )
    widths = report.blend[0].widths
    blend = _table(
        [
            "Hull spin change (end caps on)",
            (
                "time to 90 % adjustment, band over blend width "
                f"{widths[0]:g}-{widths[-1]:g}"
            ),
            "epsilon",
            "linear cap sink valid",
        ],
        [
            [
                f"{b.delta_fraction:+.0%}",
                _band(b.times),
                _num(b.max_epsilon, 2),
                "yes" if b.sink_valid else "no (order of magnitude)",
            ]
            for b in report.blend
        ],
    )
    lag = _table(
        [
            "nu_t (m2/s)",
            "End caps",
            "peak lag wind (m/s)",
            f"lag wind at {DRAG_HEIGHT:g} m (m/s)",
        ],
        [
            [
                f"{r.nu:g}",
                "on" if r.endcap else "off",
                _num(r.peak, 2),
                _num(r.wind_ref, 3),
            ]
            for r in report.lag
        ],
    )
    molecular = _table(
        [
            "Gas at floor",
            "molecular nu (m2/s)",
            "time to 90 % co-rotation (1 h ramp)",
            f"peak wind at {DRAG_HEIGHT:g} m (m/s)",
            "Ekman Reynolds number",
        ],
        [
            [
                m.label,
                _sci(m.nu),
                _duration(m.time_to_90),
                _num(m.peak_wind_ref, 1),
                f"{m.ekman_reynolds:,.0f}",
            ]
            for m in report.molecular
        ],
    )
    parts = [
        (
            f"Floor condition for every row: {floor}. Hull spun from rest to full "
            "rate; the air starts at rest."
        ),
        spin,
        (
            "Sudden change of the hull spin rate, with a Rayleigh-switched viscosity "
            "(larger where angular momentum falls outward). With end caps the sink "
            "rate scales with the new hull rate, so spin-down is the slower of the "
            "pair."
        ),
        sudden,
        (
            "Sensitivity of the time to 90 % adjustment (nu_t 1 / 10, end caps on) "
            "to the blend width of the viscosity closure, as a band over the "
            "widths tried (computed at a 150 s time step). The 50 % step has "
            "epsilon = 1, so its time is order-of-magnitude only."
        ),
        blend,
        (
            "Steady lag wind while the hull changes speed by 1 % per day "
            "(spin-up: the air lags, so the wind is negative)."
        ),
        lag,
        (
            "Molecular viscosity only. Illustrative: the laminar end-cap sink is "
            "invalid at these Ekman Reynolds numbers (turbulent layers), and the "
            "1 h ramp also exceeds epsilon = 0.3."
        ),
        molecular,
        _caveat(report, "spinup"),
    ]
    return "\n\n".join(parts)


_RENDERERS: dict[str, Callable[[AtmosphereReport], str]] = {
    "gas": render_gas,
    "equilibrium": render_equilibrium,
    "numbers": render_numbers,
    "surface": render_surface,
    "parcel": render_parcel,
    "spinup": render_spinup,
}


def render_blocks(report: AtmosphereReport) -> dict[str, str]:
    """Render every named markdown block."""
    return {name: _RENDERERS[name](report) for name in BLOCK_NAMES}


def replace_blocks(document: str, blocks: Mapping[str, str]) -> str:
    """Replace the content of every ``gen`` block in ``document``.

    Text outside the markers is untouched, and the result is stable: replacing a
    second time returns identical text.

    Raises:
        ValueError: for an unknown block name, an unclosed or stray marker, or
            nested markers.
    """
    out: list[str] = []
    cursor = 0
    open_name: str | None = None
    for match in _MARKER.finditer(document):
        token = match.group(1)
        if open_name is None:
            if token == _CLOSE:
                msg = f"closing marker without an opening marker at {match.start()}"
                raise ValueError(msg)
            open_name = token.removeprefix(_OPEN_PREFIX)
            if open_name not in blocks:
                msg = f"unknown block {open_name!r}; known: {', '.join(sorted(blocks))}"
                raise ValueError(msg)
            out.append(document[cursor : match.end()])
            out.append("\n" + blocks[open_name].rstrip("\n") + "\n")
        else:
            if token != _CLOSE:
                msg = f"block {open_name!r} is not closed before {token!r}"
                raise ValueError(msg)
            open_name = None
            cursor = match.start()
    if open_name is not None:
        msg = f"block {open_name!r} is never closed"
        raise ValueError(msg)
    out.append(document[cursor:])
    return "".join(out)
