# oneill

Fictitious forces and free flight inside a spinning O'Neill cylinder: what a
person feels walking, climbing, riding elevators, dropping things and throwing
them in a rotating habitat. Pure numpy; matplotlib only for plots.

## Install and run

```sh
uv sync                    # runtime + dev tools
uv sync --extra plot       # also matplotlib for --plot / --save-plots

uv run oneill                               # every report, 3.2 km habitat, 1 g
uv run oneill --radius 500 --sections lateral vertical
uv run oneill --rpm 0.5 --velocity 3 0 0.5  # your own velocity (spin, axial, up)
uv run oneill --selftest                    # physics checks; exit 1 on failure
uv run oneill --save-plots figs/            # four PNGs
```

`--help` lists every option (speeds, headings, stairs, elevator limits, mass,
`--spin-accel` for spin-up/down, ...). Invalid input exits with status 2 and
writes nothing to stdout.

## Conventions

- The **habitat frame** co-rotates with the cylinder; the spin vector is
  `omega * z_hat`, along the axis.
- A person's **local frame** is written `(spinward, axial, up)`. Spinward is the
  way the floor moves, axial is along the spin vector and `up` points at the
  axis. That triad is left-handed.
- **Headings** are compass-like: 0 = +axial, 90 = spinward, 180 = -axial,
  270 = antispinward.
- Apparent gravity is `centrifugal + Coriolis + Euler - a_path`; the floor,
  ladder or stair supplies `-m` times that.

## Library use

```python
from oneill import Cylinder, FloorWalk, LocalVector, analyze_path, ballistic

cyl = Cylinder.from_gravity(radius=3200, length=20000)  # 1 g, ~0.53 rpm
hist = analyze_path(cyl, FloorWalk(cyl, heading_deg=90, speed=5), [0.0])
print(hist.coriolis[0], hist.felt_weight[0], hist.tilt_deg[0])
# [ 0.  0. -0.5536] 10.368 0.0   -> walking spinward makes you heavier

drop = ballistic(cyl, 100.0, LocalVector(0.0, 0.0, 0.0))
print(drop.flight_time, drop.landing_spin)
# 4.625 -17.23   -> a 100 m drop lands 17 m antispinward of the point below
```

Paths (`FloorWalk`, `StraightPath` with a `ConstantSpeed` or `Trapezoid`
profile, `NumericPath` for any callable) give analytic position, velocity and
acceleration. `integrate_rotating` integrates a free body in the habitat frame
with an optional `extra_acc(t, p, v)` callback for drag or buoyancy, built on
the generic `rk4`.

## Layout

```
src/oneill/
  core/      Cylinder, vectors, rotating-frame accelerations, rk4   (no people, no air)
  motion/    paths, profiles, forces, ballistics, scenarios, figures, self-test
  adapters/  cli, text rendering, matplotlib plots                  (the only I/O)
```

Dependencies point inward only (`adapters -> motion -> core`); import-linter
enforces it. Domain code never prints, and matplotlib is imported only in
`adapters/plots.py` (ruff bans it elsewhere).

## Development

```sh
uv run pre-commit install        # once
uv run pre-commit run --all-files
```

Hooks run ruff (`select = ["ALL"]`), `ruff format`, `mypy --strict`,
import-linter and pytest with a 100 % line-and-branch coverage gate. Tests
check the simulator against independent physics (inertial-frame closed forms,
the Jacobi integral, RK4 against exact solutions) and use Hypothesis for the
invariants. `HYPOTHESIS_PROFILE=ci uv run pytest` runs 300 examples per
property instead of 25.
