# L3 design: `atmo/annulus2d.py`

Status: design only, for physics review before any solver code. Tests H1 (docs/research/atmosphere-literature.md Q7, Q8, Q10): does axis-invariant (2D) convection driven by three-fold strip heating build large eddies, or lock to the strips?

## 1. Equations

Domain: annulus r in [r_in, R], rotating frame, Omega along z (the axis), flow invariant in z. Streamfunction u_r = psi_theta/r, u_theta = -psi_r, zeta = -laplacian(psi). Effective gravity is centrifugal, g(r) = Omega^2 r, **outward**.

Boussinesq density rho = rho0(1 - beta T'), where T' is the potential-temperature anomaly (adiabatic background). The buoyancy force per mass is the centrifugal force on the anomaly, minus the pressure gradient absorbed into p': F = -beta T' g(r) r_hat. Warm air feels a force toward the axis ("up").

**Baroclinic torque.** For F = F_r r_hat, (curl F)_z = -(1/r) dF_r/dtheta. With F_r = -beta g T' and g/r = Omega^2:

    (curl F)_z = + beta Omega^2 dT'/dtheta

Check: T' larger at larger theta means inward flow there and outward flow at smaller theta, so u_r falls with theta and zeta = -(1/r) d u_r/dtheta > 0. The signs agree, and the test suite pins this with a two-sided heating case.

**Coriolis drops out.** In 2D incompressible flow, curl(Omega x u) = Omega(div u) - (Omega.grad)u = 0. Coriolis is a pure gradient absorbed by p'. This is the Taylor-Proudman statement of Q8.

**Optional non-Boussinesq Coriolis-buoyancy.** Keeping (rho'/rho0) on 2 Omega x u adds, after curl (using curl(phi A) with A = z x u, curl A = 0): + 2 Omega beta u.grad T'. Relative to the baroclinic term its size is 2u/(Omega r), about 3 % for 3 m/s at the rim and 1 % per m/s. It ships off by default and is a flag, with a test that turning it on changes results by about that amount.

**System**

    d zeta/dt + u.grad zeta = nu_t lap zeta + beta Omega^2 dT'/dtheta - zeta/tau_E [+ 2 Omega beta u.grad T']
    d T'/dt   + u.grad T'   = kappa_t lap T' - (<T'> / tau_T) [- (N^2 / (beta g)) u_r, default N = 0]

nu_t = 1 or 10 m^2/s (as L1), kappa_t = nu_t/Pr_t, Pr_t = 1. The N^2 term is the isothermal-reference stratification (off for the adiabatic default).

## 2. Drag

**End caps.** An Ekman layer at each cap pumps w_E = (1/2) sqrt(nu/Omega) zeta (delta_E = sqrt(2 nu / f), f = 2 Omega). With suction at both caps the interior w runs from -w_E to +w_E over L, so dw/dz = 2 w_E / L and vortex stretching gives

    d zeta/dt = 2 Omega (2 w_E / L) = -(2 sqrt(nu Omega) / L) zeta = -zeta / tau_E

with tau_E = (L/2)/sqrt(nu Omega). **The factor is 1/tau_E, not 2/tau_E**: each cap contributes sqrt(nu Omega)/L, and the two caps over the full length L sum to sqrt(nu Omega)/(L/2). This matches L1 (11.8 h at nu = 1) and Q8's 2 sqrt(nu Omega)/L. Option: turbulent C_D U/(L/2). The linear form needs Ro <~ 0.3 (as L1).

**Floor (canopy).** Acts only at r = R as a wall stress, not a body force. Tangential stress balance nu r d(u_theta/r)/dr = -C_D |u_ref| u_ref, with C_D = 0.0178 (L1 default, z_r = 10 m, u_ref interpolated at R - z_r), becomes a **Dirichlet vorticity condition** at the wall (u_r = 0 there):

    zeta_R = 2 u_theta,R / R - C_D |u_ref| u_ref / nu_t

Free slip is the C_D = 0 case; no-slip (Thom) is kept for the onset test only. The condition is lagged one step with |u_ref| frozen (Picard), as in L1. The stress term is stiff (about 0.02 /s against zeta of order 0.01), so it is applied through a Robin form in the implicit solve if the lagged version is unstable.

The m = 0 mode is advanced as the azimuthal-wind equation (the L1 operator plus Reynolds-stress forcing) because psi_0 has no unique boundary value.

## 3. Heating

Heat flux H (W/m^2) enters the air at r = R: kappa_t dT'/dr = H/(rho c_p) (positive H gives warm air near the floor, unstable). At half an atmosphere rho c_p = 603, so H = 200 gives 0.33 K m/s.

- **Pattern:** three land strips (60 degrees each) alternate with three window strips, fixed in the hull frame. Land H_L(t) = 200 sin(2 pi t/24 h) by day (H = 100-350 range, Q10), -30 by night (-10 to -60); steady option H_L = 100. Windows are H = 0 by default, with an optional radiative loss h_w.
- **Stationarity:** relaxation of the domain mean, -<T'>/tau_T with tau_T = 1 day. A net-flux-neutral mode (windows balance the land mean) is the alternative; both are tested.
- **Inner wall:** insulating (dT'/dr = 0).

## 4. Numerics

- **Grid:** rfft in theta (N_theta = 512, |m| <= 170 after 2/3 dealiasing), second-order finite volume in r on a tanh-stretched grid (N_r = 80, floor cell about 1 m, r_in = 0.1 R).
- **Poisson:** per Fourier mode a tridiagonal solve of (d_rr + d_r/r - m^2/r^2) psi_m = -zeta_m, with psi_m = 0 at both walls; one batched Thomas over modes.
- **Time stepping:** IMEX SBDF2. Diffusion and linear drag implicit (tridiagonal per mode), advection, baroclinic term and wall vorticity explicit. The Jacobian is flux-form (finite volume in r, spectral in theta), conserving the integral of zeta exactly and energy and enstrophy to discretization order.
- **Polar filter:** the high-m truncation tapers with r so the physical spacing r Delta theta stays at least 2 pi R / N_theta, otherwise the inner cell (4 m) sets the CFL.
- **CFL:** dt = 0.4 min(r Delta theta/|u_theta|, Delta r/|u_r|), about 6 s at 3 m/s.
- **Inner boundary:** free-slip insulating at r_in = 0.1 R. This is an artifact, since the real habitat has no inner wall (Zhong & Sun 2026 found it matters: strong zonal flow, Nu ~ Ra^0.1). Sensitivity runs at 0.05 R and 0.2 R are part of the deliverable.
- **Cost (numpy):** about 3 ms/step, so about 1.5e4 steps and 45 s per simulated day; ten days about 8 min. **Tests use 24 x 64 grids** (under 1 s each); long runs are a script and the CLI (`oneill atmosphere --run-l3`), never a test.

## 5. Validation tests

Each fails if its ingredient is wrong.

1. **Poisson vs analytic mode:** psi_m = sin(pi (r - r_in)/(R - r_in)) r^m-type; error falls 4x on halving h.
2. **Lamb-Oseen:** m = 0 vortex in a large domain with analytic Dirichlet boundaries; compare zeta(r, t) with Gamma exp(-r^2/4 nu t)/(4 pi nu t) to 1e-3.
3. **Axisymmetric heating gives zero flow:** theta-independent H leaves every m >= 1 mode below 1e-12 and T' equal to the 1-D conduction solution.
4. **Linear drag:** advection off, only the end-cap term: amplitude decays as exp(-t/tau_E), enstrophy at 2/tau_E, and doubling L halves the rate.
5. **Budgets:** with no T', the inviscid, drag-free run conserves energy and enstrophy to O(dt^2); with dissipation, dE/dt = buoyancy work - viscous - drag closes to 1e-3.
6. **Torque sign:** a single warm patch makes zeta of the sign derived in 1.
7. **Convection onset:** thin gap eta = 0.99, rigid walls, fixed-T walls, Ra from the mean radius. Linearized single-mode growth rate sigma(Ra) at the critical m (about 2 pi r_m / (2.016 L)), onset where sigma crosses zero. Target Ra_c = 1708 within 3 % (research correction: no E^(-4/3) law). Also check Ra_c rises as eta falls (0.8, 0.5; cf. Zhong, Wang & Sun 2024, n_c = 5 at eta = 0.5).

## 6. Outputs answering H1

- **Eddy scale vs time:** peak and centroid of the azimuthal energy spectrum for r in [0.5 R, R], l = 2 pi r/m.
- **Wind:** rms u at 10 m, and gust = mean + 3 sigma (the 2D model has no turbulent gust).
- **Heat exchange:** Nu_eff = H R/(rho c_p kappa_t Delta T), Delta T from the rim layer to the inner wall.
- **Lock or condensate:** E(m = 3)/E_total and its phase against the strips (locked), versus E(m <= 2)/E_total growing with time and l tending to the domain size (inverse cascade).
- **Arrest wavelength:** measure eps (total dissipation) and mu = 1/tau_E (plus wall drag), then compare l with (0.08-0.18)(eps/mu^3)^(1/2). For eps = 1e-4 W/kg and nu = 1, mu = 2.4e-5 gives about 7-15 km, larger than R. So friction cannot arrest the cascade inside the domain and I expect a domain-limited condensate; the run tests this.

## 7. Effort and scope cuts

About 2-3 weeks: solver core and validation (1 week), heating and diagnostics (1 week), runs and write-up (0.5-1 week). Cuts: non-Boussinesq term and stratification stay off by default; steady heating first and the day/night cycle second; one resolution plus a half-resolution check; no 3D, no moist physics, no adaptive stepping. Known weaknesses: constant nu_t, lambda about 0.03 (marginal 2D-ization), r_in artifact, and 2D by construction, so H1 can be supported but the breakdown of axis-invariance cannot be tested.
