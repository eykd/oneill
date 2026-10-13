# L3 design: `atmo/annulus2d.py`

Status: design only, for physics review before any solver code. Tests H1 (docs/research/atmosphere-literature.md Q7, Q8, Q10): does axis-invariant (2D) convection driven by three-fold strip heating build large eddies, or lock to the strips?

## 1. Equations

Domain: annulus r in [r_in, R], rotating frame, Omega along z (the axis), flow invariant in z. Streamfunction u_r = psi_theta/r, u_theta = -psi_r, zeta = -laplacian(psi). Effective gravity is centrifugal, g(r) = Omega^2 r, **outward**.

Boussinesq density rho = rho0(1 - beta T'), where T' is the potential-temperature anomaly (adiabatic background). The buoyancy force per mass is the centrifugal force on the anomaly, minus the pressure gradient absorbed into p': F = -beta T' g(r) r_hat. Warm air feels a force toward the axis ("up").

**Baroclinic torque.** For F = F_r r_hat, (curl F)_z = -(1/r) dF_r/dtheta. With F_r = -beta g T' and g/r = Omega^2:

    (curl F)_z = + beta Omega^2 dT'/dtheta

Check: T' larger at larger theta means inward flow there and outward flow at smaller theta, so u_r falls with theta and zeta = -(1/r) d u_r/dtheta > 0. The signs agree, and the test suite pins this with a two-sided heating case.

**Coriolis drops out.** In 2D incompressible flow, curl(Omega x u) = Omega(div u) - (Omega.grad)u = 0. Coriolis is a pure gradient absorbed by p'. This is the Taylor-Proudman statement of Q8.

**Coriolis-buoyancy (default ON), corrected.** An earlier draft added `+2 Omega beta u.grad T` (the density factor on `2 Omega x u`) and a K0 shift, and predicted a steady retrograde mean wind of 2-3 m/s. That was an artifact. Coriolis acts on the *mass flux* `v = rho u / rho0`, and a closed fluid conserves mass: `div v = -(d rho/dt)/rho0 = beta dT'/dt`, where `T'` is the anomaly about the area mean (total mass is fixed). Then `curl(-2 Omega z x v) = -2 Omega div v`, so the vorticity source is

    -2 Omega beta dT'/dt

built from the model's own temperature tendency (advection, diffusion, wall flux, relaxation), for the m >= 1 modes. For the theta-mean, `(1/r) d(r F)/dr = -2 Omega beta d<T'>/dt`, integrated from the inner wall (zero at the outer wall by the area-mean subtraction): `F(r) = -(2 Omega beta / r) int r' d<T'>/dt dr'`. Equivalently the mean Coriolis torque is `-2 Omega <rho u_r>/rho0`, and mass conservation gives `<rho u_r> = (beta rho0 / r) int r' d<T>/dt dr'`.

Consequences: the term **vanishes in a steady state at every radius** (there is no steady zonal wind from it), and it only moves angular momentum between relative motion and the moment of inertia `I' = -2 pi beta rho0 int T' r^3 dr`: `dL_rel/dt = -Omega dI'/dt`. The conserved quantity is `L_rel + Omega I'` (the code's `absolute_angular_momentum`), which changes only through wall and end-cap torques. The flag switches the whole term on or off; production reports both and the difference should be transient.

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

Free slip is the C_D = 0 case; no-slip (Thom) is kept for the onset test only. **As built:** the stress is applied as L1's body-force layer (-C_D |u_ref| u / z_r over the bottom z_r metres), whose total stress equals the wall-vorticity form. Its curl enters the m >= 1 vorticity equation, and its theta-mean enters the m = 0 equation, with the stress the theta-mean of C_D |u| u, not C_D |ubar| ubar. The vorticity wall condition is then plain free slip, zeta = 2 u_theta / r, which avoids the stiff Dirichlet stress condition (the Thom stability trap).

The m = 0 mode is advanced as the azimuthal-wind equation, the L1 operator plus the Reynolds-stress divergence -(1/r^2) d(r^2 <u_r u_theta>)/dr (flux form, zero at the walls), because psi_0 has no unique boundary value.

## 3. Heating

Heat flux H (W/m^2) enters the air at r = R: kappa_t dT'/dr = H/(rho c_p) (positive H gives warm air near the floor, unstable). At half an atmosphere rho c_p = 603, so H = 200 gives 0.33 K m/s.

- **Pattern:** three land strips (60 degrees each) alternate with three window strips, fixed in the hull frame. Land H_L(t) = 200 sin(2 pi t/24 h) by day (H = 100-350 range, Q10), -30 by night (-10 to -60); steady option H_L = 100. Windows are H = 0 by default, with an optional radiative loss h_w.
- **Stationarity:** the primary run is **flux-neutral**: the windows are the m = 3 heat sink, H_W = -H_L (equal strip areas), so the net wall flux is zero at all times. The uniform-sink alternative (windows 0, domain-mean relaxation with tau_T = 1 day) is a sensitivity run.
- **Inner wall:** insulating (dT'/dr = 0).

## 4. Numerics

- **Grid:** rfft in theta (N_theta = 256, |m| <= 85 after 2/3 dealiasing), second-order finite volume in r on a tanh-stretched grid (N_r = 64, first cell about 8 m, r_in = 0.1 R). A 1-2 m floor cell would make the radial CFL bind at about 0.8 s, so the floor layer is represented by the C_D drag instead.
- **Poisson:** per Fourier mode a tridiagonal solve of (d_rr + d_r/r - m^2/r^2) psi_m = -zeta_m, with psi_m = 0 at both walls; one batched Thomas over modes.
- **Time stepping:** IMEX SBDF2 (backward Euler first step). Diffusion and end-cap damping implicit (tridiagonal per mode), advection, baroclinic and Coriolis-buoyancy terms, drag and the lagged wall vorticity explicit. Advection is in advective form with centred radial differences and 2/3 dealiasing; the mean-wind Reynolds stress is in flux form so angular momentum telescopes exactly. Energy and enstrophy are conserved to discretisation order only.
- **Polar filter:** the high-m truncation tapers with r so the physical spacing r Delta theta stays at least 2 pi R / N_theta, otherwise the inner cell (4 m) sets the CFL.
- **CFL:** dt = 0.4 min(r Delta theta/|u_theta|, Delta r/|u_r|), about 6 s at 3 m/s.
- **Inner boundary:** free-slip insulating at r_in = 0.1 R. This is an artifact, since the real habitat has no inner wall (Zhong & Sun 2026 found it matters: strong zonal flow, Nu ~ Ra^0.1). Sensitivity runs at 0.05 R and 0.2 R are part of the deliverable.
- **Cost (numpy):** to be benchmarked on one step before runs are planned; the estimate is 6-10 ms/step at dt about 3 s, so 1.5-3 min per simulated day, possibly several times more. **Tests use 24 x 64 grids** (under 1 s each); long runs are a script and the CLI (`oneill atmosphere --run-l3`), never a test.

## 5. Validation tests

Each fails if its ingredient is wrong.

1. **Poisson vs analytic mode:** error falls 4x on halving h.
2. **Lamb-Oseen:** m = 0 vortex in a large domain; compare with Gamma exp(-r^2/4 nu t)/(4 pi nu t).
3. **Axisymmetric heating gives zero flow:** every m >= 1 mode stays at rounding level and T' equals the 1-D conduction solution.
4. **Linear drag:** end-cap term only: amplitude decays as exp(-t/tau_E); doubling L halves the rate.
5. **Energy:** inviscid, drag-free, no T', polar filter off: total energy conserved to discretisation order. **(a) Inviscid shear instability** exchanges energy between m = 0 and m >= 1 while the total stays constant (catches a sign error in the m = 0 Reynolds stress).
6. **(b) Axisymmetric spin-down in L3 matches L1.**
7. **(c) Angular-momentum budget:** the Reynolds torque integrates to zero; the Coriolis-buoyancy torque equals `-Omega dI'/dt`; `L_rel + Omega I'` changes only by wall and end-cap torques; and a steady state feels no Coriolis-buoyancy force at any radius.
8. **(d) Omega scaling:** changing Omega with beta Omega^2 and tau_E held fixed leaves results unchanged when Coriolis-buoyancy is off.
9. **(e) Torque sign:** the centre of a warm patch moves inward (mean u_r < 0 over the patch).
10. **(f) Heat budget:** d(int T')/dt = boundary flux - relaxation, exactly.
11. **Convection onset:** thin gap eta = 0.99, rigid walls, fixed-T walls, Ra from the mean radius, linearised single-mode growth rate at the critical m (about 2 pi r_m / (2.016 L)); target Ra_c = 1708 within 3 % (no E^(-4/3) law). Also rigid-rigid fixed-flux onset at Ra_c = 720 (Hurle, Jakeman & Pike 1967, wavenumber to 0, lowest mode), since production heats by flux with an insulating inner wall; and Ra_c rises as eta falls (cf. Zhong, Wang & Sun 2024).

## 6. Outputs answering H1

- **Eddy scale vs time:** peak and centroid of the azimuthal energy spectrum for r in [0.5 R, R], l = 2 pi r/m.
- **Wind:** rms u at 10 m, and gust = mean + 3 sigma (the 2D model has no turbulent gust).
- **Heat exchange:** Nu_eff = H R/(rho c_p kappa_t Delta T), Delta T from the rim layer to the inner wall.
- **Lock or condensate:** E(m = 3)/E_total and its phase against the strips (locked), versus E(m <= 2)/E_total growing with time and l tending to the domain size (inverse cascade).
- **Arrest wavelength:** derive eps from the forcing: the surface buoyancy flux over land is B_s = (g/T) H/(rho c_p), about 0.034 x 48/603 = 2.7e-3 m^2/s^3 for a daily-mean 48 W/m^2; over a depth of about 1 km with half the floor land, eps is about 1-3e-6 W/kg, and only part of it feeds the 2D cascade. Take mu = 1/tau_E (2.4e-5 /s at nu = 1) plus the domain-averaged floor drag (about 5e-6 /s), so mu is about 3e-5 /s. Then (0.08-0.18)(eps/mu^3)^(1/2) is about 0.5-1 km for eps = 1e-6 and 5-11 km for eps = 1e-4. I make no prediction: the run measures eps and l and decides whether l sits below or at the domain scale.
- **Mean wind:** the theta-mean wind profile ubar(r) and its sign are a headline output, with Coriolis-buoyancy ON and OFF.
- **Production set:** primary (flux-neutral windows, Coriolis-buoyancy ON, nu_t = 1); OFF; nu_t = 10; r_in = 0.05 R and 0.2 R; uniform sink. Run to statistical stationarity judged from the energy series; report wall time.

## 7. Effort and scope cuts

About 2-3 weeks: solver core and validation (1 week), heating and diagnostics (1 week), runs and write-up (0.5-1 week). Cuts: no stratification term; steady heating first and the day/night cycle second; one resolution plus a half-resolution check; no 3D, no moist physics, no adaptive stepping. Known weaknesses: constant nu_t, lambda about 0.03 (marginal 2D-ization), r_in artifact, and 2D by construction, so H1 can be supported but the breakdown of axis-invariance cannot be tested.
