<!-- gen:gas -->
| Preset | p dry (kPa) | p vapour (kPa) | M (g/mol) | R_s (J/kg/K) | c_p (J/kg/K) | gamma | mu (1e-5 Pa s) | nu (1e-5 m2/s) | rho (kg/m3) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| O'Neill 1974 at 290.0 K | 101.3 | 0.0 | 28.97 | 287.0 | 1005 | 1.400 | 1.786 | 1.467 | 1.217 |
| SP-413 at 290.0 K | 49.7 | 1.0 | 29.96 | 277.5 | 978 | 1.396 | 1.834 | 2.970 | 0.618 |

*Dry mixture ideal gas; water vapour is a passive trace carried in the SP-413 total pressure but not in R_s, c_p or the density (about 2 % bias). Viscosity is Wilke mixing of per-species Sutherland laws.*
<!-- /gen -->

<!-- gen:equilibrium -->
| Preset | phi axis (kJ/kg) | p_axis/p_rim isothermal | p_axis/p_rim adiabatic | dp isothermal (kPa) | dp adiabatic (kPa) | dT adiabatic (K) | floor lapse (K/km) | air mass (1e11 kg) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| O'Neill 1974 at 290.0 K | 15.69 | 0.828 | 0.824 | 17.4 | 17.8 | 15.6 | 9.76 | 7.14 |
| SP-413 at 290.0 K | 15.69 | 0.823 | 0.818 | 8.8 | 9.0 | 16.0 | 10.03 | 3.61 |

| Gas at floor | RH | T at cloud base (degC) | cloud base, true height (m) | cloud base, Earth-equivalent (m) |
| --- | --- | --- | --- | --- |
| O'Neill 1974 at 290.0 K | 50% | 4.1 | 1828 | 1306 |
| O'Neill 1974 at 290.0 K | 70% | 10.1 | 787 | 690 |
| SP-413 at 290.0 K | 50% | 4.1 | 1750 | 1271 |
| SP-413 at 290.0 K | 70% | 10.1 | 763 | 672 |

*Still, horizontally uniform air in hydrostatic balance; the isothermal and dry-adiabatic profiles bracket reality. Cloud base uses Bolton's formula for air with a dry-adiabatic ascent; true height converts the geopotential, Earth-equivalent height is phi / g.*
<!-- /gen -->

<!-- gen:numbers -->
| Gas at floor | 2 omega (1/s) | vs Earth | inertial period (s) | Ro (U=3 m/s, L=1 km) | Ro (D=10 m) | N (1/s) | Rossby radius (m) | shedding f (Hz) | shedding / 2 omega |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| O'Neill 1974 at 290.0 K | 0.1107 | 759x | 56.7 | 0.027 | 2.71 | 0.0182 | 525 | 0.060 | 3.4 |
| SP-413 at 290.0 K | 0.1107 | 759x | 56.7 | 0.027 | 2.71 | 0.0184 | 532 | 0.060 | 3.4 |

| Viscosity | nu (m2/s) | Ekman number | Ekman depth (m) | E^1/4 H (m) | E^1/3 H (m) | spin-up time |
| --- | --- | --- | --- | --- | --- | --- |
| molecular, O'Neill 1974 at 290.0 K | 1.47e-5 | 2.65e-12 | 0.016 | 12.8 | 1.4 | 128 d |
| molecular, SP-413 at 290.0 K | 2.97e-5 | 5.36e-12 | 0.023 | 15.2 | 1.8 | 90 d |
| eddy, nu_t = 1 | 1.00e+0 | 1.81e-7 | 4.250 | 206.2 | 56.5 | 11.8 h |
| eddy, nu_t = 10 | 1.00e+1 | 1.81e-6 | 13.440 | 366.6 | 121.8 | 3.7 h |

*Viscous scales use the half-length L/2 = 10 km. Eddy viscosities are assumed constants, not derived. N is for an isothermal atmosphere at the floor.*
<!-- /gen -->

<!-- gen:surface -->
Mean wind 3 m/s at 20 m. Unstable layer: spinward wind, heights where dU/dz > 2 omega.

| Terrain | Gas at floor | z0 (m) | d (m) | h_c (m) | u* (m/s) | gust (m/s) | gust factor | q mean (Pa) | q gust (Pa) | unstable layer (m) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| suburban 2-3 storey | O'Neill 1974 at 290.0 K | 0.5 | 5.0 | 10.0 | 0.353 | 5.49 | 1.83 | 5.5 | 18.3 | 4.6 to 13.0 |
| suburban 2-3 storey | SP-413 at 290.0 K | 0.5 | 5.0 | 10.0 | 0.353 | 5.49 | 1.83 | 2.8 | 9.3 | 4.6 to 13.0 |
| orchard | O'Neill 1974 at 290.0 K | 0.6 | 4.0 | 6.0 | 0.365 | 5.58 | 1.86 | 5.5 | 18.9 | 2.6 to 12.3 |
| orchard | SP-413 at 290.0 K | 0.6 | 4.0 | 6.0 | 0.365 | 5.58 | 1.86 | 2.8 | 9.6 | 2.6 to 12.3 |
| forest | O'Neill 1974 at 290.0 K | 1.5 | 10.0 | 15.0 | 0.633 | 7.46 | 2.49 | 5.5 | 33.9 | 8.7 to 24.3 |
| forest | SP-413 at 290.0 K | 1.5 | 10.0 | 15.0 | 0.633 | 7.46 | 2.49 | 2.8 | 17.2 | 8.7 to 24.3 |

| Terrain | Wind toward | 1 m | 2 m | 5 m | 10 m | 15 m | 20 m | 30 m |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| suburban 2-3 storey | spinward | S | S | U | U | S | S | S |
| suburban 2-3 storey | antispinward | S | S | S | S | S | S | S |
| suburban 2-3 storey | axial | N | N | N | N | N | N | N |
| orchard | spinward | S | S | U | U | S | S | S |
| orchard | antispinward | S | S | S | S | S | S | S |
| orchard | axial | N | N | N | N | N | N | N |
| forest | spinward | S | S | S | U | U | U | S |
| forest | antispinward | S | S | S | S | S | S | S |
| forest | axial | N | N | N | N | N | N | N |

U unstable (shear turbulence augmented), S stable, N neutral.

*Area-averaged canopy profile; near-wall shear around individual buildings is not represented. The canopy attenuation follows from shear continuity and is low for the suburb (0.87, against 1-3 in urban fits). The 20 m reference height lies inside the forest's roughness sublayer, so the forest's u* and gusts are overestimates.*
<!-- /gen -->

<!-- gen:parcel -->
Warm parcel released at the floor in dry-adiabatic air (O'Neill air, 290 K). Short run 300 s; deflection read at 2 m of rise; cloud base at RH 50%, computed from a long run where shown.

| dT (K) | drag k (1/s) | max rise in 300 s (m) | peak speed (m/s) | deflection at 2 m rise (m) | drift at 2 m rise (m/s) | steady rise (m/s) | steady drift (m/s) | time to cloud base |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.000 | 5.5 | 0.61 | 0.92 | 0.222 | 0.000 | 0.00 | - |
| 1 | 0.010 | 10.9 | 0.53 | 0.91 | 0.212 | 0.027 | 0.30 | - |
| 1 | 0.050 | 35.6 | 0.35 | 0.91 | 0.176 | 0.115 | 0.25 | - |
| 1 | 0.111 | 45.5 | 0.23 | 0.87 | 0.125 | 0.153 | 0.15 | 4.9 h |
| 3 | 0.000 | 16.4 | 1.82 | 0.49 | 0.222 | 0.000 | 0.00 | - |
| 3 | 0.010 | 32.4 | 1.60 | 0.49 | 0.217 | 0.082 | 0.91 | - |
| 3 | 0.050 | 105.4 | 1.06 | 0.48 | 0.197 | 0.344 | 0.76 | - |
| 3 | 0.111 | 134.4 | 0.69 | 0.48 | 0.169 | 0.458 | 0.46 | 99 min |

*A point-mass parcel with no perturbation pressure or added mass: a ceiling for compact thermals, not a prediction. Real thermals are elongated along the axis by inertial waves and pressure absorbs the Coriolis force.*
<!-- /gen -->

<!-- gen:spinup -->
Floor condition for every row: canopy drag, C_D = 0.0178 (suburb, z0 = 0.5 m). Hull spun from rest to full rate; the air starts at rest.

| Hull ramp | nu_t (m2/s) | End caps | time to 90 % co-rotation | peak wind at 10 m (m/s) | peak wind anywhere (m/s) | peak epsilon | linear cap sink valid |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 60 min | 1 | on | 23.0 h | 18.3 | 159.3 | 1.26 | no (order of magnitude) |
| 60 min | 1 | off | 14 d | 18.6 | 165.2 | 1.28 | n/a (no end caps) |
| 60 min | 10 | on | 6.8 h | 23.1 | 131.5 | 1.02 | no (order of magnitude) |
| 60 min | 10 | off | 38.0 h | 24.1 | 147.2 | 1.06 | n/a (no end caps) |
| 24.0 h | 1 | on | 37.7 h | 5.3 | 66.3 | 0.78 | no (order of magnitude) |
| 24.0 h | 1 | off | 15 d | 6.8 | 133.6 | 1.03 | n/a (no end caps) |
| 24.0 h | 10 | on | 24.3 h | 6.2 | 24.5 | 0.44 | no (order of magnitude) |
| 24.0 h | 10 | off | 2 d | 8.9 | 74.5 | 1.00 | n/a (no end caps) |

Sudden change of the hull spin rate, with a Rayleigh-switched viscosity (larger where angular momentum falls outward). With end caps the sink rate scales with the new hull rate, so spin-down is the slower of the pair.

| Hull spin change | nu_t stable / unstable | End caps | time to 90 % adjustment | epsilon | linear cap sink valid |
| --- | --- | --- | --- | --- | --- |
| -10% | 1 / 10 | on | 23.8 h | 0.11 | yes |
| +10% | 1 / 10 | on | 22.0 h | 0.09 | yes |
| -10% | 1 / 10 | off | 15 d | 0.11 | n/a (no end caps) |
| +10% | 1 / 10 | off | 15 d | 0.09 | n/a (no end caps) |

Sensitivity of the time to 90 % adjustment (nu_t 1 / 10, end caps on) to the blend width of the viscosity closure, as a band over the widths tried (computed at a 150 s time step). The 50 % step has epsilon = 1, so its time is order-of-magnitude only.

| Hull spin change (end caps on) | time to 90 % adjustment, band over blend width 0.15-0.6 | epsilon | linear cap sink valid |
| --- | --- | --- | --- |
| -10% | 22.9-23.8 h | 0.11 | yes |
| +10% | 21.5-21.9 h | 0.09 | yes |
| -50% | 24.9-27.4 h | 1.00 | no (order of magnitude) |

Steady lag wind while the hull changes speed by 1 % per day (spin-up: the air lags, so the wind is negative).

| nu_t (m2/s) | End caps | peak lag wind (m/s) | lag wind at 10 m (m/s) |
| --- | --- | --- | --- |
| 1 | on | 0.69 | -0.365 |
| 1 | off | 10.58 | -0.987 |
| 10 | on | 0.23 | -0.233 |
| 10 | off | 1.61 | -0.963 |

Molecular viscosity only. Illustrative: the laminar end-cap sink is invalid at these Ekman Reynolds numbers (turbulent layers), and the 1 h ramp also exceeds epsilon = 0.3.

| Gas at floor | molecular nu (m2/s) | time to 90 % co-rotation (1 h ramp) | peak wind at 10 m (m/s) | Ekman Reynolds number |
| --- | --- | --- | --- | --- |
| O'Neill 1974 at 290.0 K | 1.47e-5 | 292 d | 10.7 | 11,921 |
| SP-413 at 290.0 K | 2.97e-5 | 205 d | 10.8 | 8,401 |

*One-dimensional axisymmetric model with constant eddy viscosities and a canopy-drag floor (C_D from the suburban log law at 10 m); winds are quoted at that 10 m reference height, the one the drag law is written for. The end-cap sink is the linear Greenspan-Howard result, valid for epsilon = |u| / (Omega r) <= 0.3; rows with a larger peak epsilon (fast ramps, the 50 % step) are order-of-magnitude only, and turbulent end caps are not linear. Adjustment times depend on the blend width of the viscosity closure; the sensitivity band is shown, and 50 % times drift about 3 % with the production time step.*
<!-- /gen -->
