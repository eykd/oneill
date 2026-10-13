# Weather in an O'Neill cylinder

Hurricane-force winds or gentle breezes? Breezes, mostly. The hurricane exists only while someone changes the spin.

That is our answer to worldbuilding.stackexchange question 159233. The habitat is the question's: a cylinder 3.2 km (2 miles) in radius and 20 km (12 miles) long, spinning at 0.53 rpm for 1 g at the floor (one turn every 113.5 s, floor speed 177 m/s), with three window strips between three land strips of low buildings and trees.

## TL;DR

Once the air turns with the hull, the spin makes no wind; sunlight through the windows drives the weather. It is convective, with thermals over the land strips; how fast they rise lies between the rotation-free scale (≈1.8 m/s, ≈15 min to cloud base) and the parcel bound (tenths of a m/s, 1.5–5 h). Clouds form about 0.8 km (half a mile) up at 70 % humidity and 1.8 km (1.1 miles) up at 50 % (O'Neill air at 290 K).

Our 2D simulation puts the wind 10 m (33 ft) up at about 0.9 m/s (2 mph), averaged over the floor and over day and night, with local highs (statistical, mean + 3σ around the ring) near 3.5 m/s (8 mph). That is the headline case, where the windows shed heat at a steady rate. Other assumptions about where the heat leaves give 1.0 to 1.7 m/s (2 to 4 mph), with local highs of 4 to 6.6 m/s (9 to 15 mph). Afternoons are breezier and nights nearly calm, with 25 to 50 times more wind energy by day than by night. The breeze fronts at the window and land edges are not resolved and may be windier. There is no prevailing wind.

The blender exists only while the spin changes. Spin a hull up from rest in one hour and the wind at 10 m reaches 18 to 23 m/s (40 to 50 mph), a gale at a bit over half hurricane force, and the air takes 7 to 23 hours to catch up. Those figures are order of magnitude only. A 24-hour spin-up cuts the wind to about a third, and a real habitat would spin up slowly or pressurise after spin-up, so the fast case is a worst case.

Rising air is pushed spinward and falling things land antispinward, so the question's "deflected antispinward" is half right. In the 2D model, eddies settle at 1 to 2 km (0.6 to 1.2 miles), near the strip spacing, and no habitat-scale vortex forms. The more of the heat that leaves the air away from the land strips, the more energy moves into the largest eddies, and if the air sheds heat everywhere, vortices the size of the habitat appear. Where the heat goes is a design choice that shapes the weather.

## Words we use

*Spinward* is the direction the floor moves, *antispinward* the opposite, and *axial* runs along the cylinder. *Up* means toward the axis. Air *co-rotates* when it turns with the hull, which makes it still relative to the floor. The *Coriolis force* is the sideways push that anything moving feels in a spinning frame. The *Rossby number* compares a flow's own speed to the spin: well below 1, rotation dominates; well above 1, it barely matters. The *canopy* is the buildings and trees treated as one rough layer that drags on the wind. The *end caps* are the cylinder's two flat ends.

## The question

The atmosphere the question describes (½ atm, 20 % O₂, 30 % N₂) is not O'Neill's. His 1974 *Physics Today* paper specifies "the air composition as well as the atmospheric pressure are those of sea level on Earth." The half-atmosphere comes from NASA SP-413, the 1975 Stanford torus study, at 50.8 kPa. We carry both as presets, "O'Neill 1974" and "SP-413", and every table reports both.

The question says anything moving radially is deflected antispinward, which holds only for falling things. Some answers put Coriolis at "1200×" or "thousands of times" Earth's. By rotation rate it is 759×, and one turn takes 113.5 s, not the half-minute one answer gives.

## Co-rotation is the default state

Still air feels centrifugal gravity toward the floor and a pressure gradient pushing back, and the two balance as described below. Coriolis acts only on moving air, so still air feels none. Solid-body rotation is also stable: angular momentum per unit mass rises outward, so a displaced ring of air returns to where it started (Rayleigh's criterion).

So in a sealed hull wind needs a driver, either uneven heating or a change in the spin rate. The rest of this document takes them in turn.

Warm air is lighter than cold air and so carries slightly less Coriolis force, which looks like a third driver. Mass is conserved in a sealed hull, so the net torque on every ring of air vanishes once temperatures stop changing. This *Coriolis-buoyancy* only shuffles momentum while the air warms or cools.

## Coriolis: what the residents feel

Walking spinward adds about 1.1 % to your weight for each m/s, and walking antispinward removes the same, so a brisk stroll changes your weight by under 2 %. Walking along the axis changes nothing, and nobody is pushed sideways on level ground.

Climbing brings the sideways push: spinward going up and antispinward going down, about 7.8 N (1.7 lbf) for a 70 kg person on a 1 m/s ladder. The push keeps its size with height while gravity weakens (to half, halfway to the axis), so near the axis "down" tilts strongly.

Dropped objects land antispinward of the point below: 0.53 m (21 in) from a 10 m drop, 17 m (56 ft) from a 100 m tower. A vertical jump lands slightly spinward. Wind obeys the same rules, so spinward wind presses harder on the ground than antispinward wind (the Eötvös effect).

<!-- gen:numbers -->
| Gas at floor | 2 omega (1/s) | vs Earth | inertial period (s) | Ro (U=3 m/s, L=1 km) | Ro (D=10 m) | N (1/s) | Rossby radius (m) | shedding f (Hz) | 2π f_shed / 2Ω |
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

A 3 m/s breeze across a kilometre has a Rossby number of 0.027, so rotation dominates anything kilometre-scale. A 10 m building sees 2.7 and mostly ordinary wind.

## Spin-up and spin-down: the "blender"

Changing the spin is the one way to get hurricane-class winds without weather, because momentum reaches the air only slowly, through the floor and the end caps. Our spin-up model follows the wind across the radius with turbulent (eddy) viscosity, canopy drag on the floor and friction at the end caps.

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

With end caps and a suburban canopy, a one-hour spin-up from rest puts the wind at 10 m near 18 to 23 m/s (40 to 50 mph), and the air takes 7 to 23 hours to reach 90 % of the hull's rotation, depending on the eddy viscosity. A 24-hour ramp drops that to 5 to 6 m/s (11 to 14 mph). The "peak wind anywhere" column is the lagging air a few hundred metres above the floor, just above the drag layer; residents feel the 10 m column, though towers and anything flying would meet the peak. The fast-ramp numbers are order of magnitude only, because the end-cap friction model is linear and holds only while the wind is small against the floor speed.

The end caps run the clock. Without them the same spin-ups take 38 hours to 15 days. With them, the Greenspan–Howard spin-up time is 11.8 h at an eddy viscosity of 1 m²/s and 3.7 h at 10 m²/s, and adjustment times scale with the inverse square root of the new spin rate.

We expected spin-down to be much faster than spin-up. It is not. When the hull slows, the air beside the floor outruns it, which is centrifugally unstable, while spin-up is stable (Kaiser et al. 2020). But the unstable layer is thin and the end caps set the clock, so a 10 % slowdown adjusts in 23.8 h and a 10 % speedup in 22.0 h.

Slow changes are gentle. Changing the spin by 1 % per day leaves a lag wind of 0.37 m/s (under 1 mph) at 10 m.

## Pressure, temperature and clouds

The habitat's whole sky is, thermodynamically, as deep as Earth's lowest 1.6 km (1 mile). The natural height coordinate is centrifugal geopotential, Φ = Ω²(R² − r²)/2, the energy per kilogram needed to lift something from the floor to radius r. It plays the role of g·z on Earth, so Earth's formulas for pressure, lapse rate and cloud base carry over exactly when written in Φ. Floor to axis is 15.69 kJ/kg.

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

Pressure at the axis is about 82 % of floor pressure for either gas. A well-mixed column cools at 9.76 K/km (5.4 °F per 1,000 ft) at the floor, Earth's dry rate, and the rate falls to zero at the axis along with gravity, for 15.6 K (28 °F) in all.

Cloud base is where true height and Earth-equivalent height part ways. For O'Neill air at 290 K, condensation starts 1.8 km (1.1 miles) up at 50 % humidity, more than halfway to the axis, and 0.8 km (half a mile) up at 70 %. O'Neill quoted "between 1100 and 1400 meters" using Earth lapse rates; converted to true height those become 1.4 to 2.1 km, matching the table.

## Convection and thermals

Can a thermal reach the clouds in a spinning sky? Yes, but more slowly than on Earth. With 200 W/m² of surface heating and a 1 km mixed layer, Deardorff's rotation-free velocity scale gives updrafts of 1.8 m/s in O'Neill air and 2.2 m/s in the thinner SP-413 air. The convective Rossby number is about 0.02 either way, so rotation shapes every thermal.

Our parcel model releases one warm parcel at the floor into still air, with buoyancy, Coriolis and an optional drag.

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

The drag-free rows are the extreme case. Rising pushes the parcel spinward, spinward motion makes it heavier, and it loops in a circle with a 57 s period, never rising more than about 5.5 m (18 ft) per kelvin of excess warmth. A free parcel drifts spinward about 11 m/s per 100 m of rise.

Both numbers are ceilings from a model with no pressure response. A real 100 m thermal has a Rossby number of 0.01, and rotation-driven waves stretch it along the axis within about a minute into a long column. For such a column the Coriolis force is a pure pressure gradient that the column's own pressure absorbs, so it gets no net spinward push. Real thermals fall between the parcel and the column.

Drag frees even the parcel. With drag k = 2Ω, the most favourable value for this model, a 3 K parcel climbs at 0.46 m/s and reaches cloud base (50 % humidity) in 99 minutes, and a 1 K parcel in 4.9 hours. The parcel is dry, so latent heat in the cloud is not counted.

Spinning-annulus experiments agree: convection with a hot outer wall turns quasi-2D once rotation is strong enough, and the habitat is past that point (Jiang et al. 2020).

## The 2D weather simulation

Within its assumptions, the 2D simulation finds eddies near the strip spacing, a kilometre or two across, with no habitat-scale vortices unless the air sheds its heat everywhere. Between those ends, the outcome slides with where heat leaves the air.

The model solves the flow in a cross-section of the habitat, assuming nothing varies along the axis, which strong rotation favours. In 2D, turbulence moves energy to larger scales until friction stops it, and the open question (hypothesis H1) was whether that happens near the strip scale or only at the size of the habitat. Theory says drag halts the cascade at some length (Boffetta 2007), but that length depends on an energy-input rate we could pin down only to within about two orders of magnitude, so the estimate ran from below the forcing scale to larger than the habitat. That spread is why we ran the model.

Gravity is centrifugal, so warm air is pushed toward the axis, and uneven heating around the ring spins up eddies. Land strips heat the air with a midday peak of 200 W/m² (daylight mean 127) and −30 W/m² at night. In the headline run the windows shed heat at a constant 54 W/m², the land's daily-mean gain, so the floor is heat-neutral over a day but not hour by hour. An instant balance, in which the windows warm the air at night, is kept as a sensitivity case, along with steady heating (±100 W/m²) and a volumetric sink that removes heat from all the air. Each validation test is built to fail if its ingredient is wrong (`docs/design/annulus2d.md` §5).

<!-- gen:annulus -->
| Case | days | nu_t | Coriolis-buoyancy | r_in / R | n_theta | rms wind 10 m (m/s) | mean wind 10 m (m/s) | spatial peak 10 m wind, mean\|u\|+3 sigma over theta (m/s) | peak m | m=3 fraction | m<=2 fraction | energy centroid length (km) | day/night KE ratio |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| day/night, windows always cool (headline) | 10 | 1 | on | 0.1 | 256 | 0.94 | -0.02 | 3.55 | 3 | 0.19 | 0.31 | 1.82 | 31x |
| day/night, windows balance land instantly | 10 | 1 | on | 0.1 | 256 | 1.06 | +0.01 | 4.05 | 3 | 0.19 | 0.17 | 1.40 | 42x |
| constant heating, no night | 10 | 1 | on | 0.1 | 256 | 1.37 | +0.01 | 5.26 | 3 | 0.18 | 0.11 | 1.12 | 2x |
| warm-air Coriolis term off | 10 | 1 | off | 0.1 | 256 | 1.05 | -0.01 | 4.05 | 3 | 0.19 | 0.20 | 1.49 | 39x |
| 10x eddy viscosity | 10 | 10 | on | 0.1 | 256 | 1.03 | +0.01 | 3.96 | 3 | 0.35 | 0.19 | 1.92 | 36x |
| inner wall at 0.05 R | 10 | 1 | on | 0.05 | 256 | 1.06 | +0.00 | 4.06 | 3 | 0.19 | 0.20 | 1.48 | 36x |
| inner wall at 0.2 R | 10 | 1 | on | 0.2 | 256 | 1.05 | +0.01 | 4.05 | 3 | 0.18 | 0.18 | 1.47 | 51x |
| heat lost from all the air | 10 | 1 | on | 0.1 | 256 | 1.74 | -0.00 | 6.64 | 2 | 0.11 | 0.65 | 4.28 | 94x |

Second-half (days 5-10) means; 2D (r, theta) model; heat in and out through the floor strips vs volumetric sink differ because the outcome depends on where heat leaves the air. The day/night ratio is the median over whole days after day 1 of the daily maximum over minimum kinetic energy.

Resolution check (same case, different azimuthal resolution; days 2-4 means):

| Run | n_theta | rms wind 10 m (m/s) | spatial peak 10 m wind, mean\|u\|+3 sigma over theta (m/s) | m=3 fraction | m<=2 fraction | energy centroid length (km) |
| --- | --- | --- | --- | --- | --- | --- |
| 4-day check, 256 around the ring | 256 | 1.05 | 4.03 | 0.18 | 0.18 | 1.41 |
| 4-day check, 512 around the ring | 512 | 1.03 | 3.95 | 0.24 | 0.17 | 1.48 |

*Two-dimensional (axis-invariant) Boussinesq convection with constant eddy viscosity, a free-slip inner wall at the stated radius (an artifact: the habitat has no inner wall) and a canopy-drag floor. It can support or undermine the 2-D hypothesis but cannot show where axis-invariance breaks down. Energy is conserved only to a few percent when eddies saturate.*
<!-- /gen -->

A 256 vs 512 azimuthal-resolution check over days 2-4, run on the instant-balance case, agrees within 7 % for winds, dissipation, eddy size and energy; the strip-locked fraction varies more (0.18 vs 0.24) and is quoted as a range.

Whenever heat leaves through the windows, the spectrum still peaks at the three-strip forcing scale and no habitat-scale vortex forms, but the share of energy in the largest patterns (those repeating once or twice around the ring) rises as heat leaves the air away from the land strips: 0.17 with instant balance, 0.31 with a constant window sink. The share locked to the three-strip pattern itself is 18 to 35 % across these runs, and the energy-weighted eddy size is 1.1 to 1.9 km (0.7 to 1.2 miles). Ten days of steady heating never grew the large patterns either, so the night is not what keeps eddies small.

When the air loses heat everywhere and the windows are inert (the volumetric sink), the largest patterns hold 0.65 of the energy, total energy is about 45 times higher, and the flow turns condensate-like, with vortices the size of the habitat. So H1 has no yes-or-no answer. It is a continuum set by the habitat's heat sink, and how a real habitat sheds its heat is a design choice. Little else matters: switching Coriolis-buoyancy off or moving the model's artificial inner wall (0.05 R to 0.2 R) changes nothing, and ten times the eddy viscosity raises the three-strip share to 35 % with unchanged winds.

With the constant window sink the 10 m wind is about 0.9 m/s (spatial rms, day and night mean; 2 mph), with local highs (statistical, mean + 3σ around the ring) near 3.5 m/s (8 mph). Kinetic energy is about three times the instant-balance case, and 25 to 50 times higher by day than by night; the extra energy sits in larger eddies higher up, not at street level. The other heat sinks give 1.0 to 1.7 m/s (2 to 4 mph) with local highs of 4 to 6.6 m/s (9 to 15 mph), highest for the volumetric sink.

At night the floor air cools to roughly 12 K (22 °F) below the interior. That cold pool's depth and strength depend on the constant eddy diffusivity used in a stable layer, so treat it as qualitative.

There is no prevailing wind: the day-and-night mean of the wind around the ring is under 0.03 m/s at 10 m and under 0.1 m/s anywhere. As the floor cools at night, mass shifts outward and a weak antispinward (retrograde) drift of up to about 0.4 m/s develops aloft, then reverses by day. Steady heating has no daily cycle. The windiest places are the strip edges, where breeze fronts sit fixed at every window/land boundary. The step-function heating sharpens them, and the model does not resolve their peak speed.

The model is 2D by construction, with canopy and end caps entering only as drag and end-cap friction as the main brake. It cannot show where the flow starts to vary along the axis, and its winds are area and time means, not street-level gusts.

## Winds near the floor

Spinward wind should be gustier near the rooftops than antispinward wind. Our surface model profiles the wind over suburb, orchard and forest canopies. It takes a mean wind of 3 m/s (7 mph) at 20 m as an input for studying gusts and stability; the 2D winds above are gentler on average.

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

At that wind, suburban gusts reach 5.5 m/s (12 mph), a gust factor of 1.83, typical of suburbs on Earth. The ½ atm air carries half the dynamic pressure, so wind loads halve with it.

Rotation changes which shear is turbulent, as in rotating channels (Brethouwer 2017). Spinward wind is destabilised where its shear exceeds 2Ω, antispinward wind is stable at every shear, and axial wind is neutral. Above the canopy the shear soon drops below 2Ω, and inside it the averaged profile is gentle too, so the unstable layer is a slab near canopy top: 4.6 to 13 m over the suburb, 2.6 to 12.3 m over the orchard and 8.7 to 24.3 m over the forest. This is an area-averaged result that leaves out the sharper shear around single buildings, and how much gustier spinward rooftops get has never been measured.

## Ekman spirals and vortex streets

Both of the question's sub-questions have answers: the floor gets a slight turbulent veer in place of an Ekman spiral, and spinward and antispinward winds should leave different wakes behind buildings.

There is no classical Ekman spiral on the floor, because the floor is parallel to the spin, like Earth's ground at the equator, so nothing turns the wind within the floor's plane. Axial wind does pick up a weak cross-flow near the floor in rotating-channel simulations, veering 2° to 11° (Oberlack et al. 2006), in a direction not yet pinned down for the habitat. Laminar theory would give the floor Stewartson layers (Stewartson 1957), but every layer here is turbulent.

The end caps, with the spin perpendicular to them, carry true Ekman layers with spirals. The molecular layers would be 1.6 cm thick, but they go unstable above about 0.1 m/s and are fully turbulent by roughly 0.15–0.5 m/s, so whenever there is weather the end-cap layers are turbulent and the eddy-viscosity thicknesses (4 to 13 m) apply.

Vortex streets need wind past an obstacle, so they appear only in sun-driven breezes. A 10 m building in a 3 m/s wind sheds a vortex every 17 to 28 s (Okajima 1982), a frequency about 3.4 times 2Ω/2π and comparable to the 57 s period of rotation-driven oscillations. Rotation steadies wakes only at Reynolds numbers far below a building's, so shedding persists, with tilted wakes (Machicoane et al. 2018). Flow over ridges in a rotating frame separates more on one side than the other (Savino & Wu 2024), so spinward and antispinward winds over rows of buildings should leave different wakes. That is untested for habitat geometry.

## Open questions

**H1: big eddies or locked cells.** Within 2D the answer is a continuum: eddies stay near the strip spacing when heat leaves through the floor, and grow to habitat scale as heat leaves the air everywhere. What happens once the flow varies along the axis is still open, and so is how a real habitat would shed its heat.

**H2: the axial-wind veer.** Its size comes from a low-Reynolds-number simulation, and its direction in the habitat is unverified.

**H3: the floor's boundary layer.** Theory for a sidewall that drives the flow is incomplete, and every layer here is turbulent. Whether floor drag or end-cap drag dominates the long-term momentum budget is open too; the floor has six times the end caps' area.

The turbulent spin-up time rests on our own derivation, and the spin-down eddy viscosity on a large extrapolation from simulations. We found no published dynamical study of a rotating-habitat atmosphere.

## Verdict: hurricane or breeze?

A day in the habitat runs like this. Nights are nearly calm, with cool air pooled on the floor. As the land strips warm, thermals rise toward the axis, leaning spinward and climbing more slowly than on Earth, and by afternoon there are clouds about a mile up and a breeze averaging about 2 mph, with local highs near 8 mph. It blows hardest at the strip edges where land meets window, and a little gustier over the rooftops when it blows spinward. The only hurricane-class wind comes from changing the spin, and even a one-hour spin-up from rest makes a gale for most of a day, short of hurricane force.

## Appendix: what about Navier–Stokes?

Readers ask whether the unsolved mathematics of fluid flow undermines models like these. As of October 2026, a forced Navier–Stokes blowup has been claimed: on 8 September 2026 OpenAI posted an unrefereed manuscript constructing a smooth force that drives a solution to infinite speed in finite time. The Clay Institute said "the Navier-Stokes problem has apparently been settled" but named no author and awarded nothing. The unforced problem, smooth start and no force, remains open. None of this touches habitat weather. Turbulence models describe averaged statistics on cells of metres to kilometres, far above this flow's 1 to 2 mm Kolmogorov scale, and never assumed the equations stay smooth. One early critic put the physics plainly: "A real fluid cavitates (water) or shocks (air) first" (Duraiswami).

## Reproducing these numbers

Every table here is generated by the package's models. Literature values, with confidence tags and verified URLs, are in `docs/research/atmosphere-literature.md`, and the 2D model's design and validation are in `docs/design/annulus2d.md`. From the repository root:

```sh
uv run oneill atmosphere --write-doc docs/atmosphere.md
```

This rewrites only the generated blocks (each fenced by a pair of HTML comment markers) and leaves the prose alone; a test fails if the committed file is stale.

The numbers in the Coriolis section come from:

```sh
uv run oneill --sections lateral vertical ballistics
```

The 2D cases run one at a time and merge into `docs/data/l3_summary.json`, which the `annulus` block reads:

```sh
uv run oneill atmosphere --run-l3 OUT_DIR --l3-case primary --l3-days DAYS
uv run oneill atmosphere --write-doc docs/atmosphere.md
```

The cases are `primary`, `cb_off`, `nu10`, `rin005`, `rin02`, `steady` and `uniform_sink`.

## References

Full citations, confidence tags and the URL audit are in `docs/research/atmosphere-literature.md`. The sources leaned on here:

- O'Neill (1974). The colonization of space. *Physics Today* 27(9):32–40. https://doi.org/10.1063/1.3128863 · https://nss.org/the-colonization-of-space-gerard-k-o-neill-physics-today-1974/
- NASA SP-413 (1977). *Space Settlements: A Design Study*, ch. 3–4. https://nss.org/settlement/nasa/75SummerStudy/Chapt4.html · https://ntrs.nasa.gov/citations/19770014162
- The question: https://worldbuilding.stackexchange.com/questions/159233/hurricane-force-winds-or-gentle-breezes-in-an-oneill-cylinder
- Greenspan & Howard (1963). On a time-dependent motion of a rotating fluid. *JFM*. https://doi.org/10.1017/S0022112063001415
- Stewartson (1957). On almost rigid rotations. *JFM* 3:17. https://doi.org/10.1017/S0022112057000452
- Kaiser et al. (2020). *JFM* 885 A6. https://doi.org/10.1017/jfm.2019.974
- Brethouwer (2017). *JFM* 828:424. https://arxiv.org/abs/1612.00254
- Oberlack, Cabot, Pettersson-Reif & Weller (2006). *JFM* 562:383. https://doi.org/10.1017/S0022112006001121
- Jiang, Zhu, Wang, Huisman & Sun (2020). *Sci. Adv.* 6:eabb8676. https://arxiv.org/abs/2010.07189
- Zhong, Wang & Sun (2024). *JFM* 992 A16. https://arxiv.org/abs/2402.00779
- Boffetta (2007). *JFM*. https://doi.org/10.1017/S0022112007008014 · https://arxiv.org/abs/nlin/0612035
- van Kan & Alexakis (2020). *JFM* 899 A33. https://arxiv.org/abs/1912.05394
- Deardorff (1970). *J. Atmos. Sci.* https://doi.org/10.1175/1520-0469(1970)027%3C1211:CVATSF%3E2.0.CO;2
- Okajima (1982). *JFM*. https://doi.org/10.1017/S0022112082003115
- Machicoane et al. (2018). *Phys. Rev. Fluids* 3:034801. https://arxiv.org/abs/1802.07602
- Savino & Wu (2024). *JFM* 999 A51. https://arxiv.org/abs/2403.13264
- Fefferman. The Navier–Stokes Millennium Prize problem statement (Clay). https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf
- Clay Mathematics Institute (11 Sep 2026). Navier-Stokes Announcement. https://www.claymath.org/news/navier-stokes-announcement/
- OpenAI (8 Sep 2026). Finite Time Blowup for Navier–Stokes. https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- Chen & Hou (2025). *PNAS* 122(27) (review). https://doi.org/10.1073/pnas.2500940122 · Part II: https://arxiv.org/abs/2305.05660
- Constantin, Ignatova & Vicol (2026). https://arxiv.org/abs/2609.20803 · Bastounis, Circelli & Hansen (2026). https://arxiv.org/abs/2610.08144 · Duraiswami (2026). https://arxiv.org/abs/2609.17642
