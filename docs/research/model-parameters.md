# Atmosphere model parameters

Companion to `atmosphere-literature.md`. Sources there carry full URLs; short keys here.
Frame: x spinward, y axial, z up (toward axis), right-handed ⇒ **Ω⃗ = −Ω ŷ** (ŷ points along −Ω̂; see the convention note in the literature file). Ω = 0.0554 s⁻¹, R = 3.2 km,
L = 20 km.

Modules: gas · L0 (`equilibrium`/`numbers`) · L1 (`spinup`) · L2 (`parcel`) · L2b (`surface`) ·
L3 (`annulus2d`).

## Gas presets

| Parameter | Value / range | Units | Source | Conf. | Module |
|---|---|---|---|---|---|
| p_rim, `oneill_1974` | 101.3 (sea-level Earth air) | kPa | O'Neill 1974 | HIGH | gas |
| p_rim, `sp413_torus` | 50.8 dry + 1.0 H₂O | kPa | SP-413 Tab. 4-3 | HIGH | gas |
| pO₂ / pN₂ / pCO₂ (SP-413) | 22.7 / 26.6 / < 0.4 | kPa | SP-413 Tab. 4-3 | HIGH | gas |
| T (SP-413) | 20 ± 5 (criteria: 23 ± 8) | °C | SP-413 Tab. 4-3, 3-3 | HIGH | gas, L0 |
| RH (SP-413) | 50 ± 10 (Ch. 3 text: ~40) | % | SP-413 | HIGH | gas, L0 |
| ν at ½ atm, 288 K | ≈ 3.0×10⁻⁵ | m²/s | Sutherland; 2× sea level | HIGH | gas |
| ρc_p at ½ atm | ≈ 603 | J m⁻³ K⁻¹ | derived | HIGH | gas, L2b |

## L0 checks

| Parameter | Value / range | Units | Source | Conf. | Module |
|---|---|---|---|---|---|
| Ω²R | 9.82 (1 g needs Ω = 0.05536) | m/s² | derived | HIGH | L0 |
| Φ_rim→axis = Ω²R²/2 | 15.7 | kJ/kg | derived | HIGH | L0 |
| p_axis/p_rim | 0.827 isothermal (288 K); 0.822 adiabatic; O'Neill ≈ 0.825 | — | O'Neill 1974; derived | HIGH | L0 |
| Adiabatic ΔT rim→axis | 15.6 | K | derived | HIGH | L0 |
| Dry lapse rate | Ω²r/c_p: 9.8 at rim → 0 at axis | K/km | derived | HIGH | L0 |
| Cloud base (O'Neill) | 1.1–1.4 **Earth-equivalent (Φ/g)** ⇒ ≈ 1.4–2.1 true habitat height | km | O'Neill 1974; converted via Φ | HIGH | L0 (convert via Φ before testing) |
| Coriolis ratio vs Earth | 2Ω/2Ω_E ≈ 760 | — | derived | HIGH | L0 |
| Inertial period π/Ω | 56.7 | s | textbook | HIGH | L2 |

## L0 numbers, L1 spin-up

| Parameter | Value / range | Units | Source | Conf. | Module |
|---|---|---|---|---|---|
| Ekman depth δ_E = √(ν/Ω) | 0.023 (molecular ½ atm); 4.2 (ν_t = 1); 13 (ν_t = 10) | m | GH63; derived | HIGH | L0 |
| GH length L in E = ν/(ΩL²) | **half-length = 10 km** | m | GH63 | HIGH | L0, L1 |
| τ_E = L/√(νΩ) | 90 d (molecular); 37 h (ν = 0.1); 11.8 h (ν = 1); 3.7 h (ν = 10) | — | GH63; derived | HIGH form / LOW ν_t | L1 |
| Sidewall reach in τ_E | E^{1/4}L ≈ 15 m (mol.), 200 m (ν_t = 1) ≪ R ⇒ end caps control | m | GH63 §5b | HIGH | L1 |
| Turbulent end-cap spin-down | τ ~ L/(C_D ΔU), algebraic decay; ≈ 4–12 d for ΔU = 10 m/s | — | derivation (Eliassen–Lystad form, unverified) | MED | L1 |
| Rim-drag spin-down | τ ~ R/(2 C_D U) ≈ 1–2 d | — | derivation | LOW | L1 |
| C_D (turbulent wall) | 1–2.5×10⁻³ (smooth/geostrophic); log law [κ/ln(z/z₀)]² for rough | — | Rossby similarity | MED | L1 |
| Ekman instability Re_δ (δ = √(ν/Ω)) | type II 55; type I 110–118 (theory), 125 ± 5 (expt) | — | Lilly 66; Faller 63; Tatro & M-C 67 | HIGH / MED | L0, L1 |
| Ekman fully turbulent Re_δ | ~150 (usual); ~500 (Buffett 2021) | — | secondary; Buffett | MED | L1 |
| Stewartson thicknesses | H E^{1/4}, H E^{1/3} | m | Stewartson 57 | HIGH | L0 |
| Spin-down centrifugal onset | τ_c = νt/R² = 9.4 Re^{−4/3} | — | Kaiser 2020 | HIGH (DNS range) | L1 |
| Spin-down ν_t enhancement | ν_t/ν ≈ 0.01 Re^{2/3} (Re = ΩR²/ν) | — | derived from Kaiser 2020 | LOW (extrapolated 6–7 decades) | L1 |
| Rayleigh switch | unstable where d(r²v_θ)²/dr < 0 | — | Rayleigh 1917 | HIGH | L1 |

## L2b surface layer

| Parameter | Value / range | Units | Source | Conf. | Module |
|---|---|---|---|---|---|
| κ | 0.40 | — | standard | HIGH | L2b |
| z₀, suburb (h ≈ 10 m) | 0.3–0.8 (class 6–7: 0.5–1.0) | m | Wieringa; G&O 99 | MED | L2b |
| d, suburb | 3–7 (0.5–0.7h) | m | G&O 99; Macdonald 98 | MED | L2b |
| z₀, forest/orchard | 0.1h (0.05–0.13h): 1–2 forest, 0.5–1 orchard | m | Raupach 94 | HIGH rule | L2b |
| d, forest/orchard | 0.6–0.8h (2/3h) | m | Raupach 94 | HIGH | L2b |
| Canopy element C_d | 0.2 (0.1–0.3) | — | Katul 04 impl.; Shaw–Schumann | MED | L2b, L3 |
| LAI | deciduous 3–6, conifer 4–9, orchard 1.5–4 | — | standard | MED | L2b |
| Drag length 1/(C_d a) | median 21 (trees), 0.7 (crops) | m | arXiv 2411.01570 | HIGH | L2b, L3 |
| Roughness sublayer depth | 2–3h | m | Simpson 98 via secondary | MED | L2b |
| u*/U (10–15 m above d) | 0.12 suburb, 0.17 forest, 0.07 grass | — | log law | HIGH formula | L2b |
| σ_u/u* | 2.2–2.5 | — | standard | MED | L2b |
| Gust factor (3 s / 10 min) | 1.4–1.5 open; 1.7–2.0 suburb/forest | — | Wieringa 73 (unread); peak-factor estimate | MED | L2b |
| Rotational stability | Ri_rot = S(S+1), S = −2Ω/(dU/dz); unstable −1 < S < 0 | — | Brethouwer 2017 | HIGH | L2b |
| Spinward unstable layer | dU/dz > 2Ω = 0.111 s⁻¹ ⇒ z ≲ u*/0.044 (2–7 m) | m | derived | MED | L2b |
| Antispinward | S > 0, always stable | — | Brethouwer 2017 | HIGH | L2b |
| Axial wind | Bradshaw-neutral; turbulent veer 2–11°, W/U_b 3–10 % | ° | Oberlack 06 DNS (sign unverified) | MED / LOW sign | L2b (H2) |
| Channel relaminarization | Ro_b = 2Ωh/U_b ≈ 3.0 | — | Brethouwer 2017; Wallin 2013 | HIGH | L2b |
| St (building / round) | 0.12–0.14 / 0.2 | — | Okajima 82 | HIGH | L2b |
| Rotating-wake stabilization (Ω ⟂ axis) | Re_c ≈ (275 ± 25)/Ro, Ro = U/(2Ωd) | — | Machicoane 2018 | HIGH | L2b |

## L2 parcels

| Parameter | Value / range | Units | Source | Conf. | Module |
|---|---|---|---|---|---|
| Spinward drift per rise | 2ΩΔr = 11.1 per 100 m | m/s | derived | HIGH | L2 |
| Stratified oscillation | √(N² + 4Ω²) | s⁻¹ | derived | HIGH | L2 |
| Daytime sensible heat H | 100–300 veg, 150–350 suburb | W/m² | G&O 2002 | MED | L2, L3 |
| Night H | −10 to −60 (rural) | W/m² | Oke 1982 | MED | L2, L3 |
| Bowen ratio | 0.1–0.3 irrigated … 0.7–2 suburb | — | standard | MED | L2, L3 |
| w* (H = 200, z_i = 1 km, ½ atm) | ≈ 2.2 | m/s | Deardorff 70; derived | HIGH | L2 |
| Convective Rossby w*/(2Ωz_i) | ≈ 0.02 | — | derived | HIGH | L2 |

## L3 annulus

| Parameter | Value / range | Units | Source | Conf. | Module |
|---|---|---|---|---|---|
| Ra (ACRBC) | ½ω²(R_o+R_i)βΔL³/(νκ) | — | Jiang 2020 | HIGH | L3 |
| Ro⁻¹ | 2(βΔ r_m/L)^{−1/2}; habitat ≈ 19 (L = 1 km) to 34 (L = 3.2 km) for ΔT = 1 K, scaling as ΔT^{−1/2} | — | Jiang 2020; derived | HIGH / MED | L3 |
| Onset, flat end caps | non-rotating RB, Ra_c → 1708 (η → 1); n_c = 5 at η = 0.5 | — | Zhong/Wang/Sun 2024; Westerburg & Busse 2003 | HIGH | L3 test |
| Onset with sloping caps (NOT the habitat) | Ra_c ≈ 3η_P^{4/3} ∝ E^{−4/3} | — | Busse 1970 | HIGH | L3 (counter-test) |
| Quasi-2D threshold | Ro⁻¹ ≳ 10–15 | — | Jiang 2020 | HIGH | L3 |
| Roll aspect ratio | ≈ 1 gap; 4 pairs at η = 0.5 | — | Jiang 2020 | HIGH | L3 |
| Turbulent drift | prograde; ω_d = −0.2η + 0.15 (η 0.3–0.7), → 0 as η → 1 | rad per t_ff | Wang 2022 | HIGH | L3 |
| Nu scaling | Ra^0.27 (< 10¹⁰); Ra^0.40 (10¹¹–10¹²) | — | Jiang 2020, 2022 | HIGH | L3 |
| Inverse-cascade C | 6 (4–6) | — | Boffetta 2007 | HIGH | L3 |
| Arrest wavelength | (0.08–0.18)(ε/μ³)^{1/2} | m | Grianik 2004 | HIGH | L3 |
| Quadratic-drag arrest | k_a ≈ 51 C_d (C_d = c_d/h) | m⁻¹ | Grianik 2004 | HIGH | L3 |
| 2D-ization control | λ = ℓ_in^{5/3}Ω/(Hε^{1/3}) > 0.03; H = 20 km | — | van Kan & Alexakis 2020 | HIGH | L3 (H1) |
| Canopy area fraction in (r,θ) | ≈ 2h_c/R ≈ 0.0125 | — | derived | MED | L3 |
| End-cap Ekman drag | μ_E ≈ 2√(ν_tΩ)/L (laminar form) or C_D U/(L/2) | s⁻¹ | standard | LOW / MED | L3 |
| Compressional β | ≈ 0.19 density scale heights rim→axis | — | derived | MED (strength contested) | L3 |
