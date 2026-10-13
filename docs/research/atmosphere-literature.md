# O'Neill cylinder atmosphere: literature findings

Research pass for Part 2 of the simulator plan (worldbuilding.stackexchange Q159233).
Compiled 2026-10-12 from an Extensive-mode Research run (multiple explorer agents plus two
verifier agents) and arXiv lookups. Every URL below was fetched or batch-checked on
2026-10-12. Where a publisher (AMS, Wiley, T&F, PNAS, SIAM, Annual Reviews) blocks scripted
fetches with HTTP 403, the DOI was confirmed registered through Crossref or doi.org instead;
those are marked †.
URL audit 2026-10-12: 116 unique URLs were batch-checked with curl. 95 returned 200. All 20
publisher-blocked DOIs (403, or DOIs containing "()<>" that the checker mis-split) returned
responseCode 1 (registered) from the doi.org handle API. The WB.SE page returns 403 to scripts;
its content was read through the Stack Exchange API instead.

Confidence tags: **[HIGH]** read in a primary source or a standard textbook result checked by
derivation; **[MED]** single source, secondary source, or our own derivation on sourced inputs;
**[LOW]** extrapolation, assumption, or unverified.

Habitat used throughout: R = 3.2 km, L = 20 km (half-length H = 10 km), Ω = 0.0554 s⁻¹
(Ω²R = 9.82 m/s², 0.529 rpm, period 113.4 s), 2Ω = 0.111 s⁻¹, inertial period π/Ω = 56.7 s.

---

## Sign convention (used in every section)

Local floor frame: **x = spinward, y = axial, z = up (toward the axis), right-handed.** The
floor point sits at r = −R ẑ from the axis. Requiring the floor velocity Ω⃗ × r⃗ to point along
+x gives **Ω⃗ = −Ω ŷ**. In this frame the Coriolis acceleration is

  −2Ω⃗ × u⃗ = 2Ω (w x̂ − u ẑ).

- Rising air (w > 0) is pushed **spinward**.
- Spinward wind (u > 0) is pushed **toward the floor** (Eötvös effect: it is effectively heavier).
- Axial wind (v) feels **no Coriolis force** (it is parallel to Ω).

**"Axial" is ambiguous, so pin it.** In the right-handed frame above, ŷ points along **−Ω̂**, which
is why Ω⃗ = −Ωŷ. If code instead defines axial = +Ω̂ (natural with `Z_HAT` as the spin axis),
then the plan's `LocalVector(spinward, axial, up)` is a *left-handed* triple, because
(spinward, up, +Ω̂) is the right-handed one. Either convention works, but `cylinder.py`'s local
basis and the Coriolis-sign tests must state which one they use. [HIGH, derivation;
independently reproduced by three agents]

---

## Headline findings that change the plan

1. **O'Neill specified 1 atm, not ½ atm.** The ½ atm atmosphere is NASA SP-413's, and SP-413
   designs the Stanford *torus*. `gas.py` should carry both baselines. (Q1)
2. **The Bradshaw instability is height-dependent.** Spinward wind is destabilized only where
   the shear dU/dz exceeds 2Ω (the bottom ~1–7 m). Above that, spinward wind is
   rotation-*stabilized*. Antispinward wind is stabilized at every shear. The plan's
   `B = S(S−1)` paired with `S = −2Ω/(dU/dy)` is sign-inconsistent. (Q5, claim c)
3. **Convection onset for L3 is the non-rotating value.** With flat end walls, Coriolis on
   axis-parallel rolls goes entirely into pressure. Onset is Rayleigh–Bénard-like (Ra_c → 1708
   as the radius ratio → 1). It is *not* Ra_c ∝ E^(−4/3), which needs sloping end caps. (Q7)
4. **Canopy drag is not a domain-wide Rayleigh drag in the (r,θ) plane.** The canopy is a thin
   annulus on the sidewall, covering ~1% of the disc. The true 2D-drag analog is end-cap Ekman
   friction. The arrest-scale formula also needs a prefactor of 0.08–0.18. (Q8, claim d)
5. **The Navier–Stokes paragraph is out of date.** OpenAI claimed (8 Sep 2026) a forced-blowup
   proof aimed at Clay alternatives (C)/(D). Clay calls the problem "apparently settled" but has
   awarded nothing. Chen–Hou's journal paper is PNAS **2025**. (Q12, claim h)

---

## Q1. O'Neill's specified atmosphere

- **O'Neill 1974 (Physics Today 27(9):32–40)** specifies sea-level Earth air: "the air
  composition as well as the atmospheric pressure are those of sea level on Earth." [HIGH]
  - Pressure profile p = p₀ exp(−a(R²−r²)), giving p_axis/p_rim ≈ 0.825. [HIGH]
  - Cloud base: "between 1100 and 1400 meters". He used constant-g Earth lapse rates (3.1 K per
    300 m dry), so these are *Earth-equivalent* heights (Φ/g). Converting through
    Φ = Ω²(2Rh − h²)/2 puts the true habitat height at **≈ 1.4–2.1 km**. [HIGH quote; HIGH
    conversion arithmetic]
  - The plan's "about 1 km at 50 % RH" looks low even in Earth-equivalent terms. The LCL is
    ≈ 125 m × (T − T_d); at 50 % RH and 15–20 °C that is ≈ 1.3 km. L0 should compute this, not
    take it as an input. [MED]
  - Rotation: "The rotation period is two minutes", i.e. 0.5 rpm, which gives only 0.895 g at
    3200 m. Exactly 1 g needs 0.5286 rpm (Ω = 0.05536 s⁻¹). [HIGH]
  - On thermals for soaring: "sufficient atmospheric instability to provide them with lift." He
    gives no wind or Coriolis analysis. [HIGH]
  - https://doi.org/10.1063/1.3128863 · full text: https://nss.org/the-colonization-of-space-gerard-k-o-neill-physics-today-1974/
- **NASA SP-413 (1977), Table 4-3, "Habitat Atmosphere"** (Stanford torus) [HIGH]:

  | Component | Value |
  |---|---|
  | O₂ | 22.7 kPa |
  | N₂ | 26.6 kPa |
  | CO₂ | < 0.4 kPa |
  | H₂O | 1.0 kPa |
  | **Total** | **50.8 kPa** (380 mmHg; the oft-quoted 51.5 kPa is wrong) |
  | Temperature | 20 ± 5 °C |
  | Relative humidity | 50 ± 10 % |

  - Table 3-3 (physiological criteria): pO₂ 22.7 ± 9 kPa, T 23 ± 8 °C, inert gas 26.7 kPa,
    pH₂O 1.00 ± 0.33 kPa, rotation ≤ 1 rpm.
  - SP-413 wants "weather" as a desirable unpredictable element but contains no atmospheric
    dynamics. It treats Coriolis only as a vestibular issue.
  - https://nss.org/settlement/nasa/75SummerStudy/Chapt4.html · https://nss.org/settlement/nasa/75SummerStudy/Chapt3.html · https://ntrs.nasa.gov/citations/19770014162
- **The question's own premises** (read via the SE API; the page returns 403 to scripts):
  - It attributes "½ atm, 20 % O₂, 30 % N₂" to O'Neill. This is wrong.
  - It says radial movers are "deflected anti-spinward". This is half wrong: inward (rising)
    motion is deflected *spinward*.
  - Answers quoting "1200×" or "thousands of times" Earth's Coriolis overstate it. By rotation
    rate the ratio is ≈ 760×.
  - One answer says the rotation is "half a minute"; it is 113 s.
  - https://worldbuilding.stackexchange.com/questions/159233/hurricane-force-winds-or-gentle-breezes-in-an-oneill-cylinder
- **Model parameters:** two named presets, `oneill_1974` (101.3 kPa, Earth air) and
  `sp413_torus` (50.8 kPa dry + 1.0 kPa H₂O, x_O₂ ≈ 0.45 of dry gas), with T = 288–293 K.
  [HIGH]

## Q2. Greenspan–Howard spin-up with turbulent Ekman layers

- **Laminar result [HIGH].** T = (L²/νΩ)^{1/2} = E^{−1/2}/Ω, with E = ν/(ΩL²).
  - Greenspan & Howard (GH) define L as the **half-gap**: plates 2L apart. Here L = 10 km.
  - Oruba et al. build E on the full height and quote it only to order of magnitude, which hides
    a factor of 2.
  - GH give three phases: the Ekman layer forms in O(Ω⁻¹), the interior spins up in
    E^{−1/2}/Ω, and residual oscillations decay viscously in L²/ν.
  - The DOI ends in …1415; the widely circulated …1414 is a 404.
  - https://doi.org/10.1017/S0022112063001415 · text: https://empslocal.ex.ac.uk/people/staff/gv219/classics.d/Greenspan_Howard63.pdf
- **Habitat numbers** (our calculation, L = 10 km):

  | ν (m²/s) | τ_E = L/√(νΩ) |
  |---|---|
  | 3×10⁻⁵ (molecular, ½ atm) | ≈ 90 days |
  | 0.1 | 37 h |
  | 1 | 11.8 h |
  | 10 | 3.7 h |

  [HIGH formula; ν_t is LOW]
- **Sidewall versus end-cap control (GH §5b) [HIGH].** Within one spin-up time the sidewall
  influences a layer E^{1/4}L thick. Here E^{1/4}L ≈ 15 m (molecular) and ≈ 200 m (ν_t = 1),
  both ≪ R. So the *laminar* GH clock is set by the end caps.
- **Turbulent end caps [MED, derivation; no peer-reviewed turbulent spin-up study found].**
  - Replace Ekman pumping with a quadratic stress C_D|ΔU|ΔU. Then L·dΔU/dt ≈ −C_D|ΔU|ΔU.
  - Result: **τ_turb ~ L/(C_D ΔU)**, with algebraic rather than exponential decay. To leading
    order it is independent of Ω. This is the Eliassen–Lystad hurricane spin-down form, but that
    attribution is unverified.
  - With C_D ≈ 1–2.5×10⁻³ and ΔU = 10 m/s, τ ≈ 4–12 days via the end caps.
- **Rim drag may dominate [LOW].**
  - The rim area (2πRL ≈ 4×10⁸ m²) is about 6× the end-cap area (6.4×10⁷ m²).
  - A turbulent rough floor damps azimuthal relative wind on τ ~ R/(2C_D U), about 1–2 days.
  - Laminar theory treats the sidewall as passive, so this needs testing in L1.
- **Nonlinear spin-up from rest [MED].** Wedemeyer 1964, JFM 20:383: a sidewall shear front
  propagates inward; the leading scaling stays E^{−1/2}/Ω.
  https://doi.org/10.1017/S002211206400129X
- **Ekman instability [HIGH/MED].** With Re_δ = Uδ/ν and δ = (ν/Ω)^{1/2}:
  - Type II onset ≈ 55 (Lilly 1966 theory; 56 ± 2 measured by Tatro & Mollo-Christensen 1967).
  - Type I onset ≈ 110–118 in theory and 125 ± 5 in experiment (Faller 1963).
  - Fully turbulent ≈ 150, the usual citation [MED]. Buffett 2021 finds ≈ 500 for a neutral
    layer.
  - At molecular ν, δ = 2.3 cm, so Re = 125 corresponds to U ≈ 0.16 m/s. **Any real wind makes
    the end-cap layers turbulent.**
  - https://doi.org/10.1017/S0022112063000458 · https://doi.org/10.1175/1520-0469(1966)023%3C0481:OTIOEB%3E2.0.CO;2 † · https://doi.org/10.1017/S0022112067002289 · https://par.nsf.gov/biblio/10218706
- Benton & Clark 1974, "Spin-up" (review): https://doi.org/10.1146/annurev.fl.06.010174.001353 †
  (content not read).

## Q3. Stewartson layers when the sidewall drives (H3), and claim (b)

- **Claim (b) is CONFIRMED [HIGH].**
  - Ekman layer thickness scales with the *wall-normal* component of Ω, so it is undefined on a
    wall parallel to Ω.
  - On the floor, axial wind feels no Coriolis force, and spinward/antispinward wind feels only
    a wall-normal one (absorbed by pressure). Nothing turns the flow within the wall plane, so
    **there is no Ekman spiral on the floor**.
  - The end caps (Ω normal) carry true Ekman layers.
  - The straight sidewall carries the Stewartson sandwich:
    - δ_{1/4} ~ H E^{1/4} = H^{1/2}(ν/Ω)^{1/4}, carrying the azimuthal jump;
    - δ_{1/3} ~ H E^{1/3}, ageostrophic, carrying the axial return flow.
  - Stewartson 1957, JFM 3:17: https://doi.org/10.1017/S0022112057000452
- **Conflict resolved.** One agent proposed the "equatorial Ekman layer" (E^{2/5} thick, E^{1/5}
  axial extent). That structure belongs to a *curved spherical* wall passing tangent to Ω. The
  habitat floor is a straight cylinder parallel to Ω, so the Stewartson layers apply.
  Marcotte, Dormy & Soward 2016: https://arxiv.org/abs/1602.08647 [HIGH]
- **Sidewall-driven flow [MED].**
  - van Heijst 1983 (JFM 130:1) treats a differentially rotating sidewall: the end-wall meridional
    flow changes the roles of the E^{1/4} and E^{1/3} layers.
    https://doi.org/10.1017/S0022112083000932 (content via secondary source)
  - Oruba, Soward & Dormy 2017 (stress-free sidewall): fast E^{−1/2} spin-down via the rigid
    base, then slow O(R²/ν) decay through the E^{1/3} layers.
    https://arxiv.org/abs/1611.06511
  - With stress-free or co-rotating end caps, sidewall-only forcing reduces to lateral
    diffusion, τ ~ R²/ν_t. For ν_t = 1 that is ≈ 118 days. No paper treats this exact case.
- **Practical caveat for H3 [MED].** At habitat Reynolds numbers (molecular E ≈ 5×10⁻¹²), every
  layer is turbulent. Stewartson asymptotics give structure and scalings, not magnitudes. Model
  the floor as a turbulent rough-wall layer plus the rotational-stability term from Q5.

## Q4. Eddy viscosity during unstable spin-down

- **The asymmetry is CONFIRMED [HIGH].** By Rayleigh's criterion, centrifugal instability
  appears only when the outer wall decelerates relative to the fluid (spin-down). Spin-up is
  Rayleigh-stable. Kaiser et al. 2020, JFM 885 A6 (DNS, impulsive stop, Re = ΩR²/ν = 3×10³ to
  2.8×10⁴): https://doi.org/10.1017/jfm.2019.974
  - Onset at τ_c = νt_c/R² = 9.4 Re^{−4/3}.
  - Taylor–Görtler rolls saturate at Ωt ≈ 4.3 and break down at Ωt ≈ 5.3.
  - In the turbulent stage the boundary layer still grows as (νt)^{1/2}, with
    a_turb/a_lam = 0.1 Re^{1/3} (≈ 1.4–3.0 over their range).
  - Other refs: Neitzel & Davis 1981, https://doi.org/10.1017/S002211208100267X ; Mathis &
    Neitzel 1985, https://doi.org/10.1063/1.865108 †
- **Effective ν_t [LOW].** ν_t/ν ≈ (a_turb/a_lam)² ≈ 0.01 Re^{2/3}, our derivation from Kaiser's
  growth law. The habitat's ΩR²/ν ≈ 2×10¹¹ is 6–7 orders beyond the DNS, so treat this as a
  bracket only.
- **Recommended L1 closure [MED].**
  - Use a turbulent ν_t in both directions, since a rough floor is turbulent regardless.
  - Add a Rayleigh-switched enhancement when Φ = r⁻³ d(r²v_θ)²/dr < 0 (spin-down).
  - Do not make spin-up laminar: outer-wall-driven Couette flow shows subcritical transition.

## Q5. Spanwise-rotation (Bradshaw) thresholds, and claim (c)

- **Criterion, in Brethouwer's convention [HIGH].** S = −2Ω/(dU/dy).
  - S > 0 is cyclonic and damped; S < −1 is also damped.
  - −1 < S < 0 is anticyclonic and **augmented**.
  - Bradshaw number B = S(S+1), most negative (strongest destabilization) at S = −½, where
    |dU/dy| = 4Ω.
  - The plan's `B = S(S−1)`, unstable on 0 < S < 1, is the same physics with S = **+**2Ω/(dU/dy).
    Pick one convention.
  - Brethouwer 2017, JFM 828:424: https://arxiv.org/abs/1612.00254
  - Bradshaw 1969, JFM 36:177: https://doi.org/10.1017/S0022112069001583
  - Johnston, Halleen & Lezius 1972, JFM 56:533: https://doi.org/10.1017/S0022112072002502
  - Kristoffersen & Andersson 1993, JFM 256:163: https://doi.org/10.1017/S0022112093002757
  - Pedley 1969 (the same criterion as 2Ω(2Ω − U′) < 0, a narrow-gap Rayleigh limit):
    https://doi.org/10.1017/S002211206900098X
  - Tritton 1992 (cyclonic shear stabilizes; anticyclonic destabilizes at low rotation and
    restabilizes at high): https://doi.org/10.1017/S0022112092002131
  - Lezius & Johnston 1976 (laminar Re_c = 88.53 at Ro = 0.5; turbulent roll onset at
    Ro ≈ 0.022): https://doi.org/10.1017/S0022112076001171
  - Grundestam, Wallin & Johansson 2008 (full laminarization at Ro ≳ 3.0):
    https://doi.org/10.1017/S0022112007000122
  - Analogs: spinward ≡ rotating-channel pressure side ≡ Taylor–Couette with inner-cylinder
    rotation ≡ a concave Görtler wall. Brauckmann, Salewski & Eckhardt 2016:
    https://doi.org/10.1017/jfm.2015.737
  - Port Earth equatorial-boundary-layer results by vorticity sense, never by compass direction.
    "Up" points the other way relative to the axis, so Earth's destabilized case (westward)
    corresponds to the habitat's spinward. Zikanov, Slinn & Dhanak 2003:
    https://doi.org/10.1017/S0022112003006244 [MED]
- **Channel thresholds [HIGH, Brethouwer 2017]**, with Ro_b = 2Ωh/U_b:
  - Roll cells on the unstable side at Ro ≈ 0.15–0.45.
  - The unstable side self-organizes to dU/dy ≈ 2Ω (zero absolute vorticity) up to Ro ≈ 1.5.
  - The stable side is laminar-like for Ro ≥ 1.2, and the unstable side is viscous-dominated by
    Ro ≳ 2.
  - Full relaminarization at Ro ≈ 3.0 (linear theory Ro_c = 2.80 at Re = 10⁴).
- **Habitat mapping (our derivation, cross-checked two ways) [HIGH for signs].** The mean
  vorticity of U(z)x̂ is (dU/dz)ŷ, and Ω⃗ = −Ωŷ.
  - **Spinward wind** (dU/dz > 0): vorticity anti-parallel to Ω, so anticyclonic, S < 0.
    - Destabilized **only where dU/dz > 2Ω = 0.111 s⁻¹**, strongest where dU/dz ≈ 4Ω.
    - Where dU/dz < 2Ω, S < −1 and it is **stabilized**.
    - Rayleigh cross-check: L_ang = r(Ωr + u) gives dL/dr ≈ R(2Ω − dU/dz) + u, which is negative
      exactly when dU/dz > 2Ω.
  - **Antispinward wind**: vorticity parallel to Ω, so cyclonic, S > 0. **Always stabilized.**
  - **Axial wind**: vorticity ⟂ Ω, so it is Bradshaw-neutral (it is the streamwise-rotation case,
    Q6).
  - **Verdict on claim (c): CONFIRMED with a correction.** Spinward is destabilized and
    antispinward damped, but only within the near-floor layer where the shear exceeds 2Ω.
- **Height of the unstable layer [MED].** With log-law shear u*/(κz):
  - dU/dz = 2Ω at z ≈ u*/0.044 s⁻¹, i.e. 2.3 m for u* = 0.1 m/s and 6.8 m for u* = 0.3 m/s.
  - Peak destabilization sits at about half that height.
  - At canopy top the shear can be ≈ U_h/L_s ≈ 0.4 s⁻¹, giving S ≈ −0.3 (near the peak).
  - Above a few metres, spinward flow is rotation-stabilized.
- **Bulk scale [LOW].** For H = 100 m and U = 2 m/s, Ro_b ≈ 5.5, beyond the channel
  relaminarization threshold. Whole-depth rotational suppression of shear turbulence is
  possible but unproven (no confining wall; stratification present).
- **For L2b:** replace the per-direction stability flag with **Ri_rot(z) = S(S+1),
  S = −2Ω/(dU/dz)**. Add thermal N² separately.

## Q6. Streamwise rotation: axial-wind veer (H2)

- **Literature [HIGH].** Oberlack, Cabot, Pettersson-Reif & Weller 2006, JFM 562:383 (channel
  rotating about the streamwise axis): rotation induces a mean spanwise cross-flow W(y), and all
  six Reynolds stresses become non-zero. https://doi.org/10.1017/S0022112006001121 · critique of
  their §2.3 symmetry analysis (not of the DNS): https://arxiv.org/abs/1609.08155
- **Magnitudes [MED].** One agent computed these from the public TU Darmstadt DNS files
  (Re_τ = 180):

  | Ro_τ | 1 | 4 | 20 | 25 | 65 |
  |---|---|---|---|---|---|
  | peak W/U_b | 0.029 | 0.071 | 0.095 | 0.054 | 0.084 |
  | wall veer angle | 2.2° | 5.6° | 11.2° | 4.2° | 6.8° |

  - The trend is non-monotonic. The veer decays by y⁺ ~ 30–50 and reverses weakly (−1° to −2°)
    in the outer layer.
  - Data: https://www.fdy.tu-darmstadt.de/media/fachgebiet_fdy/dns_daten/Re180_Ro20.txt
- **H2 assessment [MED].**
  - There is no laminar Ekman spiral.
  - A turbulence-mediated veer of a few degrees (2–11°) and a few percent of U exists near the
    floor for axial wind.
  - **The veer direction is unverified [LOW]:** the DNS sign convention was not pinned down.
  - Extrapolating from Re_τ = 180 to atmospheric Re is [LOW].

## Q7. Supergravitational / centrifugal-buoyancy convection, and claim (f)

- **Definitions (Sun/Lohse group, ACRBC) [HIGH].**
  - Ra = ½ω²(R_o+R_i)βΔL³/(νκ), using the mean radius and gap L.
  - Ro⁻¹ = 2ωL/U_ff = 2(βΔ r_m/L)^{−1/2}, which does not depend on ω.
  - Hot outer, cold inner wall: the same sense as a heated floor.
- **Jiang, Zhu, Wang, Huisman & Sun 2020, Sci. Adv. 6:eabb8676 [HIGH].**
  - η = 0.5, water (Pr 4.3), 5–60 g (**not** 100 g), Ra 10⁶ to 6.5×10¹⁰ (**not** 10¹²).
  - Nu ∝ Ra^{0.27±0.01} below 10¹⁰; the exponent exceeds 1/3 only above 4×10¹⁰.
  - Coriolis regimes:
    - I, Ro⁻¹ ≲ 0.1–0.15: 3D, Nu independent of rotation;
    - II: Nu falls as Taylor–Proudman takes hold;
    - III, Ro⁻¹ ≳ 10–15: quasi-2D, axial velocity ≈ 0.
  - Turbulent rolls: **4 pairs, axis-parallel, aspect ratio ≈ 1**, drifting **prograde**
    (curvature asymmetry favors hot plumes).
  - https://doi.org/10.1126/sciadv.abb8676 † · https://arxiv.org/abs/2010.07189
- **Follow-ups [HIGH unless noted].**
  - Jiang, Wang, Liu & Sun 2022, PRL 129:204502: Nu ∝ Ra^{0.40±0.01} for Ra 10¹¹–3.7×10¹²
    (ultimate regime). https://doi.org/10.1103/PhysRevLett.129.204502 · https://arxiv.org/abs/2110.01224
  - Wang et al. 2022, JFM 930 A19 (radius ratio): drift rate ω_d = −0.2η + 0.15 rad per free-fall
    time for η 0.3–0.7, going to 0 as η → 1. https://doi.org/10.1017/jfm.2021.889 · https://arxiv.org/abs/2104.03709
  - Ren, Zhong, Lai & Sun, PRL 136:214002 (2026): early ultimate-state onset attributed to
    Stewartson layers. https://arxiv.org/abs/2509.14800 [MED]
  - Zhong & Sun 2026 (preprint): a no-slip outer wall with a stress-free inner wall gives a strong
    zonal flow and Nu ∼ Ra^{0.1}. This is relevant because the habitat has no inner wall.
    https://arxiv.org/abs/2601.21396 [MED]
- **Rouhi, Lohse, Marusic, Sun & Chung 2021, JFM 910 A32** (thin shell, Pr = 0.7, i.e. air):
  [HIGH]
  - A bidirectional azimuthal wind peaks at optimal Ro⁻¹ ≈ 1.0 (Ra 10⁷) and 0.6 (Ra 10¹⁰), where
    Nu is at its minimum.
  - Beyond that optimum, an axial quasi-2D cyclone forms.
  - At large Ro⁻¹ the flow is geostrophic and Nu drops.
  - https://doi.org/10.1017/jfm.2020.959
- **Onset with flat end walls: this corrects the plan [HIGH].**
  - Axis-parallel rolls satisfy Taylor–Proudman exactly. Coriolis is balanced by pressure, so
    **onset equals the non-rotating Rayleigh–Bénard value**, with no E dependence.
  - In the annulus, Ra_c → **1708** (rigid walls) as η → 1 and is somewhat higher at smaller η.
  - At η = 0.5, Ro⁻¹ = 20, Pr = 4.3, the critical mode is n_c = 5 roll pairs, with axial
    wavenumber 0.
  - Zhong, Wang & Sun 2024, JFM 992 A16: https://doi.org/10.1017/jfm.2024.543 · https://arxiv.org/abs/2402.00779
  - Westerburg & Busse 2003: https://doi.org/10.5194/npg-10-275-2003
  - Busse & Simitev review: https://arxiv.org/abs/0906.1124
  - Pitz et al. 2017 tabulates Ra_c(η): https://doi.org/10.1017/jfm.2017.451 [MED, not read]
- **Linear onset drift is retrograde, turbulent drift prograde [MED].** Kang et al. 2019, PRF 4
  043501, abstract only: https://doi.org/10.1103/PhysRevFluids.4.043501
- **Busse thermal Rossby waves [HIGH].**
  - Ra_c ≈ 3η_P^{4/3}, which is ∝ E^{−4/3}, with wavenumber ∝ E^{−1/3}.
  - These **require a β-effect from sloping end caps.**
  - The habitat's only β source is compressional: ≈ 0.19 density scale heights rim-to-axis. It is
    weak, and its importance at high supercriticality is contested.
  - Busse 1970, JFM 44:441: https://doi.org/10.1017/S0022112070001921
  - Busse & Carrigan 1974, JFM 62:579: https://doi.org/10.1017/S0022112074000814
  - Evonuk 2008: https://arxiv.org/abs/0711.4289 ; contrary view: https://arxiv.org/abs/1203.4145
- **Habitat placement [MED].** With a 1 K perturbation, Ro⁻¹ ≈ 19 over 1 km and ≈ 34 over the full
  3.2 km. That sits inside the experimental range (18–58), so expect regime III, quasi-2D.
- **Three-fold (non-uniform azimuthal) heating:** nothing found. This is an open gap [LOW].
- **L3 validation targets.**
  - Onset ≈ the non-rotating annular RB threshold for the run's η.
  - Roll aspect ratio ≈ 1.
  - Ro⁻¹ ≫ 15 implies quasi-2D.
  - Weak curvature-driven drift as η → 1.

## Q8. 2D inverse cascade arrested by drag (H1), and claim (d)

- **Coriolis in 2D [HIGH].** For flow in the plane ⟂ Ω that is invariant along Ω, ∇×(Ω×u) = 0.
  Coriolis is a gradient absorbed into pressure. Lithwick 2007: https://arxiv.org/abs/astro-ph/0702046
  Axis-invariant modes, including those with axial velocity (2D3C), are "unaffected by
  rotation" (van Kan & Alexakis 2020).
- **Friction arrest [HIGH].**
  - Scale ℓ_α ≃ ε^{1/2} α^{−3/2}, with the inverse-cascade Kolmogorov constant C ≈ 6 (range 4–6).
  - The prefactor matters: k_a = A_l(μ³/ε)^{1/2} with A_l ≈ 35–76. So the arrest *wavelength* is
    **(0.08–0.18)·(ε/μ³)^{1/2}**.
  - The bare formula overestimates the eddy size by 5–12×.
  - Boffetta 2007: https://doi.org/10.1017/S0022112007008014 · https://arxiv.org/abs/nlin/0612035
  - Grianik, Held, Smith & Vallis 2004 (Eq. 11): https://doi.org/10.1063/1.1630054 †
  - Boffetta & Ecke 2012 review: https://doi.org/10.1146/annurev-fluid-120710-101240 †
- **Quadratic drag [HIGH].** −C_d|u|u, with C_d = c_d/h, acts as r_eff = C_d·|v|_rms and selects an
  ε-independent arrest wavenumber k_a ≈ 51 C_d (DNS). (Grianik et al.)
- **Critical rotation for 2D-ization depends on aspect ratio [HIGH].**
  - van Kan & Alexakis 2020, JFM 899 A33: an inverse cascade appears for
    λ = ℓ_in^{5/3}Ω/(H ε^{1/3}) > λ_c ≈ 0.03, with a continuous transition.
    https://doi.org/10.1017/jfm.2020.443 · https://arxiv.org/abs/1912.05394
  - Seshasayanan & Alexakis 2018: condensates form below Ro_f ≈ 0.6 at H/ℓ_f ≈ 25.
    https://doi.org/10.1017/jfm.2018.106 · https://arxiv.org/abs/1707.03800
  - Deusebio et al. 2014: a taller domain weakens the inverse flux.
    https://doi.org/10.1103/PhysRevE.90.023005 · https://arxiv.org/abs/1401.5620
  - Review: Alexakis & Biferale 2018, https://arxiv.org/abs/1808.06186
- **Habitat numbers (H = 20 km axis) [LOW, assumed ε and ℓ_in].**
  - ℓ_in = 1 km, ε ~ 10⁻³ W/kg: λ ≈ 2.8, a robust inverse cascade.
  - ℓ_in = 100 m, ε ~ 10⁻²: λ ≈ 0.03, **right at threshold**.
  - So H1 holds for km-scale forcing (strip heating) but is marginal for 100 m-scale plumes.
- **The geometry corrects H1's drag picture [MED].**
  - In the (r,θ) disc, the canopy is a sidewall annulus covering ≈ 2h_c/R ≈ 1.25 % of the area.
    It is not a uniform Rayleigh drag.
  - Domain-averaged μ ≈ c_d·a·U·(2h_c/R) ≈ 2.5×10⁻⁴ s⁻¹ for illustrative values.
  - The drag acting uniformly on Taylor columns is **end-cap Ekman friction**: μ_E ≈ 2√(ν_tΩ)/L
    laminar-form [LOW], or C_D·U/(L/2) turbulent. Both are weak at L = 20 km.
  - Expect a condensate or wall-attached circulation to be likelier than a cleanly
    friction-arrested cascade.
- **Claim (d): partly CONFIRMED.** Taylor–Proudman 2D-ization and the inverse cascade are robust
  when λ ≫ 0.03. "Arrested by canopy drag" is the wrong mechanism for the bulk.

## Q9. Roughness length, displacement height, and drag coefficients

- **Davenport–Wieringa z₀ classes [HIGH]:**

  | Terrain | z₀ (m) |
  |---|---|
  | Sea | 0.0002 |
  | Mud flats | 0.005 |
  | Open grass | 0.03 |
  | Low crops | 0.10 |
  | High crops | 0.25 |
  | Parkland (class 6) | 0.5 |
  | **Suburb / forest (class 7)** | **1.0** |
  | City centre (class 8) | ≥ 2 |

  - WMO CIMO guide: https://www.eoas.ubc.ca/courses/atsc303/Instruments/wmo_guides/CIMO_Guide_2014-Met_Site_Classification.pdf
  - Wieringa 1992: https://doi.org/10.1016/0167-6105(92)90434-C
- **Suburban low-rise (h ≈ 10 m) [MED].** z₀ ≈ 0.3–0.8 m and d ≈ 3–7 m (d ≈ 0.5–0.7h,
  z₀ ≈ 0.1h, about 2× scatter).
  - Grimmond & Oke 1999: https://doi.org/10.1175/1520-0450(1999)038%3C1262:APOUAD%3E2.0.CO;2 †
  - Macdonald, Griffiths & Hall 1998 (morphometric, A ≈ 4.43, C_D ≈ 1.2):
    https://doi.org/10.1016/S1352-2310(97)00403-2
- **Forest/orchard [HIGH rule, MED exact].** d ≈ 0.6–0.8h (2/3h), z₀ ≈ 0.1h (range 0.05–0.13h).
  For h = 10–20 m: d ≈ 7–13 m and z₀ ≈ 1–2 m; orchards ≈ 0.5–1 m.
  - Raupach 1994: https://doi.org/10.1007/BF00709229
  - Raupach 1992: https://doi.org/10.1007/BF00155203
- **Canopy drag [MED].**
  - Element C_d ≈ 0.2 (range 0.1–0.3).
  - LAI: deciduous 3–6, conifer 4–9, orchard 1.5–4. a ≈ LAI/h ≈ 0.2–0.5 m⁻¹.
  - Field drag length ℓ_d = 1/(C_d a): median ≈ 21 m for trees and 0.7 m for crops.
    https://arxiv.org/abs/2411.01570 (abstract verified)
  - Canopy turbulence reviews: Finnigan 2000, https://doi.org/10.1146/annurev.fluid.32.1.519 † ;
    Raupach, Finnigan & Brunet 1996 (mixing-layer analogy, Λ_x ≈ 8.1 L_s [MED]),
    https://doi.org/10.1007/BF00120941 ; Belcher, Harman & Finnigan 2012,
    https://doi.org/10.1146/annurev-fluid-120710-101036 †
- **Surface-layer statistics [MED].**
  - κ = 0.40.
  - u*/U at 10–15 m above d: ≈ 0.12 suburb, ≈ 0.17 forest, ≈ 0.07 open grass (log law).
  - σ_u/u* ≈ 2.2–2.5.
  - Gust factor (3 s peak over 10 min mean) ≈ 1.4–1.5 open and **1.7–2.0 suburban/forest**. The
    primary gust-factor tables (ESDU, Harris & Deaves, WMO 2010) were not retrieved. Wieringa
    1973: https://doi.org/10.1007/BF01034986 (content not read).
- **Roughness sublayer ≈ 2–3h, i.e. 40–60 m over a 20 m forest [MED].** The log law and
  Monin–Obukhov similarity fail inside it.
- **Rotation matters at canopy scale here [MED, derivation].**
  - At U = 1 m/s, 2Ω/(C_d a U) ≈ 2.3 under forest (≈ 0.002 on Earth).
  - The rotation time 1/(2Ω) is comparable to the canopy eddy turnover, h/U_h ≈ 10 s.
  - No canopy study with 2Ωh/U_h ~ 1 exists.

## Q10. Day/night heating of window and land strips

- **Fluxes over land [MED].**
  - Daytime peak sensible heat flux H: forest/crops 100–300 W/m², suburban 150–350, dry/urban up
    to ~400+.
  - Night: −10 to −60 W/m² rural; urban near zero or positive because stored heat is released.
  - Bowen ratio: irrigated 0.1–0.3, deciduous 0.3–0.6, conifer 0.5–1.5, suburban 0.7–2, urban core
    1.5–5.
  - Daytime H/Q* ≈ 0.2–0.4 vegetated and 0.3–0.5 suburban. Roughly **10–35 % of incident
    shortwave becomes H**.
  - Grimmond & Oke 2002: https://doi.org/10.1175/1520-0450(2002)041%3C0792:THFIUA%3E2.0.CO;2 †
  - Oke 1982: https://doi.org/10.1002/qj.49710845502 †
- **Window strips are transparent and produce H ≈ 0. Land strips carry all the surface heating.**
  The result is three-fold azimuthal forcing [HIGH by construction].
- **Convective velocity scale [HIGH formula, Deardorff 1970].**
  - w* = (g H z_i/(ρ c_p T))^{1/3}.
  - At ½ atm, ρc_p ≈ 603 J/m³K, so the same H gives 2× the kinematic flux and w* ≈ 1.26× Earth's.
  - With H = 200 W/m² and z_i = 1 km: **w* ≈ 2.2 m/s**.
  - https://doi.org/10.1175/1520-0469(1970)027%3C1211:CVATSF%3E2.0.CO;2 †
- **Rotation control of thermals [HIGH arithmetic].**
  - Convective Rossby number w*/(2Ωz_i) ≈ 0.02.
  - A 2 m/s updraft feels a spinward Coriolis acceleration of 0.22 m/s², about 7× the buoyancy of
    a 1 K thermal.
  - Thermals will be strongly deflected spinward.
- **Thermal wind geometry [MED, derivation].**
  - An *axial* temperature gradient gives Δu_θ ≈ (ΩR/2)(ΔT/T) ≈ 3 m/s per 10 K at the rim,
    independent of L.
  - Azimuthal gradients (strips) conflict with the end caps and Taylor–Proudman, so expect
    overturning or ageostrophic responses rather than balanced jets.
  - Radial gradients are pure stratification.
  - Rotating-annulus analogy (Hide, Fowlis, Read): Θ = ΔT/T ≈ 0.035 for 10 K, Ta ≫ 10¹⁰, which is
    the irregular/geostrophic-turbulence corner. https://doi.org/10.1098/rsta.1958.0004 †
    [MED, analogy only]

## Q11. Vortex shedding with Ω perpendicular to the obstacle axis

- **Non-rotating baselines [HIGH].** St ≈ 0.2 (circular), ≈ 0.12–0.14 (square/rectangular
  buildings). Okajima 1982: https://doi.org/10.1017/S0022112082003115
  - The plan's 0.06 Hz for a 10 m building at 3 m/s uses St = 0.2. With St = 0.12 it is
    0.036 Hz, versus 2Ω/2π = 0.018 Hz.
- **Ω ⟂ obstacle axis, the closest data [HIGH].**
  - Machicoane et al. 2018, PRF 3:034801 (horizontal cylinder towed in fluid rotating about the
    vertical).
  - Rotation stabilizes the wake up to Re_c ≈ (275 ± 25)/Ro, with Ro = U/(2Ωd).
  - At low Ro the wake becomes a stationary inertial-wave pattern.
  - **No Strouhal data.**
  - https://doi.org/10.1103/PhysRevFluids.3.034801 · https://arxiv.org/abs/1802.07602
- **Habitat buildings [MED].**
  - A 10 m building at 2–5 m/s has Ro = 1.8–4.5 and Re ~ 10⁶ ≫ Re_c ≈ 100, so **shedding
    persists**.
  - The shedding period (≈ 28 s) is comparable to the inertial period (57 s): expect wake tilting
    and inertial-wave radiation, not quenching.
  - Obstacles ≥ 50 m have Ro ≈ 0.4–0.9 and are rotation-controlled.
- **Ω parallel to the obstacle axis (island wakes) does NOT transfer [HIGH].** These results
  include anticyclone dominance and St up to 0.6 (Perret et al. 2006,
  https://doi.org/10.1063/1.2179387 † ; Teinturier et al. 2010,
  https://doi.org/10.1016/j.dynatmoce.2008.10.006). They are not applicable to vertical
  obstacles in a habitat. Boyer & Davies 2000 review: https://doi.org/10.1146/annurev.fluid.32.1.165 †
- **Better analog for azimuthal wind over axis-aligned ridges or building rows [HIGH]:
  spanwise-rotation separation.**
  - Savino & Wu 2024, JFM 999 A51: the separation bubble shrinks on the anticyclonic side and
    grows on the cyclonic side, and drag drops in all rotating cases.
    https://doi.org/10.1017/jfm.2024.898 · https://arxiv.org/abs/2403.13264
  - Barri & Andersson 2010: https://doi.org/10.1017/S0022112010004696
  - Spinward and antispinward winds should give systematically different lee wakes. This is
    untested for habitat geometry.
- **Gaps:** no measured St(Ro) for Ω ⟂ obstacle axis. The content of Boyer & Kmetz 1983
  (https://doi.org/10.1080/03091928308221763 †) is unverified.

## Q12. Navier–Stokes and Euler regularity as of 2026-10-12, and claim (h)

Every quotation below is copied from the primary document, fetched on 2026-10-12. Each item is
tagged **ESTABLISHED** (peer-reviewed, or an official statement of record) or **CLAIMED**
(posted, but unrefereed and unaccepted).

### 12.1 What the Millennium Problem actually asks

The Clay problem is Fefferman's official statement. It asks "for a proof of one of the
following four statements":

- **(A) and (B):** global existence and smoothness, on R³ and on the torus R³/Z³ respectively,
  for any smooth initial data, with the force "f(x, t) to be identically zero".
- **(C) "Breakdown of Navier–Stokes solutions on R³":** "Then there exist a smooth,
  divergence-free vector field u°(x) on R³ and a smooth f(x, t) on R³ × [0, ∞), satisfying (4),
  (5), for which there exist no solutions (p, u) of (1), (2), (3), (6), (7) on R³ × [0, ∞)."
- **(D)** is the same breakdown statement on R³/Z³.

The side conditions:
- Condition (5) requires f and all of its derivatives to decay: |∂ₓ^α ∂ₜ^m f| ≤ C(1 + |x| + t)^−K.
- Condition (6) requires the solution to be C^∞.
- Condition (7) requires bounded energy, ∫|u|² dx < C.

Fefferman adds: "These problems are also open and very important for the Euler equations
(ν = 0), although the Euler equation is not on the Clay Institute's list of prize problems."
https://www.claymath.org/wp-content/uploads/2022/06/navierstokes.pdf

**Consequence [ESTABLISHED].** A construction with a *smooth, decaying external force* that
produces breakdown is a prize-eligible answer as written. The popular line "a forced blowup
isn't the Clay problem" misreads the official text. Euler results of any kind cannot win the
prize.

### 12.2 The headline claim: OpenAI, 8 September 2026 [CLAIMED]

- **Who and when.** The manuscript is "Finite Time Blowup for Navier–Stokes", 166 pages,
  author line "OPENAI". It was posted on OpenAI's CDN on 8 Sep 2026; the date is confirmed by
  the Córdoba & Martínez-Zoroa post, Petrillo & Glimm, and others. It is not on arXiv. The
  openai.com announcement page returns 403 to scripted fetches, so it was not read.
  https://cdn.openai.com/pdf/32d9f210-8b73-45e0-91bc-82a30aef8a9a/navier-stokes.pdf
- **Equation and setting.** The setting is the viscous incompressible Navier–Stokes equations
  (ν > 0) on all of R³, with no boundary, plus an external force. The result is **blowup**, not
  regularity.
- **The theorem, verbatim.** "Theorem 1.1. For every ν > 0 there exist a force
  f ∈ C_c^∞(R³ × (0, ∞); R³), a compact set K ⊂ R³, and smooth velocity and pressure fields u,
  p on R³ × [0, 1) satisfying ∂ₜu + (u·∇)u − νΔu + ∇p = f, ∇·u = 0, u(·, 0) = 0, such that
  supp u(·, t) ∪ supp p(·, t) ⊂ K for every 0 ≤ t < 1, sup_{0≤t<1} ‖u(t)‖_{L²} < ∞,
  limsup_{t↑1} ‖u(t)‖_{L∞} = ∞."
- **The claimed Clay link, verbatim.** "This establishes alternative (C) in the Millennium
  problem statement for Navier–Stokes as stated by Fefferman … Compact support also yields the
  corresponding construction on T³ = R³/Z³, establishing alternative (D) … see Corollary
  10.6."
- **Plain reading.**
  - The fluid starts at rest.
  - A smooth push confined in space and time drives it.
  - The kinetic energy stays bounded, but the peak velocity becomes infinite at t = 1.
  - The mechanism, in the paper's words: "oscillatory pulses generate a mean momentum flux that
    supplies the missing force on a collapsing background vortex".
  - It builds on the forced-Euler and IPM "cascade" constructions of Córdoba & Martínez-Zoroa.
- **A second claim on the same date [CLAIMED; that manuscript was not read].** Córdoba &
  Martínez-Zoroa describe it this way: "OpenAI's September 8, 2026 announcement concerns two
  separate claims: blow-up for ordinary Navier–Stokes with smooth forcing, and blow-up for
  unforced Euler from smooth initial data."

### 12.3 What the Clay Institute actually said [ESTABLISHED: the statement of record]

The CMI news item "Navier-Stokes Announcement" is dated 11 September 2026. The operative
passages, verbatim:

> "The problem asks about the existence and smoothness of Navier-Stokes solutions in
> 3-dimensional Euclidean space. In recent years there has been an increasing sense of
> anticipation as breakthroughs in the surrounding field (some recognised by the Clay Research
> Award) have raised hopes that the Navier-Stokes problem might soon be resolved. The
> increasing ability of new technologies to accelerate mathematical research has heightened
> this sense of anticipation. Today, CMI shares in the excitement of the global mathematical
> community as we contemplate the announcement that the Navier-Stokes problem has apparently
> been settled. We hope to see waves of new human understanding unleashed as the innovations
> behind this work are analysed and interrogated."
>
> "The rules governing the prizes describe the process for evaluating what has been achieved
> and for assigning credit. The process is deliberately unhurried, but we will provide
> updates."

https://www.claymath.org/news/navier-stokes-announcement/

**What the notice does not do:**
- It does not name OpenAI or any author.
- It does not mention forcing.
- It does not declare a solution or award a prize.
- The Clay problem page still carries the status label **"Active"**:
  https://www.claymath.org/millennium/navier-stokes-equation/

**The prize rules, verbatim.** "(i) the proposed solution must be published in a Qualifying
Outlet (see §6), and (ii) at least two years must have passed since publication, and (iii)
the proposed solution must have received general acceptance in the global mathematics
community." So no award is possible before late 2028 at the earliest.
https://www.claymath.org/millennium-problems/rules/

### 12.4 Peer review and acceptance [ESTABLISHED as of 2026-10-12]

- **No journal publication and no referee report.** The OpenAI manuscript is not on arXiv and
  has no journal record.
- **No published refutation.** No paper I found shows Theorem 1.1 is false.
- **First community responses**, all preprints:
  - **Constantin, Ignatova & Vicol** (arXiv 2609.20803, 17 Sep 2026). They *constrain* the
    construction; they do not refute it. They prove that solutions with OpenAI's two structural
    features are regular if the force is real-analytic. "As a consequence … the force can
    neither vanish identically near the singular point, nor be real analytic in the space
    variables." https://arxiv.org/abs/2609.20803
  - **Lei & Ren** (arXiv 2609.35406, 28 Sep 2026). An expository rewrite of the
    profile-construction part. They "regard OpenAI's work as a major advance". It "will not be
    submitted to any journal", so it is not a verification. https://arxiv.org/abs/2609.35406
  - **Bastounis, Circelli & Hansen** (arXiv 2610.08144, 6 Oct 2026). They challenge the Lean
    check: "we show that the formalised Lean proof does not correspond to the NL
    [natural-language] proof of blow-up of solutions to the Navier-Stokes equations."
    https://arxiv.org/abs/2610.08144 [MED: one critique, not yet answered]
  - **Duraiswami** (arXiv 2609.17642). On the physics: the construction's cone condition
    "requires radii of order 10^20" in similarity variables; "A real fluid cavitates (water) or
    shocks (air) first"; the theorem leaves "the unforced equations of practice as they were."
    https://arxiv.org/abs/2609.17642
  - **Petrillo & Glimm** (arXiv 2609.23868). They summarize: "a forced Navier-Stokes
    singularity at every fixed viscosity, statements (C) and (D) of the Clay problem … The
    unforced problem, statements (A) and (B), stands open." https://arxiv.org/abs/2609.23868
  - **Córdoba & Martínez-Zoroa** (guest post on Tao's blog, 4 Oct 2026). The forced-NS claim
    "if verified under the stated hypotheses, would establish an allowed breakdown alternative
    in the Millennium Problem formulation. It is distinct from unforced Navier–Stokes blow-up.
    At this moment we are not in a position to say anything substantive about a comparison
    between OpenAI's constructions and our cascade methods."
    https://terrytao.wordpress.com/2026/10/04/on-classical-solutions-and-singularity-formation-in-incompressible-fluids/
- **Summary.** The forced alternatives (C)/(D) are *claimed* and are being taken seriously.
  Formal acceptance is at least two years plus publication away. The unforced alternatives
  (A)/(B), the question most specialists mean by "the Navier–Stokes problem", are **open**.

### 12.5 The role of AI and computer-assisted proof

These three things are different and should not be blurred:

| Kind | What it is | Examples here |
|---|---|---|
| **Computer-assisted proof** | A human-designed argument whose numerical estimates are certified by rigorous (interval) arithmetic | Chen–Hou |
| **ML-discovered candidates** | Neural networks find approximate singular profiles; these are numerics, not proofs | DeepMind et al. 2025; Ganeshram–Duruisseaux–Anandkumar 2026 |
| **AI-generated proof plus Lean formalization** | An LLM agent system writes the proof, and an autoformalized Lean version is machine-checked | OpenAI 2026 (claimed); Alpöge–Buckmaster 2026 ("LLM-assisted") |

- **OpenAI.**
  - *What the PDF itself says:* it contains **no mention of AI, Lean, or computer assistance**
    (checked by full-text search of the PDF). The author line is simply "OPENAI".
  - *What secondary sources add* [LOW, not confirmed from a primary OpenAI source]:
    "166 pages, 616,000 lines of Lean, produced in 88 hours by some 10,000 software agents".
    That wording is Gamburd's essay, quoting the OpenAI press release.
    https://arxiv.org/abs/2609.28591
  - *The open question this raises:* Bastounis et al. argue that the Lean statement does not
    match the paper's theorem. If they are right, the Lean check certifies a different
    statement from the one in the PDF.
- **Alpöge & Buckmaster (with Coiculescu on IPM), 7 Sep 2026 [CLAIMED].**
  - The IPM paper's arXiv comment reads: "LLM-assisted and Lean formalized proof". It proves
    IPM blowup on T² with a uniformly space-time smooth force. https://arxiv.org/abs/2609.16470
  - Tao: they demonstrated forced blowup for IPM, 2D Boussinesq and **3D Euler**, "the
    arguments here are heavily AI-assisted", and they "do not quite achieve" Navier–Stokes,
    although the method "has a high likelihood of also extending to Navier-Stokes".
    https://terrytao.wordpress.com/2026/09/07/finite-time-blowup-with-smooth-forcing-term-for-the-incompressible-porous-medium-boussinesq-and-incompressible-euler-equations/
  - The construction builds on Córdoba & Martínez-Zoroa (arXiv 2410.22920, human-authored).
  - A priority and credit dispute between this group and OpenAI is reported. I have only
    secondary accounts of it and have not adjudicated it [LOW].
- **DeepMind et al., "Discovery of Unstable Singularities"** (arXiv 2509.14185, 17 Sep 2025)
  [ESTABLISHED as a preprint; numerical].
  - Machine-learning architectures plus a high-precision Gauss–Newton optimizer find "new,
    unstable self-similar solutions for the incompressible porous media equation and the 3D
    Euler equation with boundary".
  - The precision "meets the requirements for rigorous mathematical validation via
    computer-assisted proofs". It is not itself a proof, and it has no journal reference.
  - Follow-up arXiv 2511.22819 (Nov 2025) reaches round-off-level residuals for more unstable
    profiles, "bridging the gap … for computer-assisted proofs". Still numerics.
  - https://arxiv.org/abs/2509.14185 · https://arxiv.org/abs/2511.22819

### 12.6 Relation to arXiv 2609.10867 and to Chen–Hou

| Work | Equation | Domain | Forcing | Status | Relation to the Clay problem |
|---|---|---|---|---|---|
| **Chen & Hou** | 3D axisymmetric **Euler** (and 2D Boussinesq), ν = 0 | **Bounded, with a solid boundary** | None, smooth finite-energy data | **ESTABLISHED in part**: computer-assisted proof. Part II published, *Multiscale Model. Simul.* 23(1) (2025). PNAS 122(27), 27 Jun 2025, is a *review* ("we review a recent computer-assisted proof"). Part I (arXiv 2210.07191, v4 Aug 2026) has no journal reference. | Not eligible: Euler, with a boundary. It is the first rigorous smooth-data Euler blowup, but the boundary "helps confine the singularity". |
| **Ganeshram, Duruisseaux & Anandkumar**, arXiv 2609.10867 (9 Sep 2026; v2 24 Sep) | 3D **Euler** | **R³, no boundary** | None | **CLAIMED, and explicitly not a proof.** "We provide evidence of a finite-time singularity"; a PINN finds a profile at blowup rate 0.5. The companion arXiv 2609.10860 is "Conditional on rigorous certification of the estimates and constants". | Not eligible (Euler), but it targets exactly the case Chen–Hou cannot reach: free space with no boundary. |
| **OpenAI** (8 Sep 2026) | **Navier–Stokes**, ν > 0 | R³ (and T³) | **Smooth, compactly supported force** | **CLAIMED**, unrefereed, Lean-fidelity challenged | Targets Clay (C)/(D) directly. |
| **Alpöge–Buckmaster** (7 Sep 2026) | IPM, Boussinesq, 3D Euler | T², R², R³ | Smooth force | **CLAIMED**, preprints, LLM-assisted + Lean | Not eligible; Tao: "do not quite" reach NS. |

- **Chen–Hou sources.**
  - Part II: https://doi.org/10.1137/23M1580395 † · https://arxiv.org/abs/2305.05660
  - PNAS review: https://doi.org/10.1073/pnas.2500940122 † (Crossref: PNAS 122(27), 2025-06-27) ·
    https://authors.library.caltech.edu/records/40zbk-gep55
  - Part I: https://arxiv.org/abs/2210.07191
- **The 2609.10867 group's account.** Anandkumar's guest post on Tao's blog says they released
  "the day before the OpenAI's announcement". It describes their route as self-similar rather
  than the "multi-scale structure with discrete jumps in scales" of the NYU and OpenAI
  constructions.
  https://terrytao.wordpress.com/2026/09/10/stable-singularity-of-the-euler-equations-on-r3/ ·
  https://arxiv.org/abs/2609.10867 · https://arxiv.org/abs/2609.10860
- **How the pieces fit.**
  - Chen–Hou: Euler with a boundary, no forcing. Proved.
  - 2609.10867: Euler without a boundary, no forcing. Evidence plus a conditional framework.
  - Alpöge–Buckmaster: Euler, with a force. Claimed.
  - OpenAI: Navier–Stokes, with a force. Claimed.
  - Still open, and not claimed by anyone: Navier–Stokes, no force. That is Clay (A)/(B).

### 12.7 Context results [ESTABLISHED]

- Elgindi 2021, *Annals* 194(3), C^{1,α} (not C^∞) Euler blowup on R³:
  https://doi.org/10.4007/annals.2021.194.3.2
- Tao 2016, *JAMS* 29(3), blowup for an *averaged* NS. This is a barrier result: any regularity
  proof must use the fine structure of the nonlinearity. https://arxiv.org/abs/1402.0290
- Buckmaster–Vicol 2019 (*Annals*), non-uniqueness of weak solutions:
  https://arxiv.org/abs/1709.10033
- Albritton–Brué–Colombo 2022 (*Annals*), non-unique Leray–Hopf solutions with a force that is
  singular at t = 0: https://arxiv.org/abs/2112.03116
- Hou 2022/2023 (*Found. Comput. Math.*): numerical NS and Euler blowup *candidates*, not proofs.

### 12.8 What this means for an engineer modelling habitat weather

Nothing in the models changes. Three reasons:

1. **The singularities are engineered, not generic.** Every construction needs one of: a solid
   boundary with exact symmetry (Chen–Hou), a C^∞ force tuned to feed a collapsing vortex
   (OpenAI), or a profile that infinitesimal perturbations destroy (DeepMind's unstable
   families). Duraiswami finds OpenAI's mechanism needs similarity radii of order 10²⁰, and that
   real air would shock first. Habitat weather is noisy, generically forced and viscous.
2. **Closures never assumed global smoothness.** RANS, LES and the reduced-order models here
   (L0–L3) model averaged or filtered statistics at finite Reynolds number. They are calibrated
   against experiment and DNS. A blowup, if real, sits below the Kolmogorov scale and below
   the continuum hypothesis's own range of validity. The Kolmogorov scale (ν³/ε)^{1/4} here is about 1–2 mm, versus the metre-to-km
   cells we use.
3. **The answer we need is statistical.** "Gentle breezes or hurricanes" depends on energy
   input, rotation, stratification and drag, which are the parameters in
   `model-parameters.md`. Whether an exact solution can become singular in some measure-zero
   corner of data space does not bear on it.

The honest one-liner for `docs/atmosphere.md`: *"As of October 2026, a forced Navier–Stokes
blowup has been claimed (OpenAI, unrefereed) and Clay says the problem has 'apparently been
settled'. The unforced question remains open. None of this affects engineering turbulence
models, which never relied on it."*

## Claims (a)–(h): summary

| Claim | Verdict |
|---|---|
| (a) Rising parcel drifts spinward ≈ 2ΩΔr (11 m/s per 100 m) | **CONFIRMED for ballistic parcels and axisymmetric rings** [HIGH]. Angular momentum and −2Ω×v agree; the exact value is 11.26 m/s per 100 m. Sinking parcels drift antispinward. **For thermals embedded in air it is a ceiling** [HIGH principle, LOW magnitude]: a plume that is uniform along the axis feels a Coriolis force that is a pure gradient (the claim-d argument), absorbed by its pressure perturbation, so it gets no net spinward push. Compact 3D thermals fall in between. This matches the plan's "upper bound for compact 3D thermals" framing for L2. |
| (b) No Ekman spiral on the floor; Ekman layers on end caps, Stewartson on the curved wall | **CONFIRMED** [HIGH]. Not the spherical E^{2/5} equatorial layer. |
| (c) Spanwise rotation: anticyclonic destabilized, cyclonic stabilized | **CONFIRMED with correction** [HIGH]. Spinward is anticyclonic and unstable only where dU/dz > 2Ω (bottom ~1–7 m), stable above. Antispinward is always stable. Use B = S(S+1), S = −2Ω/(dU/dz). |
| (d) Taylor–Proudman 2D-ization; inverse cascade arrested by drag | **PARTLY** [MED]. Robust for km-scale forcing (λ ≫ 0.03), marginal for 100 m forcing over a 20 km axis. The bulk drag is end-cap Ekman, not canopy. Arrest wavelength ≈ 0.08–0.18 (ε/μ³)^{1/2}. |
| (e) Greenspan–Howard spin-up with turbulent Ekman layers | **CONFIRMED form** [HIGH], with L = half-length = 10 km. The turbulent version τ ~ L/(C_D ΔU) is [MED] (derivation). Rim drag may dominate [LOW]. |
| (f) Centrifugal-buoyancy convection onset and rolls | **CORRECTED** [HIGH]. Flat end caps give non-rotating onset (Ra_c → 1708 as η → 1), axis-parallel rolls with aspect ratio ≈ 1, and quasi-2D at Ro⁻¹ ≳ 15. Not E^{−4/3}. |
| (g) Literature on rotating-habitat atmospheres | **NONE** [MED, negative result]. SP-413 and O'Neill 1974 have no dynamics. Only forum posts and an unexecuted OpenFOAM proposal exist. Jensen 2024 (https://arxiv.org/abs/2408.00152) covers design limits, including axis pressure. |
| (h) NS/Euler regularity | **UPDATED** [HIGH]. Forced NS blowup has been claimed (OpenAI, 8 Sep 2026, Clay C/D), is unrefereed, and its Lean fidelity is challenged. Clay says "apparently been settled" and has awarded nothing. Unforced A/B are open. Chen–Hou (Euler with a boundary) is proved; PNAS 2025 is a review. Engineering closures are unaffected. See Q12. |

## Could not verify

- A peer-reviewed turbulent spin-up time (only our quadratic-drag derivation), and the
  Eliassen & Lystad 1977 attribution.
- A ν_t number for spin-down at habitat Re (Kaiser extrapolation only).
- The Ekman fully-turbulent threshold Re_δ ≈ 150 (secondary; Buffett finds ≈ 500).
- The direction of the axial-wind veer (Oberlack DNS sign convention).
- Strouhal versus Ro for Ω ⟂ obstacle axis; Boyer & Kmetz 1983 content.
- Primary gust-factor tables; Macdonald constants; the exact Raupach curve.
- Ra_c(η) values (Pitz 2017 not read).
- Any study of three-fold azimuthally non-uniform heating.
- How strong the compressional β is in a 0.19-scale-height annulus.
- Details of the OpenAI manuscript (Lean size, agents, hours); the separate OpenAI unforced-Euler
  manuscript (not read); the OpenAI/Alpöge–Buckmaster priority dispute; whether Chen–Hou Part I
  has been accepted anywhere.
- Kalpana One atmosphere figures; the text of *The High Frontier* (1976).

## Tooling notes

- Gemini-backed research agents failed: `GOOGLE_API_KEY not set in ~/.claude/.env`. Those
  agents fell back to Claude WebSearch, Crossref and arXiv.
- The arXiv Atom API returned HTTP 429 during the ArXiv-skill queries. arXiv IDs were verified
  through abstract pages instead.
