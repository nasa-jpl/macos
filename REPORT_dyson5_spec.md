# REPORT — dyson5: the SPECTROMETER (addendum 49)

TO, 2026-10-07.  Brief: `BRIEF_to_dyson5.md` addendum 49.  Step 0, the
paper's units, is in `REPORT_dyson5_cprime.md` (Step 0).  This report holds
steps 1–2 (step 3, throughput by band, is CC's).

**Status:** step 1 done (3k and 1.5k ladders); step 2 done (Option B solved, controlled, laddered).

## What the Dyson of record IS — Option A, not "spherical/spherical"

The 3k Dyson of record (`dyson5_size.mat`, `size:F:240`, CaF2 240 mm,
emitted by `spectrometer_rx`) carries an aspheric lens: BlockSphereOut
`Surface= Aspheric`, K = 0.0432, h⁴ = −0.479, h⁶ = −5.99.  The grating is
spherical (`Surface= Conic`, K = 0).

That is the paper's **Option A**: aspheric Dyson lens + spherical grating
(Bradley 2024, Table 2 and Sec. 3, read from the page).  Its **Option B**,
the final design, is aspheric lens + conic grating.  Addendum 49 describes
the record as "spherical lens, spherical grating"; it is not.  Step 2 is
therefore one new DOF, the grating conic, on top of the record's asphere.

Paper Table 2 (DSI with an even-asphere TMA; smile in % of a co-added
pixel, keystone in % of a pixel):

| | A design | A CBE | B design | B CBE |
|---|---|---|---|---|
| smile < 5 % | 0.3 % | 3.5 % | 1.3 % | 3.5 % |
| keystone < 10 % | 0.4 % | 5.7 % | 2.7 % | 5.8 % |
| SRF < 1.8 co-added | 1.37 | 1.51 | 1.33 | 1.51 |
| CRF < 2.8 | 1.52 | 2.21 | 1.35 | 2.04 |

"Option A is more sensitive, having resulted in similar Smile and Keystone
predictions with worse performing CRF" (after tolerancing).

## Step 1 — the tolerance ladder, 3k Dyson of record

**Method.**  `design/src/spectrometer_sens` (+ `_defaults`, `_dn`), runner
`challenges/dyson5/dyson5_sens_run('3k')`, record
`dyson5_sens_3k.{txt,mat}`.
- **One perturbation at a time** on the deck TEXT.  The lens is one rigid
  part: both passes of its face and sphere, pivoting about the slit face.
  The grating pivots about its vertex, the detector about the FPA vertex,
  and the slit moves as the launch.
- **Scored twice:** by the record's own `spectrometer_score` (7 × 7, alone),
  and END TO END through the R9 join (`dyson5_t5f`, centroid launch = slit
  filled, roll 0; the template's Dyson blocks perturbed the same way —
  same frame as the alone deck, checked).
- **Compensators:** detector FOCUS (along its normal, minimising the
  mean-square ray width) and detector CLOCK (about its normal, minimising
  smile + keystone).  Residuals are judged against the nominal at the SAME
  compensator's optimum, because the record's deck is not at that focus:
  refocusing it alone takes CRF 1.213 → 1.175.
- **Gates:** `tSpectrometerSens` (SUITE_FAST) — a zero perturbation gives
  zero rows, the focus compensator recovers a pure detector defocus, focus
  + x/y == focus, and a real perturbation moves a metric.

**Nominals.**
- Alone: smile 0.0050, keystone 0.0058, CRF 1.213, SRF 2.0238, EE 0.824 px.
- At the focus optimum: CRF 1.175, EE 0.854.
- e2e: smile 0.0222, keystone 0.0058, CRF 1.174, SRF 2.0248, EE 0.849 —
  R9's centroid row, reproduced.

**The table** (Δ in px of 18 µm per the stated amount; "focus" = the
residual after the detector refocus vs the refocused nominal):

| row | amount | Δsmile | Δkeystone | ΔCRF | ΔSRF | ΔEE | after focus ΔCRF / ΔEE | e2e ΔCRF |
|---|---|---|---|---|---|---|---|---|
| lens decenter x (along slit) | 10 µm | 0 | **+0.0099** | +0.0005 | 0 | 0 | 0 / 0 | +0.0006 |
| lens decenter y (across) | 10 µm | +0.0009 | 0 | 0 | 0 | 0 | −0.002 / +0.001 | 0 |
| lens tilt x | 10 µrad | 0 | 0 | +0.0001 | 0 | 0 | 0 / 0 | 0 |
| lens tilt y | 10 µrad | 0 | 0 | 0 | 0 | 0 | 0 / 0 | 0 |
| lens radius | dR/R 1e-4 | 0 | 0 | **+0.261** | +0.001 | **−0.145** | +0.0002 / 0 | +0.241 |
| lens index | dn 1e-5 | 0 | 0 | +0.0014 | 0 | −0.0008 | +0.0003 / 0 | +0.0005 |
| lens thickness | 10 µm | 0 | 0 | **−0.114** | 0 | +0.095 | −0.003 / +0.002 | −0.064 |
| grating decenter x | 10 µm | 0 | **+0.0099** | +0.0015 | 0 | −0.0008 | +0.0006 / 0 | +0.0002 |
| grating decenter y | 10 µm | −0.0007 | 0 | −0.0004 | 0 | 0 | +0.0005 / 0 | +0.0002 |
| grating tilt x | 10 µrad | −0.0007 | 0 | +0.0003 | 0 | 0 | +0.0001 / 0 | +0.0002 |
| grating tilt y | 10 µrad | 0 | **+0.0074** | +0.0014 | 0 | −0.0008 | +0.0006 / 0 | +0.0002 |
| grating radius | dR/R 1e-4 | 0 | +0.0003 | **+1.932** | **+0.415** | **−0.699** | −0.050 / +0.039 (SRF 0) | +1.889 |
| grating period | 1e-5 | 0 | 0 | 0 | 0 | 0 | 0 / 0 | 0 |
| grating clock | 10 µrad | 0 | **+0.0029** | 0 | 0 | 0 | 0 / 0 | 0 |
| air space (lens–grating) | 10 µm | 0 | 0 | **−0.171** | 0 | +0.176 | +0.005 / −0.002 | −0.116 |
| slit x (along) | 10 µm | 0 | 0 | +0.0015 | 0 | −0.0008 | +0.0007 / 0 | +0.0009 |
| slit y (across) | 10 µm | 0 | 0 | +0.0004 | 0 | 0 | +0.0003 / 0 | 0 |
| slit defocus | 10 µm | +0.0001 | −0.0009 | **+0.261** | +0.001 | −0.146 | +0.0001 / 0 | +0.242 |
| detector defocus | 10 µm | −0.0001 | +0.0010 | **+0.261** | +0.001 | −0.146 | 0 / 0 | +0.242 |
| detector tilt x | 10 µrad | 0 | 0 | −0.0016 | 0 | +0.0016 | +0.0006 / 0 | −0.0007 |
| detector tilt y | 10 µrad | 0 | 0 | +0.0084 | 0 | −0.0064 | +0.006 / −0.003 | +0.006 |

Alone vs e2e agree in every sign and nearly every magnitude.  The e2e focus
terms are 25–40 % smaller: the telescope's own focus position offsets part
of the spectrometer's defocus.

**What it says.**
1. **Every big row is a FOCUS term, and detector refocus removes it.**
   Lens radius, thickness, air space, slit defocus, detector defocus and
   grating radius all cost CRF/EE.  After refocus each is ≤ 0.006 px CRF,
   and SRF is 0.  The grating radius is the strongest: dR/R = 1e-4 =
   79 µm needs a 76 µm refocus.  A one-axis detector focus stage is the
   one compensator this spectrometer needs.  Detector x/y buys nothing
   (`focus + x/y` CRF 1.1746 vs `focus` 1.1748).  Every metric here is
   relative, so a translation cannot move it.
2. **Keystone is set by DECENTERS ALONG THE SLIT and grating tilt y/clock,
   and no compensator removes it.**  Detector clocking leaves those
   keystone rows unchanged: they are lateral colour across the band, not a
   rotation of the spectral axis.  Lens and grating decenter x each give
   0.99e-3 px/µm, grating tilt y 0.74e-3 px/µrad, clock 0.29e-3 px/µrad.
3. **Smile is gentle.**  Decenters and tilts across the slit give
   0.7–0.9e-4 px per µm or µrad.
4. **Period and slit translations are calibration terms.**  A period change
   moves wavelengths, a slit shift moves the spatial origin; the relative
   metrics are blind to both, by construction.

**Linearity** (Δ at 2× / (2Δ at 1×), on the rows that dominate):
- keystone vs lens decenter x: 1.10 (linear);
- CRF vs grating radius: 1.03 (linear, uncompensated);
- SRF vs grating radius: 2.61 (quadratic: a blur added to the slit floor).

Smile and keystone are max − min of a pattern, so near nominal a small
perturbation can lower them first.  On grating decenter y the smile ratio
is −0.13.  Linear tolerances below hold once a perturbation's pattern
dominates the nominal one.

**The budget** (equal RSS allocation over each metric's contributors;
allowance = target − e2e nominal):

| metric | target | allowance | contributors (sensitivity) | tolerance each |
|---|---|---|---|---|
| keystone | paper A CBE 5.7 % = 1.03 µm = 0.057 px | 0.051 px | lens dx, grating dx (0.99e-3 px/µm); grating tilt y (0.74e-3 px/µrad); clock (0.29e-3 px/µrad) | **lens dx 26 µm, grating dx 26 µm, grating tilt y 35 µrad, clock 88 µrad** |
| keystone | Joe 0.1 px | 0.094 px | same | lens dx 47 µm, grating dx 47 µm, tilt y 64 µrad, clock 162 µrad |
| smile | paper A CBE 3.5 % = 1.26 µm = 0.070 px | 0.048 px | lens dy (0.9e-4 px/µm), grating dy, grating tilt x (0.7e-4) | lens dy ~310 µm, grating dy / tilt x ~400 µm / µrad |
| smile | Joe 0.1 px | 0.078 px | same | ~500 µm / µrad |
| CRF | paper A CBE 2.21 px | 1.03 px | focus-compensated residuals ≤ 0.006 px per 10 µm / 10 µrad | not limiting with a focus stage; without one, ~0.026 px per µm of any focus-equivalent |

So the tight tolerances are the **along-slit decenters (26 µm) and the
grating tilt about y (35 µrad)** for the paper's keystone CBE.  Joe's
0.1 px doubles them.  The focus terms are free with a detector focus stage.

**Thermal, in K** (the brief's two rows):
- **Soak via the index:** dn = 1e-5 is 1 K at dn/dT ≈ −1e-5/K (CC; CaF2's
  is about −1.0 to −1.1e-5/K over the band, not in this tree).  It costs
  ΔCRF +0.0014 px alone, +0.0005 e2e, at its focus: negligible.  The
  Sellmeier scaling that imposes it has a residual dispersion of 4.6e-7
  (dn 0.989–1.035 × 1e-5 over 0.38–2.5 µm).
- **The gradient-free soak via the lens radius:** at CaF2's CTE (≈ 18.9e-6/K,
  literature) 1 K is dR/R 1.9e-5, i.e. ΔCRF ≈ +0.050 px before refocus.
- **The same soak also grows the lens thickness** by 4.5 µm/K, i.e.
  ΔCRF ≈ −0.051 px.  The two nearly cancel, to first order, on the lens
  alone.
- **The structure terms are the open part:** the aluminium air space
  (≈ 12.9 µm/K over 547 mm → ΔCRF ≈ −0.22 px/K) and the N-BK7 grating
  radius (≈ 7.1e-6/K → ΔCRF ≈ +0.14 px/K).  They are focus terms of
  opposite sign, compensable by refocus.  This matches the paper's
  "CRF limits the DSI temperature requirements".  A proper soak row, every
  CTE at once, is the next row to add.

## Step 1 — the 1.5k Dyson of record (`dyson5_sens_1k5.{txt,mat}`)

Silica 130 (`size:D:130`), the 27 mm slit, through the 1.5k telescope of
record's join (`dyson5_tA_GM_1k5_bAs.in`, roll 180, centroid launch).

**Nominals.**
- Alone: smile 0.0065, keystone 0.0090, CRF 1.026, SRF 2.0241,
  EE 0.999 px.
- e2e: smile 0.088, keystone 0.0091, CRF 1.211, SRF 2.205, EE 0.295.  That
  telescope is the earlier aspheric one (2.2 px spots), not an R9-class
  freeform, and it dominates the e2e here.

| row | amount | Δsmile | Δkeystone | ΔCRF | ΔSRF | ΔEE | after focus ΔCRF | e2e ΔCRF |
|---|---|---|---|---|---|---|---|---|
| lens decenter x | 10 µm | 0 | **+0.0182** | +0.0002 | 0 | 0 | +0.0002 | +0.0001 |
| lens decenter y | 10 µm | +0.0008 | 0 | 0 | 0 | 0 | +0.0001 | +0.0001 |
| lens tilt x / y | 10 µrad | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| lens radius | dR/R 1e-4 | 0 | 0 | +0.079 | 0 | −0.092 | +0.0001 | −0.023 |
| lens index | dn 1e-5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| lens thickness | 10 µm | 0 | 0 | +0.003 | 0 | +0.001 | 0 | +0.025 |
| grating decenter x | 10 µm | 0 | **+0.0182** | +0.0002 | 0 | 0 | +0.0002 | +0.0002 |
| grating decenter y | 10 µm | −0.0008 | 0 | 0 | 0 | +0.0008 | 0 | +0.0001 |
| grating tilt x | 10 µrad | −0.0003 | 0 | 0 | 0 | 0 | 0 | −0.0004 |
| grating tilt y | 10 µrad | 0 | **+0.0053** | +0.0001 | 0 | 0 | +0.0001 | +0.0001 |
| grating radius | dR/R 1e-4 | 0 | +0.0002 | **+0.831** | +0.030 | **−0.559** | −0.0006 | +0.294 |
| grating period | 1e-5 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| grating clock | 10 µrad | 0 | **+0.0017** | 0 | 0 | 0 | 0 | 0 |
| air space | 10 µm | 0 | 0 | +0.010 | +0.0005 | +0.001 | +0.0002 | +0.065 |
| slit x / y | 10 µm | 0 | 0 | 0 | 0 | 0 | 0 | ±0.0001 |
| slit defocus | 10 µm | +0.0001 | −0.0013 | **+0.171** | +0.0005 | −0.167 | 0 | −0.035 |
| detector defocus | 10 µm | −0.0002 | +0.0013 | **+0.171** | +0.0005 | −0.166 | 0 | −0.035 |
| detector tilt x / y | 10 µrad | 0 | 0 | +0.0001 | 0 | ±0.0008 | 0 | +0.001 |

**Against the 3k:**
- **Keystone per µm of along-slit decenter is ~2× the 3k's** (0.0182 vs
  0.0099 px per 10 µm): the smaller Dyson maps the same decenter to more
  lateral colour per pixel.
- **The grating radius and the focus terms are weaker** (CRF +0.83 vs
  +1.93 per 1e-4).
- **Every focus term is again recovered by detector refocus** (≤ 0.0006 px).
- **Keystone again has no compensator.**

**Budget (1.5k):**
- **Paper A CBE keystone (0.057 px):** allowance 0.048 → lens dx and grating
  dx **13 µm** each, grating tilt y **45 µrad**, clock **141 µrad**.
- **Joe's 0.1 px:** 25 µm, 25 µm, 86 µrad, 268 µrad.
- **Smile** stays loose (≥ 300 µm / µrad), as on the 3k.

(Linearity on this run used the default row pair, too small to read on
smile; the 3k's dominant-row linearity stands.)  Index row: residual
dispersion 8.4e-7 (silica's Sellmeier).

## Step 2 — Option B (the grating conic): `dyson5_optionB('3k')`, `dyson5_sens_3kB`

**The solve.**  `spectrometer_geom` gained `P.grat_Kc` (default 0, so every
existing chain is unchanged).  `dyson_ladder` gained rung **RB** (index 7):
R3's variables + the grating conic, warm from the record, with the size
trade's own settings (family F; 30 iterations).
- **The chain still equals the engine with a conic grating** (K 0.3,
  unsolved): centroids agree to 9e-3 px, against 2e-3 at K 0.
- **A bug found en route, fixed:** `dyson_ladder` seeded a meniscus corrector
  for every rung index ≥ 4, so the first Option B solve silently carried a
  meniscus.  Its seed was CRF 5.0 px, not the record's 1.21, and it ended at
  CRF 4.65; that run is kept as `*_meniscusbug`.  The seeding is now
  restricted to the meniscus rungs 4–6.
- **The CONTROL:** R3 re-solved from the record with no conic does not move
  (merit 0.4934, CRF 1.213).  The record is a converged Option A, so
  Option B's gain is the conic's.

| alone, engine | record (Option A) | Option B |
|---|---|---|
| smile | 0.0050 | 0.0038 |
| keystone | 0.0058 | 0.0022 |
| CRF | 1.213 | **1.029** |
| SRF | 2.0238 | 2.0235 (the floor) |
| EE | 0.824 | **0.987** |
| merit (chain) | 0.4934 | 0.3761 |

| e2e with R9 (centroid, roll 0) | record | Option B |
|---|---|---|
| smile | 0.0222 (0.40 µm) | **0.0122 (0.22 µm)** |
| keystone | 0.0058 | **0.0036** |
| CRF | 1.174 | **1.064** |
| SRF | 2.0248 | 2.0244 |
| EE | 0.849 | **0.955** |

The solved grating conic is small, K = −0.0038.  Freeing it lets the rest
move to a better basin: lens K 0.043 → 0.084, h⁴/h⁶ −0.48 / −6.0 →
−0.76 / −11.9, face offset 0.51 → 1.19 mm, block dz −0.55 → −1.62 mm,
R_g 787.2 → 792.4 mm.  Clearance 13.6 mm.  Deck
`dyson5_optionB_3k.in`.

**Its tolerance ladder** (`dyson5_sens_3kB.{txt,mat}`; e2e through a join
template = the record's up to its Slit, then Option B's Dyson):

| sensitivity (per 10 µm / 10 µrad / 1e-4) | Option A (record) | Option B |
|---|---|---|
| keystone, lens decenter x | 0.0099 | 0.0105 |
| keystone, grating decenter x | 0.0099 | 0.0105 |
| keystone, grating tilt y | 0.0074 | 0.0080 |
| keystone, grating clock | 0.0029 | 0.0042 |
| smile, lens / grating decenter y, grating tilt x | 0.0007–0.0009 | 0.0007–0.0009 |
| CRF, grating radius | +1.93 | +1.87 |
| CRF, lens radius / slit / detector defocus | +0.26 | +0.21 |
| CRF, lens thickness / air space | −0.11 / −0.17 | +0.001 / +0.012 |
| after detector refocus, every focus term | ≤ 0.006 | ≤ 0.006 |

The thickness and air-space rows differ because Option B's nominal IS at
its focus (a defocus error there is second order).  The record's nominal
sits 0.04 px of CRF off its best focus (refocusing it alone gives 1.175),
so those rows are first order for A.  Refocused, the two behave alike.

**Budget, paper A CBE keystone (0.057 px):**
- **Option B:** allowance 0.0534 → along-slit decenters **25 µm**, grating
  tilt y **33 µrad**, clock **64 µrad**.
- **Option A:** 26 µm, 35 µrad, 88 µrad.

**Reading.**  On this geometry Option B buys **nominal** margin:
- CRF −0.18 px alone and −0.11 e2e; EE +0.16 alone and +0.11 e2e; smile
  and keystone halved.
- That margin is what makes the paper's CBE CRF lower (2.04 vs 2.21 in
  Table 2).

It does **not** buy alignment insensitivity.  Keystone per µm or µrad is
6–45 % HIGHER, and the resulting tolerances are the same or a little
tighter.  The paper's "Option A is more sensitive" is not reproduced by an
alignment ladder alone.  Its CBE also rolls up surface fabrication and
thermal errors, which this ladder does not model.  So the claim is untested
here, not refuted.  Keep Option B as the spectrometer of record candidate,
for its nominal margin.
