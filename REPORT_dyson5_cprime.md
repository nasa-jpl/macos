# REPORT — dyson5 (c′): the 3k telescope from the SBG VSWIR seed

TO, 2026-10-07.  Brief: `BRIEF_to_dyson5.md` addendum 48 (with Dave's framing: a
telescope design TEMPLATE, `mmacos/templates/10_telescopes/tma_longslit/`).
Status: **steps 1–4 done.**  Step 2: conics + aspheres leave the strip far above spec, reason measured.  Step 3: the reweighted pole-frame freeform (R7) puts the telescope at the pixel floor at every field.  Step 4: e2e on R4–R7.  **Deck of record = R7** (was provisionally R4).  Open: smile 0.14 px.

## Step 1 — first order from the seed

### The model, and the one fact it forces

Unfolded thin-mirror chain along the chief: Joe's spec f 330 mm, D 183 mm, the
Fig. 4b legs scaled by 330/345, i.e. **258.3 / 245.8 / 299.4 mm**.  Three
conditions (EFL 330, focus at the slit, chief exit slope 0) on three powers.
The legs are the same in the fold (tangential) and strip (sagittal) sections,
so both sections need the SAME unfolded power per mirror.  A mirror tilted by i
gives φ_t = 2/(R_t cos i) and φ_s = 2 cos i / R_s, so each mirror must carry
**local radii at the chief with R_t/R_s = 1/cos² i** (1.33 / 1.57 / 1.07 at
30 / 37 / 15°).

- **A tilted sphere cannot do it.**  Solving all six conditions over the three
  radii and three AOIs collapses to AOI = 0.  The first order is anamorphic by
  construction.
- **An off-axis conic section can.**  Cut a conic of revolution at height h,
  with angle θ between the local normal and the parent axis.  There
  R_s = h/sin θ, R_t = R_s³/R², R_s² = R² − K h².  With θ = AOI every mirror
  comes out K = −1: three off-axis **paraboloids** give the exact first order.
  θ (and with it K) is left as a third-order DOF for step 2.

### The first-order table (`tls_first_order.txt`)

| stop | mirror | f (mm) | R_t (mm) | R_s (mm) | beam (mm) | footprint x × y (mm) |
|---|---|---:|---:|---:|---:|---|
| **M2 (record)** | M1 | +1012.5 | +2338.3 | +1753.7 | 183.0 | 239.9 × 211.3 |
| | M2 | −452.0 | −1132.0 | −722.0 | 136.3 | 136.3 × 170.7 |
| | M3 | +245.8 | +509.0 | +474.9 | 166.0 | 220.2 × 171.9 |
| M1 | M1 / M2 / M3 | −567.7 / +364.7 / +1130.8 | | | 183 / 266 / 166 | M2 309 × 333 |
| ½ way M1→M2 | M1 / M2 / M3 | −3024 / +1028 / +393.5 | | | 183 / 199 / 166 | M2 219 × 249 |

(Footprint x = across the strip at the ±4.7° edge field; y = on the surface in
the fold plane.)

- **Telecentricity is available with the stop at M2**, in closed form: M2 sits
  at M3's front focus, f3 = M2→M3 = 245.8 mm.
- **There is no intermediate focus.**  M1 is weak concave, M2 convex, and M3
  carries the power; the marginal ray crosses the axis only at the slit.  This
  is consistent with addendum 37 (a real M2–M3 focus and a pupil at infinity are
  incompatible): this family does not have that focus.
- **The entrance pupil is virtual**, 346.7 mm downstream of M1 (the engine
  computes the same point: StopPos z = 346.69 mm).  The strip angle at M2 is
  ×1.34.
- **A stop at M1 (or halfway to M2) is also telecentric, but it flips the form**
  (convex M1, concave M2) and M2 grows to ~310–330 mm (M1 stop) or ~220–250 mm
  (halfway).  The M2 stop is the small-M2 form of the paper.
- **The drawing checks out.**  The paraxial footprints match the Fig. 4b mirror
  chords (M1 ~240, M2 ~140, M3 ~205 mm).  The "~85 mm" M2 footprint in
  addendum 48 does not hold: the beam at M2 is 136 mm.

### The section in the engine (`tls_section.in`, `tls_section.txt`)

Three off-axis paraboloid sections: VptElt = parent vertex, RptElt = the chief
pole, psiElt = the parent axis.  The stop is an ELEMENT stop on M2, written as
the element's `ApStop= dx dy` = pole − vertex in M2's TElt.  The fold is in the
y-z plane and the strip runs along global x.  Engine, model 256, 41-pt grid
(1185 rays), 9 fields across ±4.7°:

| quantity | measured | spec |
|---|---|---|
| chief to slit normal, max over strip | **0.093°** (edge; 0.023° at ±2.35°) | < 0.5° |
| plate scale, local / to the edge | **329.99 / 329.87 mm** | 330 |
| image of the strip | 54.24 mm | 54 mm slit |
| admitted at the slit | **1.000** every field | |
| chief miss of the M2 pole | < 2e-13 m | 0 |
| working distance M3 → slit | 299.4 mm | > 250 mm |
| footprints (all fields) | M1 234.6 × 201.3, **M2 131.7 × 161.4**, M3 212.0 × 163.5 mm | |
| clearance (footprint + 5 mm mount; every leg vs every body) | **CLEAR, min +10.2 mm** (M3→slit leg vs M2); then input vs M2 +61.9 | > 0 |
| working F/# x / y from the extreme rays | 1.92–1.93 / 1.955 | [1.7, 1.8] |
| spot rms, as placed / best focus | 1.88–2.23 / 1.78–2.06 mm | pixel 18 µm |
| bow of the strip image (across the slit) | **445 µm** (centre vs the end-field line) | |

How to read these:

- **First order is met:** telecentric, plate scale, admitted, stop, packaging.
- **The image is not.**  Spots of ~2 mm rms on the uncorrected paraboloid
  sections.  The cone F/# read from the extreme rays carries that aberration,
  so it cannot be scored against [1.7, 1.8] until step 2: paraxially the cone
  is F/1.803 in both sections by construction.  The 445 µm bow (25 px across
  the slit at the strip ends) is the plane-symmetric distortion of the tilted
  train, also a figure-stage item.  It becomes a row in step 2, because a
  straight slit needs a straight strip image.
- The nearest clearance, +10.2 mm, is the final leg passing M2.  The paper's
  Fig. 4b has the same pinch.

Layout render: `tls_section_layout.png` (fold plane).

### Engine finding (for CC; no engine change made)

`macos_api_mod.F90:stop_info_set` refuses a stop element `iElt >= nElt-2`
(line 6333).  Every four-element telescope (M1, M2, M3, focal plane) therefore
cannot set `macos.stop(2)`.  The CLI's STOP has no such limit, and the deck's
element `ApStop=` applied at load works: StopPos is computed and the chief hits
the M2 pole to 2e-13 m.  The template works around it by loading one deck copy
per field.  The guard looks like a leftover; its comment reads
"0 < iElt < nElt-2".

## Step 2 — conics + aspheres with the cone bounded (in progress)

**The solve** (`tls_figure`) is an engine-traced Levenberg–Marquardt solve
with Jacobian scaling (the tGM settings, residuals in µm).  It varies the
GENERATOR, not deck text: every evaluation rebuilds the section, so the chief
always runs through the three poles and the stop sits exactly at the M2 pole.

Rows per solve field (0 … 4.7°, 5 fields; the section is mirror-symmetric
about the fold plane):

| row | what it measures | weight |
|---|---|---|
| SPOT | every ray, about the as-placed centroid on the slit | × sqrt(254/N) |
| PLATE | chief vs f·tan θ | × 30 |
| PLATE_Y | along-track focal length, from a 0.05° probe at the centre and edge | × 30 |
| BOW | the strict chief intercept across the slit, vs the centre field | × 30 |
| TELE | chief angle to the slit normal | 1e4 µm/rad |
| CONE | hinge outside F/1.7–1.8 (±0.005), per axis, on the extreme ray angles | 1e4 µm per F/# |
| WD | working distance (leg 3 + focus) | hinge below 250 mm |

**Why the cone rows are not the engine's `OptBeamSize=`:** that row is twice
the largest ray distance from the bundle centroid at one element
(`utilsub.F GetBeamSizeCmd`), a single isotropic footprint extent.  At the
slit it measures the spot, not the cone, and it cannot separate the two axes.

**What the first three tries established** (records kept as
`runs_figure_try*.log`):

1. **Free R_t/R_s with soft rows lose the first order.**  With the local radii
   as DOFs and only slit-axis plate rows, R1 sold the fold-plane first order
   for spots: F/2.37 across the slit, M2's R_t −1132 → −1858 mm, the slit
   moved 20 mm.  Adding an along-track plate row and a 10× cone wall was not
   enough either: the next R1 put the slit 117 mm in (the WD row checked the
   leg, not leg + focus; fixed) and the edge chiefs 2.2° off.  About 12k spot
   rows outweigh any soft first-order row.
2. **So the first order is re-derived, never penalized** (the afocal4
   doctrine).  `X.closure`: at every iterate R_t and R_s are re-solved from
   the first order of the current legs and AOIs.  EFL, back focus and
   telecentricity then hold exactly, paraxially.  Closure R1 (θ + focus):
   plate scale 329.96 / 329.98 mm (slit / across), chiefs ≤ 0.050°,
   WD 299.4 mm, but spots 138–341 µm rms.
3. **The R1 residual is ordinary third order** (Zernikes on the engine's
   rays about best focus):

   | field | rms OPD | dominant terms | after 11 terms |
   |---|---:|---|---:|
   | centre | 7.2 µm | spherical Z11 −9.9 µm, balanced by focus; coma/astig ~0 | 0.29 µm |
   | edge (4.7°) | 21.2 µm | field coma along the strip +31.7 µm; astig −15.6 / +13.7 µm; trefoil 9.3 µm | 0.36 µm |

4. **Even aspheres are the wrong basis on these sections.**  R2 (θ + focus +
   h⁴/h⁶ on all three) moved the cost 3 % (1.279e8 → 1.237e8) and the centre
   spot 153 → 138 µm; the edge did not move; the coefficients ran to 63 mm
   of "sag at the lit radius" on M1.  The poles sit 0.9–1.5 m from their
   parent axes, so an h⁴ term about the parent axis is almost entirely tilt
   and curvature at the pole.  The closed-form construction absorbs that into
   R and K, leaving almost no local 3rd/4th-order lever.  The rotationally
   symmetric asphere is a basis for sections cut near their parent axis;
   these are not.  **Step 2 therefore leaves the strip far above spec**
   (FWHM 7–16 px), and the brief's condition for step 3 holds.

5. **An unbounded layout abandons the form.**  The geom rung, with the AOIs
   and legs free (try 4, from the conic rung), dropped the AOIs from
   30 / 37 / 15° to 2.8 / 2.9 / 1.7°.  That is a near-coaxial train: tilt
   aberrations gone, spots 62–154 µm, chiefs ≤ 0.013°, but **clearance
   −70 mm** (M1's body in the M2→M3 beam) and the WD hinge on its 250 mm
   bound.  The solve has no clearance row, so nothing kept it in the
   unobscured family.  The layout DOFs are now bounded about the seed (AOI
   ±5°, legs ±15 %: `P.aoi_span_deg`, `P.leg_span`), keeping the form the
   template is for.  Clearance stays a per-rung CHECK, not a row; a rung that
   fails it is a failed rung.
6. **Bounded, the layout still does not pay** (try 5).
   - AOIs 27.6 / 32 / 10°: M2 and M3 ran to their lower bounds, the solve
     still wanting to un-tilt.
   - Legs 234.9 / 238.5 / 254.5 mm.
   - Spots 153–341 → 141–254 µm rms.
   - Clearance **−4.9 mm** (the M3→slit leg into M2).

   A failed rung.  The seed layout, the paper's, stays.  The freeform rungs
   branch from the conic rung with the layout frozen; that is where the
   paper's own freedom is.

## Step 3 — freeform (pole-frame polynomial, in progress)

The freeform is a **pole-frame polynomial** on the conic base:
`Surface= Monomial`, `MonCoef` x^j·y^(i−j) in the section frame about the
pole, normalized by the footprint half-size.  It keeps degree 3–6 and even j
(mirror symmetry about the fold plane), 12 terms per mirror.

- **Degree ≥ 3 adds nothing to value, slope or curvature at the pole**, so the
  first-order closure is untouched.  This is gated: random 50 µm terms on all
  three mirrors change the spot (1876 → 1370 µm rms), the chief stays on
  every pole to 2e-16 m and leaves exactly along the slit normal
  (`tTmaLongslit/test_the_pole_frame_freeform_is_first_order_neutral`).
- **A Zernike departure would not do this.**  Coma carries a tilt and
  spherical a curvature at the pole, so it would bend the chief and break the
  closure.

### The rung table (engine, model 256, 41-pt grid, 9 strip fields; `tls_figure.txt`)

FWHM per the record's scorer (rays ⊗ 1 px ⊗ Airy at 2.5 µm, F/1.8), in
18 µm pixels.  "x" = along the slit (the CRF direction); "y" = across it
(the ARF direction).  Bow = the strict chief intercept across the slit;
cbow = the spot centroid's.

| rung | DOFs | rms µm ctr / edge | FWHM x worst | FWHM y worst | chief max | F/# x | F/# y | bow / cbow max µm | x err edge µm | M2 foot mm | clear mm |
|---|---|---|---|---|---|---|---|---|---|---|---|
| R0 seed | — | 1876 / 2230 | 15.3 | 14.5 | 0.093° | 1.92–1.93 | 1.95–1.96 | 445 / — | −10 | 132 × 161 | +10.2 |
| R1 conic | θ + focus | 153 / 341 | 15.5 | 15.3 | 0.050° | 1.91–1.92 | 1.91 | 185 / — | −50 | 131 × 161 | +15.2 |
| R2 asph | + h⁴ h⁶ | 139 / 338 | 15.8 | 13.8 | 0.050° | 1.91–1.92 | 1.90 | 187 / — | −51 | 131 × 161 | +16.0 |
| R3 geom | θ + focus + AOI + legs | 141 / 254 | 11.6 | 14.7 | 0.037° | 1.91 | 1.91 | 7 / — | −44 | 128 × 148 | **−4.9** |
| R4 ff34 | θ + focus + pole poly deg 3–4 | 28 / 59 | **2.56** | 2.56 | 0.034° | 1.89 | 1.82 | 2.7 / 9.1 | −82 | 132 × 165 | +13.1 |
| R5 ff | θ + focus + pole poly deg 3–6 | **4.1** / 51 | 6.68 | 2.14 | 0.031° | 1.89–1.90 | 1.75 | 2.2 / 0.6 | −51 | 132 × 168 | +10.5 |
| R6 ff2 | R5 + CBOW rows, 3000 eval | 4.5 / 51 | 6.76 | 2.00 | 0.031° | 1.89–1.90 | 1.78–1.79 | 2.3 / **0.3** | −51 | 131 × 166 | +15.5 |
| **R7 ffw** | from R4: deg 3–6, along-slit spot rows ×3, outer fields ×2/×3, CBOW, 3000 eval | **1.8 / 2.3** (max 2.5) | **1.02** | **1.02** | 0.033° | 1.889–1.894 | 1.797–1.806 | 3.2 / 0.75 | −79 | 132 × 166 | +12.2 |

First order holds on every rung by construction (plate 329.96–329.99 mm,
along-track 330.0–330.1 mm, WD 299.4 mm, every ray admitted).  No rung has a
ray below F/1.7.  The cone sits at F/1.75–1.95 because Joe's D 183 / f 330 is
itself F/1.803; "a little faster than the spectrometer" needs a larger
entrance beam (the paper's 192 mm at f 345).

### End to end on R4 and R5 (`tls_e2e_R4.txt`, `tls_e2e_R5.txt`)

The (c) column is the record re-scored after the grating-aperture fix (macos 29f41da).
Joined to the 3k Dyson of record (CaF2 240, `size:F:240`) by the dyson5
engine join (`dyson5_t5f`, unchanged).  7 fields × 7 wavelengths, the grating
the stop, both rolls:

| e2e | Dyson alone (record) | R4 roll 0 / 180 | R5 roll 0 / 180 | R6 roll 0 / 180 | **R7 roll 0 / 180** | (c) record | Joe | paper |
|---|---|---|---|---|---|---|---|---|
| smile px | 0.005 | 0.732 / 0.732 | **0.090 / 0.093** | 0.149 / 0.147 | 0.140 / 0.141 | 1.64 | < 0.1 | < 0.05 |
| keystone px | 0.006 | 0.009 / 0.019 | 0.031 / 0.038 | 0.037 / 0.043 | **0.006 / 0.009** | 0.05 | < 0.1 | < 0.1 |
| CRF px | 1.21 | 2.51 / 2.50 | 6.98 / 6.93 | 7.16 / 7.22 | **1.175 / 1.255** | 4.02 | < 1.5 | < 2.8 |
| SRF px | 2.02 | 3.90 / 3.88 | 2.52 / 2.55 | 2.48 / 2.50 | **2.025 / 2.025** | 4.99 | < 1.5–2.0 | < 1.8 |
| ARF px (telescope FWHM across slit) | — | 2.56 | 2.14 | 2.01 | **1.02** | — | — | < 2.8 |
| energy in a pixel (min) | 0.82 | 0.03 | — | — | **0.853 / 0.810** | — | > 0.75 | |
| joined clearance mm | +0.54 | +0.6 | +0.6 | +0.6 | +0.6 | +0.06 | > 0 | |

- **The e2e smile is the telescope's CENTROID bow, not its chief bow.**
  R4's chief bow is 2.7 µm, but its centroid bows 9.1 µm (0.50 px): coma
  moves the centroid off the chief.  The smile read 0.73 px.  R5's centroid
  bow is 0.6 µm, and the smile is 0.09 px, passing Joe.  The spectrometer
  scores centroids, so the merit now carries a centroid-bow row (CBOW)
  beside CC's chief-intercept row.
- **CRF and SRF are the telescope's spots.**  R4's 2.6 px slit-axis spots
  take CRF from the Dyson's 1.21 to 2.5.  R5's 6.7 px slit-axis spot at
  3.5° makes CRF 7.0, while its small across-slit spots give the best SRF
  (2.5 vs the Dyson's 2.0).
- **The joined clearance is +0.6 mm:** the telescope's final leg against the
  Dyson block face, right at the slit.  This is the same pinch the record
  carries (Dyson alone +0.54, 3k (c) +0.06).  Every other telescope/Dyson
  pair clears by ≥ 13 mm, and M2 clears the final leg by 17 mm.  It is
  measured by `tls_clearance_joined` (footprint bodies from engine rays on
  the joined deck).  t5f's live clearance cannot take these sections: it
  lifts bodies onto the parent's base sphere, which the poles lie beyond.

**R6 (ff2: R5 + the CBOW rows, 3000 evaluations) stalled.**  Cost
8.4385e6 → 8.4218e6 (0.2 %).  The off-axis angles wandered at nearly
constant cost (M1 22.9 → 12.3°, K −1.66 → −5.50), so this is a flat valley.

- **The centroid bow halved but the e2e smile rose.**  0.6 → 0.3 µm in
  centroid bow, 0.09 → 0.15 px in smile.  The slit-plane centroid bow at
  fixed field angles is not the whole smile below ~0.1 px.  One difference:
  t5f launches each field on the sky direction whose chief lands on the slit
  line.  NOT yet separated.  R4 → R5 (9.1 → 0.6 µm; 0.73 → 0.09 px) still
  shows it is the large-scale driver.
- **CRF is the wall on every freeform rung past R4.**  The slit-axis spot is
  5–7 px FWHM at 2.35–4.7°.  The merit's spot rows are rms about the
  centroid, equally weighted in x and y, and R5/R6 spent them on the centre
  (4 µm rms, EiP 1.00).  More budget on the same merit does not move it.
  The next lever is the merit, not the budget.

### R7: the reweighted rung (CC 2026-10-07, from R4)

The merit change was the lever, not the budget.  Starting from R4 (the
CRF-best rung), the along-slit spot rows were weighted ×3: that is the CRF
direction, which the Dyson cannot fix, while the slit itself truncates the
across-slit width.  The outer solve fields were weighted 1 / 1 / 2 / 3 / 3,
and the CBOW rows kept.  3000 evaluations took the cost 6.66e7 → 1.08e7;
the remainder is plate + cone, spot only 1.5e5.

**The telescope then sits at the pixel floor at every field:**
- 1.8–2.5 µm rms, FWHM 1.02 × 1.02 px (the 1 px ⊗ Airy floor), EiP 1.00 —
  the paper's 2–3 µm class;
- chiefs ≤ 0.033°;
- cone F/1.889–1.894 along the slit, F/1.797–1.806 across, no ray below
  F/1.7;
- M2 131.7 × 165.9 mm, clearance +12.2 mm.

**The freeform is not a different mirror.**  Its sag departure over the
normalized disc is 0.46 / 0.55 / 1.11 mm P-V (M1 / M2 / M3), the same class
as R4's (0.24 / 0.57 / 0.85) and below the record's freeform (1.3–2.1 mm).
The figure is overstated: the disc is larger than the lit ellipse.

**End to end, the telescope is no longer what limits CRF, SRF or keystone:**
- CRF 1.175 px (at the Dyson's 1.02 floor except the edge field, 1.17);
- SRF 2.025 px = the Dyson alone's 2.024;
- keystone 0.006 px;
- EiP 0.85.

**What remains:**
- **Smile 0.14 px (fails Joe's 0.1).**  It is not the slit-plane centroid
  bow (≤ 0.75 µm = 0.04 px); parked, per CC.
- **SRF at 2.025 is the Dyson's own.**  It fails the paper's 1.8 and sits
  on Joe's 2.0; the paper's aspheric Dyson lens + conic grating (Option B)
  is the spectrometer form for that.
- **The cross-track distortion of 79 µm** (0.24 %) at the strip end is a
  ground-mapping term.

## Step 4 — end to end

On R4–R7 above.  **Deck of record: R7** (`challenges/dyson5/dyson5_cprime_3k.in`,
e2e decks `dyson5_t5f_cprime_roll{000,180}_e2e.in`; R4 was provisional until
R7 reported).

## What the seed bought over the (c) section

The (c) section put marginal rays at F/1.19 into an F/1.8 spectrometer, with
M2 at 126 mrad of freeform.  End to end it reached CRF 4.02 / SRF 4.99 /
smile 1.64 px.

The SBG VSWIR zig-zag seed (Fig. 4b) is a different family: stop at M2,
telecentric in closed form, no intermediate focus, and the Dyson hanging in
a 300 mm working distance.  Two things were then held:
- the first order, exactly: re-derived every iterate, never penalized;
- the two long-slit constraint row sets in every rung: TELE (chief angle
  to the slit normal) and CONE (a hinge on the marginal rays' F/# at the
  slit, along and across).

The figure came from a pole-frame freeform that cannot touch that first
order, with the merit pointed at the along-slit width and the outer fields.
The result:
- the telescope at the pixel floor at every field (FWHM 1.02 px, 1.8–2.5 µm
  rms);
- chiefs within 0.033°;
- the cone F/1.80–1.89 with no ray below F/1.7;
- M2 132 × 166 mm;
- e2e CRF 1.18, SRF 2.03, keystone 0.006 px.

CRF and SRF are the Dyson's own floors.  Smile, 0.14 px, is the one
requirement still open.
