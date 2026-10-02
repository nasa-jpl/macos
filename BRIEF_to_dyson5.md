# BRIEF for TO: dyson5 -- a VSWIR Dyson imaging-spectrometer challenge (2026-09-30)

From Dave via CC.  Your 2026-09-30 assessment is accepted with the corrections
and decisions below.  Run this lane on **Opus 5.5**; escalate to Fable only if
an engine gate fails against its closed form (that becomes an engine slice).

## Decisions (Dave)

1. **Our own prescription, from the spec.**  No proprietary sources: do NOT
   request Lori's CDR slides or any JPL-internal prescription.  Public papers
   (Carbon-I arXiv:2505.22545; the EMIT and CWIS design papers from Mouroulis
   and Green's group; Mouroulis & Green 2018 review) are fine for the FORM, the
   design conditions and the SPEC, never for surfaces.  dyson5 is therefore a
   design-from-spec challenge scored against the spec; say so in the record.
2. **The spec is Joe's, and it is EMIT-class:** F/1.8, 3000 x 500 px at 18 um,
   380-2500 nm, smile/keystone < 0.1 px (0.2 acceptable), SRF < 1.5-2.0 px
   FWHM, XRF < 1.5 px FWHM, plus the radiometric chain (throughput, grating
   efficiency, QE, slit loss) vs wavelength.  Record Jim's realism alongside:
   as-built SRF 2.5-3 px, 2-px slits, photon-limited, smile/keystone drive.
   The grating is the stop.  Note 3000 px x 18 um = a 54 mm slit -- larger than
   EMIT's 1280 -- so the block scale is part of the problem, not a detail.
3. **Talk is 1-4 months out.**  Tight beats fast.  Every number gated, every
   convention stated before its number, a runner from day one.

## Corrections to the assessment (verified in the engine 2026-09-30)

- "Physical-optics legs cannot pass a grating" is wrong as written.  The
  propagation chain's ray re-trace HAS a Grating branch (propsub.F ~740), so
  rays that seed and carry the grid pass a grating and its dL enters the grid
  phase.  What is absent is a WAVE-optics grating (orders, efficiency).  Keep
  grating diffraction analytic -- for that reason.
- The immersed grating math is right (`Snells_Law_Grating`: u = (na*i +
  m*lambda0/d*s)/nb, lambda0 vacuum), but how the engine assigns na/nb to a
  REFLECTOR embedded in glass is unpinned: no Rx in the corpus does it, and
  `IndRef(iElt)=CurIndRef` appears at two reflector sites (tracesub.F ~4108,
  ~4141).  Gate 1a must build the grating INSIDE a silica block and check the
  outgoing direction against sin(out) = sin(in) + m*lambda0/(n*d).  If it
  fails, STOP and report -- engine fixes are not this lane's.
- The stop wrapper's range is `0 < iElt < nElt-2` (macos_api_mod), i.e. at
  least the Return + FocalPlane tail after the grating.  Build the deck so.

## Build order (yours, kept), with gates and pre-registered nulls

1. **Engine gates** (mmacos tests, SUITE_FAST, both must have a must-PASS leg):
   - `tGratingImmersed`: reflection grating on a concave base inside a fused-
     silica block vs the closed form, 3 wavelengths, 1e-10 on the direction.
     Include the in-air case as the control (the existing pymacos grating
     fixture is air only, `IndRef= 1`).
   - `tGlassDispersion`: `GlassElt= Silica` (and CaF2 once added) at 3
     wavelengths via `set_src_wvl`, index read back vs Sellmeier, 1e-12; a
     trace-level check that the refracted angle FOLLOWS the wavelength (the
     null: a run where the index does not move is the failure to detect).
   - **CaF2 into the glass table** = an ENGINE change (macos repo:
     `macos_f90/macos_glass_list.txt` -> `tools/gen_glass_builtin.py` ->
     `glass_builtin.f90`, rebuild `makems.sh release gfortran`, then mmacos
     `make` after removing the mex).  Engine first, resources second.
   - **DO NOT rebuild libsmacos or the mmacos mex before
     `templates/40_benches/tg_psi_dm96_oap/runs/fix2x2.done` exists** (the
     descent 2x2 is running until ~17:00 today; each cell starts a fresh
     MATLAB and would pick up a different engine).  Until then, gates that
     need only Silica can proceed on the current build.
2. **`dyson_layout` closed-form seed** in `mmacos/design/src/` beside
   `offner_layout.m`: concentric Dyson, grating radius R_g = n*R_b/(n-1) with
   the slit/FPA plane at the block's flat face through the common centre,
   slit and FPA offset either side of the axis.  VERIFY the condition against
   Mouroulis & Green 2018 before pinning it, and record where the classical
   seed breaks at F/1.8 with a 54 mm slit (block size, offsets, the field
   flattener or meniscus the real ones carry).  A pre-registered null: if a
   single-block Dyson cannot meet F/1.8 at this slit length, report the
   scaling law first; do not add elements until that is on record.
   `SPECTROMETER_DESIGN_REFERENCE.md` in `optical_design/` follows the
   telescope reference's pattern -- and goes on the dev-only strip list
   (`DEV_FILES.md` + `release-exclude.txt`) like its sibling.
3. **`spectrometer_score`**: the two 2D maps (field-angle map, wavelength map)
   from `set_src_fov` x `set_src_wvl` sweeps and `macos.spot` centroids; smile
   (centroid drift along the slit at fixed lambda) and keystone (across lambda
   at fixed slit position) in PIXELS with the sign and axis convention stated
   first; SRF/XRF as the convolution chain geometric spot (x) 2-px slit image
   (x) pixel (x) Airy (analytic; ~0.3 px at 2500 nm, F/1.8), FWHM in pixels;
   the radiometric chain closed-form with the SLIT loss the one measured term
   (a physical-optics leg to the grating plane, what overfills its aperture).
   One scorer for Dyson and Offner; the Offner comes from `offner_layout`.
4. **Native optimize** over grating radius, block radius, offsets, grating
   period: the CALIB multi-field / multi-wavelength path (<= 12 FOV x 6
   lambda), merit in pixel units.  Zernike-solve doctrine applies (memory
   `project_zern_solve_doctrine`): power pinned multi-field, tilt verified.
5. **Record + slides**: `mmacos/challenges/dyson5/` in the rodgers3 pattern
   (target, the extension, one open question for the room), the Dyson-vs-
   Offner trade slide, and a live-demo candidate (five-element geometric
   deck, seconds per solve).

## Standing rules (DAVE'S, not optional)

- ONE user-editable runner from day one: `dyson5_params.m` + `dyson5_run.m`,
  every stage run THROUGH it, README "run it yourself".
- Real examples: save the `.in` and the `.mat`, never `exit(0)` in an example.
- State the convention before the number (challenges README).
- Corpus/engine facts from the ENGINE, never by text-parsing `.in` files.
- Gates fail closed: every gate ships a must-PASS leg.
- New guards warn once per run; hard errors need Dave's sign-off.
- Never edit an executing script; ONE model-1024 MATLAB on this box (yours
  are small -- keep them small).
- Commit on `dev-candidate` (both repos as needed); state SHA + branch; do
  not push -- Dave reviews.
- Write every resolved oddity into the record at resolution time.

## Deliverable of the first beat

Gates 1a/1b green (or a STOP report), `dyson_layout` with the condition
verified and the F/1.8 / 54 mm scaling on record, and the runner skeleton.
Report as a BRIEF-style note in `mmacos/challenges/dyson5/`.

## Addendum (Dave, 2026-09-30): back every ray metric with a PROPAGATION run

Every analytic / ray-trace metric (smile, keystone, SRF, XRF, the slit and
grating diffraction losses) gets a physical-optics twin: seed the field at the
slit, propagate leg by leg with the near-field (`PropType` 4/5, plane-to-plane
and sphere-to-sphere), the DFT legs (13-15, which let you set the OUTPUT
sampling -- a few um over a few-hundred-um window at the FPA per (field,
lambda) point, never the whole 54 x 9 mm FPA) and the far-field leg, and read
the COMPLEX FIELD back (`macos.complex_field` / `cfield_get`, re + im), not
just the intensity.  Phase is what carries the grating: the ray OPD is an
unwrapped LENGTH and does not wrap, but as a phase it spans thousands of waves
across the grating, so any comparison between the ray and wave pictures is
done on the complex field (or on intensity-derived quantities), never on a
wrapped phase map.  Metrics from the wave twin: PSF centroid per (field,
lambda) -> smile/keystone maps to compare with the ray centroids (expect
agreement well under 0.01 px where the PSF is symmetric; where they differ,
that IS the finding); SRF/XRF from the propagated PSF (x) 2-px slit (x) pixel,
replacing the analytic Airy term; slit loss and grating-aperture overfill from
the field at the grating plane.

**Engine gap, verified 2026-09-30, CC's lane, NOT yours:** every propagation
kernel in `propsub.F` (NFPROP, PPPROP, SFPROP, FRPROP, NFPropDFT, FFPropDFT,
SPH2PL/PL2SPH, FFPROP, and `FnCalc`) is handed `WaveBU`, the VACUUM
wavelength, with no division by the medium's index.  The inter-leg geometric
phase is right (`CumRayL` accumulates `CurIndRef*RayL`, an optical path, and
`TPL = 2pi/WaveBU`), but a diffraction leg INSIDE the silica block runs at the
wrong Fresnel number by n ~ 1.45.  In a Dyson the slit and the FPA sit at the
block's flat face, so essentially every leg is inside glass -- the wave twin
of the Dyson is blocked until the kernels take `WaveBU/n_leg` (n_leg = the
index of the medium the leg traverses, `CurIndRef` at the leg's start
element).  CC does that slice after the descent 2x2 releases the build
(`fix2x2.done`), gated by the Fresnel identity: a leg of length z in glass of
index n must equal the same leg of length z/n in vacuum, bit-for-bit, and the
pre-fix engine must FAIL it by n.

**What you do meanwhile:** build the wave-twin harness on the OFFNER first --
it is all-reflective, entirely in air, so the current engine is correct for
it today.  Same scorer, same (field, lambda) sweep, same complex-field
readback.  When the medium-aware kernels land, the Dyson runs through the
identical harness.  Two conventions to state before any number: the field's
sampling comes from `macos.dx_at` / `elt_dx_get` at the element (SI, via the
base-unit factor -- memory `project_pymacos_dx_unit_convention`), and the
grating's m-th order is the ONLY order the engine carries (the ray model gives
one diffracted direction), so grating efficiency stays in the radiometric
chain.

## Addendum 2 (2026-09-30): Mouroulis & Green 2018 is now on disk
Dave placed the PDF in `mmacos/challenges/dyson5/` (git-ignored: SPIE
copyright, local only -- cite, never commit).  CC's digest is
`NOTE_mg2018_digest.md` beside it: the Sec. 5.3 design principles (distortion
~1% px at design, >75% ensquared, degraded spots OK for uniformity, grating =
stop), the SRF/CRF/ARF definitions verbatim in form for the scorer, the
incoherent-approximation rule (Airy diameter < pixel and slit -- met at
F/1.8), the Fig. 15 CaF2 Dyson spec table as the flight-class column, and
Table 3 placing Joe's spec in the ALIS regime (3200 px, 380-2500, 7 nm).  The
paper restates neither the concentric condition nor a prescription; beat 1's
verification stands.

## Addendum 3 (2026-10-01): finding #3 fixed; OPD across a grating is MODULO LAMBDA (Dave)
`dL = Order*lambda/RuleWidth * dot(s0, rho)` landed (macos, CC); `tGratingOpl`
is the gate.  Dave's rule for every comparison the wave twin makes from here:
the OPD of a wavefront that has crossed a grating is a STAIRCASE (the groove
structure), i.e. defined only modulo lambda; the engine's ray OPL carries the
smooth order-m phase function, which equals the physical wavefront mod
lambda.  So: compare pupil OPD, ray-vs-field phase, and chain-vs-engine phase
in PHASE SPACE -- wrap to (-lambda/2, lambda/2], or compare the complex
fields / exp(i 2 pi W/lambda) -- never as unwrapped lengths across the
grating.  `tGratingOpl`'s rms on the reference sphere is fine as written only
because the fixed engine's OPD there is ~0 (no wraps); state the rule in the
record and wrap before any rms where the OPD could exceed lambda/2.  Re-run
s2w unchanged once the mex is relinked; the order -1 rows then replace the
"defect" rows.

## Addendum 4 (2026-10-01, Dave): the COMPACT variant is the path; and the deck starts now

**Decision.** Pursue the paper's compact Dyson variant (a separate mirror plus
a meniscus) as rung R4 -- not the free-radius concentric block.  The beat-3
ladder settled the fixed-scale question: R3's de-concentred block solves the
distortion (keystone 0.011 px), the face asphere is inert, and the blur is a
size statement (h^4/r^3 at the ~32 mm effective field).  Score R4 under the
IDENTICAL operands and the identical scorer, and put R3 and R4 in ONE trade
table (keystone, smile, CRF FWHM, ensquared, SRF, length, footprint, element
count, glass volume) so the trade is a single comparison.  Keep the free-radius
record as the "size alone" column of that table.

**Re-score first.** `dyson5_s2_maps.png` and every SRF/CRF number quoted in
beats 2-3 predate the chord-ruled engine (the 2.1-3.1 px SRF ramp in the s2
maps IS finding #2).  Re-run s2 and s3 on the current mex before any number
reaches a slide; beat 2c's s2w order -1 rows likewise (finding #3 is in).

**The deck.** Dave is starting `demo_session/deck_dyson.md` -> `deck_dyson.pptx`
(the gauge-deck toolchain; CC builds it; DECK_STYLE.md governs: lean main
path, one Backup divider, every result slide pairs its LAYOUT figure with its
PERFORMANCE map, figures are the tools' own output unmodified).  What the
dyson5 runner must therefore emit, stably named, one set per stage, from a
committed producer (no hand-rendered figures -- see memory
`feedback_figure_producers_rot`):
1. **Layout figures to deck standard** -- the current `dyson5_s1_layout.png`
   (whole spheres and axes in metres) is a chain sketch, not a layout.  Needed:
   a dispersion-plane section and a slit-direction section per form, drawing
   the block outline (flat face + spherical face), the grating arc, the air
   gap, the slit and FPA as marks, the meniscus/mirror where present, rays at
   380 / 1440 / 2500 nm from slit centre and both ends, axes in mm, a scale
   bar.  `spectrometer_layout_fig.m` (design/src), called by every stage that
   emits a deck: `dyson5_<stage>_layout_<form>.png`.
2. **Performance maps** per scored deck: the 2 x 4 map panel as now
   (field-angle map, smile, SRF, ensquared) plus keystone and CRF panels, with
   the convention of each axis quantity in the axis label (pixel units; sign
   and axis stated in the figure, not the caption).
3. **The ladder bars** as now, plus the R4 column, and the same chart for the
   free-radius record.
4. **The propagation twin** evidence: order-0 Airy / pupil OPD (validation),
   then per (field, lambda) PSF centroid vs ray centroid, SRF/CRF from the
   PSF, the sinc^2 slit-loss comparison -- each as a figure the deck can place.
Every figure's producer is in the tree and runs from `dyson5_run` stages;
the deck's `figs_dyson/` copies are taken from the runner's outputs by a
script (`crop_panels.py` pattern), never edited by hand.

## Addendum 5 (2026-10-01, Dave): layouts are ENGINE RENDERS of the .in file, not chain sketches
Dave: "I'd expect to see rendered optics, not abstract circles.  Where is the
.in file?"  The layout figure of record is `macos.view_rx` on the emitted
deck -- bodies, rays and labels read back from the engine.  Producer:
`challenges/dyson5/dyson5_view_figs.m` (CC, committed; two views per deck,
3-D and the dispersion plane); `demo_session/tile_views.py` composes them.
This REPLACES addendum 4's `spectrometer_layout_fig.m` ask: do not write a
separate drawing tool; call `dyson5_view_figs` from every stage that emits a
deck, so each ladder step and R4 has its render beside its score.
Two things the renders exposed in the emitter, for you to decide and record:
1. **The block does not render as glass.**  The slit sits ON the flat face,
   so the emitter makes the glass the SOURCE medium and only the spherical
   face is a Refractor; `view_rx` joins consecutive Refractor pairs into a
   solid, so the block shows as a thin face.  Emitting the flat face as a
   Refractor (slit in air at the face, cone F/1.8 in air refracting to the
   same in-glass half-angle) is more physical (the slit is bonded or
   air-spaced on real Dysons) and renders the block as a body.  Ray results
   are identical by construction; verify with tSpectrometerRx before
   switching.
2. **Labels pile up at the face** (E1/E4-E7 at z = 0): the Reference/Return/
   FocalPlane bookkeeping planes.  Name the elements in the emitter
   (`EltName=`) so the viewer's labels read "Slit", "BlockFace", "Grating",
   "FPA" instead of indices.

## Addendum 6 (2026-10-01, Dave): the layouts show OBSTRUCTED beams -- clearance is a number, apertures are declared, the Offner's slit moves
Dave, on the engine renders: "serious blockage of the beams; the Offner needs
separation to clear M2; the Dyson looks like the Offner."  Measured from the
engine's ray history on the emitted decks (`dyson5_clearance_probe.m`, one
field point at the slit's +6 mm, 1 um):

| deck | what crosses what | numbers |
|---|---|---|
| Offner `dyson5_s1_offner.in` | leg 1 (slit -> M1) at the grating's plane z = -250 mm | beam y in [-39.4, +51.4] mm; grating footprint y in [-44.6, +44.6] mm -> the grating body sits INSIDE the incoming beam |
| Offner | leg 3 (M3 -> FPA) at the same plane | beam y in [-57.5, +33.3] mm -> through the grating again |
| Offner | M1 footprint | radius 89.5 mm about y = +6; the "M1" and "M3" zones overlap (both at y ~ 0) |
| Dyson `dyson5_s3_r3.in` | nothing crosses a body | block face footprint radius 42 mm at z = 217; grating footprint radius 134 mm at z = 691 (a 270 mm wide grating of 0.70 m radius); slit at +6.0 mm, this wavelength lands at -11.95 mm on the face |

**Why the trace cannot see it:** every element has `ApType= None`, and the
sequential trace never tests a ray against an element it is not currently
traversing.  Obstruction is a CLEARANCE property, not a ray-trace property
(memory `feedback_clearance_gates`: a gate's margin is a NUMBER, not a body).

**Required, before any further rung is scored:**
1. **Declare apertures** on every element from the multi-field, multi-lambda
   footprint union plus a mount margin (ApVec from the hull, as
   `dmg_bench_clearance` does for the benches).  Then the engine vignettes
   what it should, `view_rx` draws real bodies (the Dyson block becomes a
   cylinder, not the beam's hull -- that is why it "looks like the Offner"),
   and ensquared energy is measured through the real stops.
2. **`spectrometer_clearance`:** for every leg, the ray hull of that leg
   against every element body it does NOT traverse (aperture + mount), the
   minimum clearance in mm per (leg, body), printed as a table and FAILING
   the stage when negative.  Run it from every stage that emits a deck.
3. **The Offner's slit offset.**  Clearing the grating needs the slit->M1
   and M3->FPA beams beside the grating: offset >~ grating half-width + beam
   half-width at that plane + margin ~ 45 + 45 + 10 = 100 mm, i.e. ~0.2 R
   (the classical Offner geometry, M1 used in two separated zones).  Re-pose
   the Offner at that offset and re-score; it is the all-reflective twin and
   the propagation twin's validation ground, so it must be a real layout.
   Note the Offner deck runs at F/2.8 (source cone 0.359 rad full); state
   that beside the Dyson's F/1.8 wherever the two are compared.
4. **The Dyson's slit/FPA mechanics.**  The slit mask and the FPA package
   share the block's face: slit at +6 mm, the band landing 12-21 mm away.
   Add a package model (active area 54 x 9 mm inside its carrier + window;
   a cold shield if the SWIR end needs one) and score the slit-to-package
   clearance; the review's answer when it fails is a fold prism / in-built
   reflector at the slit -- make that part of R4, not an afterthought.
5. **Element sizes in every record and on the deck:** block diameter and
   thickness, grating diameter and radius, mirror diameters, overall length;
   the trade table (addendum 4) gains those columns.
**Correction to addendum 5 item 1:** the emitter ALREADY carries the flat
face as `BlockFaceIn`/`BlockFaceOut` refractors; the block renders as the
beam's hull only because no aperture is declared.  Item 1 above fixes it.

## Addendum 7 (2026-10-01, Dave): R4 is the design of record; addendum 6 applies to IT; the telescope spec is EMIT's
**R4 accepted as the headline design** (keystone 0.0026 px, smile 0.0051,
CRF 1.33, EE 0.76 at 694 mm and 22.7 L -- the free-radius result at 63 % of
the length and 27 % of the glass).  Two things the record keeps beside it:
the solve ends on its bounds (4 mm plate, ~2 m radii) and the landscape is
multimodal -- beat 4's global search over the meniscus is where that
resolves; and the twin's diffraction ensquared energy (0.85 Offner, 0.45
Dyson) is quoted for which deck?  Run the twin ON R4 and report its EE with
diffraction next to the geometric 0.76; the deck will show the pair.

**Addendum 6 is the collisions brief -- do it on R4, not R3:** apertures from
the footprint hulls on every element (block, meniscus, grating, slit mask,
FPA), `spectrometer_clearance` per (leg, body) with the minimum clearance
in mm and a FAIL on negative, the Offner re-posed at ~0.2 R (it is the
twin's validation ground and the trade's reference, so it must be a real
layout), the slit / FPA package model with the fold-prism split when it
fails, and element sizes in the trade table.  The layout figure of record
stays the ENGINE render of the .in (`dyson5_view_figs.m`); your mm-scaled
chain sections are the 2-D complement, placed beside it, never instead.

**The telescope (Dave: "go with the EMIT parameters for now").**  Beat 5 =
the fore-optics that feed the slit, scored under the same runner.  Spec:
| item | value | from |
|---|---|---|
| orbit altitude | 420 km | EMIT (ISS) |
| ground sample distance | 60 m | EMIT |
| IFOV | 0.143 mrad | 60 m / 420 km |
| pixel | 18 um | Joe's FPA |
| focal length | 126 mm | 18 um / 0.143 mrad |
| F-number | 1.8 | the spectrometer's, matched at the slit |
| aperture | 70 mm | f / 1.8 |
| cross-track field | 24.6 deg full (3000 px) | 3000 x 0.143 mrad; swath 180 km |
| image | flat, telecentric, 54 mm slit, 2-px (36 um) slit width | the Dyson's entrance condition |
Form: the review's Sec. 5.1 (a wide-field reflective telescope with a
telecentric relay, Fig. 12; Carbon-I used a freeform TMA).  Start from the
design layer's TMA / fold builders (`macos.design.Telescope`, the fold and
centred-TMA tools -- memory `project_fold_extraction`), two-mirror or TMA
seeds, telecentric output as a merit term (chief-ray angle at the slit, in
mrad), spot size at the slit under the pixel, along-slit mapping linearity
(the telescope's share of smile/keystone), and the PUPIL MATCH: the
telescope's exit pupil must land on the spectrometer's stop (the grating).
Realizability rules stand (AOI < 15 deg, shroud fit, clearance as a number).
Deliverable of the first beat: the seed laid out and engine-rendered, its
spot and telecentricity maps over the 24.6 deg field, and the pupil-match
number -- then the spectrometer and telescope traced END TO END as one deck.

## Addendum 8 (2026-10-01): beat 4 order -- the centroid question first, with its discriminators
Beat 3c accepted (354186e): apertures in the engine's frame with the no-vignetting
gate, `spectrometer_clearance` failing on negatives, the Offner real at 0.22 R
with the classical corrections, the package model, R4's twin pair 0.747 / 0.759.

**The open finding is now the headline risk** and the deck says so on its title
slide: if the detector-seen keystone is the WAVE centroid, R4 reads up to 0.12 px
against Joe's 0.1, not 0.0026.  Settle it before any optimizer runs, and settle
it with discriminators, not a single confirmation:
1. **The amplitude-weighted ray centroid** from the engine's pupil amplitude
   (|WFElt| at the seed, or per-ray transmission), compared with the PSF
   centroid per (slit, lambda).  If it reproduces the 0.122 px, the mechanism
   is proven and the operand changes.
2. **The wavelength law is the tell.**  Amplitude weighting of the RAY
   aberration is achromatic to first order (ray spots do not scale with
   lambda; Fresnel transmission barely does), so a difference that GROWS with
   lambda points instead at a diffraction-scale effect: the PSF's width grows
   as lambda, so (a) truncation of asymmetric wings by the twin's window and
   (b) the sampling of a lambda-wide PSF on the grid both grow with lambda.
   Null test: double the twin's window and halve its sampling at 2500 nm --
   if the 0.122 px moves, it is the twin's numerics, not the design.
3. **Separate the two forms of weighting:** an unweighted PSF centroid over
   a window that holds > 99 % of the energy vs the same over the 1-px box;
   and the centroid of the engine's far-field leg vs the DFT leg at the
   same (slit, lambda).  Agreement between legs and windows isolates the
   physics from the propagator.
Report the result as measured, with the mechanism stated as proven or
refuted; the deck carries whichever keystone the detector sees.  Then: the
global meniscus search, the native optimize with that operand, R5's fold
prism under the clearance check (cold shield height as the parameter), and
beat 5's telescope at the EMIT parameters.

## Addendum 9 (2026-10-01): centroid closed; the slit-loss factor is next; deck asks
**Centroid (4d95a9a): accepted.**  Mechanism proposed, tested, refuted; numerics
ruled out by the window and pitch variants; cause found (the far-field grid's
frame, both signs on the Dyson, one on the Offner) and taken from the engine's
source frame.  That is how a finding closes.  R4's 0.0026 / 0.0051 px stand on
the deck as what the detector sees, with the diagnosis in Backup.

**Slit diffraction loss -- the factor is now the open item** (Dave asked why it
does not match sinc^2).  The record's model: a 36 um x 0 mm slit, far-field leg
to the grating plane, 255-point grid, window 11 % wider than the acceptance;
engine / closed form = 0.30 at 380 nm, 1.04 at 700, 1.32 at 1440, 1.36 at
2500.  A ratio that changes SIGN across the band is not a normalisation.  Three
tests, each a one-knob change, reported as a table before any explanation:
1. **Grid:** 255 -> 511 -> 1023 points at 380 and 2500 nm (the slit is ~61
   samples wide at 255; the far-field pitch is 0.44 to 2.9 mm).  Convergence or
   not is the first fact.
2. **Window:** acceptance-to-window 1.11 -> 2 -> 4.  Energy diffracted beyond
   the window aliases back INTO it on a periodic FFT grid and is counted as
   loss; the tail beyond the window scales as lambda, which is the right sign
   for the long-wavelength excess.
3. **The slit's length:** 0 mm -> 54 mm (2-D).  A zero-length slit spreads
   uniformly along x over the whole window; with the real length the x-pattern
   is a narrow sinc and the acceptance should be the F/1.8 CIRCLE, not the
   |y| strip -- check which the loss definition uses.
Also confirm the plotted "sinc^2 closed form" is the exact integral over the
acceptance (the record says so) and not the 1/(pi^2 u0) asymptote, which is
off at 2500 nm where the acceptance is only ~4 sidelobes wide.

**Deck asks that reach the runner:** (1) element tables are generated from
the .in files (`demo_session/rx_elt_table.py`: E#, name, type, surface, |R|,
conic, aperture radius, vertex z) -- keep `EltName=` meaningful on every
emitted deck, they are now on the slides; (2) a methods slide defines engine /
chain / twin / ladder / seed / operand in plain words -- use those words the
same way in the records; (3) `dyson5_s3_r3.in` on disk has no apertures (only
the seeds and R4 were re-emitted in 3c) -- re-emit every ladder deck with
apertures so the tables and renders agree across steps.

## Addendum 10 (2026-10-01, Dave): next steps approved; three future-work items on the record
Dave: "OK for TO next steps" -- the addendum 8/9 order stands: the slit-loss
factor tests, the global meniscus search, the native optimize with the
smile/keystone operands, R5's fold prism under the clearance gate (cold-shield
height as the parameter), then beat 5's telescope at the EMIT parameters.
Add when the design of record is stable (not before R5 and the telescope seed):
1. **A surface-by-surface tour of the prescription,** as the CTB record did
   leg by leg: for each surface in order, its role in a sentence, the ray
   footprint (centroid, size), the clearance to its neighbours, and where a
   propagation leg ends the field / intensity there; one strip of figures per
   surface from the runner, the engine render marking the surface.  It is how
   a reader who did not build the deck learns it.
2. **Spot diagrams in the Mouroulis & Green form** (their Fig. 16): a grid of
   spot diagrams, slit positions down, wavelengths across, each drawn inside
   the 18 um pixel box at a common scale, geometric from the engine's rays,
   with the propagated PSF beside it where the twin has run.  A producer in
   the tree, run from the scoring stage, stable names per deck.
3. **The telescope** (beat 5, addendum 7): the EMIT-parameter fore-optics
   feeding the slit, then the end-to-end deck.

## Addendum 11 (2026-10-01, Dave): the closure envelope -- "for which parameters should designs close?"
The run-it-yourself slide now states the EXERCISED envelope (F/1.8, 220 mm
silica, 54 mm slit, 380-2500 nm, 18 um, order 1; the scaling scan 50-500 mm
at order 0; the Offner at F/2.8) and says closure outside it is untested.
Make that a measured statement.  A sweep stage (`s4env` or similar, opt-in,
runs R4's solve from the design-of-record seed at each point, scores it,
records pass/fail against the spec with the failing metric named):
| axis | points |
|---|---|
| F-number | 1.6, 1.8, 2.0, 2.2, 2.8 |
| block radius | 150, 180, 220, 260, 300 mm |
| slit length | 30, 40, 54, 60 mm (pixel count follows) |
| pixel | 18, 30 um (the paper's) |
| glass | Silica, CaF2 |
One-axis-at-a-time from the design of record first (5 axes x ~5 points =
~25 solves, minutes each on the chain), then the two-axis corner that fails
first.  Output: a table + one figure (pass/fail map per axis pair), and the
sentence the slide needs: "designs close for F/x-y, blocks of r-s mm, slits
to t mm; the first metric to fail outside is ...".  Record bound-hitting
solves separately: a solve on its bounds is not a closed design.

## Addendum 12 (2026-10-01, Dave): polarization sensitivity, Dyson vs Offner -- future work
Mouroulis & Green credit the Dyson's low polarization sensitivity to its
near-normal incidence; quantify it against the Offner with the engine's
polarization machinery (memory `project_polarization_plan`; mmacos
`macos.jones_pupil`, `pol_maps`, `pol_zernike`; `coat_set` for coatings).
When the design of record is stable:
1. **Coatings as built:** AR on the silica faces and the meniscus, a metal
   (aluminium or silver) grating with a protective overcoat, aluminium
   mirrors on the Offner; thicknesses PHYSICAL (`coat_set`), overcoats at
   lambda/4 of a stated working wavelength (the quarter-wave trap is in the
   engine cheatsheet).
2. **Per (slit position, wavelength):** the Jones pupil of each form, its
   diattenuation and retardance (pupil mean and variation separately; the
   mean is a state change, the variation an aberration), and the one number
   the spectroscopist wants: the instrument's polarization sensitivity,
   (I_max - I_min)/(I_max + I_min) over input linear polarization angle, at
   the detector, per wavelength -- plus the Stokes-to-measured matrix if the
   remote-sensing group wants Mueller terms.
3. **The comparison figure:** sensitivity vs wavelength, Dyson and Offner on
   one axis, with the angle-of-incidence histogram of each form beside it
   (the mechanism in one picture).  Grating efficiency's s/p split is NOT in
   the engine (scalar order model) -- state it as the term the figure omits,
   or add it from a published efficiency curve as a closed-form factor.
Record both forms' numbers in the trade table as a column.

## Addendum 13 (2026-10-01, Dave): full diffraction from the slit through the grating -- future work, AFTER the telescope
The current scorer is the incoherent chain (slit (x) LSF (x) pixel), valid
here because the Airy disk (11 um at 2500 nm) is under the pixel and the slit;
Mouroulis & Green put the full treatment's correction at ~10 % of the
response functions.  Dave: quantify it, and wait for the telescope -- the
field on the slit IS the telescope's image, and the result depends on how the
telescope's cone fills the spectrometer's pupil.  When beat 5's telescope is
on record:
1. **The input:** for each scene point across the slit width (a few points
   suffice; the scene is spatially incoherent), the telescope's point image
   at the slit plane -- propagated, not the ray spot -- truncated by the slit
   mask (36 um, the mask as a hard aperture on the slit plane).
2. **The chain, coherent per source point:** slit plane -> block (medium-
   aware kernels) -> meniscus -> grating as the design-order phase (the
   engine carries one order; efficiency stays in the radiometric chain) ->
   back -> detector; one far-field or DFT leg per segment as the propagation
   twin already does.  Sum the intensities over the source points: that is
   the partially coherent slit image of a uniform scene.
3. **What to report, against the incoherent chain:** SRF and CRF FWHM per
   (slit position, wavelength), ensquared energy, the slit-truncation loss
   (the energy the slit mask removes and the energy the grating aperture
   loses -- the s2l stage's number generalised), and the pupil fill at the
   grating (over- or under-filled by the telescope's cone).  The paper's
   ~10 % is the expectation to test, not assume.
Precursor available now, if wanted before the telescope: a uniformly lit
slit (what s2l does) with the incoherent sum across its width -- a bound,
not the answer.

## Addendum 14 (2026-10-01): beat 4c's four engine items -- status
1. **CALIB SPOT-target derivative stride (the blocker): CONFIRMED and FIXED.**
   `design_optim.F funcs_app` advanced the derivative columns by `opd_size`
   (= mpts^2) where the value loop advances by `obj_size` (1 for SPOT,
   n_wf_zern for WFE_ZMODE, = opd_size for WFE -- which is why every
   Telescope optimize was fine).  Pinned on the bounds-checked CLI exactly as
   you reported (`YFIT` subscript 16385 of 30, line 785).  Both sites now
   `obj_size`.  Gate: the reproducer on the debug CLI + tSpectrometerRx's
   native leg (un-mark it when the engine of record carries the fix) +
   tDesignTelescope (the WFE path must be bit-unchanged).
2. **OptAsph slice (`smacos_compute.inc:382`): CONFIRMED and FIXED** --
   `ptbArr(ia+1:ia+n_optAsphArr(ie))`.
3. **`OptRayGrid=` heap corruption: UNDER INVESTIGATION** on the bounds-
   checked CLI with your seed deck + `OptRayGrid= 21 / 41` (the parse at
   msmacosio.inc:224 caps only at mpts; the default is nGridpts/2-1; the
   optimizer swaps `npts=opt_npts` at macos_cmd_loop.inc:582).  Will report
   the overrun site.  Until then: leave `OptRayGrid=` unset.
4. **A centroid-position operand for CALIB: QUEUED as an engine feature,**
   not a bug -- a per-(field, lambda) chief / centroid position target so
   smile and keystone become native operands.  Scoped after items 1-3 land;
   the chain's lsqnonlin path is the operand's home until then.

## Addendum 14, item 3 resolved (2026-10-01): `OptRayGrid=` was the blocker seen through another knob
With fix 1 on the engine (macos 0d257ff) your seed deck with `OptRayGrid= 21`
and `41` runs CALIB to completion through the mmacos mex and MATLAB exits
cleanly (2 s and 5 s), with no other change -- the crash-at-exit / hang /
first-iteration crash were the SPOT-stride write (row 16385 of a 30-row
array) landing on different live memory for each ray-grid size; the default
grid happened to put it somewhere benign.  The parse cap at mpts was never the
problem.  `OptRayGrid=` may be used again.  (A bounds-checked one-iteration
run of the 41-point deck is finishing as the "nothing else out of bounds"
statement; if it reports anything, that goes in a further addendum.)
Still true: item 4, the centroid-position operand, is a queued feature.
Closed: the bounds-checked CLI ran the 41-point deck through a full
optimization iteration ("Optimization iterations = 1") with no out-of-bounds
access; item 3 is item 1, nothing further.

## Addendum 15 (2026-10-01): the deck is at beats 4b-4e; two producer items; one engine note acknowledged
`demo_session/deck_dyson.md` now carries the resolved slit-loss factor (slide
19), the fold prism as traced / prescription / scored / the cold-shield sweep
(slides 20-23, R5 of record = h 0 mm per the corrected 4d), the closure
envelope (slide 24 + the run-it-yourself sentence), the meniscus search
(slide 15), the native stage as "built and gated, result pending" (slide 25),
and the CALIB stride as the fifth engine finding (Backup).  Every refreshed
figure is the runner's own PNG from the re-emitted-with-apertures records.
Two producer items, yours, not blocking: (1) `dyson5_s5_layout_r5.png`: the
two panel titles overlap at the figure's width (the R4 title is shorter and
clears) -- the deck uses the engine's y-z render instead until it is fixed;
(2) the zoom inset of the fold at the base, as your 4d brief notes.
Engine note from 4b (a FarField kernel option to zero |f| > 1/lambda so a wide
window cannot carry evanescent energy into an energy-fraction metric):
acknowledged, queued behind item 4 of addendum 14; the record's propagating
normalisation stays the convention until then.

## Addendum 16 (2026-10-01): beat 4c sections 3.5 and 3.6 -- both FIXED (macos, local)
3.5 The asphere differential step is `das_rel = 1e-3` of the coefficient
(design_optim.F); a ZERO coefficient takes the step that moves the sag at the
element's circular aperture radius by 1e-7 |Kr|; with no circular aperture the
legacy 1e-10 step remains and one line says so (so declare the aperture --
every dyson5 deck does).  Your reproducer (`dyson5_s4_r4n_seed.in` +
`OptAsph= 2 1 2` on the block face) runs 4 LM iterations in 9 s on the CLI
where it was singular at once; the zero-term variant (`OptAsph= 1 3`) runs too.
3.6 The LM failure branch no longer `stop`s: the optics go back to the last
accepted parameter vector, `rtn_flg=1`, normal cleanup -- `macos.calib()` raises an ordinary MATLAB error
(`mmacos: calib_run failed`, the engine's reason printed just before it;
`dyson_native` should try/catch it) and MATLAB lives (the no-aperture variant of your
deck exercises it: "Optimization aborted; optics restored").  Gate
`tAsphCalib` (SUITE_FAST), fixture `Rx_AsphCalib.in` (a metre paraboloid with
a spoiled h^4 term CALIB drives back to zero; the failure path with the host
alive).  `P.native_asph` can come on when the fix is on the engine of record.
Deck: slide 25 carries the native result as a confirmation (R4n = R4; the
15-px keystone row as the reason for operands) and the asphere as the next
freedom.

## Addendum 17 (2026-10-01, Dave): CLEAR before beat 5 -- the resume list
Beat 5 is a new design (the telescope) and the engine of record now carries
every fix the spectrometer beats found (macos f1d2617, pushed; resources
b3d2644).  Start it from a CLEARED context, not a compaction: the numbers
that matter are in the record files, not in the conversation, and a long
context is where withdrawn numbers (the first R5 sweep) come back as facts.
Re-read, in this order, before the first command:
1. this brief: addendum 7 (the telescope spec = EMIT: 420 km, 60 m ground
   sample, 0.143 mrad IFOV, focal length 126 mm, 70 mm aperture at F/1.8,
   24.6 deg cross-track field onto the 54 mm slit, telecentric and flat,
   exit pupil on the grating), addenda 10, 13 (full slit diffraction waits
   for the telescope), 15-16 (deck state; the asphere step and the
   calib_run error are fixed -- try/catch `macos.calib()`);
2. `challenges/dyson5/README.md` + `BRIEF_dyson5_beat4d.md` section 2 (the
   fold's detector frame `G.fpa.xhat/yhat/normal` -- the telescope's image
   must land on the SLIT frame the same way) + `dyson5_s5.txt` (R5 of record
   = NO shield; the 2 mm/1.295/0.794 numbers are withdrawn);
3. `dyson5_params.m` / `dyson5_run.m` stage list (s0-s5, s4env): the
   telescope is a new stage through the SAME runner (Dave's rule: every
   stage through the parameterized runner), scored with the same gates
   (apertures declared, clearance per leg/body, engine renders per deck,
   chain-vs-engine identity to 1e-9 m);
4. `macos_f90/CLAUDE.md` sections dated 2026-09-30 / 10-01 (glass catalog,
   medium-aware kernels, chord-ruled gratings + modulo-lambda rule,
   short-line pad, CALIB stride / asphere step / failure path) and the
   design layer's telescope tools (`macos.design.Telescope`, `add_fold`,
   `realize_apertures`, `view_layout`) before writing a new emitter.
Deliverable of beat 5: the telescope deck of record, its score at the slit
(spot inside a pixel at every field, telecentricity, field flatness), the
end-to-end deck (telescope + R4 and + R5) traced as ONE prescription and
scored by the spectrometer's scorer, the clearance gate across both, and
two deck slides (the telescope as traced + the end-to-end layout and score).

## Addendum 18 (2026-10-01, Dave): the telescope goes THROUGH the design layer and the runner
Build it with `macos.design.Telescope` (TMA layout from the f-numbers,
Seidel-seeded conics, multi-field conic optimize, freeform / asphere
refinement, `add_fold`, `realize_apertures`, `view_layout`) -- read
`templates/10_telescopes/tma_widefield/example_tma_widefield.m` first (a
two-mirror hits a wall as the field opens; the third mirror buys astigmatism:
24.6 deg at F/1.8 IS that case) and `templates/80_end_to_end/e2e/README.md`
for the stage pattern (telescope stage -> instrument stage, one parameter
file, each stage consuming the previous stage's prescription).
It runs as a NEW STAGE of `dyson5_run` with its knobs in `dyson5_params`
(Dave's rule: one parameterized runner per modeled system) -- not a
standalone template first.  When it closes, extract the telecentric
wide-field TMA into `templates/10_telescopes/` (nothing there is telecentric
today).
Operands the Telescope optimizer does NOT have (it scores WFE or spot only):
telecentricity (chief-ray angle at the slit, mrad), field flatness, and the
PUPIL MATCH (exit pupil on the grating).  Hold them on the CHAIN side with
the lsqnonlin machinery you already run for smile and keystone; the
engine-side operand is queued with the centroid operand (addendum 14 item 4).

## Addendum 18 (2026-10-01, evening): beat 5 run from the cleared context -- the telescope is laid out and traced end to end; its image is the open item
Beat 5 was run in the cleared session per addendum 17 (resources
dev-candidate, local).  On record (`challenges/dyson5/BRIEF_dyson5_beat5.md`,
`dyson5_t1.txt`, `dyson5_t2.txt`, gate `tTelescopeRx`): the shared exact
tracer (`chain_trace`, the Dyson's lifted verbatim; R4/R5 re-emit
byte-identical), the telescope chain (`telescope_geom`: a folded
three-mirror, the field as the sky line that images onto the straight slit),
the end-to-end chain (`e2e_geom`, the grating the stop), the ladder, the
engine scorers, the runner stages t1/t2, two deck slides.  Facts: R4 is
telecentric at the slit to 0.09 deg (apparent pupil 16.84 m); the PNP first
order is a one-parameter family in that limit; the coaxial section cannot be
unobscured at F/1.8 with f-scale spacings -- the chief is folded at each
mirror (32/-32/21 deg).  Result: layout clears (+0.75 mm), pupil match 9.4
mm on the grating, 100 % admitted, EFL 126.4 by the map -- and 67 px rms at
the slit (a 1.6 mm field swing): conics + even aspheres on folded mirrors do
not image at the pixel over +-12.3 deg at F/1.8.  NEXT (TO or CC): the
two-mirror modified Schwarzschild (M&G: the widest field and lowest F-number)
and freeform mirrors on the folded three-mirror (the engine's Zernike
surfaces + OptZern; the chain lacks Zernike surfaces -- add them to
chain_trace); a cheaper first check: seed the conics from the coaxial
parent's anastigmat and fold afterwards.

## Addendum 19 (2026-10-02, Dave): beat 5b -- the telescope through the offset_imager ladder (rodgers3's product)
Written for a CLEARED TO (Opus 5.5 this round): re-read addendum 17's list,
then 18, then this, then `templates/10_telescopes/offset_imager/README.md`
and `BRIEF_dyson5_beat5.md` sections 4-5 (your own t1/t2 record: the layout
closes, the image does not -- 67 px rms, a 1.6 mm field swing).

**Why.**  Beat 5 started from spheres folded at the chief by 32 deg at each
mirror and let conics + even aspheres chase the field swing: that is the
offset_imager ladder's rung 4 without rungs 1 and 3.  The ladder that closed
Mike Rodgers' 20x20 deg box offset 22 deg at F/4 (`challenges/rodgers3`) is:
r1 coaxial anastigmat parent on axis, symmetric aspheres -> r3 the field box
pushed OFF AXIS, aspheres re-solved -> r4 + tilts, decenters, radii -> r5
8th-order Zernike freeform surfaces THROUGH THE ENGINE (Surface= Zernike;
`oi_zern_seed`, the asphere->Zernike seed is NEGATED).  A pushbroom
telescope IS an offset-field imager with a thin box: 24.6 deg cross-track x
~0.3 deg along-track, pushed off axis ALONG-TRACK to clear the beams.

**Run, in this order (every step through `dyson5_run`, a stage `t3`, its
knobs in `dyson5_params`; the template is called, not copied):**
1. `oi_story(struct('EPD_m',0.070,'Fno',1.8,'box_deg',[24.6 0.3],
   'offset_deg',OFF,'clear_m',[0.005 0.005],'exit_dir',<toward the slit>, ...))`
   at OFF = 4, 6, 8, 10 deg -- the along-track offset is THE unknown at
   F/1.8 (the template's wide-offset trial found 12 deg unpackageable at
   F/4 with a 20 deg box; a strip field is easier, a fat F/1.8 beam is
   harder).  Seed `z_m1_m`, `spacings_m`, `seed_R1_m` from beat 5's T0 first
   order (R [700 125 152], t [140 76 141] mm, Petzval 0), scaled to the
   template's sign conventions (read `oi_paraxial.m`).  If the cold solve
   stalls (the README's F8 rule), `oi_walk` the box from [6 0.3] outward.
2. Add the three residual ROWS the template lacks, in `oi_solve`'s exit-wall
   pattern (weighted rows on ITERATES, never a boolean wall -- lesson 3):
   telecentricity (chief-ray angle at the slit, mrad, per field),
   field flatness (best-focus z per field against the slit plane, um), and
   the pupil match (chief miss of the grating vertex, mm, with the R4
   spectrometer appended -- beat 5's t2 machinery).  Keep rodgers3's solver
   lessons: TRUE GN on stacked PER-RAY residuals (per-field rms plateaus);
   VERTEX-RADIAL natural scales for aspheres/Zernikes; spend on FIELDS not
   iterations, odd field counts across the slit.
3. Score every rung with beat 5's gates (engine render, chain-vs-engine
   1e-9 m, the clearance gate with the spectrometer's bodies, `tTelescopeRx`)
   and the END-TO-END row (telescope + R4 through the spectrometer's scorer:
   smile / keystone / SRF / CRF / admitted fraction).  Beat 5's T3 row (67 px,
   +0.75 mm) is the baseline every rung must beat; the target is the spot
   inside the 18 um pixel at every field with telecentricity < 0.5 deg.
4. Decision rule: if no OFF in step 1 gives a packageable r3 (clearance
   gate PASS) at F/1.8, STOP and report the number -- that sends the design
   to the review's two-mirror modified Schwarzschild (beat 5c), not to more
   iterations.  If r4 packages but r5 is what buys the pixel, the Zernike
   rung runs through the engine (the chain tracer does not carry Zernike
   surfaces; `tSpectrometerRx`'s identity gate is then engine-vs-engine on
   the re-emitted deck, as rodgers3 did).
Deliverable: the ladder table (r1-r5 x spot / telecentricity / flatness /
pupil miss / clearance / end-to-end CRF), the telescope of record if one
closes, deck slides 25-26 replaced (the design-of-record rule: show the
best, not the history), `BRIEF_dyson5_beat5b.md`, the t3 stage committed.
The edit-deck sync and every push wait on Dave's word.

## Addendum 20 (2026-10-02): beat 5b step 1 read -- the IMAGE closes; the clearance test could not pass as I specified it; scan the ENVELOPE before beat 5c
Read from your `t3/dyson5_t3_off{04,06,08,10}_REPORT.md` (CC, read-only).
**1. The ladder images.**  Dense-map max strict RMS WFE over the 24.6 deg
strip at 1 um: S1 coaxial 172 nm; S3 re-solved at the offset 162 / 173 / 196
/ 216 nm at 4 / 6 / 8 / 10 deg.  That is ~0.2 waves -- of order a tenth of an
18 um pixel of blur -- against beat 5's 67 px.  The form is alive on the
image side, and the offset costs little up to 10 deg.
**2. The clearance deficit is FLAT in the offset, so the offset is not the
variable:** -53.9 mm at the ON-AXIS S1, -52 to -59 mm at every offset and
rung.  Cause: the envelope.  Spacings are 140 and 76 mm for a 70 mm beam
(your S3 side view: the incoming beam passes through M2 and M3).  The lateral
walk an along-track offset buys is spacing x tan(offset) = 140 x tan(10 deg)
= 25 mm; clearing a 70 mm beam plus margins needs ~80 mm.  No offset in 4-10
deg could pass.  **That is my brief's error** (addendum 19 told you to seed
the spacings from beat 5's T0 and called the offset THE unknown); the
decision rule fired on a test that could not pass by construction, so beat
5c is NOT yet earned.
**3. Stop the four S5 runs.**  A freeform figure cannot move a beam 55 mm;
their only product is an image number on an envelope that cannot be built.
**4. The scan that decides it: envelope scale x offset, S1-S3 only** (your
timestamps: ~28 min per case).  Lengthen the spacings through your
first-order family (the t1 / y2 knobs of `telescope_seed`, EFL 126 mm and
Petzval 0 held) so that  t1 x tan(offset) >= 80 mm :
| t1 (mm) | offset (deg) | walk (mm) |
|---|---|---|
| 300 | 15 | 80 |
| 450 | 10 | 79 |
| 600 | 8 | 84 |
| 600 | 10 | 106 |
Report per case: the clearance floor WITH ITS WORST PAIR NAMED (the template
gate reports the pair; the flat -54 mm never said which), S3 map max, the
largest mirror diameter and the overall length.  A weak, long three-mirror
with a short focal length is an inverted telephoto, "large against f" like
the Schwarzschild; ~0.6 m beside a 0.7 m spectrometer is acceptable.
**5. Decision rule, restated.**  If a case packages (floor >= +5 mm) with S3
under ~250 nm: the TMA lives -- add the three residual rows (telecentricity,
flatness, pupil match), then S4 / S5.  If the image collapses at every scale
that packages, or packaging needs more than ~0.8 m: THAT is the number that
sends the design to beat 5c.

## Addendum 21 (2026-10-02): the envelope scan will likely come back flat too -- prediction on the record, and a first-order SCREEN before any more solves
Your note that the M2-to-M3 spacing stays at 76 mm in every case is the
tell.  Addendum 20's walk (t1 x tan(offset)) clears the FRONT-end pairs
(`in->M1 x M2`, `in->M1 x M3`); two BACK-end pairs do not depend on t1 at all.
From your `dyson5_t3_off08_s3.in` (M1 z 200, M2 = stop z 60, M3 z 136.2, FP
z -54 at y +17.46 mm; stop semi-diameter 29 mm):
- **`M3->FP x M2`.**  The design is telecentric (exit chief along the axis),
  so the image cone runs from M3 to the FP centred at the IMAGE HEIGHT,
  y' = f tan(offset) = 126 mm x tan(8 deg) = 17.7 mm -- the FP's own y.  It
  passes M2's plane 114 mm from focus with radius 114 / (2 x 1.8) = 32 mm.
  Needed separation from M2: 32 + 29 = 61 mm (more with the x1.15 disks);
  available: 17.7 mm.  Deficit ~ -45 mm, set by f, the offset and
  (BFD - t2) only: 22 mm at 10 deg, 34 mm at 15 deg.
- **`M1->M2 x M3`.**  Within 76 mm of the stop the M1->M2 beam is ~76 x
  tan(chief angle) off the axis whatever t1 is, and M3's patch is the same
  order on the other side: ~15-25 mm apart against ~60 needed.
**PREDICTION (pre-registered):** all four cases of the running scan report
floors of about -40 to -50 mm with the worst pair `M3->FP x M2` or
`M1->M2 x M3`; none packages.  Let them finish (S1-S3, ~28 min): the NAMED
pair is the measurement.  If a case packages, this addendum is wrong and
addendum 20's rule 5 applies.

**Next, before any further 28-minute solve: screen at FIRST ORDER.**  All
nine template pairs can be evaluated from `telescope_seed`'s first-order
layout in seconds, engine-free: per leg the chief height and the beam
half-width (marginal + field) at each obstacle's plane, per obstacle its
half-size; clearance = centre separation - (half-widths) x 1.15.  Tabulate
over the family's real knobs -- t1 (300-600 mm), y2 (0.3-0.9: it sets t2 and
the back focal distance) and the offset (8-30 deg) -- and print, per row,
the nine clearances, t2, BFD, the mirror diameters and the length.  Check
the screen against the template gate on the four finished cases first (it
should reproduce their floors and name the same pair to a few mm); then
solve ONLY rows the screen passes.
Two things the screen should look for: (a) a back end with the focal plane
near or inside M2's plane (BFD ~ t2: the image cone is small where it passes
M2, and the slit is thin along-track, so `M2->M3 x FP` asks only
f tan(offset) > M2's radius + a few mm, i.e. ~15 deg); (b) the offset at
which `M3->FP x M2` closes for the present back end, ~atan(61/126) = 26 deg
-- which is beat 5's 32 deg folds found from the other side.
**Decision rule (third statement, now with a cheap test):** if the screen
finds NO row with all nine clearances >= +5 mm at an offset <= ~15 deg, the
telecentric three-mirror of this family does not package at F/1.8 and beat
5c (the two-mirror modified Schwarzschild) is earned -- report the screen's
best row and its binding pair as the evidence.  If rows pass, solve the best
two through S3 and report image + gate.

## Addendum 22 (2026-10-02): the stalled 600 mm cases -- no verdict on EITHER count; the t1 continuation waits behind the screen
Your read is right: a cold S1 that stalls at 20 um is a solve outside its
basin, not evidence about the form.  One step further: **a stalled case's
CLEARANCE number is not a verdict either.**  The template's pairs are cut
from traced footprints, and the back-end pairs depend on the image cone
being the F/1.8 cone at the image height; a design that does not image puts
those legs somewhere a real design would not.  So both t1 = 600 mm cases
count for nothing, whatever S3 prints; the 300 and 450 mm cases count only
if their S1 converges.
**Order of work:** (1) let 300 / 450 finish; (2) addendum 21's first-order
nine-pair screen; (3) build the t1 continuation ONLY for rows the screen
passes -- if it passes none at <= 15 deg the continuation is never needed.
When it is needed, your design is the right one (warm-start each envelope's
S1 from the converged 140 mm solution, step t1 up, the box walk applied to
the spacing).  One thing to carry into it: the seed holds R1 = 5 x t1, so M1
weakens as the envelope grows, and its conic's leverage on the wavefront
falls as 1/R1^3 (16x weaker at 3.0 m than at 1.2 m) -- K1 ran to +17.6 while
the merit did not move.  Scale K1's step with R1^3, or hold K1 and let M2 /
M3 carry the correction on the long envelopes (rodgers3 lesson 2: natural
scales must match the variable's real leverage, or the LM damping runs away).
Do not chase the 8-vs-10 deg trace difference or the smoke-test mismatch
(`pose_stop_once_`) unless a screened row depends on it.

## Addendum 23 (2026-10-02): beat 5b read -- the screen did its job; one more round (the y2 continuation) with a hard stop
Read from `BRIEF_dyson5_beat5b.md` (CC).  Accepted: the first-order screen
(`tma_screen`, validated to 3-5 mm against the engine gate, same binding
pair), its finding that t1 is NOT the knob and the back end closes at y2 =
0.3-0.4 near 12-14 deg, the engine's confirmation that y2 = 0.4 packages
(+17.5 / +10.4 mm), and the stall rule as you applied it.  My ~26 deg figure
in addendum 21 was for the y2 = 0.6 back end; the short-BFD route (the other
option in that addendum) is the one that closed.  The t1 continuation of
addendum 22 is withdrawn: no row needs it.
**Next: your y2 continuation, with these rules.**
1. Walk y2 0.6 -> 0.55 -> 0.5 -> 0.45 -> 0.4 at t1 140, S1 only, each step
   warm-started from the previous step's solved S1, the R2/R3 branch held
   (`seed_R_m`), `nsolve` 5 across the slit.
2. **Raise the iteration cap.**  "Still descending at the 12-iteration cap"
   has now appeared three times; a capped solve is a stall by another name.
   Run S1 and S3 to convergence (cap 40; stop on the LM's own criterion) and
   report the iteration count.
3. A step COUNTS when its S1 dense-map max is <= 1000 nm (your
   `tel3_s1_conv_nm`).  If a step fails to converge, halve the y2 step once;
   if it still fails, the walk ends there.
4. At y2 = 0.4 (or the last converged step that the screen passes): S3 at
   14 deg seeded FROM THAT S1, not from a sphere.  The template's S3 has no
   clearance row and has already moved a layout from -7 to -38 mm, so a
   solved S3 is accepted only if the gate still reads >= +5 mm afterwards;
   if it does not, run S4 (which carries the clearance hinge) from it and
   report that instead.
5. Report per step: S1 max / avg, iterations, the radii and conics, M3's
   rho/|R| at the box corners, the fraction of rays surviving at the
   +-12.3 deg cross-track fields (the y2 = 0.4 seed already loses 17 of 277
   there), and the exit chief error.
**HARD STOP for this form.**  Beat 5c (the two-mirror modified Schwarzschild)
is earned if ANY of these holds: the walk ends above y2 = 0.4 with no
screen-passing step; S1 at the packaging corner cannot reach 1000 nm; or S3
/ S4 at 14-15 deg ends above 250 nm with the gate satisfied.  It is also
earned by vignetting: more than 5 % of rays lost at the cross-track edge on
the solved design.  Report whichever fires with its number.  If none fires,
the three-mirror lives: add the three residual rows (telecentricity,
flatness, pupil match) and go to S4 / S5.
**The template fix (section 4 of your brief) before it is committed:**
`tOffsetImager` and `tRodgers3` (SUITE_FREEFORM) green on the DEFAULT path --
"rodgers3 untouched" has to be a gate count, not a statement -- plus one new
test: with `seed_R_m` at y2 = 0.3 the re-solved M2 stays CONVEX at every
iterate, and without it the fixed Newton start lands on the concave root
(the must-fail leg that documents the defect).
**Commit your beat-5b work locally now** (`git add` BY PATH: the dyson5
files, `tma_screen.m`, the four template files; the bench and sensitivity
files showing as modified in that tree are not yours).  Push waits on Dave.

## Addendum 24 (2026-10-02): beat 5b accepted; HOLD beat 5c -- the stop fired on a solve whose own gate is red, and the telescope's target is about to change
Read from your committed `BRIEF_dyson5_beat5b.md` sections 6-9 (0a70058).
**Accepted, and good work:** the y2 walk (every step converged on its own
test), the `oi_score` telecentric fix (the same degenerate crossing the
engine's FEX guards; your note that sections 1-5's image numbers are
superseded is right -- including the 162-216 nm I leaned on in addendum 20,
which was bought at R1 1.2 m, off the family), `hold_R1`, `seed_R_m`,
`tOiSeedBranch`.
**What the walk establishes, and it is the real finding:** in this family at
24.6 deg, packaging and imaging pull opposite ways.  The on-axis parent goes
349 -> 787 nm as y2 goes 0.60 -> 0.40, M3's reach jumps to 0.89 of its radius
at the packaging corner, and the edge field starts losing rays there.
**Why 5c is NOT yet earned.**
1. The stop fired on S3 at 14 deg: start 552 081 nm, "4 iterations (own
   stop)", map max 1 206 588 nm.  A solve that ends where it started is a
   stall (addendum 22), and it was a cold jump from on-axis to 14 deg.
2. The template's own gate for exactly that solve is RED on the committed
   tree: `tOffsetImager/test_s3_resolve_recovers`, S3 29 274 vs S2 29 024 nm
   where 2026-08-20 recorded s3/s2 = 0.71.  A design cannot be ruled out by
   a solve whose smoke test fails.  **CC is chasing it now** (engine-vs-
   MATLAB bisect: the August test commit against today's engine).  Do not
   spend time on it.
3. **The target is moving.**  CCMac's block-size trade (resources 9d3b6e2,
   round 2 in progress): with the slit shared by two 27 mm modules the Dyson
   needs no meniscus and a 100-130 mm block.  If Dave takes two modules, the
   telescope question becomes one 24.6 deg field split ahead of the slits,
   OR two telescopes of 12.3 deg -- half the cross-track field, where M1's
   width, M3's reach and the edge losses all ease.
**Until Dave rules on the modules and the S3 gate is understood:** no beat
5c work.  One cheap run IS worth doing now, because it is S1 only and S1 is
sound: the same y2 walk at a 12.3 deg cross-track box (t1 140, R1 held,
branch held, 5 x 3 solve set).  Report the same table.  If the on-axis
parent at the packaging corner is several times better at half the field,
that is the number Dave needs for the module decision.

## Addendum 25 (2026-10-02): the hold on S3 is released for the 12.3 deg telescope; the red gate is understood
**Your 12.3 deg walk is the number Dave needed:** at the packaging corner
(y2 0.40) the on-axis parent is 93.8 nm against 786.7 nm at 24.6 deg, M3's
reach falls from 0.89 to 0.57 of its radius, no edge rays are lost.  Commit
it locally (the opt-in knobs `tel3w_xtrack_deg`, `tel3w_s1_only`,
`tel3w_suffix`, the record and decks), by path.
**The `tOffsetImager` failure, bisected (CC):** engine commit 81d3308
(2026-09-08, re-traces made idempotent) is where `test_s3_resolve_recovers`
turns red; the engines of 08-22, 09-05 and the two commits just before it
pass with identical numbers.  Mechanism: before that fix every
finite-difference Jacobian column carried a few-ulp 2-cycle of the re-trace;
after it they are clean.  The test's 5-iteration S3 is ill-conditioned (LM
damping ~1e-6), so the cleaner Jacobian moves step 1 by 0.006 nm, step 2 by
7 %, and step 3 is rejected at 23.7 um where the old path reached 13.6 um.
A path, not a property: with a 15-iteration budget today's engine takes that
instance from 906 to 167 nm (s3/s2 = 0.18).  **The S3 solver is healthy at
your cap of 40.**  CC re-pins the test; nothing for you to do.
**So the full-field S3 result stands as what it is** -- a cold jump from an
on-axis parent of 787 nm to 14 deg that started at 552 um and did not move.
It is a stall, not a verdict, but the 24.6 deg three-mirror is not worth
more time: its parent is three times over the bar before any offset.
**Next: S3 at 14 deg on the 12.3 deg box, from the y2 0.40 parent (93.8 nm).**
1. Direct first: stop re-posed at the offset, R1 and branch held, 5 x 3
   solve set, cap 40, own stop.  Report the start value, the per-iteration
   trace, the dense-map max and average, clearance after S3 with its worst
   pair, exit error, M3 reach, edge rays kept.
2. **Stall test, stated in advance:** a solve that stops within 5 iterations
   having gained less than 20 % from its start, with the damping run up, is a
   STALL.  If that happens, walk the OFFSET instead: 0 -> 5 -> 10 -> 14 deg,
   S3 at each step seeded from the previous step's solved design with the
   stop re-posed (the box walk applied to the offset; your step-1 runs solved
   S3 at 4-10 deg from a carried design, so the first steps are known to
   work).  Halve a step once if it stalls.
3. **Rule.**  Dense-map max <= 250 nm with the gate >= +5 mm after the solve
   (S4 with its clearance hinge if S3 un-packages it): the three-mirror
   LIVES for a two-telescope instrument -- then the three residual rows
   (telecentricity, flatness, pupil match) and the end-to-end deck against
   the two-module, no-meniscus Dyson (CCMac's `dyson5_size_D_r130.in`, 27 mm
   slit, 130 mm silica block).  Converged ABOVE 250 nm: run S5 (freeform)
   once before any other form is considered.  The walk cannot reach 14 deg:
   report where it ends and why.
Geometry note for the record: a 12.3 deg box centred on the axis is ONE of
two telescopes, each 70 mm at F/1.8, their axes +-6.15 deg apart cross-track,
each feeding one 27 mm slit.  One 24.6 deg telescope with a field splitter
is the other reading of "two modules" and is NOT what this run tests.
