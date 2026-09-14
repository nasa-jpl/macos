# PLAN: one deck for the DM surface gauges -- configurations, layouts, parts, performance

Dave 2026-09-13: "bring all the DM surface gauge concepts together into
one new deck, comparing performance, showing layout and parts, for the
lab bench configurations."  CCMac's reflective work is not complete; TO's
PDI is begun, not complete.  This is the plan; nothing is built yet.
Written by CCL for Dave's ruling on the decisions at the end.

## 1. What the deck is for

One question: which DM surface gauge to build on the bench, and what it
takes.  Audience: the people who will build and run it (Fang's group and
JPL), so every configuration gets its layout drawn to scale from the
engine, its parts named, and its performance in the one currency every
lane already scores in (`dm_gauge_lib`: actuator-space rows on the same
96x96 DM and the same 30 nm working surface, the same seeds; photons per
measurement; the closed-loop hold; the capture range to 10%).  The three
existing decks (deck_zwfs 28 slides, deck_tg_fang, deck_pdi) stay as the
lane records; the new deck takes their figures and numbers, not their
narrative.

## 2. The configurations (six, on one front end)

All six share the TG96 front end: source, collimator L1, the 7-deg
beamsplitter, the 700 mm leg to the 96 mm DM, and the tuned pupil-imaging
tail.  They differ in what sits at the internal focus and what follows.

| # | configuration | what is added to the front end | lane | state |
|---|---|---|---|---|
| A | Twyman-Green IFO, lens rig, four-step PSI | reference flat on a PZT (or the polarization snapshot: polarizers, arm QWPs, MacNeille cube, analyzer -- the v1/v2 rigs) | CCMac | rows, loop, coatings done; layout figure needs work; parts list not written; capture range partial |
| B | Twyman-Green IFO, reflective (OAP) rig, bare Al | two OAPs (bare or protected Al) replacing the lenses in the test arm | CCMac | open-loop rows done; loop shows it does NOT hold a 2 pm walk (fold cross-talk 0.18); design number not closed |
| C | Zernike sensor, scalar dimple (readings L / I+ / S) | the etched-dimple mask substrate at the focus (VSG2 part); nothing else | CCL | complete: rows, capture range, photons, loop, systematics; layout done |
| D | vector Zernike (polarized dimple, reading V) | geometric-phase metasurface in the mask seat; quarter-wave plate + 12.7 mm MacNeille cube behind the field lens; second camera | CCL | complete but for the cube's channel leakage (unpriced); layout + decks done |
| E | point-diffraction, stepped pinhole (reading P) | pinhole mask substrate with an attenuating surround, stepped like the dimple | TO | rows, loop, step-scheme and camera-drift trades done; capture range print + parts list to do |
| F | P/SRI, waveguide reference (reading PF) | pickoff BS, reference arm (Lr1, pinhole, Lr2, folds, compensator), single-mode waveguide + photonic phase shifter, recombination BS | TO | buildable two-arm bench drawn (psri_bench); PF still synthesized, not yet run through the two decks; reference-arm drift unpriced |

## 3. The common currency and where each configuration stands

Every cell below is either a run tag (exists) or a name (who runs it).

| metric (96x96, 30 nm working surface, matrix measured on it) | A lens IFO | B OAP IFO | C ZWFS S | D vZWFS V | E PDI P | F P/SRI PF |
|---|---|---|---|---|---|---|
| 10 nm single: gain / floor | 0.99 / 4.8 pm (lens_base ladder at 30 nm; its rows are on a 16 nm base -- CCMac re-runs the rows on the 30 nm surface) | 0.996 / 2.8 (same caveat) | matbase385 0.989 / 5 | v193base 0.9935 / 4 | pdi193fbase | pdi193fbase |
| 1 nm on 47 sites | CCMac (47-site row, 30 nm base) | CCMac | 0.9993 / 4 | 0.9992 / 3 | pdi193fbase | pdi193fbase |
| dense random 10 nm | 0.988 / 686 pm (16 nm base) | 0.95 (bare Al, oap_bareAl) | 0.984 / 0.68 nm | 0.9999 / 0.33 nm | pdi193fbase | pdi193fbase |
| capture range to 10%, aging from 30 nm | single-site 120-158 nm (wrap); CCMac: 47-site + the runner's print | 60-120; CCMac | 42 nm (cap385) | 70 nm | TO (ladder run; the print exists in zwfs_run) | TO |
| capture range, matrix re-measured on the surface + photons | CCMac | CCMac | to 160 nm; photons x40 | to 160 nm; x40 | TO | TO |
| photons per measurement for 1 pm (S5 form) | CCMac (D7 gives sig_n 12 pm at 1e12 -> ~1.4e14; state it in the noise-stage form) | CCMac | 5.6e13 | 4.7e13 | 3.3e13 | TO |
| closed loop: 3 pm noise-only / 2 pm walk | 5.5e12 / 2.0e13 | 3.4e13 / never | 2.6e12 / 7.5e12 | 1.5e12 / 5.3e12 | 2.3e12 / 7.0e12 | 7.0e12 / 2.5e13 |
| thermal floor / noiseless step | 13.1 pm / 87 pm rising (8.6%) | 39 / 276 (27.6%) | 10.0 / 0.000 | 9.9 / 0.000 | ploop193 | ploop193 |
| camera drift (1/f, electron and relative forms) | CCMac (dmg_loop 'cam' exists) | CCMac | pcam193* | pcam193* | pcam193* | pcam193* |
| systematics priced | BS diattenuation (v1/v2), coatings (item B), stencil bias; PZT step error: CCMac | coating retardance band | model defocus, sampling, color | metasurface retardance (V2), arm polarization (V3); cube leakage: CCL | step error, camera drift | reference-arm drift: TO |

## 4. Layouts and parts -- what exists, what each lane owes

Rule for every layout figure (Dave 2026-09-12): drawn by the engine from
the emitted deck (`macos.view_rx`), the fold plane seen from above,
passive bookkeeping planes hidden, elements named (not E-numbers) with
leader labels, the crowded node (beamsplitter, mask seat) as a cropped
panel at full slide width, type readable at slide size.  The recipe is
`zwfs_dm96/zwfs_vlayout.m` (deck_zwfs slide 17).

| configuration | layout figure | parts list |
|---|---|---|
| shared front end | zwfs_render_rig.png (deck_zwfs slide 3) -- redo in the recipe | CCL: source, L1, BS plate + compensator, DM, L2, the tail (FL, camera); apertures and focal lengths from `zwfs_params` / `twyman_green` |
| A lens IFO | CCMac's lens_render (b6ed9fe) -- E-numbers, no reference arm; redo in the recipe with the reference arm, the PZT flat, the BS node inset (CCMac's own offer: yes to all three) | CCMac: reference flat + PZT, the polarization snapshot parts (v1 plate / v2 cube rigs), compensator |
| B OAP IFO | oap_render: the table-plane panel is edge-on (fold plane not in view); redo | CCMac: OAP1/OAP2 (off-axis distance, f, AOI, coating), folds |
| C ZWFS | done (slide 3 + slide 2 mask figure) | CCL: mask substrate (VSG2, 9 spots, etch), camera (385 px per pupil) |
| D vZWFS | done (zwfs_vlayout.png, slide 17) | CCL: metasurface (geometric-phase HWP, dimple pattern), QWP, 12.7 mm MacNeille cube, camera B |
| E PDI | pdi_layout.png / pdi_layout_tail.png exist (deck_pdi) -- check against the recipe | TO: pinhole substrate (diameter, surround transmission t, the stepping) |
| F P/SRI | psri_layout.png / psri_render.png exist (deck_pdi slides 5-6) -- check against the recipe | TO: pickoff BS, Lr1/Lr2, pinhole, folds M1/M3, compensator, waveguide chip + phase shifter, BS3 |

## 5. Dave's rulings (2026-09-13)

1. Audience: the JPL HWO WFS&C discussion group.  Title: **DM Surface
   Gauge Comparison**.
2. OAPs are the traditional (budget) implementation.  The deck must show
   why the lens route is worth its cost, if it is, and what the OAP
   front end means for the OTHER approaches -- so the ZWFS, vZWFS and
   PDI readings get run on the OAP rig too (section 8, CCMac).
3. All three IFO phase-shift forms are compared: PZT four-step,
   polarization snapshot (v1 plate, v2 cube), and the hybrid (snapshot
   for the change measurements, PZT as the calibrator).
4. TO's points, ruled here: the P/SRI models the paper's pickoff form
   (the pinhole form P stays as the stepped-pinhole reading); pinhole
   diameter of record = whichever of 2.0 lam/D (1024/193) and 1.0 lam/D
   (2048/385) performs better on the surface rows and the loop, with
   both recorded; the PDI gets its own directory `40_benches/pdi_dm96/`
   (pdi_params, pdi_run, the psri bench and decks, README, runs), the
   code staying shared through zwfs_run / dm_gauge_lib.
5. The draft stance is approved for performance.  Added: **capturing the
   DM's initial figure, ~100-200 nm WFE (50-100 nm of surface)** -- every
   configuration must show how it gets from there to the hold regime
   (section 7).
6. CCMac's renders and layouts are not of deck quality: redo in the
   recipe; CCL QAs every figure before it enters the deck.
7. CCMac and TO are tasked at the Opus level (briefs in section 8, written
   to be executed without CCL); CCL does QA, the vZWFS cube leakage, and
   the assembly.
8. Style: succinct; add slides rather than crowd one; jargon-free,
   fact-based, no fluff.  Slide rule: one figure or one table, at most
   three bullets, at most ~60 words of bullets, a footnote that fits.
9. Time is limited: the main body presents each APPROACH once, in its
   best-performing configuration; the alternatives and the stories of how
   we got there (the model correction, the readings that lost, the rigs
   that do not hold) are backup slides.  Each lane's report names the
   best configuration of its approach with the numbers that make it so.

## 6. Deck outline (deck_gauges.md; main body ~15 slides, backup the rest; DRAFT until Dave signs)

Main body -- each approach once, in its best configuration:
1. Title: DM Surface Gauge Comparison.
2. The requirement: picometers on a 96x96 DM holding a 30-60 nm surface
   in a servo; first capture of a 100-200 nm WFE figure; the four scores.
3. The shared front end: layout + parts.
4. Interferometer, best configuration (lens rig; the phase-shift form the
   numbers pick -- expected the hybrid): layout + parts.
5. Zernike sensor, best configuration (the stepped reading, matrix on
   the surface): layout + parts.
6. Vector Zernike sensor (the polarized pair, cube split): layout + parts.
7. Point-diffraction, best configuration (stepped pinhole or P/SRI, as
   the numbers pick): layout + parts.
8. Performance side by side (rows on the 30 nm surface).
9. Capture range (aging and re-measured, with the photon cost).
10. Photons for 1 pm and the closed-loop hold (3 pm vs light; fixed
    errors) -- two slides if one overflows.
11. Capturing the initial figure: each approach's route from 100-200 nm
    WFE to hold (the descent run).
12. Lenses vs OAPs: the budget question answered with the IFO on each and
    the other gauges on the OAP front end.
13. Systematics: priced and open, one line per approach.
14. Recommendation.
15. The modes, and the flow from launch to hold (section 11.1, a flow
    diagram).
16. The complex amplitude for two-DM control: which readings give it
    (section 11.2).
17. One bench, all modes: what is inserted or removed to switch
    (section 11.3, the universal-bench layout).
18-20. Future work (section 10).
21. Run it yourself.

Backup: the three IFO phase-shift forms in detail; the OAP rig's rows and
loop; the scalar Zernike readings that lost (linear, one-frame exact) and
the fold; the pinhole-diameter and P/SRI-vs-pinhole trades; drift
(camera, DM, within-measurement); the metasurface and arm terms (V2,
V3); the model correction and sampling story; the coating story (item
B); provenance.

## 7. Capturing the initial figure (new measurement)

The DM arrives with 100-200 nm WFE (50-100 nm of surface).  The capture
range slides say: with a calibration aging from 30 nm every ZWFS reading
is dead by 100 nm of surface; re-measured on the surface it holds to 160
nm at 5-40x the light; the IFO reads to its wrap (158 nm of surface,
316 nm WFE) and never folds; the P/SRI's reference does not depend on the
surface (TO's ladder 1.02 / 1.06 / 1.13 at 120 / 240 / 480 nm).  The
measurement that settles it is a DESCENT run: the loop started at a 100
nm rms surface (200 nm WFE) with the matrix measured there, gain 0.5,
recalibrated every K cycles (or once), photons per cycle at the hold
level -- does each reading converge to the hold regime, in how many
cycles, at what light?  Knobs `loop.start_rms`, `loop.recal_every` in
dmg_loop / zwfs_run (TO builds, section 8), mirrored in tg96_run
(CCMac).  Two more routes to state on the slide: capture with the IFO,
then hold with the sensor (the hybrid bench); and the longer wavelength
(the color stage: 780 nm widens every fold by 1.23x).

### 7.1 Capture results so far (2026-09-14, TO descent ladder at 193 rays; CCMac descent_lens at 1024)

The premise "capture is a wrap problem" held for one reading only.
Largest starting surface reaching 3 pm in 40 cycles, matrix at the
start, 1e13 and 1e15 alike: L 30 nm, S 30, V 60, P 60, PF 60 -- unwrap
on or off, bit-identical for S / V / P: past ~60 nm their maps carry no
2 pi residues and no gradient near pi; they return 24-28 nm whatever the
truth (70 to 270 nm).  Their focal-plane REFERENCE has collapsed; they
are blind, not folded.  Only PF's map genuinely folds (4928 residues at
70 nm), and PF from 100 nm captures only with BOTH the unwrapper and
re-calibration (neither: 64 nm; unwrap: 5.4 nm contracting; recal: 64
nm; both: 0.245 pm, 3 pm at cycle 26), at 1e13 as at 1e15 --
reference-limited, not light-limited.  The interferometer (CCMac,
unwrapper, no re-calibration) captures from 150 nm rms to 2.1 pm at 1e13
(0.21 at 1e15); 200 and 300 nm running.  Pending: P with a shutter frame
(tracks PF to four digits so far), the OAP descent.

Update 2026-09-14 (TO): the recommended configuration captures -- the
stepped pinhole P WITH A SHUTTER FRAME (the reference re-measured per
state), unwrapped, re-calibrated every 10 cycles, from a 100 nm surface
(200 nm WFE): 0.207 pm at 1e15, 2.07 pm at 1e13, 3 pm at cycle 26,
contraction 0.716 (measured, cap_state_uw_recal); at 150 nm it reaches
10 nm, not 3 pm -- the ceiling is between 100 and 150 nm.  Recalibration
cadence is not the constraint (every 2 cycles 6.5 pm, every 5 11.2, every
10 enough); the cycle count is.  So the per-state reference is what a
self-made reference needs to capture: the plain dimple and pinhole do
not, the shuttered pinhole and the P/SRI do, the interferometer does.

Within-scan DM drift (V4, TO, `loop.intra`): it HELPS the stepped
readings -- 2 pm walk at 1e15: S 2.32 -> 1.57 pm, P 2.32 -> 1.58, PF
2.33 -> 1.71; 5 pm ramp: S 10.05 -> 7.16, P 9.86 -> 6.70, PF 9.88 ->
7.30; L and V unchanged to the digit (the plumbing gate).  A scan
measures the surface at its midpoint, half a measurement closer to now:
a free half-step of prediction against a smooth drift.  Camera drift
within a scan costs 5-11 pm because it is additive bias, which a
zero-sum scheme cancels only while constant across the scan.  Same knob,
opposite signs -- the slide states it that way (the V4 framing
inverted).

Consequence for the recommendation: a self-referenced reading (the
Zernike dimple, scalar or vector; the stepped pinhole) captures to ~60
nm of surface (120 nm WFE); an externally referenced one (the
interferometer's flat, the P/SRI's waveguide) captures the whole 100-200
nm WFE range with unwrapping.  So the bench needs either the hybrid
(interferometer or P/SRI for capture, the sensor for hold), the P/SRI
alone, or the sensor with a second color for capture (10.1).  Slide 11
carries this table; slide 14's stance is updated to it.

### 7.1b The P/SRI's reference-arm walk (TO deliverable 5, 2026-09-14; slide 21 / 13)

A path-length walk between the P/SRI's arms is a piston on the retrieved
phase, the mode the actuator estimator nulls (the S10 rank-one term), so
it is nearly free in hold: 3 pm under the 2 pm DM walk at 4.8e13 /
4.9e13 / 5.0e13 photons per cycle for 1e-3 / 1e-2 / 1e-1 rad per cycle
(4% across a hundredfold range); the common-path pinhole P is
bit-identical (no such arm; below 1e13).  Scoped: for absolute field
reconstruction (the paper's use) the same walk is a direct error; a
reference arm that TILTS is not a piston and is not modeled -- the
as-built term to carry (10.3).  deck_pdi stays the lane's draft record;
the comparison deck takes the report's numbers.

### 7.2 The OAP front end under the mask gauges (CCMac fe36c66, 2026-09-14; slide 12)

Resolved: the OAP rig's seat needed a solved trim (MASK_TRIM 6.14 mm;
the root is OAP1's residual collimation defocus refocused by OAP2), not a
geometry fix -- the "8 um, 3 lambda F/D" of CCL's scan was defocus from
the scan's 0.25 mm rung spacing (0.11 mm off best focus at NA 0.069); at
best focus the blur is fold coma, 0.82 lambda F/D at 9 deg AOI, and G1 =
2e-15, G3 = 4e-16 pass.  The gauges then split by how focus-critical
their mask feature is: the stepped dimple ZWFS survives (single 10 nm
1.000 / 4 pm, grid SNR 221, capture 38 nm -- comparable to the lens rig),
the vector pair fails its fold gate (19.6 pm on the 100 nm pokes), the
pinhole fails its exactness gate (94 pm).  Deck statement: the OAP front
end carries the interferometer and the dimple ZWFS on a marginal,
~0.05 mm-alignment-critical focus, degrades the vector reading and
breaks the pinhole -- the more focus-critical the mask feature, the
worse the fold coma.  Open, not for this deck: trace-solving the OAP
seat trim inside twyman_green (the S1 recipe; a builder change that
risks the lens rig's byte-identical gate).

## 8. Tasking (Opus level; each brief stands alone)

- `BRIEF_ccmac_gauge_deck.md` -- CCMac: rows on the 30 nm surface (lens,
  OAP); capture range both ways + photons; the three phase-shift forms'
  systematics; the descent run; the ZWFS / vZWFS / PDI readings on the
  OAP rig; layouts in the recipe; parts lists; close the reflective
  design.
- `BRIEF_to_gauge_deck.md` -- TO: `pdi_dm96/`; PF through the two decks;
  the pinhole-diameter choice; capture range + photons; the descent
  knobs in dmg_loop + the within-measurement DM drift knob (V4, shared
  with the IFO's PZT form); the reference arm's drift; layouts checked
  against the recipe; parts lists; deck_pdi conclusions.
- CCL: the vZWFS cube leakage; figure QA; the combined figures; the deck.

Gate before assembly: every number on slides 14-22 has a run tag in a
committed report; every layout figure passed CCL's QA at slide size.

## 9. Decisions still open

None blocking.  The recommendation slide's wording is Dave's at sign-off.

## 10. Future work (three slides; Dave 2026-09-13)

### 10.1 Capture beyond one wave: a second color

Unwrapping resolves a wrapped differential only while neighboring pixels
differ by less than pi: at 4 px per actuator a 100-200 nm rms actuator
pattern is 0.7-1.4 rad per pixel and unwraps; a figure of more than a
wave with actuator-scale structure (quilting, print-through, a stuck
actuator) does not.  The second color removes the ambiguity without
unwrapping: two wavelengths give a synthetic wavelength lambda1 lambda2 /
|lambda1 - lambda2| -- 632.8 + 700 nm: 6.6 um, an unambiguous surface
range of 1.6 um (double pass); 632.8 + 780 nm: 3.35 um, 0.84 um -- at a
noise cost of lambda_syn / lambda (5-10x), so the coarse color pair
captures and the single color holds.  The machinery exists: the color
stage runs one physical mask at five wavelengths (its phase 1.57 rad at
632.8 nm, 1.42 at 700; the stepped depths scale as (n-1)/lambda), and a
two-color IFO is the classic form.  Measurement to add: the start-rms
ladder with a two-color coarse solve feeding the single-color loop, per
reading; the number is the largest capturable start and its light.

### 10.2 Other approaches to consider

One line each; none is modeled yet.
- Shack-Hartmann or modulated pyramid as the capture stage: range of
  many waves, no wrap, sensitivity far from picometers; hands off to the
  gauge once inside its range.
- Phase diversity (two defocused pupil images, no mask): wide range,
  iterative, common path; a candidate for capture with the existing
  camera and a translation stage.
- White-light (low-coherence) scanning in the interferometer: absolute
  surface with no wrap, slow; a one-time capture tool.
- Model-based large-figure solve: the DM's influence functions as the
  basis of a nonlinear (iterative) fit to the sensor's frames; the
  actuator-space estimator already carries the basis, the nonlinearity
  is the addition.
- Vector Zernike with a polarization camera (micro-polarizer array):
  one camera instead of the cube and two; costs a quarter of the pixels
  per image.
- Heterodyne or lock-in detection: immunity to slow drift and 1/f
  electronics, at the cost of a frequency-shifted reference (the
  interferometer and the P/SRI can carry it; the common-path sensors
  cannot).
- Direct actuator metrology (capacitive, optical) as the coarse reference
  the optical gauge is calibrated against.

### 10.3 The path to as-built performance

The model today: ideal optics, a perfect camera, photon noise, DM and
camera drift.  As-built adds, in the order they are likely to matter:
1. **The camera's throughput sets the measurement time, not the laser.**
   1e14 photons per measurement over ~1e5 pupil pixels is 1e9 electrons
   per pixel; a 1e5-electron well means 1e4 co-added frames -- at 100
   frames per second, 100 s per measurement, against the 0.13 s the
   laser needs.  Price it: well depth, frame rate, read noise per frame
   (which then enters as sqrt(frames) x read noise), bit depth
   (quantization at 1e-4 of the signal), gain nonlinearity and pixel
   response nonuniformity (a flat-field term the differential mostly
   cancels), persistence between frames (the stepped readings).
2. **Optical surface errors:** each lens, OAP, plate and the mask
   substrate with a typical figure (lambda/10 to lambda/20 PV, a 1/f^2
   spectrum) as GridData on the element; the sensor's reference core
   sees the low orders (the S7 lesson: 0.16 of defocus moved every
   number), the differential rows cancel what is fixed.  The mask
   substrate's flatness and the dimple's etch-depth error (the
   metasurface's retardance error is V2) are the ZWFS-specific ones.
3. **Alignment and stability:** the alignment sensitivities exist (D4:
   10 um, 10 urad); add the mask's centering drift on the focal spot
   (the S1 registration is the sensitivity), thermal expansion of the
   bench (a leg length per degree), source pointing and wavelength drift
   (the mask phase goes as 1/lambda), laser polarization angle (V3), and
   for the two-arm instruments the non-common-path air and mount drift
   (the P/SRI's reference walk is TO's item 5).
4. **More drift terms in the loop:** actuator hysteresis and creep,
   command quantization (14-16 bit: 0.1-0.5 nm steps as a floor), the
   influence-function error (absorbed by a matrix measured through the
   sensor, not by a modeled one), vibration within a stepped scan (the
   PZT form; `loop.intra`).
5. **The error budget in the JPL form:** per reading, fixed terms
   (calibratable, the residual after calibration), drift terms weighted
   by the servo bandwidth, and noise terms (photon, read, quantization),
   summed in quadrature to a held error at a stated light and time --
   the sensitivity-factor form the ZWFS deck's conclusions name.  The
   deck's comparison table becomes the first column of it.

## 11. Dave's notes of 2026-09-14: modes, complex amplitude, one bench

### 11.1 The modes, with a flow diagram (slide 15)

The deck names the modes it tests and draws the flow between them:

- **Mode 0, ground flat:** the ground-calibrated voltage map (~0 WFE on
  the ground) is applied; on orbit the residual is the launch, gravity-
  release and thermal change, 100-200 nm WFE.
- **Mode 1, image-based phase retrieval:** the WFS&C loop's own
  phase-retrieval (focal-plane images, no gauge) drives the WF toward
  the gauge's capture level.  It wraps too: at the half-wave level with
  high-spatial-frequency WFE -- the same failure class as the gauges'
  wrap, one wave earlier.  Its reach sets what mode 2 must capture.
- **Mode 2, capture:** an externally referenced reading (interferometer
  or P/SRI) with unwrapping, re-calibrated on the surface as it moves
  (7.1); or the sensor with a second color (10.1).  Ends at the hold
  regime (~30 nm of surface, then the matrix measured there).
- **Mode 3, closed-loop hold:** the sensor at picometers (the stepped or
  vector Zernike reading, or the pinhole), gain 0.5, the matrix
  re-measured on the held surface when the calibration ages (slide 9).
- **Recalibration events** between modes: the matrix on the current
  surface (photon cost, 7.1); the flat re-taken.

The figure is a boxes-and-arrows flow with the capture limit of each
reading written on its arrow (30 / 60 / 100 / 150 nm of surface, from
7.1), drawn by a script (graphviz `dot` -> PNG, committed with the deck
as `gauge_modes_flow.py`); CCL, at assembly.

### 11.2 The complex amplitude: intensity across the pupil as well as the WF (slide 16)

With starlight through a coronagraph, controlling DM1 and DM2 needs the
pupil's amplitude as well as its phase.  What the existing frames give,
per reading:

| reading | amplitude from the existing frames? | how |
|---|---|---|
| ZWFS linear L, exact I | no | one frame; the solve assumes the flat's amplitude |
| ZWFS stepped S | yes | its clear frame IS the pupil intensity; the three depths give the complex field E conj(b) |
| vector Zernike V | yes, with a solver change | two simultaneous images = two equations per pixel; with b iterated they give A and phi together (Doelman 2019's "phase and amplitude"); the present solver takes A as known -- extend and gate |
| pinhole P, P/SRI PF | yes | the phase-stepped solve is the complex field against a known reference (Dube 2024's "complex E field reconstruction") |
| interferometer four-step | yes | fringe modulation = the test amplitude times the reference's, the reference known |

So no separate camera is needed: every phase-stepped or two-image form
returns the complex field; the one-frame Zernike readings do not.  The
measurement to add: an amplitude aberration on the pupil (an apodizing
patch on the test-optic aperture, `macos.apodize`, 5% and 20% dips) and
each reading's recovered amplitude map against it, with the phase rows
unchanged as the gate.  CCL extends `solveV` (A and phi from the pair,
the 'amp' machinery reused) and gates it; TO gates P / PF; CCMac gates
the four-step.  On orbit the photon budget per stellar magnitude is
future work (10.3).

### 11.3 One bench, all modes (slide 17)

Yes: one layout switches between every mode by inserting or removing a
part, nothing realigned:

| switch | part | modes |
|---|---|---|
| the mask seat translates | one substrate with the etched dimples, the pinholes (with their attenuated surrounds), the metasurface and a clear window (the VSG2 nine-spot idea) | Zernike / vector Zernike / pinhole / clear (interferometer, phase retrieval) |
| reference-arm shutter | the interferometer's reference flat on its PZT stays built | interferometer on / sensors (arm shuttered) |
| quarter-wave plate in or out | the MacNeille cube and camera B stay behind the field lens; with the laser p-polarized the cube transmits (98%) to camera A in every scalar mode; the plate in makes the vector split | vector Zernike / all others |
| flip-in pickoff plate | the P/SRI's second arm (Lr1, pinhole, Lr2, folds, compensator, waveguide, BS3) on its own breadboard behind a flip-in plate | P/SRI / all others |
| in-arm quarter-wave plates in or out; analyzer | the polarization-snapshot form of the interferometer on the plate rig; the PZT form needs neither; the hybrid uses both | interferometer forms |

The v2 cemented-cube interferometer is the one form that does not switch
in (it replaces the plate splitter); it stays a backup slide.  The
universal-bench drawing: the two-arm interferometer with the vector tail
and the flip-in P/SRI arm, switchable parts marked -- TO draws it (the
`psri_bench` + `zwfs_vlayout` recipe, both theirs to combine), CCL QAs.
