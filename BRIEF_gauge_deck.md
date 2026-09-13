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
15. Run it yourself.

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
