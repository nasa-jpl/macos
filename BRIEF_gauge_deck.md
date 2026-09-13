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

## 5. Deck outline (deck_gauges.md, ~18 slides, DRAFT until Dave signs)

1. Title.
2. The requirement: a DM surface gauge at picometers, holding a real
   surface (30-60 nm rms), 96x96, in a servo -- not at null.  The four
   numbers every candidate is scored on (rows, capture range, photons,
   hold).
3. The six configurations at a glance: one small layout each + the parts
   each adds to the shared front end.
4. The shared front end: layout + parts (source to the tail).
5-10. One slide per configuration: layout large (left), parts list
   (right), its numbers row (bottom).  A and B may share a slide if the
   OAP rig is dropped as a candidate (decision 2).
11. Performance side by side: the table of section 3 in slide form.
12. Capture range (both ways) across all six.
13. Closed loop: hold vs photons, all six on one figure (the dmg_loop
   figures exist per lane; the combined figure is a runner-figs job).
14. Systematics priced per configuration, and the ones still open.
15. What each lane has not finished (honest status).
16. The recommendation (Dave's; a draft stance is in section 7).
17. Run it yourself: one runner per lane, one scoring library.
18. Provenance.

Style: American English; figures = the tools' own PNGs, crop only;
every axis quantity defined once; photons per measurement, never per
state; the pre-write gate in doc/STYLE_REPORTS.md.

## 6. Work split, in order

**CCMac (reflective + lens IFO completion):**
0. The rows on the 30 nm working surface, matrix measured on it, with the
   47-site 1 nm grid row (the ZWFS convention; the current rows are on a
   16 nm base and single-site) -- lens and OAP.
1. The runner prints the capture range (port the `zwfs_run` line: the
   largest base rms with the 47-site 10 nm gain in 0.9..1.1), and runs
   the ladder both ways (matrix aging from 30 nm; re-measured at 60 / 90
   / 120 / 160) with N(1 pm) at each surface, lens and OAP.
2. Photons for 1 pm in the S5 noise-stage form (not only the loop's
   sig_n), lens and OAP.
3. Layout figures in the recipe: lens rig with the reference arm and the
   PZT flat, OAP rig with the fold plane in view, BS-node inset, names.
4. Parts lists (section 4).
5. Close the reflective design: is the 0.18 fold cross-talk that walls
   the loop reducible (OAP AOI, coating retardance spec, the waveplate
   azimuth) or is the OAP rig open-loop-only?  One paragraph with a
   number, in REPORT_oap.md.

**TO (PDI completion):**
1. PF run through the two decks (the reference physically through
   psri_ref.in's pinhole) instead of synthesized -- the P/SRI's number of
   record.
2. Capture range for P and PF: a ladder run (the print is in zwfs_run
   already), aging and re-measured, with photons.
3. The reference arm's own drift (non-common path) in the loop -- the
   term the paper lists; a knob in dmg_loop or the PF instrument.
4. Layout figures checked against the recipe; parts lists (section 4).
5. deck_pdi conclusions stated.

**CCL (assembly + the vector sensor's last term):**
1. The cube's channel leakage for the vZWFS (the coated diagonal is in
   the decks; a crosstalk term between the two images) -- the one open
   systematics on D.
2. The shared front-end parts list; the combined capture-range and loop
   figures across lanes (runner figs jobs).
3. deck_gauges.md + build; the side-by-side table; QA renders.

Sequencing: lanes 1-2 of CCMac and TO are the numbers the comparison
slides need; their layout/parts items feed slides 4-10; CCL assembles
once the numbers land and fills the rest from the existing records.

## 7. Decisions for Dave

1. Audience and title of the deck (Fang's bench build? the JPL review?).
2. Is the reflective (OAP) IFO still a candidate, given it cannot hold a
   2 pm walk in closed loop at any photon count?  Keep as a row with the
   verdict, or drop to a backup slide.
3. Which IFO phase-shift form is the candidate: the PZT four-step, the
   polarization snapshot (v1 plate / v2 cube), or the hybrid (snapshot
   for the change measurements, PZT as the calibrator) the ZWFS deck
   recommends.
4. TO's three points: which paper/layout the P/SRI models (pickoff vs
   intermediate-focus pinhole); the pinhole diameter of record (2.0 lam/D
   = the dimple, or 1.0 at MODEL 2048); whether the PDI gets its own
   directory or stays as readings of zwfs_dm96.
5. The recommendation's stance.  Draft, for your judgment: build the
   scalar Zernike sensor first (the mask substrate is the only new part,
   and its stepped reading matches the vector pair on every closed-loop
   line at four frames); the vector sensor as the upgrade (metasurface,
   QWP, cube, second camera: two frames, the best reading, no fold, 70 nm
   capture range aging); the lens IFO with a PZT as the absolute
   calibrator (its capture range reaches its wrap near 150 nm and it
   never folds, at 2.6x the light and a 9% fixed error in hold mode);
   the P/SRI where a surface-independent reference is worth a second
   arm; the OAP rig not for hold mode.
