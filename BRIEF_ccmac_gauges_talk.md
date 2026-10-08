# BRIEF for CCMac — standing by for Dave's DM-gauge talk (Thursday 2026-10-08, 4 pm PDT)

You are starting from a clear context on Dave's Mac.  Dave runs the meeting from the
Mac and may ask you for a slide's provenance, an explanation behind a number, a figure
re-rendered, or a short live demo.  This note gives you everything to be useful in
that role without reading the whole arc first.  Read sections 1–3 before the meeting;
4–6 as needed.

## 1. The deck, and getting current

**The talk's deck is `macos/demo_session/deck_gauges_short.pptx` — FINAL** (Dave's
sign-off 2026-10-08 morning; DRAFT marks removed; export mark on).  54 slides: title +
21 main + a Backup divider + 31 backup.  Its source is `deck_gauges_short.md` with the
layout sidecar `deck_gauges_short.geo.json`; `deck_gauges_short_edit.pptx` is Dave's
editing copy (identical to the build as of the last sync; `_baseline.pptx` is the
reference the sync script diffs against).

The full deck `deck_gauges.md` / `.pptx` (the long version the short one was cut from)
is still DRAFT-marked and is NOT the talk; it is where the deeper material lives.
`deck_gauges_material.md` is the raw material behind every figure (pixel sizes, run
tags).  Figures are in `demo_session/figs/` (164 committed PNGs); the panel crops are
produced by `demo_session/crop_panels.py` from the runners' own outputs — never
re-rendered, only cropped/recomposed (Dave's rule: deck figures are the tools' own
output files).

Everything for the talk is pushed, both repos, branch `dev-candidate`:

```bash
cd ~/dev/macos && git pull            # deck_gauges_short at b69c974 or later
cd ~/dev/MACOS_resources && git pull  # runners, records, lane reports
```

If you need to rebuild or edit the deck: `cd ~/dev/macos/demo_session && python3
make_brief_slides.py deck_gauges_short.md` (python-pptx), then
`./sync_edit_deck.sh deck_gauges_short` (it DIFFS Dave's edit copy first and refuses
if he has unfolded edits; it keeps a `.bak`).  Dave's standing rules for decks:
diff the edit copy before every sync; American English; no new figures unless they are
a runner's output; numbers only from committed run records.

If you run the engine: pull brings engine commits since your last build (the grating
aperture fix, the forward root pick, spcOption, stop_info_set range — none touch the
gauge benches, but the mex must match the tree).  Rebuild per
`macos/RUNBOOK_mac_cycle3b.md` "Setup" (your own build/mex recipe on the Mac), then
`./run_mmacos_tests.sh tDmgLoop` as the smoke test (15/15 expected).

## 2. What the talk says (so you can answer "where does that come from")

**One optical bench, two front ends (lens rig and off-axis-mirror rig), four ways of
sensing the figure of a 96-actuator DM to picometers:** a Twyman-Green interferometer
(hybrid: polarization snapshot for change measurements, a piezo four-step for the
absolute-calibration leg), a scalar Zernike sensor (stepped quarter-wave dimple), a
vector Zernike sensor (geometric-phase metasurface, two circular polarizations at ±90°,
two cameras), and a point-diffraction sensor (stepped 5.3 µm pinhole with an attenuated
surround and a shutter frame).  Scored for accuracy, precision, repeatability and hold,
photon-noise driven; three jobs (measure the 30 nm working surface to pm, capture the
100–200 nm post-launch figure, hold 3 pm in a servo).

**Headline numbers (slide 2, the interferometer on both rigs):**
- accuracy: a 10 nm single-actuator change through the matrix measured on the 30 nm
  surface reads at gain **1.004** (lens) / **0.986** (mirror), floors 1.4 / 2.3 pm; the
  raw 30 nm map against the engine's field: polarization snapshot **1.3 %** (392 pm,
  gain 0.993) / **1.8 %** (551 pm, gain 0.987); the piezo four-step **0.01 pm, gain
  1.0000** on both rigs; the pupil imaging's share 0.13 % / 0.29 %.
- precision: 2.3 / 0.7 pm at 1e14 / 1e15 photons; repeatability in a servo 1.4 / 0.5 pm;
  hold 3 pm at 2.3e13 (lens) / 2.0e13 (mirror) photons per cycle, 7.1e13 / 6.6e13 under
  a 2 pm/actuator/cycle random walk.
- the three sensors (mirror rig) read within 1 % and 2 pm of the interferometer
  (slide "Performance side by side").  Vector Zernike needs 4.6e13 photons for 1 pm,
  stepped Zernike 1.2e14, pinhole 1.5e14 ("Photons for 1 pm").  Servo at gain 0.5: the
  vector Zernike holds 3 pm from 2.2e12 photons per cycle, the pinhole from 3.1e12.
- capture: from 100 nm rms only the externally referenced readings converge (the
  interferometer; the stepped pinhole if a second color is added for range).
- recommendation: hold with the vector Zernike; capture with the interferometer (or the
  stepped pinhole); the technology to develop is listed on the last main slide.

**Scope statement (asked for by Dave, on slide 2's footnote):** gain/floor rows are
single noiseless runs (systematics only); photon rows are Monte Carlo, 6–24
realizations per level; servo rows one 60-cycle run per level, rms over the last 30;
one drift realization.  Ideal optics, perfect camera, photon noise, mirror and camera
drift models only — no vibration, thermal, detector systematics (Steeves et al. 2020,
Optica 7, 1267: 1.6 pm in 4.3 s on a real interferometer, within 2× of our
photon-limited floor).

**Provenance:** the backup slide "Provenance" lists every run tag per approach; the
lane reports are `MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/
REPORT_gauge_ifo.md` (+ `REPORT_oap.md`, `REPORT_bench_realism.md`,
`REPORT_reflective.md`), `zwfs_dm96/README.md` + `deck_zwfs.md`,
`pdi_dm96/REPORT_gauge_pdi.md`; the plan and rulings `macos/BRIEF_gauge_deck.md`
sections 5–11; the pictures' raw material `demo_session/deck_gauges_material.md`.

## 3. The explanations that came up while the deck was finished (likely questions)

Each of these was measured, not argued; the runs are in the records named.

1. **"Why does the scalar Zernike's raw map look so bad?"** It is a GAIN error
   (0.78), not noise: the reading's reference amplitude is taken from the flat, so on
   the 30 nm surface the linear reading scales low; through the matrix measured on the
   surface it reads correctly.  Slide "Each gauge reading the 30 nm working surface,
   station by station" says so in its last column; `zwfs_dm96/runs/stations193_oap`.
2. **"The interferometer station figure had a ring around the estimate"** (gone in
   the final deck): the shared interference-frame mask `ctx.msk` includes a
   reference-only annulus (the reference arm's cone is wider than the test arm's); the
   figure now masks on the TEST arm's support.  The runner's own pixel-rms metrics still
   carry the annulus — an open item (PLAN_CONSOLIDATION 1.4), diluting them ~23 %;
   actuator-space rows are unaffected.  If asked: the deck's numbers in actuator space
   stand; the pixel-rms rows are conservative.
3. **"Why a hybrid interferometer?"** The polarization snapshot reads every change in
   one frame (no drift within the measurement) but carries the splitter's diattenuation
   (the 0.7–1.3 % gain); the piezo four-step is slow but reads the surface at gain
   1.0000 — so it is the absolute-calibration leg.  `fourstep_pzt` in `tg96_run.m`;
   the "51 nm fold" that once appeared on the lens rig was the four-step's SIGN, now
   set by a +push calibration.
4. **"Pupil image distortion and blur — what does it cost?"** After calibration 0.13 %
   (42 pm) on lenses, 0.29 % (92 pm) on mirrors, piston and tilt removed; 0.5 % / 0.9 %
   in the top spatial band.  The DM is the stop; the detector leg's PSF zone by zone
   from rays (`tg96_pupilsim`, `tg96_pupilq`).
5. **"How is the pinhole stepped / the dimple stepped / the snapshot implemented?"**
   All on the slides (the gauge slides have Reads / Frames / Reconstruction / Gets
   right / Gets wrong / Scored).  The snapshot: input polarizer at 45°, a double-passed
   quarter-wave plate in each arm returns the arms in orthogonal polarizations, an
   output quarter-wave plate turns their phase difference into the orientation of a
   linear state, and a polarization camera (micro-polarizer array 0/45/90/135° in 2 × 2
   super-pixels) reads the four steps in one exposure; the piezo form steps the flat in
   time by λ/8 of surface.  The stepped dimple and pinhole translate the mask seat.
6. **"Linear vs exact Zernike reading?"** Slide "Zernike sensor": the linear reading
   imprints error under a slow drift; the exact reading re-computes the reference wave
   (the first import from the literature, backup "What the Zernike-sensor literature
   offers").
7. **"Thermal ramp"**: a 5 pm/cycle low-order growth of the mirror's figure, the
   model's stand-in for thermal drift of the DM (no temperature anywhere in the model);
   the loop lags at rate ÷ (gain × response), ~10 pm regardless of photons — a cadence
   and thermal-control requirement, not a sensing one.
8. **"Two DMs with one sensor"**: slide "Controlling two DMs with one sensor: the full
   complex amplitude" cites Redding et al. 2023 (Integrated active control…, Optics +
   Photonics 2023) and 2024 (HWO WFSC, Proc. SPIE 13092).  The paper numbers were
   inferred from Dave's `Documents/dcr26` folder — if someone asks for the exact
   citation, check the PDFs there before quoting.
9. **Capture range (DAVE'S RULE):** devices do not run at null; every reading's
   working-surface range is recorded to 10 % accuracy (slide "Capture range").
10. **Calibration vs assessment (slide 19's title, Dave's wording):** calibration
    removes the systematic errors; assessment covers those it cannot.

## 4. Live demo candidates (time them BEFORE the meeting on the Mac; do not improvise)

Cheap and visual, in order of safety:
- **Layouts**: `tg96_run('stages',{'bench'})` / `zwfs_run('stages',{'bench'})` /
  `pdi_run('stages',{'bench'})` draw the benches from the parameter sheets
  (`tg96_params.m`, `zwfs_params.m`, `pdi_params.m`); the three-leg colored layout
  figures (red source→BS→flat, blue BS→DM, green BS→camera) are `dmg_leg_draw`.
- **A station figure** (one DM command through a train; the deck's row-2 figures):
  there is no separate stage — `tg96_run`'s `'figs'` stage draws it after `'bench'`
  and `'battery'` (`stations_ifo_`, the battery's lit set), and `zwfs_run` draws one
  per reading in S/V/P when `P.figs.stations` is on (`stations_fig_`).  The deck's
  own figures came from the run tags `redo96_lensstn2`, `redo96_oapstn2`
  (`tg_psi_dm96_oap/runs/`) and `stations193_oap` (`zwfs_dm96/runs/`): their `.log`
  files carry the exact command line and the wall time — reproduce from those, and
  expect minutes to tens of minutes at model 1024, not seconds.
- **The mask round trip / reference-wave gates**: each stage asserts its own gates and
  prints them — a safe thing to show "the model checks itself".

NOT live: the loop, noise and descent runs (hours each; the model-1024 runs need
11 GB and the batch wrappers run ONE at a time — `tg96_batch.sh`, `zwfs_batch.sh`);
anything at 2048/385.  The "What comes next, and how to run it" slide lists the full
command set per approach; each template README has the run-tag index.

## 5. Background reading, in priority order

1. `deck_gauges_short.md` itself (the text IS the record; read the footnotes).
2. `macos/BRIEF_gauge_deck.md` sections 5–11 (plan, rulings, what was cut and why).
3. `tg_psi_dm96_oap/REPORT_gauge_ifo.md` and `REPORT_reflective.md` (your own lane's
   reports — the bench, the 22.5° node, the reflective front end).
4. `zwfs_dm96/README.md` + `deck_zwfs.md`; `pdi_dm96/REPORT_gauge_pdi.md`.
5. `deck_gauges.md` (the long deck) for anything the short deck compresses.
6. The shared library `dm_gauge_lib/` (readings, loop, noise, capture; `dmg_loop`,
   `dmg_arm_maps`) when a "how is X computed" question goes below the runner.

## 6. Rules that apply if you touch anything

- Push nothing; commit locally with the SHA and branch stated; Dave says when.
- Never relink a shared mex or rebuild a shared tree while another MATLAB runs; on
  the Mac that is your own MATLAB only — still one model-1024 MATLAB at a time.
- Deck figures: the tools' own outputs, cropped only; fix the producer, not the PNG.
- Every number quoted comes from a committed run record; if a question needs a number
  that is not on a slide or in a report, say so rather than estimate.
- American English; Dave's rules on decks are in the agent memory on the Linux box
  (`feedback_deck_*`, `feedback_capture_range`, `feedback_deck_best_results_only`) and
  summarized above.

— CC (Linux), 2026-10-08, for Dave.
