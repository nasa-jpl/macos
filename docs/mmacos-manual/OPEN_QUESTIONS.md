# mmacos manual -- open questions and findings for Dave (TO, from 2026-10-23)

The report on Dave's return gathers this list (BRIEF_to_mmacos_manual.md, Report).

## Defaults taken so the work could start (overturn in one line)
1. **The thread's telescope = e2e** (`templates/80_end_to_end/e2e`, the 4 m f/18 folded Korsch TMA, all
   seven stages).  The brief first described it as "the 3 m f/20 TMA", which is e2e2 (design stages
   only); CC corrected the brief 2026-10-23.  e2e2 is chapter 1's pointer to the improved design flow.
2. **The forty uncalled helpers** are documented in one table (Part II), not deleted.

## Findings met while writing (no library change made; PLAN_CONSOLIDATION 1.2 if Dave wants them fixed)
- **ch. 1: the e2e template's committed stage-1 outputs are stale.**  A fresh run of the unchanged
  `s1_telescope.m` gives conics K = [-0.98972 -1.51153 -0.97231] (committed report: -0.98946 -1.50791
  -0.98712), worst +-1' WFE 0.0138 waves (committed 0.0158), M1 hole r 0.232 m (committed 0.095 m --
  the `hole_min_r_m` floor in `e2e_params.m` postdates the committed report).  The manual prints the
  live run; the template's `s1_report.txt` / `.in` / `.mat` / PNGs are from August.  Regenerating the
  committed template outputs is Dave's call (stages 2-7 consume them).
- **ch. 1: `design_report`'s "EFL" is the local plate scale at the biased field point** (chief-ray
  image displacement per field angle, +-0.1' about the 2' bias): 72.32 m, vs the engine's SYSPROP EFL
  72.14 m on the emitted deck.  Both are right for what they measure; the report's label could say
  "plate-scale EFL at the field point".

- **ch. 2: `templates/30_instruments/coro_walkthrough` no longer shows a coronagraph** (Rx_Coro_FPM.in,
  model 1024, run unchanged 2026-10-23).  The radial-contrast curve with the coronagraph is the no-mask
  curve's ring structure scaled down uniformly (~5e-8 at 0 lambda/D, dark-zone mean 1.49e-09); the
  FPM-plane panel is at the floor everywhere and the Lyot pupil is a filled disk, not the edge ring the
  README describes.  That is the signature of a focal mask covering the whole computed field (leakage
  ~1e-4 in amplitude), not of a dark hole.  Its figures are gitignored, so no committed picture shows
  what it drew before -- a producer that rotted silently, or a scaled deck that never worked; not
  diagnosed (no library change from this brief).  Chapter 2 uses `coro_planet_demo` (CoroExample, the
  CLI journal's own recipe) instead.
- **ch. 2: `coro_planet_demo`'s planet page no longer shows a planet** (defaults, model 512, run unchanged
  2026-10-23).  Asked for a planet at sky-angle cosines [+2e-5 +2e-5], the engine puts the planet's
  chief ray at x = -2.0e-5, y = +2.0e-5; "planet light, by itself" is a shifted, clipped speckle
  field filling one corner of the composite, not a PSF (RMS spot 197x the star's; peak at the grid edge,
  [-21 -43] px, centroid [+9.9 -9.9] px); the planet peak equals the star-residual peak to the digit
  (2.13e+03), and the composite grid is 91x91 where the help text's default `npix` is 64.  The
  surfaces page (pupil, FPM, Lyot edge ring, suppressed star) is sound and is what chapter 2 prints.
  The demo was pushed working (memory: the planet-injection redo); not diagnosed.  The index's T3
  coronagraph entry still names `coro_walkthrough` (above).
- **ch. 2: ImageDemo.in declares no BaseUnits**, so the engine's metre conversions use CBM = 1 and
  `first_order_properties().dpup_m` reads 380 (the deck's millimetres).  The camera template
  (PLAN_MMACOS_MANUAL, the ImageDemo note) should declare `BaseUnits= mm`.
- **ch. 3: the live pie clips 0 gap rays; the committed template run clipped 376.**  Live (s3 unchanged,
  on the live s2): 12128 source rays, 11950 pass bare and 11950 with apertures, first-fail element 1
  (the center segment's hole); committed `e2e_pie_seg_report.txt`: 12520 source rays, 12474 bare,
  12098 with apertures (376 clip), first-fail element 18.  The live footprint figure shows the gaps
  empty and every footprint inside its polygon, so the zero is no ray reaching a gap, not a missing
  aperture; why the source grid changed (12520 -> 12128; the parent 1305 -> 1185 rays) is not
  diagnosed -- the stale stage-1 outputs (above; the hole 0.095 -> 0.232 m) are the first suspect.
  Stages 4-7's committed products were built on the 376-clip deck.
- **ch. 5: e2e's committed `s5_report.txt` did not come from the committed `s5_met.m`** (PLAN_CONSOLIDATION
  1.2 (ggg)).  It says "layout held from preset (no search)" (run_met's `opt_preset`); the script passes
  none and run_met defaults to `optimize = true`.  The live search: center-hexagon class angles
  [170 190 50 70 -70 -50] deg vs the preset's [170 190 70 90 -90 -70] (wedges identical); optimized
  WEM MET-only 3.877 vs 3.707, edge+MET 1.782 vs 1.755.  The manual prints the live search, run on a
  copy of `s5_met.m` with its preset-library promotion cut (presets md5-identical before and after).
- **ch. 5: `e2e_pie_met_metric.png` drops a bar and `e2e_pie_met*_layout.png` mislabels an axis.**  The
  metric chart's log axis bottoms near 0.09, so the "edge+MET opt" bar (about 0.04 of prior by the WEM ratio; the chart does not print it) does not show; the
  hub-plane inset of the layout figures is labeled mm but plots metres (fiducials at r = 0.196).  The
  manual prints the layouts (the caption says the units) and leaves the metric chart out.  PLAN_CONSOLIDATION
  1.2 (hhh), (iii); the label is hard-coded at src/+macos/+design/met_view.m:199-200.
- **ch. 6: `run_compare`'s engine-vs-linear gate mixes two OPD references** (PLAN_CONSOLIDATION 1.2 (kkk)).
  The engine's OPD default is 'chief' since 2026-10-10; `run_sensitivities` harvests at 'mean'; the
  compare's engine truth takes the default.  Live worst w rel 0.945 / 0.990 / 0.987 (x / z / grid) vs the
  committed August 2.6e-3 / 1.8e-3 / 0.29.  Measured (ch06c): Seg 1 Tz +100 nm, engine 'chief' 188.9 nm
  rms (178.4 nm piston), 'mean' 61.98, the linear column 61.98.  Which reference the harvest uses is
  Dave's ruling.
- **ch. 7: `run_budget`'s per-group chart (`<name>_budget.png`) leaves out the rows that make its total.**
  It bars the alignment groups only; the actuator and sensing rows (2.595 + 1.281 + 1.251 + 0.092 nm of
  the 3.156 nm in its title) are not in it, and the controlled groups' bars are projection residuals
  (1e-9 .. 1e-12 nm, numerical zeros) drawn as results.  The manual prints it with a caption that says
  both.  No fix.
- **Part II: `tma_unobscured`'s closing design report reads NaN on every field ring** (step [6], after
  `add_pupil`; run unchanged 2026-10-23, the finder then the demo), while the template's own field ladder
  ([3]) is diffraction-limited to the 1' ring (0.032 waves, tilt removed).  `freeform_unobscured`'s report,
  same path, gives numbers.  Not diagnosed; the page prints the steps before [6] and says so.
- **Part II: `tma_onaxis`'s README table predates the engine it runs on** (1' bias: 0.0074 waves in the README,
  0.0059 live; the whole sweep lower by about a fifth, same shape).  The page prints the live sweep.
- **Part II: `offset_imager()` with its defaults does not reproduce the rodgers3 ladder of record past S1.**
  Run unchanged 2026-10-23 (9841 s): S1 37.8 / 25.8 / 9.2 nm (map max / avg / std), exactly the record's
  (`challenges/rodgers3/t3/r3t_REPORT.md`); then S2 40824.4 (record 303585.5), S3 41849.9 (252.0), S4 28907.4
  (113.6), S5 12415.1 (118.2).  The record is `run_t3`: the exit beam pinned to [0 0 -1] at 0.1 deg and 30
  Gauss-Newton iterations a stage, against the defaults' report-only exit and 12.  The two part at S2 already
  (a focal-plane refit alone), so those differences are not shown to be the cause.  Not diagnosed; the README
  presents `offset_imager()` as "the five-stage ladder (rodgers3 instance)".  The page prints today's run.
- **Part II: `tg_psi_dm96_oap`'s layout figure (`<tag>_vlayout.png`) mislabels its colors.**  Its title reads
  "reference arm red, test arm (DM) blue, camera leg green"; the lens rig's figure (run 2026-10-23, bench stages)
  draws both arms red.  The figure-producers-rot class; the page's caption says so.
- **Part II: `e2e6m`'s telescope stage no longer reproduces its committed report.**  Run unchanged on a copy
  2026-10-23 (layout search, s1_telescope, s2_segmentation; 1095 s): the same all-sphere layout and S0 (0.8594
  waves), then S1 inner 0.3512 (committed 0.0200), S2 worst 0.0702 (0.0355), dense 7x7 raw max 0.1059 [FAIL]
  (0.0502 [PASS]), EFL 70.579 m f/11.76 (152.327 m f/25.39).  The template's files are unchanged since
  2026-08-25 (a9787a9); the library under them is not.  Not traced.  The same shape as offset_imager's (lll):
  the first stage reproduces, the next departs.
- **Tools appendix: `design_layer_api/example_sensitivities_from_rx.m` ends in `exit(0)`** (line 49), and its README
  says all four examples do; the library's own rule is no `exit` in an example (a batch wrapper supplies it).
  `example_telescope_design.m`, the one the page runs, does not.  One line to remove, yours to allow.
  **The two rules conflict, and that is the likely origin:** `mmacos/CLAUDE.md` ("`clear mmacos` hangs
  `matlab -batch` ... and so does implicit exit") says to end every batch script with `exit(0)`, because a
  batch MATLAB that has loaded the mex can hang at an implicit exit; your rule puts the exit in the
  `-batch` command instead, never in the example.  Section 10 follows your rule.  Which one the library's
  record should state is yours.
- **Part II: refusals today** (pages print them): `zwfs_run` / `pdi_run` bare call (coat_oap 'protectedAl' on the
  lens rig, zwfs_run.m:195, (nnn)); `e2e2` s1 conic-fixture gate (seidel_seed re-targeted 2026-10-03,
  the template's gate not, (ppp)); `rc_offaxis` final clipping assert (M2 cuts a foreign beam, (rrr)).
  `tg_psi_dm` / `tg_psi_dm_v2`'s registration gates (intensity zero after load_rx with POLARIZATION in
  force, (ooo)) pass since the engine fix macos c6c83d55 (CCL, 2026-10-24); their pages print the runs.

## Questions
- Authorship / acknowledgment block in the front matter: copied from the MACOS manual as a placeholder.
- **The thread's error budget (ch. 7):** e2e has no budget stage, and its group table (an allowed sigma per
  error line, the control set, the sensing and actuator noise) is a design choice, so chapter 7 runs the
  library's e5hex1 budget (the same seven-segment class) and says how `run_budget` takes the thread's
  `s4_sens.mat`.  Do you want an e2e budget spec authored (and by whom), so chapter 7 joins the thread?
- **Part II, the READMEs the 15 un-README'd design templates lack:** the manual keeps their generated text in
  `docs/mmacos-manual/part2/<name>.md` (from the entry script's header and its run); writing them into the
  templates themselves is a library change -- yours to allow.
