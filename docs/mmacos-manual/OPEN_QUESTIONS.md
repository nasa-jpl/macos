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
- **ch. 7: `run_budget`'s per-group chart (`<name>_budget.png`) leaves out the rows that make its total.**
  It bars the alignment groups only; the actuator and sensing rows (2.595 + 1.281 + 1.251 + 0.092 nm of
  the 3.156 nm in its title) are not in it, and the controlled groups' bars are projection residuals
  (1e-9 .. 1e-12 nm, numerical zeros) drawn as results.  The manual prints it with a caption that says
  both.  No fix.

## Questions
- Authorship / acknowledgment block in the front matter: copied from the MACOS manual as a placeholder.
- **The thread's error budget (ch. 7):** e2e has no budget stage, and its group table (an allowed sigma per
  error line, the control set, the sensing and actuator noise) is a design choice, so chapter 7 runs the
  library's e5hex1 budget (the same seven-segment class) and says how `run_budget` takes the thread's
  `s4_sens.mat`.  Do you want an e2e budget spec authored (and by whom), so chapter 7 joins the thread?
