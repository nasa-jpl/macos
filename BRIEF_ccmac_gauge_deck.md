# BRIEF for CCMac (Opus): the IFO's part of the DM Surface Gauge Comparison deck

From CCL for Dave, 2026-09-13.  Plan: `BRIEF_gauge_deck.md` (read it
first; sections 2-8).  Standing rules: work on dev-candidate; pull both
repos first (your render commits 5e45851 / b6ed9fe are merged here);
one MODEL-1024 MATLAB at a time on this box (the DM-gauge batch wrappers
serialize; check `pgrep -af 'seq.*sh'` for live sequences before
launching); `exit(0)` only in batch wrappers; the loop is `dmg_loop`,
never a copy; every number in a committed report; American English;
photons "per measurement", never "per state".  Do not edit files under
`zwfs_dm96/` or `dm_gauge_lib/` (CCL and TO own them); `tg96_run` may
call them.  Push when tBench is green and say so.

## What the deck needs from the IFO lanes (lens and OAP rigs)

Every item is a number with a run tag or a figure from the tool.  Order
matters: 1-3 feed the comparison slides, 4-6 the layout slides.

### 1. The rows on the 30 nm working surface (the ZWFS convention)

Your rows are on a 16 nm random base and single-site.  Re-run, lens and
OAP (bare Al), matrix measured ON the surface: base 30 nm rms (seed 7,
the same field the ZWFS uses: `P.battery.seed_base`), rows = single 10
nm at the hold-out site, 1 nm on the 47 grid sites (every 8th actuator),
dense random 10 nm; report gain / floor / SNR as `zwfs_run` does.  If
`tg96_run` lacks the 47-site row, add it (`dmg_lit` + the grid the ZWFS
uses, `Afig(4:8:end, 4:8:end)`).

### 2. Capture range, both ways, with photons

Port the `zwfs_run` line (search "capture range to 10%": the largest
base rms with the 47-site 10 nm gain in 0.9..1.1, log-interpolated
between rungs) into `tg96_run`'s ladder.  Run, lens and OAP:
- aging: matrix measured once on the 30 nm surface, ladder 30 / 40 / 50 /
  60 / 80 / 100 / 120 / 160 / 240 / 480 nm;
- re-measured: the matrix measured on 60 / 90 / 120 / 160 nm surfaces,
  the 1 nm grid row on each;
- photons: N(1 pm) in the S5 noise-stage form (sigma ~ c/sqrt(N) fit,
  single 10 nm on the surface, matrix on the surface) at 30 / 60 / 120 /
  160 nm.  The loop's sig_n is not this number; state both.

### 3. The three phase-shift forms

The four frames are the same in the model; the forms differ in what
they get wrong.  Price each with existing machinery, lens rig, loop
stage where a loop is named:
- PZT four-step: a step-size error of 2% and 5% (the frames with the
  error, the solve with the nominal steps -- TO's `pdi.step_err`
  pattern; add the knob to tg96_params) on the surface rows and the
  loop; drift between the four frames -- the camera walk within a scan
  (`dmg_loop` 'cam', `cam_intra`) and, once TO lands it, the DM drift
  within a scan (`loop.intra`); say what a 1e-3-of-signal camera walk
  and a 2 pm-per-cycle DM walk do to the four-step in hold mode.
- Polarization snapshot: all four frames at once, so no within-scan
  drift; the systematics are polarization: the v1 plate rig's BS
  diattenuation (11.7% high, the analyzer sweep fixes it), the v2 cube's
  R_p 2.1%, the coating retardance (item B).  State each as a gain error
  or a floor, from the runs that exist (tg_psi_dm, tg_psi_dm_v2), on the
  surface rows.
- Hybrid: snapshot for the change measurements, PZT for the absolute
  calibration -- one paragraph: which errors each half removes; no new
  run unless one is needed to state a number.

### 4. Lenses vs OAPs, and the OAP front end under the other gauges

- The IFO on each (you have it): the open-loop rows, the loop rows, the
  alignment sensitivity (D4), the coating story (item B).  Close the
  reflective design in one paragraph with a number: is the 0.18 fold
  cross-talk that walls hold mode reducible (OAP AOI, coating retardance
  spec, waveplate azimuth), or is the OAP IFO open-loop-only?
- The other gauges on the OAP rig: run `zwfs_run` on your OAP test arm
  (`bench.optics 'oap'`, `bench.coat_oap 'bareAl'`, plus whatever tail
  tuning your rig needs -- `zwfs_run` passes every `P.bench.*` field to
  `twyman_green`), readings {'L','S','V','P','PF'}, stages bench +
  battery (rows on the 30 nm surface, matrix on it) + noise + loop, and
  for V the arm maps on: `mask.v_arm 'engine'` (the engine's polarized
  traces of YOUR arm: the OAPs' retardance and diattenuation per
  channel; `mask.v_laser_deg` 45 and 90).  The bench stage's gates G1
  (mask sandwich round trip) and G3 (DM-conjugate pupil) may fail on an
  OAP tail; if they do, report the numbers and stop -- do not force them.
  The deck's slide 15 is this run's rows beside the lens rig's.

### 5. The descent run (capturing the initial figure)

Once TO lands `loop.start_rms` / `loop.recal_every` in `dmg_loop`
(brief for TO, item 4), mirror them in `tg96_run`'s loop stage and run:
start 100 nm rms (200 nm WFE), matrix measured at the start, gain 0.5,
recalibrate every 10 cycles and never, photons per cycle 1e13 and 1e15,
K 60; report cycles to reach 10 nm and 3 pm, the final residual, and
whether it converged, lens and OAP, PZT form.

### 6. Layout figures and parts lists

Your renders are not of deck quality (Dave): the OAP table-plane panel
views the fold plane edge-on, elements are labeled by E-number, passive
planes crowd the nodes.  Redo in the recipe of
`zwfs_dm96/zwfs_vlayout.m` (read it): `macos.view_rx` with `'labels',
false` and `'hide'` on the Reference / FocalPlane planes, the fold plane
seen from above (`view(ax, 0, 90)`, axis equal), elements named by a
`text` with a leader line placed off the beam, type 15-17 pt in an
1800-px figure, and the crowded node (BS + compensator; the OAP folds) as
a second panel cropped to it.  Draw: the lens rig with the reference arm
and the PZT flat; the OAP rig; the v2 cube rig if its deck still emits.
Then the parts list per rig as a table (part, size / focal length / AOI,
coating, count, what it is for), from the bench parameters -- lenses,
BS plate + compensator (or the cube), reference flat + PZT, OAPs
(off-axis distance, f, AOI, coating), polarizers / QWPs / analyzer for
the snapshot forms, camera (385 px per pupil).  Put them in
`tg_psi_dm96_oap/README.md` and REPORT_oap.md; CCL folds them into the
deck.

## Report

Dave (2026-09-13): the deck's main body shows each approach ONCE, in
its best-performing configuration; everything else is backup.  Open the
report with that configuration named and the numbers that make it the
best; keep the route there for the backup slides.

One report file (REPORT_gauge_ifo.md in tg_psi_dm96_oap/), numbers first,
departures flagged, every figure and run tag listed, the three
phase-shift forms in one table, the OAP-front-end rows beside the lens
rows.  Commit the pruned run artifacts as before.
