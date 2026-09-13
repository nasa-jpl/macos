# BRIEF for TO (Opus): the PDI's part of the DM Surface Gauge Comparison deck, and two shared loop knobs

From CCL for Dave, 2026-09-13.  Plan: `BRIEF_gauge_deck.md` (read it
first; sections 2-8).  Standing rules: dev-candidate, one MODEL-1024
MATLAB at a time (the batch wrappers serialize; check `pgrep -af
'seq.*sh'` before launching a sequence); `exit(0)` only in batch
wrappers; every number in a committed report; American English; photons
"per measurement".  You may edit `dm_gauge_lib/dmg_loop.m`,
`dmg_pdi_gauge.m`, `zwfs_run.m` (the loop and PDI paths; keep the
ideal-arm and vector paths bit-identical -- the v3dev regression is G4 =
0.296 pm at model 512 / NGRID 65 / grid 256 / 0.42, `dm_use 2`, `reg.mode
record`, readings V) and `tests/tDmgLoop.m`.  Commit into the shared
tree with `git add` of your files only; `git diff` a shared file for
foreign hunks before adding it.

## Dave's rulings on your three points

1. The P/SRI models the paper's pickoff form (your `psri_bench`); the
   stepped pinhole at the intermediate focus stays as reading P.
2. Pinhole diameter of record: whichever of 2.0 lam/D (1024 / 193) and
   1.0 lam/D (2048 / 385, `param_file 'macos_param_2048.txt'`,
   `ZWFS_MEMMAX=20G`) performs better on the 30 nm-surface rows and the
   loop; record both, state the choice with its numbers.
3. The PDI gets its own directory `templates/40_benches/pdi_dm96/`:
   `pdi_params.m`, `pdi_run.m` (calling `zwfs_run`), the psri bench
   decks and figures, `README.md` (the P / PF section moves there, with a
   one-line pointer left in zwfs_dm96/README.md), `runs/`.  Code stays
   shared through `zwfs_run` / `dm_gauge_lib`; nothing is copied.

## Deliverables, in order

### 1. PF through the two decks

Run the P/SRI reading with the reference physically traced through
`psri_ref.in`'s pinhole and the test through `psri_test.in` (your NEXT),
instead of the synthesized reference: the rows on the 30 nm surface
(matrix on it: single 10 nm, 1 nm on the 47 grid sites, dense 10 nm),
photons, the loop rows.  State the difference to the synthesized PF.

### 2. Capture range and photons, P and PF

The runner prints "capture range to 10%" after every ladder now.  Run
the aging ladder (matrix at 30 nm; 30 / 40 / 50 / 60 / 80 / 100 / 120 /
160 / 240 / 480 nm, `battery.ladder_sites 'grid'`) and the re-measured
rungs (`battery.base_rms` 60 / 90 / 120 / 160 with `calib_surface
'base'`, rows {'base/grid'}), and N(1 pm) at 30 / 60 / 120 / 160 nm
(`stages` bench + noise, `noise.nstates 10.^(8:2:14)`, `noise.nreal 6`)
-- the pattern of runs/cap385*, noise193_b* in zwfs_dm96.  P and PF,
385 rays for the ladders, 193 for the photons (the ZWFS record's
choices).

### 3. The pinhole-diameter choice (ruling 2)

The 1.0 lam/D pinhole at 2048 / 385 versus 2.0 at 1024 / 193: rows on the
30 nm surface and one loop run (`loop.nph [1e13 1e15]`, `loop.drifts
{'walk'}`, `loop.steps 1e-6`) for each; choose; record both.

### 4. Two shared loop knobs in `dmg_loop` (with tDmgLoop gates)

- `loop.start_rms` (mm): the loop starts from a surface of this rms
  (the same random field as the base, scaled) instead of the set-point;
  `loop.recal_every` (cycles; 0 = never): the response matrix is
  re-measured on the current surface every that many cycles (through the
  instrument's own calibration, `calib_matrix_`, on the loop's current
  DM command).  The DESCENT run: start 100 nm rms (200 nm WFE), matrix
  measured at the start, gain 0.5, recal every 10 and never, photons per
  cycle 1e13 and 1e15, K 60; report cycles to reach 10 nm and 3 pm, the
  final residual, converged or not -- readings L, S, V, P, PF.  CCMac
  mirrors the knobs in tg96_run; keep the interface to `ins` handles
  (add `ins.recal(cmd)` returning a new `est`).
- `loop.intra` (the V4 item): the DM / thermal drift ADVANCES between
  the frames of one stepped measurement, as your `cam_intra` does for
  the camera -- frame j of nf sees the drift at (j-1)/(nf-1) of the
  cycle's step; the simultaneous readings (L, I+, V, the polarization
  snapshot IFO) see none.  Run: walk 2 pm per cycle and thermal 5 pm per
  cycle at 1e13 / 1e15 photons per cycle, readings S, P, PF (and L, V as
  the controls), K 60; report ss, bias, the fixed error.  This is the
  IFO's PZT-form drift term too (CCMac runs it there).

### 5. The reference arm's own drift (P/SRI)

The non-common path: a phase walk of the reference arm relative to the
test arm, per cycle (a knob on the PF instrument: `pdi.ref_walk`, rad per
cycle rms, a random walk), in the loop at 1e13 / 1e15 with the walk
drift; report the floor it sets versus its size (1e-3, 1e-2, 1e-1 rad
per cycle).

### 6. Layouts and parts

Check `psri_layout.png` / `psri_render.png` / `pdi_layout*.png` against
the recipe in `zwfs_dm96/zwfs_vlayout.m` (read it): `macos.view_rx`
with `'labels', false` and `'hide'` on the passive planes, the fold plane
from above, names with leader lines off the beam, 15-17 pt type in an
1800-px figure, the crowded node as a cropped second panel.  Redo what
does not meet it.  Parts lists as tables (part, size / focal length /
AOI, coating, count, purpose): the pinhole substrate (diameter, surround
transmission, the stepping), and for the P/SRI the pickoff BS, Lr1 /
Lr2, the pinhole, folds, compensator, the waveguide chip and phase
shifter, BS3, camera.

### 7. deck_pdi conclusions; README

State the conclusions from the records above; move the P / PF README
section to `pdi_dm96/README.md` with a "Run it yourself" that reproduces
every number here.

## Report

One file `pdi_dm96/REPORT_gauge_pdi.md`: numbers first, departures
flagged, run tags, figures.  Say which of 1-7 is done, and push when
tDmgLoop and the fast suite are green.
