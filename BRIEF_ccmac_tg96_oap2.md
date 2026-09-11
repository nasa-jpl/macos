# BRIEF for CCMac: tg96-oap, the last mile -- windows from the ray affine, calibration by the measured response matrix

From CCL for Dave, 2026-09-11.  Written for a CLEARED or compacted
session: everything you need is in this file and the files it names;
pull `origin/dev-candidate` (both repos) and `origin/tg96-oap`
(MACOS_resources) first.  Your report `tg_psi_dm96_oap/REPORT_oap.md`
was reviewed 2026-09-11: the feasibility result stands, the lens gate
is exact (re-run here: tBench 9/9 incl. `test_twyman_green_optics`),
and the route for the remaining item is changed -- see "The route".

## Where the branch stands (facts, verified 2026-09-11)

- `origin/tg96-oap` = `b6bb3fe`, seven commits on top of `fb6fb98`
  (2e56522 builder; cfe7a3b runner; 8bb7d8e tail retuner + lens gate;
  9098182 report; b884dcd / 4ec2e4a in-pupil poke placement; b6bb3fe
  artifacts).  It does NOT carry what landed on dev-candidate after
  fb6fb98 (TO's dwd* plotter commits; ZWFS S9/S10); it touches only
  `src/+macos/+design/twyman_green.m`, `tests/tBench.m` and
  `templates/40_benches/tg_psi_dm96_oap/`, so it merges cleanly.  Start
  by merging `origin/dev-candidate` INTO `tg96-oap` (you need the ZWFS
  S10 runner code named below).
- Validated: `twyman_green('optics','oap')`; 'lens' byte-identical;
  `tg_psi_dm96_oap/tg96_{params,run,tail,run_batch}.m`; the lens rig
  through the runner reproduces `tg_psi_dm96/tg96_report.txt` line for
  line at model 1024.
- The reflective rig IMAGES the DM as well as the lens: in-pupil
  actuators recover 129-143 nm against the lens's 131-146 (a poke at the
  exact centre reads 0 because the four-step map is referenced there --
  see "Reference" below); `macos.pupil_quality` says the OAP exit pupil
  is cleaner than the lens's (|astig| 0.144 vs 1.025); fold AOIs 5 / 9
  deg, margins printed.  The flat-DM null is 12.9 nm (lens 0.134): a
  low-order arm difference of the same-plane fold that the common tail
  cannot null -- the fold cost Dave asked to measure; it should cancel
  in the differential rows and must be reported beside them.
- OPEN: the detector -> DM pixel mapping of the folded rig (a flip, a
  rotation, and 0.8x the lens's scale; illuminated mask 18376 vs 28917
  px) defeats `tg96`'s inline two-poke registration (|corr| 0.0014),
  which gates the closure / transfer / differential table.

## The route (replaces "swap in dmg_register")

`dmg_register` searches only the eight flips and transposes; it cannot
express a rotation that is not a multiple of 90 deg or a scale, so it
would not resolve this mapping either.  Two things resolve it, and both
are already on dev-candidate:

1. **Window placement from the full ray-traced affine.**
   `dm_gauge_lib/dmg_frame.m` already fits the complete 2x3 affine
   detector-mm -> DM-mm over the traced rays
   (`Aaf = [xy_d(:,okr).' ones(nnz(okr),1)] \ xy_to(:,okr).'`) and then
   returns only `mag = sqrt(|det|)`.  Extend it to return `Aaf` and the
   detector reference (the chief ray's pixel and `dxd_mm`), invert it,
   and every actuator's lattice point (x = lat(c), y = lat(r), the
   convention `dmg_anchor` uses: `tax = xg(tc)`) maps to a detector
   pixel (u, v) directly.  Flip, rotation, scale and shift are all in
   `Aaf`; no parity search, no blob matching.  Gate: the centre of mass
   of each poked actuator's response lands within 2 px of its predicted
   (u, v) for > 99 % of lit actuators on BOTH rigs (the lens rig is the
   control: its mapping is known).
2. **Calibration by the measured response matrix** (Dave 2026-09-10;
   the ZWFS default since S10).  Poke every 8th actuator in a sparse
   grid, step the grid through its 64 offsets so every lit actuator is
   poked once, cut each response from its own +/- half-step window at
   (u, v), assemble J (detector px x lit actuators), and estimate
   actuator changes by regularized least squares on J.  No single-site
   kernel, no shift-invariance, no modal correction; registration only
   PLACES the windows, and the columns carry the response wherever it
   lands (position dependence, a real DM's irregularities).  Reference
   implementation to lift: `zwfs_dm96/zwfs_run.m` -> `calib_matrix_`
   (build) and `est_matrix_` (solve); knobs `battery.matrix_step` 8,
   `battery.matrix_lam` 1e-3 (of the median column energy),
   `battery.matrix_sign` 'same' (alternating +/- measured neutral to
   slightly worse in the model; keep the option for bench drift).
   Replace only its window placement (it uses the ZWFS bench-stage
   anchor + parity) with route 1's (u, v).  Measured on the ZWFS with
   this: single actuator on the flat 0.994 / floor 4 pm (single-site
   kernel: 0.946 / 37 pm); response 0.99-1.08 at every spatial
   frequency with no correction.

**Reference -- read before building the matrix.**  A multiplexed frame
carries whatever the instrument's map reference does to a sum of pokes.
The ZWFS cannot see piston, so its frames are mean-referenced and the
pokes' shared pedestal sat in every window; cut naively, the estimator
over-responded 2-4x below 12 cycles/aperture.  The fix in `calib_matrix_`
is two lines: subtract the frame's median over the mask before cutting,
and carry each column's own volume spread over the mask as a rank-one
term (`J'J = Jl'Jl - v v'/A`, `J'm = Jl'm - v (1'm)/A`).  YOUR map is
referenced to the chief-ray / centre pixel (why the centre actuator
reads 0): a poke's value AT the reference pixel is subtracted from the
whole map, so the centre column would read zero and its neighbours
would carry a pedestal.  Do not special-case it: mean-reference every
map over the mask -- columns and measurements alike -- and the ZWFS
code then applies verbatim (piston is the one direction neither
instrument's differential should claim).  Say in the report which
reference the four-step map carries and that you mean-referenced it.

## Deliverables, in order (stop at the timebox, never thin gates)

1. `dmg_frame` returns the full affine (+ the chief pixel); window
   placement gate above on both rigs.  Non-vacuity: windows placed by
   the old lens-tuned mapping on the OAP rig must FAIL the gate.
2. `tg96_run` gains `battery.calib_mode` 'matrix' (default) | 'kernel'
   (the record), lifted from `calib_matrix_` / `est_matrix_`; the lens
   rig through it: report single-actuator gain and floor, the 12-mode
   response, and the Stage-E rows beside the record's kernel numbers
   (0.9654 / 21 pm flat single-10 nm; 0.9203 random) -- expect
   improvement, state what you measure.
3. The OAP rig through the same runner: Stage C-E rows + break ladder
   in actuator units, all pm, side by side with the lens rig; the
   12.9 nm null reported beside them and its cancellation in the
   differential measured, not assumed.
4. D4 alignment sensitivity (OAP1/OAP2 decenter 10 um, tilt 10 urad,
   one at a time; null and single-actuator row; pm per um, pm per urad;
   what the differential leaves).
5. STRETCH: D5 coated-Al OAPs (`coat_set`, the protected-Al anchor).
6. README section + REPORT_oap.md updated (numbers first; what departed
   from this brief and why); push `tg96-oap` when gate 1 is green on
   the lens rig; Dave orders the merge.

## Constraints (standing law)

No engine work.  'lens' stays byte-identical (tBench).  Two decks, two
traces.  Same-plane folds, field lens transmissive (Dave's rulings).
Conjugates not focal lengths (pole-to-focus = F1 / F2).  Functional
stops on flat marker planes.  Commits on `tg96-oap`; the report is the
record.  One MODEL-1024 MATLAB at a time (~11 GB); a trimmed size table
in the run dir (`zwfs_dm96/macos_param_2048.txt` pattern) if memory-
bound.  matlab -batch: `exit(0)` only in the batch wrapper.

## Cold-start reads

`macos/BRIEF_ccmac_tg96_oap.md` (the original brief + its addendum on
the stencil-site bias), `tg_psi_dm96_oap/REPORT_oap.md` and README,
`tg_psi_dm96/README.md` (the lens record), `zwfs_dm96/README.md` S9-S10
bullets and `zwfs_run.m` (`calib_matrix_`, `est_matrix_`, `frames_`),
`dm_gauge_lib/README.md`, `dmg_frame.m`, `dmg_anchor.m`, memory-style
summaries in `macos/CURRENT_SLICE.md`.
