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

## Addendum 2026-09-11: the closed-loop HOLD metric -- the IFO half (deliverable 7)

Dave: on orbit the DM surface must hold to << 10 pm under frequent
remeasurement and closed-loop actuator servo; the metric that compares
the gauges in THAT mode is the steady-state hold error.  Spec:
`macos/BRIEF_loop_metric.md`.  The loop code is SHARED and already
gated: `dm_gauge_lib/dmg_loop.m` (read its header: the instrument
interface is four function handles + a lit mask) with `tests/tDmgLoop.m`
(8 gates on a synthetic linear instrument: geometric convergence at
1 - gG, the noise-only steady state sigma_n sqrt(g/(2-g)), the random-
walk and ramp laws, a biased reading converging to a non-zero surface,
seeded drift, a single-shot reference as a fixed bias, Parseval).  Do
NOT copy or re-implement the loop; run it.  The ZWFS reference
implementation is `zwfs_dm96/zwfs_run.m` -> `stage_loop_` (+ `P.loop` in
`zwfs_params.m`, `noisy_frames_`, `hold_photons_`, the loop figure in
`zwfs_run_figs.m`): lift it after the matrix calibration (deliverable 2)
is in `tg96_run`, because the loop reads the state through the measured
matrix on the working surface.

What `tg96_run` stage 'loop' has to provide to `dmg_loop`:
- `ins.measure(cmd)`: the DM at command `cmd` (nact x nact, mm) traced,
  the four fringe frames captured NOISELESS (your `frames`/four-step
  capture);
- `ins.noisy(F, nph, seed)`: photon noise at `nph` photons per MEASUREMENT
  (one DM shape measured once -- Dave 2026-09-11: say "measurement", not
  "state", which collides with state-vector controls terminology)
  split over the four frames (nph/4 each), from a `RandStream` seeded
  with `seed` (the S5 model, `zwfs_run` `noisy_frames_` verbatim in form);
- `ins.diff(F1, F0)`: the four-step differential map between two frame
  sets, mean-referenced over the mask (the reference ruling above);
- `ins.est(map)`: the matrix estimator (`est_matrix_`), on the matrix
  calibrated ON the working surface (`calib_surface 'base'`);
- `ins.lit`: the lit actuators.
Then the same run matrix as the ZWFS: readings = the four-step reading
(one row; add the two-step or single-frame variants only if the runner
has them), drifts 'walk' 2 pm per actuator per cycle and 'thermal' 5 pm
rms per cycle, noiseless steps 1 and 10 nm, photons per cycle {1e12,
1e13, 1e14, 1e15}, g 0.5, K 60, the noise-only floor at every level,
`seed 77` (the drift realization is drawn on the full actuator grid and
masked by lit, so both instruments see the SAME pattern where both are
lit), `hold_spec 3e-9` (3 pm).  Cost: K + 1 traced states per run; at
one state ~ 4 frames the record run is an hour-class background job
(`zwfs_batch.sh` pattern, memory-capped).

Report (the same tables as the ZWFS stage, so the comparison is row by
row): the step response (rho, tau, the noiseless floor), the hold error
vs photons per cycle for none / walk / thermal (ss, bias, sig_n, the
theory line), the ONE number = photons per cycle to hold 3 pm rms per
drift, and the spectrum bands of the held residual.  Lens rig first (the
record), OAP rig when its matrix exists.  Also say what the four-step
reading's 4-theta harmonic does in the loop (absolute vs differential:
the differential to the set point's frames should cancel a fixed
harmonic; measure, do not assume).
