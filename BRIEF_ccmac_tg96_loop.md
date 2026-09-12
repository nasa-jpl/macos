# BRIEF for CCMac: the closed-loop hold metric on the T-G IFO (deliverable 7), on dev-candidate

From CCL for Dave, 2026-09-12.  Dave: "Push yes. Merge yes. Brief CCMac!"
Written for a cold start.  `tg96-oap` is MERGED into MACOS_resources
`dev-candidate` (merge 0fd6786) and both repos' dev-candidates are pushed;
work on dev-candidate from here, no more branch merges.  Pull both repos
first.  Read, in this order: this file; `BRIEF_ccmac_tg96_oap3.md` (the
review + its 2026-09-12 addendum: the bare-Al / Jones-pupil question is
item B below); `BRIEF_loop_metric.md` (the spec); `zwfs_dm96/README.md`
S11 and V1 (the ZWFS record you are comparing against);
`dm_gauge_lib/dmg_loop.m` header and `tests/tDmgLoop.m`.

## What is now on origin/dev-candidate that you have not seen

- `dm_gauge_lib/dmg_loop.m` -- the ONE closed-loop hold code for both
  gauges (proportional loop through an instrument given as four handles
  + a lit mask; drift 'walk' / 'thermal' / 'step'; set-point reference;
  theory lines; divergence guard `rmax`).  `tests/tDmgLoop.m`: 9 gates on
  a synthetic instrument, in SUITE_FAST.  Do not copy or re-implement.
- `zwfs_run` stage 'loop' (`stage_loop_`, `noisy_frames_`, `hold_photons_`)
  + `P.loop` in `zwfs_params.m` + the loop figure in `zwfs_run_figs.m`:
  the reference implementation to mirror, line for line where the
  instrument allows.
- The ZWFS S11 record (`runs/loop193`, `loop385`) and the V1 vector
  reading (`runs/vloop193`): the rows your IFO row goes beside.
- Terminology (Dave): photon budgets are "per MEASUREMENT" (one DM shape
  measured once, every photon the camera detects over the pupil, summed
  over the frames that reading needs); never "per state" (collides with
  state-vector controls), not "per frame".  Decks in American English.

## A. Deliverable 7: `tg96_run` stage 'loop'

Instrument handles for `dmg_loop`, exactly as `zwfs_run` builds them:

- `ins.measure(cmd)`: the DM at command `cmd` (nact x nact, mm) traced,
  the four fringe frames captured NOISELESS (your four-step capture).
- `ins.noisy(F, nph, seed)`: shot noise at `nph` photons per MEASUREMENT
  split over the four frames (nph/4 each) from a `RandStream` seeded
  with `seed` (lift `noisy_frames_`: `I .* (1 + randn ./ sqrt(max(I /
  sum(I) * n, 1)))`).
- `ins.diff(F1, F0)`: the four-step differential map between two frame
  sets, mean-referenced over the mask; the phase DIFFERENCE wrapped
  (the ZWFS lesson from V1: absolute maps wrap at +-pi individually).
- `ins.est(map)`: the matrix estimator (`est_matrix_tg`) with the matrix
  calibrated ON the working surface (`calib_surface 'base'`, 30 nm rms,
  seed_base) -- the loop reads through the operating-point matrix.
- `ins.lit`: the lit actuators.

Run matrix = the ZWFS's (`P.loop` defaults): set point = the 30 nm
working surface; g 0.5; K 60; photons per measurement {1e12, 1e13,
1e14, 1e15}; drifts 'walk' 2 pm per actuator per cycle and 'thermal'
5 pm rms per cycle (defocus + astigmatism); noiseless steps 1 and 10 nm;
the noise-only floor at every level; `seed` 77 (the drift is drawn on
the full actuator grid and masked by lit, so both instruments see the
same pattern); reference frames 'noiseless'; `hold_spec` 3e-9 (3 pm);
`rmax` 1e-3.  Report the same four tables (`stage_loop_`): step response
(rho, tau, r(K/2), r(K)); hold error vs photons for none / walk /
thermal (ss, bias, sig_n, theory); the ONE number (photons per cycle to
hold 3 pm per drift); the spectrum bands of the held residual.  The
figure from the same code path.  Lens rig first (the record), the OAP
rig second (coated, once item B says which coating).  Cost: K + 1 states
per run, 14 runs per reading -- an hour-class background job
(`tg96_batch.sh`, memory-capped).

**The comparison row** (numbers, side by side, same seeds):

| reading | 3 pm held, noise only | 3 pm held, 2 pm walk | thermal floor | noiseless step at cycle 60 |
|---|---|---|---|---|
| ZWFS linear L (1 frame) | 2.1e12 | 7.3e12 | 27.6 pm | 1.2 pm, still falling |
| ZWFS exact one-frame I+ | diverges | diverges | diverges | 99 nm |
| ZWFS stepped S (4 frames) | 2.6e12 | 7.5e12 | 10.0 pm | 0.000 pm |
| ZWFS polarized pair V (2 frames) | 1.5e12 | 5.3e12 | 9.9 pm | 0.000 pm |
| T-G IFO four-step, lens | ? | ? | ? | ? |
| T-G IFO four-step, OAP | ? | ? | ? | ? |

Expected from S5 / S11 physics: the IFO's single-shot noise is ~8e14
photons per measurement for 1 pm (vs the ZWFS ~5e13), so its 3 pm
crossing should sit ~15x to the right; what the loop will tell us is
whether it has a FIXED error (the 4-theta harmonic in absolute mode;
the differential to the set point's frames should cancel a fixed one --
measure) and whether its low-order gain holds a ramp at rate/g.  Say
what you measure, especially where it departs from that.

## B. Before the OAP loop row: the bare-Al / Jones-pupil runs (oap3 addendum)

The D5 "0.95" is accepted as a measurement, not yet as a design number.
Two one-command runs: (1) bare aluminium on both OAPs (`coat_set`,
n = 1.373, k = 7.62 at 632.8 nm), same rows as D5; (2) the Jones pupil
of the test arm (`macos.jones_pupil` + `pol_maps`) for ideal / bare Al /
protected Al, with the fringe visibility at the band centre and at the
edge as numbers.  If bare Al carries the band, the coating's retardance
(at both OAP AOIs, with its chromatic slope) is a design parameter and
the cheaper fix is re-solving the test-arm waveplate azimuth for the
folded arm; if it fills in, the ideal-reflector idiom was the
idealization and D3 is retired.  One paragraph in REPORT_oap.md with the
numbers, then the OAP loop row with whichever coating the answer says.

## Deliverables, in order

1. `tg96_run` stage 'loop' + `P.loop` knobs; lens rig; the four tables +
   figure; the comparison row filled in.  Commit the run artifacts
   (pruned results + report + figure), as you did for 662e76a.
2. Item B (two runs + the paragraph).
3. The OAP loop row.
4. README + REPORT_oap.md updated (numbers first; departures flagged);
   `zwfs_dm96/README.md` gets NO edits from you -- put the IFO rows in
   `tg_psi_dm96_oap/README.md` and I fold the comparison into the
   campaign record and the deck.
5. Commits on dev-candidate; push when the lens loop row exists and
   tBench is green; say so in your report.

## Constraints (standing)

No engine work; 'lens' byte-identical (tBench 9/9); one MODEL-1024
MATLAB at a time; `exit(0)` only in the batch wrapper; the loop code is
`dmg_loop`, never a copy; every number in the report.
