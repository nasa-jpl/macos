# BRIEF — the per-gauge resolution sweep (TO), one `format`-class run

CC for TO, 2026-10-13, on Dave's word (deck_cophi review: "ask TO to run the resolution
sweeps").  Self-contained: read this, `macos/BRIEF_cophi_campaign.md` §7 (rulings; ruling 6 =
full actuator resolution, nothing coarser than 1:1 is a gauge), the template
`MACOS_resources/mmacos/templates/40_benches/cophi_dm96/` (README, `cophi_params.m`,
`cophi_run.m` — the `noise` and `format` stages and the `dmg_*` noise hook that bins cameras
onto the model grid), and the deck slide it feeds: `macos/demo_session/deck_cophi.md`, "How
much resolution? IFO, vector ZWFS, HET" (the table whose 1 / 1.5 / 2 / 4 rows are today
INFERRED from the blur demo and the well-bound photon rule; every record is at 3.3–4 px per
actuator).  Both repos `dev-candidate`; state SHAs when you commit.  One model-1024 MATLAB on
the box (`ps -C MATLAB` first); do not relink the shared mex.

## The question

For each gauge, what does pixel count per actuator buy and cost — where does the actuator
fit collapse, and above that is the hold / speed really indifferent?  Dave's words: "price
resolution for the IFO, vZWFS and HET, trading hold performance, speed, complexity".

## The sweep

Sampling: **1, 1.5, 2 and 4 pixels per actuator across** (plus the record's own 3.3 as the
anchor).  Gauges: IFO (Twyman-Green, the record rig), stepped ZWFS, vector ZWFS, pinhole /
P-SRI, HET with Feng's camera, HET lock-in reading (`opt.iq`).  The mirror rig of record
(`bench.optics 'oap'`, plain splitter for HET), the 30 nm working surface, seed and matrix
conventions of the record (the composite-deck rows: `redo96_oap*`, `camS_oap`, `camI_oap`,
`cophi_dm96/runs/*_oap_plain`, `runs/lockin_*`).

Mechanism — your call, but say which and why in the README: binning the model's 1024
detector grid (k = 1 → 3.3, k = 2 → 1.64, k = 3 → 1.1 px per actuator; non-integer targets
by resampling the pupil image, with the same resample applied to the matrix calibration so
registration is treated as the bench treats it) for ≤ 3.3; the tail's magnification (the
field lens / `385 px/pupil` label in the layout) or a finer source grid for 4.  Whatever the
mechanism, **the matrix is measured THROUGH the resampled reading** (blur and registration
absorbed by calibration, as on the bench) and the registration error row is ± 0.1 pitch.

Per (gauge, sampling) record: (1) the battery's three columns — 10 nm single actuator, 1 nm
on the grid, dense random 10 nm: gain / floor / error, one-step AND converged; (2) the
noise stage's pm per reading and seconds per pm on the gauge's record camera (the
polarization camera for IFO / vector ZWFS, Photonfocus for HET, C4 for the lock-in) — binning
restores read noise per sample, say how; (3) the loop's noise-only hold and the thermal-ramp
row at the record's cycle time; (4) the registration row (± 0.1 pitch: does the fit read the
wrong actuator — the blur demo's collapse signature); (5) the data rate (pixels × 16 bit ×
frame rate) and the per-pixel solve cost for the ZWFS (it solves pixel by pixel).

## Gates

- A `tResolutionSweep` (SUITE_FAST-sized subset: one gauge, two samplings) that pins the
  mechanism: the anchor sampling reproduces the record's battery row to the digit; a
  below-collapse sampling fails the single-actuator row by name.
- The record's numbers at 3.3 unchanged (the sweep adds, it does not move).
- Non-vacuity: the collapse must APPEAR in the sweep (1 px per actuator: the single-actuator
  gain or floor visibly broken) — a sweep where every row reads the same is a binning that
  never reached the reading.

## Output

`runs/ressweep_<gauge>/` with one table `ressweep.txt` (gauge × sampling × the five
quantities) and the composite figure the deck slide will carry: hold (pm) and seconds per pm
vs px per actuator, one line per gauge, the collapse band shaded.  Report to CC: SHAs, the
table, what surprised you (the inferred rows of the deck table that the records overturn).
Expected: one day of runs; if a gauge needs its own mechanism (the ZWFS's core sets its
exposure), say so and do the others first.
