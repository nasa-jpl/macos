# BRIEF: the closed-loop hold metric -- one number to compare the DM gauges on-orbit

Dave 2026-09-11: "On-orbit, the DM surface needs to remain constant to
<< 10 pm, with frequent remeasurement and closed-loop DM actuator servo
control.  How can performance in this mode be made into a metric that
will allow direct comparison between these methods?" -- then "let's
build the loop metric for both the T-G IFO and the ZWFS."  Written for
a cleared or compacted session; the ZWFS half is CCL's (this box), the
IFO half rides CCMac's `tg96_run` (branch `tg96-oap`,
`BRIEF_ccmac_tg96_oap2.md`) once its matrix calibration lands.  The
loop code itself is SHARED (one file in `dm_gauge_lib`) so both
instruments run the identical loop.

## What closed loop changes (the reasoning that sets the metric)

- A sensor GAIN error only sets the convergence rate: a loop of gain g
  corrects a fraction g*G of the residual each cycle and keeps going.
  It is not a floor.
- Every noiseless systematic measured so far is MULTIPLICATIVE (floor
  proportional to the change, sensitivity stage S2: no additive floor
  over five decades).  In closed loop the change per cycle is
  picometres, so those floors shrink with the residual.
- What survives is what is injected or biased EVERY cycle: photon noise
  (measured: ZWFS one-frame ~9e13 photons/state for 1 pm, stepped
  1.0e14, IFO 8e14 -- S8/S5), any state-dependent bias (the ZWFS
  one-frame readings' sign fold; a stepped reading's reference-intensity
  drift; the IFO's 4-theta harmonic in absolute mode), and calibration
  ageing -- which a HELD surface suppresses (the matrix measured on the
  working surface ages 7 % per 20 nm rms of surface change, S10).
- A fixed sensor bias b makes the loop converge to -b, not to zero: the
  held surface sits at the negative of the sensor's bias on the
  zero-change state.  In the noiseless model that bias is zero
  (identical states read identically); the loop test is what shows
  whether a real bias appears once noise, the fold, or the reference
  drift enter.

## The metric

**Closed-loop hold error**: the steady-state rms surface error over the
lit actuators, in pm, when the DM is driven by a specified drift and
held by a loop of gain g at cadence T with N photons per cycle, read by
the instrument and fitted through its measured response matrix.
Report it as a curve against photons per cycle, plus the time constant
(cycles to 1/e of an initial error) and the spatial spectrum of the
residual.  The ONE-NUMBER comparison: the photons per cycle each method
needs to hold the surface at a stated level (3 pm rms) against the
stated drift.  Secondary: the residual with the noise removed (the
bias floor), and the largest per-cycle disturbance each method tracks
without aliasing (dynamic range in this mode).

## The loop (shared code: `dm_gauge_lib/dmg_loop.m`)

Inputs: a DRIFT generator (actuator-command increments per cycle), a
MEASURE function (DM command -> the instrument's map, re-traced; the
frames captured noiseless), a NOISE injector (photons per state split
over the reading's frames, the S5 model), an ESTIMATOR (map -> actuator
changes; the measured matrix), loop gain g, cycles K, seeds.  Per
cycle: truth += drift; frames = measure(truth) with noise; the
reference is the PREVIOUS cycle's frames (the differential protocol,
as every row so far); a_hat = estimator(frames - ref); command -= g *
a_hat.  Record the residual (truth - target) rms over lit each cycle.
Outputs: residual history, steady-state rms (mean over the last K/2
cycles), time constant, the residual's spatial spectrum (the (p,0)
projection), and the noise-only floor (the same loop with drift = 0).

Drift models (parameters): 'walk' -- per-cycle Gaussian increments of
sigma_d per actuator (a random walk); 'thermal' -- a fixed low-order
shape (defocus + astigmatism) growing linearly at r pm per cycle;
'step' -- one disturbance of amplitude D at cycle 1 (measures the time
constant and the dynamic range).  Defaults: walk 2 pm/cycle, thermal
5 pm/cycle, step 1 nm; g = 0.5; K = 60; photons per cycle {1e12, 1e13,
1e14, 1e15}; working surface 30 nm rms (the operating point; the
matrix calibrated ON it, S10).

Gates: (G1) noiseless, drift 0, an initial error of 1 nm: the residual
falls geometrically at (1 - g*G) per cycle to < 1 pm within K/2 cycles
for every reading (the estimator is consistent); (G2) noise only, drift
0: the steady-state rms equals sigma_n * sqrt(g/(2-g)) with sigma_n
from the S5-style single-shot noise -- the loop propagates noise as
theory says; (G3) non-vacuity: a reading known to be biased on the
working surface (the ZWFS one-frame exact reading without its refined
sign map) converges to a NON-zero surface, and the report says by how
much.

## Instrument halves

- ZWFS (CCL): `zwfs_run` stage 'loop' (`P.loop.*` knobs), readings L /
  I+ / S with `calib_mode` 'matrix' on the working surface (S10:
  stepped 0.99 / 5 pm, linear 1.05 / 23 pm there); frames per state
  1 / 1 / 4.  Cost: K states re-traced per point (~0.7 s per frame at
  193 rays) -- 60 cycles x 4 photon levels x 3 readings at 193 rays is
  minutes, not hours; the 385-ray confirmation once.
- IFO (CCMac): `tg96_run` stage 'loop' with the same `dmg_loop`, the
  four-step reading through its measured matrix; 4 fringe frames per
  state; the same drift models, gains, photon levels and seeds.  The
  lens rig first (the record); the OAP rig once its rows exist.
- The comparison table: rows = instrument x reading; columns = photons
  per cycle to hold 3 pm rms (walk, thermal), time constant, bias floor,
  largest tracked step.  Both from the identical loop code, the
  identical drift realizations (seeded), the identical scoring.

## Deliverables

1. `dm_gauge_lib/dmg_loop.m` + a unit test on a synthetic linear
   instrument (a known J, Gaussian noise) that pins G1 and G2 to
   theory (the loop code must be right before either instrument runs
   through it).
2. `zwfs_run` stage 'loop' + `zwfs_params` `P.loop`; the ZWFS table at
   193 rays; the figure (residual vs cycle per photon level; the
   photons-to-hold curve).
3. README S11 bullet + `BRIEF_zwfs_campaign.md` section + deck slide
   (the one-number comparison, once the IFO half exists).
4. Hand the IFO half to CCMac as an addendum to
   `BRIEF_ccmac_tg96_oap2.md` (the loop stage on `tg96_run`, same
   `dmg_loop`, same seeds).

## Constraints

No engine work.  The runner rule (Dave, standing): every stage runs
THROUGH `zwfs_run` with its knobs in `zwfs_params`; every number in the
report; the loop code shared, never duplicated.  Sequential MODEL-1024
MATLABs.  Commit local; push on Dave's review.
