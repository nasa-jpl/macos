# CALIB: a failed trial is a REJECTED trial; zero accepted steps is 'failed' or 'stalled' (PLAN_CONSOLIDATION 3.4, 2026-10-11)

**Change.**  `funcs_app` turns an evaluation that could not be traced (`lOpsFail`: a trace inside
it lost the chief ray) into yfit = 1e100, so `mrqmin` REJECTS the trial by its own rule (lambda
x10, parameters unchanged; nls.F untouched), counted in `n_failed` beside `n_rejected`.  After the
loop, ZERO accepted steps is never 'floor': with failed trials -> **failed** (catchable error,
'no trial step was accepted -- N of K trials lost the chief ray'); without -> stop reason 5
**stalled** (the start is a stationary point -- typically already optimal; counts as converged).
A failed INITIAL evaluation still aborts (338416c).  Rulings: CC, three rounds (the first
`converged = false` for stalled amended after the four movers below).

**Matrix (tCalibMatrix cells + TIP+TILT variants, 24, mmacos, model 128) vs the pre-change build:**
all 13 well-posed cells and all BEAM cells IDENTICAL.  Movers, both intended:
| cell | before | after |
|---|---|---|
| rigid x SPOT, TIP+TILT free | failed at iteration 1 ('a trial step lost the chief ray') | failed: 'no trial step was accepted -- 6 of 10 trials lost the chief ray' (10 trials, optics restored) |
| mixed x SPOT, TIP+TILT free | failed (gaussj singular at 10) | unchanged |
A TILT-only rigid solve (CC's candidate for a stalled fixture) is 'failed' too: 1 of 10 trials
loses the chief ray.  The stalled fixture is `Rx_CalibOptimal` (the UNspoiled paraboloid, h^4
free) under SPOT: 10 trials, all evaluated, none lowers the merit -> 'stalled'.

**CALIB classes (tAsphCalib, tBeamRows, tCalib, tCalibTargets, tCalibStop, tDesignTelescope,
tSpectrometerRx) vs the 338416c build:** iterations and merits IDENTICAL everywhere; the STOP
REASON moves on four cases whose start was already optimal (no trial lowered the merit):
- tBeamRows / centroid leg (merit 0 from the start): 'round-off floor' -> 'stalled'
- tDesignTelescope, one optimize case: 'round-off floor' -> 'stalled'
- tSpectrometerRx, the two T=1 calibs: 'OptMaxItrs reached' -> 'stalled'
`converged` stays true for all four (no re-pin).  The result structs carry the reason:
`macos.calib` / pymacos `calib()` `stop_reason`, `Telescope.optimize` `.stop` and `.n_iter`.
