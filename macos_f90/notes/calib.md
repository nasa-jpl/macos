# MACOS engine notes — CALIB / design_optim: wrappers, targets, the LM stopping rule, beam rows, asphere DOFs and step

> **Read this file before touching:** `design_optim.F`, `nls.F`, `macos_ops.F`, `calib_*` api, `Opt*` keywords.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## CALIB wrappers in macos_api_mod (Phase 1a/1b/1c)
**Documented 2026-10-16:** the keyword family, defaults, per-load / SAVE status, row and unit
conventions and the stopping rule are Manual §4 Table 13 ("Design Optimization (CALIB)
Keywords"); the command and its dialog are cmdref `CALib` (+ SETC / AVAR / AFOV …), the
bindings cmdref Part II (`calib_*`).  Engine facts kept here: CALIB *is* `nls_optim_dvr`
(LM + SLSQP) in `design_optim.F`, the design layer's `native` engine (PLAN_DESIGN_LAYER §8);
8-DOF `VarElt` mask with **DOF 7 = radius, DOF 8 = conic**; `calib_run` + `calib_buffer_dims`
(Phase 1a), the setters `calib_set_var_elt` / `calib_clear_var_elts` / `calib_set_iter` /
`calib_set_tol` / `calib_set_target` (1b); `calib_clear_var_elts` carries `allocated()` guards
(1c: a no-op, not a SIGSEGV, on a fresh `macos_init_all` state); programmatic FIELD setters
(`calib_add_fov` / `calib_set_wavelens`) remain the deferred Phase-1d gap -- fields come from
the `.in` keywords.

## CALIB targets: WFE_ZMODE works, OPL removed (2026-10-10, PLAN 3.1)
**WFE_ZMODE (`OptTarget= ZWF`, api target 2) STOPPED THE HOST** on every deck
that does not declare `ZCOZernType=`: `macos_mod_init_vars` reset the Zernike-fit
basis `zcoType` to 0 on every load (after the LOAD handler had set
`ZernCircular` -- the rule-6 trap again), so `zern_fit` took its non-circular
branch with type 0, read `ZernAnnuRatio(ie_m=0)` out of bounds, the fit's SVD
failed and `nerror` STOPped ('svd_decomp: no convergence'; the ifx .so got
garbage instead and `lmlsq` failed).  Now the per-load reset is TO
`ZernCircular`.  Measured on `Rx_AsphCalib.in` (h^4 spoiled to 5e-3): modes
[4 11 12] -> 3.5e-13, [11] -> 1.4e-12, [12] leaves it (4.4e-3) -- the modes are
NOLL, 1-based (11 = primary spherical).  **OPL (target 5) had no objective**
anywhere in design_optim: removed from the enum; refused by name at the api and at the Rx
(Manual Table 13).  Gates
`tCalibTargets` (mmacos, the ZMODE leg in a SUBPROCESS) + `test_calib_targets.py`;
red on the pre-fix mex / .so.  `zern_fit`'s pupil = `opdValid` since 2026-10-12
(see the OPD-reference section's last paragraph).

## CALIB's LM stopping rule: the derivation (2026-10-11, PLAN 3.3)
*(The rule as the user sees it -- OptTol, the five stop reasons, what is printed -- is Manual
§4 Table 13 "Stopping"; this section is the engine record and the derivation.)*
**The old loop** (`nls_optim_dvr`, design_optim.F) checked convergence only at itr =
int(.21T), int(.5T), int(.71T) and only when T > 4, judged the RMS WFE whatever the target,
and exited on `itr_done .and. itr > 0.7T` -- with `itr_done` never set before the loop.
`nls_optim_dvr` carries a blanket `Save`, so `itr_done` is a STATIC that starts .false. and
then PERSISTS between calls: a T<=4 CALIB (no check ever runs) exited on whatever the
previous CALIB left.  Measured (CLI, both compilers, pre-fix and a -finit-logical=true/false
pair, identical): `tst_calib_spot4.in` 4 iterations fresh, 3 after a converged
`tst_calib_spot20.in` in the same process.  The 0.7T floor also meant no solve ever stopped
before 70% of its cap.
**The rule now (LM / NonLin only; the SVD path keeps the old schedule, `itr_done`
initialised):** accepted = lambda fell (nls.F `mrqmin_private`).  `rms = sqrt(chisq/rows)` in
the target's units.  Stop when
- `converged` (calib_stop 1): two consecutive ACCEPTED steps each with `|d rms| <= OptTol`
  (ABSOLUTE, the old WFE check's own definition; api `calib_set_tol`, default 1e-12) AND
  `chisq_prev - chisq <= 1e-6 * chisq_prev` (`ftol_rel`, a fixed constant -- MINPACK's
  ftol idea; OptTol stays the only user knob);
- `floor` (2): TEN consecutive REJECTED trials since the last accepted step (lambda up 1e10);
- `cap` (3): T; `failed` (4) on the LM failure path.  Never before itr 2.
Reported: CLI `Final LM merit (chisq) =` + `Stopped: ...` after `Optimization iterations =`;
api `calib_stats_get(n_iter, merit, stop_reason)` (stop_set's stop aim, which also calls
`nls_optim_dvr`, saves and restores them); mmacos `macos.calib` `.n_iter/.merit/.stop_reason`
+ `macos.calib_stats`; pymacos the same keys + `calib_stats()`.
**Why each number, MEASURED** (full-precision full-T chisq/lambda traces of every CALIB in
tAsphCalib / tBeamRows / tCalib / tDesignTelescope, scored against the old loop):
- **Do NOT shorten the floor.**  K=3 stopped `tDesignTelescope/base_sphere_zernike_correction`
  at itr 9 with chisq 481x the old loop's 120-iteration result (WFE 1.16e-4 vs 5.28e-6 m; the
  test still passed -- its gate is `wfe_after < wfe_before/100`): a 3-rejection run at itr 7-9
  is followed by an accepted step lowering chisq 6.2%.  The longest MID-SOLVE run of
  rejections a real decrease later followed is SIX (`convex_conic_seed_zero_and_optimizes`,
  itr 9-14, then -1.4e-6 relative); K=10 is four over it.
- **At a converged design's residual** (chisq ~2.3e-16, WFE ~1.3e-9 m) LM recovers after runs
  of 11-14 rejections for 1e-7..6e-7 RELATIVE chisq (`tma_multifield` 14 rejects then
  -5.9e-9 / -7.3e-8; `elts_subset` 11 then -5.8e-7; `2d_field_cross` 12 then -4.4e-7) -- each
  ~4e-16 m of WFE against OptTol 1e-12 m, round-off refinements.  A strict-relative reading
  of the tolerance would need K >= 15, and the floor would then save 8%.  Who wants a
  stricter stop sets OptTol, not K.
- **The relative term is load-bearing.**  With the absolute test alone, five
  tDesignTelescope cases stopped 'converged' at itr 3 with chisq +4..39% (fourth_mirror
  1.6587e-16 -> 2.2984e-16): on a nm-class merit the first accepted steps of an LM still
  ramping lambda down from 1e-3 lower chisq ~1e-3 relative = ~1e-12 m of rms and read as
  flat.  The old 0.7T schedule hid this.  rtol 1e-6 (simulated, with the absolute test only
  stopping later): worst +2.7e-6 relative chisq (~1e-15 m), iterations 804 -> 594.
- A quadratically convergent solve (Rx_AsphCalib SPOT: 1.6e-25 -> 1.6e-33 -> 5.4e-43 ->
  1.7e-43, measured at field weight 0 -- x1e10 since 2026-10-17, below) never takes a flat
  step; it stops at the FLOOR (spot20: 14), never 'converged'.
**A failed evaluation is a named failure, not a STOP (2026-10-11, PLAN 3.4).**  The matrix's
rigid x SPOT cell with the unspoiled TILT free: a centroid-referenced spot size is first-order
NULL in that DOF, J'J near-singular, and LM's first trial rotated the mirror ~76 deg; the
beam-frame chief ray (`LocalCoord`) then missed the focal plane and `MACOS_OPS` STOPPED --
killing MATLAB (rc=139).  Both `MACOS_OPS` STOPs (the SPOT beam frame, 'load Rx first') now set
`macos_ops_mod%lOpsFail` and return; `nls_optim_dvr` resets it at entry and checks it after the
initial evaluation and after every `lmlsq` call -> the existing failure path (last accepted
state, `rtn_flg=1`, stop 'failed').  Ill-posed problems END BY NAME: BEAM x shape DOFs (a zero
Jacobian column -> `gaussj` singular at n=0), SPOT x a null DOF.  **Follow-on (same day): a
failed evaluation INSIDE an LM trial is a REJECTED trial** -- `funcs_app` returns yfit = 1e100 so
`mrqmin` rejects it by its own rule (no change to nls.F), counted in `calib_n_fail` beside
`calib_n_rej`.  Stop reasons now: 1 converged, 2 floor (progress THEN ten rejections), 3 cap,
4 failed (incl. ZERO accepted steps with failed trials: 'no trial step was accepted -- N of K
trials lost the chief ray'), 5 STALLED (zero accepted, none failed: the start is a stationary
point of the merit, typically already optimal -- rtn 0, `converged` TRUE; four existing suite
cases with an already-optimal start read 'stalled', iterations and merits unchanged).  A failed INITIAL
evaluation still aborts (nothing to retreat to).  Gates `tCalibMatrix` /
`test_calib_matrix.py` (the host-survival case runs in a subprocess).
Setup: the variable packing and the configuration checks that precede both
`nls_optim_dvr` entries (the LOAD handler's `OptNomSens= Y` nominal-sensitivity run and
the CALIB command) are ONE function, `calib_setup_vars` (design_optim_mod, 2026-10-11) --
change it there, not in macos_cmd_loop.inc.
Gates: `cli_tests/cmd/CALIB_history.jou` (the same T=4 count fresh and after a T=20 CALIB, a
regex BACKREFERENCE; red on prefix9 + the old-loop side builds: 4/15/3), `tCalibStop`
(mmacos, SUITE_FAST) + `test_calib_stop.py`; the old-vs-new table is
`cli_tests/records/calib_loop_old_vs_new.md`.

## CALIB beam rows ride on ANY target; per-field targets; centroid (2026-10-03)
*(Keywords, defaults and the rows/units convention: Manual Table 13.  Below: the mechanics.)*
`design_optim.F`: the beam rows (`OptBeamDir=` chief direction at an
element, `OptBeamPos=` position, `OptBeamSize=`, reference rays) were
sized and scored ONLY under `OptTarget= BEAM`; the BEAM-only value /
derivative / linear paths reset their row offset to 1 for every field (a
multi-field beam solve scored the LAST field alone); the ZMODE derivative
stepped by `opd_size` and used `obj_size-1` modes; and the beam rows were
written one row LATE (`off+1`) after a WFE / SPOT block.  Now ONE helper
`beam_rows_` writes the rows for every target in the sizing order
(dir 3, pos 3, size 1, refray 3 x nRefRay), each (wavelength, field)
block is `obj_size + mBeamOff` rows, `sig` divides the beam rows by
`sqrt(OptBeamWt)`.  New (OptBeamStr fields + `msmacosio.inc` keywords,
reset per load in `macos_cmd_loop.inc`): `OptBeamPosFov=` (one 3-vector
row per CALIB field, in order, <= 12 -- a dyson5 smile / keystone solve
needs one per field), `OptBeamWt=`, `OptBeamCentroid= Y` (the position is
the CENTROID of the rays that pass the train at `beamPosElt`, read from
the new `RayPosAtElt(3,mRay)` that CTRACE fills at that element, instead
of the chief ray; under a BEAM target the whole beam is traced with
`GBS` to that element).  `OptBeamDir=` is unitised at parse.  API:
`calib_set_beam(kind 1/2/3, iElt, target, on)`, `calib_set_beam_pos_fov`,
`calib_set_beam_wt(wt, centroid)`; mmacos `macos.calib_set_beam` /
`calib_set_beam_pos_fov` / `calib_set_beam_wt` + Session methods.
Legacy trap kept: the `OptBeamPos=` keyword resets `OptTgtElt` to its
element.  Units trap: WFE rows are in the target's units, SPOT and beam
rows in base units -- weight accordingly.  Gate: `tBeamRows` (mmacos,
SUITE_FAST, `Rx_BeamRows.in`): direction rows on a SPOT target (pre-fix:
zero rows, the mirror never tilts), per-field position targets both met
by one FP piston, centroid vs chief on a comatic field.

## Asphere DOFs from the design layer (2026-10-04, dyson5 round 4 step 4)
`elt_asph_get(OK, coef, n, iElt)` (api; mmacos `macos.get_elt_asph`, pymacos
`get_elt_asph`) reads `AsphCoef(1:n, iElt)` -- what CALIB's `OptAsph=` DOFs
leave on an element.  `Telescope.optimize` gained `'asph_elts'` /
`'asph_terms'`: the element is emitted `Surface= Aspheric` from a zero seed
(the freeform declare-to-perturb rule), `OptAsph= n t1..tn` is written
AFTER its `VarElt=` line (the parser keys it on `isVarElt`), coefficients
are read back into `spec.elt(k).asph`.  Trap it carries: CALIB's
zero-coefficient asphere step is sag-based ONLY with a CIRCULAR aperture
(`ApType==1`), else the legacy round-off step on a metre deck -- Telescope
decks declare `ApType= None`, so the hook declares a vertex-centred circle
ENCLOSING the element's footprint for the solve (for an off-axis section
that radius runs from the PARENT vertex to the far edge of the footprint,
the scale the term acts on; it clips nothing) and removes it afterwards.
Gate `tAsphHook` (SUITE_FAST): on the dyson5 TMA parent h^4+h^6 on M1/M3
take the 3-field WFE [529 484 753] -> [55 117 308] nm vs the conic-only
control; readback; clean emit; the circle gone.

## CALIB: asphere differential step + the LM failure path (2026-10-01, dyson5 beat 4c)
*(The user rule -- declare a circular aperture or the zero-coefficient step is round-off on a
metre deck -- is Manual Table 13 "Aspheric degrees of freedom".  Below: the mechanics.)*
`nls_optim_dvr` (design_optim.F) probed every aspheric coefficient with a
FIXED step `das = 1e-10` (x 1e-5 per higher order).  On a deck in metres an
h^4 coefficient of 1e-3 probed by 1e-10 moves the sag by 1e-15 m at a 70 mm
aperture -- round-off, a zero derivative column, `gaussj: singular matrix (2)`
at the first LM step (TO, sec. 3.5).  Now the step is `das_rel = 1e-3` of
the coefficient; a ZERO coefficient (an asphere started from a conic) takes
the step that moves the sag at the element's CIRCULAR aperture radius by
1e-7 |Kr| (`ApType==1`, `ApVec(1)` = radius); with no circular aperture the
legacy absolute step remains, now scaled per ORDER of the term
(`das*1e-5**(j-1)`; the old code scaled by the term's POSITION in the
list), and ONE line says so (that branch is still round-off on metre decks
-- declare the aperture).  The same failure also
killed the HOST: the loop's `if (.not. lmlsq_success)` branch had a bare
`stop` before its own `rtn_flg=1` (the CLI happened to survive only when
the failure came elsewhere; the mex segfaulted, TO sec. 3.6), and the
"last run" branch `return`ed without deallocating.  Both now `go to 105`:
the optics are put back at the last ACCEPTED parameter vector (the
pre-optimization state when the first step failed), `rtn_flg=1`, normal
cleanup -- which is what the CLI's CALIB handler (`Optimzation aborted!`)
and the bindings' `calib_run` already expect: it answers OK=FAIL, which the
mmacos wrapper raises as an ordinary, catchable MATLAB error (`mmacos:
calib_run failed`) after the engine has printed the reason -- so
`macos.calib`'s `converged=false` branch is unreachable from MATLAB; catch
the error.
**Follow-on (2026-10-03): the failure path leaked `lmlsq`'s SAVEd arrays.**
`nls.F`'s `mrqmin_private` allocates `atry/beta/da` on the `alamda<0`
initialisation call and frees them only on the `alamda=0` final call; a
CALIB that ends on the failure path (or `gaussj` failing, `return` before
the free) never makes that call, so the NEXT CALIB in the same process
died at the `allocate` with "Attempting to allocate already allocated
variable 'atry'" -- a Fortran runtime abort, host killed.  Invisible
before 10-01 because the failure `stop`ped the process anyway; caught by
the fast suite ONLY in order (`tBeamRows` right after `tAsphCalib`'s
failure test; both classes green standalone -- an order-dependent crash
is the signature of leaked SAVE state).  Now the initialisation frees any
leftover first.  Gate: the suite order itself; `scratchpad/two_class.m`
pattern (tAsphCalib -> tBeamRows in one process).
Gates: `tAsphCalib` (mmacos, SUITE_FAST) on `Rx_AsphCalib.in` -- a metre
paraboloid with a spoiled h^4 term that CALIB must drive back to zero, and
the no-aperture/zero-term variant that must FAIL with the flag, the optics
untouched and MATLAB alive.  Measured on TO's reproducer
(`dyson5_s4_r4n_seed.in` + `OptAsph= 2 1 2`): pre-fix singular at once;
post-fix 4 iterations in 9 s; the no-aperture variant returns to the prompt.

## A field whose rays are all lost is a FAILED trial (2026-10-30, IODC 2010 outer loop, TO)
`funcs_app` built the WFE rows of every (wavelength, field) block from `OPDm` as `MACOS_OPS`
returned it.  When a trial step loses EVERY ray of a field (a surface walked past the source
plane; a wild LM step), `OPD` zeroes `OPDMat` before its ray loop, finds no passing ray, prints
'All rays are obscured or lost' and returns the 99999d32 sentinel in the RMS ONLY -- the map
stays zero, so that field's rows read ZERO, the trial's chisq FELL, and lmlsq ACCEPTED it: the
optimizer was rewarded for losing a field.  Measured on the IODC 2010 symmetric pair at R 130
(three fields, PIST on four surfaces + the image plane): fields 2-3 lost mid-solve, the solve ran
to the cap with chisq 1.8e-9 on the axis alone, the AFTER print showed 9.9999E+36 for the two
fields; the zero Jacobian columns that follow are the `gaussj: singular matrix (2)` failures the
later runs hit.  Now, after each `smacos_compute_perturbed` (the value and every derivative
probe), `wfe_rms >= 1d30` on a WFE / WFE_ZMODE target sets `lOpsFail`, and the PLAN 3.4 tail
(yfit = 1d100, `nFailTrials`) rejects the trial -- the same rule as a lost chief ray, no new
refusal.  The deck-side half of the IODC item (the CALIB chief rows started 100 mm ahead of the
SEED's pupil and the solve moved surface 1 upstream by 31 mm, past the source plane) is the
emitter's to fix and is fixed there (`iodc2010_emit`: a CALIB deck backs its source off by the
lens's length, every field's chief starts on that plane).  Gates: tCalibStop 4, tCalibMatrix 18,
tCalibTargets 3, tAsphCalib 2, tBeamRows 3, tCalib 9 (39/39 on the fixed mex); the corpus A/B
record names the movers (expected none: the branch fires only on the sentinel).

## Field weights default to 1 (2026-10-17, PLAN 1.2 (a))
`fov_weight` was 0 after every load: the LOAD handler's `fov_weight=1.0` ran before `MBFile6`,
whose `reinitialise_variables -> dopt_init_vars` zeroed it (the `LOptIfFEX` trap), and
`design_optim` reads |w| < 1e-10 as 1e-10 -- every row of a deck without `OptFOVWt=` divided by
1e5, chisq 1e-10 x the weight-1 value, the `OptTol` test on an rms 1e-5 x the target-unit rms.
The default is now 1d0 in `dopt_init_vars` (the LOAD line kept: one default, both paths).  Gate
`cli_tests/cmd/CALIB_fov_weight.jou` (tst_calib_rigid.in 7.98e-43 -> 7.98e-33; OptFOVWt= 0.5
halves it).  **LM is invariant to a uniform row scale in exact arithmetic** (alpha, beta, chisq
all scale by 1/sig^2; gaussj's test is an exact zero pivot): measured over the CALIB classes,
pymacos calib tests and the CLI CALIB journals, every solve kept its path with chisq x1e10
EXCEPT the rank-1 ones -- `mixed x SPOT` (one spot-radius residual, two DOFs) is round-off
chaotic (power-of-2 weights reproduce bit-for-bit, w = 1, 1e-3, 1e-6 each end elsewhere; only
the old 1e-10 met T=10): it now runs T=35 (tilt 6.0e-10) in tCalibMatrix / test_calib_matrix /
CALIB_SPOT_mixed.jou; the Rx_AsphCalib SPOT floor stops moved 19 -> 17 and 15 -> 18 iterations.
**Open (not changed):** the post-solve covariance call (`alamda = 0`, undamped) on a rank-1
alpha can hit the exact-zero pivot and turn a finished solve into stop 4 'failed' (the optics
kept at the last accepted step): mixed x SPOT at T = 18..22.  Decks written by the design layer
always carried `OptFOVWt=` and are unaffected.


## SAVE writes the family (2026-10-18, Dave; PLAN_CONSOLIDATION section 4 item 1)
The 26 `Opt*` / `Var*` keys are written by SAVE when they differ from their load defaults
(`PrtCalibHeader` / `PrtCalibElt`, iosub.inc) and announced at LOAD (`CalibInForceNote`); the
mechanics and the reader-counting gotchas are in `notes/rx_io.md` (the 2026-10-18 entry).  Table 13's
"none is written by SAVE" is history.  Gate `cli_tests/cmd/CALIB_save_roundtrip.jou`.
