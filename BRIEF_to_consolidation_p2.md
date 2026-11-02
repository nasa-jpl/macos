# BRIEF — PLAN_CONSOLIDATION phase 2, engine lane (TO)

CC for TO, 2026-10-10, on Dave's word.  Phase 1 is EXITED (Dave, today) — every item
of `BRIEF_to_consolidation_p1.md` and the two rulings that followed it (OPD reference,
IRIS) are pushed and gated; the record is `macos/PLAN_CONSOLIDATION.md` §1.1 / §1.1a
and the §5 exit line.  Phase 2 per the plan's §3: "Commands + optimizer".  Read the
plan first — §0 (the ten rules) is still the definition of done, **§1.1b is your
first list** (the items phase 1 carried, Dave: "keep these undone items in the plan"),
§1.3 is the second.  Then `macos_f90/CLAUDE.md` (every item below has its entry).
Work in `~/dev/macos` and `~/dev/MACOS_resources`, `dev-candidate`, tips 5362009 /
400f825 on origin (CC's bookkeeping commit 257cf70 is local; it changes only
PLAN files).

Dave's rulings that frame this phase: the feature freeze is **not required but
preferred until the dev-candidate PR is approved** — so nothing new unless an item
below needs it; the `Opt*`/trace-state SAVE question is deferred (do not touch the
writer for it); COPHI starts later, not now.

## How we work (unchanged from phase 1, with the lessons added)

ONE item in flight.  Its gate FIRST, shown red on the pre-fix binary (keep the
`build_release_prefix*/` convention; CC gates from the clean worktree
`~/dev/macos_cli_base` and never from the shared tree).  Both compilers rebuilt.  The
mex relinked only when no other MATLAB runs (say when).  Fast suite after.  Corpus
A/B before the commit when the item touches surfsub / elemsub / tracesub / propsub /
msmacosio / iosub / macos_cmd_loop (`cli_tests/run_cli_tests.sh load` is the per-
commit form; the full corpus for a root-pick or reference change); every moved deck
named in a committed record (`cli_tests/records/ab_*_moved_<sha>.txt`).  Commit
pairs local, SHA + branch stated; CC gates each commit; Dave says push.

Lessons phase 1 paid for, now rules: **measure before you guard** — grep the value's
later assignments first (a parse-time `nGridpts` refusal broke 195 tests because both
load paths already clamp); **a guard warns and counts**, a new refusal needs Dave
(precedent: a crash or a corrupted trace → a clean NAMED refusal is approved);
**user options at allocation**, never in `*_init_vars`; **a pymacos red leg needs a
pre-fix `.so`** — build it in a side tree (the gap on `test_opd_ref`; the §1.4a
generator will make this routine, until then do it by hand and say so).

## A. The carried items (plan §1.1b), in this order

1. **3c — `tst_save_keys.in` traces to NaN at the GRIN leg** (its twin
   `tst_block_comment.in` is the same deck, one item).  Diagnose before anything:
   is `GradInd 6` with the fixture's coefficients a VALID gradient-index element, or
   a fixture that asks the integrator for something it cannot do?  If the deck is
   wrong, fix the FIXTURE (it exists to exercise SAVE keys, not GRIN physics) and
   say what was wrong with it.  If the deck is right, the integrator's NaN is the
   defect: find the step that produces it (the pty CLI + one `WRITE` per step, the
   pattern that found the Return parity).  Either way the OUTCOME must be: the
   corpus carries no NaN trace, and a GRIN leg that cannot integrate says so ONCE
   and marks its rays lost (`SetRayFail`, a new `RayStat_*` only if none fits) —
   warn and count, no refusal.  Gate: the deck's `opd nElt` finite on both
   compilers (CLI load gate record); the pre-fix binary prints NaN.  Note for the
   FwdRoot counter: this deck was the one `NaN .NE. NaN` case — after the fix,
   check the `ABS(new-old) > 0` test still has a reason to exist or say it does not.
2. **3b — the Refractor branch's IndRef stamp.**  Item 3 removed the Return branch's
   stamp of the running medium onto the Rx; the Refractor branch has the same
   statement, kept because it feeds the NEXT ray's `nb` for GRIN / variable-index
   elements.  Audit: does it ever change a value the Rx declared (a SAVE after a
   trace differing from a SAVE before, on any corpus deck with a Refractor)?  The
   corpus A/B already answers the surface question — `tst_save_keys` was the only
   mover and it was the Return.  So: name the decks it CAN change (variable-index,
   GRIN, LensArray-in-glass) and measure them; if it mutates the Rx, move the state
   to a per-ray running variable like item 3 did; if it does not, close the item
   with the measurement and a comment at the statement saying why it is safe.
3. **Even-grid point-source centring (6d13d07) has no gate.**  The fix completed PR
   #70 at the `PtSource` site — a grid of `nGridPts == mdttl` wrote `RayID(.,mdttl+1)`
   — and it CHANGES the trace of even-grid point-source diffraction (the nine
   `Rx_Coro*` decks at `nGridpts= 511` clamp to the model and are the class).  The
   gate is the PROPER compare: a point-source far-field or near-field phase on an
   EVEN grid (`pymacos/tests/proper_compare/`, phase 1's `Rx_Cass_FarField` form
   with the grid forced even), post-fix agreeing at the odd-grid level
   (1e-11 peak-normalised far field; the phase-2 near-field 2.4e-14 rms), the
   pre-fix `.so` off by the half-pixel.  If the pre-fix `.so` does NOT differ, the
   centring did not reach the diffraction seed and the item is "measured, no gate
   needed" — write that down with the numbers.
4. **`Grating_example_001.in` — eight knife-edge rays.**  7303 rays pass on ifx,
   7295 on gfortran (RMS OPD 177.168 vs 177.155): rays on the concave grating's
   aperture edge whose pass/fail depends on the compiler's arithmetic.  Find the
   edge (which rays, which aperture test, how far inside/outside in ulp); then
   decide with CC: move the deck's aperture by a pixel's worth so no ray sits on it
   (the fixture is a grating-GETTER test, not an aperture test), or accept it as
   the one known compiler-sensitive deck.  Today `compare_compilers.py` tolerates
   it by name — the outcome must let that exception go.  No engine change expected.
5. **The `SegDemo.in` load TIMEOUT** — seen in every corpus run on 10-09 and never
   since, on every binary including pre-fix ones.  Not an item to chase; a WATCH:
   if any CLI run of yours sees a load take more than a few seconds on that deck,
   stop and capture it (`gdb -p` on the CLI, the stack is the finding).  The
   `pupilsim_redo_lens_*` timeouts are CC's (model-aware driver).

Item 6 — `segment_rx` writing `NSCount= 1` for single-hit segmented groups — is
resources-side and small; take it after A.2 when you have the mex relinked anyway
(gate: an emitted segmented deck round-trips through SAVE with the line present
and traces identically with and without it; the fixture `tst_zrngr_roundtrip.in`
shows the form).

## B. The optimizer (plan §1.3 → PLAN.md §3), after A

The urgent two first, both in `design_optim.F` / `smacos_compute.inc`:

1. **`WFE_ZMODE_TARGET` and `OPL_TARGET`: implement or remove** (PLAN §3.1).  Since
   the CALIB wrappers landed, `m.calib_set_target('WFE_ZMODE', …)` is a documented
   api call that today no-ops or fails silently.  Measure first: drive each target
   from mmacos on a deck where the answer is known (a Zernike-mode WFE target on the
   `Rx_AsphCalib` paraboloid; an OPL target on a two-mirror deck) and record what
   the engine does.  Then ONE of: wire the `objfun_nom` branches so the target does
   what its name says, with a gate that converges a spoiled deck (the `tAsphCalib`
   form); or remove both from the enum, the parser and the wrappers so the api
   REFUSES them by name (a refusal of a non-functional mode is the approved
   precedent class — say so in the commit).  Recommend implement for WFE_ZMODE (it
   is what a segmented-pupil solve wants) and measure OPL before deciding.  Either
   way the pymacos docstring and the mmacos help text state the resolution.
2. **`SPOT_TARGET` reports 0 iterations** (PLAN §3.1) even after the `dopt_init_vars`
   fix.  Find the early-exit gate (it is a control-flow bug, not a physics one: a
   convergence test that fires at iteration 0 on a SPOT merit, or the stride bug's
   sibling); gate = `tCalib` gains a SPOT × rigid-body case that takes ≥ 1
   iteration and converges (the pre-fix mex returns at 0).

Then §3.2 (the keyword fixes — `OptTol`, the `OptMxItrs` prefix, `OptAsph` naming,
the stray `MBFile6: Unidentified string` under a non-matching target), §3.3 (the
convergence check every iteration — today it never runs for `tot_itrs ≤ 4`; one
`nls_optim_dvr` call site; `quit` after CALIB; the `dopt_init_vars` audit), each a
small commit pair with its own gate.  §3.4, the target × DOF matrix, is written
PROGRAMMATICALLY (`tCalib` / `test_calib.py`, the 2026-06-02 re-scope): twelve cases,
one fixture deck per DOF type, every case converging a spoiled deck; CC mirrors
each as a `cli_tests/cmd/CALIB_*.jou` so the CLI surface is gated too (rule 2).  The
asphere+Zernike surface option and EPFIX (PLAN §3.1, Dave 2026-10-05) are FEATURES —
they wait under the preferred freeze unless Dave says otherwise.

## C. After B (or when Dave re-orders)

The §1.2 deferred items in the order they bite — `NSRefractor`'s null grid frame
(the one remaining piston-only GridData path), the `ExtFringe` converter or its
error, `RayFailElt = nElt+1` for obscuration — then the bindings program §1.4a/b
(the one generator for mex + f2py + veneer; the 21 unbound pymacos wrappers; the
cross-binding `tDefaults` / `test_defaults.py`; the cmdref binding column).  Not
before B: the freeze is preferred and the optimizer is the user-facing hole.

## Reporting

One item in flight; a line in `NOTE_to_ccl_consolidation_p2.md` (or this channel)
when it lands — SHA pair, gate counts red/green on both compilers, the A/B's moved
list — and CC gates the commit from the clean worktree before it goes on the push
list.  What you find that is NOT the item (phase 1 found a defect behind nearly
every item) goes into the plan's §1.1b or §1.2 as a dated line, not into the fix.
