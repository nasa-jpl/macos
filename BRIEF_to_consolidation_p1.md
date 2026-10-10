# BRIEF — PLAN_CONSOLIDATION phase 1, engine lane (TO)

CC for TO, 2026-10-09, on Dave's word ("start on PLAN_CONSOLIDATION … take it slow
and careful").  Read `macos/PLAN_CONSOLIDATION.md` first — §0 (the ten rules) is
the definition of done for everything below, §1.1 is the list, §2 is the suite
CC is building in parallel.  Then `macos_f90/CLAUDE.md` (the engine cheatsheet;
every item below has an entry there or in PLAN.md §0).  Work in `~/dev/macos`
and `~/dev/MACOS_resources`, `dev-candidate`, tips 4a24ba3 / 465fe18 or later.

## The order, one item at a time, each its own commit pair (engine + gate)

Slow and careful means: ONE item in flight; its gate written FIRST and shown red
on the current binary; both compilers rebuilt (`makems.sh release` and `makems.sh
release gfortran`); the mex relinked only in a window when no other MATLAB runs
(CC's harness runs the CLI only — never MATLAB — so CC is never the blocker; say
when you relink); the fast suite after; the corpus A/B before the commit when the
item touches the files in §2.3 (surfsub / elemsub / tracesub / propsub /
msmacosio / iosub / macos_cmd_loop).  For the A/B until CC's `cli_tests` lands,
use the FwdRoot driver: `scratchpad/ab_corpus.py <binary> decks_all.txt out.csv`
(CC's scratchpad, 2026-10-04 session; CC will copy both into `macos/cli_tests/`
today) — pre-fix binary vs post-fix, `opd nElt` RMS / P-V / nPass / lost identical
on every deck the change is not meant to move, each moved deck named in the commit.

1. **Short multi-value lines: `ZernCoef=`, `MonCoef=`, `FFZernCoef=`, and the rest
   of the bare `READ(VALUE,*)` keywords** (PLAN §0 line 270; CLAUDE.md "Short
   multi-value lines").  `ReadRealsPad` exists (elt_mod) and is wired for
   `AsphCoef=` / `AnaCoef=`; extend it to every multi-value real keyword in
   `msmacosio.inc` (list them in the commit: ZernCoef, MonCoef, FFZernCoef,
   MonZernCoef, FFCoef, GridSrfdx?, ObsVec/ApVec groups — audit, do not guess).
   Gate: `tRxShortAsph`'s sibling `tRxShortCoef` (mmacos) — a deck with a short
   `ZernCoef=` block loads, pads with zeros, prints ONE line; the pre-fix mex
   kills MATLAB (the must-fail leg is the engine abort itself: run it in a
   `matlab -batch` subprocess and assert the exit code).  CLI leg: the same deck
   in `ZGD_test_files/tst_short_coef.in` loads on both compilers.  Guards warn,
   not error.
2. **`macos.trace(26)` then `trace(27)` on the jwst zoom deck returns a wrong
   first OPD** (3.72982e-05 vs 6.85038e-06; PLAN §0 line 51).  Check the restarted-
   trace `PrevNonSeg` fix of 2026-10-03 first: reproduce on the current mex; if it
   is gone, the gate is the reproduction (`tTraceRestart` gains the zoom pair) and
   the item closes with the measurement, not a code change.  If it persists, it is
   the incremental-trace / cached-OPD path: bisect with the pty CLI (`onecall_cli.py`
   pattern: `opd 26; opd 27` vs `opd 27`) before touching code.
3. **lensarr trace-time overrun** (PLAN §0 line 42): tracing a `LensArrayIndRef=`
   deck stamps a lenslet index onto a later element's `IndRef` (`tst_save_keys.in`:
   load→save clean, load→opd→save shows `IndRef 1.0 → 1.51242597`).  Chase the
   lenslet index array bounds on the trace path (`lensarr_indexes.inc` + the
   LensArray branch in tracesub/propsub).  Gate: CLI load→opd→save on
   `tst_save_keys.in` byte-identical to load→save except the trace-state keys, both
   compilers; and `tst_save_keys` no longer traces to NaN if the overrun was the
   cause (if NaN remains, say what it is — a second item, not this one).
4. **`ifLNsrf` root pick in `RefSrf` / `ObsSrf` / `PolElt` (elemsub.F) and `IntSrf`
   (didesub.F)** still use the vertex-distance metric (PLAN §0 line 271);
   `LNsrfRoot` (surfsub) is the rule.  The three-case check in CLAUDE.md
   ("ifLNsrf root pick is AXIAL") is mandatory: a pupil Reference just ahead of a
   convex conic (`Rx_SchwarzschildEP.in`), the FEX sphere, an OAP after a
   Reference (tBench).  Gate: `tTraceRestart` / `tFwdRoot` extended with a
   Reference/Obscuring/PolElt-after-Return case that the vertex metric gets wrong
   (construct it; show it red first).  Corpus A/B mandatory (surfsub/elemsub).
5. **Load-time crashes in the extended corpus** — the FwdRoot A/B found five decks
   in `MACOS_sandbox/old_Rx/` that CRASH the CLI at load (`ape`, `mars`, `nngx`,
   `wfpc3`, `wideAngle`) and the manual's `SegDemo.in` that times out at load.
   Triage only, one line each: a Fortran runtime abort on malformed input (the
   ReadRealsPad class → folds into item 1), a validator gap, a real engine bug, or
   a dead legacy keyword.  A crash at LOAD kills the host when the engine is a mex;
   each one is either fixed under item 1 or becomes its own numbered item with a
   reproducer deck in `ZGD_test_files/`.

Waiting on Dave's rulings (PLAN_CONSOLIDATION §4), do NOT start: the IRIS real-grid
reload SIGSEGV (fix vs document) and the `LUseChfRayIfOK` global default.

## Rules that bite here

- Rule 1: a fix ships with a gate that is RED on the pre-fix binary; keep the
  pre-fix build dir (`build_release_prefix/` from the parent commit) for the
  A/B and the red leg, and say its SHA in the commit.
- Rule 5: guards warn once per run and count; a hard error needs Dave's word.
- Rule 6: nothing a command or start-up sets as a user option goes into
  `elt_mod_init_vars`.
- Rule 8: fixed-form ≤ 72 columns.
- Rule 9: commit locally, SHA + branch in your message to CC; push only on Dave's
  word.  No `--amend`.  Never relink the shared mex while another MATLAB runs.
- Report to CC after each item (the cross-session channel): the gate's red/green
  counts on both compilers, the A/B's moved-deck list, the fast-suite count.
