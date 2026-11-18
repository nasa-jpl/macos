# BRIEF to TO (macos-fd lane) — cleanup round 1: the ruled-LATER small fixes, the dead source files, the ToDo census (CC, 2026-10-22)

Dave 2026-10-22: these three are yours; CC's Fable time is short, so this brief is the hand-off.  Engine lane rules
apply throughout (`macos_f90/CLAUDE.md` + the `notes/` file for each site; every fix ships with a journal that is
RED on the pre-fix binary; the corpus A/B from the clean worktree for anything that can move a trace -- none of
these should; explicit-path commits; push only on Dave's word; plan lines in `PLAN_CONSOLIDATION.md` 1.2 as DONE).
Origin: macos 06d3f648 / resources 1ea48dc once Dave pushes CC's local stack; CC's clean worktree is at 0efc2dec.

## 1. The ruled-LATER small fixes (PLAN 1.2; each a one-branch change + a journal)
- (k) DOPD after ZABERR / ZRM reports `RMS dOPD = -1.0` (the sentinel) although two OPDs were run -- find where the
  second OPD's record is not stored for those commands; journal: `opd; zaberr; opd; dopd` reads a real number.
- (m) LSPOT draws regardless of the NONE plot type (eight `giza ... No device open` errors) -- LSPOT must honour
  the plot-type gate the other plot commands use (notes/cli.md, giza/CLAUDE.md); journal under `/null` + NONE.
- (i) RESET does not reset the plot type (STATUS still NULL after NONE then RESET; obscuration option does reset)
  -- add the plot type to RESET's list (macos_cmd_loop.inc ~346, where zcoType is reset); the cmdref REset entry.
- (v) api `save_rx` on a NEW file prints `** Unknown command = YES`: the wrapper passes `CARG(2)='YES'` as if to a
  prompt; make the api path not re-enter the dispatcher (notes/api_bindings.md: SMACOS handlers must not re-prompt).
  Gate: the pymacos / mmacos save twins' stdout no longer carries the line (a must-fail leg on the old .so).
- (dd) ImageDemo.in carries `PropType= Annulus` on its three obscurations, not a PropType name, silently accepted:
  a named refusal or a named default at load (Dave's refusal rule: it is a loader finding -- ask him which in your
  report; implement the WARN-and-count form unless he says refuse), and fix the deck (`docs/macos-manual/examples/`,
  Appendix A regenerates: `make transcripts`).
- (f) e5hex1.in: repurpose one of the existing FLAT grid files as the fixture's `flat.txt` (Dave's ruling), so the
  GPERTURB journal and FFSegDemoZern stop staging a copy; corpus A/B must read 0 movers.

## 2. Dead source files
`macos_f90/sv_surfsub.F`, `wk_surfsub.F`, `tst_zern_wf.F` each carry their own `ZerntoMon7` and are not in the
CMake lists (verify: `grep -n sv_surfsub\|wk_surfsub\|tst_zern_wf macos_f90/CMakeLists.txt makems.sh`).  If nothing
builds them, `git rm` them; if something does, say what.  A build of both compilers + `cli_tests load` proves it.

## 3. The ToDo / TODO / FIXME census (35 markers: elemsub 9, iosub 7, api 5, tracesub 3, msmacosio 3, macos_mod 3,
macos_cmd_loop 3, utilsub 2)
Yesterday's `ToDo --> FIX` in MOD was a real out-of-bounds write ((xx)).  For each marker: (a) a defect -> a 1.2
finding with the measurement, for Dave; (b) a done thing -> delete the marker; (c) a genuine future item -> move
the sentence into the relevant `notes/` file and delete the marker.  Deliverable: one commit per class, and a table
in the report (file, line, class, one line).  No fix beyond the marker itself in this item -- findings go to Dave.

Report: SHAs, gate counts, the (dd) question, and the census table.  CCMac's next pair = your tip after Dave's word.
