# BRIEF to CCMac — the pre-release checks on the Mac (Dave, 2026-10-18)

Phase 4 of `PLAN_CONSOLIDATION.md` (§3, "Release gate") needs the three suites green on a
clean clone.  Linux is CC's (ifx and gfortran, the clean worktree `~/dev/macos_cli_base`, the
Windows configure smoke).  **The Mac is yours**: the same gate on Apple Silicon, from a clean
clone, and `MAC_PORT.md` rewritten as the record of that run (today it still says "nothing
validated" in its first half and carries the 2026-07-22 hardware pass in its second).

Wait for the SHA.  CC names one pair (macos `dev-candidate` + MACOS_resources `dev-candidate`)
after each Linux gate goes green; run on that pair, never on a tip you pulled yourself.  Open
engine items are still landing (today: CALIB keys in SAVE, the output-coordinate proposal, the
numeric PropType resolution), so a run on an unnamed tip measures a moving target.

## The checks, in order (report each with its counts)

1. **Clean clone, both repos**, into a new directory (not `~/dev/MACOS_sandbox`, not your working
   tree): `git clone` + `git checkout <sha>` for each.  Nothing from your existing trees, no
   `addpath` outside the clone, no `MACOS_HOME` pointing elsewhere (the one legitimate setting is
   the clone's own `macos_f90/`, for the `macos_param.txt` lookup — MAC_PORT's load-bearing gotcha).
2. **Engine and CLI, gfortran** (ifx has no Apple Silicon arm): `source ./makems.sh release
   gfortran`, then the CLI target (`-DBUILD_MACOS=ON`, `--target macos`; XQuartz present).  Report:
   the compiler versions, `nm` showing `__macos_api_mod_MOD_` symbols in `libsmacos.a`, any warning
   class that is new to clang/gfortran-16 since July.  A build fix is a Darwin-guarded build-file
   change only; an engine source change on the Mac goes back to CC with the compiler message, not
   into the tree.
3. **The CLI suite** from the clone: `cli_tests/run_cli_tests.sh load` on `corpus_core.txt` and
   `cli_tests/run_cli_tests.sh cmd` (the pty driver sets `TERM=dumb` itself; the two `must_fail`
   journals must stay red).  Compare with CC's Linux gfortran records for the same SHA
   (`cli_tests/records/{load,cmd}_gfortran_<sha>.csv`): a deck or journal that differs only on the
   Mac is reported by name with its transcript (`<out>_transcripts/<name>.txt`).  Expect the
   formatted output to match Linux gfortran to the digit; a list-directed difference is a finding.
4. **mmacos**: `rm src/mmacos.mexa64` then `./run_mmacos_tests.sh fast` from the mmacos root
   (the Darwin Makefile arm, `mexmaca64`; `gen_mex_wrappers.py` only if `macos_api_mod.F90`
   changed since the last mex).  Report pass / fail / skipped; the July pass was 226/0/3 with
   `tRunCompare` skipped (SegMirMaker under `build_release_ifx/` cannot exist on arm64 — say so
   again rather than fix it).
5. **pymacos**: the cmake path with the Darwin arm (`FC=gfortran`, Python 3.13 venv), then the
   whole `pytest` from `pymacos/tests/` with the PROPER phases excluded (CC's Linux line is
   6924 + 59 on 2026-10-18; yours is the first full Mac run — July ran one import/trace/opd check).
6. **Graphics**: one plot journal on `/xw` under XQuartz and one on `/png` (the giza device), e.g.
   `cli_tests/cmd/PLOT_TYPES.jou` and `tools/figs/FIG56.jou`; a PNG that renders black-on-white
   with its labels is the pass.  Note finding 1.2 (t) (the multiplication sign on PNG labels) is
   known on Linux; confirm whether the Mac shows it too.
7. **Self-containment**: the whole run above from the clone with your `~/dev` trees renamed or
   off the path for the duration; anything that then breaks is a hidden dependency — report it.

## The record

- Rewrite `MAC_PORT.md`: the first half ("prep notes", "nothing validated") becomes a short
  "how to build on a Mac" that matches the files as committed; the validation section becomes
  THIS run (date, SHA pair, toolchain versions, the seven counts, the deltas from July).  Keep
  the July gotchas that still bite; drop the ones the tree now handles.
- `cli_tests/records/mac_gate_gfortran_<sha>.md`: the counts and the named differences, one line
  each, so CC can put the Mac line beside the Linux lines in PLAN_CONSOLIDATION §5.
- Commits on `dev-candidate` of each repo, explicit paths only (no `-a`), no `--amend`, no push —
  Dave pushes on review.  State the SHAs in your report.

## What CC needs back

One message: the seven checks as a table (check, count, Linux count, differs?), the named
Mac-only differences with one line each, the build-file changes you made (file, guard), and the
two record files.  A check you could not run is listed as "not run: <why>", not omitted.
