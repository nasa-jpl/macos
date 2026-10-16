# State sync for CCMac — 2026-08-19 (from Dave + CC-Linux)

Purpose: bring your local memory files current.  Everything below is
landed fact unless marked pending.  Hashes are origin truth.

## Branch state (both repos) — update any memory that says otherwise

- **`dev-candidate` is the current consolidation branch in BOTH repos,
  pushed:** macos = `485160a`, MACOS_resources = `8c172ff`.
  Each is a strict descendant of its `dev`; Andy will fast-forward
  `dev` to these tips, then run the dev→main release steps
  (MERGE_HANDOFF_ANDY.md carries a status note: steps 1–2 DONE).
- PR #67 (pol-core → dev, the vector-polarization arc) was merged by
  Andy as a SQUASH (`3b9dfc8` on macos dev).  pol-core is therefore
  content-complete on dev but its commits are not ancestors — the
  branch must never be merged again.
- **Deletable branches** (worktrees removed first; Andy or Dave
  executes): macos `pol-core`; resources `pol-core` (last commit
  2db52a4 = overcoat-chromatic ladder, merged as `90e7d19`),
  resources `pol-ifo` (merged as `e462c4f`), local
  `pol-ifo-merge-prep`.  **`develop` (Norbert's)**: fully accounted —
  the zrn_freeform arc was ported long ago (11dec8e/fa093a2, our copy
  strictly better) and the one unique piece (spot() return-code check
  cf027c5) is cherry-picked as `8c172ff` with his authorship.  Tag
  `archive/develop-final` before deletion.

## What landed on dev-candidate (beyond #67)

- **Segment-routing audit round 2 closes**: ff_hex2 regenerated (was
  the corpus's one 180°-permuted deck; both lineages now
  byte-identical), tSegMirMaker byte-identity references committed
  (were gitignored — never tracked), `assert_draw_frame_global` guards
  the three Telescope.m b.U-as-global sites (errors, orientation test
  at cos 8.1° — equality false-alarms on off-axis field points).
- **GridFile silent-flat family closed (Luis's GridData report)**:
  (1) parse-order fix — `GridFile=` before `nGridMat=` used to read a
  0×0 grid silently; now deferred until nGridMat is known, loud error
  if never set; (2) name buffer widened 24→256 chars (GridFile= AND
  AmplFile= — 25-char names used to truncate → "does not exist" →
  silent flat).
- **rodgers1 post-pupil-fix regen COMMITTED** (`bd4af5d`) — the
  2026-08-01 regen that sat dirty for 17 days.  PACKET.md Addendum-10
  rung-3/4 columns still quote pre-regen values; reconciliation is in
  the demo rewalk (below).
- **dw_dgrid 'elts' filter fixed** (was computed then discarded — all
  16 columns at 8× FD cost) + tDwDgridElts.
- **Two test re-pins after the merges, both test-side, no engine
  defect**: tPolElement's bit-exact pins were minted on the pre-#70
  pol-core engine — the merged engine shifts them by the uniform
  ColSource pupil factor (1.01577535231445 amplitude, its square in
  intensity) and lands the Malus pins on the EXACT closed forms
  [1, 3/4, 1/4, cos²(0.37)]; tViewRx's Return-label match was string-
  literal and pol-ifo's offset-label rework changed the format.

## Gates on the consolidated pair (all green)

- mmacos full suite 593/593 (7 model-size groups, gfortran, Linux).
- polfloor 15/0 (256) + 6/0 (512) including the new overcoat ladder.
- pymacos: main 6694/0 (1 skip) + PROPER-compare 26/26, ifx build.
- GMI 6/6 untouched.

## Action items for you (CCMac)

1. **Re-pull both repos to `dev-candidate`** and rebuild: engine
   (gfortran), then delete the mex to force the relink.  Your Mac
   memories that anchor on `dev` tips are stale until Andy
   fast-forwards; anchor on dev-candidate for now.
2. Your queued item stands: **dw_dx_multi per-field EP reset**
   (NOTES/TO_CCMAC_2026-08-03_1416.md — now on dev via #67; default +
   regen decisions are Dave's).
3. **Heads-up: the mmacos examples tree is being REORGANIZED this
   week** (Keysight CodeV demo ~2026-09-01, they browse the repo):
   `examples/` + `design/examples/` collapse into numbered
   `templates/10_…90_…` + `challenges/` (rodgers1/2, afocal4).  Plan:
   `mmacos/PLAN_EXAMPLES_REORG.md` (rev 2).  **Do not land work into
   the old example paths this week** — coordinate, or base on the
   `examples-reorg` branch once it exists.  The `run_dwd*` runner
   names/call forms are a preserved surface (Dave).
4. Two Linux-side operational facts worth mirroring in your memory:
   pymacos suite must run with cwd = `pymacos/tests/` (context.py
   builds the import path from `Path(".")`; wrong cwd → ~6600 bogus
   FAILs); the only pymacos venv (with PyPROPER3) lives in the
   pol-core worktree and is symlinked — it gets recreated when that
   worktree dies.

## Near-term calendar

- ~2026-09-01: Keysight (CodeV team, ex-ORA; Mike arranged) demo of
  the design layer, near-release quality.  Week 1 = reorg (TO),
  Week 2 = rewalk of every design example on the new tree (expect the
  documented post-#70 count shifts in regenerated artifacts —
  expected, not regressions).
