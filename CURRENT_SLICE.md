# Current Slice — in-flight working state

> **CC:** This file is the *ephemeral* working memory for the ONE sprint
> slice in progress.  It is the deliberate complement to the permanent
> record: `PLAN.md` / `PLAN_DESIGN_LAYER.md` hold *landed* state (sprint
> checkboxes, `CORE COMPLETE` blockquotes, the §10 Decisions ledger);
> the agent `MEMORY.md` holds durable learnings; `CLAUDE.md` (+ nested)
> holds mechanical gotchas.  THIS file holds only the half-done middle
> that compaction throws away — and it is **promoted then cleared** when
> the slice lands.  It is never a second source of truth.

> **CC (post-compaction / session resume):** read this file FIRST, then
> the plan section it anchors, then the CLAUDE.md set per the root
> directive.  If this file is at the empty template below, no slice is
> in flight — pick the next unchecked item from the active sprint.

> **CC (standing conventions, do not relax here):** work lands on
> `sls-dev` (macos) + companion `sls-dev` (MACOS_resources); every new
> function ships with a matlab.unittest test; `./run_mmacos_tests.sh
> fast` between edits, full suite pre-commit; every `matlab -batch`
> ends in `exit(0)`; each sprint tag ships a runnable worked-example.

---

## Active slice

- **Sprint / item:** GridData (SrfType-9) grid-frame engine fix + the dw/dgrid SegDemo example (mmacos sensitivities)
- **Plan anchor:** macos PLAN.md §0 (GridData closed item); MACOS_resources/mmacos/CLAUDE.md "dw/d(grid) segment example"; memory [[project-dwdgrid-seg-example]]
- **Branch / worktree:** sls-dev (macos) + sls-dev (MACOS_resources). Engine + build fixes LANDED on sls-dev AND opt-dev.
- **Definition of done (honest):**
  - [x] Engine GridSrf null-frame fix — sls-dev `03db580`, opt-dev `1b535a5` (cherry-pick), pushed
  - [x] GMI Makefile slsqplib link — sls-dev `af68528`, opt-dev `de85fa0`, pushed; GMI regression 6/6 bit-identical
  - [x] SegDemo3conic dw/dgrid images localized (GridSrfdx=0.01); generic `plot_dw_per_element` (center+multi) wired into all run_dwd*_multi examples
  - [ ] Examples committed (UNCOMMITTED — pending curation)
  - [ ] Cleanup items (Next step below)

## In-session state NOT yet committed
> Engine + GMI-build fixes are COMMITTED + PUSHED (both branches, both repos). Uncommitted = mmacos example feature work + this session's doc edits:
- `mmacos/mmacos_setup.m` (new), `mmacos/sensitivities/plot_dw_per_element.m` (new) — uncommitted
- `mmacos/sensitivities/examples/` tree (run_dwd*_multi self-contained examples + SegDemo3conic + run_dwdgrid_multi_SegDemo3; SegDemo3conic.in @ GridSrfdx=0.01) — uncommitted; HAS CRUFT to prune: stray `SD3ff.in`, experimental `run_dwd{x,z}_multi_SegDemo/` dirs
- Doc edits (uncommitted in working trees): `macos_f90/CLAUDE.md`, `GMI/CLAUDE.md`, `mmacos/CLAUDE.md`, `sensitivities/README.md`, `PLAN.md` (+ Dave's new ApStop SAVE to-do, §0), `MEMORY.md` + the project memory

## Just tried / ruled out (with why)
- GridData piston root cause: NOT EP-conjugate, NOT grid→EP offset, NOT fex-vs-sxp, NOT basis (GS vs circular) — ALL ruled out across a long debug. It WAS the engine `GridSrf` null grid frame (center-pixel-only). DON'T re-chase the wrong hypotheses.
- Decisive test = the FAITHFUL `dw_dgrid` pipeline (FreeForm == GridData at matched GridSrfdx); one-off `m.trace()`+`m.opd()` probes gave setup-dependent artifacts — use the real pipeline, not hand-rolled traces.
- `zern41em5z155em3` figure diverges at GridSrfdx=0.0071 for FreeForm too → figure-scale, NOT a GridData bug; 0.01 is the figure's design scale. Flat segments rejected (Dave: unusable in a real telescope).

## Next concrete step
- Curate + commit the mmacos examples: prune experimental `run_dwd*_SegDemo` dirs + stray `SD3ff.in`, then commit `mmacos_setup.m` + `plot_dw_per_element.m` + `sensitivities/examples/` to MACOS_resources sls-dev (feature work, NO opt-dev cherry-pick). Then deferred cleanup: co-locate the GridFile, retire the bespoke `mimg`+band-aid loop / consolidate the two SegDemo drivers, standalone G-S generator (saves GridMat).

## Open micro-questions (slice-local)
- Consolidate the two SegDemo grid drivers (generic `run_dwdgrid_multi_SegDemo3` vs bespoke GS `run_dwdgrid_multi_SegDemo`) into one, or keep both?
- GridFile co-location: copy the 917 KB `zern41em5z155em3.txt` into each example dir, or a shared fixtures path?

## Promote-on-land  →  then CLEAR this file
> Same commit as the `design-sprint-N` tag: move each item to its
> permanent home, then reset this file to the empty template below.
- [ ] PLAN checkbox(es) ticked: <which>
- [ ] `CORE COMPLETE <date>` blockquote added to the sprint: <one-line summary>
- [ ] §10 Decisions: <new "Made" entry / resolved "Open" item>, dated
- [ ] CLAUDE.md / nested gotcha captured (if a new trap was found): <where>
- [ ] agent `MEMORY.md` learning (if workflow-level): <one line>
- [ ] worked-example committed + named: <example_*.m>
- [ ] **reset CURRENT_SLICE.md to empty template**

---

## Empty template (reset state — copy over the above on land)

```
## Active slice
- Sprint / item: —
- Plan anchor: —
- Branch / worktree: sls-dev + sls-dev @ —
- Definition of done (honest): —

## In-session state NOT yet committed
—

## Just tried / ruled out (with why)
—

## Next concrete step
—

## Open micro-questions (slice-local)
—

## Promote-on-land → then CLEAR this file
- [ ] PLAN checkbox(es)
- [ ] CORE COMPLETE blockquote
- [ ] §10 Decisions entry
- [ ] CLAUDE.md / nested gotcha (if any)
- [ ] agent MEMORY.md learning (if any)
- [ ] worked-example committed
- [ ] reset this file
```
