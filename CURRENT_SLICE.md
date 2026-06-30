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

- **Sprint / item:** Per-segment GridMat generator + per-segment influence path (mmacos sensitivities) **and** build-system fixes (issue #56) — both LANDED + pushed.
- **Plan anchor:** PLAN.md §0 (Surface=Zernike-Reference TODO; ApStop SAVE; GridData closed); memories [[project-gridmat-generator]], [[project-build-fixes-issue56]]; MACOS_resources/mmacos/CLAUDE.md.
- **Branch / worktree:** sls-dev (both repos). Build-break + GMI-build fixes ALSO on opt-dev (cherry-pick).
- **Definition of done (honest):**
  - [x] Issue #56 (smacos_dvr own `Fortran_MODULE_DIRECTORY`) — sls-dev `84f0b9c`, opt-dev `cf01617`; GMI FC-stamp `5c729f0`/`c2f3c19`; 4 dead `makefile_*.sh` deleted + makegfortran GMI report + HOW_TO_COMPILE cleanup `7041420`. Pushed.
  - [x] iElt fix (dw_dx_multi / dw_dz_zernike_multi propagate canonical `out.iElt`) — sls-dev `0f071f2`, pushed
  - [x] GridMat generator (`segment_grid_basis` + `write_grid_file` + per-segment `influence` path in grid_channels/dw_dgrid[_multi] + `gen_segment_gridmat` example) — sls-dev `245af94`, pushed; validated
  - [x] PLAN.md §0 `Surface=Zernike` for `Element=Reference` engine TODO — macos sls-dev `d8dd5b8`
  - [ ] Deferred follow-ons (Next step below)

## In-session state NOT yet committed
> Everything above is COMMITTED + PUSHED. Remaining uncommitted (deferred):
- The OTHER `run_dwd*_multi` self-contained examples + `plot_dw_per_element.m` (mmacos/sensitivities) — STILL uncommitted (only the `gen_segment_gridmat` example + `mmacos_setup.m` shipped, in 245af94). Prune cruft (stray `SD3ff.in`, experimental `run_dwd{x,z}_multi_SegDemo/`) before committing.
- These compaction doc edits (CLAUDE / MEMORY / CURRENT_SLICE / PLAN).

## Just tried / ruled out (with why)
- Per-segment GridMat SINGLE-vs-PERSEG **94%** on SegDemo3conic was NOT a bug — bases ~98% congruent (3-fold tilt pattern, opposite segments identical); 94% = the **max** at ~100 mask-EDGE pixels (one mask has the pixel, the other doesn't → full poke vs zero). Per-segment matters for clipped EDGE segments. [[project-gridmat-generator]]
- Masks render as Voronoi **WEDGES** (not hexes) — correct: `ApType=None` segments → engine assigns each ray to the nearest centre. DON'T "fix" them.
- masks.png rendered **grayscale** (headless imagesc-of-binary + colormap trap) → DROPPED it; the basis montage already shows the mask shapes. [[feedback-headless-png-grayscale]]
- **DON'T push generated figures into the Rx** — Dave: the generator makes the per-segment mode BASIS for `run_dwdgrid*`; the Rx-collapse (modes×coefs in the engine) is a deferred bigger task.

## Next concrete step
- Wire the `run_dwdgrid*` examples to consume the per-segment `.mat` (grid_channels now takes the `segment_grid_basis` struct directly); curate + commit the remaining `run_dwd*_multi` examples + `plot_dw_per_element` (prune cruft, sls-dev). Deferred engine: `Surface=Zernike` for `Element=Reference` (FreeForm-segmented apertures); Rx-collapse; more Zernike types.

## Open micro-questions (slice-local)
- opt-dev cherry-pick the iElt fix (`0f071f2`)? It's a bug fix; the generator is sls-dev feature work.
- Consolidate the bespoke `run_dwdgrid_multi_SegDemo` (gs single-basis) with the new per-segment generator path?

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
