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

- **Sprint / item:** <e.g. Sprint 2A-ii — `macos.design.Telescope` 2-mirror builder>
- **Plan anchor:** <PLAN_DESIGN_LAYER.md §8 Sprint 2A-ii ; §5.2 ; §6.2>
- **Branch / worktree:** sls-dev (macos) + sls-dev (MACOS_resources) @ <short-sha>
- **Definition of done (honest):** <the checkbox is only tickable when ALL of these are true>
  - [ ] <e.g. Mirror component contract (§6.2) implemented>
  - [ ] <e.g. closed-form Kr/Kc for Cassegrain-class (§5.2) matches fixture>
  - [ ] <e.g. emitted .in passes ValidatePrescription>
  - [ ] <test(s) green in SUITE_FAST>
  - [ ] <worked-example script runs end-to-end, exit(0)>

## In-session state NOT yet committed
> What exists only in this conversation / working tree right now and would
> be lost at compaction. Keep terse; this is a ledger, not prose.
- <file:line — what changed, why, committed? Y/N>
- <decision reached this session, not yet in §10 — e.g. "power-split sign:
  using concave-M1 alternate per Open #5; revisit if layout fights it">

## Just tried / ruled out (with why)
> The expensive-to-rediscover cul-de-sacs. Stops the post-compaction redo loop.
- <approach → outcome → why abandoned. e.g. "tightened SQP tol to chase
  gradient noise → no help → fall back to patternsearch per §1.3(3), NOT
  tighter knobs">

## Next concrete step
- <the single next action, specific enough to resume cold>

## Open micro-questions (slice-local)
> Not §10 Open items (those are durable). These die when the slice lands.
- <e.g. "does e5hex1 fixture already carry a BFD>0 case or add one?">

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
