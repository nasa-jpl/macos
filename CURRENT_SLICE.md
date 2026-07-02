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

**No slice in flight** — the conforming-Reference slice landed and pushed
2026-07-02 (promoted to [[project-conforming-reference]]).  Pick the next item
from the loose ends or the plan.

### Last landed (2026-07-02) — conforming Reference (PASSIVE)
- **`Element=Reference` now accepts `Surface=Zernike`/`Aspheric`** so a
  conforming Reference CARRIES a Zernike basis definition (segment shapes) for
  GS-basis dev, but has **no effect on the light** (RefSrf unchanged; coeffs
  stored, never injected — I first built it ACTIVE = WRONG, Dave caught it,
  reverted; verified with-ref==no-ref OPD 9e-12).  Engine = 3 files
  (`EltSurfCompat` gate + 2 shared-parser fixes: `ZernModes` single-vs-wrapped
  IOSTAT read; `SrfTypeName(EltID→SrfType)` warning mislabel).  **macos sls-dev
  `c9fa767`** (pushed).  Example `e5hex2_refzern` (passivity + `make_gs_basis`
  + `run_dwdgrid_multi` split + `verifyall`) + `segment_grid_basis`/
  `grid_channels` exclude refs/non-segments: **MACOS_resources sls-dev
  `52d688d`** (pushed).  Other-session fixes: `b0d044c`.  **Rx GOTCHA: segment
  grid frame `pData/xData/yData/zData` must = clocked `pMon/xMon/yMon/zMon` or
  pokes don't localize.**  [[project-conforming-reference]]
- **NOT run:** pymacos regression (needs ifx build_release + pymacos rebuild;
  mmacos suite was green for the change).  Dave's `macos` alias (ifx
  build_release) rebuilt to current/passive.

### Landed 2026-07-01 — pointers for the next session
- **dwdgrid segmented examples** (MACOS_resources sls-dev `c07ed65` + `8beb774`,
  pushed): `run_dwdgrid_multi_multisegbasis` (per-segment `segment_grid_basis`
  struct fed as `dw_dgrid_multi` `influence`; verified 72 chan / 0.0% inter-seg
  overlap / 92px spread) + `run_dwdgrid_multi_singlesegbasis` (single shared
  `gs_zernike_segment_basis`) + generic `run_dwdgrid_multi` + `run_dwd{x,z,surf}
  _multi`; both segmented library drivers (generic RX='' flavor) +
  `plot_dw_per_element` + `gs_zernike_segment_basis` committed.  Experimental
  SegDemo scratch dirs deleted; unique FreeForm/Zernike fixtures preserved in
  `~/dev/MACOS_sandbox/segdemo_fixtures/`.  [[project-gridmat-generator]]
- **GMI build fix** (macos opt-dev `796dc51` + sls-dev `2f4948e`, pushed): opt-dev
  `makeall` `-DBUILD_GMI=ON`→`OFF` (Makefile is sole GMI builder) + MATLAB_ROOT
  autodetect (both branches, was hardcoded R2025b) — fixes Scott's GMIG.F
  `fintrf.h` fail (opt-dev double-build).  [[project-build-fixes-issue56]]
- **iElt fix** cherry-picked to MACOS_resources opt-dev (`0f071f2`→`15603a8`).
- All four branch-refs (both repos × sls-dev/opt-dev) in sync with origin.

### Open loose ends (next candidates)
- **Design layer NEXT (planned 2026-07-02): Sprint 5 — simultaneous
  focal + pupil optimization** (`PLAN_DESIGN_LAYER.md` §8 Sprint 5).
  Six slices, all MATLAB over existing XPS/`pupil_quality`/
  `check_clipping`: distortion metric → `pupil_quality_multi` →
  pupil-station clearance → stacked-residual objective wiring →
  3+1 field-mirror builder support → `design/tma_3plus1/` worked
  example (null sz_tma's +1.67 mm pupil defocus / 1.77 mm astig while
  holding WFE diffraction-limited).
- ~~Deferred engine: `Surface=Zernike` for `Element=Reference`~~ DONE
  2026-07-02 (passive), see above.  Remaining:
  Rx-collapse (modes×coefs in the engine); more Zernike types.
- (Pre-existing) `tma_onaxis` designer + the convex-secondary REVERT question
  ([[project-design-drivers]]); `test7.in` / `zmode_end` rename
  ([[project-engine-fixes-lega-shipped]]); layout realizability
  ([[project-layout-realizability]], local unpushed commits).

## In-session state NOT yet committed
—

## Just tried / ruled out (with why)
—

## Next concrete step
—

## Open micro-questions (slice-local)
—

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
