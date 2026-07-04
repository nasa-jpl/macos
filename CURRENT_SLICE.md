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

**SAVE audit element-data bucket (PLAN §0 item 2 follow-through)** —
2026-07-04, engine work COMPLETE, gates: GMI 6/6 PASS, mmacos full
suite RUNNING (mex relinked).  On green: commit macos sls-dev
(iosub.inc + macos_IO.f90 + msmacosio.inc + utilsub.F +
lensarr_indexes.inc + SAVE_KEYWORD_AUDIT.md + PLAN.md + fixture
tst_save_keys.in/tst_save_ampl.dat) — Dave steering interactively,
ask before push unless he pre-authorizes.

What shipped in the working tree (all emission gated on
value-set-in-memory; e5hex1 legacy SAVE byte-identical):
- All 18 element-data keys round-trip: Coating (thickness un-scaled
  ×IndRef/Wavelen), Grad*, Doe*+OrderHOE (new DoeTrGrating CASE),
  Ampl*, LensArrayIndRef (covers XYIndRefFile), ArrIndRef +
  ArrWaveLen (λ2..λn), SegAp*, ZernCenter/XDir/YDir/Rad,
  ZernAnnularRatio (**writer misspelled `ZernAnnualRatio` since
  birth — ratio silently lost every reload**; + FF/Mon emission),
  ZCOZernType, GridSrfOrder, lData, nMetPos/tMetElt, EdgeSensors.
- Structural: pData..zData widened to any nGridMat>0 (SrfType 9/11);
  nGridMat/GridFile/GridSrfdx emitted for grids on NON-grid SrfTypes
  (iris_dp_ZGD Conic NSReflector segments lost their whole grid).
- **lensarr_indexes.inc heap stomp fixed**: 107×107 rec table vs
  mLenslet=250 → LensArrayIndRef/XYIndRefFile parse corrupted memory
  (symptom: garbage nObs/ApType/PropType/nECoord on OTHER elements,
  SAVE "output statement overflows record"); table now sized from
  mLenslet (ndim=15) + clamp in InitLensletIndexAndCtr + count guard
  at the parse site.
- FmtD: '-0.0E+00' → '0.0E+00' (cross-compiler byte-identity).
- Verified: tst_save_keys.in load→SAVE→reload→SAVE byte-identical
  ifx AND gfortran, ifx==gfortran; iris converges after known 1-ulp
  psiElt DUNITIZE settle (pre-existing class, like ChfRayPos re-aim).
- REMAINING for Dave: Opt*/CALIB-family + trace-state singles SAVE
  policy; Lou-UpdateNotes cross-check.

Other open follow-ons: `Glass=` BK7 fixture (zero coverage); mmacos
opd-veneer clean error on the 9.9999e36 sentinel; FEX
footprint-autoswitch arm unexercised / journals unregressed; eac2 +
obscured-rays parked at PLAN_DESIGN_LAYER Sprint 5 head.  Dave: NO
opt-dev cherry-picks for now (incl. the older 60f886d).

### Previous state (2026-07-03, all pushed)
PLAN §0 items 1–5 ALL LANDED + PUSHED: #1 FEX `1139f78`, #2
SAVE/ApStop `96696fa`, #3 comments `d2da3cc`, #5 DXCALC `65a9c7f`,
#4 glass `65b2a67` (all macos sls-dev).

### Items 5 + 4 as landed (2026-07-03 pm)
- **Item 4 (glass tables)**: `tools/gen_glass_builtin.py` (new) →
  generated `macos_f90/glass_builtin.f90` (`LoadBuiltinGlass`, 200
  glasses from the canonical macos_glass_list.txt; verified
  programmatically identical, all 200×6).  All THREE loader
  fallbacks (macos_glass / smacos_glass / rl_macos_glass .inc) call
  it instead of hardcoded Air/BK7/LAK9; a found file still
  overrides.  CMakeLists: glass_builtin.f90 in COMMON_SOURCES.
  Tested both branches (file present → unchanged; file hidden →
  "using built-in catalog: 200 glasses").  NOTE: no Rx anywhere uses
  `Glass=` — zero test coverage for glass-by-name refraction; a
  BK7-refractor fixture is a worthwhile follow-on.
- Pending gate: final combined full mmacos suite RUNNING (GMI
  already all-pass; ifx+gfortran builds + mex rebuilt with both
  fixes).  On green: ready to commit as two commits (#5 utilsub
  DXCALC guard; #4 glass bake-in + codegen + CMake) — awaiting Dave.
- **Root cause (gdb, ifx build): `DXCALC` (utilsub.F) — when EVERY
  pupil slice has zero valid rays, the slice-search loop never
  assigns `iimax` (only set when `l>lmax`, and all l=0), and the
  second pass indexes `LRayOK` from `j=(iimax-1)*npts+1` =
  uninitialized garbage.**  ifx stack garbage → wild index → SIGSEGV
  (`forrtl severe(174)`, Dave's repro); gfortran's garbage happened
  benign (CLI survived) — why it looked ifx-only.  The earlier
  `dxRef=1d0` l==0 guard was necessary but not sufficient.
- Fix: `lmax=0; iimax=1` init + skip the second slice pass entirely
  when `lmax==0` (l stays 0 → existing warning + degenerate-beam
  guards take over).  Repro Rx: `scratchpad/optiix_badstop.in`
  (optiix_dp_conic with ApStop moved to 1e6 → all rays obscured at
  elt 1); driver `scratchpad/drive_macos_opd.py` (pty+gdb, answers
  the giza device prompt with /null).
- Verified: ifx CLI `OPD 55` → "All rays are obscured or lost" +
  9.9999e36 sentinel + prompt returns (was SIGSEGV); mmacos opd on
  the same Rx leaves MATLAB alive (wrapper returns non-struct on the
  sentinel — optional mmacos-veneer polish, not engine).  GMI
  all-pass.  Full mmacos suite RUNNING.  Both builds + mex rebuilt.

### Item 3 as landed (d2da3cc):
  - Capture: `GET_EQ` whole-line-'%' branch → `RxCommentCapture`
    (iosub.inc; banner-furniture filter for idempotency); storage in
    `elt_mod` (`RxCommentLine/Anchor`, cap 500, `RxCommentCtx`,
    `LRxCommentCapture` gate — MOD dialog also calls GET_EQ and must
    not capture; disarmed at both hosts' MBFile6 exits).
  - Anchor: `RxCommentCtx` = 0 in header, iElt per element block
    (set in msmacosio.inc's element DO loop).
  - Emit: `RxCommentEmit(0)` at `PrtSourceInfo` end, `(iElt)` at
    `PrtSingleEltInfo` end — comment before `iElt= k+1` is read in
    block k → end-of-block-k emission reproduces position (verified:
    "%Verify fElt, eElt for the following…" lands right before its
    ellipsoid element).
  - Rode along: `FmtD` 16→**17 sig digits** (bit-exact real reload);
    skip `BaseUnits/WaveUnits= (default-unit)` placeholder (GET_EQ
    eats parens on reload → came back empty; idempotency breaker).
- Verified: eac2 (14 comments) all survive at correct positions;
  pass2==pass3 byte-identical EXCEPT pre-existing 2-ulp `ChfRayPos`
  re-aim oscillation (ApStop-driven recompute each load — NOT a SAVE
  defect, out of scope).  ifx/gfortran CLI SAVE byte-identical.  GMI
  all-pass.  Full mmacos suite RUNNING (commit+push on green — Dave
  gave advance permission).
- **NEXT: #5 (opd all-rays-lost SIGSEGV), #4 (glass tables)**; work
  SAVE_KEYWORD_AUDIT element-data bucket (Dave review); eac2/
  obscured-rays parked at PLAN_DESIGN_LAYER Sprint 5.

### Last landed (2026-07-03) — FEX EP-radius rework + guards (§0 item 1)
- **FEX EP radius = chief-ray distance EP→(iElt+1) plane, ALWAYS** —
  whatever iElt+1 is (FP, coronagraph mask, ZWFS…): the far-field
  propagation distance for physical optics (Dave's spec).  Sign per EP
  EltID (Return −cr1dir / Reference +cr1dir).  Legacy `zp_iEm1` =
  fallback + autoswitch alternative.  Guards (all noisy): telecentric
  (parallel probe chief rays → FindCrossPt 0/0; keep station, radius =
  station→iElt+1, FLAT 1d22 fallback); beam-footprint sanity
  autoswitch; Rx-order flag (Return before EP return should usually be
  a Reference — fires on most legacy Rx, deliberate nudge).
  `tracesub.F` FEX block + decls; SXP/XPS untouched (SXP redundant now).
- **Compat pass green:** conforming double-pass Rx → legs equal by
  construction (pre-EP Return AT the focus) → round-off no-op.
  `Cassegrain.in` −1.2295 → +6.7907 == CassWithExitPupil exactly
  (right answer where legacy was degenerate).  **`eac2_7seg` −289.7 →
  +52597.6 — material shift, Dave to review.**  GMI all-pass; mmacos
  30 classes / 0 verification failures (tProperCompareCassFF crash
  pre-existing); pymacos has no FEX tests.
- **fex_sweep tool** saved at `MACOS_resources/mmacos/tools/fex_sweep/`
  (per-Rx process isolation; corpus = old_Rx sandbox + manual
  examples).  Corpus notes: ape.in trips a legacy loader parse bug
  ("Bad real number", segfaults under ifx build — pre-existing);
  several old Rx fail load (SegDemo/FFSegDemoAll/Luneberg/eac1/…)
  or the sweep's heuristic stop (j18sa/j18sc/6MST_segV3/dmt6seg).
- **Not yet exercised:** the footprint-autoswitch arm (no corpus Rx
  has a degenerate DEFAULT leg).  Journals not regressed.
- **Parked → design layer (Dave 2026-07-03):** the eac2_7seg leg
  discrepancy (−289.7 vs +52597.6; coronagraph, mask obscures all
  rays) and the "should FEX include obscured rays?" question — see
  the deferral note at PLAN_DESIGN_LAYER §8 Sprint 5 head (FEX
  centroid mode's all-obscured path builds psip from uninitialized
  CentroidSpot — latent; vertex/radius/footprint already use
  obscured rays).  Do NOT chase eac2 further in the engine for now.
- Negative-L answer (Dave's Q): ordinary train rejects L<0 (GO TO 98
  → miss); near Reference/Obscuring/Return/LensArray ifLNsrf=TRUE
  allows it and ConSrf picks the root by |L²−mpr| PROXIMITY to the
  ray's vertex distance (surfsub.F:~148) — no flow-of-light sense.
  That heuristic is WHY a too-small reference sphere fails silently.

opt-dev cherry-pick of 60f886d still pending.

### Last landed (2026-07-02 pm) — NS grid-frame fix + ns_griddata decomposition
- **Grid figures on NSReflector segments collapsed to a center-pixel
  per-segment piston**: the NS-block Reflector call sites passed the grid
  frame `pData/xData/yData/zData` indexed by `iElt` (NS-group entry elt,
  frameless → zeros) instead of `imin` (elt actually hit).  Two-line fix,
  tracesub.F:3714 + propsub.F:983.  **macos sls-dev `60f886d`** (pushed;
  **opt-dev cherry-pick NOT yet done** — it has the same slip via 1b535a5).
- **ns_griddata example redone as a 5-variant decomposition** (Luis Q&A):
  Conic/Zernike/GridData/ZrnGrData+flat/ZrnGrData variants; superposition
  3e-8 rms (real 2nd-order slope-resampling cross-term, mostly double-pass);
  flat grid bit-identical to none.  GridSrfdx=1.1 (=280/255 — grid span is
  (nGridMat−1)·dx centered on pData; Luis's 0.2 covered only ±25 of the
  140-radius segment).  **MACOS_resources sls-dev `e6f4f17`** (pushed).
- mmacos regression: 30 classes, 0 failures; tProperCompareCassFF heap
  crash pre-existing/known.  Both release builds + mmacos mex rebuilt.
- Remaining NS gap: the `NSRefractor` ROUTINE (refractive NS) still
  null-frames GridData — iris (all-reflective) unaffected.

### Landed earlier 2026-07-02 — conforming Reference (PASSIVE)
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
