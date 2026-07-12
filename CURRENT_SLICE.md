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

**2026-07-12: PLAN_DESIGN_LAYER Sprint 2D — SEGMENTATION + SENSING
(Dave's objective): segment the design flow incl. per-segment local
coordinate systems, edge sensors, laser metrology; MET truss from M2
to points around M3; future = optimize the MET configuration.**

### The frame (Dave, 2026-07-12 — supersedes §6.6 tier-3 draft merit)
Deliverable = the **linear forward model of a segmented system on one
shared per-segment DOF vector x**, for closed-loop control analysis:
- optics `w = dwdx·x + dwdz·z + dwddm·dm + w0` (channels exist)
- edge sensors `e = dedx·x + e0` (dedx = SegMirMaker `Hx`)
- laser MET `l = dldx·x + l0` (poke → METcalc → FD; NEW channel)
- simple control `dx = −pinv(dwdx)·(w − wtarg)`, w from the ESTIMATE.
**Configuration merit = post-control wavefront residual.**  With
full actuation the correctable part cancels exactly and
`w_post = dwdx·(x − x̂)`, so
`merit = E‖w_post‖² = trace(dwdx · P_δx · dwdxᵀ)`,
`P_δx = X − X·Hᵀ(H·X·Hᵀ + R)⁻¹·H·X`, `H = [dedx; dldx]`, X = prior
cov of x (deploy tolerances), R = sensor noise.  Prior wavefront cov
`W = dwdx·X·dwdxᵀ` is the baseline to beat (report the ratio).  When
actuation ⊄ states, report the uncorrectable projection
`(I − dwdx·pinv(dwdx))·w` separately.  Unobservable directions
saturate at their prior wavefront cost (no singular inverse).  MET
placement optimization (tier 3, future) minimizes this merit over
launcher/fiducial placement + beam topology — pure MATLAB on the S4
outputs.

### Decisions (Dave, 2026-07-12)
1. Dev/test parent = **e5mono import** (SegMirMaker's canonical parent,
   committed references); showcase on the design-layer 3M after.
2. **Hx normal-height edge-sensor model accepted** (relative
   surface-normal displacement at edge midpoints — piston+dihedral, no
   in-plane gap/shear); richer sensor models later = pluggable dedx
   backends behind the same Jacobian contract.
3. **MET default geometry = Stewart-platform trusses**: 6 launchers on
   EACH segment illuminating 3–6 fiducials on M2; 6 launchers around
   M3 illuminating 3–6 fiducials on M2; one measurement per
   launcher/fiducial pair = change in straight-line distance;
   ≥ as many measurements as DOFs.

### Slices
- [x] **S0 — SegMirMaker refresh + batch driver (Q7).  DONE 2026-07-12,
  commit pending fast-suite gate.**  CMakeLists repointed
  build_release_giza→build_release, npsol/lapack/blas→slsqplib;
  **`-fp-model strict` added → a rebuilt binary reproduces the
  committed test_in references BYTE-IDENTICALLY** (first non-strict
  build differed at 1 ulp; strict closed it; verified 2× fresh-dir
  sha256).  Batch mode = scripted stdin (NO control-file mode exists;
  README's "nine questions" is really ~15–17 prompts; DOF default is
  3 not the documented 6 — README fixed).  Committed answer fixtures
  `test_in/e5pie.stdin`/`e5hex2.stdin`.  MATLAB driver
  `macos.design.segmirmaker_run` (fresh scratch dir; copies parent +
  macos_param.txt + GridFile= refs; positional answers; 'Done.' gate)
  + `tSegMirMaker` (3/3: pie + hex2 byte-identity, size-args error)
  wired into SUITE_FAST.  GOTCHA: rerunning in a used dir trips the
  overwrite prompt and shifts the answer stream — always fresh dir.
- [ ] **S1 — `segment('M1',…)` splice.**  Emit parent Rx → S0 driver →
  splice .presc blocks in place of the parent elt + downstream
  renumbering (emitter owns numbering; watch ApStop/element-index
  refs) → reload → validate (ValidatePrescription + segment-center
  ray check ≈ SEGRAYTRACE's: center-ray endpoint vs RptElt, tol
  ~1e-11).  Spec carries per-segment frames; `segment_frames()`
  readback.  Segment local csys = RptElt + TElt/pMon,xMon,yMon,zMon
  triad (SEGMENT = full parent surface: same conic/FF/psi/Vpt every
  segment; Mon slot zeroed = reserved per-segment figure channel;
  DOF columns = [rot·x̂ rot·ŷ rot·ẑ | trans·x̂ ŷ ẑ] in the triad).
- [ ] **S2 — edge sensors (dedx).**  Ingest Hx.m + MeasToSeg (nMeas ×
  nState, nState=6·nSeg; row 1 = master-segment absolute piston);
  map columns onto the design-layer DOF ordering; VALIDATE with an
  independent MATLAB kinematic recompute (normhat·Tb + moment-arm
  crossmat at edge midpoints) — catches frame/order slips.
- [ ] **S3 — MET emission + engine wrappers.**  Builder `add_met(...)`:
  Stewart trusses per decision 3 (fiducial/launcher placement in the
  segment triad → global; emits nMetPos/SrfMetPos rows + tMetElt +
  metBeamFlg).  Engine (small): `met_calc` + `met_get` wrappers in
  macos_api_mod (codegen Path A → mmacos veneer `macos.met()`), and
  CAPACITY BUMP `elt_mod.F:332-336` mMetSrf 20→64, mSysMetBeam
  128→512 (19 segs + M2 + M3 = 21 surfaces; 19·6+6 = 120 beams —
  2-ring barely fits today, 3-ring doesn't).  Ship with Q8-style
  closed-form tests (gauge Δ == LOS projection of relative fiducial
  motion; orthogonal-motion null; reciprocity) = the §6.6 tier-2 gate
  for SrfMetCalc.  Document limitation: straight-line lengths, NO
  LOS/obscuration check (PLAN §4.5 ShowMetObscur deferred).
- [ ] **S4 — `dmet_dx` channel + `forward_model()` + demo.**  FD
  Jacobian via perturb→METcalc→read (dw_dx supervisor pattern; the
  API perturb path CPERTURB_PROG DOES move SrfMetPos).
  `forward_model()` = {dwdx, dedx, dldx, w0,e0,l0} on ONE x ordering
  (per-segment 6-DOF local + M2/M3 rigid).  Worked example: segmented
  e5-class primary → forward model → draw x~X, measure with noise,
  estimate, control dx=−pinv(dwdx)·ŵ → simulated RMS w_post ==
  analytic trace(dwdx·P_δx·dwdxᵀ).  That demo = sprint acceptance.

### Engine facts (2026-07-12 surveys; file:line current that day)
- MET: `nMetPos`/`tMetElt`/`metBeamFlg` parsed msmacosio.inc:1937-1975
  (nMetPos MUST precede tMetElt or parser STOPs); global-frame points
  `SrfMetPos(3,48,20)`; ONE command `METcalc` (3-char min-match,
  macos_cmd_loop.inc:2511) → `SrfMetCalc` utilsub.F:1810 = **pure
  straight-line Euclidean distance** source-point→target-point;
  flat output `metMeasBuf(1:nMetMeas)` (GMI output #9, pflg(27)).
  No LoadStack entry, no api_mod wrapper, no mmacos surface yet.
- SrfMetPos moves under CPERTURB + CPERTURB_GRP + **CPERTURB_PROG
  (funcsub.F:372-380 — the API/mmacos/pymacos path, so FD works)**;
  CPRead / CPERTURB_2 / LnkEltCPERTURB do NOT (PLAN §4.5 gaps — avoid
  those paths for met work).
- Engine `EdgeSensors` keyword = parsed+stored+SAVEd but **dead** (no
  consumer; PLAN §4.5 recommends not building on it) — dedx comes
  from SegMirMaker Hx per §6.6 tier 1.  SAVE round-trips all met
  keys since 662e86e (nMetPos precedes tMetElt in the writer).
- SegMirMaker: standalone fixed-form Fortran, own CMake vs
  build_release; loads the parent through the REAL engine loader
  (needs macos_param.txt + GridFile data in cwd);
  `SEGRAYTRACE` = engine-side QA command (macos_cmd_loop.inc:1391).

## In-session state NOT yet committed
- MACOS_resources: segmirmaker CMakeLists/makesegmirmaker.sh/README/
  CLAUDE.md edits + test_in/*.stdin + mmacos segmirmaker_run.m +
  tSegMirMaker.m + run_mmacos_tests.sh SUITE_FAST entry.  Fast suite
  running as the pre-commit gate → commit on green.
- macos: this CURRENT_SLICE rewrite only (no engine change yet;
  engine work starts in S3).

## Just tried / ruled out (with why)
- Byte-parity chase vs the retired giza-era binary — RESOLVED, not
  ruled out: `-fp-model strict` reproduces the committed references
  exactly; no reference refresh needed.
- Reusing a scratch dir for repeat runs — ruled out (overwrite prompt
  shifts the positional answer stream; driver always makes fresh).

## Parked threads (carry, do not lose)
- **eac5 thread (Dave 2026-07-07):** reflective e5 back end =
  `~/dev/tst_dir/eac5mono.in` m3 (R=0.328 m, lMon 20 mm) + m4
  (R=0.407 m, lMon 30 mm) compact mini-relay ~0.2 m past the f/21
  intermediate focus — the reference for the freeform_unobscured
  "+2" (supersedes single-M4).  File breakage: CRLF + tabs-in-TElt
  (l.146) + `/* */` comment blocks (l.151+); Phase-1 validator
  rejects ("blank line inside TElt block").  NEXT: does the parser
  know `/*`? clean a copy, load/trace, extract m3/m4 conjugates.
- **±1′ design fork (Dave's call pending):** (a) ANSI-45 mode-depth
  test via one Jacobian assembly + column-space projection (~15 min);
  (b) pupil-aware re-layout (Sprint 5); (c) re-scope field/λ.  The
  shipped freeform_unobscured 3M numbers rest on pathological
  surfaces (see [[project_zern_solve_doctrine]]); re-solve under the
  doctrine needs Dave's OK.  Probe scripts:
  `~/dev/MACOS_sandbox/freeform_4m_probes/`.
- Others' actions on the AppScan response: Luis closes #61 + deletes
  ft-dev (both repos); IT rescans + dispositions the 2 conftest FPs.

## Next concrete step
Fast suite green → commit S0 (MACOS_resources + this file on macos)
→ S1 `segment()` splice: decide splice host (Telescope method vs
System.from_rx path — likely both share one private splice helper),
renumbering rules for downstream element-index keywords, and the
segment-center ray check via mmacos trace.

## Open micro-questions (slice-local)
- S1: which element-index-bearing keywords must renumber on splice
  besides nElt (ApStop index, tMetElt targets, OptTgtElt, RefElt,
  LnkElt/group refs?) — enumerate from msmacosio.inc when building.
- S3: met-point placement API — accept explicit point lists AND the
  Stewart preset; where do M2-fiducial defaults sit (radius on M2
  face)?
- S4: x ordering — segments first (6/seg, triad order rot|trans) then
  M2/M3 rigid?  Must match dw_dx channel DOF naming.

## Promote-on-land  →  then CLEAR this file
> Same commit as the `design-sprint-N` tag: move each item to its
> permanent home, then reset this file to the empty template below.
- [ ] PLAN_DESIGN_LAYER Sprint 2D checkboxes ticked (Q7, segment(),
      post-splice validation) + NEW: S2/S3/S4 sensing items recorded
- [ ] `CORE COMPLETE <date>` blockquote added to Sprint 2D
- [ ] §10 Decisions: Stewart MET default; post-control-residual merit
      (supersedes tier-3 min-σ(H)); Hx sensor scope — all dated
- [ ] PLAN.md §4.5: met wrapper/capacity items ticked when S3 lands
- [ ] CLAUDE.md / nested gotcha captured (segmirmaker batch gotchas →
      segmirmaker/CLAUDE.md — DONE in S0)
- [ ] agent MEMORY.md: sprint-2d project memory
- [ ] worked-example committed + named (S4 demo)
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
