# MACOS engine core — subsystem cheatsheet

> **Post-compaction / resume:** if you are resuming engine-core work: `macos_api_mod.F90`, IO/parse (`msmacosio.inc`/`iosub.inc`/`macosio.F`), `elt_mod.F`, `macos_cmd_loop.inc`, `funcsub.F`, `surfsub.F`, `mathsub.F` (includes the mhist/readline CLI sub-prompt cache), re-read
> THIS file — nested CLAUDE.md files are NOT auto-injected after
> compaction; they reload only when CC next reads a file in this
> directory. Root rules, conventions, and the subsystem index live in
> `../../CLAUDE.md`.

*Sections below are lifted verbatim from the former monolithic root
CLAUDE.md. Move text, don't paraphrase — engine gotchas are exact.*

---

## §0 hygiene cluster (latent bugs closed on opt-dev)
- **`macos_realloc` after `macos_init_all`.**  Code paths that call
  `macos_init_all(modelSize)` directly (smacos_dvr, pymacos init,
  mmacos init) must set `macos_realloc = .true.` afterward so the
  first `SMACOS()` call re-allocates its scratch buffers
  (`L1, R1, R2, PertVec, DrawEltVec, DV1, DV2, D2, CD1, CD2, DWF`)
  for the new `mdttl`.  Interactive macos.F gets this for free via
  the main command-loop path.  Without it, model_size transitions
  (512 ↔ 1024) silently corrupt the heap on the first plot.
- **DFOURN scratch growth.**  `sunsub.F:DFOURN` had a SAVE'd
  `first_entry` LOGICAL that allocated the scratch `DATA(:)` array
  once and never grew it.  Pymacos's model_size transitions blow
  past the original `SzData`.  Fixed: re-allocate when
  `.not.allocated(DATA) .or. size(DATA) < SzData`.
- **`psiElt` renormalization in `CPERTURB_PROG` (funcsub.F).**  After
  `Q · psi_in → psi_out` the result drifts off unit-norm for large
  rotations (≥ a few mrad cumulative).  Divide by
  `sqrt(sum(D1**2))` after the DMPROD.  Caught by Src vs all-optics-
  group equivalence test.
- **`dopt_init_vars` sentinel defaults.**  `OptTgtElt=-1` and
  `OptAlg=NonLin` initialized so prescriptions that omit them parse
  cleanly through the new constrained-optim path (previously
  uninitialized → silent garbage element index).
- **`MBFile6` parser noise.**  Unknown-key catch-all for `OptBeam=`
  added to `msmacosio.inc` so older prescriptions don't trip on it.
- **ConSrf near-tangent fallback (`surfsub.F`).**  When
  `k2 = b*b - 4*a*c < TOL_TANGENT * b*b` the quadratic is effectively
  linear; fall back to `L = -b/(2*a)` instead of `sqrt(k2)` to avoid
  loss-of-precision NaN.
- **`smacos_compute.inc` slice overrun.**  Five `1:mZern` → `1:mZernModes`
  fixes (mZern≠mZernModes when both Zern and FF aspheres are active).
  Cherry-picked back to release-candidate.

## CALIB wrappers in macos_api_mod (Phase 1a/1b/1c)
- **CALIB *is* the native multi-field design optimizer** (`nls_optim_dvr`
  LM + SLSQP/NPSOL in `design_optim.F`): multi-field × multi-λ,
  FoV-weighted least-squares over per-element DOFs (8-DOF `VarElt` mask;
  **DOF 7 = radius, DOF 8 = conic**), Zernike modes, aspheric coeffs;
  targets WFE/SPOT/WFE_ZMODE (≤12 FOV × 6 λ).  `calib_run` returns
  per-(FOV,λ) WFE before/after.  This is the **`native` engine** for the
  design layer's multi-field `optimize()` (PLAN_DESIGN_LAYER §8 Sprint
  2B).  FOV set + vars + target are configured via the `.in` Opt*
  keywords today; programmatic field setters (`calib_add_fov`/
  `calib_set_wavelens`) are the deferred Phase-1d gap.
- `calib_run` + `calib_buffer_dims` (Phase 1a): drive `'CALIB'`
  end-to-end via SMACOS and expose result dims to pymacos/mmacos.
- Programmatic setters (Phase 1b): `calib_set_var_elt`,
  `calib_clear_var_elts`, `calib_set_iter`, `calib_set_tol`,
  `calib_set_target` — let bindings configure a calibration run
  without round-tripping through a `.in` prescription edit.
- `calib_clear_var_elts` carries defensive `allocated()` guards
  before deallocating optional arrays (Phase 1c) — silently no-ops
  rather than SIGSEGVing on a fresh `macos_init_all` state.
- pymacos: `m.calib(...)` plus `m.calib_set_*(...)` setters in
  `macos.py`.  mmacos: `+macos/calib*.m` plus Session methods.


## Current work: SrfType_FreeForm = 14
Composite surface: conic + Mon monomial + FF monomial + grid data.
- FreeFormSrf() and helpers (MonomialEval, MonomialD2, ConicDeltaS)
  now in surfsub.F as module procedures inside MODULE surfsub.
- MonSrf, GridSrf, MonGridSrf are thin wrappers calling FreeFormSrf.
- SetFreeFormFlags(ie) in elt_mod.F sets ifMon/ifFF/ifGridTerm at
  prescription-read time based on lMon, lFF, nGridMat sentinels.
- New elt_mod.F arrays: lFF(mElt), pFF/xFF/yFF/zFF(3,mElt),
  FFCoef(mFFCoef,mElt), ifMon/ifFF/ifGridTerm(mElt),
  FFZernCoef(mFFCoef,mElt), FFZernTypeL(mElt),
  MonZernCoef(mMonCoef,mElt), MonZernTypeL(mElt).

## FreeForm user input design
- Mon component: user enters MonZernCoef/MonZernTypeL (Zernike).
  Converted to MonCoef at trace time via ZerntoMon in tracesub/propsub/srtrace.
- FF component: user enters FFZernCoef/FFZernTypeL (Zernike).
  Converted to FFCoef at trace time via ZerntoMon.
- Both support sparse input (nXxxZernCoef/XxxZernModes/XxxZernCoef).
- Interactive dialog: mode-index + value pair loop (0=done).
- Prescription output: MonZernType/MonZernCoef + FFZernType/FFZernCoef (dense).
- SrfType 4 (Monomial) and 12 (MonGrData): still use direct MonCoef input.


## SPOT BEAM reference frame: obscured chief ray (tracesub.F LocalCoord)
- `SPOT ... beam` (and `WINDOW/PLOCATE ... beam`) build the beam coordinate
  frame from a **chief-ray trace** in `LocalCoord` (tracesub.F ~2612): a
  nominal chief ray + two differential rays (x/y translate/rotate), each via
  `CRTrace` (single ray, `SetSourceRayGrid(nMinPts=3)`).  Each `CRTrace`
  result was gated on `IF (nBadRays.NE.0)` → "chief ray becomes undefined" →
  abort.  **The on-axis chief ray of any centrally-obscured telescope is
  VIGNETTED**, so `nBadRays≥1` from obscuration alone, and beam-frame SPOT
  failed at **every** element of every obscured Rx (reproduced identically in
  the CLI; TOUT/TELT work because they use the stored `Tout`/`TElt` frame, no
  chief-ray trace).  Obscuration sets only the flux flag (`L1`/`LRayPass`),
  NOT the geometric intersection `RayPos` (CTRACE comment ~4189), so the frame
  is well-defined.  Fix: the 5 `LocalCoord` failure checks now read
  `nBadRays.NE.0 .AND. RayStatus(1).NE.RayStat_Obscured` — tolerate a
  vignetted-but-geometrically-valid chief ray, fail only on a genuine
  geometric miss/bracket.  Do **not** try to fix this by toggling `iObsOpt`:
  the obscuration `nBadRays` bump (tracesub.F ~4364) is unconditional
  (`obs_set 0` does not suppress it), and it MUST stay that way so WARN's
  per-category loss accounting fires.  mmacos regression: `tSpot` (obscured
  Cass, beam ray count == tout).

## Per-ray status tracking
- RayStat_* constants in elt_mod.F: OK(0), Obscured(1), Miss(2), Bracket(3), MaxIter(4), Undef(5).
- Per-ray arrays: RayStatus(mRay), RayFailElt(mRay), RayFailMsg(mRay) — allocated in elt_mod.
- `LOGICAL FUNCTION SetRayFail(iRay, iStat, iElt, cMsg)` records first failure only (avoids
  overwriting root cause) and returns .TRUE. iff it actually recorded a status this call.
  All six `nBadRays = nBadRays + 1` sites in tracesub.F and propsub.F are gated on the return:
  `IF (SetRayFail(...)) nBadRays = nBadRays + 1`. This prevents an obscured ray being counted
  once per propagation pass — the bug that produced the "26076 rays were lost / Obscured: 4346"
  diagnostic where 4346 × N_elements gave the inflated headline number.
- LZPFailed module variable in MODULE surfsub: set by *ZPSolve on bracket/max-iter failure,
  checked by surface routine callers (IF (LZPFailed) GO TO 98), reset after recording status.
- All STOP statements in surfsub.F converted: AZPSolve, GSZPSolve, GSZPB, SFFZPSolve, UDSZPSolve.
- All PAUSE statements converted: DZPSolve (didesub.F), ZPSolve (tracesub.F).
- Status recorded at nBadRays increment points in CTRACE (tracesub.F) and CPROPAGATE (propsub.F).
- Obscuration recorded at end of iRay loop (L1/LRayPass check).
- RayStatus/nZPFailMsg reset at trace start (alongside nBadRays=0).
- WARN subroutine enhanced with per-category breakdown (miss/obscured/bracket/other).
- Message throttling: nZPFailMsg counter with mZPFailMsg=20 threshold in MODULE surfsub.
  First 20 bracket/iter messages print; rest suppressed. WARN prints suppression summary.

## PERTURB notes
- 5 routines perturb coordinate frames: CPERTURB, CPRead, CPERTURB_GRP (funcsub.F),
  CPERTURB_2 (macos_ops.F), LnkEltCPERTURB (lnk_pert.inc).
- All now handle pFF/xFF/yFF/zFF for SrfType_FreeForm (position translates+rotates
  relative to RptElt; orientation vectors rotate only).
- pMon condition: funcsub.F uses SrfType>=4 .AND. !=10 .AND. !=11 (correct for all
  Mon-using types including FreeForm). macos_ops.F and lnk_pert.inc were <=9,
  now fixed to match funcsub.F pattern.
- pData condition: includes FreeForm so grid coord frame perturbs when nGridMat>0.
- pData condition bug fixed: was SrfType<=13 (fired for all 1-13),
  now ==12 .OR. ==13 .OR. ==SrfType_FreeForm.

## FEX EP-radius rework (2026-07-03) + SXP command (Set eXit Pupil)
**FEX now defaults to the EP→next-element radius** (Dave's spec): the
EP Return radius is ALWAYS the chief-ray distance from the EP
(CrossPt) to the **iElt+1 plane — whatever iElt+1 is** (FP,
coronagraph mask, ZWFS…) — because that is the far-field propagation
distance for physical optics.  No element-type scan.  Sign per EP
EltID (Return reverses the beam → −cr1dir; Reference passes →
+cr1dir).  Legacy `zp_iEm1` = fallback (no iElt+1 / degenerate plane)
and the footprint-autoswitch alternative.  Guards, all noisy:
(1) **telecentric** — parallel probe chief rays make FindCrossPt
divide 0/0 (Cp2→0, NO guard there); FEX detects sin<1d-12, keeps the
element station, radius = station→iElt+1 plane, FLAT 1d22 fallback;
(2) **beam-footprint sanity** — a reference sphere smaller than the
beam footprint at the EP guarantees k2<0 "surface miss" for marginal
rays (the SegDemo3 failure); autoswitches to the other leg if usable;
(3) **Rx-order flag** — a Return immediately preceding the EP return
usually marks an intermediate focus that should be a passive
Reference (pattern: Reference@FP, Return@EP, Return@FP); fires on
most legacy Rx (corpus predates the convention — deliberate nudge).
**Compatibility (fex_sweep 2026-07-03):** conforming double-pass Rx
have the pre-EP Return AT the focus, so both legs are equal by
construction → round-off-level no-op (e5hex1, 6MST, iris, j18*,
dmt6mono, CoroExample, HOEExample, lst3zern, CassWithExitPupil).
Divergent: manual `Cassegrain.in` (EP Reference right after SM)
legacy −1.2295 → new +6.7907 == exactly CassWithExitPupil's value —
the rework computes the right answer where legacy was degenerate;
`eac2_7seg` −289.7 → +52597.6.  GMI regression all-pass.  Sweep tool:
`MACOS_resources/mmacos/tools/fex_sweep/`.  Negative-L background
(why small spheres fail SILENTLY): ordinary train rejects L<0 (GO TO
98 → miss); near Reference/Obscuring/Return/LensArray
`ifLNsrf=.TRUE.` allows it and ConSrf picks the quadratic root by
|L²−mpr| PROXIMITY to the ray's current vertex distance
(surfsub.F:~148) — no flow-of-light sense.

SXP (retained, now largely redundant with FEX): FEX clone whose EP
radius is the chief-ray distance EP→FP (assumed iElt+1) — the
geometry FEX now uses by default.  Captures FP-Tz perturbations
through the EP radius and catches EP shifts when upstream optics are
perturbed (use case: sensitivity / dw/dx loops).
- `SUBROUTINE SXP` in `tracesub.F` — cloned from FEX; after
  `FindCrossPt` replaces zp with `t_FP = -((FP_vpt - CrossPt) · psi_FP)
  / (cr1dir · psi_FP)`, i.e. the chief-ray distance from CrossPt to
  the FP plane (chief ray reverses at EP, so post-EP direction is
  `-cr1dir`).  Falls back to `zp_iEm1` if iElt+1 is out of range.
- Dispatch added in `macos_cmd_loop.inc` alongside FEXIT (uses label
  8073 to avoid the existing 73 collision) — reachable from both
  interactive macos AND SMACOS, since SMACOS `#include`s
  `macos_cmd_loop.inc` from `smacos.F:210`.
- `LoadStack` entry in `smacosutil.F` — one IARG (iElt) + one CARG
  (YES/NO), same signature as FEXIT.
- pymacos wrapper `m.sxp(mode=1)` (in MACOS_resources/pymacos/dr-dev2).

Pre-existing dispatch puzzle worth knowing: `macos_ops.F:MACOS_OPS`
is NOT the SMACOS dispatcher — it's only called from the inner
optimization loop (`smacos_compute.inc`).  SMACOS proper dispatches
via `#include "macos_cmd_loop.inc"` from `smacos.F`.  When adding
new commands callable via SMACOS, add them to `macos_cmd_loop.inc`
and to `smacosutil.F`'s LoadStack (so the args get on the STACK
that IACCEPT_S reads in SMACOS mode), NOT to `macos_ops.F`.


## mathsub.F modernization
- `DMPROD(A,B,C,NA,NB,NC)` in `mathsub.F` is now a thin wrapper over
  `MATMUL(B,C)` with explicit `A(NA,NC) / B(NA,NB) / C(NB,NC)` shapes
  and `IMPLICIT NONE`. 287 call sites unchanged — the wrapper preserves
  the old flat-indexing contract (caller's physical leading dim == NA,
  which all existing call sites satisfy since the dominant shapes are
  3×3, 3×1, 1×3). Do not inline `MATMUL` at call sites: many callers
  pass flat scratch buffers (`D1(100)`, `D2(100)`) where MATMUL would
  force awkward `RESHAPE` noise, and keeping one chokepoint means a
  future DGEMM swap is a single edit. Aliasing is undefined for both
  old and new (old code zeros `A(NAPTR)` before accumulating, which
  corrupts B or C if aliased), so no regression.

## Reference
- macos_f90/Archive/surfsub_old.F : pre-FreeForm surfsub; reference for SGSrf, GSZPB, GSZPSolve
- macos_f90/Archive/ : obsolete source files (sourcsub variants, test programs, etc.)

## Name collision warnings
- macos_mod exports scalar LOGICAL ifGrid (ray-grid-established flag, used everywhere).
  Do not reuse ifGrid for any new per-element logical in elt_mod — rename to avoid
  ambiguity (e.g., ifGridTerm for the FreeForm grid-component flag).
- Inside MODULE surfsub, do not redeclare sibling functions in a local REAL*8 list.
  That overrides host association and produces an undefined reference at link time.
  Just call module-sibling functions directly without redeclaring them.

## elemsub.F call interface notes
- FindSrf, Reflector, and Refractor signatures now include FreeForm arrays
  as explicit args: pData,xData,yData,zData (grid coord frame) and
  lFF,pFF,xFF,yFF,zFF,FFCoef.
- All call sites updated: tracesub.F (3 Reflector + 1 FindSrf + 1 Refractor),
  propsub.F (3 Reflector + 1 FindSrf + 1 Refractor), srtrace.F (2 Reflector).
- FindSrf uses `use elt_mod, only:` with explicit list; Reflector and Refractor
  use only: with mMonCoef,mAnaCoef,mFFCoef,SrfType_FreeForm.
- Refractor also has MonGrData (12/13) dispatch calling MonGridSrf.

## GridSrf (SrfType-9) passed a NULL grid frame — CLOSED (sls-dev 03db580, opt-dev 1b535a5)
- `GridSrf` (the SrfType-9 thin wrapper over `FreeFormSrf` in surfsub.F)
  forwarded an all-zero grid coordinate frame (`pData_null/xData_null/
  yData_null/zData_null = 0`) to FreeFormSrf.  FreeFormSrf maps each ray into
  grid-pixel space via `xi=(xData·rhom)/dAct+diCtr`, `yj=(yData·rhom)/...`;
  with xData=yData=0 those collapse to the **center pixel for every ray**, so a
  GridData grid contributed `GridMat(center)` uniformly = a pure **piston** and
  the spatial figure was discarded.  Symptom: a GridData segment grid-figure
  poke imaged as a per-segment piston (whole-pupil staircase), not localized.
- Introduced by the FreeForm refactor (when GridSrf became a FreeFormSrf
  wrapper) — SrfType-9 grid figures were silently piston-only since.
  `MonGridSrf` (12/13) and `FreeFormSrf` (14) always forwarded the real frame,
  so FreeForm/MonGrData segments imaged clean (this is why the symptom looked
  like "FreeForm vs GridData" — it's the SrfType-9 dispatch dropping the frame).
- Fix: thread the element's real `pData/xData/yData/zData` through `GridSrf`
  (signature + decls + the FreeFormSrf call) and its FindSrf/Reflector/Refractor
  call sites in elemsub.F.  GMI regression 6/6 bit-identical (its prescriptions
  use FreeForm/Conic, not SrfType-9), confirming the change is isolated.
- **NS follow-on CLOSED (sls-dev 60f886d, 2026-07-02):** the frame-threading
  fix above used `iElt` at all ten call sites — correct at the eight
  sequential ones, WRONG at the two non-sequential Reflector sites
  (tracesub.F:3714, propsub.F:983) where the intersected element is `imin`
  and `iElt` is the frameless NS-group entry (vectors all zero → every ray
  → center pixel → per-segment piston).  Caught by the `ns_griddata`
  5-variant decomposition (mmacos examples; Luis's iris_dp NSReflector
  segments).  NSReflector elements dispatch through `Reflector` in the NS
  block, NOT the `NSRefractor` routine.  **Not yet cherry-picked to opt-dev**
  (which has the same slip via 1b535a5).
- **Deferred:** the `NSRefractor` ROUTINE (refractive non-sequential only)
  doesn't carry the element grid frame — its `GridSrf` call still passes a
  local null frame; that GridData path stays piston-only until threaded.
  Noted in-line at the NSRefractor decl.
- **GridSrfdx semantics (user-facing):** grid span = `(nGridMat−1)·GridSrfdx`
  base units, centered on `pData`; rays outside the span get `fh=0` (no
  figure).  256-grid across a 280-diameter segment → `GridSrfdx = 280/255
  ≈ 1.1`.  Too-small dx = figure only near segment center (looks piston-ish).

## Prescription I/O notes (msmacosio.inc / iosub.inc)
- The parser has two chains: label 50 (header/global) and label 61 (per-element).
  FreeForm keywords must be in BOTH chains. The element chain (61) now includes:
  lData, lFF, pFF, xFF, yFF, zFF, FFZernType, nFFZernCoef, FFZernModes, FFZernCoef,
  MonZernType, nMonZernCoef, MonZernModes, MonZernCoef.
- nFFZernCoef/FFZernModes/nMonZernCoef/MonZernModes reset per-element at loop start
  (before LOHIN2 call).
- Output (PrtSingleEltInfo in iosub.inc): pData/xData/yData/zData written for
  FreeForm only when nGridMat > 0.
- ChkDf2: FreeForm with nGridMat=0 silently satisfies pData/xData/yData/zData/lData
  (no spurious "Default used" warnings).

## ZernType parsing (ParseZernType in elt_mod.F)
- All eight ZernType / FFZernType / MonZernType parser sites — three in
  msmacosio.inc (load) and three MOD-dialog sites in macosio.F — go through the
  shared `ParseZernType(value, typeL, found, annuRatio, hasRatio)` helper in
  elt_mod.F. Don't add a new local copy; extend the helper.
- Norm-prefixed names (NormANSI, NormBornWolf, NormFringe, NormHex, NormNoll,
  NormAnnularNoll) match by 7 chars; un-normalized names (ANSI, BornWolf,
  Fringe, Noll, ExtFringe) by 3. NormAnnularNoll (typeL=9) parses an annular
  ratio from the trailing token.
- Historical bug: a "Noll → NormNoll" rewrite shim existed at one site only,
  and the un-normalized loop ran indices 1..3 (omitting 10=Noll, 11=ExtFringe),
  so `MOD ZernType=Noll` silently became NormNoll. Both fixed in the helper.
- Legacy alias: "Malacara" (pre-modernization name for ANSI) maps to typeL=1
  via a `Mala`-prefix check at the bottom of the un-normalized branch. Older
  JPL Rx files still use it; without the alias the parser
  silently dropped the field and ZernTypeL stayed 0 — which then silently
  broke the Zernike-apply trace path (see next section).

## Zernike-apply trace dispatch requires ZernTypeL ≠ 0
- `propsub.F:230-249` (and the parallel `srtrace.F:145-153`) dispatches
  `ZerntoMon1/2/3/4/6/7` based on `ZernTypeL(iElt)`. The chain has NO `ELSE`
  branch — if ZernTypeL is 0 (the elt_mod.F:944 initial value), the IF/ELSEIF
  falls through, ZerntoMon is never called, MonCoef stays zero, and the
  surface intersection (which uses MonCoef, not ZernCoef) sees no perturbation.
  Result: silent zero-response.
- ChkDf2 defaults ZernTypeL to 1 (ANSI) for elements whose parse path visited
  the ZernType default block. **NSReflector+Conic elements don't visit that
  block** — their ZernTypeL stays 0 unless explicitly set in the Rx.
- GMI workaround at `GMI.F:1066`: when forcing `SrfType(iElt)=8` to apply a
  Zernike perturbation to an arbitrary element, also force `ZernTypeL=1` if
  currently 0. Without this, GMI's pzern channel silently no-ops when applied
  to a Conic-typed element. Caught by `test03_zern_response_e2e_pie` (had been
  masked under ifx because the same broken path also corrupted memory and
  triggered an exit-time SIGSEGV — gfortran exposed the actual zero-response
  bug cleanly).
- Latent gap (not yet fixed): the propsub/srtrace IF/ELSEIF chains don't
  handle ZernTypeL=10 (Noll) or =11 (ExtFringe). A user Rx with
  `ZernType= Noll` parses correctly and sets ZernTypeL=10, then the trace
  dispatch silently no-ops. Add ELSE-with-error or extend the chain when this
  surfaces. (We didn't extend it now because no current test exercises it.)

## Conforming Reference: Surface=Zernike/Aspheric is PASSIVE (sls-dev c9fa767)
- `Element=Reference` accepts `Surface=Zernike` (8) and `Surface=Aspheric` (3)
  so a reference can CARRY a Zernike basis definition (modes/coeffs/lMon/
  aperture = "segment shapes") for GS-basis dev.  Gate: `EltSurfCompat`
  (iosub.inc) allows SrfType 3/8 for EltID=3 (was ≤2 or =9).
- **It has NO effect on the light.**  A Reference passes rays straight through
  (RefSrf sets rout=ihat; total path to the next real elt unchanged).  RefSrf
  is UNCHANGED — the Zernike/aspheric coeffs are parsed and stored but NEVER
  injected.  (Do NOT mirror `Reference`+`GridData`'s `ifOPDModGrid` `L=L+figure`
  — that's a deliberate ACTIVE phase grid, the exception.  I did that first;
  it broke Dave's "references have no effect" invariant.  Verified passive:
  with-ref == no-ref OPD to 9e-12.)
- Two SHARED-parser fixes shipped alongside (CLI/mmacos/pymacos share
  msmacosio.inc): (1) **`ZernModes=` now reads single-line OR wrapped**
  (`READ(VALUE,*,IOSTAT=ios)` all modes; only if short do the Grp-per-line
  read) — >Grp modes on one line previously pulled phantom continuation rows →
  ate the next keyword → "Bad integer".  (2) the "Invalid Element/Surface
  Combination" warning printed `SrfTypeName(EltID)` — fixed to
  `SrfTypeName(SrfType)` (2 sites msmacosio.inc + 1 macosio.F).
- Rx GOTCHA (mmacos-side, not engine): a grid/FreeForm segment's grid frame
  `pData/xData/yData/zData` must equal its clocked `pMon/xMon/yMon/zMon` or a
  per-segment poke won't localize (dW→central dot).

## Deferred PolyObsVec/Poly3DObsVec projection (msmacosio.inc + iosub.inc)
- For polygon obscuration, the 3D-vertex forms (`PolyObsVec=`/`Poly3DObsVec=`)
  used to require xObs to be parsed BEFORE the polygon line, because
  `SetCvxPolyObsVtx` projects the 3D vertices to element-local 2D at parse
  time and needs xObs/yObs to do so. Wrong order produced a hard error.
- Now the order doesn't matter. On a 3D-vertex polygon line:
  - If `xObs_FLG` is true, project immediately (legacy behaviour).
  - Else, stash raw 3D vertices into `PolyObsVtx3DPending(3, mPolySide, mObs,
    mElt)` and the count into `nPolyObsVtx3DPending(mObs, mElt)` (both in
    elt_mod), defer the projection.
- xObs parser, after setting xObs_FLG, walks the pending-list for the
  current element and runs `SetCvxPolyObsVtx + SetCvxPolyObsBound` on any
  stashed vertices.
- ChkDf2 (iosub.inc): after the existing psiElt-based xObs default kicks in
  for elements that need it, runs the same pending-resolution pass. If
  xObs is still unavailable for an element with pending vertices, ELEM_OK
  is set false with a targeted message naming iElt and iObs.
- Error message text now references both accepted keyword names —
  `PolyObsVec/Poly3DObsVec` — not the never-implemented `PolyObs3DVec`.
- Note: `PolyApVec=`/`Poly3DApVec=` have the same xObs-ordering wrinkle
  but were left as-is (still error on bad order). Same deferral pattern
  would extend to them if needed.

## MOD command: empty-value crash guard (macosio.F MOD_LOH)
- `mod ngridpts = 33` (spaces around `=`) used to crash with
  `forrtl: severe (24): end-of-file during read, unit -5`. CACCEPT in
  MOD_LOH only consumes the first whitespace-separated token, so VALUE
  ends up empty and the subsequent `READ(VALUE,*) nGridpts` (which has
  no ERR=/END= label) hits EOF on the empty string.
- Guard added right after GET_EQ: if `LEN_TRIM(VALUE) == 0` and the
  command isn't HELP, print a clear warning ("assignment must be
  \"<var>=<value>\" with no spaces around \"=\"") and drain remaining
  buffered tokens (`read_len(pstack)=0, read_cur(pstack)=0`) so the
  leftover `=` and `33` tokens don't cascade into more bogus commands.
- Smaller scope than fixing the parser to slurp multi-token assignments.
  Users get one targeted message and prompt continues. `mod ngridpts=33`
  (no spaces) still works as before.
- Also dropped the cosmetic `WRITE(*,*)` blank line in the QUIT branch
  of MOD_LOH — left a stray empty line between the MOD prompt and
  `MACOS>` on every exit.

## Prescription validator (validate_prescription.F90)
- Phase-1 pre-validator: `validate_prescription_mod%ValidatePrescription
  (filename, ios, msg)` runs before MBFile6 opens the .in file. Pure character
  I/O: catches missing values after `Key=` (inline or continuation), bad keys,
  file-not-found, and open errors. Also catches blank lines inside multi-row
  continuation blocks (e.g. a stray blank in the middle of a `Tout=` matrix),
  which trip msmacosio.inc's fixed-format READ even when every individual
  line is well-formed. Returns ios=0 on success; ios/=0 with a
  bare `msg` (no filename prefix — caller adds context).
- Per-key empty-value exceptions (KeyEq helper inside the module): `EltName`,
  `BaseUnits`, `WaveUnits` may have empty values. Older prescriptions leave
  these blank and the parser tolerates it. Add more keys to the `.OR.` chain
  in `ValidatePrescription` if Norbert flags additional ones.
- Wired into MBFile6 in macosio.F (interactive: re-prompt on failure via
  `GO TO 43`) and smacosio.F (SMACOS: clean abort with LOAD_SUCCESS=.FALSE.,
  `GO TO 99`). Both paths run the validator only when iopt=1 (existing file
  being loaded).
- Top-level `VALIDATE` command added in macos_cmd_loop.inc — same subroutine,
  no state change, just prints `<file>: OK` or `<file>: <msg>`.
- Phase 2 (deferred): IOSTAT-armoring of msmacosio.inc parser internals;
  enum-value validation against the type tables.

## CLI sub-prompts via mhist cache (macosio.F + mhist.c)
The old fragile pattern was: Fortran ACCEPT routine WRITE-s a prompt
(no trailing newline), then calls READ_LOH which goes through MHIST →
`readline(empty)` for the input. readline couldn't see the prompt and
mis-managed cursor state, producing the recurring "MACOS> the new
element data?" overwrite class of bugs.

Replaced with a small cache in mhist.c:
- `sub_prompt_buf[256]` (file-scope, OUTSIDE the READLINE_LIBRARY
  guard so non-readline builds still link)
- `set_sub_prompt_(s, slen)` / `clear_sub_prompt_()` — Fortran-callable
- `mhist_` for sub-prompts (mp[0]==' ') passes `sub_prompt_buf` to
  readline, so readline renders the prompt itself and owns cursor
  state authoritatively. Setter appends one trailing space to match
  the legacy `' ',A,' [',A,']: '` format the old WRITE produced.

CACCEPT, IACCEPT, DACCEPT, RACCEPT in macosio.F (NOT in smacosio.F —
that one has no readline) now:
1. Build the full prompt string from TEXT + `[default]:`.
2. If `read_len(pstack) > 0` (token will come from the existing
   buffer, no MHIST call), WRITE the prompt directly with format
   `(A,$)` — matches legacy behavior for multi-token lines like
   `stop elt 1 0,0`.
3. Else `CALL set_sub_prompt(...)` then `CALL READ_LOH(...)` and
   `CALL clear_sub_prompt` afterwards.

Use `LEN_TRIM(prompt_str)` for the slice length, NOT `ICLEN` —
`ICLEN` returns the length of the first *non-blank* prefix (length
of the first word), which gives 4 for ` Enter...` etc. That bug ate
several debug cycles; don't repeat.

When adding new ACCEPT-style prompt routines, follow the same
pattern. New top-level commands that need their own prompt should
either reuse the existing ACCEPT routines or push their prompt
through `set_sub_prompt` directly.

## Test prescriptions (ZGD_test_files/)
- All FreeForm test prescriptions use nGridpts=11 (gives 89 rays for Circular grid).
  Reason: model 256 has bRay=500 (BUILD limit) and mDrawRay=101 (nDrawElt array size);
  nGridpts=11 gives 89 rays which is under both limits.
  Grid-using prescriptions (tst_FF_g/mg/fg/refl) also need model 256 for mGridMat=256.
- Obscured rays now increment nBadRays so WARN fires: tracesub.F and propsub.F both
  have `nBadRays=nBadRays+1` inside the `IF (.NOT.L1(iRay))` / `IF (.NOT.LRayPass(iRay))`
  blocks that call SetRayFail(RayStat_Obscured).

## Regression tests via pymacos (~/dev/MACOS_resources/pymacos/)
- CodeV-comparison tests in `tests/test_api_rx_grating.py` and `tests/test_masks.py`
  cover the geometric / ray-trace paths (6601 tests, all passing).
- PROPER-comparison tests in `tests/proper_compare/` cover the physical-optics paths
  (INT/PIX/DFT-propagation) that CodeV can't validate. Organised by phase:
  - **Phase 1** (`test_cass_ff*.py`, results in `results_phase1/`): far-field
    image-plane PSF comparison on `Rx_Cass_FarField.in`. Nominal + SM Tx/Ty/Tz
    perturbations.  Peak-normalised agreement at 1e-11 with OPD pass-through.
  - **Phase 2** (`test_coro_nfprop.py`, results in `results_phase2/`):
    near-field plane-to-plane propagation between Elt 2 and Elt 3 of
    `Rx_Coro.in` (HCIT-style coronagraph, simplified-conic version).
    Sum-normalised agreement at 2.4e-14 RMS, 4.8e-13 max -- at
    double-precision FFT round-off (needs the runtime `dx_at()` query
    described below + odd `nGridpts=511`).
  - **Phase 3 + 4 + 5** (`test_coro_nfprop_phase3.py`, results in
    `results_phase3/`): six additional Coro steps:
      - 3a NFPlane 5→6 (sampling-limited 3.7e-8)
      - 3b sphere→plane 8→9 (4.2e-11)
      - 4a pupil→pupil 8→10 through focus (5.9e-13)
      - 4b NFPlane 13→14 (7.4e-5; localized post-focus diffraction ring)
      - 5.1 ExitPupil 20→FocalPlane 21, no mask (2.3e-9)
      - 5.2 same with FPM=400 um + Lyot=14 mm coronagraph
        (2.6e-6 Strehl-norm against a peak suppressed by 3.2 million
        relative to the un-coronagraphed baseline).
  - **Phase 6a** (`test_coro_apodizer.py` + `apodizer.py` + new
    `pymacos.apodize` wrapper): pupil apodisation via a user-supplied
    NxN amplitude mask, multiplied into macos's WFElt and PROPER's
    wavefront via the same numpy array.  Apodiser construction uses
    KxK super-sampling for sub-pixel area-weighted aperture edges so
    low-N results converge to high-N as the same physical shape.
    First test (Gaussian-edge pupil at Elt 5, NFPlane to Elt 6):
    4.0e-8 max, 1.1e-9 RMS (same sampling-limited regime as 3a,
    confirms the wrapper integrates correctly with downstream
    propagation).
  - **Contrast scoring** (`test_coro_contrast_curve.py` +
    `contrast.py`): radial dark-zone contrast vs lambda/D at the
    science focal plane.  lambda/D derived empirically from the
    no-mask PSF's first Airy null (1.22 lambda/D with a fractional-
    depth guard).  Phase 5.2 dark-zone contrast hits ~3e-10 in the
    7-10 lambda/D range; bright outer-ring artefact at ~15 lambda/D
    (Lyot edge diffraction).  Baseline for Phase 6b/c apodised
    designs to score against.
- Run via `./run_proper_tests.sh` at the pymacos root.  It rebuilds pymacosf90
  (if needed) and invokes each phase in its own pytest process to dodge a
  pymacos state-leak across model_size transitions (512 ↔ 1024).  Artefacts
  go into per-phase `results_phaseN/` directories.
- The pymacos<->PROPER coupling needs three reconciliations to reach the
  observed agreement:
  - **Aperture match.** PROPER takes amplitude DIRECTLY from macos's mask
    via `prop_multiply`, NOT from an analytical circular_aperture +
    obscurations model.  Mismatched apertures put light where macos has
    zeroed it out, breaking the wavefront tilt and halving the apparent
    PSF shift under perturbation.
  - **Sign flip.** macos OPD sign convention is opposite to PROPER's
    `prop_add_phase` input.  Default `opd_sign_flip=True` reconciles them.
  - **Normalisation choice.** `compare_and_record` takes `norm_kind=`
    `'peak'` (Strehl-norm, default; right for image-plane PSFs) or `'sum'`
    (flux-norm; right for pupil-plane / near-field intensities).  Peak-norm
    inflates the residual on flat-top NF PSFs (where peak position inside
    a uniform region is noise-dominated); sum-norm gives the physically-
    meaningful flux-conservation precision.
- Centroid-based (not peak-based) alignment in `compare_and_record`:
  intensity-weighted center of mass is robust for both sharp Airy peaks
  (Phase 1) and flat-top NF pillars (Phase 2).

## macos_api_mod: language-neutral wrapper layer for pymacos + mmacos

`macos_f90/macos_api_mod.F90` (compiled into `libsmacos.a`) is the
shared SMACOS-call wrapper that both `pymacos` (Python) and `mmacos`
(MATLAB) bind into.  Lives in `MODULE macos_api_mod`; holds the SMACOS
call-line buffer (`command`, `CARG/DARG/IARG/LARG/RARG`), the shared
output arrays (`PixArray`, `OPDMat`, `RaySpot`, `OPD`, `SPOT`, `PIX`,
`USER`), the package-level state flags (`firstEntry`, `rxLoaded`), and
~96 wrapper subroutines covering init, load_rx, save_rx, modified_rx,
n_elt, trace_rays, opd_val, spot_cmd/get, int_cmd/get, cfield_cmd/get,
cfield_apodize, elt_dx_get, base_unit_to_metres, perturb_elt,
prb_elt(_grp), the elt_* Rx-inspection/modification family, the
source-side set_src_*/get_src_*, the elt_grp_* family, stop_info_*,
xp_get/set/fnd, sxp_fnd, ors_run, srs_run, plus internal helpers
(`SystemCheck`, `EltRangeChk`, `StatusChk1`, `translateSurfaceID`,
`translateEltID`, `checkSurfaceID`, `checkEltID`).

`ray_status_get(OK, RayStatus_, RayFailElt_, nRays)` (added 2026-06-12,
Q2) exposes `elt_mod`'s per-ray `RayStatus(:)` / `RayFailElt(:)`
(`RayStat_OK/Obscured/Miss/Bracket/MaxIter/Undef`) verbatim from the
last trace — the per-category complement to `ray_info_get`'s binary
LRayOK/LRayPass.  mmacos veneer: `get_ray_status(N)`.  Backs the design
layer's ray-loss guard (`PLAN_DESIGN_LAYER.md` §1.3.4).

The MATLAB **design layer** (`macos.design.System`,
`MACOS_resources/mmacos/src/+macos/+design/`) is the first heavy consumer of
this wrapper surface — import (`from_rx`) → `sensitivities` (harvests the
Phase 7 `dw_dx` channels) → `vary`/`evaluate`/`optimize`.  It adds no
Fortran; it's a pure MATLAB layer over these routines.  Note: SLSQP
objects live in `libslsqplib.a`, a separate archive from `libsmacos.a`
— any binding that links libsmacos must also link slsqplib (the SLSQP
back end references `slsqp_`).

Promotion history: until §5.2 the file lived in
`MACOS_resources/pymacos/src/cmake/source/macos_api_mod.F90` and was
compiled separately by each binding's own CMake/Make build.  Now it
lives in `macos_f90/` and the binding builds just link `libsmacos.a`
+ pull the `.mod` from `mod_smacos/`.  pymacos's `CMakeLists.txt`
takes `MACOS_BUILD_DIR` as a cache var (default
`~/dev/macos/build_release`); mmacos's Makefile takes the same.

gfortran cleanliness: the file is portable across ifx and gfortran
as of §5.2.  Notable porting fixes — replaced `(/'m','cm',...'/)`
implicit-length char arrays with explicit `character(len=4)` forms;
converted all `do concurrent` to plain `DO` loops (some inside
`BLOCK` for index-decl scoping); 50× `LOGICAL == PASS/FAIL` sites
swapped to `.eqv.`; 16× INTEGER-comparison sites switched to integer
equality (`var == 0` / `var /= 0`) where the variable was INTEGER;
LOGICAL-array `ifGlobal` tests rewritten in `prb_elt_grp`; bare
LOGICAL `filter` test in `set_src_csys`.  pymacos pytest 6601/6601
+ PROPER-compare green under both compilers; mmacos 11/11 under
both.


## pymacos intensity() / complex_field() / dx_at() / apodize() wrappers
> These span `macos_api_mod.F90` + `pymacos_f2py.f90` + `macos.py`.
> The Python-facing detail already lives in
> `../MACOS_resources/pymacos/CLAUDE.md` — keep it there (single source
> of truth) and reference it; do not duplicate. The buffer-wrapper
> *template* (spot_cmd/spot_get pattern) belongs with macos_api_mod above.

## Coronagraph mask element type (silent failure mode)
Circular obscurations declared in a prescription (`nObs`, `ObsType=Circle`,
`ObsVec=...`) are **only applied to the diffraction-grid wavefront when the
element is `Element= Obscuring`**.  If the element is `Element= Reference`
(or any non-Obscuring type), macos still parses the obscuration metadata
and applies it to geometric rays during ray tracing, but the
diffraction-grid `WFElt` array sails through untouched -- a hard-edge
"mask" that's invisible to the diffraction propagation.

Caught while building the Phase 5 coronagraph test
(`tests/proper_compare/test_coro_nfprop_phase3.py`).  Original
`Rx_Coro_FPM.in` had Elt 9 (the FPM) as `Element= Reference` with a
132 um circular obscuration -- it looked like a working coronagraph
in the prescription, but `complex_field(9)` and `intensity(9)` were
byte-identical with and without the mask.  Suppression at the science
focal plane was ~17%, all of which came from flux trimming at the
Lyot stop (Elt 14) -- not the FPM.  After changing to
`Element= Obscuring`, the FPM finally bites: on-axis suppression at
Elt 21 became factor ~1e5 from FPM alone, factor ~3e6 with the Lyot.

The right model to follow is `MACOS_resources/docs/macos-manual/examples/CoroExample.in`
Elt 6 (the working CoroMask).  Any new coronagraph test prescription
should use `Element= Obscuring` for the mask element.

## IACCEPT_S reprompt-loop sweep (macos_cmd_loop.inc, 2026-07-16)
Interactive commands that validated an element/ray id and jumped BACK
to their ACCEPT prompt on failure (`GO TO <prompt label>`) hang or
misbehave under SMACOS batch: IACCEPT_S reads the exhausted stack and
returns the DEFAULT, so a sticky default (`IACCEPT_S(iElt, iElt, ...)`
— the OPD case) spins forever, and a fixed default (1, iCurWFElt,
nElt−1) silently substitutes the wrong element.  ALL such validation
reprompts now abort with `'<CMD>: command aborted.'` + `GO TO 1`:
BUILD, ORS, SRS (both prompts), FEXIT, SXP, XPS, SPOT, MVAR, DVAR,
GPERTURB, LPERTURB, and AVAR — whose old SMACOS branch called `stop`,
killing the host process when the engine is a mex/library.
Deliberately KEPT as loops: CTRACE / SEGRAYTRACE "0=quit" multi-entry
prompts (default 0 self-terminates on an exhausted stack), the
`#ifdef MACOS_CMD` MRESET size prompt (interactive-only), and the
Zernike-range prompt (sets valid defaults before looping, so it
self-heals).  When adding a NEW command prompt, never re-enter the
ACCEPT on validation failure — abort to label 1.

## LOG command cleanup (macos_cmd_loop.inc)
- The `LOG` interactive command (log10-intensity wavefront display, not
  transcript logging -- which is `JOU`/`JOURNAL`) previously dumped an
  always-on debug print (`*****B4 SrfOut: mdttl=`) plus a 5 MB ASCII file
  (`IntLog.txt`) to the cwd on every invocation, gated by an `#if 1`.
  Both removed; diagnostic block gated `#if 0`.  `LOG` syntax and 3-char
  min-match unchanged -- existing `.jou` scripts using `log <iElt>` keep
  working.

## Polarization physics (PLAN_POLARIZATION, Phase 0/1)
See `PLAN_POLARIZATION.md` (root) + `POLARIZATION_PHASE0_AUDIT.md` for the
full audit.  Quick map of what the engine has and the conventions to hold.

**What exists (all gated on `ifPol`):**
- **Per-ray vector E-field** `RayE(3,mRay)` (`Complex*16`, `elt_mod.F`) --
  the polarization state carrier, s/p-decomposed at each surface in
  `Reflector`/`Refractor` (`elemsub.F`).  Source state is `Ex0/Ey0`
  (`src_mod`), mapped onto each ray's local x̂/ŷ in `sourcsub.F`/`ssrcray.inc`.
- **Fresnel + recursive multilayer thin-film coatings** with complex index
  in `Reflector`/`Refractor` (`elemsub.F:432-547`).  `Reflector`'s TP/TS
  transmittance sub-blocks are `if(.false.)` dead code (mirrors return R only).
- **Vector diffraction** = the 3 Cartesian components Ex/Ey/Ez propagated as
  independent scalar fields, gated `ifVecDif3`.  Since **Phase 3a Tranche 1**
  this covers EVERY leg (see the Phase 3a section below); before it, only the
  far-field sphere→plane leg (PropType 3, via `PFFPROP`) was vectorized.
- CLI: `POLARIZATION`/`NOPOLARIZATION` (sets `ifPol`; enables `ifVecDif3`
  when `mWF≥3` -- all stock model sizes have `mWF=3`), `VECTOR`/`SCALAR`.
  `POLARIZATION` is SMACOS-dispatchable (LoadStack packs `Ex0/Ey0` as
  `DARG(1..4)`, `smacosutil.F:160`).

**Stubs / gaps (do not assume these work):**
- `RfPolarizerElt(14)`/`TrPolarizerElt(15)` are **name-table-only** -- no
  trace dispatch anywhere.  `JmatElt(2,2,mElt)` is allocated/zeroed and
  otherwise dead.  No Jones/Stokes/Mueller math, no waveplate, no VVC.
- `srtrace.F`'s `ifPol=.false.` is a **local `parameter` in the dead
  `SRTRACE_Test` driver only** (its caller is under `#if 0`).  The production
  single-ray paths (`SRTrace`/`CRTrace`, tracesub.F) thread `ifPol` as an
  argument -- so there is nothing to "lift" here.

**Conventions (pinned in Phase 0 -- assert in tests, do not relegislate):**
| Convention | Value |
|---|---|
| Time-harmonic | `exp(+iωt)`; spatial propagator `exp(−ikz)`.  Derived from `elemsub.F:387` (`C1=exp(−i·2π·L·N/λ)`, phase decreases as `L` grows) + the coating recursion `elemsub.F:512-516`, consistent with the 2026-07-25 IFO finding (field phase advances as OPL shortens) and the pymacos↔PROPER `opd_sign_flip=True`. |
| Absorbing index | `N = n − iκ`, κ>0 = loss (as stored in `IndRefArr`/`ExtincArr`, applied `DCMPLX(n,−κ)`). |
| Jones storage basis | Linear (x,y); circular via a unitary change of basis (Phase 2 decision). |
| Coating thickness | Rx `Coating=` layer thickness is **waves at parse-time `Wavelen`**, converted to physical (`·Wavelen/IndRef`) at load (`msmacosio.inc:2660`); the trace applies phase at the *current* λ so broadband sweeps are already correct.  `Coating=` must follow `IndRef=` (boundary media snapshot).  **`coat_set` takes PHYSICAL thickness** and sidesteps all of this. |

**Two coating subsystems** (do not conflate): Model A = `Coating`/`EltCoat`/
`IndRefArr`/`ExtincArr`/`EltCoatThk` (complex index, drives the polarization
path, `coat_set`/`coat_get` target this).  Model B = `nCoatElt`/`CoatIndxElt`/
`CoatThkElt` + `AirGap` (real index, non-sequential refractive path only).
Extend A, leave B (Dave 2026-07-25).

**Phase-1 API surface** (`macos_api_mod.F90`, `beam_set`/`beam_get` template):
`pol_set`/`pol_get`, `vecdif_set`, `coat_set`/`coat_get`, `rayfield_get`.
`Ex0Ey0=` Rx keyword now parses (4 reals `ExRe ExIm EyRe EyIm`) + round-trips
through SAVE (STATE only; on/off is API/CLI).  See `SAVE_KEYWORD_AUDIT.md`.
`pol_set`/`vecdif_set` call `modified_rx` (Phase-2 fix): the POLARIZATION
command changes trace-relevant state (`Ex0/Ey0` seed `RayE` at source-grid
setup) but does NOT reset the cached trace -- without the dirty, a pol-state
change + re-trace harvests the PREVIOUS state's `RayE` (verified stale).

**Coated-branch fixes (Phase 2, this box).**  Two latent bugs in the
`ifPol`+`nCoat/=0` recursion, both in `Reflector` AND `Refractor`
(`elemsub.F`), both invisible to intensity-level tests and fatal to Jones
work; found by the Fresnel-analytic gate (`tJonesPupil`), which now pins
them at 1e-12:
1. **Incident medium**: the recursion read `nb_arr(0)` = the parser's
   `IndRefArr(0,iElt) = IndRef(iElt-1)` -- for a coated mirror FOLLOWING
   another mirror that slot holds the previous mirror's conductor-idiom
   substrate (`Extinc=1e22`), i.e. light modeled as arriving from inside a
   perfect conductor.  Fixed: the coated branch now uses `na,kxa`
   (= `CurIndRef/CurExtinc`, the medium the ray actually travels in --
   what the uncoated branch always used).  The stored slot 0 is now unused
   by the trace; the parser convention is unchanged.
2. **Signed incident cosine**: `ccfb_arr(0)=DDOTC(ihat,Nhat)` is NEGATIVE
   when the normal faces the beam (psi convention), which turns every
   interface coefficient into its RECIPROCAL (1/r): |R|>1, s/p roles
   swapped, retardance sign flipped -- while |D| survives (the sneaky
   part).  Fixed with `DABS` (mirrors the uncoated branch's `ccfa`).
   Diagnostic signature if it ever returns: measured/analytic RS/RP ratio
   = (RP/RS)^2 exactly.

**Jones input basis (engine launch frames, `ssrcray.inc`).**  Collimated
sources launch every ray with `E = S*(Ex0*xGrid + Ey0*yGrid)` -- the
source-frame pair, UNIFORM over the grid.  Point sources use a PER-RAY
frame: `yray = unit(RayDir x xGrid)`, `xray = yray x RayDir` (reduces to
the global frame on the chief).  `S` = flux normalization (~1/sqrt(nRays))
-- a common real scalar carried by the Jones pupil.  ColSource re-
orthogonalizes the Rx frame as `z=+-Chf; y=unit(z x x); x=y x z`.  The
perfect-conductor mirror idiom (`IndRef=1, Extinc=1e22`) gives RP=RS=-1
+O(1e-22): polarization-neutral, which makes any stock conductor Rx a
unitarity gate for free.

**Phase-2 binding layer** (mmacos `+macos/jones_pupil.m` + `pol_maps.m`,
pymacos `macos.py` same names): two-trace Jones pupil (double-pole default
basis / local-sp / global) + closed-form 2x2 Pauli polar decomposition
(D/Dvec, ret/retvec in [0,pi] w/ near-pi ambiguity FLAG, T, phase; pupil
mean and variation reported SEPARATELY -- the mean absorbs the system's
geometric rotation and is a state change, not an aberration).  D/T are
singular-value invariants (basis-independent); ret is exactly what
double-pole makes artifact-free (local-sp inflates ret var ~10-250x --
asserted in tests as the documented artifact).  Gates: `tJonesPupil`
(mmacos, incl. the Bench fold Fresnel gate) + `test_jones_pupil.py`
(pymacos, ifx-linked = the standing ifx smoke).

**Phase 2b low-order expansion** (`pol_zernike`, 2026-07-26, binding-side
only -- no engine change).  Expands the Dvec/retvec maps onto ANSI
Zernikes so results are comparable with the polarization-aberration
literature term-by-term.  The engine-relevant fact: an on-axis
rotationally symmetric two-mirror system MUST reduce to pure
**polarization astigmatism** -- astig0 in Pauli s1, astig45 in s2, equal
magnitude, no circular component, no defocus, and a rho^2 radial law
whose on-axis extrapolation vanishes (measured: forbidden terms at 1e-15
of the astigmatic term, D(rho=0)/D(rho=1) ~ 1e-4).  This is now the
sharpest cheap check that the coating/s-p machinery and the pupil
reference frame are both right: a broken frame breaks the mode pattern,
not just the map.  Residual asymmetries are DISCRETIZATION, verified by
scaling (astig-pair mismatch 1.9e-7 at model 128 -> 5.8e-8 at 256; the
symmetry-breaking magnitude term is quadrafoil-X, aligned with the pixel
grid's own axes, while quadrafoil-Y stays at 1e-17).

## Phase 3a Tranche 1 -- vector propagation across the whole chain

All in `propsub.F` (`CPROPAGATE` + two new module helpers).  Vector mode
now closes a multi-leg chain; the coronagraph case (pupil→FPM→Lyot→focal)
works, the VVC acceptance test is unblocked.

**What changed.**
- **Every leg loops the SAME scalar kernel over the 3 component planes.**
  One `kWF1..kWF2` range is computed once before the leg dispatch
  (`1..3` under `ifVecDif3`, else `iWF..iWF`), so every `DO kWF=kWF1,kWF2`
  degenerates to the original single call in scalar mode -- those paths
  stay **bit-identical**.  Covered: `NFPROP` (PropType 2, 5), `PPPROP`
  (4, 6), `SFPROP` (7, 8), `SPH2PL`/`PL2SPH` (10, 11 -- NOT in the plan's
  list, added for coverage completeness), `FRPROP` (12), `NFPropDFT`
  (13, 14), `FFPropDFT` (15).  Leg dx/z bookkeeping stays OUTSIDE the loop;
  the `DWF` scratch is reusable across planes.
- **`PFFPROP` retired** (3a.1(1)).  It was a bare per-component FFT that
  omitted the Fresnel output factors `FFPROP` applies via `applyfac2`
  (`1/(i·λ·dz)·dx1²` + the output quadratic phase).  The far-field leg now
  runs `FFPROP` per plane.  The routine survives in the file, uncalled,
  with a deprecation banner.  **A/B (Rx_Cass_FarField, model 128):** vector
  total power was `8.9377e-01`, is now `1.8155e+06` == the scalar total
  exactly (Parseval).  Factor 2.03e6 in intensity, 1.43e3 in amplitude.
- **Assembly: seed once, then update** (3a.1(2)).  Both `ifPol` branches
  used to RELOAD the grid from `RayE` at every physical-leg assembly,
  erasing earlier legs' diffraction (the non-pol branch always MULTIPLIED
  by the incremental geometric phase).  A SAVE'd `LWFSeeded`, reset with
  `CumLStart` at trace start (`iStartElt==0`), now seeds on the first
  assembly and multiplies thereafter, advancing `CumLStart` identically in
  both branches.
- **Vignetting at the seed** (found while implementing, not in the plan).
  `RayE` carries NO aperture clipping -- the surface routines report it via
  `LRayTrans` and never zero `Evec` -- so a seed straight from `RayE`
  resurrected rays the ray-side zeroing had already extinguished.  Both
  polarized branches now gate the seed on `LRayPass(iRay).OR.(iRay.EQ.1)`
  (the `.OR.` mirrors the `iRay.GT.1` guard the zero sites carry for the
  chief ray).  Consequence worth knowing: **polarization-ON / vector-OFF is
  now BIT-IDENTICAL to polarization-OFF** -- it was wrong by 21% after one
  leg and 38% after two.
- **Masks are 3-plane.**  A 0/1 transmittance is a diagonal Jones `t·I`.
  `FFObscure` is called per plane under `ifVecDif3`, each call handed a
  FRESH copy of `xObs` (it re-orthogonalizes its `xGrid` argument IN PLACE,
  so reusing the mutated vector could move an edge pixel by a ULP between
  components); the scalar call is untouched.  The 13 ray-side clip sites
  and the 2 taper sites go through new module helpers `WFZeroPt` /
  `WFScalePt`, which are the original single-plane statement in scalar
  mode.

**PHASE CONVENTION AT THE SEED -- read before "fixing" it.**  The seed is
`WFElt(i,j,k)=RayE(k,iRay)` with **no** phase correction, and that is a
MEASURED result.  Reading `elemsub.F:395` alone (`S1=-TWOPI*L/lambda`,
`C1=exp(i*S1*Na)` applied to `Evec` at every surface) suggests `RayE`
advances phase OPPOSITE to `WFElt`, whose non-pol assembly uses
`exp(+i*TPL*Cphi)` and whose kernels carry `exp(-i*pi*lambda*dz*f**2)`
(paraxial forward `exp(+ikz)`).  Two independent experiments say `RayE` is
already in `WFElt`'s convention: (a) seeding `RayE*exp(i*c*TPL*CumRayL)`
and scanning `c` over `{-2,-1,0,1,2}` on a tilted Cassegrain moves the
far-field centroid by exactly `(c+1)x` the scalar tilt shift, so `c=0`
reproduces the scalar pupil phase; (b) on `Rx_VecChain.in` the `c=0` seed
reproduces the scalar intensity to 4.5e-16 at every leg for x-polarized,
45° and circular input.  **RESOLVED by the Fable review (2026-07-26):**
no bridge at the seed is CORRECT, but the reason is structural, not "RayE
happens to share WFElt's convention".  The Return-leg sign flip
(`IF (ifReturnElt(iElt)) RayL(iRay)=-RayL(iRay)`, applied AFTER the
surface call) makes the two phase bookkeepings diverge by construction:
`CumRayL` SUBTRACTS return legs (physical OPD bookkeeping, equalized at a
proper reference sphere) while `RayE`'s per-surface
`C1=exp(-i*2pi*L*N/lambda)` phases ADD them (`L` pre-negation).  So the
RayE-vs-CumRayL phase relation is **TRAIN-DEPENDENT**: on the
Return-terminated Cass FF the assembled vector EP field matches the
scalar field's convention (measured: circular-concentration same=0.994 /
conjugate=0.002; signed far-field tilt response equal at +0.94x), while
plain trace-to-detector `RayE` measures CONJUGATE to the OPD map (slope
-0.9995, corr -0.9995).  There is NO universal convention bridge to
apply -- correctness rests on the behavioral gates (exact on the
mask-type/normal-incidence class Tranche 1 claims), and **Tranche 2's
J_run must carry phases explicitly relative to the CumRayL bookkeeping**
rather than assuming either sign.  The tests are the authority; re-run
`tVecChain` before changing anything here.  Note the scalar-pol branch
cannot serve as the vector seed regardless: it rebuilds from `|RayE|` and
throws away the per-component phases.
Debug gotcha that cost this review an hour: obscured rays carry
`RayE=0`, and `ATAN2(0,0)=0` reads as "zero phase, zero Ez" -- gate ANY
per-ray probe on `LRayPass` before interpreting phases.

**Gate prescription.**  `mmacos/tests/Rx/Rx_VecChain.in` (copied to
`pymacos/tests/Rx/`) -- collimated on-axis source, flat normal-incidence
uncoated planes, TWO bracketed near-field legs (`NFPlane`→`Geometric`, the
`Rx_Coro.in` idiom; consecutive same-PropType elements MERGE into one leg,
which silently defeats a multi-leg test), aperture + central obscuration.
The ray E-field direction is a CONSTANT unit vector, so `E_k = e_k·u(x,y)`
and the vector sum `Σ|E_k|²` must equal the scalar intensity EXACTLY.  On a
real off-normal train it cannot: `Rx_Cass_FarField` carries `|Ez|/|Ex| ≈
8.8e-2` at the exit pupil and vector/scalar legitimately differ ~2.6e-3.
**Run the gate at 45° and circular input, not just x-pol** -- with an x-only
source all the energy sits in plane 1, which the OLD single-plane
propagator carried correctly, so an x-pol-only gate passes VACUOUSLY.

**Constraint (document, don't work around):** vector mode repurposes the
`mWF=3` planes as Ex/Ey/Ez of ONE wavefront -- no multi-WF / COMPOSE
concurrently.

**`cfield_plane_get` (2026-07-26) -- reaching the component planes.**
`cfield_get` only ever returned `iEltToiWF(iElt)`, so in vector mode the
per-component field was unreachable from the bindings: callers could see
the summed intensity and never how Ex/Ey/Ez made it.  New sibling
`cfield_plane_get(OK, RE, IM, N, iElt, iPlane)`: `iPlane=0` is
`cfield_get` exactly (bit-identical, and that is what the bindings pass by
default); `iPlane=1..3` are Ex/Ey/Ez and are REFUSED unless `ifVecDif3` --
in scalar mode plane k is an unrelated wavefront, not a component, and
handing it back would look plausible and be wrong.  Surfaced as
`macos.complex_field(srf,'plane',k)` (mex cmd gained an optional 4th arg)
and `pymacos.complex_field(srf, plane=k)`.

This is what CLOSED the Tranche-1 attribution, and it corrected it: the
vector/scalar difference on an off-normal train is NOT simply the
out-of-plane intensity (that guess is off by ~2x).  Two mechanisms --
(1) the scalar seed `|RayE|` puts ALL the power, including the
out-of-plane part, into ONE propagating plane while the vector run leaves
only f=0.9979 in Ex (a near-pure rescale, 1-corr = 4e-8), and (2) Ey/Ez
diffract to their own pattern.  `Iv ~ f*Is + Iy + Iz` takes 2.56e-3 down
to 2.90e-4.  Pinned by `tVecChain/test_vector_scalar_difference_decomposition`.

**Still open -- Tranche 2.** Between physical legs the grid field is
advanced by a scalar phase, exact only when the intervening elements are
non-polarizing (Obscuring / Reference / FocalPlane).  A COATED or
reflecting surface BETWEEN legs (IFO recomb→detector with folds under
`DO_NEARFIELD`; VVC layouts with an OAP between masks) needs the per-ray
running Jones `J_run(3,3,mRay)` design in PLAN_POLARIZATION §3a.3.

**Gates:** `tVecChain` (mmacos, SUITE_FAST) + `test_vec_chain.py`
(pymacos, ifx-linked).  Both are non-vacuous -- checked 2026-07-26 against
the pre-fix engine, which fails them at 0.21..0.38 relative error and
mis-states total power by 4-7%.

