# MACOS engine notes — Surfaces and traces: FreeForm, grids, Zernike dispatch, root picks, gratings, media, GRIN, the far-field cut

> **Read this file before touching:** `surfsub.F`, `elemsub.F`, `tracesub.F`, `propsub.F`, `srtrace.F`, `elt_mod.F`.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

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


## GRIN: only model 5 integrates; any other GradInd= LOSES its rays by name (2026-10-10)
*(Documented: Manual §4 Refractors (gradient-index paragraph).  Below: the engine record.)*
The Refractor's gradient-index branch (`elemsub.F`, `IsGrinRef`) implements ONLY
`GradInd= 5` (radial n(r) = n0 + n2 r^2 + n4 r^4, `modelCoef(1)` = the R-K step `dt`,
`getIndex_grin5` in utilsub.F).  Every other model fell into `Else ! other lens models --
to implement` and built `pr` / `rhat` from UNINITIALISED `pcrnt` / `Tcrnt`: NaN OPD
(`tst_save_keys.in` and its twin `tst_block_comment.in`, GradInd 6, SAVE-keys fixtures "not
meant to trace").  Now: `LROK=.FALSE.` + surfsub `LGrinUnsup`, so CTRACE / CPROPAGATE (all
four classification sites) record `RayStat_Undef` 'GRIN model not implemented', counted in
`nGrinUnsup` (reset per trace, WARN line), ONE note per run.  The fixture is unchanged (its
GRIN keys still round-trip through SAVE).  Gate `cli_tests/cmd/GRIN_unsupported.jou` (red
pre-fix: NaN).  Model 5 itself cannot be loaded or saved today -- PLAN_CONSOLIDATION §1.2
(the impossible `<=2 .AND. >=5` GradCoef branch, the bare `return`, the SAVE writer, the
absolute-position depth test).  The FwdRoot counter's `ABS(new-old) > 0` form stays: no NaN
deck remains in the corpus, but a NaN ray from any other path would otherwise count as a flip.

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
*(User-facing parts documented: Manual §4 Interpolated Surfaces (grid span, oversize grid declined, grid frame = Zernike frame), §4 "Editing .In-Files" (NSCount), §5 per-ray status.  Below: the engine record.)*
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
- **NSRefractor CLOSED (2026-10-11):** the `NSRefractor` ROUTINE (refractive non-sequential)
  had kept a LOCAL null grid frame for its `GridSrf` call, so a GridData figure there imaged as a
  piston.  The element's frame (of `imin`, the element hit) is now an argument, threaded from both
  call sites; ZGD twins `tst_nsref_{seq,ns}.in` trace and propagate bit-identically (pre-fix the
  NS twin read OPD exactly 0).  Gates `cli_tests/cmd/NSRefractor_grid.jou`,
  pymacos `test_nsrefractor_grid.py`.
- **GridSrfdx semantics (user-facing):** grid span = `(nGridMat−1)·GridSrfdx`
  base units, centered on `pData`; rays outside the span get `fh=0` (no
  figure).  256-grid across a 280-diameter segment → `GridSrfdx = 280/255
  ≈ 1.1`.  Too-small dx = figure only near segment center (looks piston-ish).
- **Grid pixel-index guard is INT-overflow safe (2026-10-10).** The four
  grid-term sites (`SFFSrf`/`FreeFormSrf`/`SGSrf`/`NGSrf`) reject a
  non-finite / `≥2e9` `xi`/`yj` in floating point BEFORE `IDFLOOR` → clean
  off-grid miss (`fh=0`), counted (`nGridIdxOvf`, api `grid_idx_ovf_get`,
  WARN).  Without it, a ray whose surface-solve bracket diverges (`SFFZPB`
  on a huge `L`) overflows `i0=IDFLOOR(xi); i1=i0+1` (INT_MAX→INT_MIN wrap),
  the `i1>nGridMat` test reads false, and `GridMat` is indexed out of bounds
  → SIGSEGV.  Finite in-range indices are bit-identical to the old path.
  **ARM-only, by hardware (TO, measured 2026-10-10):** the INT_MAX
  saturation is ARM's (`fcvtzs`).  On x86_64 an out-of-range double ->
  INT gives INT_MIN (`cvttsd2si` "integer indefinite"; `floor(1d10)` =
  -2147483648 under gfortran AND ifx), `i0+1` does not wrap, and the OLD
  `i0<1` test already caught it -- `tst_zrngr_overflow.in` traces
  bit-identically to its finite twin on the pre-fix Linux engine, release
  or debug.  So the red is Apple Silicon only (CCMac's crash_opd); a
  Linux gate of B can only be a no-host-kill + counted regression gate.
  (The gfortran DEBUG build then aborted at start-up in math_mod LZERO on
  a legacy `A(1)` dummy under -fcheck=bounds; fixed the same day -- see
  "gfortran debug build RUNS" below.)
- **gfortran debug build RUNS (2026-10-10).**  `source ./makems.sh debug
  gfortran` is the only build that bounds-checks every line, and it could
  not reach the prompt: F77 `A(1)` dummies ("any length") read as extent 1
  under `-fcheck=bounds`.  All 16 converted to assumed-size `A(*)` /
  `X(LD,*)` -- math_mod LZERO, DSMPROD, DNEGATE, DEQUATE (ChkDf1, every
  load), rdEQUATE, DUNITVEC, DMPRDB, DMPRDC, DDOTLP, CAXPY, DSVDC, CSCAL,
  DROT, CNRM2, CSWAP; sunsub RFOURN `DATA(1)` (DZERO/IZERO were already
  `A(*)`).  Same meaning: release A/B 406 decks, both compilers, 0
  differences.  ifx's `-check bounds` skips a last extent of 1, which is
  why only gfortran saw it.  **Rule: never declare a dummy `A(1)`; write
  `A(*)`.**  Running, it found two REAL defects on the core corpus:
  (1) CTRACE's NS candidate test read `SrfType(iminPrev)` with `iminPrev =
  0` (`.OR.` does not short-circuit) -- guarded, bit-identical by
  construction; (2) **`PtSource` placed an EVEN grid one row high** --
  `jwf = j+(mdttl-nGridPts)/2+1`, so a grid of `nGridPts == mdttl` wrote
  `RayID(.,mdttl+1)` (the nine `Rx_Coro*` decks, `nGridpts= 511`, which the
  load clamps to the model: the long-standing "Too many grid points.
  Resetting npts to" reset in BOTH macosio.F and smacosio.F, so an
  oversized Rx grid never wrote far out of bounds -- a parse-time refusal
  tried first broke 195 fast-suite tests and was reverted).  PR #70 had
  fixed `ColSource` only; now both use `(mdttl-nGridPts+mod(nGridPts,2))/2`
  (identical for odd grids and for `iwf`): a trace change for EVEN-grid
  POINT-SOURCE diffraction -- measured as a pure ONE-ROW translation of the
  seeded field (Rx_Coro, n=510: pre = post rolled by 1 in j; n=511 bit-identical),
  gate `pymacos/tests/test_ptsource_even_centring.py`.  (3) **a grid DECLARED larger than `mGridMat`**
  (`Rx_Coro_DM`'s `nGridMat= 1024`, its data installed by the bindings at
  run time) was traced as a phantom and indexed `GridMat` past
  `mGridMat`: the grid term now also requires `nGridMat <= mGridMat`
  (inert otherwise), a grid FILE that cannot fit is refused by name
  (`GridTooBig`, elt_mod: at the Rx `GridInit` sites and as a backstop
  inside it), and both api grid setters name their refusal.  (The "defined
  data" condition was NOT added to the trace: the api grid setters never
  set `ifGridDataDefined`, so it would switch off every api-installed
  grid; an in-bounds phantom reads the zeroed `GridMat` and is inert.)
  Debug core corpus at model 256: 114 of 114 ok, both IRIS fixtures ok.
- **`NSCount=` is the non-sequential HIT BUDGET and round-trips through SAVE
  (2026-10-10).** `NSCnt(iElt)` caps the surface hits a ray may take inside
  an `Element= NSReflector` group (tracesub.F ~4499): absent/0 = UNLIMITED
  (CornerCube's retroreflection), 1 = single-hit segmented groups (the right
  value for a grid-surfaced NS group, where an unbounded search diverges the
  bracket solver — the IRIS `save_rx` crash).  It is authoring intent, NOT
  structurally derivable; `PrtSingleEltInfo` (iosub.inc) now emits it when
  declared (`NSCnt>0`).  Was DROPPED on SAVE pre-2026-10-10.  Gates:
  `ZGD_test_files/tst_zrngr_roundtrip.in` (SAVE keeps NSCount) +
  `tst_zrngr_overflow.in` (tiny GridSrfdx → the overflow guard).

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

## Two Zernike refusals by name (PLAN_CONSOLIDATION 1.2 (tt) + (vv), Dave 2026-10-20)
- `NormAnnularNoll` (typeL 9) parses (ratio into `ZernAnnuRatio`) and dispatches to
  `ZerntoMon7`, but that converter never reads the ratio (`azdef.inc` is a fixed table):
  the surface traces as circular NormNoll at any ratio, modes 46-66 apply nothing
  (macos-fd's probe for the manual's Appendix E).  REFUSED wherever the type is set:
  the five surface-type parse sites in msmacosio.inc (ZernType / MonZernType /
  FFZernType; `ZernTypeUnsupported` in elt_mod.F), the three MOD sites, and the api
  `elt_srf_zrn_type` / `elt_srf_zrn_set`.  The SAVE writer's `ZernAnnularRatio` line
  (iosub.inc) is now unreachable and left in place for the day the converter reads r.
  The header key `ZCOZernType= NormAnnular` (Zernike OUTPUT fitting, zern_wf.F) is NOT
  refused -- a separate path with the same converter, not yet ruled on.
- `ExtFringe` indices 47-49 and 58-64 (radial order 11-14) are accepted by the parser,
  but `ZerntoMon_ZFE` maps Extended Fringe onto the 66-mode ANSI set, which ends at
  order 10: those ten coefficients applied NO surface.  A non-zero coefficient there is
  refused by mode: ChkDf2 on the three coefficient arrays (`ZernCoefInert`), the three
  MOD coefficient sites (`ZernModeInert`), the api coefficient setters (`ZernModeInertV`).
  Zero coefficients there still load (a dense 66-value deck is the common case).
- (xx) the same day: a coefficient PAST THE TYPE'S MODE COUNT (NormHex 46+, Fringe 38+) is
  refused by mode at every face -- ChkDf2 on the three arrays (`ZernCoefPastMax`), the MOD
  coefficient sites (index bounded by the array first: the old 'ToDo' no-op let the READ
  write past it), and `FFZernModes` / `MonZernModes` bounded by their arrays at load (the
  remap `Coef(Modes(:n),iElt) = ...` wrote past the array with no visible symptom in the
  release build -- the pre-fix binary loads `ZERN_mode_count.jou`'s mode-500 deck silently).  The api setters already refused past the count.
- Gate `cli_tests/cmd/ZERN_refused.jou` (red on the pre-fix binary).  The three SAVE-keys
  fixtures (tst_save_keys, tst_block_comment, pymacos Rx_SaveKeys) moved from
  `NormAnnularNoll 0.3` to `NormNoll` -- identical trace by the measurement above.

## NormHex opened (PLAN_CONSOLIDATION 1.2 (uu), Dave 2026-10-20: repair)
- `ZernType_MaxMode(NormHex)` was 0 ('implemented ???'): every `nZernCoef=` on a NormHex
  surface was refused at load ("outside 1 .. 0"), though `ZerntoMon4` (surfsub.F, table
  `hzdef.inc`) is called by all three dispatch chains for Zern / FF / Mon.  Now 45.
- VALIDATED FIRST (`cli_tests/tools/hzdef_check.py`, textbook, never the engine): the 45
  `hzdef.inc` polynomials are the Noll circle Zernikes Gram-Schmidt-orthonormalized over the
  unit regular hexagon, area-mean inner product (Mahajan & Dai, JOSA A 24, 2994 (2007));
  worst mode 38, 8e-12 of its largest coefficient; Mahajan & Dai's closed forms of modes
  4-11 to 1e-14.  ORIENTATION: vertices on +-xMon (the vertex-on-yMon set differs in 29
  modes).  `hzi(i,j,1)` is the x exponent (MonCoef's first index).
- PRECISION BUG found by the probe, fixed: the table's 313 literals carried no exponent
  letter, so they compiled as REAL(4) -- traced sag 1.8e-6 relative off the construction at
  mode 38 (the float32-rounded table predicts 1.5e-6, median 1e-7 both).  All now `d0`:
  traced sag = construction to 4e-12.  TRAP met doing it: line 480 had a BLANK inside a
  literal (`38.7129633223604929 8068692574968100787442`) -- fixed form ignores blanks, so it
  compiled right, but a `d0` inserted before the blank made it an exponent (5 m sag error
  at mode 28).  The check now parses as the compiler does (blanks dropped, d exponents
  required, finite values asserted).  `fzdef.inc` / `azdef.inc` have no long literals.
- The api setters take NormHex modes 1-45 only (`elt_srf_zrn_coef`, `_mon_zrn_coef`,
  `_ff_zrn_coef`, `elt_srf_zrn_set`, `_mode_set`, `validate_Zrn_Data`), the twin of
  Fringe's 37.  NOT closed (a finding, not a new refusal): a Mon/FF Zernike coefficient
  past the type's MaxMode (`nMonZernCoef` / `nFFZernCoef` are bounded by mMonCoef /
  mFFCoef only) and a MOD `ZernCoef(k,i)` past it are accepted and converted to nothing,
  Fringe (38+) and NormHex (46+) alike.  MOD MEASURED on tst_zern_normhex.in:
  `ZernCoef(46,1)=1e-4` is taken with no message, and with mode 10 zeroed the SRT ray at
  (0.5,0) meets z = 0 exactly (macosio.F's `ToDo --> FIX` range check).  Mon/FF read from
  the code (msmacosio.inc `nMonZernCoef` <= mMonCoef), not measured.
- Gates: `cli_tests/cmd/ZERN_normhex.jou` (SRT sag at (0.5,0) = Mahajan & Dai mode 10 to
  every printed digit; SAVE writes `ZernType= NormHex` and the reload traces the same;
  mode 46 refused sparse and dense; red on the pre-fix binary), pymacos
  `test_zernike_orderings.py` (all 45 modes vs the construction to 1e-9, the vertex-on-y
  must-fail leg, Table E.2 to its six figures).

## Zernike-apply trace dispatch requires ZernTypeL ≠ 0
*(User rule documented: Manual §4 Zernike Surfaces (ZernType default ANSI; the names; ExtFringe).  Below: the dispatch record.)*
- `propsub.F:230-249` (and the parallel `srtrace.F:145-153`) dispatches
  `ZerntoMon1/2/3/4/6/7` based on `ZernTypeL(iElt)`. The chain has NO `ELSE`
  branch — if ZernTypeL is 0 (the elt_mod.F:944 initial value), the IF/ELSEIF
  falls through, ZerntoMon is never called, MonCoef stays zero, and the
  surface intersection (which uses MonCoef, not ZernCoef) sees no perturbation.
  Result: silent zero-response.
- ChkDf2 defaults ZernTypeL to 1 (ANSI) for elements whose parse path visited
  the ZernType default block. **NSReflector+Conic elements don't visit that
  block** — their ZernTypeL stays 0 unless explicitly set in the Rx.
  **And SAVE then wrote an unreloadable deck (fixed 2026-10-10):** the
  writer's `Srf_Zernike` emitted a BLANK `ZernType=` plus `nZernCoef= 1`
  for such an element (sized by `ZernType_MaxMode(0)`, out of the array),
  and the parser's `nZernCoef=` range check read the same index 0 and
  `GO TO 99` with NO message -- a silent refusal of the reload
  (`tst_zrngr_roundtrip.in`; every binary back to 8cf28c3).  Now: no
  type declared, no Zernike block written (`lMon` still is; the Mon/FF
  writers likewise).
  **PLAN item 7 (Dave, 2026-10-10): EVERY Zernike-carrying surface
  (`Zernike`, `ZrnGrData`) defaults to `ZernType= ANSI`** when it declares
  none, whatever the element type -- at the end of ChkDf2, and at
  `nZernCoef=` in the parser (which prints `Default used for ZernType(  N)`,
  ChkDf2's 112 form).  MEASURED: a plain `Zernike` surface on an NSReflector
  ALREADY got the old default; the real gap was `ZrnGrData`, whose default
  is pre-empted at `Surface=` (the FLG is set for any surface but
  `SrfType_Zernike`), so its `ZernCoef` NEVER reached the trace.  The
  `nZernCoef=`-without-a-type refusal now applies only to a surface that
  carries no Zernike terms (`tst_nzern_conic.in`); `tst_nzern_no_type.in`
  (ZrnGrData) loads with the default.  Gate `cli_tests/cmd/
  ZernType_default.jou` on `tst_zerndef_{none,ansi,fringe}.in` (a
  ZrnGrData NS primary, no grid, mode 9 = 2e-7): none == ansi to the printed
  digit, fringe differs; PRE-change the none twin traced as the bare conic
  (its coefficient silently dropped).  Corpus: one mover, `xx.in`
  (load_fail -> ok) -- an OLD corrupted SAVE whose `ZernType=` is garbage
  bytes (`ZernTypeNameL(0)` read out of bounds by the pre-fix writer); it
  now takes the default with zero coefficients.  No corpus SAVE changed.
- GMI workaround at `GMI.F:1066`: when forcing `SrfType(iElt)=8` to apply a
  Zernike perturbation to an arbitrary element, also force `ZernTypeL=1` if
  currently 0. Without this, GMI's pzern channel silently no-ops when applied
  to a Conic-typed element. Caught by `test03_zern_response_e2e_pie` (had been
  masked under ifx because the same broken path also corrupted memory and
  triggered an exit-time SIGSEGV — gfortran exposed the actual zero-response
  bug cleanly).
- Noll IS handled now (verified 2026-09-14): all three trace chains dispatch
  `ZernTypeL=10` (Noll) / `ZernType_NormNoll` to `ZerntoMon6` —
  `propsub.F:267-269`, `tracesub.F:3413-3415`, `srtrace.F:154-156` (and the
  parallel FF/Mon blocks). `NormAnnularNoll` (9) → `ZerntoMon7`. So a user Rx
  with `ZernType= Noll` traces correctly; the earlier "silent no-op" note was
  stale (and inconsistent with the same section's own `ZerntoMon1/2/3/4/6/7`
  dispatch list).
- ExtFringe (`ZernTypeL=11`) CLOSED 2026-10-12: the converter `ZerntoMon_ZFE`
  (surfsub) had existed all along -- no dispatch chain called it, so an ExtFringe
  surface traced as its bare conic (one 'no ZerntoMon implementation' line per
  trace).  Now wired into all nine chains (Zern / FF / Mon x CTRACE / CPROPAGATE /
  SRTrace).  Its table was checked against the textbook rule (group p = (n+|m|)/2,
  |m| descending, cos before sin) -- 56 of 56 entries -- and equals Fringe's
  `ZerntoMon3` on modes 1-36; ExtFringe 37 is (6,6), Fringe 37 the (12,0) of CODE V.
  Inputs stop at 66 (`mZernModes`).  No corpus deck declares it.  Gates
  `cli_tests/cmd/ZernType_extfringe.jou`, pymacos `test_zern_extfringe.py`.

## Conforming Reference: Surface=Zernike/Aspheric is PASSIVE (sls-dev c9fa767)
*(Documented: Manual §4 Reference Surfaces (passive basis carrier; grid Reference the active exception; masks must be Obscuring) and Sparse Zernike input (wrapped modes line).  Below: the engine record.)*
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

## Physical-optics kernels are MEDIUM-AWARE (2026-09-30, dyson5)
*(Documented: Manual §6 "Setting up Diffraction Propagations" (media paragraph).  Below: the mechanism and the gate.)*
Every kernel call in `propsub.F`'s leg dispatch (NFPROP, PPPROP, SFPROP,
FRPROP, NFPropDFT, FFPropDFT, SPH2PL/PL2SPH, FFPROP, `FnCalc`) was handed
`WaveBU`, the VACUUM wavelength in base units, whatever medium the leg ran
in.  The inter-leg geometric phase was always right (`CumRayL` accumulates
`CurIndRef*RayL`, an OPTICAL path, against `TPL=2pi/WaveBU`), but a leg
INSIDE a medium ran at the wrong Fresnel number by n -- invisible to every
PROPER comparison because those are mirror trains in vacuum; fatal for a
Dyson, where the slit and FPA sit on the silica block's face.  Now:
`LegIndRef` = `CurIndRef` captured when `nEstart` advances (the medium the
NEXT leg traverses; `StartIndRef` for the first leg), and the dispatch
computes `WaveMed = WaveBU/LegIndRef` once per leg for all 21 kernel /
FnCalc sites (a `WRITE` announces any leg with n /= 1).  A leg is assumed
HOMOGENEOUS -- a leg that crosses a refracting face was never something the
kernels modelled.  Gate = an identity: a leg of z in index n equals the same
leg of z/n in vacuum (`pymacos/tests/test_prop_medium.py` on
`Rx_PropMedium_{glass,vac}.in`, twins of `Rx_VecChain.in`; intensities
compared -- the fields differ by a global piston).  Measured: pre-fix the
glass twin matched the UNSCALED vacuum deck to 1.1e-15 and missed the scaled
one by 74%; post-fix glass == scaled to ~1e-13 and != unscaled by 71-73%.
Vector mode inherits it (same `WaveMed` per component plane).

## Gratings on CURVED surfaces are CHORD-ruled (2026-09-30, dyson5 finding #2)
*(Documented: Manual §4 Reflective Gratings (the chord convention, OPD modulo λ, grating apertures).  Below: the engine record.)*
`Snells_Law_Grating` (elemsub.F, serves Grating / TrGrating / DoeTrGrating)
built the local grating vector from the UNITISED projection of the rule
direction onto the local tangent plane, i.e. a groove period constant ALONG
THE SURFACE.  A straight-ruled concave grating -- the element's own name, the
classical ruled grating, the CODE V / Zemax convention -- has equidistant
parallel groove PLANES (normal s0 = unit(h1HOE) in the vertex plane, spacing
RuleWidth) cutting the surface, so the local grating vector is
`(m lambda/d) * (s0 - (s0.N) N)` UN-normalised: |G| falls as sqrt(1-(s0.N)^2)
where the surface tilts.  The difference is a spectral blur proportional to
lambda and uniform over the slit -- TO's chain measured 2.8 px (Offner) /
3.3 px (Dyson) rms at 2500 nm against 0.003 / 0.04 px for the chord model --
that no concentric design can correct; it set every engine SRF in dyson5
beat 2.  FLAT gratings are unchanged (N = psi, the projection IS s0).  Fixed:
`shat = RuleDir - dot(RuleDir,Nhat)*Nhat`, no unitise.  Why nothing caught
it: the pymacos grating tests check parameter GETTERS only, and the one
concave fixture (`Grating_example_001.in`, Kr -6 m) is never ray-compared.
Gates: `tGratingImmersed` (mmacos; its closed form is now the chord model
and the pre-fix engine fails it by up to 0.5% of the kick on a 100 mm /
20 mm fixture) and `pymacos/tests/test_grating_chord.py` (the vector
grating equation per ray, chord vs surface model as the must-fail leg).
Holographic / variable-line-space gratings are a DIFFERENT model and are
not what `h1HOE`/`RuleWidth` describe.
**Finding #3 (2026-10-01, TO's propagation twin): the grating's OPL jump
had the same defect.**  `dL = -(na i - nb r).rho_prj` projected the hit
vector into the LOCAL tangent plane, = `(m lambda/d)(s0.rho_prj)`; the groove
count of equidistant planes is `(m lambda/d)(s0.rho)` with rho from the
VERTEX along the fixed ruling direction.  Difference `(m lambda/d)(s0.N)
(rho.N)` ~ rho^3/2R^2: cubic, zero on a flat grating; 12 waves on a 45 mm
footprint at R = 250 mm while the rays converged to 0.05 um -- rays and
path lengths disagreed.  Fixed: `dL = Order*lambda/RuleWidth*dot(s0,rho)`.
Gate `tGratingOpl` (TO; order 0 vs order -1 pupil OPD on a reference sphere
about the chief's focus).  **DAVE'S RULE: across a grating the OPD is
defined MODULO LAMBDA.**  The physical wavefront is the groove staircase;
the engine's smooth order-m phase function equals it mod lambda (it IS the
order-m Fourier term).  Compare OPDs across a grating in PHASE (wrapped, or
the complex field), never as unwrapped lengths; the ray-side scorer never
uses path lengths and is unaffected.

## ifLNsrf root pick is AXIAL, not radial; restarted traces reset PrevNonSeg (2026-10-03)
*(Documented: Manual §4 Conic Surfaces (the sheet rule, the RptElt pick on Return/Reference).  Below: the engine record.)*
TO's R2c Schwarzschild (`tests/Rx/Rx_SchwarzschildEP.in`): the EP
Reference (elt 1, the stop) sits 5.3 mm ahead of a CONVEX hyperboloid.
`macos.trace(5)` in ONE call gave a 93 mm rms spot; `trace(1)..trace(5)`
gave 10.6 um and matched the exact chain per ray.  Two defects:
1. **Root pick.**  The five base-conic routines (ConSrf / AsphSrf /
   AsphGridSrf / FreeFormSrf / UDSrf) choose between the quadratic's two
   roots; on a surface ON or right AFTER a Reference / Return (`ifLNsrf`,
   negative L allowed -- the FEX exit-pupil sphere traced BACK from the
   focus needs it) they used `|L^2 - mpr|` proximity with `mpr = |pin-pv|^2`,
   the ray's distance to the VERTEX including its LATERAL height.  A ray
   90 mm off axis has mpr ~ h^2 = 8.3e-3 and roots -0.107 / +0.011: the
   sheet BEHIND the ray "wins".  Now `LNsrfRoot` (module function) picks
   the root whose HIT POINT is nearest the element's REFERENCE point
   `prot` = `RptElt` (`d^2 = L^2 - 2Lq`, `q = ihat.(prot-pin)`): the vertex
   for an ordinary element, the POLE for an off-axis section.  M1 here ->
   +0.011; the FEX sphere about the focus (roots +-R, vertex behind) ->
   -R, a strict choice where the old metric relied on an exact TIE
   (mpr = R^2); a 90-deg OAP takes the side its pole is on.  **A first cut
   used the axial distance to the VERTEX plane, `s = ihat.(pv-pin)`: right
   on both of those, and WRONG on tBench's 90-deg OAPs (the ray runs
   parallel to the parent's vertex plane, the metric is degenerate, the
   chief went 600 mm off the pole).**  Any future rule here must be checked
   on all three: a pupil Reference just ahead of a convex conic, the FEX
   sphere, an OAP after a Reference.  Measured A/B (pre-fix HEAD worktree
   vs post-fix, gfortran): FEX radius and OPD bit-identical on e5hex1 /
   Rx_Cass_FarField / jwst zoom (stop 25) / SegDemo3conic.  Forward trains
   never reach this branch (ifLNsrf false -> the positive-root branch,
   untouched).  `LNSFlowL` (NS probes) is a third branch, also untouched.
2. **Stale `PrevNonSeg`.**  CTRACE advances `PrevNonSeg` element by
   element and never reset it between rays; a trace RESTARTED at
   `iStartElt > 0` (what `OPD` does when a trace exists and
   `iCurRayElt <= iEndElt` -- the API's `trace(ie)` loop) left every ray
   after the first seeing the previous ray's LAST element as "previous",
   so `ifLNsrf` was FALSE for them and the positive root was forced: the
   stepwise answer was right BY ACCIDENT.  Reset per ray at the ray-loop
   head.  `trace(nElt)` and the stepwise loop are now bit-identical.
Gate: `tTraceRestart` (mmacos, SUITE_FAST; pre-fix 9.33e-2 m, CLI RMS OPD
8.87e-2 m).  Diagnosis pattern worth keeping: a pty-driven CLI
(`onecall_cli.py` pattern) bisected `opd 1; opd 5` (right) vs `opd 2; opd 5`
(wrong) in minutes, then one `WRITE` in AsphSrf per ray showed the two
paths picking the SAME root at M1 with DIFFERENT `ifLNsrf` -- the
ray-by-ray diff found the divergence at ray 2, not ray 1 (the chief is
on axis, h = 0, and never sees the metric's flaw).

## Forward root pick is the VERTEX SHEET, not the vertex distance (2026-10-04)
*(Documented: Manual §4 Conic Surfaces.  Below: the engine record and the corpus A/B.)*
The 10-03 entry above closed the `ifLNsrf` branch and left the FORWARD
branch alone ("forward trains never reach this").  They do, through the
other defect of the same metric.  The `.NOT.ifLNsrf` pick in ConSrf /
AsphSrf / AsphGridSrf / FreeFormSrf / UDSrf (+ IntSrf in didesub.F and
the RefSrf / ObsSrf / PolElt / ReturnSrf conics in elemsub.F) chose the
root by `|L^2 - mpr|`, `mpr = |pin - pv|^2` -- the ray's distance to the
VERTEX, lateral height included, with no notion of which SHEET of the
conic the hit is on.  dyson5 stage B (TO, addendum 38): on the 3k
eccentric TMA section at -3 deg / 180 mm, `Telescope.optimize` reported
mm-scale per-field WFE (8.5 / 6.8 / 3.9 mm) with every ray passing while
the trace's 4-sigma-cut spot said 400 um.  CALIB was FAITHFUL: the plain
`OPD` at the FP is 7.8 mm rms on that deck because M3 (K = -15.8,
Kr = -0.1617 m, sheets 2a = 21.9 mm apart) had 40 of 1185 rays on its FAR
sheet (dz = +22.0 mm = 2a, at h 7-8 mm from the vertex), 44 mm of extra
path and a 10 mm spot tail -- all "passing" (the far sheet is a perfectly
good intersection).  The bimodal OPD histogram (90 % within 1.5 mm, 10 %
at +42 mm) is the tell; an outlier-cut spot metric hides it completely.
**Rule (`FwdRoot`, surfsub.F):** the real surface is the sheet the vertex
is on.  With `p` the hit relative to the vertex and `z = psi.p`, the conic
is `|p|^2 + 2 Kr z + Kc z^2 = 0` and the normal's axial component
`Nvec.psi = Kr + (1+Kc) z` has the SIGN OF Kr on the vertex sheet and the
opposite sign on the other (they meet where the normal is perpendicular
to the axis).  A forward root (`L > 0`) on the vertex sheet is taken,
first crossing on a tie; with none, the legacy pick is kept to the bit.
**"First forward crossing" alone is WRONG, measured on the same deck:**
M2's far sheet is a bowl opening toward M1, 37 mm IN FRONT of the mirror,
and 125 rays cross it before the real surface -- min-positive put them on
it and lost 350 rays at M3 (RMS OPD 0.34 m).  So the sheet decides, not
the order.  Measured post-fix: elements 1-2 bit-identical, 40 flips at
M3, FP OPD 7.84e-3 -> 5.79e-4 m, CALIB initial WFE [8.47 0.24 3.89 6.77
..] -> [0.58 0.24 0.29 0.50 ..] mm.  **Corpus A/B (406 decks, CLI, pre
vs post, `opd nElt`): 398 load on both, (RMS, P-V, nPass, lost) identical
to 10 digits on every one; the only flip count on a legacy deck was
`tst_save_keys.in`, a NaN trace on both engines -- `NaN .NE. NaN` is
TRUE, so the counter now tests `ABS(new-old) > 0`.**  Every pick that
differs from legacy is counted in `nFwdRootFlip` (reset per trace with
`nZPFailMsg`; api `fwd_root_flips_get`; WARN line) and ONE note prints
per RUN (a per-trace note floods CALIB).  A nonzero count on a legacy
deck means that deck was tracing onto the wrong sheet.  Gate: `tFwdRoot`
(mmacos, SUITE_FAST, `Rx_TwoSheetTMA.in` = the stage-B section as
emitted): vertex-sheet residual < 1e-6 (pre-fix 2.2e-2), FP OPD < 1 mm,
no spot tail, M2 bit-identical, the Schwarzschild twin unmoved, CALIB
sees every field; the pre-fix mex fails 3 of 5 by exactly those numbers.
The OPD at a FP the chiefs hit at 31 deg is still 0.5 mm-class on a
400 um spot -- that is the tilted-plane OPL metric (rodgers1 A.1), the
reason `add_pupil` exists; it is now a REAL number, not a wrong one.

## A Grating never vignetted by its aperture (2026-10-06, dyson5 3k join)
*(Documented: Manual §4 Reflective Gratings (last sentence).  Below: the mechanism and the gate rule.)*
`elemsub.F`'s four grating routines (`Grating`, `TrGrating`, `FzpTrGrating_`,
`DoeTrGrating`) computed the aperture / obscuration verdict with
`ChkRayTrans` and then OVERWROTE `LRT` on the next statement with the
interpolated-surface check (`LRT = .not.(SrfType.EQ.Interpolated .AND.
IERROR.NE.0)`), so a ray outside a grating's `ApType= Circular` aperture
(or inside its obscurations) PASSED; only its E-field was zeroed.
`Reflector` / `Refractor` keep the verdict.  Found because the engine's own
render of the dyson5 3k end-to-end deck drew the grating smaller than the
beam: 251 of 1185 center-field rays hit it 154-213 mm from its axis against
a 154 mm aperture, every one "passing" -- the (c) freeform telescope's
marginal rays leave at up to F/1.2 into an F/1.8 spectrometer.  Fixed: the
two conditions are ANDed, `LRT = .TRUE.` when no aperture is declared.
Measured after: the 1.5k row moves by 0.01 px (its 183 mm sky beam overfills
the footprint+5 mm grating by 1 %, admits 0.988), the 3k row's smile 2.05 ->
1.64 px, CRF 4.22 -> 4.02, admits 0.953; every spectrometer-alone dyson5
record and TO's Offner rows are unchanged (their apertures are their own
footprints + 5 mm).  Gate `tGratingAperture` (mmacos, SUITE_FAST): pass
flags == the geometric inside-the-circle set ray by ray, the clipped rays
carry `RayStat_Obscured` and stay clipped to the FP, the no-aperture control
passes every ray; the pre-fix engine passes all 89.  Lesson, same as
2026-09-30's groove law: the grating tests checked parameter GETTERS, the
chord law, and the OPL jump -- nothing ever put a ray OUTSIDE a grating's
aperture.  Any element type that vignettes needs one ray outside the
aperture in its gate.

## Far-field evanescent cut, opt-in (2026-10-03, dyson5 addendum 15)
*(Documented: Manual §6 Far-Field (the evanescent paragraph); cmdref Part II `ffcut`.  Below: the mechanism and the gates.)*
A far-field leg (`FFPROP` / `FFPropDFT`) maps spatial frequency f to the
output coordinate `x = lambda*dz*f`, i.e. `sin(theta) = x/dz`.  The
output window is `lambda*dz/dx1` wide, so a pupil sampled finer than
`lambda/2` (a 20 um pinhole at 1 um on a 64-point grid; any heavily
zero-padded fine pupil) carries pixels with `|x| > dz` -- `|f| > 1/lambda`,
EVANESCENT, no propagating energy -- which an energy-fraction metric over
the window then counts as light (0.64 % on `Rx_FarFieldPinhole.in`).
`dft_mod%ifFFEvanCut` (default OFF -- the record's propagating
normalisation stays the convention) makes both kernels zero every output
pixel with `x^2 + y^2 > dz^2` right after `applyfac2` (`FFEvanCut`, same
index convention: pixel `n/2+1` is x = 0) and print ONE line per call
with the count (`nFFEvanCutPix`).  API `ffcut_set(on)` / `ffcut_get(on,
nPix)` -- session state, not reset by a load, dirties the cached
propagation; mmacos `macos.ffcut`, pymacos `pymacos.ffcut`.  Gates
`tFFCut` (mmacos, SUITE_FAST) + `test_ffcut.py` (pymacos): the pinhole
(pixels zeroed, the inside bit-identical, the total lower) AND
`Rx_Cass_FarField` as the must-not-change twin (dx1 >> lambda: zero
pixels, bit-identical) -- a cut that zeroed everything passes the
pinhole leg alone.  Fixture trap met on the way: a Return sphere must FACE
the incoming light (`psi` against the arrival direction) with its centre
along `+psi`, exactly the Cass deck's EP idiom; the other way round every
ray is lost and the field is silently all zero.

