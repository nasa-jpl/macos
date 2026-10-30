# MACOS engine notes — Engine hygiene: the closed latent-bug cluster, per-ray status, mathsub, name collisions

> **Read this file before touching:** any engine file; `SetRayFail`, `RayStatus`, `DMPROD`.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

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
- **`Get_Values` over-read (iosub.inc, 2026-09-09).**  The value
  tokenizer scanned to `MacosCharLen` (256) over a `MacosValLen` (220)
  buffer -- 36 bytes past the argument -- so `ArrWaveLen=`/`ArrIndRef=`
  lines (its only callers) parsed NONDETERMINISTICALLY: a non-blank byte
  beyond the buffer became a token and `Read(tok,*)` died "Bad real
  number in item 1" (tst_save_keys.in loaded 3 of 5 times).  Bounded by
  `LEN(ValBuf)`, token/array capped, trailing token completed; gate =
  20 consecutive loads on both compilers.  Same class as the lensarr
  heap stomp: fixed-length buffers scanned by a DIFFERENT constant.
- **`smacos_compute.inc` slice overrun.**  Five `1:mZern` → `1:mZernModes`
  fixes (mZern≠mZernModes when both Zern and FF aspheres are active).
  Cherry-picked back to release-candidate.
- **`RFOURN` (sunsub.F) takes its COMPLEX*8 array as interleaved REAL*4**
  (2026-10-18, PLAN_CONSOLIDATION 1.2 (p)).  It is Numerical Recipes'
  `fourn`: `IP1=2*NPREV`, `DATA(I3+1)` = the imaginary part.  019a2981
  (2023) declared `COMPLEX*8 DATA`, so every index addressed a whole
  complex element and the transform ran to 2*NTOT -- writing past the
  caller's array: ATMOS screens were checkerboards of seed-dependent
  amplitude (1e5x), ifx non-deterministic under a fixed SEED, every ray
  lost at a Reflector; LPI's linear-model PSF the same.  Now `REAL*4
  DATA`.  Rule: an NR-style FFT's DATA is REAL, whatever the caller
  passes -- never "tidy" it to the caller's complex type.  (DFOURN's NR
  branch converts through its own REAL*8 buffer for the same reason.)
- **A local named like a use-associated symbol** (`sp` vs `kinds`' single-
  precision kind, found in AddAtmos 2026-10-18): ifx accepted it silently,
  gfortran refused ("conflicts with symbol from module").  Build gfortran
  before trusting an ifx-only compile of new declarations.

## Per-ray status tracking
*(Documented: Manual §5 "Per-ray status reporting" (statuses, the clipping element, restarts, the bindings' getter).  Below: the engine record.)*
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
- Obscuration recorded AT THE CLIP (2026-10-12): CTRACE stamps
  `RayStat_Obscured` with the current element at the end of the element
  iteration in which `L1` went false (CPROPAGATE, element-major, always
  did).  It used to stamp after the element loop with the loop's exit value
  `iEndElt+1`, and never for a ray that later MISSED (GO TO 14 skips the
  stamp: such a ray read Miss at the later element; it now reads Obscured
  at the clip, as CPROPAGATE always said).  The end-of-ray stamp remains as
  the fallback.  WARN prints `at element N: count` under `Obscured:`.
- RayStatus/RayFailElt/RayFailMsg reset at a trace from the SOURCE only
  (`iStartElt==0`); a RESTARTED trace keeps them and sets `nBadRays` =
  rays already lost (it used to reset them and re-record every earlier
  loss at the restart's first element).  nZPFailMsg etc. still reset.
- "Obscured" is NOT "geometrically OK": the LocalCoord / FFP-PFP chief-ray
  tolerances test `LRayOK(1) .AND. RayStatus(1)==Obscured`.
  Gates `tRayFailElt`, `tGratingAperture`, `cli_tests/cmd/WARN_obscured_elt.jou`.
- WARN subroutine enhanced with per-category breakdown (miss/obscured/bracket/other).
- Message throttling: nZPFailMsg counter with mZPFailMsg=20 threshold in MODULE surfsub.
  First 20 bracket/iter messages print; rest suppressed. WARN prints suppression summary.

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


## 2026-10-19 (b): the manual's generated dialogs as a defect finder
- Debug prints found by regenerating every dialog and removed: `** In DBINOUT` (sunsub.F), `*** CPERTURB src`
  (funcsub.F), `Debug: VALUE` (macosio.F MOD Flux=), `Debug: VAR_NAM` (MOD, earlier the same day), FEX's
  `CrossPt =` / `zp_iEm1 =` (tracesub.F).  Rule: a `print*` that names a variable is a debug line; the cmd
  journals carry EXPECT-NOTs for each.
- `BINOUT` / `DBINOUT` (sunsub.F): a direct-access RECL is `INQUIRE(IOLENGTH=)` of one row, never a byte or word
  count -- gfortran counts bytes, ifx words (PLAN 1.2 (jj)).
- `CTRACE`'s composite-partials chain (DO 11) steps over an element the build did not trace (`nTr`); `ifRayElt`
  is reset per ray, `PrevElt(0)=0` (1.2 (ll)); a refused BUILD still traces (1.2 (mm)); CENTER's count is
  `nStopIter` (elt_mod, set by ChiefRayAiming) (1.2 (kk)).

## Future items from the ToDo census (2026-10-22, BRIEF_to_cleanup_round1)
- `ObjInt` is declared in macos_mod.F but allocated in utilsub.F (the reset in macos_mod_init_vars
  leaves it alone).  A tidy would declare it beside its allocation.
- macos_mod_init_vars initializes the CMACOS buffers in an `#ifdef CMACOS` block; its comment
  asks for SMACOS's initialization to be combined with it.  A refactor, not a fix.
