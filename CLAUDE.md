# MACOS Source Tree


NASA/JPL optical ray tracing code. Legacy Fortran, some files date to the 1980s.
Fixed-form source: .F files use the C preprocessor, .f files do not.

## Build
CMake-based. All scripts live in `macos/` and accept `[giza|pgplot] [debug|release] [gfortran]`
in any order (defaults: giza, release, ifx). Each combination gets its own build directory.

| Script | Targets |
|---|---|
| `source ./makejoint.sh` | macos + smacos + smacos_dvr + GMI (all four) |
| `source ./makems.sh` | macos + libsmacos.a |
| `source ./makesd.sh` | smacos_dvr (builds macos/smacos too if needed) |
| `source ./makegmi.sh` | GMI.mexa64 (requires macos+smacos built first) |
| `source ./makegfortran.sh` | macos + smacos + smacos_dvr via gfortran |

Build directory naming: `build_{release|debug}_{giza|pgplot}[_gfortran]`

- CMakeLists.txt: top-level `CMakeLists.txt` (macos_f90 sources), `macos_f90/npsol/CMakeLists.txt`.
- CMakePresets.json: debug and release presets for VS Code CMake Tools integration.
- C compiler must be gcc (not icx) — legacy C files use implicit function declarations.
- smacos_dvr re-compiles macos_mod.F with -DCMACOS (smacos_lib's copy lacks CMACOS-only symbols like ifPGColor).
- GMI.mexa64 is built via the standalone `MACOS_resources/GMI/Makefile` (not cmake) — more robust across MATLAB versions. makejoint.sh and makegmi.sh both call it with `MACOS_BUILD_DIR` pointing to the cmake build tree.

## Key files for current work
- macos_f90/elt_mod.F        : per-element data arrays and SrfType constants
- macos_f90/surfsub.F        : all surface intersection routines
- macos_f90/elemsub.F        : calls surface routines during ray trace (FindSrf, Reflector, Refractor)
- macos_f90/param_mod.F      : array dimension parameters (mElt, mGridMat, etc.)
- macos_f90/iosub.inc        : ChkDf2 defaults, PrtSingleEltInfo output (included by macosio.F)
- macos_f90/msmacosio.inc    : prescription file reader (included by macosio.F and smacosio.F)
- macos_f90/macosio.F        : interactive UI dialog for element entry
- macos_f90/tracesub.F       : ray trace loop (Reflector + Refractor + FindSrf call sites)
- macos_f90/propsub.F        : propagation (Reflector + Refractor + FindSrf call sites)
- macos_f90/srtrace.F        : single-ray trace (Reflector call sites)
- macos_f90/funcsub.F        : CPERTURB, CPRead, CPERTURB_GRP (perturbation routines)
- macos_f90/macos_ops.F      : CPERTURB_2 (perturbation, macos_ops_mod)
- macos_f90/lnk_pert.inc     : LnkEltCPERTURB (linked-element perturbation)

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


## Per-ray status tracking
- RayStat_* constants in elt_mod.F: OK(0), Obscured(1), Miss(2), Bracket(3), MaxIter(4), Undef(5).
- Per-ray arrays: RayStatus(mRay), RayFailElt(mRay), RayFailMsg(mRay) — allocated in elt_mod.
- SetRayFail(iRay, iStat, iElt, cMsg) records first failure only (avoids overwriting root cause).
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

## Display polarity + plot annotations
- `plot_mode_mod.F` (registered early in `COMMON_SOURCES`) holds session
  state:
  - `ifAstroMode` (default `.TRUE.` = astronomy / negative polarity,
    large→dark — matches legacy PGPlot users' expectation)
  - Bottom-of-plot annotation: `annotOn`, `annotLine1`, `annotLine2`
  - `wedgeLabel` (default `'pixel value'`) — text under the color wedge
  - `ClearAnnot()` resets all three (annotOn + both lines + wedgeLabel)
    after the draw routine emits them.
- `IMGMODE` command toggles polarity interactively (prompts `NEG|POS`,
  accepts `ASTRO|CONV` as synonyms). Must be dispatched BEFORE the
  `LCMP(command,'IMG',3)` branch in `macos_cmd_loop.inc` — placed at
  ~line 4362 for that reason (otherwise "IMGMODE" is absorbed by the
  3-char 'IMG' match and treated as a plot command requiring `ifLoad`).

### Raster inversion — the key Giza gotcha
Giza's `giza_render_gray` (backing `PGGRAY`) internally calls
`giza_set_colour_table_gray()` and forcibly resets its own ramp,
**ignoring any `PGSCR` setup the caller does** (see
`macos_f90/giza/src/giza-render.c:368`). You cannot use `PGGRAY` for
polarity control. The fix in `pgplotsub.F:GRAY`:
- Render the image via `PGIMAG` + a 2-point `PGCTAB` we build ourselves
  (both the color and gray paths).
- For gray: `CTAB_R/G/B = [1.0, 0.0]` (white→black) under astro mode,
  `[0.0, 1.0]` under POS.
- For color: walk `CInd` in reverse under astro mode (palette flipped
  end-to-end; hue order preserved).

### The second Giza gotcha: `PGCTAB` clobbers CI 0 and 1
`giza_set_colour_table` (called by `PGCTAB`) ends with
`_giza_set_range_from_colour_table(_giza_colour_index_min,
_giza_colour_index_max)` — that range starts at **0**. A 2-point ramp
written with `L=[0,1]` overwrites CI 0 and CI 1 with the ramp endpoints,
clobbering the bg=white / fg=black that axes, tick labels, wedge axis,
and `PGMTXT` rely on. Symptoms (all debugged this way): axis-label text
invisible on first plot, annotation visible but axes ticks gone on a
second plot, etc.
- Fix: **reserve CI 0=white and CI 1=black** via explicit `PGSCR` at
  start of GRAY, then `PGSCIR(16, 255)` pushes the image ramp into CI
  16..255 so `PGCTAB` can't touch CI 0/1. Always use `PGSCI(1)` (never
  `PGSCI(ICILO)`) for axes/box/labels/wedge/annotation.
- Re-assert `PGSCI(1)` before `PGWEDG` and before `PGMTXT` — the
  intervening `PGCTAB`/`PGIMAG` may shift the current CI.

### Non-raster plots — always black ink on white
`CONTOUR`, `SLICE`, `WIRE`, `SPOTDIAG`, `PLOTCOL` all set
`PGSCR(0,1,1,1); PGSCR(1,0,0,0); PGSCI(1)` at entry — overrides any
`IMGMODE` setup. Harmless under PGPlot (same as its default).

### Labels, titles, wedge
- `PGLABEL(' ', ' ', CTITLE)` in GRAY suppresses the placeholder
  `'X-Axis'` / `'Y-Axis'` arg strings (the 76 callers that pass those
  literals are untouched — GRAY just ignores them). Only the title
  survives above the plot.
- `PGMTXT('B', 3.5, 0.5, 0.5, line)` centers the annotation under
  the x-axis (coord=0.5 is middle, fjust=0.5 is center-justify).
- Tick numeric labels come from `PGENV(...,axis=0)` → internally
  `PGBOX('BCNST', ...)` — the `N` option prints numeric values at
  each major tick. Keep as-is; for OPD the tick values are projected
  dimensions on the reference surface (usually the OPD reference,
  e.g. Elt 9 in SegDemo3) which is sometimes meaningful.
- Titles:
  - OPD → `'OPD, Surface N'` (no `File=` suffix)
  - INT / PIX → `'<kind>, Surface N, File=<filnam>'` (via `SrfOut`
    and `PixOut` in `utilsub.F`; `File=` kept because the in-file
    stretch variant — `LOG10 Intensity` etc. — can distinguish runs)
- Wedge labels (`wedgeLabel` in `plot_mode_mod`):
  - OPD with BaseUnits → `'OPD (<cUnit>)'` (e.g., `OPD (m)`)
  - OPD without BaseUnits → `'OPD'`
  - INT → `'Intensity'`
  - PIX → `'Pixel value'`
  - Others → default `'pixel value'`

### Annotation lines (bottom of plot)
Command handler sets `annotLine1/2` + `annotOn=.TRUE.` before the draw
call. `GRAY` and `SPOTDIAG` emit via `PGMTXT` after `PGLABEL`, then
`ClearAnnot()`.
- OPD → `OPD=<rms> <cUnit> RMS, <pv> <cUnit> P-V` (from `RMSWFE`,
  `WFEPV`, `cUnit`; falls back to no-unit form if `.NOT.BaseUnits_FLG`)
- SPOT → `RMS spot radius=<r> <cUnit>` — computed from `RaySpot`
  about the centroid in the SPOT handler (not a pre-existing var).
- INT → `Peak pixel=<MaxInt>` (from `Ca2Int` output).
- PIX → `Peak pixel=<maxval(PixArray)>`.
- Phase/Amplitude/etc. intentionally left unannotated.

## Giza plot improvements (macos_f90/giza/)
- Supersampling: off-screen image surface is 2× (GIZA_XW_SUPERSAMPLE) in giza-driver-xw.c;
  `cairo_surface_set_device_scale(surface, 2, 2)` keeps MACOS logical coordinates at 1×
  while rendering at 2× physical pixels.  Downsample to 1× pixmap on flush, then XCopyArea
  or cairo-scale to the window (handles user resize).  Text/font antialiasing: GRAY hinting
  enabled in giza-drivers.c; image data uses CAIRO_FILTER_NEAREST in giza-render.c
  (sharp pixel boundaries for OPD/PIX/INT raster arrays).
- Idle window responsiveness: mhist.c sets rl_event_hook to `_macos_rl_event_hook` which
  (via weak `giza_process_events`) calls `_giza_process_events_xw` — drains ConfigureNotify,
  Expose, and WM_DELETE_WINDOW while readline blocks at the prompt.
- Close button (WM_DELETE_WINDOW): unmaps the window and sets XW[id].window_hidden=1;
  the next `_giza_flush_device_xw` re-maps it.  Does NOT call `giza_close_device()` —
  MACOS tracks its own `ifPlot` flag and skips `GRAINI`/`PGBEGIN` after first plot.
  Closing the device would strand MACOS with spam of "No device open" errors.
- `giza_set_colour_representation` (giza-colour-index.c) no longer requires an open
  device — `colourIndex[]` is static global, matching PGPLOT's pgscr semantics.
  `giza_set_colour_index(ci)` is still guarded (uses Cairo context).

## Intel RPATH (self-contained ifx binary + smacos_dvr + GMI mex)
- Top-level CMakeLists.txt embeds RPATH for ifx builds so `/opt/intel/oneapi/...` libs
  are found without sourcing setvars.sh. Applies to macos executable, smacos_dvr, and
  GMI.mexa64: `BUILD_RPATH/INSTALL_RPATH` = `${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib`,
  `-Wl,--disable-new-dtags` forces DT_RPATH (transitive) over DT_RUNPATH (direct-only).
  libimf → libintlc transitive dep requires DT_RPATH.
- For gfortran builds: `CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES` (auto-set by CMake) used
  as RPATH so libgfortran/libquadmath are found without a module file.

## smacos_dvr graphics initialization
- smacos_dvr.F calls `SMACOS()` directly; unlike the interactive macos.F command loop
  it never opens a Giza device. Any plot-producing SMACOS command (OPD, SPOT, INT, ...)
  then emits "No device open" for every Giza call. Fix: after `macos_init_all(modelSize)`,
  set `nPgPanel = 1` and call `GRAINI` (opens device via `PGBEGIN(0,'?',1,1)`).
- `macos_init_all` does NOT include `macos_init.inc`, so `nPgPanel` is not defaulted to 1
  — it stays at 0 and `GRAINI`'s panel-layout if/elseif (1/2/3/4) silently skips
  `PGBEGIN`. Consumers outside the interactive path must set `nPgPanel` themselves.

## Giza GRAY CR/CG/CB array bounds
- macos_f90/pgplotsub.F GRAY subroutine caches default color representations for
  restoration when exiting color mode. It calls `PGQCIR(ICILO,ICIHI)` then loops
  `Do J=ICILO,ICIHI; CALL PGQCR(J,CR(J),CG(J),CB(J))`. Under Giza, `PGQCIR` returns
  ICILO=0 (GIZA_COLOUR_INDEX_MIN), whereas classic PGPLOT starts at 1 — so CR/CG/CB
  must be declared 0-based: `REAL, Save :: CR(0:511),CG(0:511),CB(0:511)`. Same total
  size (512), but bounds now cover index 0. Otherwise the first-entry gray image crashes
  on CR(0) out-of-bounds write (hard to debug inside giza since graphics libs are
  typically stripped).

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

## Debug builds
- `source ./makems.sh debug`    — macos + smacos, -O0 -check all (ifx)
- `source ./makejoint.sh debug` — all four targets, debug
- CMake debug uses -check all,noarg_temp_created (suppresses harmless array temporary warnings).
- VS Code: F5 launches debug session; builds automatically via preLaunchTask.
  launch.json points to build_debug_giza/macos and build_debug_giza/smacos_dvr.
  tasks.json runs cmake -S/-B with debug flags. Requires ms-vscode.cpptools extension.
  `debug.allowBreakpointsEverywhere: true` in settings.json enables breakpoints in .F files.
- smacos_dvr: compiles smacos_dvr.F with -DCMACOS, links against smacos_lib.a
  plus macosio.o/pgplotsub.o/macos_vars_mod.o/macos_mod.o from MACOS_OBJS.

## Build notes
- Never leave a surfsub.f (lowercase .f) in macos_f90/ alongside surfsub.F — cmake
  may pick it up and try to build it as a standalone executable.
- Use `source ./makems.sh` for macos_f90-only changes (faster than makejoint.sh).
  cmake recompiles only changed files; delete the build directory for a clean rebuild.
- macos.F line 59: `use propsub_mod128` was a stale name; fixed to `use propsub_mod`.
- GMI (Matlab mex): built via standalone MACOS_resources/GMI/Makefile (invoked by
  makejoint.sh and makegmi.sh). MATLAB auto-detected under /usr/local/MATLAB/R*.
  ifx 2025.3 quirk: -C flag implies -fsanitize=memory; use -check all instead.
- jGridSrf mapping: tracesub.F, propsub.F, srtrace.F use nGridMat(iElt).GT.0
  (not SrfType checks) so all grid-using surfaces get the correct GridMat slot.

## Test prescriptions (ZGD_test_files/)
- All FreeForm test prescriptions use nGridpts=11 (gives 89 rays for Circular grid).
  Reason: model 256 has bRay=500 (BUILD limit) and mDrawRay=101 (nDrawElt array size);
  nGridpts=11 gives 89 rays which is under both limits.
  Grid-using prescriptions (tst_FF_g/mg/fg/refl) also need model 256 for mGridMat=256.
- Obscured rays now increment nBadRays so WARN fires: tracesub.F and propsub.F both
  have `nBadRays=nBadRays+1` inside the `IF (.NOT.L1(iRay))` / `IF (.NOT.LRayPass(iRay))`
  blocks that call SetRayFail(RayStat_Obscured).

## Conventions (new code)
- IMPLICIT NONE throughout
- Sequential DO loops only — no DO CONCURRENT
- DBLE() not FLOAT() for real conversion
- WRITE(*,*) + STOP instead of PAUSE
- Fixed-form column layout: code in cols 7-72, continuations in col 6
- Compiler flag -132 is set; lines may extend to col 132, but keep to 72 for readability
- `!-->` / `!<--` markers bracket added code blocks — follow the same pattern
- BLOCK construct is fine for local variable scoping inside ELSE IF branches

## Whitespace / editing notes
- Source files use mixed tabs and spaces for indentation (legacy; do not reformat).
  Fixed-form is column-sensitive so bulk retabbing risks introducing subtle bugs and
  destroys git blame. Leave indentation as-is.
- When the Edit tool fails due to whitespace ambiguity (tab vs spaces mismatch),
  use a Python script with exact string replacement instead:
    python3 -c "
    with open('file.F','r') as f: c=f.read()
    c=c.replace(old,new)
    with open('file.F','w') as f: f.write(c)
    "