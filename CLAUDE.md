# MACOS Source Tree


NASA/JPL optical ray tracing code. Legacy Fortran, some files date to the 1980s.
Fixed-form source: .F files use the C preprocessor, .f files do not.

## Build
- Full build (npsol, pgplot, readline, fitsio, macos, smacos): source ./makealldcr.sh
- macos/smacos only (macos_f90 changes): source ./makeMSdcr.sh
- GMI mex (Matlab interface): source ./makeGMIdcr.sh
- All scripts use /opt/intel/oneapi/compiler/latest/ for portability.

### CMake build (alternative)
- Build script: source ./makeCMdcr.sh [debug] [gfortran]
  - `source ./makeCMdcr.sh`                — Release ifx
  - `source ./makeCMdcr.sh debug`          — Debug ifx (-O0 -check all)
  - `source ./makeCMdcr.sh gfortran`       — Release gfortran
  - `source ./makeCMdcr.sh debug gfortran` — Debug gfortran
- Each combination gets its own build directory (build_release_ifx, build_debug_gfortran, etc.).
- CMakeLists.txt files: macos_f90/CMakeLists.txt (top-level), macos_f90/npsol/CMakeLists.txt.
- CMakePresets.json: debug and release presets for VS Code CMake Tools integration.
- Targets: macos (executable), smacos_lib (static library), smacos_dvr (executable, -DBUILD_SMACOS_DVR=ON), GMI (mex, -DBUILD_GMI=ON).
- pgplot and fitsio are pre-built externals — build them first with their own scripts before running CMake.
- C compiler must be gcc (not icx) — legacy C files use implicit function declarations.
- smacos_dvr re-compiles macos_mod.F with -DCMACOS (smacos_lib's copy lacks CMACOS-only symbols like ifPGColor).

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

## Next tasks
1. Verify surfsub.F and elt_mod.F compile cleanly -- done
2. Add UI input for lFF, pFF, xFF, yFF, zFF, FFZernCoef, MonZernCoef -- done
3. Prescription file output -- done
4. Add FreeFormSrf calls in elemsub.F (FindSrf + Reflector) -- done
5. Fix OLD command surface type loop (was DO i=1,13; now mSrfType) -- done
6. Fix OLD prescription reader: FreeForm keywords in element chain -- done
7. Fix pData/xData/yData/zData output and ChkDf2 defaults for FreeForm -- done
8. Fix UI dialog: pData/xData/yData/zData prompts for FreeForm+grid -- done
9. Reorder NEW dialog: FF-first (FFZernType→FFZernCoef→lFF→Mon→Grid→coords) -- done
10. Reorder prescription output: FF-first -- done
11. Fix jGridSrf mapping for grid-using surfaces (was only SrfType 9/11) -- done
12. Add FreeForm to Refractor in elemsub.F + update call sites -- done
13. Fix makeGMIdcr.sh: remove -fsanitize=memory, replace -C with -check all, remove eval -- done
14. SHOW UI display during input -- deferred
15. Fix FreeFormSrf normal: zc must use total sag (fh+fhFF), not just Mon sag -- done
16. Add pFF/xFF/yFF/zFF to all PERTURB routines -- done
17. Fix pData perturbation condition (was <=13, now ==12/==13/==FreeForm) -- done
18. Add MODIFY handlers for nFFZernCoef, FFZernModes, nMonZernCoef, MonZernModes -- done
19. Archive/remove obsolete source files -- done
20. CMake build system (CMakeLists.txt, presets, VS Code integration, makeCMdcr.sh) -- done
21. Per-ray status tracking + convert STOP/PAUSE to graceful failures -- done
22. Fix short filename crash in OLD command (macosio.F substring bounds) -- done
23. Giza plot sharpness (supersampling + HiDPI via cairo_surface_set_device_scale) -- done
24. Giza window resize (no blank on resize; 2× render → 1× pixmap → scaled copy to window) -- done
25. Readline arrow keys for CMake build (-DREADLINE_LIBRARY in top-level CMakeLists.txt) -- done
26. Embed Intel RPATH so ifx-built macos runs without sourcing setvars.sh -- done
27. Giza close-button handling: unmap window instead of closing device -- done
28. Linefeed after MODIFY Q exit so next MACOS prompt doesn't overwrite -- done
29. Embed Intel RPATH in GMI.mexa64 (CMakeLists.txt) so MATLAB loads it without setvars.sh -- done

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

## Intel RPATH (self-contained ifx binary + GMI mex)
- Top-level CMakeLists.txt embeds RPATH for ifx builds so `/opt/intel/oneapi/...` libs
  are found without sourcing setvars.sh. Applies to both macos executable and GMI.mexa64:
  `BUILD_RPATH/INSTALL_RPATH` = `${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib`,
  `-Wl,--disable-new-dtags` forces DT_RPATH (transitive) over DT_RUNPATH (direct-only).
  libimf → libintlc transitive dep requires DT_RPATH.
- For gfortran builds: `CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES` (auto-set by CMake) used
  as RPATH so libgfortran/libquadmath are found without a module file.

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
- Makefile: source ./makeMSdcr.sh debug   # macos + smacos with -O0 -check all
- CMake: source ./makeCMdcr.sh debug      # macos + smacos + smacos_dvr
- Makefile uses `OPT ?= -Ofast`; scripts export `OPT="-O0 -check all"` for debug.
- CMake debug uses -check all,noarg_temp_created (suppresses harmless array temporary warnings).
- VS Code: F5 launches debug session; builds automatically via preLaunchTask.
  launch.json points to build_debug_ifx/macos and build_debug_ifx/smacos_dvr.
  tasks.json runs cmake --preset debug. Requires ms-vscode.cpptools extension.
  `debug.allowBreakpointsEverywhere: true` in settings.json enables breakpoints in .F files.
- smacos_dvr: compiles smacos_dvr.F with -DCMACOS, links against smacos_lib.a
  plus macosio.o/pgplotsub.o/macos_vars_mod.o/macos_mod.o from MACOS_OBJS.

## Build notes
- Never leave a surfsub.f (lowercase .f) in macos_f90/ alongside surfsub.F — make
  picks it up via implicit rules and tries to build it as a standalone executable.
- Full link requires environment from makealldcr.sh (macossrc_dir, intel64_lib, etc.).
  Running make macos directly will compile but fail at link due to missing paths.
- Use makeMSdcr.sh for macos_f90-only changes (faster than makealldcr.sh).
  makeMSdcr.sh runs `make clean-macos` then `make macos` — this removes the old binary
  first. If the link fails (e.g. missing propsub_mod128), the binary is gone; use CMake
  instead: source ./makeCMdcr.sh debug (from macos/, not macos_f90/).
- macos.F line 59: `use propsub_mod128` was a stale name; fixed to `use propsub_mod`.
- GMI (Matlab mex): built separately via makeGMIdcr.sh.
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