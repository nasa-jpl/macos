# MACOS Source Tree


NASA/JPL optical ray tracing code. Legacy Fortran, some files date to the 1980s.
Fixed-form source: .F files use the C preprocessor, .f files do not.

## Build
- Full build (npsol, pgplot, readline, fitsio, macos, smacos): source ./makealldcr.sh
- macos/smacos only (macos_f90 changes): source ./makeMSdcr.sh
- GMI mex (Matlab interface): source ./makeGMIdcr.sh
- All scripts use /opt/intel/oneapi/compiler/latest/ for portability.

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
- All build scripts accept an optional `debug` argument:
  source ./makeMSdcr.sh debug   # macos + smacos with -O0 -check all
  source ./makeSDdcr.sh debug   # smacos_dvr with -O0 -check all
- Makefile uses `OPT ?= -Ofast`; scripts export `OPT="-O0 -check all"` for debug.
- VS Code launch configs in macos_f90/.vscode/launch.json:
  "Debug MACOS" (interactive macos) and "Debug SMACOS_DVR" (smacos driver).
  Requires ms-vscode.cpptools extension.
- smacos_dvr: compiles smacos_dvr.F with -DCMACOS, links against smacos_lib.a
  plus macosio.o/pgplotsub.o/macos_vars_mod.o/macos_mod.o from MACOS_OBJS.
  Requires prior make macos + make smacos.

## Build notes
- Never leave a surfsub.f (lowercase .f) in macos_f90/ alongside surfsub.F — make
  picks it up via implicit rules and tries to build it as a standalone executable.
- Full link requires environment from makealldcr.sh (macossrc_dir, intel64_lib, etc.).
  Running make macos directly will compile but fail at link due to missing paths.
- Use makeMSdcr.sh for macos_f90-only changes (faster than makealldcr.sh).
- GMI (Matlab mex): built separately via makeGMIdcr.sh.
  ifx 2025.3 quirk: -C flag implies -fsanitize=memory; use -check all instead.
- jGridSrf mapping: tracesub.F, propsub.F, srtrace.F use nGridMat(iElt).GT.0
  (not SrfType checks) so all grid-using surfaces get the correct GridMat slot.

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