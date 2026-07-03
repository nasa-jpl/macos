# SAVE keyword audit — parser-accepted keys the SAVE path does not write

Generated 2026-07-03 for PLAN §0 item 2 ("Add `ApStop` to the SAVE
path; review the other not-currently-saved variables").  Method: every
`LCMP(VAR_NAM,'…')` keyword in `msmacosio.inc` (201 keys) diffed
against every `PrintFmt`/`PrintFmtArray` name written by the SAVE
writers in `iosub.inc` (102 keys), prefix-matched both ways.
Regenerate with the snippet at the bottom.

**Status: 101 keys not written.**  Categorized below for review —
most are fine by design; the "element data" bucket is where real
load→SAVE data loss may hide.

## Fixed 2026-07-03 (this slice)
| Key | Form | Written by |
|---|---|---|
| `ApStop` (header, 3-vector `StopPos`) | `IF (ifStopSet .AND. .NOT.EltStopSet)` | `PrtSourceInfo` |
| `ApStop` (element-bound `dx dy [auto]`) | `IF (EltStopSet .AND. iElt==StopElt)` | `PrtSingleEltInfo` |
| `OPDRefRayLen` | `IF (OPDRefRayLen_FLG)` | `PrtSourceInfo` |
| `RxNoStopSet` | `IF (RxNoStopSet_Flg)` | `PrtSourceInfo` |

Also fixed en route: **`macos_IO.f90` built a malformed runtime format
`(A17,'=',,…A)` (double comma) at six `cmd` assembly sites — ifx
silently tolerated it, gfortran hard-errored ("Unexpected element ','
in format"), so every gfortran-build SAVE (CLI and `save_rx`) crashed.**
ifx/gfortran SAVE output now byte-identical (e5hex1 diff).

## By-design canonical alternatives — no action
The parser accepts these as *input* conveniences; SAVE writes the
canonical form instead: `fElt`, `eElt` (→ `KrElt`/`KcElt`), `EltType`
(→ `Element`), `CopyElt`, `DuplElt`, `LINK` (block macros expand at
load), `JZLouChfRayDir/Pos`, `JZLouWavelen` (debug aliases).

## Run/debug/output toggles — probably skip (session state, not Rx)
`DumpRayPosH`, `SaveRayPosH`, `SaveVis3dDat`, `SaveOPDMap`,
`PgplotImage`, `PrtEltRots`, `ShowMetData`, `drawXGrid`, `drawYGrid`,
`Regrid`, `NSCount`(?).

## Optimizer / CALIB configuration — design decision needed
Whole `Opt*` family (25 keys): `OptAlg OptAsph OptBeam OptBeamDir
OptBeamPos OptBeamRefRayDir OptBeamSize OptChfRayDir OptChfRayPos
OptFEX OptFOVWt OptMaxItrs OptNomSens OptPinvFile OptRayGrid
OptSavePinv OptSpotSize OptSrcRpt OptSvdSvCut OptTarget OptTgt
OptTgtElt OptTgtWF OptUseSavedPinv OptWaveLen OptWFElt OptZern
OptZernModes` plus `VarElt VarZern VarCons VarSrc CalBWK`.
These configure CALIB runs; whether a SAVE'd Rx should carry them is
a policy call (they round-trip a *workflow*, not the optical system).

## Trace / OPD reference state — review individually (like ApStop)
`FEXCentroid` (`Rx_FEXCentrFlg`), `UseChfRay4OPD` (`LUseChfRayIfOK`),
`RayTgtElt`, `DpElt`, `Ex0Ey0`, `PolSrc`, `PolBeam`, `UDBeam`,
`UDSrcProf`, `xFrame`/`yFrame`/`zFrame` (source local frame — flagged
`SrcLF_FLG`), `ObjFile`.

## Element data — POTENTIAL LOAD→SAVE DATA LOSS, verify per key
`Coating` (nCoat written? verify), `GradInd GradCoef GradLensZ`
(gradient-index), `DoePhase DoeWL DoeFl` (DOE gratings),
`AmplFile AmplSrfdx nAmplMat` (amplitude grids), `XYIndRefFile`,
`LensArrayIndRef ArrIndRef ArrWaveLen`, `SegApType SegApVec`
(segment apertures), `ZernCenter ZernRad ZernXDir ZernYDir
ZernAnnularRatio ZCOZernType` (Zernike frame overrides — note
`ZernType= NormAnnularNoll <ratio>` may lose the ratio on SAVE),
`GridSrfOrder`, `lData` (grid-frame length — `pData/xData/yData/zData`
are written for FreeForm; `lData` is not), `PolyApVec PolyObsVec
Poly3DObsVec` (2D projections are saved as ObsVec; 3D forms need
`Save3DApVec=Y` — see Lou item 47), `EdgeSensors nMetPos tMetElt
MrEltGrp` (metrology).

## Regeneration snippet
```python
import re
parsed = {m.group(1) for m in re.finditer(r"LCMP\(VAR_NAM,\s*'([^']+)'", open('msmacosio.inc').read())}
written = {m.group(1) for m in re.finditer(r"PrintFmt(?:Array)?\(\s*'([^']+)'", open('iosub.inc').read())}
cov = lambda k: any(w.lower().startswith(k.lower()) or k.lower().startswith(w.lower()) for w in written)
print('\n'.join(sorted((k for k in parsed if not cov(k)), key=str.lower)))
```
