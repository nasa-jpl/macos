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
`RayTgtElt`, `DpElt`, `PolSrc`, `PolBeam`, `UDBeam`,
`UDSrcProf`, `xFrame`/`yFrame`/`zFrame` (source local frame — flagged
`SrcLF_FLG`), `ObjFile`.

**`Ex0Ey0` — RESOLVED (polarization Phase 1, 2026-07-25).** The source
polarization state (`Ex0`,`Ey0`) now parses (`msmacosio.inc`, four reals
`ExRe ExIm EyRe EyIm`) and round-trips through SAVE (`PrtSourceInfo`,
after `Flux=`), gated on a non-zero state so polarization-free Rx SAVE
output stays byte-identical.  Carries STATE only; polarization on/off is
an API/CLI concern (`pol_set` / `POLARIZATION`), not an Rx keyword.
`PolSrc`/`PolBeam` remain unrevived (commented block in `msmacosio.inc`).

## Element data — RESOLVED 2026-07-04 (all keys now round-trip)
Every key below is now written by SAVE, gated on the value actually
being set in memory (prescriptions that never used the key produce
byte-identical output — verified against the 2026-07-03 e5hex1
baseline).  Fixture: `ZGD_test_files/tst_save_keys.in` (+
`tst_save_ampl.dat`) exercises all of them; load→SAVE→reload→SAVE
is byte-identical under ifx AND gfortran, and ifx==gfortran.

| Key(s) | Gate | Notes |
|---|---|---|
| `Coating` | `EltCoat>0` | thickness UN-scaled by ×IndRef/Wavelen (inverts parse-time scaling) |
| `GradInd GradCoef GradLensZ` | `IsGrinRefElt` | coef count mirrors parse (type 1/6→3, 7→5) |
| `DoeWL DoePhase DoeFl` + `OrderHOE` | new `CASE (DoeTrGratingElt)`, `DoeId>0` | DOE had NO writer at all |
| `nAmplMat AmplFile AmplSrfdx` | `ifAmplSrf` | nAmplMat precedes AmplFile (AmplInit needs it) |
| `LensArrayIndRef` (also covers `XYIndRefFile` input) | `IsVarIndRefElt` | canonical explicit-value form; file name not retained |
| `ArrIndRef` | `nOptWavelen>1 .AND. ANY(/=0)` | one line (parser reads from VALUE) |
| `ArrWaveLen` (header) | `nOptWavelen>1` | writes λ2..λn only — `Wavelen=` restores slot 1 on reload |
| `SegApType SegApVec` (header) | `SegApType==6` | after SegCoord (parser ordering) |
| `ZernCenter ZernXDir ZernYDir ZernRad` | `zernUsrOpt` | after VptElt so the default→override order holds |
| `ZernAnnularRatio` | typeL==NormAnnularNoll (Zern OR FF/Mon) | **was written misspelled `ZernAnnualRatio` since the writer was born — silently lost on every reload** |
| `ZCOZernType` (header) | `zcoType/=0` | NormAnnular / NormHex |
| `GridSrfOrder` | `nGridMat>0 .AND. ==3` | default 1 implied |
| `lData` | `/=0` | parse-and-store |
| `pData/xData/yData/zData` | widened to ANY `nGridMat>0` | SrfType 9/11 use the frame in the trace (GridSrf fix) but SAVE only wrote it for 12/13/FF |
| `nGridMat GridFile GridSrfdx` on NON-grid SrfTypes | `nGridMat>0` outside grid-surface CASEs | **Conic NSReflector segments with GridData figures (iris_dp_ZGD) lost their whole grid on SAVE** |
| `nMetPos tMetElt` | `nMetPos>0` | nMetPos precedes tMetElt (parser STOPs otherwise) |
| `EdgeSensors` | `EdgeSensor>0` | 9-value row per sensor |

No action (already canonical): `MrEltGrp` (saved as `EltGrp`),
`PolyApVec PolyObsVec Poly3DObsVec` (2D projections saved as
ObsVec/ApVec vertices; 3D aperture form via `Save3DApVec=Y`; no 3D
re-save for obscurations — 2D is equivalent).

Rode along:
- **`lensarr_indexes.inc` heap stomp**: the rectangular lenslet
  table hardwired ndim=107 → 11449 lenslets vs `mLenslet=250` —
  `RecLensletIdx(500)` overflowed by 22k ints the moment any Rx
  used `LensArrayIndRef=`/`XYIndRefFile=` (zero corpus coverage;
  broken since the 107×107 table replaced the 5×5).  Now sized
  from mLenslet (ndim=15), plus a clamp in InitLensletIndexAndCtr
  and a count guard at the LensArrayIndRef parse site.
- **FmtD negative-zero normalization**: computed frame defaults can
  carry −0.0 under one compiler and +0.0 under the other; `-0.0E+00`
  now prints as `0.0E+00` so ifx/gfortran SAVEs stay byte-identical.
- Re-save settle: `psiElt` is re-DUNITIZEd on every load, so a
  17-digit unit vector can drift 1 ulp on the FIRST re-save and is
  stable from then on (same pre-existing class as the ChfRayPos
  re-aim under ApStop; not a SAVE defect).

## Regeneration snippet
```python
import re
parsed = {m.group(1) for m in re.finditer(r"LCMP\(VAR_NAM,\s*'([^']+)'", open('msmacosio.inc').read())}
written = {m.group(1) for m in re.finditer(r"PrintFmt(?:Array)?\(\s*'([^']+)'", open('iosub.inc').read())}
cov = lambda k: any(w.lower().startswith(k.lower()) or k.lower().startswith(w.lower()) for w in written)
print('\n'.join(sorted((k for k in parsed if not cov(k)), key=str.lower)))
```
