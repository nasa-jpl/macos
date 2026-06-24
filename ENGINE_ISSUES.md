# MACOS engine issues — found during Leg-A grid-data / FreeForm MonZern work

Logged 2026-06-24 (sls-dev).  Status: 🔴 open · 🟡 worked-around, engine fix pending · 🟢 fixed.

These surfaced while building the `dw_dgrid` (GMI `pgrid`) sensitivity channel and
its grid-vs-MonZern cross-check (`mmacos/examples/sensitivities/e5hex1`).  The
cross-check is the diagnostic: a Zernike applied as **grid data** must equal the
same Zernike applied as a **MonZern coefficient** (identical sag), so any
disagreement is an engine convention/bug.

---

## Resolution status (2026-06-24)
- **E1** 🟢 fixed (`iZernTypeForNorm`) + VERIFIED — MonZern RMS now matches the grid.
- **E2** 🟢 guard prints a loud rejection (was a silent return).
- **E3** 🟢 fixed — `MonZernCoef` + `FFZernCoef` try the whole line first (both chains); 7-coef-on-one-line now loads.
- **E4** 🟢 documented — `GridMat` i=+x, j=+y in `elt_srf_grid_data_add` + `NGSrf`.
- **E5** 🟢 already resolved — all 3 dispatch files have ELSE+warn + Noll (`ZerntoMon6`); the CLAUDE.md note was stale.
- **E6** 🟢 correct-by-design — Data & Mon frame defaults are identical (`VptElt`/`psiElt`) AND independently overridable; differences are allowed.
- **M1** 🔴 mmacos — rename `dw_dz_zernike` `n_zcoef` (END mode) → `zmode_end`.

All engine fixes built clean (gfortran) + mex relinked; E1/E3 verified via MATLAB.

---

## E1 🟢 MonZern / FF RMS normalization keyed off the wrong ZernType
**Where:** `surfsub.F:3381` (inside `ZerntoMon1`)
```fortran
IF (ANY(ZernType_rms_normalised == ZernTypeL(ie))) zc = zc * NORM_RMS_PARAM_ANSI
```
**Symptom:** `MonZernType=NormANSI` (and `FFZernType=Norm*`) come through
`ZerntoMon1` **un-normalised** — the √ RMS factor is skipped.  grid-vs-MonZern
RMS/P-V differ by exactly √6 (astig) / √3 (focus) / √8 (coma).  Shapes match
(cosine 0.998); only the scale is wrong.
**Cause:** `ZerntoMon1` is shared by all three coefficient paths but hardcodes
`ZernTypeL(ie)`.  For the Mon coefs the governing type is `MonZernTypeL(ie)`, for
the FF coefs `FFZernTypeL(ie)`.  The dispatch (tracesub/propsub/srtrace) already
*selects the variant* off the correct type — only the norm gate is wrong.
**Fix (this commit):** module integer `iZernTypeForNorm` in MODULE surfsub, set by
each dispatch path (element-Zern→`ZernTypeL`, FF→`FFZernTypeL`, Mon→`MonZernTypeL`)
just before the `ZerntoMon*` call; `ZerntoMon1` gates on it (−1 sentinel falls back
to `ZernTypeL(ie)`) and resets it.  Norm stays applied *inside* `ZerntoMon1`
(post-permutation, ANSI order — moving it out would mis-order it for
BornWolf/Noll/Fringe).  9 set-sites (3 paths × tracesub/propsub/srtrace), no
signature change, no external callers.

## E2 🟡 Grid-data add overflows the GridMat slot when nGridMat > mGridMat
**Where:** `macos_api_mod.F90` `elt_srf_grid_data_add`; root: `GridMat(mGridMat,
mGridMat,mGridSrf)` is model-tied (`mGridMat`=128/256/1024 from macos_param.txt).
**Symptom:** model-128 (`mGridMat`=128) + `nGridMat`=256 writes `GridMat(1:256,…)`
past the 128-slot → corrupts the neighbouring element's grid (a poke on one
segment leaks into others).
**Fix:** guard `Nx>mGridMat → return` ADDED (held; engine rebuild + verify pending).
TODO: make it a loud error not a silent return; consider guarding the engine's own
grid-set paths too.

## E3 🔴 >6 modes on one `MonZernModes=` line fails to load
**Symptom:** `MonZernModes= 4 5 6 7 8 9 10` (7 modes) → `load_rx failed`.
**Cause:** TBD — likely the fixed-format row width in `msmacosio.inc`
(`MonZernModes`/`MonZernCoef` reader), or a `MonZernCoef` value-format edge.
**Fix:** TBD — investigate the MonZern* parser; either support continuation rows
or emit a clear "max N modes per line" error.

## E4 🟡 GridMat array orientation (i=+x, j=+y) undocumented → coma transpose trap
**Where:** `surfsub.F` `NGSrf` (~2823): first array index maps to +x, second to +y.
**Symptom:** a MATLAB `meshgrid` map (rows=y) handed to `elt_grid_add` is x↔y
transposed vs the analytic Mon; EVEN modes survive, ODD (coma/trefoil) orthogonalize.
**Fix:** worked around in mmacos (`zernike_grid_basis` uses `ndgrid`).  TODO: add a
one-line convention comment at `NGSrf` and in `elt_grid_add` so future callers know.

## E5 🔴 Zernike-apply dispatch silent no-op when ZernTypeL = 0 (pre-existing)
**Where:** `propsub.F:230-249`, `srtrace.F:144-158` (and the parallel chains) — the
IF/ELSEIF over `ZernTypeL` has no ELSE, so `ZernTypeL=0` (uninitialised) silently
skips `ZerntoMon`, leaving MonCoef zero → silent zero-response.  Also unhandled:
`ExtFringe` (11).
**Fix:** add an ELSE that warns (and/or default `ZernTypeL` to ANSI at parse time
for paths that miss the ChkDf2 default).

## E6 🔴 Grid (Data) frame default for FreeForm segments
**Symptom:** stock `e5hex1.in` carries `pData/xData/yData/zData` at the PM vertex
(global) for every segment, so a per-segment grid Zernike mis-centers / mis-rotates;
had to set the grid frame = the MonZern frame per segment by hand.
**Question:** should `ChkDf2` default the Data frame to the element's own
frame (`pMon/xMon/…` or `RptElt`) for a segmented FreeForm, rather than leaving it
global?  TODO: decide + either fix the default or document that the Data frame must
be set per segment.

---

## Non-engine (mmacos)

## M1  `dw_dz_zernike` `n_zcoef` is the END mode, not a count
**Where:** `mmacos/src/+macos/dw_dz_zernike.m:71` `target_modes = zmode_start:n_zcoef`.
**Symptom:** `n_zcoef=6` gives modes 4:6 (no coma columns) — misleading name; cost
real debugging time on the coma cross-check.
**Fix:** rename to `zmode_end` (or change to `zmode_start + n_zcoef - 1`) and update
callers/docstring.
