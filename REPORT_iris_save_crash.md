# The IRIS save_rx -> reload SIGSEGV: a phantom-grid SAVE bug (2026-09-09)

CCMac's round-2 item 3, root-caused and fixed.  Off the sens-core merge
path (the merge stands); a pre-existing crash in the SAVE writer since
the July element-data bucket work (macos `662e86e`).

## Symptom

`iris_dp_ZGD`: load -> (any or no stop) -> `save_rx` -> reload SIGSEGVs
in `SFFSrf -> FreeFormSrf -> MonGridSrf -> FindSrf -> CTRACE`.  The
ORIGINAL deck loads and traces (12760 rays); only the round-trip
corrupts it.  Stop-independent (CCMac 3a), so not a stop bug.

## Root cause (two independent defects, both from the nGridMat>0 idiom)

IRIS declares 38 elements with `nGridMat= 99` and `GridFile= None` -- a
grid slot with NO data loaded (`ifGridDataDefined` false).  Two SAVE
sites and one trace site keyed on `nGridMat > 0` alone, which is true
for these phantom grids:

1. **SAVE over-emits the grid frame.**  `iosub.inc`'s element-data
   bucket wrote `pData/xData/yData/zData` (and `GridSrfdx`) for ANY
   `nGridMat>0` element.  The 38 phantom elements gained a frame they
   never had in the source; `GridSrfdx` was written as its default 0.
2. **The trace activates the grid term on a zero pitch.**  In
   `surfsub.F` (`SFFSrf` and `FreeFormSrf`), `ifGridTerm = nGridMat>0`.
   With a frame present and pitch 0, the grid index
   `xi = (xData.rhom)/dAct` divides by zero; the NaN index reads
   `GridMat` out of bounds -> SIGSEGV (macOS) / NaN OPD (Linux).

The original deck survives because it carries NO frame on those 38, so
the frame defaults harmlessly; the round-trip is what materializes the
crashing frame.

## Fix

- **`surfsub.F` (both grid-term sites):** `ifGridTerm` now also
  requires the pitch positive (`GridSrfdx > 0` / `dAct > 0`).  A zero
  pitch cannot index a grid, so a declared-but-empty grid is inert --
  the belt.
- **`iosub.inc` (SAVE writer, 3 sites):** the grid frame, `GridSrfdx`
  and the frame vectors are emitted only for a DEFINED grid
  (`GridDefinedElt(iElt)` = `nGridMat>0 .AND. ifGridDataDefined(slot)`)
  -- the suspenders.  `nGridMat`/`GridFile` are still preserved (they
  are real in-memory state); a blank `GridFile` is written as the
  `none` sentinel GridInit already accepts, so the SAVE re-loads
  (a blank value trips the prescription validator).

## Gate (synthetic; the real IRIS deck is JPL-private -- CCMac confirms)

`ZGD_test_files` decks + two synthetic phantoms (a FreeForm element with
`nGridMat= 99` and `GridFile= None`, matching the IRIS idiom):

| deck | class | load->save->reload->OPD | SAVE vs pre-fix |
|---|---|---|---|
| tst_FF_fg / mg / g | FreeForm + real grid | OPD unchanged | byte-identical |
| FFSegDemoData | 6 FreeForm segments + grids | OPD unchanged | byte-identical |
| SegDemo3data | 6 grid segments | OPD unchanged | byte-identical |
| nsA_conicgrid | grid on a Conic NSReflector (July case) | OPD unchanged | byte-identical |
| nsB_ffgrid | FreeForm+Zernike+grid NSReflector | OPD unchanged | byte-identical |
| phantom2 (nGridMat=99, GridFile=None) | the IRIS idiom | reloads clean; OPD == grid-free twin | frame + GridSrfdx dropped |
| phantom (nGridMat=99, no GridFile) | blank-name variant | reloads clean (was validator-refused) | frame dropped; GridFile= none |

Real grids keep their frame (the July fix stands); the phantom grids
lose the frame that crashed on reload and trace inert.  Both compilers
(ifx + gfortran); SAVE byte-identical across the two.  A no-trace
load->save of every deck is byte-identical pre/post fix.

## NOT this bug (found alongside, pre-existing, separate)

Tracing `tst_save_keys.in` (a broken all-keys fixture that traces to
NaN, and the only local deck with `LensArrayIndRef=`) leaves a
`LensArray` index value (1.51242597) stamped on the NEXT element's
`IndRef` -- a TRACE-time lensarr overrun (the `lensarr_indexes.inc`
heap-stomp class the July `662e86e` note fixed at parse time but not on
this path).  Trace-triggered, not SAVE-triggered: load->save is clean,
load->opd->save shows it, on the same binary.  Independent of this fix;
filed as its own PLAN section 0 item.
