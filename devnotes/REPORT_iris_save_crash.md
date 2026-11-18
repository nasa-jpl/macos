# The save_rx grid-SAVE crashes: a phantom-grid bug FIXED, a real-grid bug STILL OPEN (2026-09-09)

CCMac's round-2 item 3.  Off the sens-core merge path (the merge stands);
pre-existing in the SAVE writer since the July element-data bucket
(macos `662e86e`).

> **STATUS after CCMac round 3 (`verify/REPORT_round3.md`).**  This report
> fixes ONE of two defects.  The **phantom-grid** side is FIXED and gated
> (below): the 38 `nGridMat= 99` / `GridFile= none` elements no longer
> gain a frame (pData/GridSrfdx count 44 -> 6).  But the **real
> ZrnGrData grids** on IRIS (iElt **17/19/21/35/37/39**) STILL SIGSEGV on
> reload after `save_rx` -- same stack -- and `cda178e` does NOT cover
> them.  CCMac's attribution (eliminated: not a pre-existing grid-trace
> bug -- the original with the grid activated traces fine, 12737 survive;
> not the missing file; not the added `lData`).  The trigger: `save_rx`
> CLEANS the original `GridFile=` TAB-comment that had silently disabled
> those grids (the GridFile tab bug), so reload now actually loads the
> grid and hits a corrupted rewrite of the real-grid block.  OPEN,
> engine-side, best chased with a debug build on the IRIS deck.
> Candidates for that session (unverified): the `GridSrfOrder= 3` bicubic
> edge stencil overrunning `GridMat` near the grid rim; the rewritten
> `GridFile` name/path; `nGridMat` vs the file's actual dimensions; the
> pData/xData frame written for a ZrnGrData(13) element.  Note the public
> `SegDemo3data` (6 grids) round-trips clean, so the ZrnGrData save path
> is not universally broken -- the trigger is IRIS-specific.

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

---

# CLOSURE (CCMac, 2026-10-10) -- two engine defects, both fixed

The real-grid half (above, OPEN since 2026-09-09) is closed.  Reproduced
on the current engine (dev-candidate `d69b408`+, gfortran) and resolved as
**two independent defects that chain**, not the "corrupted rewrite" the
round-3 note guessed.

## Mechanism (measured end-to-end)

1. **Trigger -- `save_rx` DROPS `NSCount=`.**  The SAVE writer
   (`PrtSingleEltInfo`, `iosub.inc`) never emitted `NSCount=`.  The IRIS
   deck declares it on all 21 non-sequential segments (`NSCnt=1`, a
   per-group HIT BUDGET: one surface hit per ray in the NS group).  On
   reload the budget is 0 = UNLIMITED, so the non-sequential search
   (`tracesub.F` ~4499) never terminates the group; a ray over-searches,
   the composite-grid surface-solve bracket (`SFFZPB`, `surfsub.F`)
   diverges and the ray parameter `L` -- hence the grid pixel index
   `xi=(xData.rhom)/dAct` -- grows geometrically without bound.
2. **Crash -- INT32 overflow defeats the grid-index guard.**  The four
   grid-term sites (`SFFSrf`, `FreeFormSrf`, `SGSrf`, `NGSrf`) tested
   `i0<1 .OR. i1>nGridMat` with `i0=IDFLOOR(xi)`, `i1=i0+1`.  For
   `xi>2^31`, `IDFLOOR` saturates `i0` to `INT_MAX` and `i1=i0+1` WRAPS to
   `INT_MIN`, so the test reads all-false, the guard passes, and
   `GridMat(i0,j0)` is indexed ~2e9 out of bounds -> SIGSEGV in `INTNORM`
   (`mathsub.F:738`).  lldb stack: `INTNORM <- SFFSrf <- SFFZPSolve <-
   FreeFormSrf`.  (On Linux the same OOB read may land in mapped memory
   and return NaN instead of faulting -- same defect, platform-dependent
   symptom.)

**One-variable proof (my side; the IRIS deck/numbers stay JPL-private):**
`iris_clean_gridfile.in` (original with the GridFile tab-comment stripped,
grids active, NSCount present) traces 12737 survivors and round-trips
clean; the SAME deck with `NSCount` removed by hand crashes identically
(same SFFZPB divergence, same SIGSEGV).  So neither the grid data, the
frame, the file, nor the `lData` the writer adds is the trigger -- the
dropped `NSCount` is.

## Fix (Dave's rulings, 2026-10-10)

- **A -- round-trip `NSCount` (iosub.inc).**  `PrtSingleEltInfo` emits
  `NSCount=` whenever the element declares it (`NSCnt>0`), mirroring the
  `Link=` precedent.  **Emit in SAVE only; do NOT auto-derive** -- the
  budget is authoring intent (0/absent = UNLIMITED, which CornerCube's
  repeated bounces require; 1 = single-hit segmented groups), not derivable
  from structure, and defaulting absent->1 would alter every public NS deck.
  The misleading `tracesub.F` comment ("NSCount should NOT be defined in
  Rx") is replaced with the correct semantics.
- **B -- harden the grid-index guard (surfsub.F, all four sites).**  Reject
  a non-finite / `>=2e9` `xi`/`yj` in floating point BEFORE `IDFLOOR`; the
  ray is then a clean off-grid miss (`fh=0`) and the solve reports a normal
  bracket failure instead of crashing.  A finite in-range index is
  bit-identical to the old path (normal decks do not move).  Each rejected
  sample is counted (`nGridIdxOvf`, api `grid_idx_ovf_get`, WARN line); the
  first per run prints one note.  Rule-5 shape: warns and counts, does not
  error.

With B alone the IRIS deck no longer crashes (rays over-search, are lost
and COUNTED, trace completes) -- the host-killer is contained even when a
deck's `NSCount` was already lost.  With A the trace is also CORRECT (the
budget is restored on round-trip; 12737 survivors, NSCount preserved).

## Public gates (two focused fixtures; `ZGD_test_files/`)

The combined "`NSCount`-drop -> reload SIGSEGV" crash is NOT reproducible
in a public deck -- it needs IRIS's clocked/decentered NS tiles presenting
grazing candidates, and the public NS templates do not (segment tooling
emits sequential `Element= Segment`, which never over-searches; CornerCube
genuinely needs unlimited bounces).  Each defect is instead gated by its
own fixture, both must-fail on the pre-fix binary:

- `tst_zrngr_roundtrip.in` (+ `tst_zrngr_grid.txt`): a Cassegrain with a
  single-member NS `ZrnGrData` primary declaring `NSCount= 1` (built from
  the public `Rx_Cass_NS` geometry, flat 64x64 grid, nPass 32168).  **A1
  gate:** load -> SAVE; the saved deck must contain `NSCount= 1`.  Pre-fix
  SAVE drops it (must-fail); post-fix preserves it.
- `tst_zrngr_overflow.in`: identical but `GridSrfdx= 1.0E-10`, so
  `xi ~ 1e10` for ordinary rays.  **B gate:** pre-fix SIGSEGV (`crash_opd`
  in the CLI load gate); post-fix clean miss, `grid_idx_ovf_get > 0`,
  trace completes (nPass 32168).

**Verified gfortran, pre (`macos_prefix` @ d69b408) vs post:** A1 pre = save
lacks NSCount / post = present; B pre = crash_opd / post = ok.  **Remainder
(Linux/MATLAB side):** the ifx leg of both gates; `mmacos/tests/
tZrnGrRoundTrip.m` + the pymacos twin; the full corpus A/B (surfsub.F +
iosub.inc are shared surfaces, rule 3) -- expected 0 moved since the guard
is bit-identical for finite in-range indices and the SAVE change only adds
`NSCount=` to decks that declare it (no public deck does).

## Follow-up (resources, not this item)

An NS group containing a grid surface SHOULD declare `NSCount` -- belongs in
the segmentation emitter (`segment_rx`) as a default when it writes a
grid-surfaced NS group.  A resources-side follow-up.
