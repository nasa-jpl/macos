# BRIEF for CCMac: round 3 -- confirm the two fixes on IRIS, then the remaining decks

From CCL for Dave, 2026-09-09.  Follows `BRIEF_ccmac_jpl_round2.md` (all
three of your round-2 items are closed there; this is what to do next).
Cold-start: everything needed is here.

## 0. What landed from your round-2 findings (both pushed)

- macos `dev-candidate` `c3f3738`, engine tip **`cda178e`**: your item 3
  (save_rx -> reload SIGSEGV) is FIXED.  Root cause exactly as you
  measured: the 38 `nGridMat= 99` / `GridFile= None` elements gained a
  `pData..zData` frame on SAVE that indexed an unallocated GridMat at zero
  pitch on reload.  Two guards -- the trace ignores a grid term with no
  positive pitch, and SAVE writes the frame/`GridSrfdx` only for a grid
  that has data.  Gated on 7 public real-grid decks (byte-identical) and
  2 synthetic phantoms; `macos/REPORT_iris_save_crash.md`.
- resources `dev-candidate` **`e36db7c`**: your item 2 (reset_xp=true
  empty canvas) -- the emptyOPD warning now RE-READS the engine at the
  read surface and names the cause: "the ENGINE read is NOT empty (N rays
  pass, M non-zero samples, max|w| ...): canvas/mask fault, report it"
  vs "K rays pass: vignetting".  The canvas-construction fix itself
  needs the IRIS deck (below).

Rebuild the engine at `cda178e` or later (`source ./makems.sh release
gfortran`), pull resources to `e36db7c` or later, relink the mex
(`rm src/mmacos.mexa64` first -- the runner does not relink on an
engine-only change).

## 1. Confirm item 3 on the real deck (one MATLAB, numbers only)

```matlab
nE = macos.load_rx(rx);  macos.stop(24);
macos.save_rx('iris_rt_stop.in');   nE2 = macos.load_rx('iris_rt_stop.in');
t = macos.trace(nE2);  rs = macos.get_ray_status(t.nRays);
fprintf('stopped round-trip: rays %d, lost %d\n', t.nRays, nnz(rs.status>=2));
```
and the same with no stop (`iris_rt_nostop.in`).  Expected: both reload
and trace with the ORIGINAL's ray count (12756 stopped / 12760 not); no
SIGSEGV.  Also confirm the saved deck's 38 phantom elements now carry
only `nGridMat= 99` + `GridFile= None` (no `pData..zData`, no
`GridSrfdx`), and the 6 real grids still carry theirs.  If anything
differs from that: STOP-AND-REPORT with the element numbers.

## 2. Item 2: run the instrumented supervisor on IRIS

Same reset_xp=true harvest as round 1 (stop 24, on-axis field).  The
warning now prints which case it is.  Send the warning text verbatim
(numbers only) plus, if it says "canvas/mask fault", the three values
the single-DOF front sees at that field: `size(sf.w_nom_2d)`,
`nnz(sf.w_nom_2d)`, and `numel(sf.w_nom_vec)` -- from `dw_dx` called
directly with the same options the multi driver passes (`reload_rx`
false, `exit_pupil_elt` as the multi resolved it).  That localizes
whether the scatter (`m2v`) or the read (`local_wf`) drops the map.

## 3. Then the remaining IRIS decks -- reduced protocol

Per round-2 section 4: intake card, Test-3 pupil-read check, and a
note on whether reset_xp=true produces rows (item 2's numbers wherever
it does not).  No PRE/POST pair; OLD/NEW engine axis optional.  Dave
points you at the next deck.

## 4. For the record

- If any IRIS deck carries `LensArrayIndRef=` and a later element's
  `IndRef` reads as a lenslet index after a trace, that is a separate
  pre-existing lensarr TRACE-time overrun (PLAN section 0), not your
  finding and not the merge's.
- OPTIIX stop-less run: skipped per round 2 section 1 -- stands.
