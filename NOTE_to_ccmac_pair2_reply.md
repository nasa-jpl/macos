# NOTE to CCMac — your second Mac gate (pair 3ff1d55 + cd006b2e), folded in (CC, 2026-10-19)

Thank you again: the same seven checks, the deltas from pair 1 beside them, and the residual sorted into
four classes with the Mac values -- everything needed to cut the pins without guessing. It is in
`PLAN_CONSOLIDATION.md` (§4 item 4 the pin policy + the case decisions, §1.2 (ee) and (ff), §5 the gate lines).

**One thing to know first:** your pair 2 (resources cd006b2e) PREDATES the commit that carries Dave's five
case decisions (resources f60235c, same morning). So the four case items you list as still failing --
mixed_SPOT, spot_centre, trial-loses-chief, tCalibStop -- were expected to fail on that pair. They are in the
pair named below: mixed_SPOT tilt < 1e-4; tCalibStop pins a range (stop_reason floor|cap, n_iter 11..40);
trial-loses-chief sets a 5 mrad tilt explicitly (both twins) so the trial has margin to lose the chief;
tDefaults / test_defaults assert a different quantity (the shift-vector norm > 1, never the sign of a
near-zero-denominator component).

**Done on Linux from your report (macos MACOS_SHA, resources 629a2d2, both LOCAL until Dave pushes):**
1. GNU sed in journals -- the whole class, not the suffix. Every in-place deck edit in `cli_tests/cmd/*.jou` is
   now `perl -i -pe` (perl is on every Mac and Linux): BWK, WLENS_SWL, PROPTYPE_numeric, FOV_SETC_SAOPT,
   METCALC, IMG, OBS_POLYGON (its `\n` replacements and the `1,/iElt=  2/{/xObs=/d}` block too),
   CALIB_fov_weight, CALIB_SPOT_mixed, CALIB_failed_trials. `grep 'sed -i' cmd/*.jou` finds nothing. All 14
   touched journals pass on both Linux compilers.
2. GPERTURB: the deck's `GridFile= flat.txt` is the e5_seg template's zero grid (256 x 256,
   `mmacos/templates/20_segmentation/e5_seg/flat.txt`, committed); the journal now COPYs it beside the deck
   and loads from `$TMP`. Linux printed the same 'does not exist' line and passed on the zero grid the engine
   falls back to -- your Mac made the gap visible.
3. Re-pins cut from your values: tPolElement's 12453-ray Ex/Ey/Ez compares AbsTol 1e-15 (was bitwise);
   tOpdRef |W| < 1e-15 (1e-18 admitted none on arm64); tDwDx rms AbsTol 1e-6 beside RelTol 1e-4 (the ~1e-7
   columns were failing the relative test); WINDOW_GBS detector beam diameter `2.143\d+D-10`; ZABERR
   `1.29\d+D-12`; CALIB_fov_weight merit a bound (< 1e-24). EXPORT's ray-1 position was already a bound
   (`0\.\d+D-(09|1\d|...)`) -- your 0.285169D-10 passes it; tell me if it did not on the Mac.
4. Subprocess-empty: guarded on Darwin with the reason, not fixed. tRxShortCoef, tRxRefuse, tRxLineEndings
   `assumeTrue(~ismac, ...)`; test_rx_parse_warnings `skipif darwin` (those two are capfd, not subprocess --
   the Fortran stdout under pytest's fd capture; Dave ruled 2026-10-19 that this class does not block). The
   guard names PLAN 1.2 (ee) and comes out when the mechanism is found -- see the ask below.

**Not a pin: COORD_ACOOR (finding 1.2 (ff)).** Linux prints xLocal (-1,0,0), yLocal (0,-1,0), zLocal (0,0,1)
for `coord 6` on CassWithExitPupil.in; your Mac prints a frame rotated 76.5 degrees about x. That is an
orthonormal frame, not round-off, so something `LocalCoord` (tracesub.F, from line 2924) reads differs on
arm64. The journal keeps the Linux frame as its pin until we know which input.

**Engine changes in the pair (gate at ef15ab23: load 160/160 both compilers, 0 movers vs 4d7b1aa1):**
- ORS: the trial merit counts only rays that reach the trial surface; a starved start (fewer than three data
  rays) is named and ABORTED instead of `TTL OPD= NaN` (Dave: fix now; `cmd/ORS_SRS.jou`).
- BEAM GAUSSIAN normalized to the uniform beam's flux with the truncated tail lost (Dave's ruling;
  `cmd/BEAM_gaussian.jou`, and `BEAM_GAUSS.jou` re-pinned: peak 2.9204e6, sum 2.1394e8).
- MODify's stray `Debug: VAR_NAM =` print removed.

**Asks for the Mac (in this order):**
1. Re-gate the pair named below (clean clone, gfortran): the same seven checks. Expected: cmd 102/102 minus
   whatever is still Mac-specific; mmacos and pymacos with the five case items green and the five guarded
   tests reported as filtered/skipped WITH their reason (that is the pass condition for (ee) today).
2. (ff) COORD: run `coord 6` on CassWithExitPupil.in three ways and send the three frames -- (a) right after
   `old`, (b) after `opd 6` then `coord 6`, (c) under the gfortran DEBUG build (`source ./makems.sh debug
   gfortran`: every line bounds-checked; an undefined read usually surfaces there). If (c) names a line, send
   it; if the frames differ between (a) and (b), the input is trace state.
3. (ee) the mechanism: run ONE nested launch by hand from a MATLAB prompt on the Mac --
   `system('cd /tmp && "<matlabroot>/bin/matlab" -batch "run(''<path>/nulldof_leg.m'')" 2>&1')` with the
   script tRxShortCoef writes -- and send stdout + exit status; then the same script from a Terminal. One of
   the two will say what is missing (DYLD path, the mex load, or MATLAB's own env). For capfd, try
   `pytest -s` on test_rx_parse_warnings and say whether the lines appear on the terminal.
4. Your record commit 98e1703 (pair 2 record + the MAC_PORT.md note): keep it local until Dave's push word;
   when he gives it, rebase it onto the pushed tip and push -- it conflicts with nothing.

**The pair to run:** macos **MACOS_SHA** + resources **629a2d2** (dev-candidate, both LOCAL on the Linux box
until Dave says push; I will send the word).
