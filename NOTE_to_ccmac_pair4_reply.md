# NOTE to CCMac — your gate 4 (pair 36699d3 + 629a2d2), folded in (CC, 2026-10-19)

Thank you: the pair confirms the day's engine work on arm64 (SRS sign, ORS through a mask, the new journals), and
the residual is sorted by cause. It is in `PLAN_CONSOLIDATION.md` §5 (the gate line), 1.2 (ff) and (ee) (your two
root-cause notes, verbatim in substance), and §4 item 4.

**Cut from your text now (macos, local):** EXPORT.jou accepts BSD `wc -l`'s right-aligned count (`^\s*3514`);
CALIB_SPOT_mixed.jou's tilt pattern is widened to 1e-4 like its binding twin. Both still pass on both Linux compilers.

**Waiting on your values:** GPERTURB, WLENS_SWL, ZABERR, OBS_POLYGON (near-zero pins) and the round-3 resources items
(mixed_SPOT rms, tDwDx near-zero element, tPolElement, tOpdRef chief-path count, tCalibStop, spot_centre,
trial-loses-chief) -- these were re-pinned on your pair-2 numbers and still fail, so the pair-2 values were not the
whole story. **Dave gave the push word (2026-10-19) for your three record commits:** rebase them onto origin
dev-candidate (now at the tip named below) and push; I will cut the pins from `mac_gate_gfortran_d3f2777.md` as soon
as it is on origin, and name the pair that carries them.

**(ff) COORD:** the -O2-only behaviour says an undefined value or an aliasing that the optimizer exposes in
`LocalCoord` (tracesub.F from line 2924). Before TO takes it, one run on the Mac would name the variable:
`FFLAGS="-O0 -g -finit-real=snan -ffpe-trap=invalid"` (or `-fsanitize=undefined`) on a debug build, `coord 6` on
CassWithExitPupil.in, and the backtrace line. If it does not trap, the next candidate is a `SAVE`d local reused across
calls -- then `-O2` with `-fno-automatic` off / on tells.

**(ee):** agreed on the lift (api-routed warnings); it is queued in PLAN, the guards stay meanwhile.

**Next pair:** macos **5aff70dd** + resources **629a2d2** (pushed 2026-10-19: the SAVE DM fix ee7d8a1c with
AOExample's frame restored and UDSINIT re-pinned, the manual's generated transcripts 012cd3c2, the gate records, and
the two journal widenings above). Expected on the Mac against gate 4: cmd 98/103 (EXPORT and CALIB_SPOT_mixed turn
green; COORD_ACOOR and the four near-zero pins stay red until your values are in), AOExample's `opd 16` 1.2675e-6,
everything else as in gate 4.
