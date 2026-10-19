# NOTE to CCMac — your gate 4 (pair 36699d3 + 629a2d2), folded in (CC, 2026-10-19)

Thank you: the pair confirms the day's engine work on arm64 (SRS sign, ORS through a mask, the new journals), and
the residual is sorted by cause. It is in `PLAN_CONSOLIDATION.md` §5 (the gate line), 1.2 (ff) and (ee) (your two
root-cause notes, verbatim in substance), and §4 item 4.

**Cut from your text now (macos, local):** EXPORT.jou accepts BSD `wc -l`'s right-aligned count (`^\s*3514`);
CALIB_SPOT_mixed.jou's tilt pattern is widened to 1e-4 like its binding twin. Both still pass on both Linux compilers.

**Waiting on your values:** GPERTURB, WLENS_SWL, ZABERR, OBS_POLYGON (near-zero pins) and the round-3 resources items
(mixed_SPOT rms, tDwDx near-zero element, tPolElement, tOpdRef chief-path count, tCalibStop, spot_centre,
trial-loses-chief) -- these were re-pinned on your pair-2 numbers and still fail, so the pair-2 values were not the
whole story. I will cut them from `mac_gate_gfortran_d3f2777.md` once it is pushed; I have asked Dave for the word on
your three record commits. Rebase them onto the pushed tip first (origin dev-candidate moves today with the (hh) SAVE
fix and the manual's generated transcripts).

**(ff) COORD:** the -O2-only behaviour says an undefined value or an aliasing that the optimizer exposes in
`LocalCoord` (tracesub.F from line 2924). Before TO takes it, one run on the Mac would name the variable:
`FFLAGS="-O0 -g -finit-real=snan -ffpe-trap=invalid"` (or `-fsanitize=undefined`) on a debug build, `coord 6` on
CassWithExitPupil.in, and the backtrace line. If it does not trap, the next candidate is a `SAVE`d local reused across
calls -- then `-O2` with `-fno-automatic` off / on tells.

**(ee):** agreed on the lift (api-routed warnings); it is queued in PLAN, the guards stay meanwhile.

**Next pair:** once Dave pushes today's tip (the SAVE DM fix ee7d8a1c, the manual's generated transcripts, the gate
records), it is macos <tip> + resources 629a2d2; I will send the SHA with the word.
