# Draft reply to Luis — the dw/dsurf centre-channel speckle (for Dave to send)

Luis —

We chased the speckle in the `dwdsurf` centre channel and it is not a
response: it is the finite-difference noise floor, and it is under *every*
column, not just element 4's.  Two measurements settle it.  First, a **zero**
poke — set the parameter to its nominal value, re-trace, do it again, and
divide the two traces' difference by 2·delta exactly as the Jacobian does —
reproduces element 4's column to 0.1% (rms 1.193e-06 either way, |corr| 0.9992);
and element 5's column, on the 2060 pixels *outside* its own segment footprint
where its poke provably cannot reach, carries the same 1.1916e-06.  Second, the
reason there is anything to divide: on this deck at its nominal field the
re-traced OPD is not idempotent — ten identical traces alternate in a strict
2-cycle, differing by up to 6 ulp of the 24 459 mm accumulated path on 379 of
2184 rays.  That is why the pattern is *fixed* for a given (field, zoom), which
is what made it look physical.  It is engine-side (the interactive CLI gives
the same two values bit-for-bit), it does not come from your three header lines
(the committed fixture behaves identically, and `UseChfRay4OPD= Y` actually
*helps* — it confines the noise to the 379 rays that really move, where the
aperture-mean reference smears a common shift over all 2184), and it appears
only in the five centre-field blocks: the twenty off-axis blocks are exactly
bit-reproducible and their floor is identically zero.  Element 4 is simply the
only channel whose real response never clears the floor, so it is the only one
where you see nothing else.

The fix on your side is one number: **raise the `dw_dsurf` delta from 1e-6 to
1e-4.**  The floor scales exactly as 1/delta (1.19335e-06 → 1.19335e-08, same
mantissa) while the live columns do not move — the entire difference between
the 1e-6 and 1e-4 Jacobians *is* the 1e-6 run's own floor, for all six live
channels we checked.  Two things worth knowing while you are in there: (1)
`flag_zero_norm_channels` does **not** currently drop element 4 on this rung —
its threshold is 1e-6 of the median live block and elt 4 sits at 2.6e-06, five
decades down rather than six — so its two noise columns are in your Jacobian
today; at delta 1e-4 the ratio falls to 2.6e-08 and the existing flag catches
it with no code change.  (2) `find_powered_elts` now returns elements 4–24, so
`dwdsurf` harvests 42 channels on this deck, not the 4 the committed report and
README still quote — those artifacts predate the Segment-powered ruling.

Full numbers, scripts and the scope checks are in
`macos/REPORT_sens_noise_center.md`; the resolution is also recorded in the
`zoom_5x5` README.  What alternates *inside* the trace we did not chase — that
is an open engine question for Dave and CCL, and the harvest is usable in the
meantime.
