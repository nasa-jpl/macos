To: Kent (cc Jeff, Fang)
Subject: Re: Follow up on iterative ZWFS phase reconstruction

Kent — thanks, this is exactly what we needed, and the good news is that it is what our DM
gauge model already runs.

Scalar: our iterated reading is your iterative arcsin step for step — A from the off-dimple
image, the reference wave b (amplitude and phase) re-propagated from the current estimate
through the dimple's amplitude mask, the arcsin with the π/4 offset and β added, iterated.  We
write it as a cosine form so a general dimple phase and amplitude work; for the quarter-wave
dimple it is your formula, same branch.  Near a flat DM it converges in 3–4 iterations
(Δφ ~ 1e-6 rad), and the map error against the model's own field is ~1e-6 rad.

Two things we measured that your diagram does not carry, which you may find useful:

1. The fold.  On our 30 nm rms working surface (0.54 rad rms of phase) about 8 % of the lit
   pixels sit beyond the arcsin's principal branch (φ − β outside −π/4 … 3π/4).  There the
   iteration stalls — Δφ plateaus at ~0.025 rad and the map error is 0.25 rad — because the
   iteration cannot choose the branch.  We resolve it with a per-pixel branch prior from the
   phase-stepped frames of the same state; with it the error is 3e-4 rad.

2. The piston null.  What remains after the branch is a pure piston, because the intensity is
   invariant under a common phase on E and b; an iterated b cannot recover it (a modeled b
   can).  For a DM hold that is harmless; for an absolute read it has to come from elsewhere.

Vector: your atan2 of the ±π/2 pair is our vector reading in the ideal-metasurface case —
exact per pixel, no fold, so the range question becomes ordinary unwrapping as you say.  We
also carry the real metasurface's retardance error and leak and the arm's polarization
aberration per channel, which is where the vector reading's error budget actually lives on
our bench.

We will add your stopping rule (Δφ threshold rather than a fixed iteration count) so the bench
and the model stop the same way, and cite the recipe on the gauge rows.

Dave
