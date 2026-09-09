# Draft reply to Luis -- the dw/dsurf centre-channel speckle (for Dave to send)

Luis --

We chased the speckle you saw in the dwdsurf centre channel.  It was not a
response; it was the finite-difference noise floor, and it has now been
removed at its source in the engine.

What it was.  A zero poke -- re-trace the nominal deck twice and divide the
difference by 2 delta exactly as the Jacobian does -- reproduced element 4's
column to 0.1 %, and element 5's column carried the same 1.19e-6 outside its
own segment footprint, where its poke cannot reach.  The reason there was
anything to divide: on this deck at its nominal field, ten identical traces
alternated in a strict 2-cycle, up to 6 ulp of the 24 459 mm accumulated
path on ~2200 rays.  That determinism is why the pattern was fixed for a
given field and zoom, which is what made it look physical.  Element 4 is the
virtual centre segment that only passes the chief-ray sliver, so it is the
one channel whose real response never clears the floor; every other channel
had the same floor under its signal.

The cause was the source-frame re-orthogonalisation the engine performs at
every grid setup.  On a frame like this deck's, that operation has no
floating-point fixed point: fed the frame it produced, it returns a
neighbour one ulp away, and the two alternate on every trace.  Every ray
launch moved with it.  The fix keeps the incoming frame when the recomputed
one differs by round-off only, and applies the same dead band to the
chief-ray aim under an object-space stop, which had been oscillating
prescription saves by 2 ulp.  Measured after: ten traces bit-identical on
this deck with the header stop, with `stop 25`, and with no stop; the
zero-poke difference is zero on zero rays; element 4's dw/dsurf column is
exactly zero; element 5's columns are unchanged.  Single traces never
change, so nothing you have computed once is affected; only repeated traces
now agree.

What to do on your side: pull macos dev-candidate at f2fd6ea or later
(engine commit 81d3308), rebuild, relink the mex, and rerun.  Keep the
dw_dsurf delta at 1e-6; we deliberately did not raise it, since a larger
step invites nonlinearity and can obscure rays.  With the floor gone,
element 4's columns are zero and the existing zero-norm flag drops them.
Two things worth knowing while you are in there: your three header lines
were not the cause (the committed fixture behaved identically), and the
committed zoom_5x5 dwdsurf artifacts are stale -- the harvest now returns
42 channels on this deck, not the 4 the report still quotes; they predate
the ruling that made Segment elements powered.

Full numbers are in macos/REPORT_sens_noise_center.md (with the engine
addendum) and the resolution is recorded in the zoom_5x5 README.

Dave
