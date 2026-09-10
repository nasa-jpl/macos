# Reply to Luis (round 4) -- FINAL 2026-09-10

Luis,

All three are real and fixed on dev-candidate.

1. run_sensitivities never forwarded 'orient' or 'sign' (nor 'elts' to
   the dwdsurf channel).  It now takes 'orient', 'sign', 'opd_ref' and
   passes all of them, plus 'elts', to every channel; the report header
   names the conventions in use.

2. The flat value on the other segments is the OPD reference, not the
   poke.  Every map is referenced to the whole-aperture mean path, so
   poking one segment shifts that mean and every other segment reads
   the same constant, -(N_k/N) x mean of the poked response.  On your
   jwst_ote_designc deck (Seg2, Kr/Kc, orient xy, centre field) that is
   5.4% / 4.5% of the poked segment's rms on all 17 others; under the
   chief-ray reference they read exactly zero.  Neither is wrong -- the
   two columns are the same data differing by one constant, and PTT
   removal is a third convention (a tilt instead of a flat offset).
   The tools all take 'opd_ref' now and re-apply it after every reload:

       run_sensitivities(RX, ..., 'orient','xy', 'opd_ref','chief')

3. The streaks in your xy pictures were the per-element page plotter:
   under 'orient','xy' it rebuilt its pixel index from the transposed
   map while the rows stayed in raw order.  Fixed; the xy page is now
   the raw page transposed exactly.

Tests for all three fail on the old code and pass now.  Slides attached.

Dave

(attach: demo_session/deck_dwdsurf_options.pptx)
