# Reply to Luis (run_sensitivities 'orient' + the single-segment residual) -- FINAL 2026-09-09

Luis,

Both items are real, and both are fixed on dev-candidate for your next
pull.

1. run_sensitivities did not forward 'orient' (or 'sign') to the dw_d*
   tools, so no run_dwd* user could choose the OPD orientation.  It now
   takes 'orient' ('raw' | 'xy'), 'sign' ('opl' | 'wavefront') and
   'opd_ref' (see 2), passes them to all four channels, and prints the
   conventions in use at the top of the report.  It also now forwards
   'elts' to the dwdsurf channel (it did for the other three, not that
   one), so a one-segment Kr/Kc request stays one segment.

2. The residual on the other segments is the OPD reference, not the
   poke.  Every map is referenced to the whole-aperture mean path length
   (the engine default), so poking one segment shifts that mean and
   every other segment reads the same constant, -(N_k/N) times the mean
   of the poked segment's response.  Measured on your jwst_ote_designc
   deck through run_sensitivities (Seg2 poked, Kr and Kc, orient xy,
   PTT removal off, centre field): the other 17 segments carry 4.4e-4
   (Kr) and 1.7e-2 (Kc) per unit parameter, 5.4% and 4.5% of the poked
   segment's rms, identical on all 17 to 1e-15.  Under the chief-ray
   reference they read exactly zero.  That reference has existed since
   your August fix (macos.opd_ref), but the sensitivity tools never set
   it and a reload resets it; they all take 'opd_ref' now and re-apply
   it after every reload.  So:

       run_sensitivities(RX, ..., 'orient','xy', 'opd_ref','chief')

   'surf_remove_ptt' is not a substitute: it fits a global piston/tip/
   tilt to the whole column, which the poked segment biases, so the
   other segments trade the flat offset for a tilt.

   On the JWST deck every real segment is fully local under the chief
   reference, because the chief ray sits on the virtual centre segment.
   The one case it does not localise is a deck whose chief ray sits on
   a real segment (e5hex1's centre segment): poking that segment moves
   the reference, and the others then read about 5% of that segment's
   rms.  The clean fix for every column is a fixed nominal reference
   length, an engine branch that already exists and needs one api
   wrapper; I have the call on that and on flipping the default.

3. The pictures you sent are a third thing: under 'orient','xy' the
   per-element centre-field page (the *_pages/*_center.png files)
   rebuilt its pixel index from the transposed nominal map while the
   Jacobian rows stayed in the raw order, so one segment's poke smeared
   into diagonal streaks.  The raw-orientation page was always clean.
   Fixed (sensitivities/per_field_indx.m); the xy page is now the raw
   page transposed, exactly.

The tests that pin all three (tOpdRef, tRunSensitivities) fail on the
old code and pass now.  Slides showing the Seg2 column under each option
are attached.

Dave

(attach: demo_session/deck_dwdsurf_options.pptx)
