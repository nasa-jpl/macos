# DRAFT reply to Luis (run_sensitivities 'orient' + the single-segment residual)

Luis,

Both items are real, and both are fixed on dev-candidate for your next
pull.

1. run_sensitivities did not forward 'orient' (or 'sign') to the
   dw_d* tools, so no run_dwd* user could choose the OPD orientation.
   It now takes 'orient' ('raw' | 'xy'), 'sign' ('opl' | 'wavefront')
   and 'opd_ref' (see 2), passes them to all four channels, and prints
   the conventions in use at the top of the report.

2. The residual on the other segments is the OPD REFERENCE, not the
   poke.  Every map is referenced to the whole-aperture mean path
   length (the engine default), so poking one segment shifts that mean
   and every other segment reads the same constant, -(N_k/N) times the
   mean of the poked segment's response.  We measured your case through
   the driver (e5hex1, segment 2, Kr and Kc, orient xy, PTT removal off):
   the other six segments carry 4.3e-4 (Kr) and 2.0e-2 (Kc) per unit
   parameter, 15% and 12% of the poked segment's rms, identical on all
   six to 3e-16.  Under the chief-ray reference they read exactly zero.
   That reference existed since your August fix (macos.opd_ref) but the
   sensitivity tools never set it and a reload resets it; they all take
   'opd_ref' now and re-apply it after every reload.  So:

       run_sensitivities(RX, ..., 'orient','xy', 'opd_ref','chief')

   'surf_remove_ptt' is not a substitute: it fits a global piston/tip/
   tilt to the whole column, which the poked segment biases.

   One case the chief reference does not localise: the chief ray's own
   segment (the centre one), because its reference moves with the poke;
   there the others read a constant 5% of that segment's rms (Kr).
   The clean fix for every column is a fixed nominal reference length
   (an engine branch that already exists, one api wrapper away); Dave
   has the call on that and on flipping the default.

The test that pins this (tOpdRef) fails on the old code and passes now.

-- Dave
