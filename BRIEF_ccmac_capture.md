# NOTE for CCMac (Opus): on your status of 2026-09-13, and the descent

From CCL for Dave.  Follow-up to `BRIEF_ccmac_gauge_deck.md`.

1. **Accepted:** lens_deck rows (0.9911 / 0.9896 / 0.9887 on the 30 nm
   surface), the aging capture range 45 nm (the 47-site metric -- so the
   interferometer's "never fails" was the single-site ladder; state both
   numbers in the report, they answer different questions), photons
   2.5e14 at 30 nm, the report framing, parts lists, and the figure code
   in the recipe.  Commit the main-body increment when oap_deck lands.

2. **"Re-measured collapses" needs its mechanism before it goes in the
   deck.**  On the ZWFS the matrix re-measured on 60 / 90 / 120 / 160 nm
   surfaces holds the gain within 5% (runs/cap385_b60 .. b160 in
   zwfs_dm96); the four-step reads phase linearly, so a collapse there is
   unexpected.  Report the numbers per rung (gain, floor, correlation)
   and check two things: (a) that the re-measured calibration's own
   differentials and the row's differential are the WRAPPED phase
   difference (`fsdiff_`) everywhere, including inside `est_matrix_tg`'s
   calibration pokes -- an absolute map at 120 nm rms (2.4 rad rms of
   phase) wraps, a 10 nm differential does not; (b) the window-placement
   gate (D1) on the larger surfaces.  If it is the wrap, TO's unwrapper
   (below) is the fix; if it is something else, say what.

3. **Capture is a wrap problem (TO's finding, accepted).**  From a 100 nm
   rms start the differential to the set point is ~2 rad rms, so every
   reading's wrapped difference is wrapped -- the interferometer's
   four-step included.  TO is adding `dm_gauge_lib/dmg_unwrap.m` (2-D
   least-squares unwrapping on the lit mask, residue count returned) and
   a `battery.unwrap` knob, then the start-rms ladder both ways
   (`BRIEF_to_capture.md`).  When it lands: mirror it in `tg96_run`'s
   differential path (`fsdiff_` -> unwrap when `loop.start_rms` is set or
   the knob is on), gate that the record reproduces bit-for-bit with it
   off, and run the same start-rms ladder for the four-step (lens and
   OAP): start 30 / 60 / 100 / 150 / 200 / 300 nm, matrix at the start,
   recal every 10 and never, 1e13 / 1e15 photons per cycle, K 60 --
   largest converging start, cycles to 10 nm and 3 pm, residues at the
   first cycle.

4. TO pushes the shared knobs (`loop.start_rms`, `loop.recal_every`,
   `loop.intra`, `ins.recal`, `ins.measure(cmd, aux)`) now; item 5 and
   the within-scan drift unblock on origin.  Backup runs (PZT step error,
   camera walk, zwfs_run on the OAP arm) after the main-body increment,
   as you planned.
