# NOTE for CCMac (Opus): the mask gauges on the OAP rig -- the focuser is not diffraction-limited at the seat

From CCL for Dave, 2026-09-14.  Answers your question: not an issue
tracker; this note plus your report line.

1. **Item 3 accepted:** the four-step is robust to all three sequential
   systematics in hold (step error common-mode; camera walk 13% at 25%
   within-scan; DM within-scan 1.6e13 vs 1.8e13).  Your conclusion goes
   in the deck as stated: the hybrid's value is the absolute
   calibration, not drift immunity.  OAP descent as backup: fine.

2. **The mask gauges on the OAP rig -- what the ZWFS runner finds.**
   The runner is ready (bench.coat_oap wired and pushed, a360faf;
   bench.MASK_TRIM passes through to twyman_green).  A seat-trim scan on
   your OAP rig (zwfs_params + twyman_green 'optics' 'oap', model 512,
   65 rays, the flat DM), ray spot rms at the FocalMask and the focal
   field's peak/sum (the gauge's gate is >= 0.01; the lens rig passes):
   trim -5.58 (the lens value): 0.80 mm / 0.0000; 0: 0.42 mm / 0.0001;
   +5: 79 um / 0.0000; +6.00: 10.6 um / 0.0013; **+6.25: 8.0 um /
   0.0021**; +6.5: 24 um / 0.0002; +8: 127 um / 0.0000.  The true focus
   is ~+6.25 mm from the seed (the lens rig's -5.58 corrects a thick
   lens's thin-lens seed; the parabola has a different seed error), and
   at best focus the residual ray blur is 8 um rms = 3 lambda F/D
   (lambda F/D is 2.6 um), with the focal field five times less
   concentrated than the gate wants.  The focuser as emitted is not
   diffraction-limited at the mask seat.  The interferometer never
   needed it (its seat is empty; the DM is imaged), so it never showed.

3. **What to check, on your side of the builder** (the 'oap' path of
   twyman_green / add_focuser and the rig's OAP2 geometry): OAP2's conic
   and off-axis distance against the collimated return beam's actual
   direction and offset into it (a 3 lambda F/D blur is coma or astig of
   an OAP fed off its design axis -- the D4 sensitivities are the
   handle); the FocalMask element's orientation along the reflected
   chief; then the seat trim solved on the trace (the ZWFS S1 recipe:
   find the ray focus, then the diffraction focus, MASK_TRIM as the
   knob), and the mask sandwich spheres centered on it.  Gate: the flat
   DM's focal spot at the seat with peak/sum >= 0.01 and the ray blur
   under 1 lambda F/D; then `zwfs_run` bench stage G1 / G3 pass and the
   item-4 rows follow.  If the OAP rig cannot make a clean focus at the
   seat by design, say so with the number -- that is itself a deck
   statement about the OAP front end under the mask gauges.
