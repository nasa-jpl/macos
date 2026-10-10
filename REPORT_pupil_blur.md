# REPORT — `pupil_blur_demo`: the gauge's numbers, and the box kernel

TO, 2026-10-09, for Dave and CC (BRIEF_to_pupil_blur.md).
- **Tool:** `MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/pupil_blur_demo.m`.
- **Record:** `runs/pupil_blur_demo/pupil_blur_demo_report.txt`, plus `_curve.png` and `_maps.png`.
- **Gate:** `tPupilBlurDemo` (SUITE_FAST, 5 tests).

Every number below is from the record.

**The floor, before and after.**
- **Before:** 8.1 % (checkerboard) / 3.6 % (random 30 nm).  This was the demo's
  estimator.  It solved for all 9216 commands against a lit-only map, which leaves a free
  ring outside the lit set, and it interpolated the map at the actuator sites.
- **Fixes:** the sites are now grid nodes (dx = pitch/4), the unknowns are the lit
  actuators only, and the blur is applied exactly as its MTF on the FFT grid.  The old
  pixel kernel was a delta below dx/4, so the leg's σ could not be represented at all.
- **After:** the σ = 0 floor is regularization plus noise, and it moves with both.
  Kernel, checkerboard: 3.51 / 0.16 / 0.062 % at λ 5e-2 / 1e-2 / 1e-3 with 20 pm, and
  0.0015 % noise-free at λ 1e-3.  The all-unknowns control stays at 5.8–7.2 % whatever λ is.
- **The record's estimator** (the measured matrix, λ_m 1e-3) has its own floor: 5.96 %
  (checkerboard) / 1.03 % (random), noise-free identical.  That is the regularization bias
  at the actuator Nyquist.  The bench reports the same roll-off in its own Stage D:
  `lensuw2` reads the (96,96) mode at gain 0.966.  It is a property of λ_m, not of blur.
- **Correction (2026-10-09, Dave via CC): that floor is one-step SHRINKAGE, not an
  accuracy limit.**  A single regularized read is a shrunk estimate; the gauge's servo
  applies it repeatedly and drives it to the plain least-squares answer (CC's iterate_
  in `pupil_blur_demo`, resources 0b2fae8).  So the 306 pm (1 %) at λ_m 1e-3 is not a
  bias in the hold rows.  It belongs in the error budget as a convergence-rate item
  (steps to settle, shrinkage per step about λ_m).  The photon floor is the 0.9 pm at
  1e14 photons from `pupil_blur_lam_m`.

**The legs' σ, and the cross-check against the deck's 0.13 / 0.29 %.**
- **The legs' σ:** the pupilsim redo records' stage-2 Nyquist gain (min over u, v) is
  0.9999 (lens) and 0.9994 (mirror).  The matching Gaussian 1/e radius is **0.0064 pitch
  (lens) and 0.0156 pitch (mirror)**.  CC's 0.011 is the same mirror number in the
  standard-deviation form.
- **The blur costs 0.0040 % / 0.024 % of the random 30 nm surface, uncalibrated.**  This
  is tg96_pupilsim's own metric (map − true, lit pupil, piston and tilt removed).
- **That is 3 % / 8 % of the deck's share: a FINDING, more than 3× apart.**  The
  pupil-imaging share is not blur.  A Gaussian would need 0.036 / 0.054 mm to reach
  0.13 / 0.29 %.  That is the deck's own stated "blur 0.05–0.14 mm", which the records'
  gain lines do not support.
- **What else is in the share:**
  - The records' recovered pokes are shifted by up to 0.0004 mm (lens) and 0.0001 mm
    (mirror).  A registration shift that size costs 0.047 % / 0.012 %: a third of the
    lens share.
  - The mirror's largest error is in the lowest band (0–0.06 cyc/mm, 0.50 % of the surface
    there), which no blur produces.
  - The rest is not identified here.
- **The calibrated reads at the legs' σ** sit within 0.1 % of their σ = 0 floors (gated).
  At the legs' σ, blur costs a calibrated read nothing measurable.

**The box-kernel line (COPHI's resolution half).**  One reading per detector cell: the cell
average, sampled on the cell grid, 20 pm per element.  The calibrated read uses the
cell-averaged matrix.
- A 0.5- or 1-pitch cell keeps the checkerboard (9.8 / 9.0 %) and costs 1.5 / 1.6 % on
  the 30 nm surface.
- A 1.5-pitch cell (0.58 readings per lit actuator) loses it (98 %) and costs 72 %.
- 2- and 3-pitch cells cost 86 / 93 %.
- So a photodiode array needs at least one element per actuator.  Between 1 and 1.5 pitch
  the read collapses, because there are then fewer readings than lit actuators.

**What is approximate.**
- The blur acts on the phase map: the small-phase limit (30 nm ≈ 0.1 wave).  The camera
  blurs intensity.
- The legs are a cos φ gain with sin φ cross-talk and a distortion.  They are matched to a
  Gaussian at the Nyquist MTF only.
- The lit set is 2852 actuators (0.85 × 0.74 R_ap), against the benches' 5072 / 6948.
- The box cells align with the actuator boundaries, and the box's naive read is the cell
  centre's point sample.
- **Two gates deviate from the brief**, because the record's λ_m sets a 1 % floor on the
  random surface (see the gate's header):
  - the "< 0.1 % at σ = 0" gate runs noise-free at λ 1e-3;
  - "calibrated < 1 % at the leg" is gated as "the leg adds < 0.1 %".

## Lines to paste

- **For `REPORT_gauge_ifo.md`:** pupil_blur_demo (2026-10-09, `runs/pupil_blur_demo`)
  puts the built legs at σ = 0.006 / 0.016 pitch.  It measures the blur at 0.004 % /
  0.024 % of the 30 nm surface, uncalibrated: 3–8 % of the 0.13 / 0.29 % pupil-imaging
  share.  So the share is mostly not blur: a third of the lens share is a registration
  residual of the records' size, and the mirror's sits in the lowest spatial band.  At the
  legs' σ a calibrated read is unchanged to < 0.1 %.
- **For the composite deck's pupil-imaging backup ("When blur is a concern, by eye"):**
  Blur only rolls off the DM's highest spatial frequencies.  The built legs sit at a
  hundredth of a pitch, where it costs a few hundredths of a percent before calibration
  and nothing measurable after.  An uncalibrated read starts to suffer at 0.1–0.2 pitch; a
  calibrated one holds to about 0.4 pitch.  The 0.13 / 0.29 % we quote for pupil imaging
  is therefore mostly not blur: it is registration and low-order error.
