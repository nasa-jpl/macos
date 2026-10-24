# Review of `pupil_blur_demo` (CCMac, resources 7a1edb2) — CC, 2026-10-09

For CCMac and Dave.  Short verdict first, then the measurements behind it, then
what I would change before the curve goes near a deck.

## Verdict

The demo's SHAPE is right and it is the thing Fang asked for: a known surface,
a swept blur, the actuator commands recovered naive vs calibrated, the lesson
that blur is an MTF roll-off a calibrated read removes until MTF at the
actuator Nyquist nears the reconstruction's floor.  Keep it.

Its NUMBERS are not yet the gauge's, in three ways that a skeptic would find:

1. **The "blur-free floor" (8.1 % checker / 3.6 % random) is the demo's own
   estimator, not regularization and not the gauge.**  I swept λ from 0.05 to
   1e-4 with and without the 20 pm noise: the floor moves 8.1 → 6.5 % (checker)
   and 3.6 → 3.57 % (random) and then stops; noise-free is identical to the
   digit.  Refining the map grid 2× and 4× (dx 0.28 → 0.14 → 0.07 mm) moves it
   6.5 → 6.0 → 5.9 %.  What it is: `dmg_act_fit` solves for ALL 96 × 96
   commands while the map is masked to the lit set, so the ring of actuators
   just OUTSIDE the lit set is free and takes 13.6 nm rms (36 nm max) of
   command on the checkerboard where the truth is 0, and the lit edge ring
   carries the complement.  Restricting the error to lit actuators more than
   2 pitches inside the edge: 1.8 %; solving for LIT unknowns only (the form
   the bench's record estimator uses): 1.80 %.  The remaining 1.8 % is the
   linear interpolation of the map at the actuator sites (the sites sit
   between grid points; the checkerboard's curvature there differs from the
   single poke the stencil was sampled from).  So the gauge-relevant floor of
   this demo is ~0, and the curve's left end should say so.
2. **"The gauges' OWN reconstruction" is the kernel mode, not the record's.**
   `dmg_stencil` + `dmg_act_fit` is `P.battery.calib_mode = 'kernel'` (the S3
   flavor, λ 0.05 relative to the stencil peak).  Every number in the deck is
   `'matrix'`: the measured response matrix in DETECTOR pixels (`calib_matrix_`
   / `est_matrix_`, one column per lit actuator, λ = 1e-3 of the median column
   energy, lit unknowns only).  Its floors are the deck's 1.4 / 2.3 pm on a
   10 nm change — 0.01–0.02 %, 300× below the demo's floor.  In the matrix
   form "calibrated" is automatic (the blurred columns ARE the matrix), which
   is the cleanest way to say to Fang why calibration absorbs the blur.
3. **The real leg's point on the axis is σ ≈ 0.011 pitch, not "~0".**
   `runs/pupilsim_redo_oap/pupilsim_redo_oap_report.txt`: the leg's phase gain
   at the actuator Nyquist (0.5 cyc/mm) is 0.9994 worst, 0.9999 mean; pokes
   recover width 0.999.  A Gaussian kernel with MTF(Nyquist) = 0.9994 has
   σ = sqrt(−ln 0.9994 / (2π² f²)) = 0.011 mm = 0.011 pitch.  Mark it.  (The
   leg's "blur" is a phase gain cos φ with amplitude cross-talk sin φ, not a
   Gaussian; the equivalence is at the Nyquist MTF only — say so in the
   caption.)  With the floor fixed, the NAIVE curve at σ = 0.011 pitch is the
   number that must agree with the deck's raw-map pupil-imaging share
   (0.13 % lens / 0.29 % mirror) — that is the cross-check that ties the demo
   to the record, and it cannot be made while the floor is 3–8 %.

Smaller points:

- The demo's lit set is 2852 actuators (0.85 × 0.74 R_ap); the bench's are
  5072 (lens) and 6948 (mirror) from the traced illumination.  Either use the
  bench's rule or state the demo's.
- The blur is applied to the PHASE map.  That is the small-phase limit (the
  30 nm surface is 0.1 wave of wavefront): the camera blurs intensity, and
  for the four-step the phase recovered from blurred fringes equals the
  blurred phase only when the fringe is near null.  True here; say it, since
  Fang's question is about the camera.
- `pcg` at 400 iterations past σ ≈ 0.8 pitch: immaterial, as CCMac says.

## Where it stands with respect to the gauge work

- It changes no number in the deck of record.  The deck's pupil-imaging
  numbers come from `tg96_pupilsim` and the raw-map decomposition; the demo
  is the plain-language companion and, once the floor is fixed and the leg's
  σ is on the axis, a backup slide for the composite deck ("When blur is a
  concern, by eye") beside the pupil-imaging slide.
- It is also the seed of a COPHI-campaign tool: the photodiode-array form of
  COPHI is a BOX kernel of one detector cell (the cell average) at the
  actuator pitch or coarser — add `'kernel','box'` with the cell size swept
  and the same curves answer the speed-vs-resolution trade's resolution half.
- Per-candidate error budgets (follow-up 1): the blur row is one cell; the
  demo's calibrated-minus-naive difference at the leg's σ is the "cost of not
  calibrating," which is the number a budget wants.

## Suggested edits, in order

1. Solve for lit unknowns only (mask the unknowns as well as the map), and
   sample the map at the actuator sites exactly (build the surface and the
   unit influence on a grid that contains the sites, or evaluate the
   Gaussians at the sites analytically — the demo owns the model, so it can).
   Expect the σ = 0 error to fall to the noise term (20 pm / √(pixels per
   actuator) ≈ 0.02 % of 30 nm).  Re-run the λ sweep to show the floor is
   now λ's and noise's, which is the claim the text makes.
2. Add the leg's σ = 0.011 pitch as a labelled point, from the pupilsim
   record, and quote the naive error there against the deck's 0.13 / 0.29 %.
3. Optionally the matrix estimator (blurred columns, λ 1e-3 of column energy)
   as the "calibrated" curve, so the demo's calibrated read IS the record's.
4. Use the bench's lit count or state 2852.

Scripts I ran (scratchpad, not committed; the numbers above are the record of
them): the λ × noise sweep at σ = 0/0.1/0.2/0.4 pitch, the grid refinement at
λ = 1e-3, and the edge decomposition with a lit-only `lsqr` solve.  ~40 s a run.
