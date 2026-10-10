# BRIEF — `pupil_blur_demo`: make its numbers the gauge's, then add the box kernel

For **TO**, written 2026-10-09 by CC for Dave.  You start from scratch; this
brief is self-contained.  Work in `~/dev/MACOS_resources`, `dev-candidate`,
tip **7a1edb2 or later** (CCMac's demo), with `~/dev/macos` at `a1e84c4` or
later.  Commit locally with SHA + branch stated; push only on Dave's word.
Do not relink the shared mex or rebuild a shared tree while another MATLAB
runs (this task needs no engine, no mex — pure MATLAB on `dm_gauge_lib`).

## Why this exists

Dave gave a DM-surface-gauge talk on 2026-10-08 (deck
`macos/demo_session/deck_gauges_short.md`, FINAL; the long form
`deck_gauges.md`).  Four gauges (Twyman-Green interferometer, scalar and vector
Zernike sensors, point-diffraction sensor) read a 96 × 96-actuator DM to
picometers on one modeled bench, two front ends (lens rig, off-axis-mirror
rig).  Fang Shi is skeptical of the deck's pupil-imaging analysis (the DM is
imaged onto the camera; the claim is that the leg's blur and distortion cost
0.13 % (lens) / 0.29 % (mirror) of the raw 30 nm map after calibration, from
`tg96_pupilsim`, an engine-traced, three-stage tool).  Dave wanted a version a
skeptic can check by eye.  CCMac wrote it: `templates/40_benches/
tg_psi_dm96_oap/pupil_blur_demo.m` (engine-free, ~40 s): build a known DM
surface from known commands (Gaussian influence functions), blur it with a
Gaussian kernel of swept 1/e radius σ, add read noise, recover the commands
with `dmg_stencil` + `dmg_act_fit` naive (blur ignored) and calibrated (blur
in the stencil), on a ±50 nm checkerboard (actuator Nyquist) and a random
30 nm surface; plot command error vs σ.  Records in `runs/pupil_blur_demo/`;
CCMac's hand-off `macos/NOTE_to_ccl_pupil_blur_demo.md`.

CC reviewed it (`macos/NOTE_to_ccmac_pupil_blur_review.md`, read it in full —
it holds the measurements).  The shape is right; the numbers are not yet the
gauge's.  CCMac's budget is spent; you pick it up.

## What is wrong, measured

1. The "blur-free floor" (σ = 0: checker 8.1 %, random 3.6 %) is the demo's
   estimator, not the gauge.  λ swept 0.05 → 1e-4 with and without noise:
   8.1 → 6.5 %, then flat, noise-free identical.  Grid refined 4×: 5.9 %.
   Cause: `dmg_act_fit` solves for ALL 9216 commands while the map is masked
   to the lit set; the ring outside the lit set is free and takes 13.6 nm rms
   (36 nm max) of command where the truth is 0.  Lit-only unknowns: 1.80 %.
   The remaining 1.8 % is `interp2 'linear'` at the actuator sites (sites sit
   between grid points).  The gauge-relevant floor is ~0.
2. "The gauges' own reconstruction" is the bench's `calib_mode = 'kernel'`
   (the S3 flavor).  Every deck number is `'matrix'`: the measured response
   matrix in detector pixels (`calib_matrix_` / `est_matrix_tg` in
   `tg96_run.m` ~565–740, `calib_matrix_` / `est_matrix_` in `zwfs_run.m`
   ~966–1091), one column per lit actuator, `P.battery.matrix_lam = 1e-3` of
   the median column energy, lit unknowns only.  Its floors are 1.4 / 2.3 pm on
   a 10 nm change (0.01–0.02 %).
3. The real leg is σ ≈ 0.011 pitch, not "~0": `runs/pupilsim_redo_oap/
   pupilsim_redo_oap_report.txt` gives the phase gain at the actuator Nyquist
   (0.5 cyc/mm) 0.9994 worst / 0.9999 mean, poke width 0.999; a Gaussian with
   MTF(f_Nyq) = 0.9994 has σ = sqrt(−ln 0.9994 / (2π² f_Nyq²)) = 0.011 mm.
   The demo must show that point, and its NAIVE error there is the number
   that has to agree with the deck's 0.13 / 0.29 %.
4. The demo's lit set (2852, 0.85 × 0.74 R_ap) is not the bench's (5072 lens
   / 6948 mirror, from the traced illumination).

## The task, in order (each step its own commit)

**Step 1 — the floor.**  Solve for lit unknowns only and sample exactly at the
actuator sites (the demo owns the model: evaluate the Gaussians at the sites
analytically, or put the sites on the grid).  Keep `dmg_act_fit` itself as it
is unless you find the bench's kernel mode wants the same fix — if so, that is
a library change: gate it (`tDmgLoop` must stay 15/15) and say so.  Expect
σ = 0 error at the noise term (20 pm / √(pixels per actuator) ≈ 0.02 % of
30 nm).  Re-run the λ × noise sweep (0.05, 1e-2, 1e-3, 1e-4 × 20 / 0 pm) and
print it in the report: the floor must now move with λ and noise, which is
the claim the demo's text makes.

**Step 2 — the record's estimator as the calibrated curve.**  Add
`'estimator','matrix'`: columns = the blurred unit influence at each lit site
(the measured matrix, which is why calibration absorbs blur automatically),
λ = 1e-3 of the median column energy, the same solve as `est_matrix_tg`.  Keep
`'kernel'` as the option.  Both curves in the figure; the report says which
is the deck's.

**Step 3 — the leg on the axis.**  Read σ for both rigs from the pupilsim
records (`pupilsim_redo_lens`, `pupilsim_redo_oap`: the Nyquist gain line;
convert as above), mark them on the curve, and print the naive and
calibrated errors at those σ.  **The cross-check that closes Fang's
question:** the naive error at the leg's σ on the random 30 nm surface,
against the deck's raw-map pupil-imaging share 0.13 % (lens) / 0.29 %
(mirror) — same order, or explain the difference (the leg is a cos φ gain
with amplitude cross-talk, not a Gaussian; the equivalence is at the Nyquist
MTF only; say that in the caption).  A disagreement of more than ~3× is a
finding, not a rounding — report it, do not tune to it.

**Step 4 — the geometry and the caption.**  Use the bench's lit rule
(`dmg_lit`, the traced cone) or state 2852; state that the blur is applied to
the phase map, valid in the small-phase limit (0.1 wave here) — the camera
blurs intensity.

**Step 5 — the box kernel (for the COPHI campaign).**  `'kernel','box'` with
the cell size swept in pitch units (0.5, 1, 1.5, 2, 3): the cell average of a
detector element.  Same curves.  This is the resolution half of the
speed-vs-resolution trade for Feng Zhao's photodiode-array COPHI (see
`macos/NOTE_cophi_first_look.md` §3) — the result is a line "a cell of N
pitches loses the checkerboard / costs X % on the 30 nm surface, calibrated."

**Gate.**  `mmacos/tests/tPupilBlurDemo.m` (SUITE_FAST if under ~60 s; else
its own class): σ = 0 error < 0.1 % on both patterns; the calibrated error at
the leg's σ < 1 % (no fitted value); the naive checker error at σ = 0.4 pitch
> 10× its σ = 0 value (so the demo is not vacuous); the must-fail leg = the
all-unknowns solve, asserted to give > 5 % at σ = 0 (the defect, kept as the
negative control).

## Deliverables

- The tool, its README entry (tg_psi_dm96_oap/README.md "Fang Shi's question"
  block and the Files table — CCMac's lines, updated), `runs/pupil_blur_demo/`
  regenerated (report + the two figures; the λ sweep and the box-kernel rows
  in the report).
- `macos/REPORT_pupil_blur.md`: half a page — the floor before/after, the
  leg's σ per rig and the cross-check against 0.13 / 0.29 %, the box-kernel
  line, what is approximate (phase-map blur, Gaussian vs cos φ).
- One line each for `REPORT_gauge_ifo.md` (the lane report) and for the
  composite deck's pupil-imaging backup, ready to paste.

## Rules (the ones that have bitten this campaign)

- Numbers in any report come from a committed run record (run it, commit the
  record, quote the record).
- Figures are the tool's own output; no re-rendering.
- Guards warn, not error; nothing in `dm_gauge_lib` changes without a gate.
- American English; "calibration removes systematic errors" is Dave's wording.
- Do not estimate the leg's σ by eye from a PSF plot; take it from the gain
  line in the pupilsim record as above.
