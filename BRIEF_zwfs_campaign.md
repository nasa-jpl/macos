# BRIEF: ZWFS campaign — the Zernike sensor against the PSI gauge, same DM truth

Dave 2026-09-04: "plan the ZWFS campaign."  Motivation on record (Stage
E ruling, tg96): "how well a deviation is measured — that's really the
challenge for this instrument, and it's likely what the ZWFS will do
better."  This brief is the plan; execution waits on Dave's steer at
the decision points below.

## Objective

**Ultimate target (Dave 2026-09-04): measurement error ~ 1 pm.**  For
calibration: the IFO's map-space single-poke error is 49 pm and the
sensor's 286 pm at its starved camera — both far above target.  The
road to pm class runs through the differential protocol (common
systematics cancel), actuator-space fitting (Dave's scoring ruling),
sampling that satisfies the ILLUMINATED-pupil budget, and likely
averaging / reconstructors beyond linear.  Every stage now reports
errors in pm so the distance to target stays visible.

Model a Zernike wavefront sensor measuring the SAME 96×96 Xinetics DM
truth as the TG96 polarization-PSI gauge, score it with the SAME
battery, and answer one question with a table: **does the ZWFS beat
the PSI gauge on the differential benchmark** — a 10 nm
single-actuator deviation read to 0.021 nm rms, a 10 nm rms random
deviation read to 3.67 nm (37%, instrument-transfer-limited), both
base-independent (tg96_report.txt, run 10)?

Secondary question the head-to-head sets up for free: where does each
lose?  Expected shape — ZWFS wins on sensitivity/simplicity (one frame,
no arms, no polarization train), PSI wins on dynamic range (ZWFS
response is linear only for small phase; the 30 nm working state is
safe, but somewhere it folds — find the break scale).

## What already exists (the plan is mostly assembly)

- **The mask, physically parameterized:**
  `templates/40_benches/vsg_wip/vsg2_params.m` §9 carries the real
  VSG2 ZWFS hardware — transmissive etched fused-silica substrate
  (Thorlabs W4101FT1), etch 346.2 nm → phase 2π(n−1)d/λ ≈ π/2 at
  632.8 nm (the classic quarter-wave dimple), the 9-spot table of
  dimple diameters (spot 9 = 1.06 λ/D default; 1.22, 2.0, 3.0
  alternates), a leakage parameter.  Source: "VSG2 Zernike Wavefront
  Sensor Update -v2.pptx".
- **The mask machinery:** `bench_ctb/ctb_mask_phase.m` already builds
  circular focal-plane phase masks (Roddier / dual-zone) with
  supersampled gray edges, centred on the FFT DC pixel, applied to the
  complex field via `macos.apodize_complex` — the exact family.  A
  ZWFS kind is a small extension (arbitrary phase disk: the Roddier
  form with φ = π/2 instead of π).  The CTB core-pixel supersampling
  practice is the answer to "a 1 λ/D dimple spans few focal pixels."
- **The DM truth + battery + doctrine:** tg_psi_dm96 — the same
  influence-function grid truth, the clearance-solve and
  sampling-budget stages (Dave's design rule), the two-poke
  registration doctrine (4 DOF classes; symmetric targets banned), the
  12-mode transfer battery, and the differential protocol, all as
  committed code with printed gates.
- **Pupil-reimage machinery:** the ZWFS train is pupil → focus (mask)
  → pupil (detector); the CTB chain exercises exactly that leg
  (pupil→FPM→Lyot), and the Bench builder has add_pupil/relay pieces.

## What is genuinely new

1. **The sensing train:** a Bench — collimated 96 mm beam onto the DM,
   focusing leg to the mask plane, reimaging leg to a pupil image on
   the detector.  Laid out with the SAME Stage A/A2 discipline:
   clearances with real bodies and printed margins; sampling budget
   asserting (a) detector ≥2× actuator Nyquist on the reimaged pupil
   (385-px class, as TG96), and (b) the NEW interface — focal-plane
   grid resolution at the dimple (grid pitch vs 1.06 λ/D, gray-edge
   supersampling per CTB practice).
2. **The reconstructor:** intensity → phase.  R1 = small-phase linear
   inversion about the model-computed reference wave b(x) (b is exact
   here: trace the flat-DM system with and without the mask).  R2 =
   the N'Diaye exact quadratic inversion (extends dynamic range; the
   comparison of R1/R2 IS part of the dynamic-range answer).  Height
   convention pinned in S1: single reflection off the DM → h = φλ/4π,
   same factor-2 as PSI.
3. **The battery harness adapted:** no arms, no analyzer sweep, no
   four-step — ONE frame per measurement (plus the stored reference
   frames).  The scoring code is reused unchanged.

## Stages (each gated, numbers printed, report is the record)

- **S1 — mask + response gate.**  Build the dimple mask (vsg2_params
  numbers); verify the single-frame intensity against the analytic
  ZWFS response for known nm-scale tilt/defocus; pin the sign and
  height conventions.  Gate: measured response = b(x)-based analytic
  prediction; a deliberately wrong-sign reconstruction FAILS (non-
  vacuity).
- **S2 — bench + registration.**  Lay out the train at 96 mm
  (clearance + sampling stages, margins printed); re-run the two-poke
  registration for this camera (doctrine says the flip/transpose and
  sign are DECK-dependent — never inherit from TG96).  Gate: |corr| ≥
  0.8, runner-up separation ≥ 0.3, as tg96.
- **S3 — battery.**  Null, piston gain, single actuator at 150 nm,
  the same 12-mode transfer curve, held-out random.  Deliverable: the
  side-by-side table vs run-10 PSI (same modes, same truth files).
- **S4 — the differential head-to-head.**  The same four rows
  (flat / 30 nm working × single-actuator / random 10 nm), plus the
  dynamic-range axis: grow the working state until the ZWFS breaks;
  report the break scale and what R2 buys over R1.  This is the
  campaign's headline table.
- **S5 — trades (steer-dependent).**  Dimple diameter across the real
  spot table; etch-phase error; leakage; chromaticity (the etch is
  fixed glass — phase slides with λ, machinery exists in the dual-zone
  kind).  Only if the head-to-head motivates them.

## Scoring ruling (Dave 2026-09-04, applies to BOTH instruments)

**Score by how well DM STATE is recovered, in actuator space** — not
by map-pixel residuals.  Either (a) model-based sensing: the estimator
fits actuator commands directly through the DM actuation model
(influence functions), or (b) fit the DM model to the recovered
higher-resolution wavefront — then score by how well actuator CHANGES
are measured.  Consequences:
- The TG96 run-10 battery/differential numbers are MAP-space; a
  rescore (Stage E′, actuator-space) is queued so the PSI benchmark
  is stated in the same currency before the head-to-head.  The
  actuator fit absorbs instrument roll-off inside the DM's own band —
  this supersedes the "apply measured transfer first" note.
- ZWFS S3/S4 score in actuator space from the start.  The fit uses
  the same influence-function forward model that BUILDS the truth
  grids, applied through the two-poke registration affine.

## Multi-mask phase stepping (Dave 2026-09-04: "masks with varying depth")

Frames through dimples of several etch DEPTHS (φ = π/2, π, 3π/2 at one
diameter) plus the clear frame solve the per-pixel field EXACTLY
(linear 3×3: I_k = I₀ + |c_k|²·|E_b|² + 2Re(c_k·E_b·conj(E₀))) — no
small-phase assumption, no sign ambiguity, per-pixel range ±π and
unwrappable.  The sequential cousin of the phase-2 metasurface (which
delivers two phases simultaneously).  Cost: 4 frames per measurement
vs 1 (still no polarization train; vs the IFO's 6 traces).  Hardware
implication: a mask substrate carrying spots of several DEPTHS — the
VSG2 part varies diameter at one depth.  Stage S2b (zwfs_s2b.m);
sensitivity rerun with the stepped retrieval follows it.

## Sensitivity stage (Dave 2026-09-04)

How small a change is detectable, and how accurately: single pokes
and grid pokes, against flat and a ~30 nm background, amplitudes
10 nm → 0.1 pm, differential protocol, actuator space; detection =
SNR ≥ 5 over the unpoked-actuator floor.  Objective: sensitivity at
1 pm or below.  Noiseless model → the measured floor is systematic +
numerics; photon noise is a later budgeted stage.  Twin scripts
tg96_sens.m / zwfs_sens.m.

## Multi-color stage (Dave 2026-09-08, evening: "try running both systems
   at multiple colors, maybe the combination will help with some poor SNR
   regions") — RUN, both instruments

`zwfs_s6color.m` / `tg96_s6color.m` + `dm_gauge_lib/dmg_color_comb.m`
(multi-channel Wiener on the actuator lattice, a_hat = Σ_k G_k A_k /
(Σ_k G_k² + β²)).  Five colors 480/532/632.8/700/780 nm; only the deck
header `Wavelen=` changes; every calibration redone per color; 96×96;
noiseless.  **ZWFS: YES.**  The oscillatory transfer null (near 30
cyc/ap at 632.8) migrates as 1/λ with the dimple's angular size (40 /
36 / 24 / 20 cyc/ap at 480 / 532 / 700 / 780), so the five-color
transfer never drops below 0.99 (best single 0.91).  Rows, linear
reading, 632.8 → comb: hold-out 202 → 143 pm; dense random 16.4 → 7.2
nm; single-on-30nm-base floor 720 → 224 pm (SNR 9 → 35); the one
UNDETECTED scenario (grid on base) SNR 1.46 → 3.43 — better, still under
5.  Stepped: base rows the same way; the flat hold-out slightly WORSE
(287 → 313) because the stepped systematic is larger at the other
colors and the equal-weight combiner inherits it — weight by a measured
per-color systematic before combining stepped readings.  **IFO: NO
lever** — transfer and every row identical across colors to 3 digits
(92 pm / 4.17 nm / SNR 210 / 19.8 at each λ and combined).  A
diffraction roll-off would have moved 1.6× between 480 and 780; it moved
< 0.3%, so the IFO's high-f deficit is GEOMETRIC (the null-tuned tail's
conjugate + 0.136 mm distortion) — the joint tail objective is its
lever, not the source.  Cost: K× frames; at equal total light the noise
part is neutral, the systematic part improves.  Records: both READMEs
(S6 sections), `zwfs_s6color_report.txt` / `tg96_s6color_report.txt`,
figure `zwfs_dm96/zwfs_s6color.png`.

## Decision points — RULED (Dave 2026-09-04)

1. **Scale:** 96×96 rig; may mask down to 16×16 (1 mm actuators) to
   speed development; **real work at 48×48 AND 96×96** — the battery
   and differential tables run at both.
2. **Mask form:** single dimple first.  **Phase 2 (after the scalar
   system is built and tested): a polarizing METASURFACE producing
   TWO separate phase images** — the vector-ZWFS form (opposite
   dimple phase per polarization → two simultaneous pupil images,
   phase-diverse; kills the sign ambiguity and extends range).
   Engine precedent: the vector-diffraction 3-plane chain + the CTB
   vector-vortex per-plane mask machinery (ctb_mask_vvc).
3. **Optics class:** lens train first.
4. **Location:** `templates/40_benches/zwfs_dm96` — approved.

## Cost estimate

S1–S3 ≈ one working session (the registration saga is doctrine now,
not discovery); S4 is cheap once S3 stands (frames are single traces).
Model 1024 runs at tg96-demonstrated runtimes (battery ~5 min-class
per stage on this box).

## Records

Campaign dir README per the pattern; every gate prints; failure
reports preserved; this brief's resolutions written back at resolution
time.  Fold into the Fang deck only on Dave's ask (the deck already
names the ZWFS comparison as planned).
