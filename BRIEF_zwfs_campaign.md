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

## S7 — model correction + the iterated-reference reading (2026-09-09, RUN)

Executing the literature scan's first task (`REPORT_zwfs_lit_scan.md`)
exposed a MODEL DEFECT first: the `twyman_green` 'nf' sandwich emitted
the exit reference sphere with zElt/Kr = 0.6·D_MASK_FL (23.86 mm)
against the entrance sphere's 352.7 mm.  The engine's SPH2PL leg applies
a focal quadratic factor S ∝ (Z2−Z1)·Z1/Z2 and PL2SPH is a plain FFT, so
the unmasked round trip was a Fresnel DEFOCUS of the reimaged pupil by
z_eff = 4.86 m (entrance-sphere scale), not the identity the ctb_dcr.in
precedent gets with equal radii.  Measured on the legacy deck: round
trip 0.159 with a FLAT DM; 29% rms detector amplitude modulation under
the 30 nm state (a phase-only state must give 0); the ringed poke
kernel (raw peak 0.27); the oscillatory transfer null near 30 cyc/ap =
a Talbot null (predicted 34).  **Every ZWFS number in S1–S6 was taken
on that defocused sensor; the IFO twin is all-geometric and untouched.**
Fix: 'nf' now emits the SYMMETRIC sandwich (round trip 1.8e-15);
'nf_legacy' reproduces the old emission byte-for-byte (tBench gate);
S1–S6 stand as the legacy-model record.

The reading itself (`dmg_zwfs_gauge` measI/reconI): per-pixel exact
solve with the reference wave re-propagated through the FFT surrogate
of the mask model (validated against the engine's own Eb at 2e-15), one
frame.  Gated: with the oracle b and the true branch the solve is exact
to 3e-14 on every pixel; its two residuals are the quarter-wave
sensor's per-pixel BRANCH (7.8% of pixels beyond the fold on the 30 nm
base) and PISTON (the intensity is invariant under a common phase on E
and b — the sensor's piston null).  'I+' = the same one frame plus a
branch prior from a ONE-TIME stepped retrieval of the working state,
REFINED by re-solving that retrieval with the iterated |b|² (two passes
reach the true branch on 99.99% of pixels; the plain stepped prior
misses 3%, enough to sign-flip a single-actuator differential whose
footprint sits beyond the fold — the 48×48 case).  Map space, piston
removed: the 30 nm working state (0.54 rad rms) is read from ONE frame
to 2.8e-4 rad rms with the refined prior (0.034 plain prior; 0.24–0.26
for the un-primed readings); the flat 20 nm hold-out to 7.5e-7 rad
(linear 8.5e-4).  Results, actuator
space (legacy in brackets): model correction alone, linear reading —
single-on-base floor 744 → 67 pm (SNR 9 → 80), grid-on-base SNR 1.46 →
14 (the "undetected" scenario detected by EVERY reading), dense random
42 → 9.6 nm.  I+ on the 30 nm state: floor 13 pm / SNR 584 (96×96), 24
pm / 466 (48×48); grid-on-base SNR 34 / 115; and on the break-scale
ladder **I+ holds to 60 nm rms working state (gain 0.85–0.91 at 96×96,
1.00–1.07 at 48×48) where the four-frame stepped reading falls to 0.56
/ 0.52 and the un-primed one-frame readings cliff between 30 and 40 nm
(wrong-branch pixels re-propagated into b)**.  The stepped reading keeps
dense random (3.3 nm vs I+ 5.8).  Spec: grid-on-base SNR ≥ 5 from one
frame MET; hold-out raw gain within 3% MET at 48×48 (0.997), 0.90 at
96×96 = the NGRID-193 dev grid sampling a 1 mm actuator at ~2 px (not
the reading, not the regularization).  Record: `zwfs_dm96/README.md`
S7 bullet + banner, `zwfs_s7iter_report.txt`, `zwfs_s7iter.png`.  Open:
NGRID 385; S6 color re-run on the corrected model; S5 noise pricing of
I+; deck fold (deck_zwfs is on the legacy model).

## S8 -- the runner, NGRID 385 / model 2048, colour re-run, noise of I+ (2026-09-10, RUN)

Dave: "work down the open list, your sequence; report each item as it
arrives; a PARAMETERIZED RUNNER users can modify and rerun without AI --
keep updating it as we go" (standing rule for all build tasks).
Delivered `zwfs_params.m` + `zwfs_run.m` (+ `zwfs_run_figs`,
`zwfs_run_batch`, `zwfs_batch.sh`) in `zwfs_dm96/`; stages bench /
battery / colour / noise / figs; readings L F I I+ S; README "Run it
yourself".  Equivalence gate: defaults reproduce the S7 record (64
row/ladder lines, 8 differ in the last digit).  Every item below ran
THROUGH it (`runs/<tag>/`).

- **NGRID 385:** `ng385` (1024, spot 2.0; dimple 3.96 px FAILS the 6-px
  line), `ng385s3` (spot 3.0, 5.94 px), `m2048` (MODEL 2048 via a
  trimmed size table, `macos_param_2048.txt`, dropped into the run dir
  where `find_macos_file` looks first: mGridSrf 200->4, mpts->512,
  mElt->64, mGridMat UP to 512 for the 384 DM grid; 7.92 px AND 5.03
  px/actuator = fully compliant; 32.5 min, <4 GB).  The hold-out raw gain
  moves 0.900->0.935 (96x96) with NGRID and is then identical at 1024 and
  2048 and at spot 3.0; every actuator-space row at 385 agrees between
  1024 and 2048 to 3 digits although the reference-wave profile
  (|Eb|/|E0| vs radius, new bench diagnostic) differs 2.5% between the
  4-px and 8-px dimples.  Sampling trade (Dave): mask px per lam/D =
  fill*MODEL/NGRID, detector px/actuator ~ NGRID -- opposite ways; only
  MODEL buys both.  Spot 3.0 = deeper dimple-passband dip only; spot 2.0
  stays.  Multi-site ladder (47 grid sites, `battery.ladder_sites`
  'grid'): I+ holds to 40 nm rms at 96x96 and 50 at 48x48; the record's
  "60 nm" was one site.
- **S6 colour on the corrected model:** not a lever -- combination min
  transfer 0.991 vs best single 0.962; rows neutral (I+/S) or worse (L,
  the 480 nm channel goes negative on the 30 nm base); 780 nm is the best
  single colour = a RANGE lever, not null-filling.
- **S5 noise of I+:** N(1 pm) per state L 5.4e13 / F 3.7e13 / I 6.4e13 /
  I+ 8.8e13 (prior noise costs 5%) / S 1.0e14 -- ~25x cheaper than the
  defocused-model pricing; floors converge to the battery's systematics.

Open: kernel measured AT the hold-out site (the remaining 6.5%); deck
fold (item 5, on Dave's steer).

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
