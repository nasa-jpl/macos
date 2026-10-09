# BRIEF — the COPHI campaign: a fifth DM gauge on the same bench, then the composite deck

CC for Dave, 2026-10-09, on receipt of Feng Zhao's charts
(`~/dev/MACOS_sandbox/DM_gauge/COPHI_for_DM_gauge.pptx`, "HWO Coronagraph
Technology Proposed Implementation Plan (FY26) — DM Drift Sensor", 2026-05-19).
Plan, not yet work: Dave rules on §6, then it runs.  Lanes proposed in §5.

## 1. What Feng's deck adds to the first look

(`NOTE_cophi_first_look.md` was written from two captured slides; the deck has
twelve.)  New in it:

- **The purpose is HOLD.**  "DMs operate open-loop and drift (material physics;
  'floating electrodes' from manufacturing; temperature coefficient ~2 %); it is
  the tallest tent pole in contrast stability."  The sensor runs each DM in a
  LOCAL control loop, to "drastically improve DM stability, relax DM
  manufacturing tolerance, improve DM yield."  That is our job 3 (hold 3 pm in a
  servo) and the capture job only as far as the loop's own range; it is NOT the
  picometer absolute-figure job.  The comparison must be scored on hold first.
- **Two concepts side by side:** the Zernike WFS with a dedicated laser (Ruane et
  al. 2020, demonstrated at DST, in Fang's VSG plan) and COPHI (SIM Diffraction
  Testbed and TOM3) — Feng's plan carries COPHI.  So the composite deck's
  comparison is exactly the trade the project is making.
- **The algorithm as he will implement it:** I = I0 (1 + γ cos(φ + 2π f_h t));
  four frames at α = 0, π/2, π, 3π/2 set by the AOM drivers, camera triggered from
  the same clock; φ = atan2(I4 − I2, I1 − I3); surface H = ½ (φ/2π) λ (double
  pass).  Backup: the generic N-frame least-squares solution (Zygo's four-step is
  the N = 4 case) — so the N-bucket form is in scope and is the known cure for
  step error.
- **Implementation:** as in the first look — one single-frequency laser, 1×2
  coupler, AOM-1 at 80 MHz, AOM-2 at 80 MHz + 7.5 Hz from two synthesizers on a
  10 MHz Rb standard; the ×4 multipliers mixed to 30 Hz for the camera trigger;
  7.5 Hz fringe, 30 Hz frames; measurement beam fiber → lens → DM; LO fiber →
  splitter → camera; "DM surface figure error as a function of time."
- A frame of a 1 Hz heterodyne fringe movie (SuperZygo 2008) as the illustration.

Not in the deck: the camera (format, well, rate beyond "30 fps"), any photon or
noise budget, the LO arm's optics, error terms, a performance number.  Every
number in the comparison will be ours; the survey (`NOTE_detectors_survey.md`)
supplies the cameras.

## 2. What the campaign answers

1. **COPHI scored on the same bench, same rows, same surface as the four gauges**
   (accuracy gain / floor on the 30 nm surface with the matrix on it; raw map vs
   the engine; precision vs photons; repeatability; hold at 3 pm; capture range;
   drift within a measurement) — one more column in every table of the deck.
2. **Its own error terms:** the LO arm's wavefront and drift (tilt, defocus per
   cycle); differential fiber phase within the four frames; trigger jitter;
   quarter-cycle integration; N = 4 vs N = 5–8 buckets; the beat-to-frame ratio.
3. **The camera as the variable:** photons per second at full well for each
   surveyed camera → pm per reading at its cadence; the servo at that cadence;
   Feng's 1k²/30 Hz baseline as the "as proposed" row and the best streamable
   true-global-shutter camera as the "available" row.
4. **Speed versus resolution:** formats 1k², 512², 240², 128² at their rates, the
   read as the cell average when the format is coarser than the actuator pitch
   (the box kernel; `BRIEF_to_pupil_blur.md` step 5), and the photodiode-array
   form (N channels at 10 kHz, 49 pm per zero crossing at 633 nm) — what each
   recovers of the 96 × 96 lattice and how fast.
5. **The comparison Feng's deck poses:** COPHI vs the Zernike sensor with a
   dedicated laser, for the hold job, on the mirror rig.

## 3. The template: `templates/40_benches/cophi_dm96/`

Same shape as the three benches (one parameter sheet, one runner, stages, run
tags, README, a gate), on `dm_gauge_lib`.  Reuse, do not re-derive:

- **The bench** is `tg_psi_dm96_oap`'s (lens rig and mirror rig, the collimated
  source, the splitter node at 22.5°, the imaging tail, the camera) with the
  reference flat REPLACED by an **LO arm**: a second fiber source and collimator
  entering the recombiner.  Precedent for a separately traced reference arm:
  `macos.design.psri_bench` (the P/SRI's two-deck reference).  The LO's
  wavefront is the engine's trace of that arm; its drift is a per-cycle tilt /
  defocus knob.
- **The reading class 'H'** (heterodyne) in `dm_gauge_lib/dmg_cophi_gauge.m`,
  after `dmg_ifo_gauge`: frames I_k = ∫ |E_t + e^{iα(t)} E_LO|² dt over each
  quarter cycle (the integrating bucket: contrast 0.90), α from the beat with a
  trigger-jitter term and a differential-fiber-phase walk within the scan;
  four-bucket and N-bucket least squares (Feng's backup slides); then the
  record's matrix estimator (`calib_mode 'matrix'`).  Same polarization in both
  arms (no analyzer, no QWP) — the model should show the polarization gain is
  absent by construction, not assume it.
- **Stages:** `bench` (layout, clearance, the LO arm placed), `battery` (the
  readings on the flat and the 30 nm surface; the matrix on both), `noise`
  (photons per pm; the camera table), `loop` (hold at 3 pm; walk, thermal ramp,
  camera walk; the LO-drift row), `capture` (range; temporal unwrapping as the
  free extra), `format` (the speed-vs-resolution sweep: format × rate × cell
  average; the photodiode-array row), `figs` (the station figure in the deck's
  format: command, frames, raw map vs engine, reconstruction, residual).
- **Parameter sheet:** `cophi_params.m` — beat frequency, frames per beat, N
  buckets, camera (format, pitch, well, rate, read noise, from the survey, by
  name), LO power ratio, fiber-phase walk, trigger jitter, LO drift per cycle,
  the photodiode-array option (N sites, phasemeter noise per sample, sample
  rate).  Every number cited in the deck comes from a stage record.
- **Gate:** `tCophi` (SUITE_FAST if under a minute, else its own class): the
  four-bucket on a known phase (exact steps, no noise) recovers it to 1e-12;
  the N-bucket matches; the quarter-cycle integration gives contrast 0.900
  ± 0.001; a 90° step error of 1° gives the textbook second-harmonic ripple
  (must-fail leg: the four-bucket with the step error asserted > the N = 5
  value); photon scaling 1/√N over two decades.

## 4. The decks

**`deck_cophi.md` / `.pptx`** (demo_session, `make_brief_slides.py`), parallel
to the gauge slides' structure — Reads / Frames / Reconstruction / Gets right /
Gets wrong / Scored — so it drops into the composite deck unchanged:

1. Title — COPHI as the fifth gauge, scored on the same bench.
2. Feng's proposal in one slide (his layout, his algorithm, his purpose: hold).
3. What it is in our terms (the piezo four-step, continuous, frequency-set
   steps, a fiber LO; what transfers from the deck, what is new).
4. The bench with the LO arm (the engine's render, both rigs).
5. The reading, station by station (the deck's station-figure format).
6. Performance side by side — the five gauges on the 30 nm surface.
7. Photons for 1 pm and the camera table (survey cameras → pm per reading;
   Feng's baseline vs available).
8. Hold: the servo at the camera's cadence; LO drift; fiber phase; jitter.
9. Capture range, and the temporal unwrapping COPHI gets free.
10. Speed versus resolution: the format sweep and the photodiode-array form.
11. COPHI vs the Zernike sensor with a dedicated laser, for the hold job.
12. Recommendations; what to develop.  Backup: Feng's N-frame algebra; the
    survey table with its caveats; TOM3's numbers; the blur demo's box-kernel
    curve; run tags.

**The composite deck** (`deck_gauges_v2.md`, after `deck_cophi`): the 10-08
short deck with the COPHI column in every comparison, the blur backup (TO's
corrected demo), per-candidate error budgets (follow-up 1, from the "From this
model to as-built" backup), Kent's ZWFS reconstruction if it arrives (3), and
the Feng-posed comparison as a main slide.  Same rules: figures from the tools,
numbers from records, American English, Dave's sign-off before anything leaves.

### 4a. Carried into the composite deck from the blur work (TO, 2026-10-09, `REPORT_pupil_blur.md`)

- The pupil-imaging slide's sentence is corrected: the 0.13 / 0.29 % share is the
  leg's whole imaging error; blur is 3–8 % of it (legs at 0.006 / 0.016 pitch); the
  rest is registration (a third of the lens share) and low-order error (mirror).  The
  "blur 0.05–0.14 mm" figure was pupilq's tilt-band ray walk, not the coherent PSF.
- The record's matrix estimator (λ_m 1e-3) rolls off the actuator Nyquist by 3.4 %
  (the (96,96) mode, Stage D) and biases a dense random surface by ~1 % rms: a row in
  every candidate's error budget, and λ_m a knob to trade against noise.
- The box kernel settles the photodiode-array question harder than §7's assessment:
  between 1 and 1.5 pitch per element the actuator-space fit COLLAPSES (fewer readings
  than lit actuators), not merely loses the checkerboard.  A coarse array must be read
  in a low-order (modal) basis, never in actuator space — the `format` stage's
  photodiode row is a modal reconstruction of the drift, scored on the drift's
  spectrum, or it is nothing.

## 5. Lanes and order

- **TO:** `BRIEF_to_pupil_blur.md` first (in flight; its step 5 box kernel feeds
  §2.4), then the template §3 stages bench → battery → noise → loop → capture →
  format, each stage a commit with its record.
- **CC:** this plan; the reading-class design review before TO codes it (the
  integrating bucket, the step-error algebra, the LO model); `deck_cophi` as
  the records land; the composite deck; the detector numbers verified against
  the vendors' sheets before they enter a slide.
- **CCMac:** held in reserve (budget); if used, the format-sweep stage on the
  Mac once the template exists.
- **Dave:** the rulings below; the note to Feng after `deck_cophi` is signed.

Rough size: the template's bench/battery/noise/loop stages are the tg96 ones
with a new arm and a new reading — days, not weeks, on the Linux box at model
1024 (one MATLAB at a time, the batch wrapper); the format sweep and the
photodiode form are the new modeling and the part that "deserves careful
modeling on its own."  Record every stage's wall time in its log; the
comparison must not be rushed to the deck.

## 6. Rulings needed from Dave

1. **The LO arm's form:** a fiber tip + collimator traced by the engine (Feng's
   drawing), with its drift as tilt/defocus per cycle — or keep the record's
   reference flat and treat COPHI as a reading on the existing bench (faster;
   loses the LO-drift row).  Recommendation: the traced arm; it is the one new
   error term the deck needs.
2. **Wavelength:** our bench is 632.8 nm (HeNe); TOM3's COPHI was 532 nm
   (doubled Nd:YAG); Feng's deck does not say.  Recommendation: 632.8 nm for
   the comparison (same bench, same photon rows); the λ dependence is one
   line — at a fixed photon count the surface error scales as λ (phase noise
   is λ-free, surface = λ·phase/4π), so 532 nm reads 0.84× the picometers of
   632.8 nm for the same photons — stated, not run.
3. **The camera rows:** "as proposed" = 1k² / 60 ke / 30 Hz; "available" = the
   survey's best streamable true-global-shutter 1k-class part (Photonfocus
   1024² 200 ke 150 fps, or pco.edge 5.5 in GS mode) — which two, and whether
   the sCMOS rolling-shutter cameras get a line explaining why they are out.
4. **Scope of the photodiode-array form** now (one row in §2.4) or deferred to
   its own campaign.
5. **Feature freeze (PLAN_CONSOLIDATION rule 10):** this campaign is new tool
   work in resources outside the engine; confirm it is the named exception
   alongside the pupil-blur brief, and whether phase 1 of the consolidation
   starts in parallel (CC + TO split) or after.

## 7. Rulings (Dave, 2026-10-09)

1. **LO arm: TRACED, as a drop-in to the current bench** to the degree possible —
   the reference flat's seat becomes the LO's entry; the rest of the bench, the
   tail, the camera and the records stay.
2. **Wavelength: 632.8 nm is fine for the comparison.  Feng may prefer ~1.3 µm
   (SIM MET's 1319 nm Nd:YAG class).**  Carry λ as a parameter (`P.LAM`) and run
   the record at 632.8; a 1.3 µm row needs an InGaAs camera (survey: C-RED 2/3,
   640 × 512 at 600 fps, 1.4 Me⁻ well at low gain — 2.8e14 photons/s, but the
   surface error per photon scales with λ, ×2.1) and the bench's refractive
   tail refocused — a stage, not a toggle.  Capture range doubles with λ.
3. **Cameras: several options, the three named are fine (Feng's 1k²/60 ke/30 Hz
   baseline; Photonfocus 1024² 200 ke 150 fps; pco.edge 5.5 in true GS mode),
   and READ NOISE enters every comparison** — the noise model per camera is
   shot + read (σ_read per pixel per frame, four frames), with the full-well
   ceiling setting photons per frame; the photons-for-1-pm and hold rows are
   per camera.
4. **Photodiode array: assess implementability first** (CC, below) before it is a
   stage.
5. **Order: CC finishes the Dyson deck now; TO is on the blur brief; then
   re-evaluate** (the campaign's start, the freeze exception, phase 1).

### On 4 — is a photodiode array remotely implementable? (CC's assessment)

At the ACTUATOR PITCH, no: 96 × 96 = 9216 channels, each a photodiode, a
transimpedance stage and a demodulator — TOM3 had 43 detectors and multiplexed
them ten at a time.  As a COARSE array, yes, with custom electronics of a known
kind: silicon photodiode arrays exist to 16 × 16 (Hamamatsu S-series) and 2-D
APD/lidar arrays to 32 × 32 with per-pixel readout; 256–1024 channels digitized
at ~100 kS/s for a 10 kHz beat is ~1e8 samples/s, an FPGA board set, with the
phase per channel by digital I/Q demodulation rather than TOM3's 128 MHz
counters.  What it buys: no well-depth ceiling — a photodiode takes the whole
beam (1 mW is 3e15 photons/s), three orders above any camera's photons per
second — so picometers per second on every channel, continuously, at kHz
cadence.  What it costs: a 32 × 32 array over the 96 × 96 lattice is a 3-pitch
cell — the box-kernel regime where the checkerboard (the single-actuator
scale) is gone and only the low and mid orders are recovered (TO's step 5 will
put a number on it).  That matters for Feng's own motivation: "floating
electrodes" are SINGLE-actuator failures, which a coarse array cannot localize,
while the thermal drift (~2 %/K, low order) is exactly what it reads well.  So
the realistic form is a HYBRID — the camera four-bucket for the full lattice at
its cadence, a coarse photodiode array for the low-order drift at kHz — and
that is the trade the `format` stage should quantify: cell size × cadence ×
photons against the drift spectrum.  Verdict: implementable as a coarse
companion sensor, not as the full-resolution gauge; model it as a row, decide
after the number.
