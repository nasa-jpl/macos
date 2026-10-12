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
  every candidate's error budget, and λ_m a knob to trade against noise.  Measured
  (TO, `pupil_blur_lam_m`): the bias is 306 / 31.6 / 3.2 pm at λ_m 1e-3 / 1e-4 / 1e-5
  on the 30 nm surface, photon-independent; the noise floor plateaus at 0.90 pm
  (1e14 photons) / 2.85 (1e13) — the deck's class.  The row's number is the bias at
  the λ_m the bench's Stage E reg sweep re-gates on its measured matrix.
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

6. **Full actuator resolution, no low-order forms (Dave 2026-10-11):** "Some modes
   — commanding EFC probes, for instance — will be at full actuator resolution, so
   we should not spend time investigating lower-order implementations."  The
   `format` stage's speed-vs-resolution sweep is REPLACED: the second sensor class
   is a FULL-LATTICE, well-free (continuous or lock-in) reading at ≥ 1 element per
   actuator — the lock-in camera (heliCam class), the ToF demodulating sensors and
   the Lincoln-Lab-class 128 × 128 coherent arrays are the candidates "if credible
   in the 5–10 year time frame."  Model it as a row ASSUMING a sensor at 1:1 or
   better with the actuator pitch: the question is cadence and photons per
   element per second, not resolution.  A bounded survey of those three classes
   runs in parallel (`NOTE_lockin_arrays_survey.md`, a web survey to be verified);
   the coarse-array assessment below stands as the record of why coarse forms are
   not pursued.

### On 4 — is a photodiode array remotely implementable? (CC's assessment, 2026-10-09; superseded in scope by ruling 6)

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

## 8. Stage records (as they land; the numbers are the records')

- **bench — DONE 2026-10-11 (TO, resources 362155c; gate `tCophiBench` 5/5, own class, model 1024):**
  the LO arm is `twyman_green('ref_arm','lo')` (the default `'flat'` byte-identical on four
  configs): a fiber tip + collimator with its powered face AT the flat's seat, both rigs.
  Registration PASS: chief 0.000 px on both; LO coverage 1.000; the LO overfills the test
  pupil by 34 px (radius ratio 1.215), so sub-overfill offsets are fill, not registration.
  Flat-DM null: lens 0.000 nm; mirror rig 0.782 nm rms after piston (tilt 0.647, defocus
  0.170 — the LO collimator's wavefront against the test arm, frozen into the null as the IFO
  does; its MOTION is the LO-drift row).  Contrast g 0.982 at lo_ratio 1.47 (= 2√r/(1+r),
  consistent).  Records `cophi_dm96/runs/bench_{lens,oap}/*_cophi_report.txt`.
  **Two findings:** (1) the geometric tail places the detector field BY RAY INDEX, not by
  camera position — the LO's own collimator makes the per-index null 681 nm where the
  by-position null is 0.45 nm; the LO is placed on the test pixels by position
  (`lo_on_test_`); the vector-mode field phase measures −kW on both decks (the known
  conjugate).  (2) The tg96 sheet's mirror-rig default (`POL_IN` collimated) puts the record's
  Ex = Ey launch CROSSED with the input polarizer (1e-6 of the power); the record runs pass
  `POL_IN 'source'` and never met it; COPHI applies the record's settings and solves both
  launches.  **CC review note:** the drop-in keeps the record's polarizing train (input
  polarizer, analyzer at −43.85°; the LO passes the analyzer at 1.000, the test arm at 0.265
  of its open power) — COPHI's own layout has none (same polarization both arms, a plain
  splitter), so the `noise` stage runs the photon budget for BOTH trains and the deck says
  which it prices.

## 9. Two considerations added by Dave (2026-10-11) — COPHI now, every sensor in the composite

### 9.1 Electronics architecture and complexity

A column in the composite deck's comparison, and for COPHI a slide of its own.  Per
sensor, the same rows: light source and its stabilization; modulation hardware;
detector and readout; timing / reference; real-time processing and data rate; loop
latency (sensor delay in the hold servo); parts count and single-point failures;
power / thermal; TRL.  The first cut, from the layouts as drawn:

| | TG96 PSI (record) | ZWFS (laser) | PDI / P-SRI | **COPHI (camera)** | **COPHI (lock-in)** |
|---|---|---|---|---|---|
| source | HeNe, isolator | dedicated laser | laser | single-frequency laser, 1×2 coupler, PM fiber | same |
| modulation | piezo + HV driver (4 steps) | none (dimple fixed); mask stepper for V | pinhole stepper / shutter; photonic phase shifter (P-SRI) | **2 AOMs + 2 RF synthesizers + Rb standard + ×4 multipliers + mixer → 30 Hz trigger** | same AOM pair, beat at kHz |
| detector | sCMOS camera | camera | camera | 1k² global-shutter camera + grabber | lock-in camera (I/Q out) |
| timing | frame trigger ↔ piezo settle | frame clock | frame clock | **camera trigger locked to the beat** | reference input = the beat |
| processing | atan2 per pixel per 4 frames | reconstruction (modal / matrix) | atan2 / fit | atan2 per pixel per beat (1k² at 7.5 Hz: trivial) | on-chip; I/Q at kHz |
| data rate | 4 × 1k² × 16 b × 7.5 Hz ≈ 0.5 Gb/s | one frame per reading | frames | same as PSI | 2 × 10 b × 0.27 Mpx × 1 kHz ≈ 5 Gb/s if streamed; or decimated on-chip |
| loop latency | 4 frames + settle | 1 frame | 1–5 frames | 4 frames (one beat, 133 ms) | one I/Q period (ms) |
| new error terms from electronics | piezo nonlinearity, settle | none | stepper repeatability | **synthesizer phase noise → beat jitter (`P.cophi.jitter`); AOM drive ripple → contrast; trigger jitter** | reference phase noise; ADC bits |
| parts that are new to the bench | HV driver | mask | pinhole + stage | AOMs, RF chain, Rb reference, fiber LO, collimator | + the lock-in camera |

What the model must carry: the synthesizer / reference phase noise mapped to
`P.cophi.jitter` (rad of beat phase per frame — a 1e-11 Allan-class Rb reference at 80
MHz contributes nothing; the AOM drivers' DDS phase noise is the term), AOM drive
amplitude ripple as a contrast modulation, the camera trigger's jitter, and the sensor
DELAY in `dmg_loop` (today's loop steps per reading; a reading that takes one beat is one
cycle of latency at 7.5 Hz — the servo bandwidth row).  Score in the deck: parts count,
the number of active control loops the sensor needs to run (COPHI: laser frequency,
two RF chains, fiber polarization/phase), and the failure the loop cannot survive.

### 9.2 Multi-color (a tunable laser) for capture-range extension

Single-wavelength phase has an unambiguous SURFACE range of ±λ/4 (double pass: OPD
λ/2) — ±158 nm at 632.8 nm, ±330 nm at 1319 nm.  The HOLD job (3 pm) never needs
more; CAPTURE of the DM's initial 100–200 nm rms figure (peaks 3–4× rms) wraps at
633 nm, and a dropout mid-run loses the fringe order.  Two cures, both to be priced:

- **Temporal unwrapping** (COPHI gets it free: a reading every beat, the surface
  cannot move λ/4 between readings) — covers a running loop, not a cold start.
- **Two wavelengths → a synthetic wavelength** Λ = λ₁λ₂/|λ₁ − λ₂|: the coarse
  phase at Λ picks the fringe order of the fine phase at λ.  Range ±Λ/4; the coarse
  measurement's noise is amplified by Λ/λ, and it must resolve the fine fringe
  (σ_coarse × Λ/λ < λ/4) — so Δλ is a trade: 633 ± 0.5 nm gives Λ = 0.8 mm (range
  ±200 µm, amplification 1265×, needs σ_fine < 0.125 nm/1265 ≈ 0.1 pm — too far); a
  two-step ladder (Δλ 20 nm → Λ 20 µm, then 0.5 nm) is the standard answer
  (SIM/MSTAR-class absolute metrology did exactly this with heterodyne).  Hardware:
  a tunable source at 633 (ECDL, a few nm, or a frequency comb line) is expensive;
  at **1310–1550 nm the telecom tunable lasers are cheap, agile and PM-fiber native**
  — which strengthens ruling 2's 1.3 µm row: it buys an InGaAs camera but also a
  practical tunable source and twice the single-wavelength range.
- For COPHI the two colors can run SIMULTANEOUSLY on one camera with two AOM pairs at
  two beats (frequency-multiplexed heterodyne); the lock-in camera reads one reference
  at a time (sequential colors).  For the ZWFS the sensor's own range is the dimple's
  linear range (~λ/10), so a wavelength ladder extends it only as a sequence of masks
  — a different trade; for the PSI/PDI the two-wavelength form applies as for COPHI.

Model rows: `capture` stage gains the two-wavelength ladder (Δλ sweep: range,
amplification, the fine-noise requirement at each rung, the rungs needed to capture
200 nm rms from cold), priced in photons and seconds; the 1.3 µm stage carries the
tunable-source argument.  Composite deck: a "capture from cold" row per sensor —
temporal, two-color, or neither.
- **battery — DONE 2026-10-11 (TO, resources 5d2383c; `tCophiBench` 5/5 rerun, `tCophi` 8/8):**
  the H reading runs through `tg96_run`'s own place + deck stages (opt-in
  `P.reading_ctx = 'cophi_ctx'`; absent hooks leave the record's path unchanged).  The I
  column is the SAME surface as the record (`redo96_<optics>ph`: tail md5, lit 6948, hold-out
  (62,62), 112-site grid, seed 7, the 30 nm surface, the 64-state matrix on it); the deck's
  pre-cap headline (7500 lit, 120 sites) is printed beside it, labelled.  Mirror-rig placement
  reproduces the record's D1 gate to the digit.  **Headline (gain / err pm / floor pm), H vs I
  on the same surface:** lens single 0.9947 / 1.7 / 2 vs 1.0001 / 1.4 / 1; grid 0.9964 / 2.4 / 2
  vs 1.0046 / 2.2 / 2; rand10 0.9967 / 180 vs 1.0053 / 179; mirror single 0.9941 / 1.4 / 1 vs
  0.9832 / 2.7 / 2; grid 0.9949 / 1.6 / 1 vs 0.9838 / 2.9 / 2; rand10 0.9950 / 118 vs
  0.9844 / 223.  Capture to 10 %: 480 nm+ on both rigs; the re-measured matrix's gain 0.996 /
  0.995 at 60–160 nm.  **OPEN: H sits 0.5 % below unity on BOTH rigs, every row, ladder and
  re-measured matrices** (the brief expected 1.000 ± the matrix bias); candidates the LO
  placement interpolation, the mask (LO disc vs interference support), the frozen-null form —
  TO isolates it before `noise` (CC's discriminators in the reply of 2026-10-11).  Also flagged:
  the deck stage's "flat-DM null" rms (6.8 nm on H) is taken over the whole mask including the
  no-test-light annulus; the record's 59 nm line may carry the same dilution — a composite-deck
  check.  Records `cophi_dm96/runs/battery_{lens,oap}/`.
  **RESOLVED 2026-10-11 (TO, resources fd4124d, `runs/gainprobe_lens/`): the 0.5 % is the
  ESTIMATOR's Tikhonov shrinkage at the record's λ_m = 1e-3, not the reading.**  Five
  one-run discriminators on the lens rig: (1) a synthetic 1.000 nm poke read at the pixel
  level against the engine's own field change: slope 1.000000 — the reading is EXACT; (2) LO
  := the test arm's own flat field, no placement: unchanged; (3) scalar trace (contrast 1.0000):
  identical — not the vector sum or the polarizers; (4) the test-beam mask: unchanged;
  (5) λ_m = 1e-5: 1.0007 / 1.0022 / 1.0024.  So an exact reading reads 0.995 through the
  record's estimator — the brief's "± the matrix bias", now a number.  **Corollary for every
  gauge deck:** the record's I column runs the same estimator and λ_m yet reads 1.000–1.005
  (lens) and 0.983 (mirror): its READING carries +0.5–1 % and −1.1 % against an exact one — the
  snapshot's polarization gain, which the matrix does not absorb and which the lens-rig
  "1.000" had been hiding by cancellation.  **Dave's rule (2026-10-11): regularization is
  driven to convergence, not dropped** — the converged estimator (the blur demo's `iterate_`
  form, as a `dmg_act_fit` mode) becomes the accuracy-row estimator for EVERY gauge in the
  composite deck, and the record's polarization gain is then reported in the open.  Order:
  `noise` first (a 0.5 % gain is nothing to a noise row), then the converged estimator as a
  cross-gauge item re-scoring the battery rows for H and I.  Diagnostic knobs kept, opt-in:
  `cophi.diag_lo` / `.scalar` / `.diag_mask`, `place.lit_force`.
- **noise — DONE 2026-10-11 (TO, resources c24d731; records `runs/noise_{lens,oap}[_plain]/`):**
  the deck's 30 nm matrix estimator through `cophi_ctx`'s noise hook; each camera mapped onto
  the model grid by BINNING (well × k², read × k; frame mean at 0.5 × well), 12 draws; the
  per-pixel σ_φ measured / `tCophi`'s formula = 1.000 ± 0.002 on every camera, rig and train.
  **Metric:** pm per reading = the PISTON-FREE robust floor (1.48 MAD over the unpoked lit set,
  median of draws) — the raw single-actuator σ is mostly a common lit-set piston (165 pm raw
  vs 76 piston-free; piston is unobservable and inert on a DM) and the per-actuator distribution
  is heavy-tailed (p99 ~4–6× the median; the plain std moved ±30 % between benches, the robust
  floor 1–3 %).  **The record's S5 "photons for 1 pm" uses the RAW hold-out σ (piston
  included) — the composite deck must put every gauge on the piston-free metric (CC).**
  **Headline, pm per reading, lens / mirror (seconds per pm):** Feng's baseline 1k / 60 ke /
  30 fps 57.7 / 53.7 (443 / 385 s — ~7 min per pm, inside the first look's 6–9); Photonfocus
  200 ke / 150 fps 39.5 / 36.5 (42 / 36 s); pco.edge GS 40.2 / 38.3 (129 / 117 s); C-RED 2
  19.6 / 18.7 at 632.8 nm, a placeholder for the 1.3 µm stage.  **Trains:** pm equal (the
  camera is well-bound); the record's polarizing train costs light only (test throughput 0.271
  vs 0.408, ×1.5 laser) — at ~0.2 µW on the camera a mW laser covers either, so Feng's plain
  splitter is the one to price.  **Electronics, one reading, linear in each term:** jitter
  1e-3 rad/frame ≈ 19 pm (1e-2 ≈ 195); AOM ripple 1e-3 ≈ 25 pm; fiber walk 0.01 rad/frame ≈
  66 pm (0.05 ≈ 330).  Budgets to stay under the 55 pm photon floor — **CORRECTED by the loop stage (TO 2026-10-12): the noise
  stage's 12-draw fiber row sat 2.3× LOW (a reading's whole walk is three random numbers; the
  seed set was low); exact per reading per unit: walk 1.62e4 pm per rad/frame, jitter 1.96e4
  pm/rad, ripple 3.1e4 pm/unit → jitter < 2.9e-3, ripple < 1.8e-3, differential fiber walk
  < 3.5e-3 rad per frame (0.35 nm of OPD per frame, ~10.6 nm/s, ~1.1 mK/s per metre of
  uncommon fiber at k_T 1e-5/K, to verify) — the fiber is THE term**: the engineering requirement is the differential phase stability of the
  two PM fibers between the coupler and the collimators (short, co-routed, thermally tied; or
  an active fiber-phase servo).
- **converged estimator — DONE 2026-10-11 (TO, resources "dm_gauge_lib: the converged estimator" + cba9b6e the fiber budget in engineering units):** `dm_gauge_lib/dmg_converge.m` (the shared `iterate_` form), `tg96_run battery.est_mode 'onestep'` (default = the record) | `'converged'` | `'both'`, `dmg_act_fit opts.converge`; gate `tConvergedEst` 4/4 (red leg: one-step < 0.999 on an exact reading; converged 1.0000 ± 1e-4 on a matrix and a stencil forward model; the noise trade MEASURED: ×2.1 white-noise floor at a zero-update stop, 2 steps when stopped at the one-step floor).  Re-score noise-free, both rigs, 5–7 steps to ≤ 0.001 pm; the one-step rows reproduce `battery_*` and `redo96_*ph` to the digit.  **Converged gains (single / grid / rand10): H lens 1.0008 / 1.0022 / 1.0024, H mirror 1.0000 / 1.0005 / 1.0005; I lens 1.0064 / 1.0106 / 1.0112, I mirror 0.9887 / 0.9892 / 0.9896 → the record interferometer's polarization gain in the open: I/H = 1.0056–1.0088 (lens), 0.9887–0.9891 (mirror).**  Open, small: H lens keeps +0.1–0.2 % when converged (a lens-rig matrix-model term; first candidate the box window truncating the influence wings; the mirror rig is exact).  **Pinned (C5), for the ZWFS re-score: `dmg_act_fit`'s unknowns are ALL actuators with lit-masked data, so with an unlit border the converged gain PLATEAUS (0.962 at a w = 1.0 stencil); CC's ruling 2026-10-11: the estimator's unknown set = the CONTROL set (lit by the record's `lit_margin` rule, Dave's descent-stall ruling) for every gauge — an unobservable actuator does not belong among the unknowns.**  Fiber budget in engineering units: 0.81 nm of differential OPD per frame, 24 nm/s, ≈ 2.4 mK/s per metre of uncommon fiber at k_T = 1e-5/K (coefficient to verify).  Records `cophi_dm96/runs/conv_H_{lens,oap}/`, `tg_psi_dm96_oap/runs/conv_I_{lens,oap}ph/`.
- **loop — DONE 2026-10-12 (TO, resources 25853f5; + 5a2e3f9 the unknown-set ruling, `tConvergedEst` C5 1.000000; + f7d37d7 `dmg_loop opt.delay`, `tDmgLoop` 17/17):** plain splitter, the record's knobs (g 0.5, K 60, seed 77), both estimators, M readings averaged per cycle to reach 3 pm.  **Hold steady state, pm, lens / mirror:** baseline camera (M 120, a cycle every 16 s): noise 2.75 / 2.68, walk 3.60 / 3.54, thermal 10.3 (lag rate/gG — every gauge), camera drift 2.85 / 2.77 (immune), one-cycle delay 3.65 / 3.55 (the delay costs ~0.9 pm), **LO drift 0.1 nrad tilt + 1 pm defocus per cycle: 29.6 / 29.4**; **Photonfocus (M 60, a cycle every 1.6 s): noise 2.60 / 2.53 — 3 pm held 10× faster**; fiber at the corrected budget 4.21 / 4.12 piston-free (= the √2 linearity prediction), at 3× the budget 9.9.  One-step = converged to 0.1 pm on every row.  **Two deck findings:** (1) the LO is a term the loop cannot see past — it reads LO drift as DM drift and the DM follows it; the LO pointing must hold ~10 prad per cycle (~10 pm of fiber-tip lateral motion at F1) or be referenced (the brief's 1 µrad would be ~24 nm of DM tilt per cycle); (2) the fiber walk is 30–40 % piston — `dmg_loop` now reports the piston-free `ss_pf` / `bias_pf` / `sig_n_pf` beside the record metric, identical on every non-fiber row.  Records `runs/loop_{lens,oap}_plain/`, `runs/loopfib_*/`.  Remaining: capture, format, figs.
- **capture / format / figs — DONE 2026-10-12 (TO, resources 30b124e):**
  **Capture:** truth = the engine's per-pixel OPD unwrapped (H = −ΔOPD/2 to 3e-7 nm noise-free).
  Two colours 632.8 + 652.8 nm (Λ 20.65 µm, range ±5.2 µm) frequency-multiplexed on the SAME
  8 frames (7.5 / 3.75 Hz beats, each read by `dmg_cophi_recon` at its own buckets; crosstalk
  4e-16 rad; λ₂ a re-trace of both decks): **zero wrong fringe orders on the whole capture ladder
  to 480 nm rms, noise-free and on both cameras** (0.15 / 0.10 nm per pixel); single colour wraps
  from 40 nm rms (74 % wrong at 480).  **So 200 nm rms from cold is captured in ONE multiplexed
  reading (0.27 s): the §9.2 ladder needs one rung for this job**; any Δλ from 1 to ~100 nm
  resolves (fine requirement 1.7 nm at Δλ 20 vs ~0.1 nm per pixel).  1.3 µm row analytic (range
  ×2.1).  Caveat: the bench lenses are dispersion-free (fixed n 1.5), so the re-trace cannot show
  the tail's chromatic refocus — ~0.5 mm on the focuser for fused silica at Δλ 20 nm (to verify);
  its static part freezes into each colour's null; the pupil-image sharpness effect is open.
  **Temporal unwrapping** holds at 0.9 × λ/4 per reading, fails at 1.1 ×: 1.19 µm/s (baseline)
  / 5.9 µm/s (Photonfocus) at the peak pixel, ~0.24 / 1.2 µm/s rms for a 30 nm-class shape.
  **Format (the lock-in class, heliCam C4 sheet: 512 sites = 5.45 per pitch, 500 ke well, 1 kHz
  I/Q):** WELL-CAPPED at 1 kHz from 1e9 e/s — 34 pm per I/Q sample, ~1.07 pm in 1 s, **~1.1 s
  per pm = 40× the Photonfocus, 400× the baseline**; the well-free rows beside (0.29 s per pm at
  1e9 e/s, 29 ms at 1e10) — "well-free" is a device property to VERIFY (a faster I/Q rate lifts
  the cap); the sheet has no phasemeter field, added in quadrature if it appears.  **Figs:**
  mirror rig, plain splitter, noise-free, deck format: the 30 nm surface reads 31.95 nm rms,
  read − engine 0.31 pm rms; the flat reads 0.  **Fix found en route:** the plain MIRROR train
  inherited the record's `oap.SRC_TRIM` (−0.37 mm), solved WITH the polarizer plate in the
  diverging leg — without the plate OAP1 was fed off focus (an 84 nm defocus null); now
  `SRC_TRIM 0` for the plain mirror rig: null 0.746 nm; capture identical; noise floor 54.1 →
  53.2 pm; walk budget 3.4e-3 → 3.6e-3.  The plain-mirror LOOP rows ran mis-fed; differential,
  agreeing with the lens rig to 1–3 % — rerun queued (3.4 h, background).  Records
  `runs/capture_*_plain/`, `format_lens_plain/`, `figs_oap_plain/`.  **All brief stages are in.**
  Open: `tCophiBench`'s poke / LO-tilt / no-interference / well = 0 legs (§3).
- **gates closed / wording corrected — 2026-10-12 (TO, resources f6cd069):** `tCophiBench` 9/9 (own class, ~15 min): B5 a 1 nm poke through the matrix = 0.9949 one-step / 1.0010 converged; B6 a 1 nrad LO tilt reads as a DM SURFACE tilt of 1.008 × θ/2 (0.23 pm left after the plane) — the matrix takes LO tilt as DM tilt one-for-one (0.467 in COMMAND units, the influence overlap); R2 / R3: no interference (an incoherent LO) and a dark camera (well 0) are refused BY NAME, `cophi:noFringe` — a contrast floor in the H reading itself (median amp/dc over the pupil < `P.cophi.contrast_floor`, default 0.05), with a must-pass control; `dmg_cophi_*` untouched.  **Format wording, corrected (TO):** the surveyed lock-in, even WELL-CAPPED at its 1 kHz I/Q rate, is 1.1 s per pm — **already ~40× faster than the Photonfocus (42 s) and ~400× Feng's baseline (7 min)**; a verified well-free device (continuous, or I/Q ≥ 10 kHz) would add ~4–40×.  The vendor question for Heliotis (through Dave): per-pixel phase noise vs flux at t_c = 100 ms; does the FWC bound each internal 1 MSa/s sample or the 5.6 kSa/s output.  The plain-mirror loop rerun on the fixed feed (SRC_TRIM 0, fiber rows at the corrected budget) is running detached (~3.4 h); the mis-fed rows moved to `runs/oap_plain_misfed/`.
- **mirror-rig loop on the FIXED feed — DONE 2026-10-12 (TO, resources "cophi_dm96: the plain-mirror loop rerun on the FIXED feed"; `runs/loop_oap_plain/`, 206 min):** every row within ~0.02 pm of the mis-fed run (baseline noise 2.67, walk 3.53, thermal 10.34, cam 2.77, delay 3.54, LO 29.35; Photonfocus 2.52); the fiber rows at the corrected budget, piston-free: baseline 3.95 pm (9.16 at 3×), Photonfocus 4.84 (12.65 at 3×); the 8e-3 / 2.4e-2 rows exist on the lens rig only, marked.  The mis-fed runs kept in `runs/oap_plain_misfed/`.  deck_cophi's hold table updated.
