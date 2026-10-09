# COPHI as a DM gauge — a first look (CC for Dave, 2026-10-09)

Sources: Goullioud, Lindensmith, Hahn, "Results from SIM's Thermo-Opto-Mechanical
(TOM3) Testbed," Proc. SPIE 6268, 626824 (2006) (`MACOS_sandbox/DM_gauge/
TOM3_COPHI_Renaud.pdf`); Feng Zhao's two slides captured in `gauges notes.pptx`
("Concepts" and "Implementation Example using 4 Frame Method"); the gauge deck of
2026-10-08 for every number of ours quoted below.  Feng's own slides are still to come;
this is a read of the two captures.

## 1. What is proposed

Feng's layout: one single-frequency laser, a 1×2 coupler, two AOMs driven from two
synthesizers on a 10 MHz Rb standard at 80 MHz and 80 MHz + 7.5 Hz.  The measurement
beam (F + 80 MHz) comes out of a fiber, is collimated by a lens onto the DM and
returns; the local oscillator (F + 80 MHz + 7.5 Hz) comes out of a second fiber and
is combined with the return at a beamsplitter in front of a camera that images the
DM.  Every pixel sees a 7.5 Hz beat; the camera runs at 30 Hz = 4 frames per beat,
triggered from the same standard (the 320 MHz / 320 MHz + 30 Hz mix gives the 30 Hz
trigger), and the four frames give the phase per pixel by the four-bucket formula.
The output is "DM surface figure error as a function of time"; the rate Dave quoted
is 1 Hz (about 7 beats per reading).

The TOM3 COPHI is the same interferometer with a different detector: 43 photodiodes
sampling the pupil, a 10 kHz beat (80 vs 80.01 MHz), and a VME phasemeter counting a
128 MHz clock between zero crossings — 532 nm × 10 kHz / 128 MHz = 40 pm per zero
crossing, averaged on board to sub-picometer "after a few seconds."  Phase is read
relative to the central detector, so the two fibers' and the laser's phase noise are
common mode.  Their pupil camera was an alignment aid run with the beat slowed to
~1 Hz and 360 frames, "a resolution similar to the Zygo."  Results they report:
compressor OPD vs temperature 4 nm/K measured against 4.98 predicted; guide
interferometer 99 pm wide angle / 5.4 pm narrow angle after SIM's chop-and-average
processing; a 240 nm P-V wavefront change for an 8.4 K soak matching the FEM to the
mount print-through.

## 2. What it is, in the terms of our deck

**It is our Twyman-Green interferometer's absolute-calibration leg — the piezo
four-step — with the phase steps supplied by a frequency offset instead of a
piezo, run continuously, and with the reference flat replaced by a second fiber.**
Everything the deck says about the interferometer transfers unchanged:

- The reconstruction is the same: atan2(I2 − I4, I1 − I3) per pixel, then the
  regularized fit to the response matrix measured through the gauge, in actuator
  units, piston removed.
- The photon cost is the same: the reference light is not signal, so 7e14 (mirror
  rig) to 1e15 (lens rig) photons per 1 pm reading on the 30 nm surface; 2.0–2.3e13
  photons per cycle to hold 3 pm.  Setting the LO level equal to the DM return gives
  unit fringe contrast, exactly our flat-arm case.  The integrating-bucket loss of a
  quarter-cycle exposure is a contrast factor sin(π/4)/(π/4) = 0.90, no phase bias.
- Capture is the same: an external reference reads any figure, 480 nm per pixel and
  beyond with unwrapping; from 200 nm of wavefront with unwrapping alone.  Plus one
  thing for free — the phase at each pixel is tracked continuously in time, so
  unwrapping across readings is trivial (this is how TOM3 followed hours of drift).
- The pupil imaging is the same: the DM is imaged on the camera through a lens and a
  splitter.  Fang's concern about blur and distortion applies to COPHI exactly as to
  our bench; it does not escape it.  (CCMac's blur demo answers it for both.)
- The drift rows are the same: a measurement spans 133 ms (four frames); DM and camera
  drift within a scan behave as the "Drift within a measurement" backup says; the
  thermal-ramp lag of the servo is rate ÷ (gain × response) whatever the gauge.

What COPHI changes, for the better:

- **No polarization gain.**  Both beams are the same linear polarization from PM
  fibers; the splitter's diattenuation changes contrast, not phase.  The 0.7–1.3 %
  raw gain of our snapshot (the polarization-split arms through the plate splitter)
  is absent by construction — COPHI reads at the gain our piezo leg reads at, 1.0000.
- **Exact, linear phase steps.**  The step is 90° set by two synthesizers on a Rb
  standard with the camera trigger derived from the same clock; there is no piezo
  calibration, nonlinearity or hysteresis, and no moving part.  The classic
  phase-step-error ripple (second harmonic of the fringe phase) is reduced to the
  timing jitter of the trigger.
- **Continuous.**  A map every 133 ms with no scan to command; the stepped readings'
  "mirror drifting across the scan" effects are the same size but the cadence is set
  by the beat, not by a mechanism.
- **Intensity-independent per pixel**, as every four-step is: DM reflectivity
  variations and the LO's amplitude profile do not enter the phase (unlike the
  scalar Zernike's amplitude reference, the 0.78 gain in our deck).

What it adds, as error terms we have not modeled:

- **The LO's own wavefront and its drift.**  Our reference arm is a flat; COPHI's is
  a fiber tip and a collimator, whose aberration is static (calibrated with the DM at
  null) but whose tilt/defocus drift with the fiber tip's position and the
  collimator's temperature.  Piston is immune (removed), but a slowly tilting or
  defocusing LO is a low-order error the actuator fit cannot distinguish from the DM.
  The same class as "a tilting reference is outside the model" in the P/SRI row — a
  row to add.
- **Within-measurement path drift between the two fibers.**  Differential fiber phase
  noise (acoustic, thermal) over the 133 ms is a phase-step error common to every
  pixel; with the fringe nulled (the 30 nm surface is 0.1 wave of wavefront) the
  resulting sin(2φ) ripple is nearly uniform over the pupil, i.e. mostly piston, with
  a residual proportional to the surface — a gain-like term.  Cheap to model with
  the existing intra-scan machinery; PM fiber in vacuum made it negligible for TOM3.
- **Timing jitter of the 30 Hz trigger** relative to the beat — the same term, from
  the electronics side.  A 5-frame (Hariharan) algorithm at 37.5 Hz would make both
  second order if either turns out to matter.

## 3. The number that decides it: the camera, not the heterodyne

The heterodyne buys COPHI nothing in photons.  At 30 Hz the photon budget per second
is pixels × well depth × frame rate.  With the deck's convention (about 1e6 lit
pixels at 60 ke full well, the "1e14 photons is 1600 frames" line):

- 6e10 photons per frame, **1.9e12 per second**.
- A 1 Hz reading at full well is 1.9e12 photons: scaling the deck's 2.3 pm at 1e14
  (lens rig; 2.1 pm mirror rig) as 1/√N gives **15–17 pm per reading at 1 Hz**,
  photon-limited, before any bench error.
- 1 pm needs 7e14–1e15 photons: **6–9 minutes** of 30 Hz frames.
- Holding 3 pm needs 2.0–2.3e13 photons per cycle: **11–12 s per cycle**, a 0.1 Hz
  servo; at a 1 Hz servo the hold floor is about 10 pm.

These are the same numbers as for our interferometer with the same camera — the
deck's point that the camera's well depth, not the laser, sets the time.  TOM3's
40 pm per zero crossing at 10 kHz and sub-pm in seconds came from PHOTODIODES: no
well-depth ceiling, mW on each of 43 detectors, 1e4 phase samples per second.  So:

- **COPHI-with-a-camera at 1 Hz is a 15–20 pm gauge** on a 1k × 1k, 60 ke, 30 Hz
  camera; a deeper well or a faster camera scales it as √(well × rate) — a 100 ke
  sCMOS at 100 Hz would bring 1 Hz readings to about 7 pm and 1 pm to about a minute.
- **COPHI-with-photodiodes is the version that is fast**, and it is worth a look
  for THIS problem specifically: the unknown is the 96 × 96 actuator lattice (9216
  values, ~5000–7000 lit), not 1e6 pixels — and that is the trade.  Photodiodes with
  a phasemeter at 10 kHz give a pm-class reading per second per channel with no
  well limit (633 nm × 10 kHz / 128 MHz = 49 pm per crossing, 0.5 pm per second if
  white), but one channel per actuator is thousands of phasemeter channels (TOM3
  multiplexed 43 into 10); a coarser array trades resolution for speed, and every
  element averages the phase over its cell — a fixed, known linear operator the
  response matrix absorbs, and exactly where Fang's blur question becomes concrete
  (a coarse sampler IS the blur).  Speed versus resolution is the axis to model:
  the camera at the pixel end, the phasemeter array at the other, and what the
  actuator fit recovers from each.  (Dave 2026-10-09: this deserves careful
  modeling on its own — a campaign of its own when Feng's deck arrives: a template
  `cophi_dm96` beside the three benches, `deck_cophi`, then a composite deck.)

## 4. What the model can say quickly (if Dave wants a row in the comparison)

The bench and the reading exist; the increment is small:

1. **A `COPHI` reading in `tg96_run`**: the existing four-step with (a) the reference
   arm replaced by an LO source at the recombiner with its own wavefront, (b) the
   step from the beat with a trigger-jitter parameter, (c) the quarter-cycle
   integration, (d) differential fiber phase as an intra-scan drift.  Scored on the
   same 30 nm surface with the matrix on it, the same photon rows, the same servo
   rows — one column beside the interferometer's.  Expected: gain 1.0000 like the
   piezo leg, the same photons per pm, and the LO-drift row as the new number.
2. **The camera budget line** on the photons slide: pm per reading at 1 Hz vs well
   depth and frame rate, for every gauge in the deck (it is the same line for all).
3. **A photodiode-array variant** as a reading class: cell-averaged phase at N sites,
   phasemeter noise per sample, 10 kHz; its matrix on the surface; its blur operator
   printed next to CCMac's blur-kernel demo.

Not modeled and not claimed: the AOM/RF chain, fiber polarization drift in PM fiber,
vacuum.  TOM3's results are in a 3.3 m vacuum chamber on a Zerodur sub-bench.

## 5. The rest of the follow-up list, placed

- **1. Thorough error analysis per candidate** — the deck's backup "From this model to
  as-built performance" is the skeleton (surface, alignment, stability, polarization
  errors priced per reading); COPHI joins it with the LO-drift and timing terms above.
- **2. Pupil-imaging blur demo** — CCMac (in progress): convolve a blur kernel with a
  DM surface, reconstruct the actuator state, find the kernel width where it fails.
  The photodiode-array COPHI is the extreme case of the same question.
- **3. Kent's ZWFS reconstruction experience** — fits the "What the Zernike-sensor
  literature offers" backup (our first import was re-computing the reference wave);
  ask for his reconstructor and score it on the same bench.
- **4. COPHI** — this note; Feng's slides when they arrive.

## Addendum, 2026-10-09 later: Feng's full deck

`MACOS_sandbox/DM_gauge/COPHI_for_DM_gauge.pptx` (12 slides, 2026-05-19) adds the
PURPOSE — a DM drift sensor for the HWO coronagraph: DMs drift open-loop
(material physics, "floating electrodes", ~2 % temperature coefficient), "the
tallest tent pole in contrast stability"; COPHI runs each DM in a local loop —
i.e. our HOLD job; the Zernike WFS with a dedicated laser (Ruane 2020, DST,
Fang's VSG plan) as the alternative concept; the four-frame algebra and the
generic N-frame least squares (Zygo's four-step = N 4) as backup; the surface
H = ½ (φ/2π) λ (double pass).  No camera spec beyond 30 fps, no photon or error
budget.  The campaign plan is `BRIEF_cophi_campaign.md`.
