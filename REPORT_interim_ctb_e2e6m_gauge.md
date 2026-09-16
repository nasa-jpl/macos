# Interim report: the coronagraph testbed model, the 6 m end-to-end model, and the DM surface gauge

DRAFT skeleton, CCL for Dave, started 2026-09-15; due Friday 2026-09-18.
Audience: the JPL HWO WFS&C group (the CTB story Dave sends out).  Every
number will carry its run tag; the three source decks are
`bench_ctb/deck_ctb.md` (v5, 15 + 9 backup), `e2e6m_r2/deck_e2e6m_r2.md`
(19 + 9), `demo_session/deck_gauges.md` (45, DRAFT).  Style:
`doc/STYLE_REPORTS.md` section 5 gate before the .docx/.pptx.

## Plan to Friday

| day | what lands |
|---|---|
| Tue 09-15 | this skeleton; the CTB DM-model fix ordered (TO, `BRIEF_to_gauge_close` 0b); deck_ctb DRAFT banner |
| Wed 09-16 | sections 1-2 drafted from the decks and reports (CTB, e2e6m, the overlap); section 3 from the gauge deck |
| Thu 09-17 | TO's re-scored CTB numbers in (noon); section 1's DM/EFC numbers replaced; section 4 (the gauge inside the coronagraph) with the corrected separability; figures chosen (the tools' own) |
| Fri 09-18 | style gate; .docx; Dave's read; push on his word |

## 1. The coronagraph testbed model (CTB)

- 1.1 What it is: the 8-OAP / 2-DM bench in two layers (geometry: `example_ctb`, staged optimization; diffraction: `ctb_prop_layout`, the compact and the station-to-station decks, the exit-pupil-sphere quartet at every mask plane).  Validation against PROPER (pitch ratio 1.0000 / corr 1.000000 / centroid 0.000 px on the through-focus leg).
- 1.2 Coronagraph performance: mask sweep, the 2.70 lam/D null, Lyot 0.50; the six mask families head-to-head on one annulus (APLC 2.1e-10 at 27%, BLC 2.7e-8, hard 2.5e-7, vortex-matched 2.9e-7, R&R 3.2e-6, dual-zone 6.4e-6); the vortex against the Lyot stop.
- 1.3 The loop: EFC on the engine itself (hard 2.9e-7 -> 8.1e-9; vortex 1.7e-8 -> 6.8e-15); polarization (1.1e-15 residual, a state change); bandwidth (mono 8e-13 -> 20% 5.4e-11); the vector vortex verdict (leak uncorrectable but optically removable).  **The DM model was questioned (2026-09-15) and confirmed by the traced footprint (09-16): the beam at the DMs is 21.24 mm and the 32 x 32 / 0.67 mm lattice spans it; every loop number stands.  What the probe found instead: the generator reads the engine's point-source Aperture (a full cone angle) as a half-angle, so the bench carries half the sheet's intended fill (47% of the DM's clear aperture, not 95%) -- self-consistent, and Dave's call whether the story keeps the bench as built (recommended) or regenerates at the intended beam.**
- 1.4 Hand-offs: the phase-factor export (18 stations, PROPER consumes it; per-leg replay check), the pure-PROPER run.
- 1.5 Open, in dependency order (deck slide 14): as-built surfaces, time-series drift, FALCO as the DM driver, validation against a dataset.

## 2. The 6 m end-to-end model (e2e6m, round 2) and its overlap with the CTB

- 2.1 What it is: a diffraction-limited unobscured 6 m telescope (0.0473 waves rms across the field at 500 nm), 19-segment primary with physical apertures, the 8-mirror / 2-DM relay in metres spliced on, six coronagraph families on one train, an imager on the same shroud; sensitivities (dwdx / dwdz / dwdgrid), the error budget closing engine vs model to 0.35%, the metrology truss (114 gauges + 252 edge sensors), the restart ladder, JWST-class drift held at 2.0e-9 in closed loop.
- 2.2 The overlap with the CTB, stated as what is SHARED and what is NOT:
  - shared: the Bench primitives (add_oap, the DM as a grid surface with influence functions, add_reference markers), the propagation recipe (the sphere-bracketed mask quartet, NF1/NF2), the mask library and its generators, the EFC driver, the contrast scorer, the phase-export format;
  - not shared: the pupil (segmented 6 m vs a 45 mm circular DM stop), the DM pitch (1.48 mm on a 47.5 mm beam vs 1.34 mm on 42.75 -- both 32 across), the packaging (a 3-D shroud vs a planar table), the drift model (telescope + segments vs bench thermal), and the DM-model slip (the e2e6m DMs are sized from their beam; the CTB's were not -- check and state).
  - what the CTB validates for e2e6m: the propagation and mask machinery against PROPER; what e2e6m adds: the telescope, the segments, the truss, the time series.
- 2.3 What each is for: the CTB is the lab-facing model (as-built data, phase export to external users); e2e6m is the mission-facing one (error budget, drift, hold).

## 3. The DM surface gauge

### 3.1 The question and the bench

The gauge has three jobs on a 96 x 96 deformable mirror (1 mm pitch, 96 mm
across): measure its surface to picometers, capture its post-launch shape
(100 to 200 nm of wavefront) into a servo's reach, and hold it there.  Every
test runs on a random 30 nm rms working surface with the response matrix
measured on that surface, because a flight gauge never operates at null.
Scores are the recovered change divided by the true change (gain), the rms
of the actuators that were not changed (floor, pm), the photons per
measurement summed over all frames, the largest surface at which a 10 nm
change still reads within 10% (capture range), and the photons per cycle
that hold 3 pm rms in a 60-cycle loop at gain 0.5 (hold).  Source:
`demo_session/deck_gauges.md` slides 2, 18-23 (run tags on each slide).

One bench carries four ways of sensing: a filtered HeNe at 632.8 nm, a
collimator, a plate splitter at 22.5 degrees, a 700 mm leg to the mirror,
a focuser to an internal focus at F/4.2, a field lens and a camera at the
pupil image.  The bench of record moved from the model's 7-degree splitter,
which is not buildable (eight of nine node parts sat in another beam), to
22.5 degrees on 2026-09-15; every part clears by at least 38 mm and the
sensors' rows did not move (run tag gate22_193).  The collimator and focuser
are either 103 mm lenses or off-axis parabolas; both layouts are traced by
the engine and drawn by the tools.  Real thicknesses are in the model
(10 mm splitter and compensator, 4 mm lens edges, 2-3 mm substrates on
every polarizing element and mask): the sensors read the same numbers to
the digit (sub22, thk22); the interferometer's flat-mirror null rises from
0.13 to 20 nm because a 10 mm plate at 22.5 degrees shears the transmitted
beam 1.39 mm, and its differential rows under that are being re-measured
(the Mac run item4bseq).

### 3.2 The four approaches and the record

| reading (frames per measurement) | 10 nm on one actuator: gain / floor | 1 nm on 47 sites | dense random 10 nm | hold, 3 pm under a 2 pm per-cycle walk |
|---|---|---|---|---|
| Twyman-Green interferometer, lens rig (4) | 0.991 / 2 pm | 0.990 / 1 pm | 0.990 / 168 pm | 2.0e13 photons per cycle |
| interferometer, redesigned mirror rig (4) | 0.992 / 2.4 pm | -- | 0.991 / 161 pm | 1.7e13 |
| Zernike dimple, stepped (4) | 0.989 / 5 pm | 0.999 / 4 pm | 0.984 / 0.68 nm | 7.5e12 |
| vector Zernike, polarized dimple (2, at once) | 0.9935 / 4 pm | 0.9992 / 3 pm | 0.9999 / 0.33 nm | 5.3e12 |
| stepped pinhole with a shutter frame (6) | 0.9935 / 4 pm | 0.9992 / 3 pm | 0.9985 / 338 pm | 7.0e12 |
| P/SRI, reference arm traced (4) | 0.9935 / 2 pm | 0.9924 / 1 pm | 0.9926 / 144 pm | 1.5e13 |

Run tags lens_deck, oapifo2, oapifol2, matbase, v193base, pdi193fbase,
pfdeck, vloop193, ploop193, loop193.  The four best readings are within 1%
of each other and within 2 pm on the floor.  The interferometer's dense
error (161-168 pm against 0.33-0.68 nm) is its reference arm: it does not
make its reference from the beam.  The one-frame Zernike readings lose (the
linear one is biased 4-5% and floors at 21 pm; the exact one misreads the
sign past the quarter-wave fold and diverges in a servo).

Capture range separates the readings by what they reference.  With the
calibration left as measured at 30 nm, the self-referenced sensors hold 10%
accuracy only to 36-70 nm of surface, because their reference wave is made
from the beam and collapses as the surface grows; the interferometer and the
pinhole with its shutter frame hold to 480 nm and beyond (cap385, cap385p,
lens_deck).  Re-measuring the matrix on the surface restores the sensors'
gain to 160 nm at the price of light, up to forty times the 1 pm photon cost
(cap385_b60-b160).  Capturing the initial figure is a wrapping problem
before it is a noise problem: every phase reading wraps at +/-158 nm of
surface.  With unwrapping alone the interferometer takes 300 nm of surface
to 3 pm in 43 cycles on the lens rig (descent_lens) and 200 nm in 19 on the
redesigned mirror rig (oapdesc2); the pinhole with its shutter frame and the
P/SRI capture 100 nm with unwrapping and a matrix re-measured every 10
cycles; the self-referenced sensors capture only from 60 nm.

The servo holds 3 pm from 1.5e12 photons per cycle (vector Zernike, noise
only) to 5.5e12 (interferometer), 5.3e12 to 2.0e13 under a 2 pm per-cycle
random walk, and every reading with no fixed error sits at the same 10 pm
floor under a 5 pm per-cycle thermal ramp.  That floor is the loop's lag,
rate divided by gain, not a light budget: on the mirror rig 9.2 of its
10 pm lies below 4 cycles per aperture, so a low-order loop at higher gain
removes it (oapifol2).  Both drift floors are the textbook formulas and the
engine reproduces them to 1%.

### 3.3 Systematics priced, and what is not yet modeled

Every systematic that could be modeled was either calibrated by the matrix
measured through the sensor or priced in light (deck slide 27): the
metasurface's retardance error (750 pm uncalibrated, the ideal rows after
a three-number calibration), the arm's polarization aberration (6 nm per
radian of channel phase; nothing through the matrix to 0.1 rad), the
analyzer's leak (3.5 pm), piezo step error (common mode in the servo),
camera drift within a scan (helps the stepped readings when it is the
mirror, hurts when it is the camera), the P/SRI's own reference-arm walk
(a piston, which the estimator nulls).  On the redesigned mirror rig the
vector reading's raw fold-gate error is 634 pm on bare aluminum, from a
1.08 amplitude imbalance between its two circular channels; a quarter-wave
overcoat halves it and calibrating the two channels removes it (0.054 pm,
tags oapsens22, vqw22, vmap22).  The snapshot interferometer's polarization
at the built 22.5-degree splitter rotates the test arm 1.56 degrees and
biases the gain 0.9% (7.5 degrees and 11.7% at the record's 45 degrees),
with 1.4e-3 left after a measured matrix (aoi_lens22, aoi_oap22).

Not yet modeled, stated on the slides: the camera (a 6.5 um sCMOS binned 3
to the 385 modeled pixels across the traced 7.5 mm pupil image; well depth
sets the measurement time, not the laser), the photonic phase shifter of
the two-arm pinhole form (its phase steps are ideal increments with a
step-error knob only), as-built optics, and a second color for capture
beyond one wave.

### 3.4 The reflective front end, and the recommendation as it stands

The record's off-axis-parabola rig lost the interferometer's servo and
broke the focus-critical sensors, and the cause was read as fold coma.  It
was a 25 mm conjugate error: the builder fed the collimator 25 mm inside
its focus (926 urad residual, blur linear in the fold angle), which a lens
hides in its tuned figures and a parabola cannot; the input polarizer also
sat 10 mm past the collimator's pole, inside the incoming cone.  Fed at its
focus and re-solved at the 22.5-degree splitter (OAP1 20 degrees, OAP2 25
degrees), the mirror rig clears every part by 33 mm, focuses on the mask
seat with no adjustment, nulls 0.029 nm with a flat mirror against the lens
rig's 0.13, and reads, holds and captures at the lens rig's numbers (the
table above; cross-talk 0.006-0.02 against the record's 0.22-0.48).  The
pinhole recovers (0.27 pm against 94), the vector pair passes calibrated.
Tags fold1-4, conj, zseat, oap22d, oapifo2, oapifol2, oapdesc2, oapsens22.

The recommendation as the deck states it: hold with the vector Zernike
sensor (3 pm from 1.5e12 photons per cycle, two frames at once, no fold,
no fixed error, every modeled systematic calibrated through it); capture
with the interferometer, or with the stepped pinhole's shutter frame if a
second arm is not wanted; build one bench that switches (the mask seat
translates, the reference arm is shuttered, the quarter-wave plate goes in
or out).  The front end is either: the mirror rig matches the lens rig on
every measured line and buys any color and no glass in the beam; it costs
alignment tolerance (10 urad of mirror tilt moves the null 95-103 nm) and a
polarization calibration of the vector reading.  The record's stance was
lenses; the measured trade awaits a ruling.

## 4. The gauge inside the coronagraph

`NOTES_gauge_in_coronagraph.md` and the clearance tool
`demo_session/gauge_in_coro_clearance.py`, which walks the chief ray
through each committed deck so the off-axis parabolas sit at their poles,
then scans a gauge beam about each DM's normal against every body and
science beam.

### 4.1 Packaging

A face-on gauge is blocked in both packages: the coronagraph's DM fold
(5-6 degrees) puts OAP1's body and the other DM on the normal (CTB: 54 mm
of interference at the beam, 9 mm at DM2; the flight relay: 72 mm).  The
CTB is planar, so above and below the table is empty: a gauge in from
above at 12 degrees and out below clears every body by 26 mm with its
first optic 300 mm from the DM (15 degrees: 40 mm).  In the flight relay no
angle to 30 degrees clears, because the OAPs sit 200-250 mm from the DMs;
an external gauge there is a 25-35 degree periscope the relay would be
designed around.

### 4.2 The options

Seven were listed: a face-on gauge (needs the DM fold opened and a long
DM-to-OAP leg); an out-of-plane oblique ZWFS or pinhole gauge, one per DM,
each measuring its mirror alone (fits the CTB; sensitivity 2 cos theta and
a 3.5% elliptical footprint, both absorbed by the on-surface matrix); a
gauge injected at the source and read at the first focus after the DMs; the
coronagraph's own low-order sensor and probing; a dichroic in the DM leg
(ground only); one periscope head time-shared; and opening the DM fold as a
relay design knob.

The direction Dave set on 2026-09-15 is the third option taken to its end:
a dichroic at the apodizer pupil pulls out-of-band light to a vector
Zernike or pinhole reading of the coronagraph's input field, amplitude and
phase, and the two DMs servo that field to the one recorded when the dark
hole was dug.  That holds the hole against everything upstream of the
pickoff, telescope misalignment and figure drift included, not only the
DMs.  The gauge deck's hold reading and currencies transfer to it directly.

### 4.3 The bounds

Four numbers bound it.  What the pickoff cannot see: the apodizer, FPM,
Lyot and their alignment stay with the dark-hole probing, which then runs
far less often.  Two DMs from one conjugate: amplitude and phase are two
real maps and two DMs are two real maps, so the control is determined, but
DM2's amplitude authority at DM1's plane is the Fresnel conversion
sin(pi lambda z / L^2), weak at low order (DM2 stroke per unit amplitude
correction on the CTB: 130x at 2 cycles across the beam, 21x at 5, 5x at
10, 2x at the actuator Nyquist), so low-order amplitude drift is expensive
and the servo regularizes as EFC does.  Chromaticity: surface phase
transfers out of band as OPD; amplitude made by out-of-pupil phase scales
with wavelength, so the servo holds a model-propagated target.  Photons on
a star, from the deck's vector-pair servo cost (1.5e12 photons per cycle
for 3 pm), a 6 m aperture, a 100 nm out-of-band slice and 25% throughput:
10 pm every 12 s on a V=2 star, every 200 s on V=5, every 52 min on V=8.

Separability of the two DMs from one conjugate, on the traced beams: on the
CTB (0.67 mm pitch, 21.2 mm beam, 500 mm apart) DM2's figure converts 47%
to amplitude at the actuator Nyquist and 93% at one cycle per actuator, so
the two mirrors separate at the top of the controllable band; in the flight
relay (1.48 mm pitch, 47.5 mm beam, 400 mm apart) the conversion is 7% at
Nyquist and the two DMs do not separate in band.  Per-DM knowledge there
needs a second conjugate or differential attribution through the
multiplexed matrix.

### 4.4 What is queued

The field servo on the CTB deck in three steps (BRIEF_to_gauge_close item
7; TO's REPORT_field_servo.md): the reading at the apodizer conjugate,
gated to reproduce the engine's complex field to the vector pair's record
(0.016 / 0.128 pm through 5 / 20% amplitude dips); the multiplexed matrix
over both DMs' actuators and the measured DM1-DM2 cross-talk against the
predicted 16.5-cycle crossover; the servo under an upstream drift (a 10
urad tilt of OAP1) with the dark-hole contrast scored through the EFC chain
before, after, and after the hold.

## 5. Where the three meet, and what is next

One engine, one Bench builder, one DM doctrine, one propagation recipe, one contrast scorer; the gauge's readings become the coronagraph's input-field sensor; e2e6m carries the telescope side.  Next: the CTB re-score and its aberration arc; the field-servo model; the gauge deck's sign-off; migration of the gauge work into templates.

## Sources
deck_ctb.md / CTB_PROP_STATUS.md / README (bench_ctb); deck_e2e6m_r2.md / e2e6m_r2_LOG.md / README (e2e6m, e2e6m_r2); deck_gauges.md / BRIEF_gauge_deck.md / the three lane reports (REPORT_gauge_ifo, REPORT_reflective, REPORT_gauge_pdi, zwfs_dm96 README); NOTES_gauge_in_coronagraph.md; memory project_ctb_diffraction / project_e2e6m / project_tg96_gauge.
