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

### 1.1 What it is

The CTB model is an all-reflective coronagraph bench: eight off-axis parabolas (OAPs), two deformable mirrors (DMs), an apodizer and a Lyot stop at pupil images, a focal-plane mask and a field stop at focus images, and a camera (CTB slide 1).  DM1 is the aperture stop; the beam on the DMs is 21.2 mm in diameter (slide 1; see 1.3).  The bench is geometrically diffraction-limited at 0.0014 waves (slide 1).

The model has two layers over the same optics (bench_ctb/README.md, layer table).  The geometry layer, `example_ctb.m`, places the optics by staged optimization: each pupil or focus is solved on its own plane in light order with the upstream optics frozen, so no conjugate can trade against another (README, "Why staged optimization").  The diffraction layer, `ctb_prop_layout.m`, emits two propagation prescriptions (slide 1): a compact model of 31 elements (one plane-to-plane leg DM1 to DM2, a four-surface mask block at each focus, far field to the camera) and a station-to-station model of 44 elements in which every leg between optics is propagated.  The four-surface block at each mask plane is a flat return, an exit-pupil sphere carrying the first half-propagation, the mask plane carrying the second half, and the same sphere again; the two sphere distances agree to all digits, so the block is transparent when no mask is applied (slide 1).

Validation is against MATLAB PROPER (Krist 2007).  The through-focus half-propagation from the exit-pupil sphere to the mask plane, the one step no earlier PROPER campaign covered, reproduces PROPER at matched sampling: focal pixel pitch ratio 1.0000, peak-normalized correlation 1.000000, centroid offset 0.000 px (slide 2).  The compact and station-to-station models agree on the bare image to 0.9989 correlation; with the coronagraph in, they differ by 1.76x, the mirror-to-mirror diffraction the compact model omits (slide 2).

### 1.2 Coronagraph performance

With a hard occulter, centering the mask builders on the focus pixel, sampling the camera at 4.0 pixels per λ/D (model 1024) and a mask sweep bought 16x: an interior null at an occulter radius of 2.70 λ/D with a Lyot stop of 0.50 (25% throughput) took the dark-zone mean from 4.6e-6 to 2.9e-7 (slide 3).

Six mask families were then scored on one grid, one annulus and one normalization (slide 5).  Static, before any DM control: apodized-pupil Lyot 2.1e-10 at 27% throughput; vortex with matched Lyot 1.4e-8 at 81%; band-limited (4th order) 2.7e-8 at 36%; hard occulter 2.5e-7 at 25%; Roddier π-mask 2.4e-6 at 81%; dual-zone phase 6.8e-6 at 81% (slide 5 table).  (The vortex's earlier 2.9e-7 was a sampling artifact of its singular core, cured by pixel-averaging; backup, "The vortex core".)  With the Lyot fraction as the only dial, a charge-4 vortex under-runs every fixed design at every throughput: 8.8e-11 at a 0.60 stop (36%) against the apodized Lyot's 2.1e-10 at 27%, with no apodizer to fabricate (slide 6).

### 1.3 The loop

Every number here is measured on the design train: no as-built surface errors, no drift, and sensing that reads the model's complex field directly (slide 9 footnote; slide 14).  Each DM is a 32 x 32 actuator lattice at 0.67 mm pitch (the 21.2 mm beam / 32) with Gaussian influence functions, 880 actuators inside the beam (slide 9).  The control matrix is measured: 1760 pokes of 2 nm through the full masked chain, 11 minutes (slide 9).  Electric-field conjugation (EFC, Give'on 2007) takes the hard-occulter chain 36x in 19 iterations, within 2x of the matrix's linear bound of 4.5e-9; DM1 alone stalls at 1.3e-7 (slide 9).  Re-measuring the matrix about the corrected state moves the hard chain to a floor set by occulter-edge diffraction and the vortex chain (charge 4, Lyot 0.60) to double-precision roundoff at half-nanometer strokes (slide 10; Table 1).

Polarization does not set the floor.  With quarter-wave MgF2 over aluminum on all ten mirrors, the ray-traced Jones pupil has 3.4 mrad mean retardance with 8.7 µrad rms variation; the mean is a state change the loop absorbs.  Closed loop reaches 5.8e-13, the scalar result, with an uncontrollable residual of 1.1e-15 at every bandwidth (slide 11).  Bandwidth is the floor that grows, 65x from monochromatic to a 20% band, with control wavelengths at a constant 2.5% spacing (slide 12; Table 1).  The vector vortex verdict (slide 13): the zero-order plate's retardance leak (1e-3 of the starlight at 5%, 1e-2 at 20%) is uncorrectable by the DMs, since its sign flips across band center, but a circular-polarizer sandwich removes it optically and per-wavelength control returns the sandwich to the scalar floor (7.9e-12 at 5%); a crossed-linear sandwich goes two to three decades deeper at the price of eight planet blind spots.

The DM model.  Questioned on 2026-09-15 (was `ctb_dm.m`'s 21.3 mm a radius used as a diameter?), it was confirmed on 2026-09-16 by a traced footprint: `ctb_beam_probe` on `ctb_dcr.in` (model 512, 50618 rays) measures 21.2444 mm at DM1 and 21.2451 mm at DM2, so the lattice spans the beam and every loop number stands (REPORT_field_servo.md section 0; CTB_PROP_STATUS.md, 2026-09-16).  The probe found a different slip: the generator sets the source cone as a half-angle where the engine takes the full cone angle (`sourcsub.F`, `A = Aperture/2`), so the intended 42.75 mm beam lands as 21.2 mm and the DM's 22.5 mm clear radius is filled to 47%, not 95% (README, "Source model").  Everything downstream is sized on the beam that exists; regeneration at the intended beam is a pending decision (README, 2026-09-16 note).

### 1.4 Hand-offs

One self-describing export carries the 44-element model: 18 stations (complex field, amplitude, OPD, pitch, in meters), 17 legs, 4 reference spheres and 18 phase screens, conventions stamped inside and orientation measured (slide 7).  Replayed leg by leg in PROPER, focus stations reproduce at correlation 1.000000 and gated pupils at 0.9998 or better (backup, "Per-leg replay check").  A pure-PROPER script reading only the export reproduces the bare image at 1.000000 and a dark zone at 1.4e-8, within the 2x gate above the shipped 2.9e-7 (slide 8).  A single continuous PROPER beam cannot reproduce the model (pitch ratio 0.71): every intermediate focus is sampled at the exit-pupil pitch (slide 8).

### 1.5 Open

In dependency order (slide 14): realistic sensing (pairwise probing from camera images, then FALCO on the same DMs); as-built surface maps, then drifts as a time series against the delivered loop; validation against a named testbed dataset.  `ctb_study` reruns slides 9 to 13 at other parameters with one call (backup).

Table 1.  Closed-loop results per chain (CTB slide 10) and the bandwidth ladder on the vortex chain (CTB slide 12).  N = 512, design train, perfect sensing; strokes at the floor 9.9 / 8.6 nm rms (hard chain), 0.49 / 0.55 nm rms (vortex).

| chain | band (control colors) | static | fixed matrix | re-measured matrix | polarization floor | source |
|---|---|---|---|---|---|---|
| hard occulter, Lyot 0.50 | mono (1) | 2.9e-7 | 8.1e-9 | 3.8e-9 | not run | slide 10 |
| vortex charge 4, Lyot 0.60 | mono (1) | 1.7e-8 | 5.8e-13 | 6.8e-15 | 1.1e-15 | slides 10, 11 |
| vortex charge 4, Lyot 0.60 | 5% (3) | not stated | 9.4e-12 | not run | 1.1e-15 | slide 12 |
| vortex charge 4, Lyot 0.60 | 10% (5) | 2.2e-8 | 2.5e-11 | 2.0e-11 | 1.1e-15 | slides 5, 12 |
| vortex charge 4, Lyot 0.60 | 20% (9) | not stated | 5.4e-11 | not run | 1.1e-15 | slide 12 |

## 2. The 6 m end-to-end model (e2e6m, round 2) and its overlap with the CTB

### 2.1 What it is

Round 2 of e2e6m is one model from mirror figure to a held dark zone (e2e6m slide 19).  The telescope is a 6 m unobscured three-mirror design: 0.0473 waves rms worst case over a ±0.35 arcmin field at 500 nm against the 0.071-wave diffraction limit; 7.450 m in an 8 m shroud; f/25.39 against a requested f/12 to 20 (slide 1).  The primary is 19 hexagonal segments, 1.2 m flat to flat with 25 mm gaps, each with its own polygonal aperture; 983 of 985 rays survive, and a 10 nm displacement of one segment moves the wavefront 19.91 nm over that segment's 52 rays and zero over the other 930 (slide 2).

The relay is the CTB topology in meters: an OAP collimator to a 47 mm pupil, seven 1:1 relays, DM1 and DM2 0.15 m apart at the collimated pupil, spliced onto the telescope as 37 elements of one prescription (slide 3; LOG, 2026-08-26 R1).  A deployable pick-off feeds an imager at 0.0042 waves rms; both legs fit the shroud at 7.451 m (slide 4).  Six coronagraph families run behind a circular stop at the apodizer plane; pre-control the apodized Lyot leads at 4.37e-7 with 9% throughput, and the vortex pays the segment gaps about 40x (slide 9).  The gap cost against a monolithic twin is 1298x (slide 8); closed loop the apodized Lyot reaches 1.1e-7, within 2x of its linear bound (slide 10).

The sensitivity model has 192 rigid-body, 152 figure and 114 influence channels over five fields; engine versus model closes to 0.35% worst over 18 freedom pairs (slide 14).  The metrology truss, 114 laser gauges and 252 edge sensors, leaves 0.86 nm of wavefront per nm of gauge noise; its finite-difference check closes at 0% (slide 15).  The restart ladder (EFC, relinearize, restart) takes the 1.10 m DM-spacing train from 1.2e-6 to 1.13e-9 in 10 rounds over 4.6 h, with the linear-achievable substrate at 2.0e-11 to 3.8e-11 (slide 16).  Under the JWST-class drift of Table 2 the open loop decays 4.5 decades in a day; segment control (BLUE + ridge, gain 0.5) plus a guarded EFC hold of one damped step per hour keeps 2.0e-9 all day (slide 18).

### 2.2 The overlap with the CTB

Shared, code and recipe:

- the Bench primitives `add_oap`, `add_mirror` and `add_reference` (bench_ctb/README.md, topology; LOG R1);
- the DM as an influence-function grid surface, `ctb_dm` and `ctb_dm_rx` (CTB slide 9; e2e6m slide 7, the same 20 nm poke gate);
- the sphere-bracketed mask quartet and the NF1/NF2 propagation recipe, seeded on the DM1-to-DM2 leg (CTB slide 1; LOG R1);
- the mask library and its generators (CTB slide 5; e2e6m slide 9);
- the EFC driver and the contrast scorer, `ctb_chain` and `ctb_efc` pointed at the e2e6m deck (LOG R3);
- the phase-export format (CTB slide 7).

Not shared:

- the pupil: a segmented 6 m primary against a 45 mm circular DM stop carrying a 21 mm beam (Table 2);
- the DM pitch: 1.49 mm on a 47.5 mm beam against 0.67 mm on 21.2 mm, both 32 across (Table 2);
- the packaging: a deployed 3-D shroud fit (e2e6m slide 4) against a planar table (CTB backup, "Bare-optics agreement");
- the drift model: a telescope-plus-segment time series (e2e6m slide 18) against none yet on the CTB (CTB slide 14).

The DM-sizing question was checked.  The e2e6m lattice is sized from its traced beam: `r1_dm.m` line 57 sets `beam_d = 2 * 0.023771`, commented "measured pupil at the DMs (r1 gate)", the radius `r1_seg_report.txt` records at DM1 and DM2 on a 60 mm clear aperture.  The CTB lattice was also sized from a traced footprint (`ctb_dm.m` default 21.3 mm; REPORT_field_servo.md section 0).  The difference is not in the lattices: the CTB beam is half what its generator intended (1.3); the e2e6m beam is the 47 mm the design specifies.

What the CTB validates for e2e6m is the propagation and mask machinery against PROPER (1.1, 1.4); what e2e6m adds is the telescope, the segments, the truss and the time series.

### 2.3 What each model is for

The CTB is the lab-facing model: cross-checked against PROPER, its next inputs are as-built surface maps and a testbed dataset (CTB slide 14), and its phase export lets an external PROPER user run the same planes with no macos (CTB slides 7 and 8).  The e2e6m model is the mission-facing one: the error budget (e2e6m slide 14), the metrology (slide 15), the drift and the hold (slide 18), what a flight-like train does rather than what a bench measures.

Table 2.  The two models side by side.

| parameter | CTB | e2e6m round 2 |
|---|---|---|
| aperture | 45 mm circular DM stop (22.5 mm clear radius), DM1 is the stop (CTB slide 1; README "Source model") | 6 m, 19 hexagonal segments, 1.2 m flat to flat, 25 mm gaps (e2e6m slide 2) |
| beam at the DMs | 21.2 mm diameter, 47% of the DM radius (CTB slide 1; README "Source model") | 47.5 mm diameter on a 60 mm DM, traced (e2e6m slide 7; `r1_seg_report.txt`) |
| DM pitch, actuators across | 0.67 mm, 32 x 32, 880 in the beam (CTB slide 9) | 1.49 mm, 32 x 32, 880 in the beam (`r1_dm.m`; e2e6m slide 7) |
| mask families run | six: apodized Lyot, vortex, band-limited, hard occulter, Roddier, dual-zone (CTB slide 5) | six: classical Lyot, apodized Lyot, APLC-as-implemented, band-limited, vortex charge 4 and 6 (e2e6m slide 9) |
| deepest contrast | 6.8e-15, vortex charge 4, mono, re-measured matrix, design train (CTB slide 10) | 1.13e-9, apodized Lyot, restart ladder at 1.10 m spacing (e2e6m slide 16) |
| control loop | EFC on a measured matrix, perfect sensing, relinearize once (CTB slides 9 and 10) | EFC restart ladder plus segment control (BLUE + ridge) and a guarded EFC hold (e2e6m slides 16 and 18) |
| drift model | none yet; queued after sensing (CTB slide 14) | 10 nm/hr correlated ramp + 0.5 nm/step walk, 24 h at 30-minute frames; open loop 1.1e-9 to 3.9e-5, held 2.0e-9 (e2e6m slide 18) |

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
| the three sensors on the redesigned mirror rig | -- | -- | -- | vector 5.7e12, pinhole 7.4e12, stepped 7.9e12 (photons for 1 pm on the surface: 7.3e13 / 8.5e13 / 5.8e13; oapnoise22) |
| Zernike dimple, stepped (4) | 0.989 / 5 pm | 0.999 / 4 pm | 0.984 / 0.68 nm | 7.5e12 |
| vector Zernike, polarized dimple (2, at once) | 0.9935 / 4 pm | 0.9992 / 3 pm | 0.9999 / 0.33 nm | 5.3e12 |
| stepped pinhole with a shutter frame (6) | 0.9935 / 4 pm | 0.9992 / 3 pm | 0.9985 / 338 pm | 7.0e12 |
| P/SRI, reference arm traced (4) | 0.9935 / 2 pm | 0.9924 / 1 pm | 0.9926 / 144 pm | 1.5e13 |

Run tags lens_deck, oapifo2, oapifol2, oaploop22, matbase, v193base,
pdi193fbase, pfdeck, vloop193, ploop193, loop193.  The four best readings are within 1%
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

The three efforts run on one engine and one set of tools, and that is the
point of reporting them together.  The CTB and e2e6m share the Bench
builder, the DM as an influence-function surface, the sphere-bracketed
propagation recipe validated against PROPER, the mask library, the EFC
driver and the contrast scorer; the gauge deck shares the same DM doctrine
and its readings are the coronagraph's input-field sensor in the field
servo now being modeled.  Each lane is driven by one parameter sheet and
one runner, every number carries its run tag, and every deck figure is the
tool's own output.  Three things the last week established across the
lanes: a document can be read wrongly and only a trace settles it (the OAP
poles, the CTB beam, the tail tuner); calibrate on the working surface,
never at null; and a servo's floors are analytic, so the design levers can
be read off without another run.

Next, in order: the interferometer's rows with real plate thicknesses and
the wrap-arithmetic verdict (running); the detector-leg tuner's fix; the
coronagraph field servo in three steps on the CTB deck; the CTB regenerated
at the sheet's intended beam once the gauge queue is drained; the gauge
deck's sign-off and the migration of its tools into the templates tree; and
the CTB roadmap's aberration and drift arcs, which the field servo is built
to answer.

## Sources
deck_ctb.md / CTB_PROP_STATUS.md / README (bench_ctb); deck_e2e6m_r2.md / e2e6m_r2_LOG.md / README (e2e6m, e2e6m_r2); deck_gauges.md / BRIEF_gauge_deck.md / the three lane reports (REPORT_gauge_ifo, REPORT_reflective, REPORT_gauge_pdi, zwfs_dm96 README); NOTES_gauge_in_coronagraph.md; memory project_ctb_diffraction / project_e2e6m / project_tg96_gauge.

## Notes on sources (sections 1-2)
1. Outline 1.2 lists the vortex-matched static at 2.9e-7, Roddier at 3.2e-6 and dual-zone at 6.4e-6.  Deck v5 slide 5 has vortex with matched Lyot 1.4e-8 (the pixel-averaged core; 2.9e-7 was the direct-sampled artifact), Roddier π-mask 2.4e-6, dual-zone 6.8e-6.  The deck values are used.
2. Outline 2.2 gives the CTB DM pitch as 1.34 mm on a 42.75 mm beam.  The traced beam is 21.2 mm and the pitch 0.67 mm (CTB slide 9; README "Source model"); 42.75 mm is the generator's intended beam, never delivered.  The measured values are used.
3. Outline 2.2 states the e2e6m pitch as 1.48 mm.  `r1_dm.m` gives 2 x 23.771 / 32 = 1.486 mm, which rounds to 1.49 mm; 1.48 is the truncated value.  The text uses 1.49 mm.
4. Outline 2.2 asserts "the e2e6m DMs are sized from their beam; the CTB's were not."  Both lattices are sized from traced footprints (`r1_dm.m` line 57; `ctb_dm.m` default 21.3 mm per REPORT_field_servo.md section 0).  The CTB slip is in the generator's source cone (half the intended beam), not in the DM lattice.  The text says so.
5. Outline 2.1 quotes the closure as "0.35%"; the deck prints "worst relative error 0.0035" (e2e6m slide 14).  Same number, stated as a percentage in the text.
6. The task wording gives the e2e6m beam as "47.5 mm" and e2e6m slide 7 says "47 mm"; `r1_dm_report.txt` prints 47.5 mm (2 x 23.771 = 47.54).  The table quotes 47.5 mm with the report as source.
7. The vortex mono closed-loop floor is 5.8e-13 on CTB slide 11 (the scalar/polarization run) and 8.3e-13 on slide 12 (the bandwidth ladder's mono point); both are roundoff-class and Table 1 carries both with their slides.
