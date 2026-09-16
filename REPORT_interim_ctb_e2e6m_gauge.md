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

- 3.1 The three jobs (measure to picometers, capture 100-200 nm of wavefront, hold) and the four approaches on one bench (Twyman-Green interferometer, stepped Zernike dimple, vector Zernike, stepped pinhole / P-SRI); the bench of record at 22.5 deg; the reflective front end redesigned (no fold coma: a 25 mm conjugate error; seat diffraction-limited at zero trim; rows 0.99 / 2.4 pm).
- 3.2 The record on the 30 nm working surface with the matrix measured on it (side-by-side table); capture range (self-referenced 36-70 nm; external / shutter 480 nm+); photons per measurement; the servo (3 pm from 1.5e12 photons per cycle); the descent.
- 3.3 Systematics priced and calibrated (metasurface, arm polarization, analyzer, complex amplitude V5); what is not yet modeled (substrates, camera, the photonic phase shifter -- "magic" steps).
- 3.4 Recommendation as it stands (hold with the vector Zernike, capture with the interferometer or the pinhole's shutter frame, lenses on the record with the reflective rig re-measured) -- with the reflective verdict pending TO's vector rows.

## 4. The gauge inside the coronagraph (NOTES_gauge_in_coronagraph.md)

- 4.1 Packaging: face-on blocked in both packages (OAP1 + the other DM on the normal); out-of-plane 15 deg clears the CTB by 31 mm; nothing to 30 deg clears the space relay (200-250 mm DM-OAP legs).  The clearance tool.
- 4.2 The options A-G, one line each; the pupil-dichroic field servo (C+) as the direction: holds the coronagraph's input field, amplitude and phase, against everything upstream of the pickoff.
- 4.3 The bounds: what the pickoff cannot see; DM2's Fresnel amplitude authority (stroke x530 at 2 cycles, x8 at Nyquist); the out-of-band transfer; stellar photons (10 pm per 200 s on V=5); the two DMs are NOT separable in-band from one conjugate in either package (12% / 7% conversion at Nyquist) -> a second conjugate or differential attribution.
- 4.4 What is queued (brief item 7): the reading at the apodizer conjugate, the two-DM matrix, the servo under an upstream drift with contrast scored through the EFC chain.

## 5. Where the three meet, and what is next

One engine, one Bench builder, one DM doctrine, one propagation recipe, one contrast scorer; the gauge's readings become the coronagraph's input-field sensor; e2e6m carries the telescope side.  Next: the CTB re-score and its aberration arc; the field-servo model; the gauge deck's sign-off; migration of the gauge work into templates.

## Sources
deck_ctb.md / CTB_PROP_STATUS.md / README (bench_ctb); deck_e2e6m_r2.md / e2e6m_r2_LOG.md / README (e2e6m, e2e6m_r2); deck_gauges.md / BRIEF_gauge_deck.md / the three lane reports (REPORT_gauge_ifo, REPORT_reflective, REPORT_gauge_pdi, zwfs_dm96 README); NOTES_gauge_in_coronagraph.md; memory project_ctb_diffraction / project_e2e6m / project_tg96_gauge.
