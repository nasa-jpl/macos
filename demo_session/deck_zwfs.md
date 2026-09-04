<!--
deck_zwfs.md — the Zernike wavefront sensor for the DM gauge, and how
it will be compared to the Twyman–Green IFO.  DRAFT — pending Dave's
sign-off; the builder suppresses the export mark on DRAFT decks.
Build: python3 make_brief_slides.py deck_zwfs.md
Status: covers the campaign through S1 (mask + response gates green);
battery, differential head-to-head and the vector step are pending and
say so.  Sources: templates/40_benches/zwfs_dm96 (S1, d6791d1),
vsg_wip/vsg2_params.m §9 (hardware), tg_psi_dm96 run 10 (IFO numbers),
BRIEF_zwfs_campaign.md (plan + rulings 2026-09-04).
-->

# A Zernike sensor for the DM gauge
The same 96 mm bench with the reference arm removed: a quarter-wave dimple at focus turns one camera frame into a wavefront measurement — set against the phase-shifting interferometer on the same deformable mirror
D. C. Redding, with Claude Code.
September 2026.
DRAFT — pending review.  Campaign in progress: sensor model verified; battery and head-to-head pending.

## The idea: the beam interferes with its own core | An etched dimple at the focus phase-shifts the heart of the beam; that light spreads back over the pupil as a built-in reference
::: left
- **One transparent plate, one tiny etched spot.**  At the focus, the core of the beam passes through the dimple and picks up a quarter-wave phase shift; everything outside misses it.  The shifted core light spreads back across the whole pupil image and interferes with the rest — a common-path interferometer with no second arm.
- **Pupil brightness becomes a phase map.**  Near a quarter wave of shift, the camera intensity is close to linear in the wavefront: one frame per measurement, no moving parts, no polarization optics.
- **The two prices, both already measured or scheduled:** the reference is made from the beam itself, so the slowest shapes partially cancel against themselves (defocus reads at 0.32 of its true size — measured); and the linear reading holds only for small departures — the range limit is a scheduled measurement.
::: right
![The focal spot with the dimple footprint on its core (left); the mask itself — a 1.57 rad phase disk, edges gray from area-weighted supersampling (right).](figs/zwfs_mask_fig.png){h=2.5}
~ The dimple takes 28.0% of the light — the encircled energy a 1.06 λF/D disk should take from a focused spot.

## The setup: the interferometer's test arm, alone | Same source, same splitter, same 700 mm DM leg, same lenses — the reference arm and every polarization element simply removed
::: left
- **Reused whole:** the TG96 bench solved for the 96×96, 1 mm-pitch DM — splitter at 7°, clearance-solved legs, the retuned detector optics.  The dimple sits at the internal focus of the existing detector leg, in the mask seat the bench already carried; the camera still views the DM image.
- **Removed:** the reference flat and its leg; the polarizer, four arm waveplates, output waveplate and analyzer.  Nothing else moves.
- **Why that discipline matters:** when the two instruments are scored against each other, the only difference is the sensing principle — not the glass, not the geometry.
::: right
![The traced train: source and splitter at left, the 96 mm DM on its 700 mm leg, and the detector tail with the dimple at its internal focus.](figs/zwfs_layout.png){h=2.7}
~ Two model findings on the way, recorded in the campaign README: the mask must sit at the real focus (−5.58 mm from the textbook seed — the lens is the optimized one), and representing a µm-scale mask at focus needs a diffraction bracket around that plane (reference-sphere legs, the coronagraph-bench idiom).

## The mask is real hardware, modeled at its own numbers | 346.2 nm of etch in fused silica = 1.571 rad at 632.8 nm; a 9-spot substrate, one spot in the beam at a time
::: full
- **The part:** a transmissive fused-silica substrate carrying a 3×3 array of etched dimples of graded diameter — the VSG2 wavefront-sensor mask.  Etch depth 346.2 nm → phase 2π(n−1)d/λ = 1.571 rad ≈ a quarter wave at 632.8 nm.  Translating the substrate selects a spot; only one is ever in the beam.
| spot diameter (λF/D) | 1.06 (default) | 1.22 | 2.0 | 3.0 |
| here, at F/4.2 | 2.79 µm | 3.22 µm | 5.3 µm | 7.9 µm |
- **In the model:** the etched disk is a complex mask multiplied into the propagated field at the focal plane — edges area-weighted at 8× supersampling, the disk centered on the measured spot (the same alignment the real substrate gets by translation).  An exact identity gates it: the masked field equals the sum of the direct and dimple-diffracted parts to 15 digits, so the mask provably acts on the light.
- **Spot size is also a sampling lever:** a larger spot relaxes the focal-plane grid a high-resolution pupil leaves behind — the 96×96 runs use the 3.0 λF/D spot; the trade across the real spot table is a scheduled study.
~ Hardware numbers from the VSG2 ZWFS upgrade deck, carried in vsg2_params.m; substrate class Thorlabs W4101FT1.

## First light in the model | A known 8 nm ripple commanded on the DM comes back from ONE camera frame at gain 0.93 — and reads as nothing when the dimple is removed
::: left
- **The response chain, gated:** dimple resolved at 6.3 px on the focal grid; masked frame = direct + diffracted parts at 3×10⁻¹⁶; core fraction 0.280.
- **The measurement:** an 8 nm radial ripple on the DM, reconstructed from a single frame against the flat-state reference maps — gain 0.932, residual 0.45 nm, with the camera-to-DM mapping taken from traced rays (magnification 10.16, anamorphism 0.00%).
- **The control:** the same DM state read without the dimple recovers gain −0.07 — plain pupil imaging is phase-blind; the signal is the dimple's.
- **The measured instrument property:** defocus reads at 0.32 — the self-reference cancellation.  The sensor's response will roll off at low spatial frequency, the mirror image of the interferometer.
::: right
![Commanded figure (left) against the single-frame recovery (right): same rings, same scale.](figs/zwfs_response.png){h=2.6}
~ All numbers from zwfs_s1_report.txt (five gates, 20 s a run, development resolution); full-resolution battery pending.

## Side by side with the Twyman–Green | Mirror-image instruments: the interferometer dims at the finest patterns, the Zernike sensor at the slowest — and it pays for its one-frame economy with range
::: full
| | Twyman–Green PSI (measured) | Zernike sensor (this campaign) |
| reference beam | the second arm's flat | made from the beam's own focal core |
| one measurement costs | 6 traces (3 per arm), 4 fringe frames | 1 camera frame (+1 stored flat set) |
| polarization hardware | polarizer, 5 waveplates, analyzer | none |
| moving parts | none (polarization-stepped) | none |
| null, nothing aligned | 0.134 nm | to be measured |
| response rolls off at | fine patterns: 0.50 at the 96×96 checkerboard | slow patterns: defocus 0.32 (measured) |
| range in one frame set | half a wave, and unwrappable | small departures only; the fold point is a scheduled measurement |
| a 10 nm actuator change, differenced | 0.021 nm | the head-to-head |
- **Reading:** the interferometer measures any surface and pays at high spatial frequency; the Zernike sensor measures small changes cheaply — one frame, no polarization train — and pays at low spatial frequency and in range.  Neither table column is a winner until both are scored on the same question.
~ PSI column: tg_psi_dm96 run 10.  Differential entries are map-space; both instruments will be restated in actuator space before the head-to-head (next slide).

## The head-to-head, and how it will be scored | Same DM truth, same battery — and the score is how well an actuator change is recovered, through the DM's own actuation model
::: full
- **The metric (the design ruling):** measure a working state, apply a known deviation, measure again, difference the measurements — then fit the DM actuation model to the result and score the recovered actuator changes, not map pixels.  The actuator fit absorbs each instrument's roll-off inside the DM's own band, so the comparison is fair to both.
- **The interferometer's differential floor stands at** 0.021 nm on a 10 nm single-actuator change and 3.67 nm on a 10 nm random pattern (map space, base-independent) — both to be restated in actuator space before the sensor is scored against them.
- **The campaign from here:** camera registration by the two-poke doctrine at small strokes → the same 12-pattern battery at 48×48 and 96×96 → the differential head-to-head, including the working-state size where the sensor's linear reading folds.
- **Then the vector step:** a polarizing metasurface in place of the etched dimple gives the two polarizations opposite phase shifts — two pupil images in one frame, shifted opposite ways, so the sign ambiguity disappears and the usable range grows.  The engine's vector-diffraction and per-polarization mask machinery is already in place from the coronagraph work.
~ Registration pokes drop to ~20 nm for the sensor: a 150 nm poke is 3 rad of phase — outside the linear reading.

# Backup

## Provenance and records | Every number re-derives from a committed script
::: full
- Campaign: templates/40_benches/zwfs_dm96 — zwfs_s1.m (five gates: sampling, superposition, reference wave, response, no-dimple control), zwfs_mask.m, zwfs_s1_figs.m (the figures here), zwfs_s1_report.txt, README with the findings chain; plan and rulings in macos/BRIEF_zwfs_campaign.md.
- Findings earned in S1, recorded so they are not re-derived: the stock detector leg is ray-traced (geometric) at the mask plane — a λ/D-scale mask needs the reference-sphere diffraction bracket (now a builder option, default off, existing decks bit-identical, suites green); the mask must sit at the optimized lens's real focus (−5.58 mm from the thin-lens seed, found by a peak scan the gate re-runs); the camera-to-DM frame must come from traced rays — a support-area estimate was 25% off and read a fine test pattern as "no response."
- Sign convention pinned by gate: surface height = −φ·λ/4π on this deck (single reflection doubles height).
- Hardware source: "VSG2 Zernike Wavefront Sensor Update -v2" deck, parameters carried in templates/40_benches/vsg_wip/vsg2_params.m §9.
- Interferometer comparison column: templates/90_polarization/tg_psi_dm96, run 10 (tg96_report.txt).
