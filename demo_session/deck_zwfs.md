<!--
deck_zwfs.md — the Zernike wavefront sensor for the DM gauge, and how
it compares to the Twyman–Green IFO.  DRAFT — pending Dave's sign-off;
the builder suppresses the export mark on DRAFT decks.
Build: python3 make_brief_slides.py deck_zwfs.md
2026-09-10 fold 3 (S7 + S8): the sensor-model correction (the S1–S6
sensor was Fresnel-defocused), the five readings incl. the iterated-
reference exact reading with the refined base prior, the sampling trade
closed at MODEL 2048, the colour and photon re-runs on the corrected
model, the multi-site working-state range, and the parameterized runner.
Legacy-model slides (S1 poke/defocus, the S6 colour claim) moved to
Backup with their attribution.  Figures are the tools' own PNGs
(zwfs_s7iter.png; zwfs_run outputs in zwfs_dm96/runs/<tag>/), placed
unmodified except autocrop of one panel (crop_zwfs_rec193full_color_L).
2026-09-09 fold 2: the literature scan.  Sources: templates/40_benches/
zwfs_dm96 (README S7/S8 bullets, runs/), macos/BRIEF_zwfs_campaign.md,
tg_psi_dm96 run 10 + S4/S5 (IFO numbers).
-->

# A Zernike sensor for the DM gauge
The same 96 mm bench with the reference arm removed: a quarter-wave dimple at focus turns one camera frame into a wavefront measurement — set against the phase-shifting interferometer on the same deformable mirror
D. C. Redding, with Claude Code.
September 2026.
DRAFT — pending review.  Campaign status: sensor model corrected (the early sensor was defocused); five readings measured on the corrected model, sampling closed at a 2048 grid, colour and photon pricing re-run; one parameterized runner reproduces every number; vector step pending.

## The idea: the beam interferes with its own core | An etched dimple at the focus phase-shifts the heart of the beam; that light spreads back over the pupil as a built-in reference
::: left
- **One transparent plate, one tiny etched spot.**  At the focus, the core of the beam passes through the dimple and picks up a quarter-wave phase shift; everything outside misses it.  The shifted core light spreads back across the whole pupil image and interferes with the rest — a common-path interferometer with no second arm.
- **Pupil brightness becomes a phase map.**  Near a quarter wave of shift, the camera intensity is close to linear in the wavefront: one frame per measurement, no moving parts, no polarization optics.  The exact per-pixel inversion (slide 6) removes the small-phase limit.
- **The two prices:** the sensor cannot see piston and dims toward the very lowest spatial frequencies — the reference is made from the beam itself.  On the corrected model the dip sits at 0.5–1 cycle per aperture (gain 0.5 at 1 cycle); defocus and everything finer read near 1.
::: right
![The focal spot with the dimple footprint on its core (left); the mask itself — a 1.57 rad phase disk, edges gray from area-weighted supersampling (right).](figs/zwfs_mask_fig.png){h=2.5}
~ The dimple takes 28.0% of the light — the encircled energy a 1.06 λF/D disk should take from a focused spot.

## The setup: the interferometer's test arm, alone | Same source, same splitter, same 700 mm DM leg, same lenses — the reference arm and every polarization element simply removed
::: left
- **Reused whole:** the TG96 bench solved for the 96×96, 1 mm-pitch DM — splitter at 7°, clearance-solved legs, the retuned detector optics.  The dimple sits at the internal focus of the existing detector leg, in the mask seat the bench already carried; the camera still views the DM image.
- **Removed:** the reference flat and its leg; the polarizer, four arm waveplates, output waveplate and analyzer.  Nothing else moves.
- **Why that discipline matters:** when the two instruments are scored against each other, the only difference is the sensing principle — not the glass, not the geometry.
::: right
![The train as traced: source and splitter at left, the 96 mm DM on its 700 mm leg, the detector tail converging to the dimple's focus, and the camera at the pupil image.](figs/zwfs_render_rig.png){h=2.3}
~ Representing a µm-scale mask at focus needs a diffraction bracket around that plane: two reference spheres, one before and one after the mask.  The two must carry the same radius — the finding on slide 5.

## The mask is real hardware, modeled at its own numbers | 346.2 nm of etch in fused silica = 1.571 rad at 632.8 nm; a 9-spot substrate, one spot in the beam at a time
::: full
- **The part:** a transmissive fused-silica substrate carrying a 3×3 array of etched dimples of graded diameter — the VSG2 wavefront-sensor mask.  Etch depth 346.2 nm → phase 2π(n−1)d/λ = 1.571 rad ≈ a quarter wave at 632.8 nm.  Translating the substrate selects a spot; only one is ever in the beam.
| spot diameter (λF/D) | 1.06 (default) | 1.22 | 2.0 (used) | 3.0 |
| here, at F/4.2 | 2.79 µm | 3.22 µm | 5.3 µm | 7.9 µm |
- **In the model:** the etched disk is a complex mask multiplied into the propagated field at the focal plane — edges area-weighted at 8× supersampling, the disk centered on the measured spot (the same alignment the real substrate gets by translation).  An exact identity checks it: the masked field equals the sum of the direct and dimple-diffracted parts to 15 digits, so the mask provably acts on the light.
- **Spot size selects the low band:** measured on the corrected model, the 3.0 λF/D spot deepens the dimple's own dip below 2 cycles per aperture (gain 0.22 against 0.50 at 2.0 λF/D) and changes nothing above it; 2.0 λF/D is the spot of record.
~ Hardware numbers from the VSG2 ZWFS upgrade deck, carried in vsg2_params.m; substrate class Thorlabs W4101FT1.

## First light in the model | A known 8 nm ripple commanded on the DM comes back from one camera frame at gain 0.93 — and reads as nothing when the dimple is removed
::: left
- **The response chain, checked:** dimple resolved at 6.3 px on the focal grid; masked frame = direct + diffracted parts at 3×10⁻¹⁶; core fraction 0.280.
- **The measurement:** an 8 nm radial ripple on the DM, reconstructed from a single frame against the flat-state reference maps — gain 0.932, residual 0.45 nm, with the camera-to-DM mapping taken from traced rays (magnification 10.16, anamorphism 0.00%).
- **The control:** the same DM state read without the dimple recovers gain −0.07 — plain pupil imaging is phase-blind; the signal is the dimple's.
- **Low orders are less trouble than presumed:** defocus reads at 0.986 at full resolution.  The self-reference blindness lives at piston and tilt.
::: right
![Commanded figure (left) against the single-frame recovery (right): same rings, same scale.](figs/zwfs_response.png){h=2.6}
~ All numbers from zwfs_s1_report.txt (five checks, 20 s a run, development resolution) — taken on the sensor model later found defocused (next slide); the five checks are structural and stand.

## The sensor model was defocused; corrected, the ringing goes away | The mask bracket's exit sphere carried the wrong radius — a 4.9 m Fresnel defocus of the pupil image under every S1–S6 number
::: left
- **The defect:** the two reference spheres bracketing the mask must share one radius so the unmasked round trip is the identity.  The emitted exit sphere sat at 23.9 mm against the entrance sphere's 352.7 mm; the engine's sphere-to-plane step then applied a focal quadratic phase — a Fresnel defocus of the pupil image by 4.86 m.  Found building the exact reconstructor, which needs the pupil amplitude to be state-independent.
- **The symptoms it explains:** round trip 0.159 on a flat DM (now 1.8×10⁻¹⁵); a 30 nm DM state modulated the pupil amplitude by 29% rms (now 4×10⁻¹⁶); the single-actuator kernel rang (peak 0.27, ring −0.16; now 0.75, ring-free); a transfer null near 30 cycles per aperture — a Talbot null, predicted at 34.
- **What the correction alone buys** (linear reading, 96×96, same frames): single-actuator floor on the 30 nm working surface 744 → 67 pm (SNR 9 → 80); the grid case the sensor could not detect goes from SNR 1.5 to 14; dense random 42 → 9.6 nm.
::: right
![Left: the corrected model's modal transfer for the three reading classes against the legacy linear record (dashed) with its null; right: the single-actuator differential gain as the working surface grows.](figs/zwfs_s7iter.png){h=3.0}
~ The legacy emission is kept as a builder option so the S1–S6 record reproduces byte for byte; the interferometer twin traces the same tail geometrically and was never affected.  Records: zwfs_s7iter_report.txt, README S7.

## Five readings from the same camera frames | The exact inversion with an iterated reference wave and a one-time branch prior reads a 10 nm change to 26 pm from one frame
::: full
- **The readings:** linear (L, one frame, small-phase); exact with a frozen reference wave (F, one frame); exact with the reference wave re-propagated from the estimate through the mask model, five passes (I, one frame); I with the working surface's branch resolved once by a phase-stepped retrieval of the base, refined with the iterated reference (I+, one frame per measurement after four frames once); phase-stepped through three etch depths (S, four frames).  The reference-wave surrogate matches the engine's own field to 2×10⁻¹⁵.
| 10 nm change on the 30 nm working surface, 96×96 | L | F | I | I+ | S |
| gain (raw / corrected) | 0.41 / 0.40 | 0.91 / 0.94 | 1.01 / 1.04 | 0.91 / 0.95 | 0.75 / 0.89 |
| floor over the unpoked actuators, pm (raw / corrected) | 60 / 57 | 46 / 52 | 34 / 46 | 26 / 40 | 23 / 42 |
| detection SNR (raw) | 68 | 196 | 298 | 355 | 323 |
| 47 actuators changed 1 nm on the same surface: SNR (raw) | 13.6 | 12.4 | 12.6 | 37.3 | 38.3 |
- **The branch is the one physical limit:** with the true branch and the engine's own reference wave the solve is exact to 3×10⁻¹⁴; on the 30 nm surface 7.8% of pixels sit beyond the quarter-wave sensor's fold, and the refined prior finds them on 99.99% of pixels where the plain stepped prior misses 3% — enough to flip the sign of a single actuator.
~ Compliant configuration (385 rays across, 2048 grid; slide 7); actuator-space scoring through the DM's influence model; "corrected" = after the measured modal transfer; interferometer on the same row: 0.92 / 46 pm.  Records: runs/m2048/m2048_report.txt.

## Sampling: the two interfaces pull opposite ways | Dimple pixels at the mask go as grid/rays, camera pixels per actuator as rays — only a larger grid buys both, and a 2048 grid fits the box
::: left
| rays across / grid | dimple, px | px per actuator | hold-out gain | I+ on the 30 nm surface |
| 193 / 1024 (record) | 7.92 | 2.5 | 0.900 | 0.79, 13 pm, SNR 584 |
| 385 / 1024 | 3.96 | 5.0 | 0.935 | 0.91, 25 pm, SNR 359 |
| 385 / 2048 (compliant) | 7.92 | 5.0 | 0.935 | 0.91, 26 pm, SNR 355 |
- **The rule:** mask-plane pixels per λ/D = 0.74 × grid / rays; camera pixels per actuator ∝ rays.  Halving the ray count doubles focal resolution and halves pupil sampling.
- **What it settled:** the hold-out gain is set by the ray count and not by the dimple's sampling — identical at the 4-px and 8-px dimple, identical at spot 3.0.  The remaining 6.5% is neither; the suspect is the kernel's variation between the centre, where it is measured, and the held-out site.
- **The 2048 grid runs on a 30 GB machine** through a trimmed engine size table dropped into the run directory (grid-surface slots 200 → 4 saves 1.7 GB): 32 min, under 4 GB.
::: right
![The compliant configuration (385 rays, 2048 grid): modal transfer per reading class (left) and the working-surface ladder per reading (right).](figs/zwfs_m2048_battery_n96.png){h=2.9}
~ Hold-out = a single 20 nm actuator away from the calibration site, raw gain through the measured kernel.  I+ row values raw.  Records: runs/rec193, runs/ng385, runs/ng385s3, runs/m2048.

## Working-surface range, measured on 47 sites at once | The one-frame reading with the base prior holds to 40 nm rms on 1 mm actuators and 50 nm on 2 mm; the four-frame reading to 50–60 nm
::: full
| 10 nm changes on 47 grid sites; gain (SNR) vs working-surface rms | 30 nm | 40 nm | 50 nm | 60 nm |
| 96×96, 1 mm pitch — I+ (one frame) | 0.81 (26) | 0.85 (30) | 0.76 (9) | 0.31 (1) |
| 96×96 — S (four frames) | 0.92 (29) | 0.84 (27) | 0.83 (17) | 0.66 (10) |
| 48×48, 2 mm pitch — I+ | 0.89 (56) | 0.80 (41) | 0.84 (40) | 0.83 (3) |
| 48×48 — S | 0.89 (50) | 0.83 (44) | 0.99 (16) | 0.55 (11) |
- **Why the earlier "holds to 60 nm" moved:** it was one actuator at one site.  On 47 sites the reading fails where a subset of sites falls beyond the fold, and the floor jumps from 0.3 to 3 nm — the collapse is visible as the floor, not as the gain.  Past 120 nm rms every reading aliases (the stepped retrieval itself wraps).
- **Doctrine, updated:** the one-frame reading with the refined base prior owns working surfaces to 40 nm rms (1 mm actuators) or 50 nm (2 mm); the four-frame stepped reading extends to 50–60 nm and owns dense commands; the interferometer owns everything beyond (it does not fold to 480 nm).
~ Floor here = rms of the unpoked actuators under 47 simultaneous 10 nm pokes (their crosstalk), a different quantity from the single-site floor on slide 6.  Compliant configuration; record: runs/m2048.

## Colour is no longer a lever | With the transfer null gone, five colours lift the combined transfer 3% instead of 3×, and the combination is neutral or worse on the rows
::: left
![Modal transfer of the linear reading at five wavelengths on the corrected model: no null, only the dimple's low-order dip that moves with the spot's size in λ/D; the five-colour combination (orange) sits at 0.99–1.](figs/crop_zwfs_rec193full_color_L.png){h=3.4}
::: right
- **What moved:** the "blind band near 30 cycles per aperture that migrates with wavelength" was the defocused sensor's Talbot null.  On the corrected model the five single-colour transfers bottom at 0.92–0.97 after the Wiener step and the combination at 0.991.
- **On the rows** (30 nm working surface, single 10 nm change): I+ at 632.8 nm SNR 269, five-colour combination 183, best pair 700+780 nm 255; the 47-site case 32 → 32.  The linear reading gets worse under combination (84 → 34): at 480 nm the dimple is 2.1 rad and the small-phase reading turns negative on this surface, and an equal-weight combiner inherits it.
- **What stays chromatic is range:** 780 nm is the best single colour on nearly every row — a smaller phase per nanometre of height and a 1.6 λ/D dimple.  The interferometer is unchanged at every colour (its deck).
~ 632.8 / 480 / 532 / 700 / 780 nm through one physical mask, every calibration redone per colour; record configuration (193 rays, 1024 grid); runs/rec193full.

## Photons: 10¹⁴ per DM state reaches 1 pm | Every reading prices within a factor of three, 25× cheaper than the defocused sensor's pricing; the branch prior's own noise costs 5%
::: left
![Noise on the single-actuator estimate against photons per DM state, per reading; the dotted line is the 1 pm target.](figs/zwfs_rec193full_noise.png){h=3.4}
::: right
| reading | photons per state for 1 pm |
| linear (L) | 5.4×10¹³ |
| exact, frozen reference (F) | 3.7×10¹³ |
| exact, iterated reference (I) | 6.4×10¹³ |
| I + refined base prior (I+) | 8.8×10¹³ (8.4×10¹³ with a noiseless prior) |
| phase-stepped (S) | 1.0×10¹⁴ |
| interferometer, four-step | 8×10¹⁴ |
- **Method:** noiseless frames captured once per state, shot noise added numerically, eight realizations per point; the photon budget per state is shared by a reading's frames (one, or four).  Each reading's high-photon floor converges to its systematic floor on slide 6.
- **Reading:** 10¹⁴ photons at 633 nm is 30 µJ — the 1 pm budget is systematics and gain stability, not light.
~ Scenario: the head-to-head row (single 10 nm change on the 30 nm working surface, 96×96, record configuration).  Read noise, drift and calibration noise are the next terms once a use case fixes them.  Records: runs/rec193full; interferometer tg96_s5noise_report.txt.

## Run it yourself | One parameter sheet and one entry point reproduce every number here; the defaults reproduce the S7 record to the last printed digit
::: full
- **Two files:** zwfs_params.m returns every knob at its value of record — engine grid and ray count, wavelength, the bench's lengths and tuned figures, the mask's etch and spot, the sampling budget, registration, the DM configurations, the battery's amplitudes and seeds, the colour and noise stages; zwfs_run.m runs the selected stages and writes the report, a .mat and the figures into a run folder named by its tag.
| call | what it does |
| zwfs_run | the record: bench checks, battery, figures |
| zwfs_run('tag','ng385', 'NGRID',385) | the 1 Mpix-class camera |
| zwfs_run('MODEL',2048, 'NGRID',385, 'param_file','macos_param_2048.txt') | the compliant configuration on a 30 GB machine |
| zwfs_run('stages',{'battery','color','noise','figs'}) | everything on this deck's later slides |
| ./zwfs_batch.sh TAG "the same arguments" | headless, memory-capped, logged |
- **Checks that run every time:** the round trip through the mask bracket (must be the identity), the reference-wave surrogate against the engine's field, the pupil amplitude's independence of the DM state, the sampling budget (dimple pixels; pixels per actuator), and the parity/sign search of the camera-to-DM registration.
~ Equivalence check: 64 row and ladder lines against the S7 report, 8 differ in the last printed digit (the actuator fit's solver tolerance).  README "Run it yourself" carries the parameter table.

## Side by side with the Twyman–Green | The interferometer dims at fine patterns and never folds; the sensor reads a change from one frame at the same floor and folds past 40–60 nm
::: full
| | Twyman–Green PSI (measured) | Zernike sensor (measured) |
| reference beam | the second arm's flat | made from the beam's own focal core |
| one measurement costs | 6 traces (3 per arm), 4 fringe frames | 1 camera frame (I+; base 4 frames once) |
| polarization hardware | polarizer, 5 waveplates, analyzer | none |
| moving parts | none (polarization-stepped) | none |
| null, nothing aligned | 0.134 nm | 0 by construction (the flat is the reference) |
| response rolls off at | fine patterns: 0.50 at the 96×96 checkerboard | 0.5–1 cycle per aperture (0.5 at 1); piston unseen |
| working-surface range | to 480 nm rms and beyond | 40 nm rms one-frame (1 mm actuators), 50–60 nm four-frame |
| a 10 nm actuator change on a 30 nm surface | 0.92, floor 46 pm | 0.91, floor 26 pm (I+, one frame) |
| photons per state for 1 pm | 8×10¹⁴ | 0.4–1×10¹⁴ |
- **Reading:** the interferometer measures any surface and pays at high spatial frequency; the Zernike sensor measures small changes on a small working surface more cheaply — one frame, no polarization train, a floor on par or better — and pays at the lowest spatial frequencies and in range.
~ PSI column: tg_psi_dm96 run 10 and S4/S5; sensor column: the compliant configuration (slide 7) except photons (record configuration).  Both in actuator space through the same estimator code (dm_gauge_lib).

## What the Zernike-sensor literature offers, in priority order | Twelve papers read; the first import — an iterated reference wave — is now built and measured (slide 6)
::: full
| # | paper | what it gives this model |
| 1 | Ruane, Wallace, Steeves et al. 2020, JATIS 6, 045005 | Exact per-pixel reconstruction from the pupil amplitude, the reference wave b and the masked image; a sensitivity factor per pixel; a four-term systematic budget; 1 pm reached at 4.4×10⁵ frames; no loss to 25 % bandwidth |
| 2 | Doelman, Fagginger Auer, Escuti, Snik 2019, Opt. Lett. 44, 17 | Vector Zernike: ±π/2 on opposite circular polarizations gives two pupil images, an exact phase-and-amplitude solution, an iterated b, achromatic to 100 % bandwidth; the leakage terms are written down |
| 3 | N'Diaye, Vigan, Dohlen et al. 2016, A&A 592, A79 (ZELDA II) | Second-order reconstruction with b from the mask model; 1.06 λ/D at π/2; range −0.14 to +0.36 λ; the sensitivity factor drifts 10 % day to day on SPHERE |
| 4 | Chambouleyron, Cissé, Salama, Haffert, Wallace et al. 2024, SPIE 13097 | A phase-shifted pair (±π/2) yields cos φ and sin φ; iterative reconstructors reach about 1 rad rms; the bench limits were polarization crosstalk and defocus between the two pupils |
| 5 | Haffert 2024, A&A 683, A113 | Iterative nonlinear reconstructors reach machine precision below 0.25 rad rms; phase sorting extends to 0.75 rad; adding wavelengths to 1.4 rad rms |
| 6 | Darcis, Haffert, Chambouleyron, Doelman et al. 2025, A&A | Multi-wavelength gradient descent through the forward model: dynamic range for the scalar sensor, photon robustness across a wider band, two-wavelength unwrapping of petal errors |
~ Priority = payoff against the campaign's measured weak results.  Synopses and links: macos/REPORT_zwfs_lit_scan.md.

## Demonstrations, segmented mirrors, and the in-house precedent | Picometers by alternating DM states and averaging; segment piston by model iteration, underestimated below 50 % Strehl
::: full
| # | paper | what it gives this model |
| 7 | Steeves, Wallace, Kettenbeil, Jewell 2020, Optica 7, 1267 | 1.6 pm repeatability in 4.3 s by alternating flat and waffle DM states and averaging the differences; a reconstruction robust at the highest spatial frequencies |
| 8 | Wallace, Rao, Jensen-Clem, Serabyn 2011, SPIE 8126 | All-reflective phase-shifting Zernike interferometer: a dynamic, arbitrary core phase shift read in four steps gives phase and amplitude; low sensitivity to vibration, polarization and wavelength |
| 9 | Moore and Redding 2018, SPIE 10698 | Nonlinear polychromatic physical-optics reconstruction for picometer differential metrology on LUVOIR; the in-house precedent for the reconstructor above |
| 10 | Keck vector-Zernike segment control 2024 (arXiv 2404.08728); Wallace et al. 2022 (arXiv 2205.02241) | Segment piston by model-based iteration; 11 nm rms piston uncertainty; underestimation 2 to 4× below 50 % Strehl; fabricated shifts 0.30π and 0.68π instead of ±0.5π |
| 11 | HiCAT mid-order Zernike sensor 2024 (arXiv 2409.03411) | Per-segment piston, tip and tilt by interaction matrix; a minimal step of 125 ± 31 pm at SNR 4; 14-bit DM quantization reads as steps in the response |
| 12 | Shi et al. 2015, SPIE 9605 (Roman low-order sensor) | A reflective dimple on the focal-plane mask senses Z2 to Z11 from the rejected starlight: the spatially filtered form of the same sensor |
~ Also read: PIAA-ZWFS (2026, arXiv 2606.28136), lossless pupil apodization that closes the gap to the fundamental sensitivity limit by 10× — a design lever, not a model change.

## Conclusions: what the corrected model and the exact reading settled | The model correction was the bigger lever; the primed iterated reading is the best one-frame reading; colour and photons are not the budget
::: full
- **The sensor model must keep the pupil DM-conjugate:** the diffraction bracket's two spheres share one radius, and the runner checks the round trip every time.  The correction alone took the single-actuator floor from 744 to 67 pm and made the undetected grid case detectable.
- **The exact reading with an iterated reference wave is the import that paid** (paper 1): 26 pm from one frame on a 30 nm working surface, SNR 355 — with the branch resolved once by four stepped frames of the base, refined with the iterated reference.  The plain stepped prior misses 3% of pixels and can flip a single actuator's sign.
- **Sampling is closed:** camera pixels per actuator (the ray count) set the hold-out gain; the dimple's own sampling did not move an actuator-space number at 4 px against 8 px.  The remaining 6.5% of hold-out gain is the next thing to measure — the kernel at the held-out site.
- **Colour and photons are not levers:** the chromatic null was the defocus; 10¹⁴ photons per state reaches 1 pm.  The budget is systematics: the 25–60 pm floors and gain stability.
~ Runner and records: templates/40_benches/zwfs_dm96 (README S7/S8; runs/rec193, ng385, ng385s3, m2048, rec193full).

## Conclusions: the follow-on imports, in order | The JPL budget form, model-matched DM calibration, the vector sensor priced before it is built — and the interferometer's own lever
::: full
- **The JPL budget form:** a per-pixel sensitivity factor from the measured fields predicts the photon-noise stage analytically, and the four systematic terms (pupil calibration at ½, reference wave at 1, dimple depth, initial phase) become the rows of the error budget.
- **Calibrate the DM through the sensor's own image model:** fit actuator positions and gains by matching sensor images of poke grids to a simulation that includes the propagation distances — the out-of-conjugate case, and the named suspect for the hold-out gain.  The interferometer's detector sits off the DM conjugate after the null-tuned tail, so this is its lever too.
- **Colour as a joint fit, not a combination:** the equal-weight combination inherits a bad channel; a joint per-wavelength fit through the exact reading would use 700–780 nm for range and 632.8 nm for the calibration of record.
- **Model the vector sensor before the metasurface is built:** the exact two-image solution and its leakage terms (retardance offsets, splitter rotation) run on the engine's polarization machinery; the bench limits at Keck and SEAL were defocus and crosstalk between the two pupil images, both priceable here.
- **The interferometer, for balance:** a PZT phase shift buys no model gain — its floors are geometric and colour-independent — but on hardware an absolute phase scale and freedom from polarization systematics; the recommended form is a hybrid: polarization snapshot for the differential measurements, a PZT on the reference flat as calibrator.  Its reflective (OAP) variant is in build.
~ Segmented-mirror lessons carried forward: the sensor reads segment piston but not global piston, DM quantization appears as steps in the response, and piston is underestimated below 50 % Strehl (Keck) — all three enter the e5-class segmented gauges.

# Backup

## Before the correction: the first poke and colour results | Kept for the record; every number on this slide is from the defocused sensor
::: full
![One actuator pushed 20 nm on the defocused sensor: the sensed blob is broadened and dimmed — gain 0.445, 0.29 nm rms error; the ring is the 4.9 m Fresnel defocus.](figs/zwfs_poke_triptych.png){h=1.85}
![Defocus at 8 nm amplitude on the same sensor: sensed matches applied at gain 0.986, 0.16 nm rms — low orders were never the trouble.](figs/zwfs_defocus_triptych.png){h=1.85}
- **The sampling sweep then read the spot as a band-select lever** (poke gain 0.45 / 0.59 / −0.38 at spots 2.0 / 3.0 / 1.06) and the kernel calibration recovered a 0.90 hold-out gain; both were measured on the ringed kernel.  On the corrected model the raw kernel peak is 0.75 (linear) and 0.93 (exact) and the spot only moves the low-order dip.
- **The colour claim** ("five colours fill the blind band; single-actuator floor 720 → 224 pm") was the Talbot null moving as 1/λ.  On the corrected model the null is absent and colour is not a lever (main slide 9).
~ Records: zwfs_s1 … zwfs_s6color reports and zwfs_sweep_report.txt; the legacy emission ('nf_legacy') reproduces them byte for byte.

## Provenance and records | Every number re-derives from a committed script
::: full
- Campaign: templates/40_benches/zwfs_dm96 — the runner (zwfs_params.m, zwfs_run.m, zwfs_run_figs.m, zwfs_run_batch.m, zwfs_batch.sh, macos_param_2048.txt) and its runs (runs/rec193 = the S7 equivalence check; ng385, ng385s3, m2048 = the sampling configurations; rec193full = colour + noise); the stage scripts zwfs_s1 … zwfs_s7iter as the historical record; README with the findings chain (S1–S8); plan and rulings in macos/BRIEF_zwfs_campaign.md.
- Findings earned in S1, recorded so they are not re-derived: the stock detector leg is ray-traced (geometric) at the mask plane — a λ/D-scale mask needs the reference-sphere diffraction bracket (a builder option, default off, existing decks bit-identical, suites green); the mask must sit at the optimized lens's real focus (−5.58 mm from the thin-lens seed); the camera-to-DM frame must come from traced rays — a support-area estimate was 25% off.
- The S7 correction: the bracket's exit sphere now carries the entrance sphere's radius ('nf'); 'nf_legacy' reproduces the S1–S6 emission; the interferometer twin is geometric at the mask and unaffected.
- Sign convention checked by a negative control: surface height = −φ·λ/4π on this deck (single reflection doubles height).
- Scoring: actuator space through the DM influence model, one library for both instruments (dm_gauge_lib: registration, actuator fit, modal correction, both measurement factories).
- Hardware source: "VSG2 Zernike Wavefront Sensor Update -v2" deck, parameters carried in templates/40_benches/vsg_wip/vsg2_params.m §9.
- Interferometer comparison column: templates/40_benches/tg_psi_dm96, run 10, S4, S5 (tg96_report.txt, tg96_s4_report.txt, tg96_s5noise_report.txt).
- Literature scan (2026-09-09): macos/REPORT_zwfs_lit_scan.md — twelve papers, ranked, with the PZT verdict for the interferometer.
