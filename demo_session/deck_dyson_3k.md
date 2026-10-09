<!--
deck_dyson_3k.md -- the long-slit telescope template (tma_longslit), established
on the 3k Dyson module (R9 feeding the 240 mm CaF2 Dyson of record) and applied
to the 1.5k (R5 ffc, the CLEAR row added on the way).  Parallels
deck_dyson_record.md (2026-10-06, sent to Jim and Joe).  Dave 2026-10-08: cast as
establishing a template; TO's 1.5k as its application, deviations noted as they occurred.
DRAFT -- pending Dave's sign-off.  Build: python3 make_brief_slides.py deck_dyson_3k.md
Sources: macos/REPORT_dyson5_cprime.md; templates/10_telescopes/tma_longslit/README.md;
challenges/dyson5/dyson5_cprime_3k.in + dyson5_t5f_cprime_*; BRIEF_to_dyson5.md
addenda 47-49; the SBG VSWIR paper (Bradley et al., ICSO 2024, Proc. SPIE 13699).
Figures are the runners' own PNGs; the figs_dyson/*_rc.png copies are recomposed at the panel level (recompose.py: trimmed; the 7-panel maps strip in two rows), never re-rendered.
-->

# A Long-Slit Spectrometer Design Challenge
A VSWIR imaging spectrometer for a 180 km swath at 30 m pixels, designed end to end by ray trace: a Dyson spectrometer and a freeform three-mirror telescope per module, two module sizes, both at the pixel — smile 0.02, CRF 1.2 px, 0.85–0.97 of the energy in one pixel — and the telescope flow made a reusable template.
D. C. Redding, with Claude Code.
October 2026.  Working record; no proprietary prescription is used.
DRAFT — pending review.
~ The challenge was posed by Joe Green and Jim McGuire (JPL) as a test of MACOS and its design layer; this is its second round.  The first (2026-10-06) brought the spectrometers to the pixel and left the larger telescope short; their answers and the project's published design set this round's course.

## Contents | The challenge and the approach; the pair of 3k modules; the template's default run — seed, telescope of record, the module scored, how the ladder was set, the two smile conventions, the cone; the spectrometer; the template applied to the 1.5k; the comparison; the questions
::: left
- **3** The specification, and how a design is scored against it
- **4** The challenge: two module sizes, a spectrometer and a telescope in each
- **5** The approach: the 3k telescope solved as a template, the template run on the 1.5k
- **6** Two 3k modules fill the swath
- **7** The two strips on the ground, and the two focal planes
- **8** The seed: the SBG form, first order
- **9** The telescope of record: freeform three-mirror, R9
- **10** The 3k module, scored: spectrometer, telescope, and the two together
- **11** The 3k module, end to end
- **12** How the optimization "ladder" was set: six rungs, three lessons
::: right
- **13** Smile is a convention: the launch, not the optics
- **14** The cone, the spec and the paper's beam
- **15** The spectrometer: toleranced, and a conic grating for margin
- **16** Throughput by band, and the blaze
- **17** The template applied: the 1.5k module, with the deviations as they occurred
- **18** Two of 3k or four of 1.5k, updated
- **19** Questions for Jim and Joe
- **20** Next steps
- **21** Run it yourself
- **22** Backup: the first order by the numbers; the rungs per field; the chief-tilt diagnostic; the tolerance ladder by the numbers; the engine finding; records

## The specification, and how a design is scored against it | The challenge's numbers beside the project's published design (the SBG VSWIR paper), the metrics defined once, the paper's bounds converted to this work's pixel
::: left
| item | the challenge (this work) | the SBG VSWIR paper |
|---|---|---|
| F-number | 1.8 | 1.8 |
| aperture; focal length | 183 mm; 330 mm | 192 mm; 345 mm |
| pixels | 18 µm; 3000 cross-track, 500 spectral | 18 µm; 3072 × 512, 2 co-added spectrally |
| slit | 54 mm; 2 pixels (36 µm) | 54.45 mm; 36 µm |
| band | 380–2500 nm | 380–2500 nm |
| smile, keystone | < 0.1 px (1.8 µm) | 1.8 µm each (5 % of a co-added px; 10 % of a px) = 0.10 of our pixel |
| SRF FWHM | 1.5–2.0 px (27–36 µm) | < 1.8 co-added px (64.8 µm) |
| CRF FWHM | < 1.5 px (27 µm) | < 2.8 px (50.4 µm) |
| ARF FWHM (telescope, across the slit) | – | < 2.8 px (50.4 µm) |
| field per module | 9.4° (the 3k module); 4.7° (the 1.5k) | 8.9° |
::: right
- **Smile:** centroid drift along the slit at fixed wavelength; **keystone:** drift across wavelength at fixed slit position; both from ray centroids, in pixels.
- **SRF, CRF:** slit ⊗ line-spread function ⊗ pixel ⊗ Airy, FWHM (Mouroulis & Green 2018).  With a 2-pixel slit a perfect spectrometer gives SRF = 2.010–2.023 px over the band (slit, pixel and diffraction): the floor.
- **Energy in one pixel:** the fraction of a point image inside one 18 µm pixel; design rule 0.75.
- **End to end** is one prescription from the sky to the detector, the grating its stop, scored by the spectrometer's own scorer.  Smile end to end depends on how each field is launched onto the slit — slide 13 — so it is given under both conventions.
- **The paper's SRF is in co-added pixels** (36 µm); its bounds are converted to micrometers in every table here.
~ Every number is a MACOS ray trace of a committed prescription; no coatings, no tolerances.  Paper: Bradley et al., "Finalized optical design of the SBG VSWIR Wide Swath Imaging Spectrometer," ICSO 2024, Proc. SPIE 13699, 1369945 (Tables 1–3).

## The challenge: two module sizes, a spectrometer and a telescope in each | 6000 cross-track pixels over a 180 km swath, 380–2500 nm, F/1.8, at the pixel in smile, keystone and response width; two 3k modules or four 1.5k; every number a MACOS ray trace of a committed prescription
::: left
- **The instrument:** a push-broom imaging spectrometer — a straight slit imaged on the ground, swept along the track; the spectrometer disperses the slit across a 500-pixel band at 18 µm pixels.  The challenge fixes F/1.8, a 330 mm focal length, a 183 mm aperture, and the image-quality bounds of slide 3.
- **The two module sizes.**  A **3k module** has a 3000-pixel, 54 mm slit and needs a 240 mm CaF2 Dyson block; two of them cover the swath.  A **1.5k module** has a 1500-pixel, 27 mm slit and a 130 mm fused-silica block; four are needed.  Which to build is the trade the challenge asks for.
- **Each module is a Dyson spectrometer and a telescope.**  The Dyson (a concentric block with a concave grating on its curved face) was designed first and brought to the pixel in the first round; the telescope images a strip of sky onto the slit, telecentric there, with no marginal ray faster than the spectrometer accepts.
::: right
- **Where the first round left it (2026-10-06):** both Dysons at the pixel; the 1.5k telescope (three off-axis aspheres) near it; the 3k telescope not — CRF 4 px, 4–10 % of the energy in a pixel, marginal rays at F/1.2 into an F/1.8 spectrometer.
- **The answers that set this round's course:** two modules, not four — integration, test and calibration are the cost, not mass or glass, and 3k × 0.5k detectors are a standard format where 1.5k × 0.5k is not; design the telescope separately, telecentric at the slit, a little faster than the spectrometer, no F/1.2 rays; the detector window is an order-sorting filter in air; tailor the coatings and the grating's blaze to the photon-starved long end.
- **The project's published design** (Bradley et al., ICSO 2024): a freeform zig-zag three-mirror telescope at f 345 mm, 192 mm, 8.9°, with a small convex M2 at the stop and no intermediate focus, optimized alone for small spots along the slit; a Dyson with an aspheric CaF2 lens and a conic grating (their Option B); its bounds in micrometers (smile and keystone 1.8 µm; SRF 64.8 µm; CRF and ARF 50.4 µm).
~ The challenge and the first-round answers are recorded in `BRIEF_to_dyson5.md` (addenda 47–49); the first round's deck is `deck_dyson_record` (the spectrometers' design record).  The paper: Tables 1–3, Figs. 1 and 4.

## The approach: the 3k telescope solved as a template, the template run on the 1.5k | The published form became a parameterized design flow — one parameter file, four stages, the constraint rows a slit-fed spectrometer needs — whose default run is the 3k module; run at the 1.5k spec it reached the pixel in five rungs with one new row; the slides that follow are those two runs
::: left
- **The 3k telescope (slides 8–14):** the published layout digitized and scaled to the challenge's focal length is the seed; the first order is closed-form telecentric with the stop at M2; three off-axis conic sections carry it; a pole-frame polynomial of degree 3–6 on each mirror brings every field of the 9.4° strip to a 1.02-pixel spot with the cone bounded; joined to the 240 mm CaF2 Dyson it images at the pixel end to end.  Six rungs taught the weights and rows the flow now carries; the last tenth of a pixel of "smile" was the scorer's launch convention, not the optics.
- **Made a template** (`templates/10_telescopes/tma_longslit/`): the spec, the seed legs and incidences, the bounds, the weights, the field set, the ladder and the spectrometer to join live in one parameter file; four stages — `first_order` (closed form), `section` (the conics on the chief, measured per field), `figure` (the ladder, each rung an engine-traced solve warm from the last), `e2e` (the best clear rung joined to the spectrometer under both launches) — each writing its table, record and prescription; the deck of record pinned by a test that re-scores it by name.
::: right
- **The rows that make it a long-slit front end:** spot size along the slit (weighted); plate scale; the chief's and the centroid's across-slit bow; the chief angle to the slit normal (telecentric); a hinge on the marginal F/# outside the bound, per axis (the cone); the chief − centroid offset at the slit; working distance; and, since the 1.5k run, a clearance wall.  The first order is re-derived at every iterate, never penalized, so focal length, back focus and telecentricity hold exactly under the figure.
- **The template on the 1.5k (slide 17):** the same seed, rows and weights at a ±2.35° strip and a 27 mm slit, joined to the 130 mm silica Dyson: 1.02 px at every field in five rungs; one deviation on the way — the first rung at the floor put the M3 → slit beam 0.3 mm inside M2, and a clearance wall row bought +2.4 mm with the floor held; end to end smile 0.009, CRF 1.03, 0.97 of the energy in a pixel.
- **Then (slides 15–18):** the spectrometer toleranced and its Option B assessed, throughput by band, and the two-of-3k trade restated with both modules at the pixel.
~ `BRIEF_to_dyson5.md` addenda 48–49; `REPORT_dyson5_cprime.md`; the template README.  The dyson5 join (`dyson5_t5f`) is the spectrometer side's scorer; the template calls it, it does not own it.

## Two 3k modules fill the swath | Side by side across the track at a 324 mm pitch, each canted ±4.68° to its own strip: no body touches another and no body sits in the other's sky beam; the pair is 0.75 × 0.60 × 1.36 m, 609 L, with 39 kg of optics
::: full
![The two modules as the engine traces them, each the 3k module of record moved rigidly: cants ±4.68° about the axis normal to the track and the line of sight.  Left: 3-D; right: looking along the cant axis, light entering from below.  The rays the grating stops are red.](figs_dyson/dyson5_block2_3k_views_rc.png){h=3.7}
| item | value |
|---|---|
| pitch across the track; worst margins | 324 mm (module 319 mm + 20 mm gap, in 5 mm steps); body to body +24 mm, body in the other's 183 mm sky beam +36 mm |
| envelope of the pair | 746 × 598 × 1363 mm, 609 L (the four 1.5k modules: 965 × 411 × 758 mm, 300 L) |
| optics per module; for the swath | 19.5 kg (block 14.3, three mirrors 2.9, grating 2.3); 39.0 kg (the four 1.5k modules: 14.6 kg) |
~ Bodies are the engine's surface hits + 5 mm aperture + 10 mm mount; mirrors and grating as 10 mm blanks (Zerodur, silica), the CaF2 block from `dyson5_jim_3a.txt`.  Not included: structure, detectors and their cooling, electronics, baffles.  Each module's internal clearances are the record's by construction.  Record `dyson5_block2_3k.txt`; tool `dyson5_block4.m` (`nmod` 2).

## The two strips on the ground, and the two focal planes | Each module's slit admits ±4.691°, 3000 pixels; canted ±4.68° the strips overlap by 10 pixels at the join and cover 18.7°, 181 km from 550 km, with no along-track offset; each detector holds its strip across 3000 spatial pixels and the band across 500
::: full
![Top: the swath, cross-track angle across, ground distance from nadir in gray, the dark band the 10-pixel overlap, the arrow the along-track sweep.  Bottom: the two 3000 × 500 detectors under their strips, the lines the slit's image at the seven scored wavelengths from the engine's centroids (0.24 px/nm; 380 nm at +v, 2500 nm at −v; the image is inverted along the slit).](figs_dyson/dyson5_block2_3k_swath_rc.png){h=4.4}
- **One join instead of three.**  The plate scale puts 3000 pixels over 9.381°; abutting at 9.4° would leave a 6-pixel gap, so the cants give 10 pixels of overlap.  **No along-track offset:** both lines of sight lie in one cross-track plane, so the strips image the same ground line at the same time.  **Unique pixels:** 5990 over 18.73°, 30.1 m at nadir; the detector is the paper's 3072 × 512.
~ Record `dyson5_block2_3k.txt` (strip edges per module, overlap, swath); the admitted strip from `dyson5_t5f_cprime_centroid_roll000.txt`.

## The seed: the SBG form, first order | The published layout digitized and scaled to the challenge's focal length; with the stop at M2 the first order is telecentric in closed form, has no intermediate focus, and gives 299 mm of working distance for the Dyson
::: left
![The telescope section as the engine traces it (meters): light enters from below; M1 (E1) at the right, M2 (E2) at the top, M3 (E3) at the left, the slit (E4) at the far left — a zig-zag, not a Korsch fold.](figs_dyson/dyson5_cprime_3k_viewyz_rc.png){h=4.6}
::: right
- **The seed:** the paper's Fig. 4b, read off its 200 mm scale bar and scaled by 330/345 — legs M1→M2 258 mm, M2→M3 246 mm, M3→slit 300 mm; chief incidence 30°, 37°, 15°; M2 the small mirror at the stop.
- **Two facts the first order forces.**  (1) **A tilted sphere cannot do it:** the legs are the same in and across the fold plane, so each mirror needs local radii at the chief with R_t / R_s = 1/cos² i (1.33, 1.57, 1.07); an off-axis conic section does this natively, and with the parent axis at the incidence angle all three are off-axis paraboloids — the seed.  (2) **The stop at M2 makes telecentricity closed-form:** M2 sits at M3's front focus.  M1 weak concave (f +1012 mm), M2 convex (−452), M3 f +246 = the M2→M3 leg; the marginal ray never crosses the axis between mirrors, so there is no intermediate image — the Korsch with a real M2–M3 focus cannot be telecentric with the stop at M2, and this family has no such focus.
- **In the engine, before any figure:** chiefs within 0.093° of the slit normal over ±4.7°, plate scale 329.99 mm, every field admitted, M2 footprint 136 mm, clearance +10 mm, working distance 299 mm; not imaging (2 mm spots).
~ Template `templates/10_telescopes/tma_longslit/` (parameters, runner, demo, README, test); first order `tls_first_order.txt`; `BRIEF_to_dyson5.md` addendum 48.  Beams: 183 mm at M1, 136 at M2, 166 at M3.

## The telescope of record: freeform three-mirror, R9 | The template's default run.  Off-axis conic sections with a pole-frame polynomial of degree 3–6 on each mirror: FWHM 1.02 × 1.02 pixels at every field of the strip, chiefs within 0.033° of the slit normal, marginal rays between F/1.80 and F/1.89, clearance +12 mm
::: full
| field (±) | spot rms (µm) | FWHM along / across the slit (px) | energy in a pixel | chief to the slit normal | F/# along / across | chief − centroid across the slit (µm) |
|---|---|---|---|---|---|---|
| 0 | 1.8 | 1.02 / 1.02 | 1.00 | 0.000° | 1.894 / 1.806 | 0.00 |
| 1.175° | 2.2 | 1.02 / 1.02 | 1.00 | 0.002° | 1.894 / 1.803 | – |
| 2.35° | 2.5 | 1.02 / 1.02 | 1.00 | 0.008° | 1.893 / 1.801 | – |
| 3.525° | 2.0 | 1.02 / 1.02 | 1.00 | 0.018° | 1.892 / 1.798 | – |
| 4.7° | 2.3 | 1.02 / 1.02 | 1.00 | 0.033° | 1.889 / 1.797 | 0.03 (R7: 2.43) |
::: left
- **Geometry:** the seed's legs and incidences kept; the layout stays the form (unbounded, the incidence angles collapse to a coaxial train that self-obscures).  Plate scale 329.9 mm local, 330.1 along track.  Footprints M1 235 × 202, M2 132 × 166, M3 215 × 173 mm; freeform departure 0.46 / 0.55 / 1.11 mm P-V over the lit patch (the 1.5k asphere telescope: 2.1 mm on M1).  Clearance +12.1 mm; working distance 299.4 mm.
- **Figure:** a pole-frame polynomial (`Surface= Monomial`, degree 3–6, even in x); degree ≥ 3 leaves value, slope and curvature at the pole untouched, so the first order survives the figure.  Even aspheres about the parent axis were the wrong basis on these sections (the poles sit 0.9–1.5 m off axis).
::: right
- **The constraint rows that make it a long-slit front end, in every rung:** the chief angle to the slit normal per field (telecentric); a hinge on the marginal rays' F/# at the slit, along and across, outside [1.7, 1.8] (the paper's F/# anamorphicity); and, in R9, the chief − centroid offset across the slit (slide 13).
- **Not done:** no tolerancing, no coatings, no stray-light or mirror sizing beyond the footprints.
~ Prescription `challenges/dyson5/dyson5_cprime_3k.in` (= `tls_R9_ffo.in`); engine numbers, model 256, nine fields over ±4.7°, mirror-symmetric.  The gate `tTmaLongslit` re-scores it by name.  `REPORT_dyson5_cprime.md`.

## The 3k module, scored: spectrometer, telescope, and the two together | Both meet the challenge's specification; together they image at the pixel over the strip; SRF is the two-pixel slit's floor
::: full
| criterion | challenge | paper | spectrometer alone (240 mm CaF2 Dyson) | telescope alone (R9) | together, end to end |
|---|---|---|---|---|---|
| smile | < 0.1 px | < 1.8 µm (0.10 px) | 0.005 px | – | 0.022 px slit-filled (0.4 µm); 0.016 px point-source shift |
| keystone | < 0.1 px | < 1.8 µm (0.10 px) | 0.006 px | – | 0.006 px (0.1 µm) |
| CRF FWHM | < 1.5 px | < 50.4 µm | 1.21 px | 1.02 px along the slit | 1.17 px (21 µm); 1.26 at roll 180° |
| SRF FWHM | 1.5–2.0 px | < 64.8 µm | 2.024 px (2-px slit) | 1.02 px across the slit (ARF 18.4 µm) | 2.025 px (36.5 µm), slit-limited |
| energy in one 18 µm pixel | > 0.75 | – | 0.82 | 1.00 | 0.85 (0.81 at roll 180°) |
| image blur (rms, as placed) | – | 2–3 µm (paper) | – | 1.8–2.5 µm | – |
| plate scale | 330 mm | 345 | – | 329.9 / 330.1 mm | 330.1 mm |
| chief rays at the slit | telecentric | telecentric | – | within 0.033° | – |
| marginal rays at the slit | – | F/# anamorphicity bounded | – | F/1.80–1.89, none below 1.7 | grating admits 98.9 % |
| clearance, worst pair | > 0 | – | +0.54 mm (slit vs block face) | +12.1 mm | +0.6 mm (the spectrometer's pair) |
| glass, mass | – | – | 4.5 L, 14.3 kg | mirrors not yet sized | – |
~ Spectrometer: `dyson5_jim_3a.txt` (size:F:240); telescope: `dyson5_cprime_3k.in`, `REPORT_dyson5_cprime.md`; together: `dyson5_t5f_cprime_centroid_roll000.txt` / `_roll180`, `_chief_roll000` (the point-source shift), 7 fields × 7 wavelengths.  No coatings, no tolerances in any column.  SRF 2.025 against a floor of 2.010–2.023 (slit, pixel, diffraction): the Dyson adds at most 0.0015 px, 0.03 µm.

## The 3k module, end to end | The freeform telescope feeding the 240 mm CaF2 Dyson: smile 0.02, keystone 0.006, CRF 1.17, SRF 2.025 px at every slit position and wavelength; 98.9 % of the light through the grating; clearance +0.6 mm, the Dyson's own slit-to-block-face margin
::: left
![As the engine traces it (meters), the dispersion plane: the telescope's M1–M3 (E1–E3) below, the slit (E4) on the block's face at the origin, the 240 mm CaF2 block, the grating (E7) at the top, the detector (E9) beside the slit.  The 1 % of rays the grating stops are drawn in red.](figs_dyson/dyson5_t5f_cprime_e2e_viewyz_rc.png){h=4.6}
::: right
![Per slit position and wavelength (7 × 7): field angle, keystone, smile, the two response widths and the energy in one pixel, uniform over the strip.  The runner's own figure.](figs_dyson/dyson5_t5f_cprime_centroid_roll000_maps_rc.png){h=3.4}
- **CRF 1.02 px at the inner fields and 1.17 at the strip ends; SRF 2.02 everywhere; smile 0.022 → 0.018 px from 380 to 2500 nm; keystone 0.006 px at the ends, 0 at the center.**  Against the (c) solve of the last deck (smile 1.64, CRF 4.02, SRF 4.99, 4–10 % in a pixel), every line moved to the pixel.
~ Record `dyson5_t5f_cprime_centroid_roll000.txt` (each field's centroid on the slit line); the chief-launched twin `_chief_roll000.txt` gives the point-source shift 0.016 px; roll 180° about the chief: CRF 1.26, energy 0.81, the rest unchanged.  Prescription `dyson5_t5f_cprime_centroid_roll000_e2e.in` (3 telescope + 8 slit/Dyson elements).

## How the optimization "ladder" was set: six rungs, three lessons | The 3k run fixed the template's weights and rows: equal weights spent the budget on the strip center; a constraint set weighted below one percent of the merit is not a test; the last tenth of a pixel of "smile" was the launch
::: full
| rung | what changed | outer-field FWHM along the slit (px) | e2e smile, point-source (px) | keystone (px) | CRF (px) | SRF (px) |
|---|---|---|---|---|---|---|
| R4 | conics + even aspheres on the seed layout | 5–7 | 0.73 | 0.009 | 2.5 | 3.9 |
| R5 | pole-frame freeform, centroid-bow rows added | 5–7 | 0.09 | 0.03 | 7.0 | 2.5 |
| R6 | R5 + 3000 more evaluations | 5–7 | 0.15 | 0.04 | 7.2 | 2.5 |
| R7 | warm from R4: along-slit rows ×3, outer fields ×2–3, bow rows kept | **1.02** | 0.14 | 0.006 | **1.18** | 2.025 |
| R8 | R7 + telecentric rows ×30 (a test of the chief angle) | 1.02 | 0.14 | 0.005 | 1.18 | 2.025 |
| **R9** | R7 + chief − centroid rows across the slit ×1000 | **1.02** | **0.016** | 0.006 | **1.17** | 2.025 |
- **Equal weights were the CRF wall.**  R5 and R6 weighted every field and both axes equally on rms spot size and spent it on the strip center (4 µm) while the outer half sat at 5–7 px along the slit; more budget on the same merit moved nothing.  Weighting the along-slit rows and the outer fields, from R4, took every field to the Airy floor in one rung.
- **A row set below one percent of the merit is not a test.**  R8's telecentric rows at ×30 were 0.4 % of the merit; the solver never worked on them and R8 is R7 to the third digit.  R9's offset rows at ×30 would have been 0.07 %; at ×1000 they dominated, and the solve converged in 507 evaluations.
- **The remaining 0.14 px of "smile" was the launch convention** (slide 13), and the row that removed its cause cost nothing at the pixel.  **What the template carries forward from this:** R9's weights and rows (along-slit ×3, outer fields ×2–3, the bow rows, OFF ×1000) from the first rung, and a ladder without the detours — slide 17 runs it that way.
~ `REPORT_dyson5_cprime.md` (every rung and the four solves thrown out); the template's README, "What the ladder taught."  Each rung: 3000 evaluations of the engine's own trace (R9: 507), model 256.

## Smile is a convention: the launch, not the optics | The join lands each field's chief ray on the slit line and scores the centroid, so the telescope's coma — the chief-minus-centroid offset across the slit — read as smile; the Dyson fed with the telescope's chief angles does not move
::: left
| Dyson alone, each slit point's chief tilted by R7's angles | smile | keystone | CRF | SRF |
|---|---|---|---|---|
| base (the record) | 0.0050 | 0.0058 | 1.213 | 2.024 |
| along-track component (−0.51 mrad at the strip end) | 0.0049 | 0.0061 | 1.213 | 2.024 |
| cross-track component (±0.27 mrad) | 0.0048 | 0.0063 | 1.213 | 2.024 |
| along-track × 2 | 0.0049 | 0.0064 | 1.213 | 2.024 |
| both | 0.0048 | 0.0066 | 1.213 | 2.024 |
- **The chief angle is not the cause:** keystone moves (the hook is live), smile does not, even at twice the angle.
::: right
| rung | predicted from the telescope's offset | e2e smile, chief launch |
|---|---|---|
| R4 | 0.655 | 0.732 |
| R5 | 0.086 | 0.090 |
| R6 | 0.122 | 0.149 |
| R7 | 0.135 | 0.140 |
| R8 | 0.133 | 0.138 |
- **The offset predicts every rung:** on R7 the chief sits 0 / 0.2 / 0.7 / 1.4 / 2.4 µm from the centroid across the slit at 0 / 1.2 / 2.4 / 3.5 / 4.7°; the spread of that offset is the e2e smile.
- **Two conventions, both reported.**  *Slit-filled* (the paper's: its telescope-fed smile equals its spectrometer-alone smile): each field's centroid on the slit line — the spectrometer's smile, 0.022 px.  *Point-source*: the chief on the slit line — where a star's spectrum lands, 0.016 px on R9 (0.14 on R7).  R9's row on the offset closed the second without touching the first.
~ Diagnostic `tls_dyson_chief_tilt` (the spectrometer alone, the launch tilted per field); the join scorer `dyson5_t5f` now carries both launches (`tel5f_launch` chief | centroid), the record's rows reproduce under the first.

## The cone, the spec and the paper's beam | The marginal rays leave the telescope between F/1.80 and F/1.89 in the two directions, none below F/1.7; the along-slit F/1.89 is Joe's own aperture at Joe's focal length
::: left
- **Why it is a row:** the last deck's 3k solve let its marginal rays run to F/1.19 into an F/1.8 spectrometer; the grating, sized for F/1.8, passed 95–98 % of the light, and nothing in the merit knew.  Jim: "I don't see why some rays correspond to an f/1.2 cone."  The paper constrains the F/# anamorphicity.
- **The row:** a hinge on the working F/# of the marginal rays at the slit, along and across, outside [1.7, 1.8], in every rung.  R9: F/1.889–1.894 along the slit, F/1.797–1.806 across; the grating admits 98.9 % at every field.
::: right
- **Along the slit it sits above 1.8, and that is the specification:** D 183 mm at f 330 mm is F/1.803; a cone a little faster than the spectrometer's — Jim's "a little faster" — needs a larger entrance beam, the paper's 192 mm.  A spec question for Jim and Joe, not a solve question; the telescope holds Joe's numbers.
- **Across the slit** the cone is F/1.80: the slit truncates the image in that direction, and the across-slit FWHM (the paper's ARF) is 1.02 px against its 2.8 px bound.
~ `REPORT_dyson5_cprime.md`, R9 per-field table; F/# from the marginal rays as 1/(2 sin u).

## The spectrometer: toleranced, and a conic grating for margin | SRF sits at the two-pixel slit's floor; every large tolerance is a focus term the detector's focus removes; keystone sets the alignment tolerances at 26 µm and 35 µrad; a conic on the grating takes the energy in a pixel from 0.82 to 0.99
::: left
- **SRF is not a miss.**  Slit ⊗ LSF ⊗ pixel ⊗ Airy for a perfect spectrometer is 2.010 px at 380 nm and 2.023 at 2500 nm; the Dyson of record sits 0.0001–0.0015 px above that at every wavelength, alone and behind the telescope.  In the paper's units 36.5 µm against 64.8.  Joe's "1.5–2.0" is met at its floor; 1.5 would need a narrower slit.
- **The Dyson of record is the paper's Option A:** 240 mm CaF2, an aspheric lens (K 0.043 with h⁴, h⁶) and a spherical grating, no meniscus; smile 0.005, keystone 0.006, CRF 1.21, 0.82 of the energy in a pixel; 4.5 L, 14.3 kg.
- **Toleranced, 21 rows, alone and through the join:** every large sensitivity is a focus term — lens radius and thickness, the air space, slit and detector defocus, grating radius (the strongest: dR/R 10⁻⁴ costs CRF +1.9 px uncompensated, and a 76 µm detector refocus removes it to 0.006).  One-axis detector focus is the compensator; x/y adds nothing.  Smile is loose (≥ 300 µm or µrad).
::: right
- **Keystone sets the tolerances, and no compensator helps** (it is lateral color).  For the paper's as-built keystone (1.03 µm = 0.057 px): lens and grating decenter along the slit 26 µm, grating tilt 35 µrad, clocking 88 µrad; for Joe's 0.1 px: 47 µm, 47 µm, 64 µrad, 162 µrad.  The 1.5k Dyson is twice as sensitive (13 µm).  Linear where it matters (keystone vs decenter 1.10; CRF vs grating radius 1.03).
- **Thermal:** 1 K of CaF2 index is negligible (CRF +0.001); CaF2's radius and thickness terms cancel (+0.050 / −0.051 px per K); the open terms are the aluminum air space (−0.22 px per K) and the N-BK7 grating radius (+0.14) — focus terms, refocus-compensable, which is the paper's "CRF limits the spectrometer's temperature requirements."
- **Option B — a conic on the grating (K −0.004) with the lens re-solved:** alone CRF 1.03 (1.21), energy in a pixel 0.99 (0.82), keystone 0.002 (0.006); behind R9: CRF 1.06, energy 0.955, smile 0.012; SRF the floor in both.  Its keystone sensitivities are 6–45 % higher than A's (clocking most): **it buys nominal margin, not alignment insensitivity** — the paper's preference for B on tolerances is not reproduced by alignment alone (its as-built numbers include fabrication and thermal, not modeled here).  Recommended as the spectrometer candidate for its margin.
~ `REPORT_dyson5_spec.md` (the ladder, 3k and 1.5k, Option B and its ladder); `spectrometer_sens`, gate `tSpectrometerSens`.  The compensator residuals are given for detector focus alone and focus + x/y.  A full every-CTE soak is the next row.  Throughput by band: next slide.

## Throughput by band, and the blaze | Six air-glass crossings pass 0.82 uncoated in every band; a single quarter-wave MgF2 centered at 2.2 µm lifts the long band to 0.87 and costs nothing elsewhere, where a two-layer tuned to 1 µm trades the band ends for the middle; the grating's blaze is the lever that moves photons to the long end
::: full
| module, coating on all six crossings | 380–700 nm | 700–1300 nm | 1300–2500 nm | 380–2500 mean |
|---|---|---|---|---|
| 3k CaF2, uncoated | 0.817 | 0.822 | 0.827 | 0.822 |
| 3k CaF2, MgF2 quarter wave at 2.2 µm | 0.852 | 0.844 | **0.869** | 0.855 |
| 3k CaF2, MgF2 / Al₂O₃ quarter-quarter at 2.2 µm | 0.797 | 0.801 | **0.880** | 0.826 |
| 3k CaF2, MgF2 / Al₂O₃ quarter-quarter at 1.0 µm | 0.749 | **0.920** | 0.736 | 0.802 |
| 1.5k silica, uncoated / MgF2 at 2.2 µm | 0.807 / 0.851 | 0.814 / 0.840 | 0.821 / 0.872 | 0.814 / 0.854 |
| grating, scalar blaze at 1.0 / 1.4 / 1.8 µm | 0.14 / 0.02 / 0.01 | 0.89 / 0.53 / 0.17 | 0.48 / 0.78 / 0.89 | – |
- **The crossings:** the slit face in, the block's convex face out to the grating's air gap and back in, the exit face — four on the block — and two on the detector's order-sorting filter, which sits in air (Jim: it cannot be bonded to CaF2).  Each entry is the mean over the band of the product of the six transmittances, the block's index from its Sellmeier equation at every wavelength.
- **Where to spend the coating:** on a 1.43–1.45 substrate a single MgF2 layer is never a good match, but centered at 2.2 µm it helps most where the photons are scarcest (+0.04–0.05 in the long band) and still helps at the short end; the two-layer tuned to 1 µm is the better coating in the middle of the band (0.92) and the worse one at both ends.  Jim's rule, in numbers.
- **The blaze does the real work at the long end:** a 1.8 µm blaze puts 0.89 of first-order efficiency in the long band against 0.48 for a 1 µm blaze, at the price of the short band, where the photons are plentiful.  The scalar estimate says where the light goes, not the absolute efficiency; a vector grating calculation is the next step there.
~ Record `dyson5_throughput_band.txt` (`dyson5_throughput_band.m`, `thinfilm_rt`); route 2 of the record = the 1.0 µm two-layer row.  No detector QE, slit loss or mirrors.

## The template applied: the 1.5k module, with the deviations as they occurred | The same seed, rows and weights at the 1.5k spec (±2.35°, 27 mm slit) reach the pixel in five rungs; one deviation forced a new row — the first rung at the floor entered M2's body by 0.3 mm, and a clearance wall bought +2.4 mm with the floor held
::: full
| rung | DOFs | spot rms (µm) | worst FWHM (px) | min energy in a pixel | clearance (mm) |
|---|---|---|---|---|---|
| R0 seed | first order | 1876–1968 | 15.6 | 0.00 | +10.1 |
| R1 conic | θ, slit focus | 151–187 | 10.8 | 0.00 | +17.1 |
| R2 asphere | + h⁴, h⁶ | 57–152 | 14.7 | 0.00 | +13.9 |
| R3 freeform 3–4 | pole-frame polynomial | 17–21 | 1.92 | 0.22 | +1.2 |
| R4 freeform 3–6 | degree 3–6 | 0.73–0.98 | **1.02** | 1.00 | **−0.3** |
| **R5 ffc** | R4 + the CLEAR wall | 0.76–0.96 | **1.02** | 1.00 | **+2.4** |
::: left
- **The same:** f, D, F/1.8; the seed's legs and incidences; the rows and weights the 3k ladder settled (slide 12), from R1 on — no equal-weight detour, no test rung.  **Changed by parameter only:** the strip, the slit, and the spectrometer joined (the 130 mm silica Dyson of record).
- **Deviation 1, expected:** aspheres are not enough at 1.5k either (R2 57–152 µm); the pixel needs the freeform, as on the 3k.
- **Deviation 2, the new row:** R4 reached 1.02 px at every field with its M3 → slit leg 0.3 mm inside M2's body (the 3k never came near: +12 mm).  CLEAR scores the clearance rule every iterate, a hinge to 2 mm at merit scale (weight 3000): R5 clears by +2.4 mm and keeps the floor.  The row is now in the template.
::: right
- **Deviation 3, an absence:** the two launches agree to 0.001 px.  With the OFF rows in from the first rung the chief − centroid offset that gave the 3k its 0.14 px never appeared; the 3k's R9 step was not needed.
- **Telescope alone (R5):** plate 329.98 / 329.75 mm, chiefs within 0.008° of the slit normal, F/1.89 along, 1.80 across.
- **End to end, roll 0 / 180:** smile 0.009 / 0.011, keystone 0.009 / 0.012, CRF 1.03, SRF 2.024 (the floor), energy in a pixel 0.97 / 0.99, grating admits 0.987, clearance +0.6 mm.  Joe's spec met under both launches except SRF, at the floor; the paper's five pass.  Against the three-asphere 1.5k (smile 0.08 / 0.50 point-source, SRF 2.21, energy 0.30–0.39) every line improves but keystone at roll 180: 0.012 against 0.009, a tenth of the bound.
~ `tma_longslit_1k5.m`; `tls1k5_figure.txt`, `tls1k5_e2e.txt`, `dyson5_cprime_1k5.in` (pinned by `tTmaLongslit`), `dyson5_rescore_1k5_launch.txt`; model 256, nine fields.

## Two of 3k or four of 1.5k, updated | Both modules now image at the pixel end to end; the trade is labor, detectors and integration, not optics or glass — and Jim's and the paper's answer is two
::: full
| | 2 × CaF2, 3k, 54 mm slit | 4 × fused silica, 1.5k, 27 mm slit |
|---|---|---|
| block radius = thickness | 240 mm | 130 mm |
| spectrometer alone: smile / keystone / CRF / SRF / energy in a pixel | 0.005 / 0.006 / 1.21 / 2.02 px / 0.82 | 0.007 / 0.009 / 1.03 / 2.02 px / 1.00 |
| telescope | freeform three-mirror of the SBG form, the template's default run (R9): 1.02 px at every field | the template at the 1.5k spec (R5 ffc, slide 17): 1.02 px at every field |
| end to end: smile (slit-filled / point-source) / keystone / CRF / SRF | 0.02 / 0.02 / 0.006 / 1.17 / 2.025 px | 0.009 / 0.009 / 0.009 / 1.03 / 2.024 px |
| end to end: energy in one pixel, worst field | 0.85 | 0.97 |
| light admitted by the grating | 98.9 % per field | 98.7 % per field |
| clearance, worst pair | +0.6 mm (slit package vs block face) | +0.6 mm (the same pair); telescope alone +2.4 mm |
| glass for 6000 pixels | 9.0 L, 28.5 kg CaF2 (12.4 L crystal to carve) | 3.0 L, 6.6 kg fused silica |
| the block of modules: envelope; optics | 0.75 × 0.60 × 1.36 m, 609 L; 39 kg | 1.01 × 0.53 × 0.92 m, 494 L; 21.4 kg (the three-asphere record's block: 300 L, 14.6 kg — the zig-zag's M2 and M3 are full aperture) |
| telescopes / spectrometers / detectors | 2 / 2 / 2 × (3072 × 512, a standard format) | 4 / 4 / 4 × (1.5k × 0.5k, NRE) |
| integration, test, calibration | two of each | four of each (Jim: the cost that dominates) |
| materials | "a tiny fraction of the instrument cost" (Jim) | the same |
~ Records `dyson5_jim_3a.txt`, `dyson5_t5f_cprime_centroid_roll000.txt`, `dyson5_block2_3k.txt`; `tls1k5_e2e.txt`, `dyson5_block4_cprime1k5.txt`.  The three-asphere 1.5k of the last deck is superseded (slide 17); its smaller block is the fallback.  The air gap at the window stands on both.

## Questions for Jim and Joe | Three the ray trace cannot settle
- **The entrance beam.**  Joe's 183 mm at f 330 mm is F/1.803, so the telescope's cone along the slit is F/1.89 and cannot be "a little faster than the spectrometer" without a larger beam; the paper's is 192 mm.  Is the aperture Joe's number or the paper's?
- **The slit and SRF.**  With a 2-pixel slit the SRF floor is 2.01–2.02 px (36 µm) with the pixel and diffraction, where this design sits to 0.03 µm; the paper's bound is 64.8 µm.  Is the 1.5–2.0 px range read against a 2-pixel slit, or does the project want a narrower slit at the throughput's expense?
- **The spectrometer's form.**  Our Dyson of record is the paper's Option A and meets the specification; its Option B (the conic grating) buys energy in a pixel 0.82 → 0.99 here but not alignment insensitivity.  What drove B in the build — fabrication, thermal, or the margin?
~ Still open from the last round: the detector window's gap (now answered: an order-sorting filter, not bonded) and the AR by band.

## Next steps | The spectrometer's tolerances and form; the telescope's tolerances and sizing; the record
::: left
- **Spectrometer:** the every-CTE thermal soak row; fabrication terms (figure, the grating's period and groove errors) in the ladder; a vector grating efficiency to replace the scalar blaze; Option B as the candidate, its keystone tolerances (25 µm, 33 µrad) stated beside its margin.
- **Telescope:** the same ladder on R9 and the 1.5k R5 (mirror decenter, tilt, figure; the slit's position); the mirrors sized with mounts and baffles (the blocks on slides 6 and 18 carry 10 mm blanks); the zig-zag's full-aperture M2 and M3 against the three-asphere 1.5k's smaller mirrors, if the four-module form returns.
::: right
- **The record:** the template `tma_longslit` is general (one parameter file: f, D, F/#, strip, slit, the cone bound, working distance, the seed legs); its default run is this module, `tma_longslit_1k5` the other.  A note to Jim and Joe with this deck, superseding "the 3k module is not there yet."
- **In the engine, from this round:** the api's element stop now takes the CLI's range (it refused the secondary of a four-element telescope); both launch conventions in the join scorer; the spectrometer's chief-tilt diagnostic.
~ `BRIEF_to_dyson5.md` addenda 48–49.

## Run it yourself | One parameter file, one runner, every stage through it; the telescope is a template any long-slit spectrometer can use
::: full
```
cd <path>/MACOS_resources/mmacos                   % the mmacos folder of the repository (the engine is built already)
matlab
>> run('mmacos_setup.m');
>> addpath('templates/10_telescopes/tma_longslit');  addpath('challenges/dyson5');
>> P = tma_longslit_params;                        % f, D, F/#, strip, slit, pixel, the telecentric and cone bounds, the seed legs
>> OUT = tma_longslit_run({'first_order','section','figure','score','e2e'}, P);   % the ladder to the deck of record, then the join
>> OUT = dyson5_run(struct('stages', {{'t5f'}}, 'tel5f_deck', 'dyson5_cprime_3k.in'));   % the 3k module end to end, both launches
>> OUT = tma_longslit_1k5;                          % the same ladder at the 1.5k spec (strip 4.7 deg, 27 mm slit), joined to the 130 mm Dyson
```
- **Changeable parameters live in two files:** `tma_longslit_params.m` for the telescope (spec, seed, rows and weights, field set), `dyson5_params.m` for the Dyson (F-number, pixels, band, slit, block radius, glass, order, the join).  Every record in this deck is a stage output under `challenges/dyson5/` or the template's folder.
- **Checks:** `./run_mmacos_tests.sh tTmaLongslit tSpectrometerRx tTelescopeRx tStopReload` and the grating, glass, propagation and root-pick classes, all in the fast suite (566 pass, 0 fail, 2026-10-07).  `tTmaLongslit` re-scores the deck of record by name.

## Backup

## The first order by the numbers | The seed's legs at f 330 mm, the powers the telecentric closure gives, the beams, and where the alternatives fail
::: left
| item | value |
|---|---|
| legs M1→M2, M2→M3, M3→slit | 258.3, 245.8, 299.4 mm (paper: 270, 257, 313 at f 345) |
| chief incidence M1, M2, M3 | 30°, 37°, 15° |
| focal lengths M1, M2, M3 | +1012, −452, +246 mm (M3's = the M2→M3 leg) |
| entrance pupil | virtual, 347 mm behind M1 |
| beams at M1, M2, M3 | 183, 136, 166 mm |
| intermediate focus | none (the marginal ray never crosses the axis) |
| chief to the slit normal, section as emitted | ≤ 0.093° over ±4.7° |
| plate scale as emitted | 329.99 mm local, 329.87 at the edge |
::: right
- **Stop at M1, or between M1 and M2:** also telecentric, but it flips the form and M2 grows to 310 / 220–250 mm.
- **Tilted spheres:** cannot meet the first order in both sections (R_t / R_s = 1.33, 1.57, 1.07 needed); off-axis conic sections can — paraboloids when the parent axis is at the incidence angle.
- **The paper's drawing left one thing open** — whether the beam passes through a focus between M2 and M3 — and the numbers settled it: it does not.
~ `tls_first_order.txt`; template README, "The two facts the first order forces."

## The rungs per field | What each rung's telescope looked like along the slit, and what the join made of it
::: full
| rung | DOFs | spot FWHM along the slit, fields 0 → 4.7° (px) | chief − centroid across the slit at 4.7° (µm) | e2e CRF / smile (point-source) |
|---|---|---|---|---|
| R4 | conics, θ, focus, h⁴ h⁶ | 1–2 center, 5–7 outer | 9.1 | 2.5 / 0.73 |
| R5 | + pole-frame freeform 3–6, centroid-bow rows | 1 center, 5–7 outer | 1.1 | 7.0 / 0.09 |
| R6 | R5, +3000 evaluations | as R5 | 1.5 | 7.2 / 0.15 |
| R7 | from R4: along-slit ×3, outer fields ×2–3 | 1.02 everywhere | 2.43 | 1.18 / 0.14 |
| R8 | R7 + telecentric ×30 | 1.02 | 2.4 | 1.18 / 0.14 |
| R9 | R7 + offset rows ×1000 | 1.02 | 0.03 | 1.17 / 0.016 |
- The four solves thrown out (unbounded layout, free R_t / R_s, closure variants, a resume defect) are in the report; none changed the form.
~ `REPORT_dyson5_cprime.md`; `tma_longslit/runs_*.log`.

## The chief-tilt diagnostic | How the smile's cause was isolated without a solve
::: full
- **The hypothesis:** the telescope's chief rays reach the slit up to 0.51 mrad from its normal in the dispersion direction, even in field and roughly quadratic — smile's shape — and a Dyson's spectral centroid moves with the input chief angle.
- **Why it could not be tested by hand on the joined deck:** each field has two sky angles but needs three conditions (land on the slit line, two chief-angle components), and the join re-aims every launch through the grating, so a telecentric error appears there as a pupil shift.
- **The test:** the spectrometer alone, each slit point's input chief tilted by the telescope's measured angles, scored as the record is (slide 13's left table).  Smile 0.0050 → 0.0049; keystone moves, so the hook is live.  The hypothesis is dead; the chief − centroid offset model (right table) predicts all five rungs.
- **Why R8 said nothing:** its telecentric rows were 0.4 % of the merit.
~ `tls_dyson_chief_tilt.m` (wraps the launch; the library untouched), its `.mat`; the template README.

## The tolerance ladder by the numbers | One perturbation at a time on the 3k Dyson of record, scored alone and through the R9 join; the compensator is the detector's focus
::: full
| perturbation (unit) | CRF, uncompensated (px) | CRF after detector refocus | keystone (px) | the tolerance for the paper's keystone CBE / for Joe's 0.1 px |
|---|---|---|---|---|
| lens decenter along the slit, 10 µm | – | – | 0.0099 | 26 µm / 47 µm |
| grating decenter along the slit, 10 µm | – | – | 0.0099 | 26 µm / 47 µm |
| grating tilt about the slit, 10 µrad | – | – | 0.0074 | 35 µrad / 64 µrad |
| grating clocking, 10 µrad | – | – | 0.0029 | 88 µrad / 162 µrad |
| grating radius dR/R 10⁻⁴ | +1.93 | +0.006 | – | focus-compensated |
| lens radius, thickness; air space; slit, detector defocus | focus terms | ≤ 0.006 | – | focus-compensated |
| CaF2 index, 10⁻⁵ (1 K) | +0.0014 | – | – | negligible |
| CaF2 radius + thickness per K (CTE) | +0.050 − 0.051 | – | – | cancel |
| aluminum air space per K; N-BK7 grating radius per K | −0.22; +0.14 | focus terms | – | refocus-compensated |
- **Keystone is lateral color: no compensator touches it.**  Smile is loose everywhere (≥ 300 µm or µrad).  Linearity at 2×: keystone vs decenter 1.10, CRF vs grating radius 1.03, SRF 2.6 (quadratic, and at the floor).
- **The 1.5k Dyson** has the same pattern at twice the keystone sensitivity (1.82 × 10⁻³ px per µm: 13 µm for the paper's CBE, 25 µm for Joe).  **Option B** raises the keystone sensitivities 6–45 % (clocking most): 25 µm, 33 µrad, 64 µrad.
~ `REPORT_dyson5_spec.md`; 21 rows; the paper's Option A as-built keystone 5.7 % of a pixel = 1.03 µm (Table 2).  Compensator residuals for focus alone and focus + x/y are both in the report (x/y adds nothing: the metrics are relative).

## An engine finding from this round | The api refused an element stop on the secondary of a four-element telescope
::: full
- **`stop_info_set` refused `iElt ≥ nElt−2` and `nElt ≤ 3`** — the original wrapper's range, with no engine reason — so `macos.stop(2)` failed on every M1-M2-M3-FP deck while the prescription keyword `ApStop=` and the CLI's STOP (1..nElt) accepted it.  The template worked around it with one deck per field until the fix.  Now the CLI's range; gate `tStopReload/test_stop_accepts_the_secondary_of_a_four_element_deck`.
- **Also from this round, in the tools:** both launch conventions in the join scorer (`dyson5_t5f`, `tel5f_launch`), the per-field chief-angle components in every rung table, the spectrometer's chief-tilt diagnostic.
~ macos `a09cd59`; resources `9845c2d`.

## Records behind each slide | Every number traces to a runner stage and its record file
::: full
| slide | record |
|---|---|
| 3, 4, 5 | `BRIEF_to_dyson5.md` addenda 47–49; `templates/10_telescopes/tma_longslit/` (README, `tma_longslit_params.m`, `tma_longslit_run.m`); `REPORT_dyson5_cprime.md` |
| 6, 7 | `dyson5_block2_3k.txt`, `_views.png`, `_swath.png` (tool `dyson5_block4.m`, nmod 2) |
| 8 | `tma_longslit/tls_first_order.txt`, `tls_section_layout.png`; `dyson5_cprime_3k_viewyz.png` |
| 9 | `challenges/dyson5/dyson5_cprime_3k.in`; `REPORT_dyson5_cprime.md` (R9); `tTmaLongslit` |
| 10, 11 | `dyson5_jim_3a.txt`; `dyson5_t5f_cprime_centroid_roll000.txt`, `_roll180.txt`, `_chief_roll000.txt`; `_e2e.in`, `_maps.png`, `_e2e_viewyz.png` |
| 12, 13 | `REPORT_dyson5_cprime.md` (R4–R9); `tls_dyson_chief_tilt.mat`; `tls_e2e_both_conventions.mat` |
| 14 | R9 per-field table (`REPORT_dyson5_cprime.md`) |
| 15 | `dyson5_jim_3a.txt`; `REPORT_dyson5_spec.md`; `spectrometer_sens`, `tSpectrometerSens`; `BRIEF_to_dyson5.md` addendum 49 |
| 16 | `dyson5_throughput_band.txt` (tool `dyson5_throughput_band.m`); `dyson5_jim_3b.txt` route 2 |
| 17 | `tma_longslit_1k5.m`; `tls1k5_figure.txt`, `tls1k5_e2e.txt`, `tls1k5_R5_ffc.in` = `dyson5_cprime_1k5.in`; `dyson5_rescore_1k5_launch.txt` |
| 18 | `dyson5_jim_3a.txt`, `dyson5_block2_3k.txt`, `dyson5_block4_cprime1k5.txt`; the superseded three-asphere 1.5k: `dyson5_rescore_1k5_launch.txt`, `dyson5_block4.txt` |
~ All under `MACOS_resources/mmacos/challenges/dyson5/` and `templates/10_telescopes/tma_longslit/`, branch dev-candidate; the reports under `macos/`.  Public: nasa-jpl/macos, nasa-jpl/MACOS_resources.
