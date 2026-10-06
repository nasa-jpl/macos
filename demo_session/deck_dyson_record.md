<!--
deck_dyson_record.md — dyson5: the instrument of record (telescope + Dyson) at Jim's numbers,
the short deck for Jim and Joe (2026-10-06).  The full working record stays in deck_dyson.md.
Build: python3 make_brief_slides.py deck_dyson_record.md
Figures: the runner's own PNGs, copied unmodified into figs_dyson/; the *_rc.png copies are
recomposed at the panel level only (recompose.py), the *_views.png are the engine's two
renders tiled (tile_views.py).  Producers: dyson5_t5f.m (maps), dyson5_view_figs.m (views),
spectrometer_maps_fig.m, dyson5_size_fig.m, dyson5_trade.m.
Sources: challenges/dyson5/README.md, BRIEF_dyson5_jim.md (3a), BRIEF_dyson5_size.md,
BRIEF_dyson5_tma.md, BRIEF_to_dyson5.md addenda 42–45, records dyson5_t5f_GM_*.txt.
Style: doc/DECK_STYLE.md + doc/STYLE_REPORTS.md (lean main path, one Backup divider,
layout + map per result slide, figures unmodified, conventions stated once).
-->

# A VSWIR Dyson Spectrometer and Its Telescope, Designed From the Spec
An EMIT-class push-broom spectrometer (F/1.8, 18 µm pixels, 380–2500 nm) at Jim's ground sample (30 m from 550 km), designed with MACOS from the public specification and scored on smile, keystone and the response functions, telescope and spectrometer together.
D. C. Redding, with Claude Code.
October 2026.  Working record; no proprietary prescription is used.
~ DRAFT for review, 2026-10-06.  The 1.5k module (1500 cross-track pixels, a 27 mm slit) has a telescope and a spectrometer of record that image at the pixel end to end: smile 0.51, keystone 0.01, CRF 1.22, SRF 2.20 pixels, every clearance positive; four such modules carry 6000 pixels in 6.6 kg of fused silica.  The 3k module's telescope is not yet at the pixel (CRF 4.2 pixels), so the two-module CaF2 option waits on it.

## Contents | The target, the instrument of record, the comparison for Jim, the open module, the questions
::: left
- **3** The target and how it is scored
- **4** The 1.5k module, end to end
- **5** The telescope of record
- **6** The spectrometer of record for the 1.5k module
- **7** Two of 3k in CaF2 or four of 1.5k in silica
- **8** The 3k module: where it stands
::: right
- **9** Questions for Jim
- **10** Next steps
- **11** Run it yourself
- **12** Backup: how the telescope reached the pixel; the Dyson design ladder; the 54 mm design; the block-size trade; the methods and the checks; twelve engine findings; records

## The target and how it is scored | Joe's specification at Jim's ground sample; the metrics stated once, here
::: left
| item | value |
|---|---|
| F-number (air-equivalent, image space) | 1.8 |
| pixels | 18 µm; 3000 or 2 × 1500 cross-track, 500 spectral |
| slit | 54 mm (3k module) or 27 mm (1.5k module); 2 pixels wide |
| band | 380–2500 nm, 4.24 nm per pixel |
| smile, keystone | < 0.1 pixel |
| SRF FWHM; CRF FWHM | < 1.5–2.0 pixels; < 1.5 pixels |
| ground sample (Jim) | 30 m from 550 km: f = 330 mm, 183 mm aperture at F/1.8 |
| cross-track field | 9.4° per 3k module, 4.7° per 1.5k module |
::: right
- **Smile:** centroid drift along the slit at fixed wavelength; **keystone:** drift across wavelength at fixed slit position; both in pixels from ray centroids.
- **SRF, CRF:** slit ⊗ line-spread function ⊗ pixel ⊗ Airy, FWHM in pixels (Mouroulis & Green 2018); the Airy diameter at 2500 nm is 11 µm, under the pixel, so the incoherent chain is valid.
- **Energy in one pixel:** the fraction of a point image inside one 18 µm pixel; the literature's design rule is 0.75.
- **A module** is one telescope, one spectrometer and one detector; **end to end** is one prescription from the sky to the detector, the grating its stop, scored by the spectrometer's own scorer.
- **As-built realism (Jim):** SRF 2.5–3 pixels, slits 2 pixels wide, photon-limited systems; smile and keystone drive the design.
~ Every number is a MACOS ray trace of a committed prescription, checked ray by ray against an independent solver (Backup); no coatings, no tolerances.  Public comparison point: Carbon-I (arXiv:2505.22545), F/2.2.

## The 1.5k module, end to end | Three off-axis aspheres feeding a 130 mm silica Dyson: smile 0.51, keystone 0.01, CRF 1.22, SRF 2.20 px over the whole strip; every field admitted; clearances +39 mm inside the telescope, +0.38 mm overall
::: left
![As the engine traces it (meters): the telescope's M1–M3 (E1–E3) at the lower right, the slit (E4) on the block's face, the 130 mm block and the grating (E7) above it, the detector (E11) beside the slit.  Left: 3-D; right: the dispersion plane.](figs_dyson/dyson5_t5f_GM_1k5_bAs_e2e_views.png){h=4.5}
::: right
![Per slit position and wavelength (7 × 7): field angle, keystone, smile, the two response widths and the energy in one pixel.  The response widths rise only at the strip ends.](figs_dyson/dyson5_t5f_GM_1k5_bAs_maps_rc.png){h=4.5}
~ Energy in one pixel end to end: 0.78–0.89 over the central third of the strip, 0.65 at ±9 mm, 0.30–0.48 at the strip ends (±13.5 mm); the spectrometer alone puts 1.00 in the pixel, the rest is the telescope's 8–14 µm blur.  Plate scale 327 mm at the center, 331 at the edge (330 specified).  Record `dyson5_t5f_GM_1k5_bAs.txt`; prescription `dyson5_t5f_GM_1k5_bAs_e2e.in`, 11 elements.

## The telescope of record | f = 330 mm, 183 mm aperture, F/1.8: an unobscured section of a three-mirror parent, M3 ahead of the intermediate focus so the exit pupil sits at infinity, as the Dyson requires
::: full
| mirror | R (mm) | K | h⁴, h⁶ coefficients (m⁻³, m⁻⁵) | departure from the best-fit sphere over the lit patch (mm) | vertex z (mm) |
|---|---|---|---|---|---|
| M1 (concave) | 395.6 | −1.672 | +0.944, −1.010 | 2.1 | −5.8 |
| M2 (convex) | 109.8 | −2.836 | +70.6, +1300 | 0.3 | −152.4 |
| M3 (concave) | 159.3 | −0.517 | −1.66, +2306 | 0.01 | −59.5 |
| focal plane | flat | – | – | – | −125.7 |
::: left
- **Geometry:** the beam passes 190 mm from the parent axis, the strip is biased 4° and the Dyson is rolled 180° about the chief ray; from the first-order layout M2 is tilted 0.95° and moved 6 mm, M3 tilted 6.6° and moved 24 mm, the focal plane 8.5 mm along the axis.  No hole, no central obstruction.
- **Image as placed:** 8 µm rms at the strip center, 14 µm at the ±2.3° edge (0.4 / 0.8 pixel); every field admitted; the strip's chief rays within 0.5° of each other at the slit.
::: right
![The telescope alone as the engine traces it (meters): light enters from above; M1 is the large mirror at the right, M2 the small one at the lower left, M3 above M2, the focal plane (E4) below M1.](figs_dyson/dyson5_tA_GM_1k5_bAs_views.png){h=3.4}
~ Radii are |R|; vertex z along the parent axis in the prescription's frame; the h⁴/h⁶ terms are about each mirror's own vertex and axis.  Prescription `dyson5_tA_GM_1k5_bAs.in`; solve record `dyson5_tA_GM_1k5_bAs.txt`.  How it was reached: Backup.

## The spectrometer of record for the 1.5k module | A 130 mm fused-silica block with no meniscus: CRF 1.03 px, all the energy in one pixel, 1.6 kg, four air-glass crossings
::: left
- **The block follows the slit.**  In the Dyson form the flat face sits at the center of curvature, so the block's thickness is its radius, and the corner blur grows as h⁴/r³ with the slit's half-length h.  Halving the slit from 54 to 27 mm lets the block fall from 220 mm to 130 mm and drops the meniscus.
- **The 27 mm module, at 130 mm:** CRF 1.03 px, SRF 2.02 px, smile 0.007, keystone 0.009 px, 1.00 of the energy in one pixel; grating 171 mm across on a 0.41 m radius; slit plane to grating 414 mm; 0.75 L, 1.6 kg edged; uncoated throughput 0.87 over four crossings.  It tolerates 1.5 mm of air at the slit and detector for free (CRF 1.04 px).
- **The 54 mm single module, for comparison:** silica needs the 220 mm block and a 4 mm meniscus (CRF 1.33 px, 0.76 in one pixel, eight crossings, 8.2 kg); CaF2 at 240 mm needs none (CRF 1.21 px, 0.82, 14.3 kg).
::: right
![Per block radius: CRF (left), energy in one pixel (middle), mass of the edged block (right); solid lines with the meniscus, dashed without; a filled marker meets the specification, a black ring also matches the 220 mm design's own scores.  The 27 mm slit (blue) holds the pixel down to 130 mm in silica without a meniscus.](figs_dyson/dyson5_size_rc.png)
~ Record `dyson5_size.txt` (48 engine-scored points, each a full re-solve); the working-distance scan and the throughput routes in `dyson5_jim_3b.txt`.  The 54 mm designs and their ladder: Backup.

## Two of 3k in CaF2 or four of 1.5k in silica | On the spectrometer side four small silica modules win on glass; with the telescope included, the 1.5k module images at the pixel end to end and the 3k module does not yet
::: full
| | 2 × CaF2, 3k, 54 mm slit | 4 × fused silica, 1.5k, 27 mm slit |
|---|---|---|
| block radius = thickness | 240 mm | 130 mm |
| grating diameter; slit plane to grating | 318 mm; 787 mm | 171 mm; 414 mm |
| spectrometer alone: CRF; energy in one pixel | 1.21 px; 0.82 | 1.03 px; 1.00 |
| spectrometer alone: smile / keystone | 0.005 / 0.006 px | 0.007 / 0.009 px |
| air-glass crossings; uncoated throughput | 4; 0.88 | 4; 0.87 |
| glass per block, edged | 4.5 L, 14.3 kg | 0.75 L, 1.6 kg |
| glass for 6000 pixels | 9.0 L, 28.5 kg | 3.0 L, 6.6 kg |
| single-crystal CaF2 to carve from | 6.2 L per block, 12.4 L in all | none: fused silica is a melt |
| telescopes / gratings / detectors | 2 / 2 / 2 × 3k | 4 / 4 / 4 × 1.5k |
| telescope of record | freeform three-mirror, not at the pixel (open) | three off-axis aspheres (slide 5) |
| end to end: smile / keystone / CRF / SRF | 2.05 / 0.05 / 4.22 / 4.90 px | 0.51 / 0.01 / 1.22 / 2.20 px |
| end to end: energy in one pixel | 0.05–0.09 | 0.30–0.89 (0.8 over the central third) |
| clearance, worst pair | +0.06 mm (M2→M3 beam against the block face) | +0.38 mm (detector package against the block face) |
~ Carve = a cylinder of clear aperture + 20 mm by thickness + 20 mm; CaF2 is not priced here.  A cold detector adds a dewar window to every row alike; cementing it to the block keeps the count at four.  The 3k telescope is the freeform solve of 2026-10-05, scored before its next layout step (slide 8).  Records `dyson5_jim_3a.txt`, `dyson5_t5f_GM_1k5_bAs.txt`, `dyson5_t5f_GM_3k_c.txt`.

## The 3k module: where it stands | The same parent with a 9.4° strip needs freeform mirrors and is not yet at the pixel: CRF 2.3–4.2 px, 5–9 % of the energy in one pixel, clearance +0.06 mm
::: left
![The 3k telescope through the 240 mm CaF2 Dyson as the engine traces it (meters): the same arrangement as the 1.5k module at twice the strip.](figs_dyson/dyson5_t5f_GM_3k_c_e2e_views.png){h=3.0}
- **What the solve gave:** Zernike freeform terms on all three mirrors (M2 carries 1.7 mm of departure at 126 mrad of slope) with the same merit and geometry freedoms as the 1.5k: smile 2.05, keystone 0.05, CRF 4.22, SRF 4.90 px; 98 % of the light admitted at the strip ends; the beam between M2 and M3 passes the block's face by 0.06 mm.
- **Why it stalls:** on this section the figure terms and the mirror positions work against each other (the freeform norm grows as the geometry moves), which points at the stop's position.  **Next:** the stop at M2 instead of M1; an asphere-only solve from the same seed; the field weighting across the strip as a project ruling.
::: right
![Per slit position and wavelength: the CRF is 2.3–3.7 px off center and 4.0–4.2 px at the strip center, the signature of a solve that traded the center for the ends under equal field weights.](figs_dyson/dyson5_t5f_GM_3k_c_maps_rc.png){h=4.3}
~ Record `dyson5_t5f_GM_3k_c.txt`; prescription `dyson5_t5f_GM_3k_c_e2e.in`.  The 3k solve ran to the same evaluation budget as the 1.5k and is not converged in the optimizer's sense.

## Questions for Jim | Two decisions the ray trace cannot make
- **Two modules or four.**  The 1.5k module now has a telescope of record at the pixel end to end (slides 4–5); the 3k module does not yet (slide 8).  At the same 6000 cross-track pixels, four telescopes, spectrometers and detectors in 6.6 kg of fused silica against two in 28.5 kg of CaF2: is four a packaging the project would entertain?
- **The CaF2 volume.**  Two 240 mm blocks need 12.4 L of single-crystal CaF2 to carve from (6.2 L each, clear aperture + 20 mm each way; 19.6 L at the 300 mm blocks that put all the energy in a pixel).  What does that cost against two to four fused-silica blocks of 0.75 L?
~ Also open from the last round: whether the detector window is cemented to the block in the builds (six crossings and 0.81 uncoated, against four and 0.87).

## Next steps | The 3k telescope, the cold window, the radiometric chain
::: left
- **The 3k telescope:** the stop at M2; an asphere-only solve from the freeform seed; a merit written on the response-function width for the final step; the field weighting across the strip as a ruling.
- **The 1.5k module's package:** the cold window cemented as the block's exit face with the shield as the detector housing behind it, under the same clearance check; the 1.5 mm of working distance the 130 mm module tolerates is where it goes.
- **The radiometric chain against the band:** throughput, grating efficiency, detector quantum efficiency and the slit's diffraction loss (0.3 % at 380 nm to 1.8 % at 2500 nm past the grating, before the telescope's cone).
::: right
- **In the engine:** a fixed-station reference sphere for the wavefront merit on telecentric designs; an asphere-plus-Zernike surface in one declaration; a short-line guard for every multi-value prescription keyword.
- **Later, once the design of record is stable:** full diffraction from the slit through the grating to the detector (the slit's truncation of the telescope's image, about 10 % on the response functions in the literature); polarization sensitivity of the Dyson against the Offner; a surface-by-surface tour of the prescription.

## Run it yourself | One parameter file, one runner, every stage through it
::: full
```
cd <path>/MACOS_resources/mmacos                   % the mmacos folder of the repository (the engine is built already)
matlab                                            % start MATLAB here
>> run('mmacos_setup.m');                          % puts the engine and the design tools on the path
>> addpath('challenges/dyson5');
>> OUT = dyson5_run(struct('stages', {{'s1','s2','s3'}}));   % emit the Dyson decks, score them, run the design ladder
>> OUT = dyson5_run(struct('stages', {{'t5f'}}));           % the telescope + Dyson end to end from a telescope deck
>> OUT = dyson5_run(struct('Fno', 2.0, 'block_r_m', 0.25, 'stages', {{'s1','s2','s3'}}));  % your own instance
```
- **Changeable design parameters live in one file, `dyson5_params.m`:** F-number, pixel count and pitch, band, slit width, block radius, glass, diffraction order, the ladder's steps and weights, the aperture and mount margins, the detector package, the telescope deck and the Dyson it feeds.  Every record in this deck is a stage output under `challenges/dyson5/`.
- **Checks:** `./run_mmacos_tests.sh tSpectrometerRx tTelescopeRx` and the grating, glass, propagation and root-pick classes, all in the fast suite (536 pass, 0 fail, 2026-10-05).

## Backup
How the telescope reached the pixel, the Dyson design ladder, the superseded 54 mm designs, the methods and checks, the engine findings, the records.

## How the telescope reached the pixel at the strip edge | The defect was cross-track astigmatism growing as the square of the field; figure terms alone flattened the strip by trade, the merit and the geometry removed it
::: left
- **The pupil is a first-order condition.**  The Dyson accepts chief rays within 0.09° of its own, so the telescope's exit pupil must sit near infinity.  A three-mirror with a real intermediate focus between M2 and M3 cannot be telecentric; with M3 ahead of that focus the chiefs sit within 0.5° across the strip and every field is admitted.
- **It was not focus or the mirrors' power split:** along the strip the along-track focus is flat, the cross-track focus curves 2–4 mm and the two line foci at the edge sit 2–5 mm apart; no secondary magnification flattens both, and a field lens has no leverage at the slit.
- **Figure terms alone trade the center for the edge:** Zernike departures on all three mirrors took the 1.5k edge 291 → 170 µm while the center went 29 → 73 µm.  A control with the wavefront merit written per ray gave the same trade: the merit, not the surfaces, set it.
- **Two levers removed it, neither a figure term:** the merit written on the per-ray transverse spot at the slit with the per-field image positions (at F/1.8 with tens of µm of path error the wavefront rms flattens the field by trading the center), and M2/M3 tilt, decenter and the focus as variables.  With both, the freeform terms collapsed to a re-fitted set of off-axis aspheres, re-solved as aspheres: the record.
::: right
| step (1.5k, same seed, same budget) | as-placed spot, center / edge | end to end: smile / keystone / CRF / SRF (px) |
|---|---|---|
| conics + aspheres, wavefront merit | 37 / 343 µm | 3.26 / 0.63 / 15.5 / 5.9 |
| + Zernike freeform, wavefront merit | 99 / 219 | 0.82 / 0.33 / 7.7 / 9.2 |
| same freedoms, transverse (spot) merit | 18 / 27 | 1.54 / 0.03 / 2.34 / 4.46 |
| + M2/M3 tilt, decenter, focus | 6 / 13 | 0.60 / 0.02 / 1.49 / 2.30 |
| **three off-axis aspheres, same merit and geometry** | **8 / 14** | **0.51 / 0.01 / 1.22 / 2.20** |
~ Records `dyson5_tA_EP_*`, `dyson5_tA_pz_*`, `dyson5_tA_FF_*`, `dyson5_tA_GM_*`; the end-to-end rows through the engine-placed join.  Solves run to a fixed evaluation budget; the emission is re-solved, not re-fitted (a 1.7 mrad slope residual is a millimeter of blur over 0.3 m).  The coaxial three-mirror and the two-mirror modified Schwarzschild were tried first and do not image a 4.7° strip at F/1.8 packaged: `deck_dyson.md` slide 27.

## The Dyson design ladder at the 54 mm slit, by the numbers | Keystone falls 37-fold at a fixed 220 mm block; the corner blur is bought only by the meniscus or by size
::: full
| step | what moves | keystone (px) | smile (px) | CRF FWHM (px) | energy in 1 px | length (mm) | grating footprint (mm) | glass (L) |
|---|---|---|---|---|---|---|---|---|
| R0 | concentric seed | 0.095 | 0.006 | 2.28 | 0.44 | 708 | 275 | 22.2 |
| R1 | grating-radius factor, face offset | 0.038 | 0.003 | 2.10 | 0.48 | 704 | 274 | 22.3 |
| R2 | + conic and h⁴, h⁶ terms on the block face | 0.037 | 0.004 | 2.13 | 0.47 | 704 | 274 | 22.3 |
| R3 | + block center off the grating's, along the dispersion | 0.011 | 0.008 | 2.10 | 0.48 | 698 | 272 | 22.1 |
| R4a | meniscus corrector alone | 0.063 | 0.012 | 2.49 | 0.33 | 698 | 270 | 22.3 |
| **R4** | **meniscus with all variables (the 54 mm design of record)** | **0.0026** | 0.005 | **1.33** | **0.76** | 694 | 269 | 22.7 |
| R5 | + slit plate and fold prism: the detector out of the slit's plane | 0.0085 | 0.005 | 1.27 | 0.70 | 694 | 269 | 19.2 |
| size alone | block radius free (341 mm) | 0.024 | 0.001 | 1.34 | 0.74 | 1099 | 427 | 82.8 |
~ Each step solved on the exact chain with pixel-unit smile and keystone in the merit, then emitted and scored in the engine.  SRF FWHM is 2.02–2.05 px on every step: the 2-pixel-slit floor.  The concentric condition R_g = n·r/(n − 1) is the blur minimum by exact trace; a single concentric block meets a quarter-pixel corner blur only at r ≥ 213 mm, which is why every flight form adds an asphere, a meniscus or an air gap.  Record `dyson5_s3_trade.txt`.

## The 54 mm design of record: the meniscus-corrected Dyson (step R4) | Keystone 0.0026 px, smile 0.005 px, CRF 1.33 px, 0.76 of the energy in one pixel at a 220 mm silica block; superseded for the 1.5k module by the plain 130 mm block
::: left
- **Eleven surfaces:** the block's flat and spherical faces, a 4 mm meniscus in the air gap, the grating, and the same faces on the way back; block radius 220 mm, grating 269 mm across on a 0.69 m radius, slit plane to grating 694 mm.
- **Decentering the block solves the distortion; the meniscus buys the blur.**  Keystone 0.095 → 0.011 px comes from the block's center moving off the grating's along the dispersion; the corner blur does not move on any step without the meniscus or a 341 mm block.  A global search over the meniscus (12 starts) found no better basin on the detector.
- **Checked by propagation:** the point-spread function's centroid equals the ray centroid to 0.002 px at every slit position and wavelength.
- **A cold shield cannot be bought in this form:** every millimeter of air before the detector costs about 0.5 px of CRF; a fold prism takes the detector out of the slit's plane at CRF 1.27 px and 0.70 in one pixel.
::: right
![Engine-scored maps for R4 on 18 µm pixels, per slit position and wavelength; each convention on its colorbar, maxima in the labels.](figs_dyson/dyson5_s3_maps_r4_rc.png){h=4.4}
~ Prescription `dyson5_s3_r4.in`; records `dyson5_s3.txt`, `dyson5_s3_r4global.txt`, `dyson5_s2w.txt`, `dyson5_s5.txt` (the fold prism and the cold-shield sweep).  The full ladder, the fold prism and the closure envelope: `deck_dyson.md`.

## The block-size trade, by the numbers | Two modules with no meniscus: a 130 mm silica block at 1.6 kg images at the pixel with half the air-glass crossings of the 54 mm design
::: full
| configuration | block radius and thickness (mm) | edged mass (kg) | air-glass crossings, uncoated throughput | CRF (px) | energy in 1 px | smile / keystone (px) |
|---|---|---|---|---|---|---|
| silica, one 54 mm slit, meniscus (R4) | 220 | 8.2 | 8, 0.76 | 1.33 | 0.76 | 0.005 / 0.003 |
| CaF2, one 54 mm slit, no meniscus | 240 | 14.3 | 4, 0.88 | 1.21 | 0.82 | 0.005 / 0.006 |
| **silica, two 27 mm slits, no meniscus** | **130** | **1.6 each** | **4, 0.87** | **1.03** | **1.00** | 0.007 / 0.009 |
| silica, two 27 mm slits, smallest matching R4 | 100 | 0.9 each | 4, 0.87 | 1.25 | 0.80 | 0.008 / 0.012 |
| CaF2, two 27 mm slits, smallest matching R4 | 80 | 0.9 each | 4, 0.88 | 1.02 | 1.00 | 0.002 / 0.007 |
- **The meniscus was paying for the long slit.**  With one 54 mm slit, silica closes only with the meniscus; without one, a single module closes only in CaF2 at 240 mm.  With two modules the plain block matches R4 from 220 mm down to 100 mm in silica and 80 mm in CaF2, with four crossings where R4 has eight and no thin part to make, mount or shake.  Halving the glass path also halves the wavefront error from index inhomogeneity and thermal gradients (0.4 waves each at 220 mm).
- **Throughput routes measured on the 130 mm module:** cementing the dewar window to the block, 6 → 4 crossings, 0.81 → 0.87; a two-layer broadband anti-reflection coating, +4 % (a 6:1 band is wide for a simple stack); the two together about 0.91.
~ Record `dyson5_size.txt` (CCMac, 2026-10-02): 48 engine-scored points; throughput is the uncoated Fresnel product at 1 µm, normal incidence.  Routes: `dyson5_jim_3b.txt`.

## The methods, and the checks behind every number | Ray tracing scores the design, wave propagation checks it, and two independent traces agree ray by ray
::: left
- **Ray tracing (MACOS):** rays leave a point on the slit at one wavelength, refract through the block, diffract at the grating into the design order and land on the detector; the centroid per slit position and wavelength gives smile and keystone, the spread convolved with the slit and the pixel gives the response functions.  End to end, the rays leave the sky and the grating is the stop.
- **Wave propagation:** the same prescription carried as a complex wavefront to the detector; its point-spread function's centroid equals the ray centroid to 0.002 px at every slit position and wavelength, so a detector sees the ray numbers.  Across a grating the comparison is made in phase (modulo λ).
- **Two programs, one prescription:** the design solver (MATLAB) lays out each form and writes the prescription; MACOS traces it on its own; every ray is re-traced by the solver and the two agree to 1e-9 m at every surface (1e-12 m with the aspheric block face, where a solver that ignores the asphere misses by 0.26 mm).  For the freeform telescope decks, which the solver cannot model, the end-to-end join is placed from engine traces and reproduces the solver's join on an aspheric deck to every printed digit.
::: right
- **The optimizer:** the design layer's native multi-field solve (MACOS CALIB) for the conics and aspheres; a Levenberg–Marquardt solve on the engine's own per-ray spot residuals, with Jacobian scaling, for the freeform and geometry steps.  Solves run to a fixed evaluation budget.
- **Clearance:** every design passes a check of each beam leg against every body (mirror blanks on their best-fit spheres, the block, the slit mask, the detector package with a 5 mm carrier); the number quoted is the worst pair's margin.
- **Slit diffraction loss:** a 36 µm slit at F/1.8 loses 0.3 % of the light past the grating at 380 nm and 1.8 % at 2500 nm, by the engine's far-field leg and by the closed sinc² form alike.
~ Regression tests: the spectrometer and telescope prescriptions, the grating, glass, propagation and root-pick classes, in the MACOS fast suite (536 pass, 0 fail, 2026-10-05).

## Eight engine findings from the spectrometer, each with a regression test | Found by the design's own checks between 2026-09-30 and 10-03, fixed in the engine the same day
::: full
| finding | symptom | fix |
|---|---|---|
| glass names ignored at load | a `GlassElt= Silica` element traced as air | the glass catalog is blanked once at allocation, not on every load |
| vacuum wavelength inside glass | a propagation leg in a medium ran at the wrong Fresnel number by n | every kernel takes λ/n of the leg's medium |
| groove period along the surface | 2.8 / 3.3 px rms spectral blur at 2500 nm (Offner / Dyson) | the local grating vector is the projection of the ruling direction, un-normalized (chord-ruled) |
| grating path-length jump | 4 waves rms of pupil error where the rays converged to 0.05 µm | the jump is the groove count from the vertex along the ruling direction |
| optimizer derivative stride | a multi-field spot-size solve overwrote memory on its second field | advance by the objective's size, as the value loop does |
| ray behind the surface after a pupil reference | a two-mirror with its pupil reference 5 mm ahead of a convex primary scattered every ray in a one-call trace | the conic root nearest the element's reference point; a restarted trace resets its previous element per ray |
| beam conditions only under a beam target | chief direction and position could not be scored with a spot or wavefront target | the beam rows ride on any target; a position target per field; the centroid as the position |
| optimizer asphere step and failure exit | a fixed 1e-10 probe step: a zero derivative column, a singular matrix, and a `stop` that ended the host process | the step is 1e-3 of the coefficient (sag-based for a zero one); the failure restores the optics and returns a flag |
~ Tests: `tGlassDispersion`, `tPropMedium`, `tGratingImmersed` + `test_grating_chord`, `tGratingOpl`, `tSpectrometerRx`, `tTraceRestart`, `tBeamRows`, `tAsphCalib`.  Flat gratings are unchanged by the grating fixes; the wavefront-target optimizer, which every telescope design used, was never affected by the fifth.

## Four more findings from the telescope's eccentric section | Found 2026-10-04 to 10-05, each invisible on a coaxial design, each fixed the same day
::: full
| finding | symptom | fix |
|---|---|---|
| ray on the far sheet of a hyperboloid | 40 of 1185 rays on the section's M3 took the second sheet: 44 mm of extra path, every ray reported as passing | the root on the sheet the vertex is on; 406 prescriptions pre/post identical |
| optimizer aperture sized at one field | the strip-edge fields lost every ray on M3 and were dropped by the solve | the aperture encloses the footprint over every field of the solve |
| stop aimed at the parent vertex | on an eccentric section the pupil sphere was placed 190 mm from the beam | the stop at the prescription's object-space stop point |
| Zernike origin at the parent vertex | freeform terms on a section evaluated 3.8 radii outside the unit disc | the origin is the element's pole |
~ Tests: `tFwdRoot` (5), `tAsphHook` (3), `tStopApStop` (4), `tFreeformPole` (3).  Each of the four passed every existing check (rays passed, solves converged, decks loaded); what exposed them was a section whose beam is 190 mm from its parent's vertex.

## Records behind each slide | Every number traces to a runner stage and its record file
::: full
| slide | runner stage | record | tool |
|---|---|---|---|
| 1.5k module end to end; the 3k module | t5f on `dyson5_tA_GM_1k5_bAs.in` / `dyson5_tA_GM_3k_c.in` | `dyson5_t5f_GM_1k5_bAs.txt`, `dyson5_t5f_GM_3k_c.txt`, `*_e2e.in`, `*_maps.png` | `dyson5_t5f.m` (the engine-placed join), `spectrometer_score.m`, `spectrometer_clearance.m` |
| the telescope of record | `dyson5_tma_step5.m` (the geometry + merit steps) | `dyson5_tA_GM_1k5_bAs.txt`, `dyson5_tA_GM_conicfit.txt`, `BRIEF_dyson5_tma.md`, `BRIEF_to_dyson5.md` addenda 42–45 | `macos.design.tma_layout` 'telecentric', `Telescope.optimize`, the per-ray spot solver |
| the spectrometer of record; block-size trade | `dyson5_size_trade.m` | `dyson5_size.txt`, `dyson5_size.png`, `dyson5_size_*.in`, `BRIEF_dyson5_size.md` | continuation walks over `dyson_ladder.m` |
| comparison for Jim | `dyson5_jim.m` ('3a', '3b') | `dyson5_jim_3a.txt`, `dyson5_jim_3b.txt`, `BRIEF_dyson5_jim.md` | re-cut of `dyson5_size.mat`; `macos.design.thinfilm_rt` |
| design ladder, R4, trade | s3, s3free, s3b | `dyson5_s3.txt`, `dyson5_s3free.txt`, `dyson5_s3_trade.txt`, `dyson5_s3_layout_r*.png`, `dyson5_s3_maps_r*.png` | `dyson_ladder.m`, `dyson5_trade.m` |
| fold prism, cold shield | s5 | `dyson5_s5.txt`, `dyson5_s5_r5_h*.in` | `spectrometer_geom.m` (form `dyson_fold`) |
| chain vs engine; propagation twin; slit loss | gate; s2w; s2l | `tests/tSpectrometerRx.m`, `tTelescopeRx.m`; `dyson5_s2w.txt`; `dyson5_s2l.txt` | per-ray re-trace; `spectrometer_wave.m`; far-field leg vs sinc² |
| engine renders | `dyson5_view_figs.m` | `*_view3d.png`, `*_viewyz.png` | `macos.view_rx` on the prescription |
~ Conventions: block index n(silica, 1 µm) = 1.450417; radii quoted as |R|; the engine stores `KrElt` as −|R| and the center of curvature at the vertex plus |R| along the surface normal.  Reference document: `optical_design/SPECTROMETER_DESIGN_REFERENCE.md`.
