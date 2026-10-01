<!--
deck_dyson.md — dyson5: a VSWIR Dyson imaging spectrometer designed from the spec.
DRAFT — the deck grows with the challenge; every number is from a committed
runner stage in mmacos/challenges/dyson5 (records dyson5_s*.txt).
Build: python3 make_brief_slides.py deck_dyson.md
Figures: the runner's own PNGs, copied unmodified into figs_dyson/ (no crops yet).
Sources: challenges/dyson5/README.md + BRIEF_dyson5_beat{1,2,2c,3}.md (TO),
NOTE_mg2018_digest.md (CC), BRIEF_to_dyson5.md (the lane brief).
Style: doc/DECK_STYLE.md (lean main path, one Backup divider, layout + map per
result slide, figures unmodified, conventions stated once).
-->

# A VSWIR Dyson Imaging Spectrometer, Designed From the Spec
An EMIT-class push-broom spectrometer (F/1.8, 3000 × 500 pixels at 18 µm, 380–2500 nm) designed with MACOS from the public specification alone, scored on smile, keystone and the response functions, and checked by physical-optics propagation.
D. C. Redding, with Claude Code.
October 2026.  Working record; no proprietary prescription is used.
DRAFT — in progress.  Every number is from a committed, parameterized runner (mmacos, `dyson5_run`); records and figures are the runner's own output.  Design of record: the compact variant (step R4): ray-centroid keystone 0.0026 px, smile 0.005 px, CRF 1.33 px, 76 % of the energy in one pixel, at a 220 mm block.  Open: the propagated-PSF centroid differs from the ray centroid by up to 0.12 px on R4, the detector-seen keystone until it is resolved.

## Contents | The target, the forms, the design ladder, the evidence
::: left
- **3** The target and how it is scored
- **4** The Dyson as traced by the engine
- **5** The Offner as traced by the engine
- **6** The concentric Dyson and its scaling law
- **7** The exact chain and the engine agree
- **8** The seed, scored
- **9** The departure ladder
- **10** The departure ladder, by the numbers
::: right
- **11** What each departure buys
- **12** The compact variant, as traced by the engine
- **13** The compact variant in section
- **14** The compact variant, scored
- **15** The trade at a glance
- **16** The propagation twin
- **17** Slit diffraction loss
- **18** Next steps
- **19** Run it yourself
- **20** Backup

## The target and how it is scored | Joe's specification, EMIT-class; the metrics stated once, here
::: left
| item | value |
|---|---|
| F-number (air-equivalent, image space) | 1.8 |
| focal plane | 3000 spatial × 500 spectral pixels, 18 µm |
| slit | 54 mm; spectral height 9 mm |
| band | 380–2500 nm, 4.24 nm per pixel |
| smile, keystone | < 0.1 pixel |
| SRF FWHM | < 1.5–2.0 pixels |
| CRF FWHM | < 1.5 pixels |
| radiometric chain | throughput, grating efficiency, QE, slit loss vs wavelength |
::: right
- **Smile:** centroid drift along the slit at fixed wavelength; **keystone:** drift across wavelength at fixed slit position; both in pixels from ray centroids.
- **SRF, CRF:** slit ⊗ line-spread function ⊗ pixel ⊗ Airy, FWHM in pixels (Mouroulis & Green 2018 definitions); the Airy diameter at 2500 nm is 11 µm, under the 18 µm pixel, so the incoherent chain is valid.
- **Design rules adopted from the literature:** distortion to ~1 % of a pixel at design, > 75 % of the energy in one pixel, degraded spots accepted for uniformity, the grating is the stop.
- **Realism (Jim):** as-built SRF is 2.5–3 pixels, slits are 2 pixels wide, the systems are photon-limited; smile and keystone drive the design.
~ The 54 mm slit is 2.3× EMIT's; Table 3 of the review places this spec in the ALIS class (3200 spatial pixels, 380–2500 nm).  Public comparison point: Carbon-I (arXiv:2505.22545), F/2.2, 2040–2380 nm.

## The Dyson as traced by the engine | Step R3 of the ladder with every aperture declared: the beams clear every body; the grating is 270 mm wide on a 0.70 m radius and the slit sits 18 mm from the band
- **The prescription:** `challenges/dyson5/dyson5_s3_r3.in`, seven elements; the block's spherical face (E2), the grating (E3, radius 0.71 m at the Dyson condition for a 220 mm block), and the focal plane on the flat face.  Bodies, rays and labels are read back from the engine.
- **Sizes from the traced footprints:** block face 85 mm across at the sphere (220 mm radius, 217 mm thick), grating 270 mm across (0.70 m radius) at 0.69 m from the face; slit at +6 mm, this wavelength lands at −12 mm on the same face.
- **Every surface carries an aperture** cut from the multi-field, multi-wavelength footprint plus 5 mm, and a clearance check scores every beam leg against every body it does not traverse: the worst pair on this deck is the returning beam against the detector carrier, +1.16 mm, with no cold shield.
![Left: the traced bundle and the solid bodies in 3-D.  Right: the dispersion plane (Y-Z).  Positions in metres, global frame.](figs_dyson/dyson5_s3_r3_views.png)

## The Offner as traced by the engine | The all-reflective reference at F/2.8, re-posed at 0.22 R: the two mirror zones with the grating between the beams, clearing its mount by 10.2 mm
- **The prescription:** `challenges/dyson5/dyson5_s1_offner.in`, five elements; the concave mirror (E1, used twice), the convex grating (E2) at the stop, the focal plane at the slit's conjugate.
- **At the first slit offset (6 mm) the grating sat inside the slit-to-mirror beam;** at 0.20 R the beams still clipped its mount by 4 mm; at 0.22 R they clear by 10.2 mm, and the concentric form there is astigmatic (15 px), so the classical corrections are applied under the ladder's operands: convex radius × 1.0034, second concave zone × 0.951 with a 0.3 mm offset.
- **As corrected:** keystone 0.028 px, CRF 1.20 px, SRF 3.69 px, 23 % of the energy in one pixel.  It stays in the deck as the reference form and as the ground where the physical-optics twin is validated (order 0: pupil OPD 8e-11 m, 85 % in one pixel).
![Left: the traced bundle and the solid bodies in 3-D.  Right: the dispersion plane (Y-Z).  Positions in metres, global frame.](figs_dyson/dyson5_s1_offner_views.png)

## The concentric Dyson and its scaling law | The classical condition verified by exact trace; at this slit the concentric block needs r ≥ 213 mm
- **The Dyson condition** R_g = n·r / (n − 1) is the blur minimum by exact 3-D trace: 0.67 µm at the condition against 4–5 µm at ±5 %, image at −h, residual fifth order (h⁴·¹⁷, r⁻³·¹⁹).
- **The scaling answer, before any element was added:** a single concentric block meets a quarter-pixel corner blur only at r ≥ 213 mm, about twice flight scale, which is why the flight forms add an asphere, a meniscus or an air gap.
![Rms blur at the slit corner (h = 28.2 mm) of a concentric silica Dyson at F/1.8 versus block radius, exact trace at order 0: blur ∝ r⁻³·¹⁹; the quarter-pixel budget is met at r = 213 mm (grating radius 687 mm).](figs_dyson/dyson5_s0_scaling.png)

## The exact chain and the engine agree | Every engine ray re-traced by the chain from the engine's own launch directions lands within 1e-9 m
::: left
- **The check:** the engine's chief ray and every ray of its grid are re-traced through the chain's surfaces from the engine's launch directions; agreement is per ray, every surface, 1e-9 m; the dispersed band spans the focal plane to 1 %.
- **With the aspheric block face (ladder step R2):** agreement 1e-12 m per ray; a sphere-only chain misses by 0.26 mm, so the check is not vacuous.
- **The physical-optics twin:** a far-field terminal on a reference sphere upstream of the focal plane, re-posed per slit position and wavelength; on the Offner at order 0 the pupil OPD is 8e-11 m and the Airy spot puts 94 % of its energy in one 18 µm pixel.
::: right
- **What the checks found in the engine, now fixed** (details in Backup): glass names were ignored at load; propagation inside glass used the vacuum wavelength; the grating groove period was held constant along the curved surface instead of along the chord; the grating's path-length jump had the matching defect.  Each has a regression test written from the physics, not from the engine.
~ Regression tests: `tSpectrometerRx` (3), `tGratingOpl` (2), `tGratingImmersed` (4), `tGlassDispersion` (3), `tPropMedium` (3); mmacos fast suite 501 pass, 0 fail (2026-10-01).

## The seed, scored | The concentric seed at a 220 mm block meets smile; keystone and the corner blur are the work
- **Smile is met at the seed** (below 0.01 pixel on both forms); the field-angle map is linear along the slit as designed.
- **SRF sits at 2.02 pixels on both forms**, the floor set by the 2-pixel slit; re-scored on the engine after the grating fixes (the earlier ramp with wavelength was the groove-period defect, not the design).
![Scorer maps for the two seeds at 220 mm (rows: Dyson, Offner): field-angle map, smile (pixels), SRF FWHM (pixels), ensquared energy in one pixel, versus slit position and wavelength.](figs_dyson/dyson5_s2_maps_rc.png)

## The departure ladder | Each step solved on the exact chain with pixel-unit smile and keystone in the merit from the first pass, then emitted and scored in the engine
![Engine-scored maxima per ladder step at the seed's 220 mm block (left: smile, keystone, CRF FWHM in pixels; right: minimum ensquared energy in one pixel against the 0.75 design rule).](figs_dyson/dyson5_s3_ladder.png)

## The departure ladder, by the numbers | Keystone falls 37-fold at fixed scale; the meniscus pays only with the block offset and grating radius moving with it
::: full
| step | what moves | keystone (px) | smile (px) | CRF FWHM (px) | ensquared, 1 px | length (mm) | grating footprint (mm) | glass (L) |
|---|---|---|---|---|---|---|---|---|
| R0 | concentric seed | 0.095 | 0.006 | 2.28 | 0.44 | 708 | 275 | 22.2 |
| R1 | grating-radius factor, face offset | 0.038 | 0.003 | 2.10 | 0.48 | 704 | 274 | 22.3 |
| R2 | + conic and h⁴, h⁶ terms on the block face | 0.037 | 0.004 | 2.13 | 0.47 | 704 | 274 | 22.3 |
| R3 | + block centre off the grating's, along the dispersion | 0.011 | 0.008 | 2.10 | 0.48 | 698 | 272 | 22.1 |
| R4a | meniscus corrector alone | 0.063 | 0.012 | 2.49 | 0.33 | 698 | 270 | 22.3 |
| **R4** | **meniscus with all variables (the compact variant)** | **0.0026** | 0.005 | **1.33** | **0.76** | 694 | 269 | 22.7 |
| size alone | block radius free (341 mm) | 0.024 | 0.001 | 1.34 | 0.74 | 1099 | 427 | 82.8 |
~ SRF FWHM is 2.02–2.05 px on every step: the 2-pixel-slit floor.  Length = slit plane to grating vertex; footprint = diameter of the grating hits over slit centre and ends × band edges; glass = block cap + meniscus.  Record `dyson5_s3_trade.txt`.

## What each departure buys | De-concentring the block solves the distortion; the face asphere is inert; the blur is bought only by size or by the compact form
::: left
- **Keystone:** 0.095 → 0.011 pixels at fixed scale, nine times inside the 0.1-pixel specification and at the literature's 1 % design rule; smile stays below 0.01 pixel throughout.
- **The block-face asphere does nothing here:** step R2 converges to 10 nm sags, and a one-dimensional scan about R1 is a steep bowl centred at zero.  An axisymmetric figure acts on every field point alike and cannot cancel a residual that grows as h⁴ along the slit.
::: right
- **The blur does not move on any step, and the reason is on record:** the dispersed image sits 12 mm from the order-0 image, so the slit corner's effective field height is about 32 mm, and the h⁴/r³ law predicts the measured 0.4 px rms.
- **Two routes buy it:** size alone (0.74 ensquared at a 341 mm block, 1.1 m long, 83 L of glass) or the compact variant: a meniscus corrector in the air gap, solved together with the block offset and the grating radius, reaches the same at the 220 mm block — 63 % of the length, 27 % of the glass.
![The free-radius ladder: ensquared energy reaches the 0.75 rule only as the block grows to 341 mm.](figs_dyson/dyson5_s3free_ladder.png)

## The compact variant, as traced by the engine | Step R4: the meniscus corrector in the air gap between the block and the grating; the design of record at a 220 mm block
- **The prescription:** `challenges/dyson5/dyson5_s3_r4.in`, nine elements: the block's flat face, its spherical face, the meniscus (two surfaces), the grating, and the return through the meniscus and block to the focal plane.
- **Bodies, rays and labels read back from the engine,** apertures declared on every surface; the worst clearance on this deck is the returning beam against the detector carrier, +1.19 mm with no cold shield.  Any shield height fails the check, so a fold prism at the slit is the next step's first item.
![Left: the traced bundle and the solid bodies in 3-D.  Right: the dispersion plane (Y-Z).  Positions in metres, global frame.](figs_dyson/dyson5_s3_r4_views.png)

## The compact variant in section | Dispersion plane and slit direction, to scale: block, meniscus, grating; slit and focal plane on the block's face
![Left: the dispersion plane, slit centre, rays at 380 / 1440 / 2500 nm (blue / green / red).  Right: the slit direction, slit centre and ends at 1440 nm; the focal plane in magenta.  Block shaded; 100 mm scale bar.](figs_dyson/dyson5_s3_layout_r4_rc.png)
- **Sizes:** block radius 220 mm (217 mm thick, 85 mm face), meniscus 4 mm thick at radii near 2 m, grating 269 mm across on a 0.69 m radius; slit plane to grating vertex 694 mm.
- **The solve ends on its bounds** (plate thickness and radii), and three solves from the same seed found a multimodal landscape; the global search over the meniscus is the next step.

## The compact variant, scored | Ray-centroid keystone 0.0026 px, smile 0.005 px, CRF 1.33 px, SRF 2.03 px; ensquared energy 0.76 geometric and 0.75 with diffraction
![Engine-scored maps for R4 on 18 µm pixels, per slit position and wavelength: field-angle map, keystone, smile, SRF FWHM, CRF FWHM, ensquared energy; each convention on its colorbar, maxima in the labels.](figs_dyson/dyson5_s3_maps_r4_rc.png)
- **Open:** the propagated PSF's centroid parts from the ray centroid by up to 0.12 px on this design (0.001 px on the seeds, 0.03 px on R3), growing with wavelength, with the path-length check green.  A detector measures the intensity centroid, so until this is resolved the keystone that counts against the 0.1 px specification is the wave number.

## The trade at a glance | Distortion, response functions and ensquared energy per step, both records; the specification and the literature rules as lines
![Left: keystone and smile maxima (log scale) against the 1 % rule.  Middle: CRF and SRF FWHM against the 1.5 px specification and the 2-px-slit floor.  Right: minimum ensquared energy in one pixel against the 0.75 rule.](figs_dyson/dyson5_s3_trade_rc.png)

## The propagation twin | A far-field terminal on a reference sphere, re-posed per slit position and wavelength, agrees with the ray centroids to 0.0013 px
- **What it checks:** the ray-side scorer's centroids and response functions against physical-optics propagation through the same decks; across a grating the comparison is made in phase (modulo λ), never as unwrapped lengths.
![Left: the Offner at order 0, PSF (log scale) with the 18 µm pixel box, 85 % ensquared.  Middle: wave-minus-ray centroid over the slit and band, Offner and Dyson seeds.  Right: SRF and CRF from the propagated PSF at the slit centre.](figs_dyson/dyson5_s2w_twin_rc.png)

## Slit diffraction loss | A 36 µm slit at F/1.8: the engine's far-field leg against the sinc² closed form
![Loss past the grating versus wavelength: engine far-field leg (points) against the sinc² form (line).](figs_dyson/dyson5_s2l_slitloss.png)
- **Under 2.5 % over the band** by both; the engine runs above the closed form at the long end (2.48 % against 1.83 % at 2500 nm) and below it at the short end, and that factor is an open item in the record.

## Next steps | Clearance on the design of record, the global meniscus search, and the telescope that feeds the slit
::: left
- **R5, the slit-to-detector split:** the returning beam clears the detector carrier by 1.2 mm with no cold shield; a fold prism at the slit (the literature's answer) is the next layout step, scored under the same clearance check.
- **Beat 4, first:** settle the wave-versus-ray centroid difference (amplitude-weighted ray centroid from the engine's pupil amplitude; the twin's window and sampling as the null), then the global search over the meniscus and the native multi-wavelength optimizer with the centroid that a detector sees as the operand.
::: right
- **Beat 5, the telescope:** EMIT parameters — 420 km, 60 m ground sample, 0.143 mrad IFOV, focal length 126 mm, 70 mm aperture at F/1.8, 24.6° cross-track field onto the 54 mm slit, telecentric and flat; the telescope's exit pupil on the grating.  Then spectrometer and telescope traced end to end as one deck.
- **Then the physical-optics twin on R4** and the radiometric chain against the band.

## Run it yourself | One parameter file, one runner, every stage through it
::: full
```
run('<path>/mmacos/mmacos_setup.m');
addpath('<path>/mmacos/challenges/dyson5');
OUT = dyson5_run();                                      % stage s0: the scaling law (engine-free)
OUT = dyson5_run(struct('stages', {{'s1','s2','s3'}}));  % emit, score, ladder
```
- **Knobs** live in `dyson5_params.m` (F-number, pixel, band, slit, block radius, ladder options); every record in this deck is a stage output under `challenges/dyson5/`.
- **Engine checks:** `./run_mmacos_tests.sh tSpectrometerRx` and the four grating / glass / propagation classes, all in the fast suite.

## Backup
Diagnostics, the engine findings, and the conventions behind the main path.

## Backup: four engine findings, each with a regression test | Found by the challenge's checks on 2026-09-30 and 10-01, fixed in the engine the same day
::: full
| finding | symptom | fix | test |
|---|---|---|---|
| glass names ignored at load | a `GlassElt= Silica` element traced as air in every binding | the glass catalog is blanked once at allocation, not on every load | `tGlassDispersion` (3) |
| vacuum wavelength inside glass | a propagation leg in a medium ran at the wrong Fresnel number by n | every kernel takes λ/n of the leg's medium | `tPropMedium` (3): z in n ≡ z/n in vacuum |
| groove period along the surface | 2.8 / 3.3 px rms spectral blur at 2500 nm (Offner / Dyson) | the local grating vector is the un-normalised projection of the ruling direction (chord-ruled) | `tGratingImmersed` (4), `test_grating_chord` |
| grating path-length jump | 4 waves rms of pupil OPD where the rays converged to 0.05 µm | the jump is the groove count from the vertex along the ruling direction | `tGratingOpl` (2) |
~ Flat gratings are unchanged by the last two.  A fifth item, a short `AsphCoef=` line killing the host process, now pads with zero and warns (`tRxShortAsph`).

## Backup: records and tools behind each slide | Every number traces to a runner stage and its record file
::: full
| slide | runner stage | record | tool |
|---|---|---|---|
| scaling law | s0 (engine-free) | `dyson5_s0_scaling.txt` | `design/src/dyson_layout.m`, `dyson_scaling.m` |
| as traced (both forms) | `dyson5_view_figs.m` | `*_view3d.png`, `*_viewyz.png` | `macos.view_rx` on the deck |
| seed scored | s1 emit, s2 score | `dyson5_s1.txt`, `dyson5_s2.txt`, `dyson5_s2_maps.png` | `spectrometer_geom.m`, `spectrometer_rx.m`, `spectrometer_score.m` |
| chain vs engine | gate | `tests/tSpectrometerRx.m` | per-ray re-trace, 1e-9 m |
| propagation twin | s2w (opt-in, model 512) | `dyson5_s2w.txt` | `spectrometer_wave.m` |
| departure ladder, R4, trade | s3, s3free, s3b | `dyson5_s3.txt`, `dyson5_s3free.txt`, `dyson5_s3_trade.txt`, `dyson5_s3_layout_r*.png`, `dyson5_s3_maps_r*.png` | `dyson_ladder.m`, `dyson5_trade.m` |
| slit loss | s2l | `dyson5_s2l.txt`, `dyson5_s2l_slitloss.png` | far-field leg vs sinc² |
~ Conventions behind the scaling law: block index n(silica, 1 µm) = 1.450417; in-glass marginal half-angle u = asin(1 / 2nF); field height h from the Dyson axis on the flat face; the slit corner is h = hypot(27 mm, 8 mm).

## Backup: conventions pinned on the way | Engine facts the chain depends on, stated once
::: left
- **`ChfRayPos`** is where rays start and must lie before the first surface; at load the engine folds `zSource` in once, after which it is the physical source.
- **`macos.stop`** aims immediately and its first pass is one step short on this deck (3.6 mrad): declare the stop first, then write the exact chief ray.
- **A `Return` coincident with the focal plane drops the rays**; use a `Reference` upstream.
::: right
- **An immersed grating carries its glass on its own element**: the grating branch takes the exit index from the element as written.
- **`ray_info_get`** returns the outgoing direction at the traced element.
- **Centre of curvature** is the vertex plus |Kr| along `psiElt`; `KrElt` is stored as −|R|.
~ Reference document: `optical_design/SPECTROMETER_DESIGN_REFERENCE.md` (forms, conditions, metrics, references).
