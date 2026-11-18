# BRIEF: pupil image quality of the Twyman-Green interferometer (Fang Shi's request)

CCL for Dave, 2026-09-16.  Plan only; nothing built.  Owner: CCL builds the
stage (a day), runs on the Mac (cycle 3) or this box when TO's field-servo
runs leave it free.

## 1. The question, and the formulation: the DM is the stop

Fang asks how well the interferometer's camera sees the DM: is a poked
actuator imaged sharply, in the right place, at the right size, over the
whole 96 mm?  The thing being imaged is the DM, so the DM is the aperture
stop and the entrance pupil (Dave's ruling 2026-09-16), and the camera
plane is its exit pupil.  The rodgers2 form then applies without change: a
collimated source launched at the DM with the DM's clear aperture as the
stop, tilted about the DM by a field angle theta, and two things scored:

- **at the camera, the pupil image itself:** the exit-pupil surface from
  `macos.xps` + `macos.pupil_quality` (the index-matched exit rays at two
  fields cross at the image of that DM zone; the cloud's fit gives the
  image's defocus and astigmatism in mm, its low-order shape); the
  crossing position of every zone against its DM coordinate times the
  magnification (the pupil distortion map, in DM millimeters, against the
  runner's single global affine); and the crossings' spread over the band
  of field angles (the pupil blur per zone, in DM millimeters, against the
  1 mm pitch and the 0.2 mm detector pixel).  This is the direct picture
  Fang will recognize, and it needs no second source form.
- **at the focal plane (the mask seat), the rodgers2 set per tilt:** spot,
  strict wavefront, centroid against focal length times angle, pupil
  wander.  A tilt theta at the DM is a spatial frequency theta / lambda on
  its surface, so this is the pupil image's transfer, frequency by
  frequency, and the edge-only rule applies: quote the worst tilt in the
  actuator band.

## 2. Numbers that size the scan (from the bench of record, lens rig)

| quantity | value | source |
|---|---|---|
| DM, beam | 96 x 96 at 1 mm; beam 103 mm | tg96_params |
| focuser L2 | f 429 mm, F/4.2 at the internal focus | parts slide |
| field lens at the focus | f 43 mm, 21 mm clear | parts slide |
| pupil image | 7.8 mm (lens rig), 8.0 (mirror rig); 385 modeled px; ray magnification 9.88 / 10.10 DM-mm per detector-mm | lensuw2, oaploop22 |
| actuator Nyquist as a tilt | lambda / (2 x 1 mm) = 3.2e-4 rad = 65 arcsec | -- |
| the field lens's angular acceptance | 10.5 / 429 = 0.024 rad = 77 x Nyquist | -- |
| detector sampling | 5.0 px per actuator (lens), 4.7 (mirrors) | lensuw2, oapifo2 |

So the tilt field that matters is +/-3.2e-4 rad (the actuator band) with
+/-1e-3 as margin; the tail's geometric field is 75 times wider, so the
scan is well inside the field lens.  At the tail's numerical aperture
(0.024) the diffraction limit at the DM is 13 um, so the pupil-image
resolution is geometric, not diffractive: the score is the tail's
aberration with the DM as the object, in DM millimeters against the 1 mm
pitch and the 0.2 mm detector pixel.

## 3. The stage: `tg96_run` stage 'pupilq' (sheet-driven, per Dave's rule)

`P.pupilq` knobs: `tilts` (rad, default `[0 0.5 1 2 3.2 10]*1e-4` on two
azimuths), `zones` (the DM grid the crossings are reported on, default the
lit 96 x 96 binned to 12 x 12), `rig` ('lens' | 'oap'; both run), `tail`
(the tail of record per rig: the tuned lens tail; the geometric seed on the
mirrors).

1. **The deck.**  The runner's emitted test-arm deck with everything
   upstream of the DM removed and the source replaced by a collimated beam
   launched at the DM's station along its normal (zSource 1e22, Aperture =
   the DM's clear aperture, 96 mm; `ApStop` at the DM so the DM IS the
   stop), the DM flat (the pupil imaging is about the tail, not the
   surface), the reference arm shuttered, the camera's pupil return as the
   exit-pupil element for XPS.  One deck per rig, emitted and committed
   (`runs/pupilq_<rig>/pupilq_<rig>.in`).
2. **At the camera:** `macos.xps` at the pupil return for each tilt pair
   about the nominal; `macos.pupil_quality` for the surface; the crossing
   cloud mapped to DM coordinates through the ray history (each ray's DM
   hit is its zone), giving the distortion map (crossing vs magnification
   times zone, um on the DM) and the blur per zone (crossing spread over
   the tilt band).  One figure: distortion arrows + blur ellipses over the
   DM, the runner's own.
3. **At the focal plane:** per tilt, the spot (rms radius, um and lambda
   F/D), the wavefront referenced to the chief-tied sphere with piston and
   tilt removed (nm rms; low-order Zernikes), the centroid against
   F2 x theta.  One figure: the three against theta, both azimuths.
4. **Both rigs**, and on the lens rig both tails (tuned and seed), so the
   tuner's effect on the pupil image is a number.
5. **Report section + the two figures** into `REPORT_bench_realism.md`
   (section 6, "Pupil image quality"); the deck gets one slide in the
   interferometer block: the distortion/blur map with the worst-tilt table.

## 4. What the answer looks like, and what would be a finding

Expected on the lens rig: the tail is focused by a numerical solve on a
single-actuator poke's sharpness, so the pupil PSF should be well under a
pixel (0.2 mm DM) at the center; the question is the edge of the 96 mm and
the distortion, which the runner's affine currently absorbs as one global
linear map.  A distortion residual above ~0.1 mm on the DM would mean the
affine is hiding a per-actuator registration error that a real bench must
calibrate (it would appear as the "lateral misregistration" TO found in
the lens rig's station figure, 62 nm against the engine, which this stage
would explain).  On the mirror rig with the seed tail, expect the same or
better at the center (its null is 0.029 nm) and a larger edge distortion
from the 25-degree fold.  Both results are one table for Fang.

## 5. Cost

Building the stage: a day of CCL (deck surgery, the collimated source with
the DM as stop, XPS per tilt pair, the zone mapping, two figures).
Running: 12 tilts x 2 azimuths x 2 rigs of single traces plus the XPS
pairs at model 512 (minutes); no battery, no loop.  The Mac can run it (cycle 3) the moment
the stage is committed.

## 6. Built 2026-09-16 evening: `tg96_pupilq.m` (tg_psi_dm96_oap), and three things the first passes settled

1. **A collimated source AT the DM measures a bench that does not exist.**
   On an ideal collimated beam the lens rig's detector leg lands 16 lam F/D
   of blur (1.07 um rms of wavefront) at the seat, because the leg was tuned
   to the collimator's actual beam, fed 25 mm inside its focus.  So the
   stage keeps the bench as built: the point source of record, the DM
   declared the STOP (`macos.stop`), and the field as a lateral shift of the
   source at the collimator's focus, which is a tilt about the DM by shift
   over the collimator's conjugate (832 mm on the lens rig).
2. **The seat marker in the interferometer's deck is not at the focus**
   (352 um rms of geometric spot there, 6.5 um of pure defocus in the
   wavefront): the runner parks it where the tuned tail wants the mask
   seat.  The stage finds the best focus from the rays themselves along the
   chief (pure geometry) and scores the spot there; the wavefront is scored
   with piston, tilt and focus removed by a fit -- the chief-tied sphere
   metric to first order.  The engine's exit-pupil sphere (the
   Return/Return/plane idiom) refused to load through the mex without a
   message; not chased tonight.
3. **The crossing cloud at the camera needs no engine command:** two
   traces a small field step apart, the closest point between each ray's
   two exit lines, mapped to the ray's DM coordinate from the trace at the
   DM.  Distortion, blur and the pupil surface follow from that cloud.

## 7. 2026-09-17: the simulation (CCL for Dave), and three things to rule on

Built and run: `tg96_pupilsim.m` (stages: the leg's coherent PSF per DM zone
from the rays; the DM field through those PSFs with the four-step readout; the
compromise detector plane), `tg96_pupil_engine.m` (the engine's propagation
through reference surfaces in the .in file), run-it-yourself
`tg96_pupil_batch.sh both`, sheet block `P.pupil`.  Report:
`tg_psi_dm96_oap/REPORT_bench_realism.md` section 7; deck slide 13.

**Findings.**
1. **The beam of record was the source cone, not the DM** (77 mm on the lens
   rig, 82 on the mirrors, on a 96 mm DM; nothing clips a ray).  Per your
   ruling the baffle is opened and the DM carries a 48 mm aperture: sheets
   changed (`tg96_params`, `zwfs_params`: R_BAFFLE 12.5 -> 18, D_LENS 60 -> 66,
   R_TO_AP 30 -> 28).  Every deck emitted from now on has the 96 mm beam;
   the record's rows were made on the 77 / 82 mm beam and stand as labeled
   until re-run.  **Rule:** re-run the record's rows on the 96 mm beam now
   (the Mac has cycles) or after the interim report.
2. **Every DM mode is observable on both rigs.**  Lens rig as tuned: Nyquist
   gain 0.99 center to 0.95 worst at the edge, half-Nyquist and below within
   0.3%, a single actuator 0.99, the 30 nm working surface back to 1.2 nm.
   Mirror rig: 0.997 worst, 0.23 nm.  Amplitude cross-talk up to a third of
   the phase modulation at the lens rig's edge, which the four-step ignores.
3. **The lens rig's detector is not at the pupil image**: 2.6 mm ahead of the
   on-axis image, 4-6 mm ahead of the edge zones' (the null-tuned tail is
   blind to pupil defocus).  Moving it 4.3 mm downstream takes the Nyquist gain
   to 0.993 worst and the working surface to 0.4 nm; the residual is the field
   lens's own astigmatism (1 mm split at the edge).  The mirror rig's geometric
   seed sits at the image (within 0.4-1.3 mm).  **Rule:** add the image position
   (from the rays) to the tail tuner's objective, or trim DET_TRIM + 4.3 mm on
   the lens rig and re-read the rows.
4. **The engine's plane-to-plane check is built but cannot yet check this
   leg.**  A geometric-to-physical hand-off on a reference sphere is exact only
   at a pupil conjugate, and the detector leg has none before the detector
   (the DM's image by the focuser is virtual, 876 mm past the focus); the seed
   on any upstream sphere is the undiffracted DM pattern -- the effect itself.
   Measured, not argued (the seeded phase matches the OPD to 0.3 waves; the
   propagated flat pupil is a clean disc; the modes come out wrong under every
   zElt convention).  The faithful form is the CTB's station-to-station chain
   (every leg from the DM propagated, no hand-back to rays), half a day with the
   zElt conventions validated leg by leg.  **Rule:** now, or after Friday's
   report.  Until then the zone-PSF model is the result of record, with the
   paraxial standby chain as its check.
   **Later the same day:** the station-to-station chain was built
   (`tg96_pupil_s2s.m`).  On the lens rig of record it cannot work: that
   rig's collimated space carries 41 waves of curvature (the collimator fed
   25 mm inside its focus), the rays walk off the propagation grid between
   the NFPlane legs, and the per-index lens phase lands on the wrong pixels
   (1.5 waves of spurious aberration, a 10x focal spot).  The mirror rig of
   record is fed at the focus; the chain is being validated there.
   On the mirror rig (truly collimated) the chain runs but is not yet right:
   the quartet's far-side pupil reads 20% small against the rays and the
   readout is not a gain map under any convention.  Next checks are listed in
   REPORT_bench_realism section 7 (the focal spot after the entrance sphere,
   the far sphere's pitch, the exit step alone with a 10 mm defocus).  Half a
   day more; the zone-PSF model stands as the result of record meanwhile.

## 8. Improving the pupil image (Dave's ask, 2026-09-17 afternoon): assessed

`tg96_pupil_options.m`, five variants of the lens rig (REPORT_bench_realism 7.1):

- **The tail geometry is the whole story.**  The null tuner moved the field
  lens to 39.8 mm past the focus and bent it; the SEED tail (10.8 mm, the
  mirror rig's geometry) images the DM to within 1 mm, astigmatism 0.08 mm,
  Nyquist gain 0.9992 as built and 0.9999 with a 0.7 mm detector move, the
  30 nm surface to 0.06 nm, distortion 0.003 mm (vs 0.27).  Cost: the
  flat-DM null 9.1 nm vs 0.134 tuned -- a fixed pattern the reference frame
  removes.
- **Detector move:** yes; second-order once the tail is right (0.7 mm).
- **Field lens / flattener at the detector:** no -- the tuned tail's bowl
  is a 5 mm-radius surface in image space (a -3.5 mm flattener), and the
  seed geometry removes the bowl instead.
- **True collimation:** worth doing for itself (halves the distortion,
  needed for the physical-optics chain and for any re-tune), not for the
  gain: 0.999 with the seed tail either way.  The record's source sits 14 mm
  inside the hyperbolic collimator's focus.
- **Rule:** lens rig -> seed field-lens station held, re-tune the conic +
  trim with the image surface in the objective; mirror rig -> the 0.6 mm
  move; collimator fixed before the next tune; rows re-run on the corrected
  tail + the 96 mm beam together.
