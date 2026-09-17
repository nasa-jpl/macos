# BRIEF: pupil image quality of the Twyman-Green interferometer (Fang Shi's request)

CCL for Dave, 2026-09-16.  Plan only; nothing built.  Owner: CCL builds the
stage (a day), runs on the Mac (cycle 3) or this box when TO's field-servo
runs leave it free.

## 1. The question, and two ways to ask it

Fang asks how well the interferometer's camera sees the DM: is a poked
actuator imaged sharply, in the right place, at the right size, over the
whole 96 mm?  The detector leg (focuser L2, the internal focus, the field
lens, the camera at the pupil image) is the imaging system; the DM is its
object; the aperture at the internal focus is its stop.  Two equivalent
formulations, and they use the same deck:

- **(B) Dave's, the rodgers2 form:** a collimated source at the DM, the
  entrance pupil AT the DM, tilted about it by a field angle theta, scored
  at the interferometer's focal plane (the mask seat).  A tilt theta at the
  DM is a spatial frequency f = theta / lambda on the DM surface, so the
  focal-plane image quality over the field of tilts IS the transfer of the
  pupil image, frequency by frequency: a spot that blurs or walks at angle
  theta is a DM frequency f that the camera sees attenuated or shifted.
- **(A) the direct form:** point sources ON the DM (a grid of DM
  positions), each radiating into the cone the tail accepts, imaged at the
  camera: the point-spread function of the pupil image in DM millimeters,
  its distortion against the runner's ray affine, and its wavefront.

(B) is what the rodgers2 tools score without change; (A) is the picture
Fang will recognize.  Run both from one stage; report (B) as the
assessment and (A) as its illustration.

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
azimuths), `dm_grid` (points on the DM for form A, default 5 x 5 over the
lit 96 mm), `na_stop` ('fieldlens' | a radius in mm at the focus), `rig`
('lens' | 'oap'; both run), `tail` (the tail of record per rig: the tuned
lens tail; the geometric seed on the mirrors).

1. **The deck.**  The runner's emitted test-arm deck with everything
   upstream of the DM removed and the source replaced: form B, a collimated
   beam launched at the DM's station along its normal (zSource 1e22,
   Aperture = the 103 mm beam), the DM as element 1 (its EP), the field
   angle set per point by the chief-ray direction; form A, a point source
   at each DM grid position with the cone the field lens accepts.  The DM
   flat (the pupil imaging is about the tail, not the surface).  Reference
   arm shuttered.  One deck per rig, emitted and committed (`runs/pupilq_
   <rig>/pupilq_<rig>.in`).
2. **Form B scores at the focal plane (the mask seat), per tilt, the
   rodgers2 set:** the spot (rms radius, um, and in lambda F/D); the
   wavefront referenced to the chief-tied sphere with piston and tilt
   removed (the strict metric; nm rms; low-order Zernikes); the centroid
   against F2 x theta (distortion, um and as a fraction); the exit-pupil
   surface from `macos.xps` + `macos.pupil_quality` at the camera's pupil
   return (defocus and astigmatism of the DM's image, mm) and its walk
   with theta (pupil wander, um).  The edge-only re-score rule applies:
   quote the worst tilt in the actuator band, not the mean.
3. **Form A scores at the camera, per DM point:** the spot in DM
   millimeters through the affine (the pupil-image PSF; against the 1 mm
   pitch and the 0.2 mm pixel), the centroid against the affine's
   prediction (pupil distortion map, um on the DM), the wavefront per DM
   point (nm rms).  One map figure: distortion arrows + spot ellipses over
   the DM, the runner's own figure.
4. **Both rigs**, and for the lens rig both tails (the tuned tail of
   record and the geometric seed) so the tuner's effect on the pupil
   image is a number.
5. **Report section + two figures** into `REPORT_bench_realism.md` (a new
   section 6, "Pupil image quality"), the deck gets one slide in the
   interferometer block: the form-A map with the form-B worst-tilt table.

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

Building the stage: a day of CCL (deck surgery, the two source forms, the
scoring, the figure).  Running: form B is 12 tilts x 2 azimuths x 2 rigs of
single traces at model 512 (minutes); form A is 25 point sources x 2 rigs
(minutes); no battery, no loop.  The Mac can run it (cycle 3) the moment
the stage is committed.
