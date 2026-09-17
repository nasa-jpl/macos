# BRIEF (draft): the Twyman-Green analyses redone on the improved bench, the physical-optics chain on it, the field servo parked until then

CCL for Dave's review, then TO.  2026-09-17.  Companion to `BRIEF_to_gauge_close.md`
(items 1-7) and `BRIEF_pupil_quality_tg.md` (sections 7-8).  Rules of section 0 of the
gauge-close brief hold (one model-1024 MATLAB on this box, kill by PID, never edit a
running sequence, every rms knob in mm, push on Dave's word).

## 0. Why, in one paragraph

The interferometer record (rows, capture, photons, servo, descent, station figures,
the reflective front end) was made on a bench with three defects that today's pupil
work exposed and quantified: the beam was the source cone, 77 mm on the lens rig and
82 on the mirrors, on a 96 mm DM (the outer actuator rings unlit); the lens rig's
collimator is fed 14 mm inside the focus of the hyperbola it should be, so its
"collimated" space carries 41 waves of curvature; and the lens rig's null-tuned tail
moved the field lens to its own focal length behind the focus, where the DM's image is
a bowl 2.6 to 6 mm off the camera (actuator-Nyquist gain 0.954 at the edge, distortion
0.27 mm, the 30 nm working surface read to 1.2 nm).  The seed tail station images the
DM flat (0.999, 0.003 mm, 0.06 nm), the mirror rig never left it, and true collimation
is what the engine's station-to-station chain needs.  The record's numbers are not
wrong for the bench they describe; they describe a bench we are no longer going to
build.  Redo them on the improved one, with the pupil-image stage and the
physical-optics chain as standard stages, then return to the coronagraph field servo.

## 1. The improved bench, exactly (all knobs in the sheets; nothing hand-edited in a deck)

Already in the sheets (resources 02929cc, 2026-09-17): `tg96_params` and `zwfs_params`
carry `R_BAFFLE` 12.5 -> 18, `D_LENS` 60 -> 66, `R_TO_AP` 30 -> 28 (the DM aperture =
the 96 mm actuator footprint), `P.clear.beam_r` 56.  The DM clips 11% of the rays and
is the stop in fact.

To add (work package A):

1. **Lens rig collimated for real.**  `SRC_AT_FOCUS` extended to `optics 'lens'` in
   `twyman_green`: the source at the collimator's focus, found by minimizing the
   exit-ray angular spread (the hyperbolic seed `L1_Kc = -n^2 = -2.25` leaves 6e-5 rad
   rms, ten times better than the record's 5.8e-4; `tg96_pupil_options` has the
   search, `coll_spread_`).  Then re-solve the two lens conics on the collimated
   geometry (`L1_Kc`, `L2_Kc`: the record's -0.5829 / -0.5826 were solved on the
   misfed beam; the exact forms for these plano orientations are Cartesian ovals, so
   the conic is a solve, not a formula).  Gate: exit-ray spread after L1 < 1e-4 rad
   rms; the focal spot after L2 at the marker < 1 um rms on the flat DM; the mask
   marker within 0.5 mm of the ray focus (today's marker is 5.5 mm off).
2. **The tail held at the seed station.**  `tg96_tail` gets `hold_fl_station true`:
   `D_MASK_FL` fixed at the seed (6.277 x s = 10.8 mm past the focus), only `FL_Kc` and
   `DET_TRIM` free; and a new objective term `image_surface`: the zone-image surface
   from a differential trace (what `tg96_pupilsim` stage 1 measures: rms and max of
   each zone's image offset from the detector plane, in mm), weighted so a 0.3 mm rms
   bowl costs as much as 1 nm of null.  The detector then sits at the image surface's
   mean (this is the `DET_TRIM` the objective finds; on the seed geometry it is
   0.7 mm from the thin-lens conjugate).  Winner gate unchanged (item 3).  Gate:
   `tg96_pupilsim` on the tuned deck reads Nyquist gain >= 0.998 worst, distortion
   < 0.01 mm rms, band-edge phase < 0.06 rad max; the flat-DM null is REPORTED, not
   gated (the seed tail's 9 nm is a fixed pattern the reference frame removes).
3. **Mirror rig:** `DET_TRIM` -0.6 mm (its image surface's mean); nothing else.  Gate
   as in 2.
4. **Both rigs re-emitted** (`tg96_run` stage B) and the pupil-image stage run on the
   emitted decks (`tg96_pupil_batch.sh both`) before anything else runs on them.  The
   two reports are the gate record.

Half a day.  Commit the builder and sheet changes before the first run of package C.

## 2. The physical-optics chain on the improved bench (work package B)

`tg96_pupil_s2s.m` is the station-to-station form as the CTB emits it (four collimated
NFPlane pairs at the chief pierces, the through-focus quartet with an entrance sphere
centered on the true focus and a far-side sphere at the field lens, the lens per
index, a sphere concentric with the exit beam, the scaled step to the detector).  It
cannot run on the lens rig of record (the 41 waves of curvature walk the rays off the
grid between the NF legs: 1.5 waves of spurious aberration at the entrance sphere, a
10x focal spot) and on the mirror rig it runs but is not yet right.  On the improved
bench (truly collimated), validate it leg by leg, in this order, each with a number:

1. **Collimated legs.**  After the four NFPlane pairs, the flat pupil's field phase at
   the entrance sphere must equal the ray OPD there to < 0.02 wave rms (today: 1.5
   waves on the lens rig).  A Nyquist sinusoid on the DM must reach the entrance
   sphere with its amplitude within 1% of the DM's (the 655 mm of Fresnel diffraction
   is real but the imaging undoes it; the check is that nothing else happened).
2. **The quartet.**  The focal-plane field's spot after the entrance sphere must be
   the Airy pattern of the pupil (first zero at 1.22 lam F/D = 3.3 um; today the
   mirror rig gives a 34 um blob).  The far-side sphere's pupil radius must match
   the rays' radius there to 2% (today 20% small: `pupils2s_oap`, 1.12 vs 1.40 mm).
   Suspects, in order: the pitch label the engine writes after NF2 (dx from the rays
   vs the scaled propagation pitch), the far sphere's `zElt` sign for the asymmetric
   case (the CTB's quartets are symmetric), the equivalent-radius estimator on a
   Fresnel-ringed edge.
3. **The exit step.**  The scaled sphere-to-sphere step's `zElt` convention (the one
   thing no deck of record uses): the four candidates in `tg96_pupil_s2s` against a
   10 mm known detector defocus, whose zone-PSF prediction is a Nyquist gain of
   0.95 mean (lens) / 0.85 (mirrors); the reversed sign gives ~0.7.  Then the flat
   pupil's disc must match the rays' (4.9 mm on the lens rig) to 2%.
4. **The readout.**  The gain-by-radius curves for the 2, 4, 8, 16 mm sinusoids, the
   six pokes and the working surface against `tg96_pupilsim` (the zone-PSF model of
   record).  Agreement within 0.5% in gain settles the model; disagreement is a
   finding about one of them, to be run down before either is quoted.
5. **Then the same chain on the mirror rig**, whose exit beam DIVERGES after the field
   lens (10.8 mm past the focus, inside its focal length): the exit sphere's center is
   upstream (`tg96_pupil_s2s` orients it by the measured crossing sign).

Half a day if 2 yields; a day if the quartet's pitch needs the engine read.  Report:
`REPORT_bench_realism.md` section 7 gets the validated chain; the deck's pupil slides
get the engine's curves beside the zone-PSF ones.  Once validated, the same chain is
the full physical model for the mask sensors' legs (their masks sit at the quartet's
plane), which is where it earns its keep.

## 3. The interferometer analyses redone (work package C): what, on which rig, gates

Every number the deck quotes for the interferometer, on both rigs, on the improved
bench, through the runner (`tg96_run` stages; `tg96_batch.sh` on this box, the Mac for
the mirror rig at model 1024).  In the order the deck reads:

| deck slide | record run tags | redo | gate against the record |
|---|---|---|---|
| Twyman-Green, four phase steps (rows on 30 nm) | lens_deck, oap_deck / oapifo2 | stage DECK both rigs | rows >= record's 0.99; the on-surface single row and its pm |
| Interferometer, station by station | stnoap (stnlens open: 62 nm misregistration) | both rigs, `stn` stage | flat 0.00 pm; 30 nm surface within 2% of the engine field; the lens rig's 62 nm must close or be explained (the pupil-image bowl was the size this predicted; the seed tail may close it) |
| pupil image quality + the DM's modes | pupilq_*, pupilsim_* | both rigs (package A gates) | Nyquist >= 0.998 worst |
| performance side by side | lens_deck, oap_deck | included in DECK | -- |
| capture range + its photon price | cap385 ladders, oapcap22 (Mac, running) | ladder both rigs (`wrap` stage on, `battery.unwrap`) | the beyond-fold fraction is the meter; unwrapped capture the number of record |
| photons for 1 pm | noise ladders, oapnoise22 | both rigs | -- |
| servo (loop) | loop_lens, oapifol2, oaploop22 | both rigs; the sensors' loops on the mirror rig too (oaploop22's S/V/P) | 3 pm from the record's photons; thermal floor ~10 pm |
| descent | descent_lens, oapdesc2 | both rigs, unwrap on | 3 pm by cycle ~20 from 100/200 nm |
| the reflective front end, measured | oapifo2, oapifol2, oapdesc2, oapsens22, vqw22, vmap22, oapuw2, lensuw2, oaploop22 | the mirror rig's block, on the 96 mm beam + the 0.6 mm move | the lens-vs-mirror rows (D1 100%, tail null, capture) |
| D1 window placement, D4 alignment | lens/oap `place`, `d4` | both rigs | D1 within 2 px >= record; D4 sensitivities |
| realism: plates, camera, polarization at the built angle | thk22, sub22, aoi_lens22, aoi_oap22 | both rigs | rows hold on the plates; the built-angle bias |
| the vector pair on the mirror rig | G4 rows (vqw22, vmap22) | mirror rig | the calibrated 0.054 pm |

Two days of runs with the Mac carrying the mirror rig (a cycle-3 runbook, CCL writes it
from this brief once package A lands).  What changes in the record, expected: the
lens rig's rows and station figure improve (the pupil image), the pupil image is 9.8
mm not 7.8 (the camera's binning and the 385-px sampling budget re-derived: `cam.bin`),
52 lit sites become the full 96 x 96 footprint (the sensors' 47 too), the capture
ladders shift with the lit set, and nothing about the mirror rig's rows should move
beyond the beam change.  The interim report of Friday 2026-09-18 goes out on the
CURRENT record, labeled as such; the redo lands in the next revision.

## 4. The coronagraph field servo (item 7 steps 1-3): parked until package C lands

State at parking: step 1's prescription corrected twice onto the measured pupil
(dichroic at Apodizer_Pst, 300 mm lens F/9.4, 2 lam/D dimple 11.8 um); steps 2-3 not
started.  It resumes exactly there, after the redo's deck update, unless Dave says
otherwise.  Reason for the order: the redo changes the bench every other deck slide
stands on; the field servo stands on the CTB, which is unchanged, and loses nothing by
waiting.

## 5. Order, and what the Mac does

1. A (bench, CCL + TO, half a day) -> commit -> the pupil stage on both emitted decks.
2. B (the physical-optics chain, TO, half a day to a day) in parallel with the first
   C runs on the Mac (mirror rig: DECK, station, ladder, noise, loop, descent -- the
   cycle-3 runbook).
3. C on this box (lens rig, one model-1024 MATLAB at a time).
4. Deck + reports + brief closed; then item 7.

## 6. Questions for Dave before this goes to TO

1. The tail objective's weight between the null and the image surface (proposed: a
   0.3 mm rms bowl = 1 nm of null).  Or drop the null from the objective and gate it
   only.
2. Keep the lens rig's L2 conic re-solve inside package A, or accept the record's
   tuned L2 with the collimated beam and let the tail absorb the 5 mm focus shift
   (cheaper, but the mask marker then stays off the focus).
3. Camera binning on the 9.8 mm pupil image: re-derive `cam.bin` for the 385-px
   budget (bin 4 -> the sampling changes), or hold the pixel count and accept the
   larger image.
4. The mirror rig's mask sensors (oapsens385, queued on the Mac on the old beam): let
   it run and label, or pull it into the redo.
