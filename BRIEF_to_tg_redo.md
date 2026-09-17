# BRIEF (draft): the Twyman-Green analyses redone on the improved bench, the physical-optics chain on it, the field servo parked until then

CCL for Dave's review, then TO.  2026-09-17.  Companion to `BRIEF_to_gauge_close.md`
(items 1-7) and `BRIEF_pupil_quality_tg.md` (sections 7-8).  Rules of section 0 of the
gauge-close brief hold (one model-1024 MATLAB on this box, kill by PID, never edit a
running sequence, every rms knob in mm, push on Dave's word).

## Read-up first (TO, after a clear or a compact)

Clear rather than compact: this brief and the slice carry everything this work needs,
and a fresh context reads them faster than a compacted one re-derives them.  Then, in
this order, before touching anything:

1. `macos/CLAUDE.md` (the root rules) and, when you touch the engine's propagators for
   package B, `macos/macos_f90/CLAUDE.md` (the nested one is not re-injected after a
   compact).  Memory: `MEMORY.md`, then `project_pupil_quality_tg` (the pupil work and
   the engine hand-off lesson), `feedback_sequenced_batch_jobs` (one model-1024 MATLAB
   on this box; wrappers wait), `feedback_running_script_edit`,
   `feedback_shell_self_kill` (kill by PID, never `pkill` MATLAB by name),
   `project_reflective_front_end`, `reference_gauge_flow_ifo`.
2. `macos/CURRENT_SLICE.md`, the two blocks dated 2026-09-17 (the pupil simulation, the
   stop change, the station-to-station chain's state) and 2026-09-16 evening.
3. This brief, whole.  Then `BRIEF_to_gauge_close.md` section 0 (the rules) and section
   7 (the field servo's state, where it resumes), and `BRIEF_pupil_quality_tg.md`
   sections 7-8 (the day's findings and Dave's rulings).
4. `MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/REPORT_bench_realism.md`
   sections 6, 7 and 7.1 (the numbers this brief quotes, with their run tags), and the
   headers of `tg96_pupilsim.m`, `tg96_pupil_options.m`, `tg96_pupil_s2s.m`,
   `tg96_pupil_engine.m` (what each does and what it found; `tg96_pupil_engine` is the
   single-sphere hand-off that cannot work here, kept for the record).
5. The sheets: `tg96_params.m` (the bench block with today's stop change and its
   comment, `P.pupil`, `P.bench.tail_*`, `P.bench.SRC_AT_FOCUS`, `coat_oap`, the
   realism knobs `PLATE_SUB` / `MASK_SUB` / `BS_T` / `EDGE_MARGIN`) and `zwfs_params.m`
   (the same bench block); `tg96_tail.m`'s header (the tuner and its winner gate: the
   objective you are about to change).
6. The runs to have open: `runs/pupilsim_lens`, `runs/pupilsim_oap` (the simulation of
   record: `<tag>_report.txt` and the four figures), `runs/pupil_options/
   pupil_options_report.txt` (the five variants), `runs/pupils2s_oap/
   pupils2s_oap_report.txt` and `runs/pupils2s_lens/...` (the chain's state, the
   numbers in section 2), `runs/pupileng_lens` (the hand-off that failed, with why).
7. Git: both repos on `dev-candidate`; the day's commits are LOCAL and unpushed
   (resources 02929cc, 50358d4, dad95cf; macos c13d4c5 .. 1a4f23a); `git log
   --oneline -12` in each before you start, `git status --short` for what is tracked
   in `runs/` (deck, report, figures; never the .mat).  Push only on Dave's word.
8. Before running anything: `ps -C MATLAB` (the Mac's cycle-2 jobs are on the Mac, not
   here; anything here is yours or Dave's), and the sequencer rule: never edit a bash
   sequence that is executing.

Verify before you build on a number: every figure in this brief has a run tag; re-read
the tag's report rather than this brief if they disagree.

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

## 1b. Substrates: decided and baked in before anything is emitted (work package A0)

Today the decks of record carry IDEAL polarizing elements and no mask plate:
`PLATE_SUB = []`, `MASK_SUB = []` in both sheets.  The substrates were PRICED as
variants (sub22 / thk22 / item4bseq: splitter and compensator 10 mm, 2 mm fused-silica
plates under the five polarizing elements, a 2 mm mask plate, 4 mm lens edges: rows
hold at 1.00-1.01, the flat null 0.134 -> 20 nm re-tuned, the reading attenuated 13%
by floor, not gain) and then left out of the defaults.  The redo bakes the decided set
into the sheet defaults, so every emitted deck, render, station figure and parts list
carries it.  Proposed set (Dave to confirm each line):

| part | today's default | proposed default | where it shows |
|---|---|---|---|
| splitter and compensator plates | `BS_T` 1.5 x s = 2.6 mm | 5.833 x s = 10.0 mm (the realism set) | Rx faces, render, node clearance (the 1.39 mm shear), parts list |
| input polarizer, both arm QWPs, output QWP, analyzer | ideal, zero thickness | `PLATE_SUB [1.4585 2.0]`: 2 mm fused silica each (two faces around each ideal element; stations downstream shift t/2 per plate, the tail absorbs) | Rx faces, render, station figure, parts list |
| mask plate (the sensors' seat) | none | `MASK_SUB [1.4585 2.0]`, faces ahead of the sandwich's entrance sphere inside the gap (W040 + the t(1-1/n) focus shift, measurable) | sensors' Rx, render, parts list |
| singlet edges (L1, L2, field lens) | `EDGE_MARGIN` 2.0 mm | 4.0 mm (a 113 mm singlet's edge) | Rx thickness, render |
| OAP coating (mirror rig) | `coat_oap` per run ('none' for geometric rows; bareAl / protectedAl for the vector pair) | protectedAl (MgF2 quarter wave at 632.8 nm) as the default; 'none' only for the geometric equivalence gate | Rx coatings, the polarization rows |
| lens AR coatings, camera window, DM window | not modeled | stay out (state it on the parts slide): the scalar record has no Fresnel loss; the vector mode carries the uncoated faces' loss (T 0.66) | -- |

Package A0 is a sheet change plus one run of the pupil stage and the DECK stage on
each rig to confirm the rows still hold on the substrates with the tail re-tuned on
them (the tail's re-tune in package A runs WITH the substrates in).  The renders
(`<tag>_render.png`, view_rx bodies on their real sag and thickness), the station
figures and the deck's parts slides (6-8) are regenerated from the emitted decks;
CCL carries the parts slides.  Half a day inside package A.

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

1. A0 + A (substrates decided and baked in, then the bench: CCL + TO, a day) -> commit -> the pupil stage on both emitted decks.
2. B (the physical-optics chain, TO, half a day to a day) in parallel with the first
   C runs on the Mac (mirror rig: DECK, station, ladder, noise, loop, descent -- the
   cycle-3 runbook).
3. C on this box (lens rig, one model-1024 MATLAB at a time).
4. Deck + reports + brief closed; then item 7.

## 6. Dave's rulings (2026-09-17 afternoon) -- this brief is now TO's

0. **The substrate set of section 1b: YES, as proposed** (10 mm splitter and
   compensator, 2 mm fused-silica plates under the five polarizing elements, the 2 mm
   mask plate, 4 mm lens edges, protected aluminum on the parabolas by default).
1. **The tail objective: whichever gives the best performance AS AN INTERFEROMETER.**
   Not a weight between the null and the bowl: the objective is the reading itself.
   Use the working-surface error the pupil stage measures (stage 2 of
   `tg96_pupilsim`: the recovered minus true 30 nm surface, piston and tilt out; 1.2 nm
   on the tuned tail, 0.06 on the seed) as the objective, or its stage-1 proxy where
   the full stage is too slow inside a tune (the band-edge phase rms over the pupil,
   40 s per evaluation; the two track each other in every case run today).  The null
   is reported beside it, never optimized.  The winner gate (the single-actuator row
   read through the ray affine) stays as the acceptance test.
2. **Keep the L2 conic re-solve** in package A.
3. **Camera: a camera with larger pixels.**  On the 9.8 mm pupil image the 385-px
   budget wants a 25 um pixel.  Two ways, both to be priced by the runner's camera
   report (`P.cam`): the sCMOS of record binned 4 (6.5 um -> 26 um, 377 px across the
   image, what the record already does at 7.8 mm), or an unbinned large-pixel camera
   (24 um class, 408 px).  Report both; the deck's camera line carries the chosen one.
4. **The Mac's cycle-2 results are pushed** (oapcap22 and what followed); CCL folds
   them into the deck as the record on the old beam, labeled; the redo supersedes
   them.

The original questions, for the record:

0. The substrate set of section 1b, line by line (10 mm splitter and compensator, 2 mm
   fused-silica plates under the five polarizing elements, the 2 mm mask plate, 4 mm
   lens edges, protected aluminum on the parabolas by default).

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

## 7. Mac cycle 3a landed (2026-09-17, resources b13cd38): the baseline, and what it says

Substrates in, the 96 mm beam, the seed tail, the record's collimator (`sub96_lens`,
`sub96_oap`, rows on the 30 nm surface, matrix on it):

| | lens rig | mirror rig |
|---|---|---|
| flat-DM null | 92.4 nm rms | 26.3 nm |
| 10 nm on one actuator | 1.003 / 1 pm / SNR 10645 | 0.986 / 2.2 pm / 6472 |
| 1 nm on the grid (112 / 120 sites) | 0.754 / 61.5 pm / 135 | 0.981 / 4.5 pm / 309 |
| capture to 10%, single / grid | 480+ / outside at the first rung | 480+ / 480+ |
| pupil image | 9.54 mm (341 px) | 11.07 mm (384 px) |
| clearance | compensator 16 mm vs 25 (1 of 6 node parts) | 1 of 5 node parts |

Reading for package A: the seed tail on the MISFED collimator with the plates in gives a
92 nm null (a third of a wave of wavefront); one actuator still reads to 1 pm but the
multi-site row at 1 nm falls to 0.75 through the fold.  That is the case for the
collimator fix and the reading-objective tune, in one number; the tuner's objective
(the recovered surface) penalizes a wrapping null by construction, so no separate
null term is needed -- but package A's gate must include the grid row at 1 nm on the
30 nm surface (>= 0.98) beside the pupil numbers.  The mirror rig, fed at its focus,
holds 0.98-0.99 with a 26 nm fixed-pattern null: its package-A change is the DET_TRIM
alone, as planned.  Package A0 detail: the compensator moves down the DM leg
(D_BS_CMP) to restore the 25 mm margin against the 56 mm beam.  Jobs C (pupil stage on
these decks) and D (the sensors at 385 on the mirror rig) follow on the Mac.

Jobs C and D landed (resources d4159da).  **Pupil stage on the sub96 decks** (seed tail,
substrates, 96 mm beam): lens rig -- image surface within +-0.3 mm of the camera,
band-edge phase 0.017 rad max, Nyquist gain 0.9999, distortion 0.008 mm rms, working
surface 0.075 nm; mirror rig -- within -0.6..+0.4 mm, 0.030 rad, 0.9996, 0.72 mm (the
parabola pair's mapping), 0.10 nm.  So the pupil imaging is solved on both rigs by the
seed station alone; the lens rig's 92 nm null is the mask plate's spherical aberration
plus the misfed collimator, which package A addresses.  **Sensors at 385 on the mirror
rig with the substrates** (`sub96_oapsens`, MASK_TRIM 0.632): stepped-Zernike rows
0.993 / 5 pm (single) and 0.995 / 4 pm (grid), against 0.991 / 4 and 0.999 / 4 without
the plate; capture to 10% S 35 / V 51 / P 42 nm (36 / 59 / 52 before): the mask plate
costs the self-referenced sensors 15-20% of capture and no gain; the pinhole's 100 nm
poke gate 3.5 pm (0.2 before), still passing; the vector pair's uncalibrated error
408 pm (633 before).  Package A0 stands as decided; the plate's spherical aberration is
part of the bench and the tail tune runs with it in.

## 8. Package A passed (TO, 2026-09-17 13:30; resources d1a3b57, local) -- and the gate list revised

Both rigs collimated, re-tuned on the reading, re-emitted (tags `redo_lens`, `redo_oap`;
REPORT_bench_realism section 8).  Gate record from the pupil stage on the emitted
decks at 129 rays across (the tunes ran at 65): Nyquist gain worst 1.0000 (lens) /
0.9994 (mirrors); band-edge phase 0.000 / 0.035 rad; amplitude cross-talk 0.005 / 0.036;
zone image surface flat to 0.001 mm (lens), tilt 0.43 and defocus -0.45 (mirrors); the
30 nm working surface recovered to 0.042 nm (lens; the record's tuned tail 1.22) and
0.092 (mirrors); single pokes 0.9998 peak, 0.999 width, < 0.4 um centroid shift.
**The distortion gate (< 0.01 mm rms) is withdrawn** (TO): the mirror rig's 0.72 mm is
the parabola pair's mapping, carried through every record it set while reading 0.997,
and the matrix absorbs it.  Package A's gate is therefore: gain >= 0.998 worst, band-edge
phase < 0.06 rad max, the working surface < 0.15 nm -- and, added from cycle 3a, the
DECK stage's grid row at 1 nm on the 30 nm surface >= 0.98 (package C's first item, on
the Mac).  Follow-ons TO queued (`redotail`): the seed tail on the same bench (what the
tune bought), the own-cone pupil check, the null decomposition.  Next: package B here;
package C on the Mac (cycle 3b, once these commits are pushed).

**Mac cycle 3b landed (resources 4f91561): the redo bench's rows.**  Lens rig
(`redo96_lens`, the redo tail): null 59.3 nm (the seed tail on the misfed collimator
gave 92), single row 1.004 / 1.4 pm / SNR 7741, grid row at 1 nm **0.993** / 10 pm / 246
(was 0.754): the package-A gate's grid row passes.  Mirror rig (`redo96_oap`): null
26.5, rows 0.986 / 0.981, unchanged from 3a by design.  Capture 480+ on both.  Pupil
image 11.3 / 11.1 mm (384 px).  The null left on both rigs is the mask plate's
spherical aberration in the F/4.2 beam (2 mm: about 0.1 wave rms), a fixed pattern the
reference frame removes and the reading-tuned tail does not chase; the null
decomposition prices it.  **For TO: the sensor sheet's seat scan failed on the mirror
rig** (`redo96_oapsens`: `MASK_TRIM re-scanned: -5.9285 mm (mask-plane peak/sum
2.288e-06)`, then the run died): the scan's start or range is the lens record's
(-5.6), and on the mirror rig with the plate the seat is +0.632 (cycle 3a passed
there).  The Mac reruns with the explicit value; the scan needs a start at the
plate's shift and a range that brackets both rigs.

## 7. PACKAGE A: DONE (TO, 2026-09-17).  The gates, the two calls for Dave, and what moved that was not on the list

Everything below is in `REPORT_bench_realism.md` sections 8.1-8.9 with its run
tag; `tg_psi_dm96_oap/runs/harvest_redo.sh` prints the record under those tags.
Resources commits `23131f6 .. ca70b0f` (the last one local, the rest pushed by
another lane's rebase); macos `5704c61`, `319e6d2`.

### The gates

| gate (section 1) | lens rig | mirror rig |
|---|---|---|
| exit-ray spread after the collimator < 1e-4 rad rms | **6.3e-09** PASS | **1.4e-08** PASS |
| focal spot at the seat < 1 um rms | **0.17 um** PASS | **0.09 um** PASS |
| mask marker within 0.5 mm of the ray focus | **0.000** PASS | **-0.0001** PASS |
| Nyquist gain >= 0.998 worst | **1.0000** PASS | **0.9994** PASS |
| band-edge phase < 0.06 rad max | **0.000** PASS | **0.035** PASS |
| distortion < 0.01 mm rms | 0.041 FAIL | 0.722 FAIL |

**The distortion figure is withdrawn** (8.4): it was read off a row measured on
the UNCOLLIMATED bench, where the same table already showed collimation RAISING
distortion; it is a residual AFTER the registration affine (4 % of an actuator
pitch on the lens rig); and the mirror rig has carried 0.72 mm through every
record it ever set while reading 0.997.

**The number that matters:** the bench reads its own 30 nm working surface to
**42 pm** (lens) and **92 pm** (mirror), against the record's 1.22 nm and
0.23 nm.  Single pokes 0.9998 peak, 0.999 width, under 0.4 um of centroid
shift.  The lens rig's amplitude cross-talk is 0.005 per unit phase, where the
record's was up to a third at the edge.

### Two calls for Dave

1. **The field lens's conic: -7.77 or the seed's -2.11?**  The tuned conic buys
   a factor of **2.3 in the reading** (97 -> 42 pm on the 30 nm surface, 25x of
   it at the lowest spatial band) and costs 1.8x in the distortion residual.
   Over the 1.2 mm the beam actually uses it is a 0.20 um departure from the
   sphere (0.15 um of it new) -- so the part to specify is a mild asphere over a
   2.4 mm clear aperture, and the 12 mm blank is 10x oversized for this beam
   whatever the figure.  `bench.tail_from_mat false` selects the seed.
2. **The compensator's clearance: +10.4 mm against a 25 mm spec** (8.6a), and
   it is the STOP RULING, not the plates -- a builder plate's scored radius is
   the beam plus 5 mm, so opening the beam from 51.4 to 59 mm costs margin
   twice; the record cleared by +25.6.  `D_BS_CMP` 200 -> ~225 mm physical gives
   +28.1.  A layout decision with a parts-list consequence, so it is yours.

### What moved that was not on the list

- **The lens rig's collimator was the wrong lens, not just misfed** (8.1-8.2):
  `L1_Kr` 236.866 is `(n-1) x (F1 - zsource)` exactly, principal plane and all.
  The radius is re-solved with the conic; the CONIC barely moves (-0.583016 vs
  -0.5829).
- **The mask seat belongs to the FOCUSER** (8.5): one global `MASK_TRIM` put
  the lens rig's 1.23 mm on the mirror rig, 0.69 mm of it wrong -- outside this
  brief's own 0.5 mm gate.  `P.oap.*` now overrides `P.bench.*` per optics.
- **The input polarizer's substrate defocuses the mirror rig's parabola**
  (8.5): `POL_IN 'source'` puts a 2 mm plate in the DIVERGING leg, worth 0.70
  waves; solved out with `P.oap.SRC_TRIM`.  The seat then lands on the mask
  plate's own `t(1-1/n)` to three figures.
- **The 59 nm flat-DM null is the BEAM** (8.9), 44 of it: lighting the whole
  96 mm DM takes the null from 9.8 to 53.6 nm.  The substrates cost 39 nm on the
  misfed bench and 5.6 on the collimated one -- a plate in a genuinely
  collimated beam is pure path.  It is a fixed pattern the reference frame
  removes (the gate record above was measured with it in place), but the raw
  four-step folds at +-158 nm, so `battery.unwrap` is the first thing to try if
  a package-C row table comes back speckled.
- **The lens rig's station residual is a FOLD, and this brief's hypothesis for
  it is dead** (8.6b): the pupil-image bowl is gone and the residual barely
  moved (62 -> 51 nm, against the mirror rig's 0.46).  The figure shows isolated
  pixels thrown by ~lambda/2.  Package C item, with the probe named in 8.6b.
- **A clearance guard was firing on its own parts** (8.6a): substrate faces are
  now grouped with the element they bracket.  Five rows at -94 to -102 mm were
  bookkeeping.
- **The cycle-3a runbook needs one edit** before it runs (noted at its top):
  job D pins `'bench.MASK_TRIM',0`, which overrides the zwfs sheet's `'scan'`
  and now seats the ZWFS DIMPLE -- which, unlike the interferometer's marker,
  IS the optic -- 0.63 mm off focus.

### One thing to settle before package C

CCL has tracked the redo tails under their run tags (`redo_lens_tail.mat`,
`redo_oap_tail.mat`) for the Mac's cycle 3b.  `tg96_run` looks up
`<tag>_tail.mat` first and `<optics>_tail.mat` second, so a run tagged
`redo_lens` finds them and a package-C run under any other tag falls back to the
RECORD's `lens_tail.mat` / `oap_tail.mat`.  Either package C runs under those
tags, or the redo tails are copied onto the canonical names.  Left alone rather
than overwritten under another lane's feet.
