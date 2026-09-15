# BRIEF for CCMac (Opus): the bench must be buildable -- node clearance, substrates, a real camera

From CCL for Dave, 2026-09-15.  Dave reviewed the deck's layouts: "the
collimator lens, input polarizer, compensator, analyzer and the BS
interfere severely, blocking the various beams.  This is not buildable."
And: the polarizers, quarter-wave plates, analyzer and the masks need
substrates of real thickness; the cameras are drawn small; the beam at
the camera must match a real pixel pitch.  Priority: the deck.  Standing
rules as before (dev-candidate; one model-1024 MATLAB at a time; every
number in a committed report; American English).

## 0. Work so that a dead session loses nothing (Dave, 2026-09-15)

A session can run out of tokens partway through.  Structure the work so
the benefit of everything already done survives that, and a successor
session (yours or anyone's) resumes from the record, not from memory:

1. **Order the items cheapest-complete-first and land each one whole**
   (code, gate, report section, run dir) before starting the next --
   TO's practice on the PDI lane.  Section 6 below is the order.
2. **Commit after every item, locally, into the shared tree** (`git add`
   your files only; a foreign-hunk check before adding a shared file).
   The commit is the checkpoint; Dave pushes.  Never hold a day's work
   uncommitted for a single big commit.
3. **Every run longer than a few minutes is a detached batch job**
   (`tg96_batch.sh` / `zwfs_batch.sh`: nohup + the lock, the log in
   `runs/<tag>.log`, the exit code recorded), never a foreground call
   in the session.  A detached run finishes after the session dies; its
   report and .mat are harvested by whoever comes next.  Chain runs in a
   sequence script (`runs/<name>seq.sh`) so the whole ladder proceeds
   unattended.
4. **Write the report as you go, numbers first**: each item appends its
   own section with its run tag the moment its run lands.  Keep a status
   table at the TOP of the report (item | state: done / running (tag,
   started, expected end) / not started) -- the successor reads that
   table first, then `runs/*.log` for the exit codes, and continues from
   the first item not done.
5. **Keep the context small**: grep the reports and tail the logs rather
   than reading them whole; do not paste long outputs into the session;
   let the runner's gates do the checking.
6. **If you see the budget ending**: commit what is committed-able,
   update the status table with what is running and where its output
   goes, and stop cleanly -- do not start a new item.

With this, a session cut at any point costs one uncommitted edit; the
runs keep going and the record shows exactly where to resume.

## 1. What is wrong, measured

`dm_gauge_lib/dmg_bench_clearance.m` (CCL, 2026-09-15): builds the
polarizing TG96 with the record's bench block plus overrides, traces
both arms, and for every physical part reports the smallest clearance to
any beam it is not in -- lateral distance of the beam's chief where it
crosses the part's plane, minus the beam radius (the 103 mm DM aperture)
minus the part's radius and an 8 mm mount.  Records of one part (a
plate's faces, its two passes) are grouped by name stem.  It also draws
the bench from above with a node panel (`zwfs_dm96/bench_bs7.png`,
`bench_bs22.png`, `bench_bs30.png`).

At the record's 7 deg, 8 of 9 node parts sit in another beam: output
QWP -109 mm, reference QWP -106, analyzer -106, compensator -74,
test-arm QWP -66, L2 -61, input polarizer -54, L1 -50.  Cause: the
Stage-A solve in `tg96_run` clears the three END bodies (DM, reference
flat, camera) with a 700 mm leg cap and prefers the smallest angle; it
never looks at the node.  Two legs 2*AOI apart clear a part at distance
d only if d sin(2 AOI) >= 103 + 8 (mm); at 7 deg that is 459 mm, and the
node parts sit at 17-257 mm.

Scan, with the output QWP and analyzer moved from 17 / 27 mm behind the
splitter to 160 / 170 mm (`D_RECOMB` 150, `D_RC_L2` 55: L2 stays at 207
mm, so the tuned tail is untouched) and nothing else moved:

| splitter | worst part (mm) | next |
|---|---|---|
| 7 deg | -109 (output QWP) | eight parts negative |
| 15 deg | compensator -27, output QWP -17, analyzer -11 | L2 +10 |
| 22.5 deg | compensator +16 | output QWP +55, analyzer +65, L2 +100, L1 +153 |
| 30 deg | compensator +57 | output QWP +178, analyzer +196, L2 +257, L1 +350 |

The one negative left at 22.5 / 30 is the reference QWP's "In" record,
which the builder places 25 mm behind the splitter; physically it is one
plate at the flat's end (the "Out" record).  Same for the test-arm QWP
(its "In" record sits after the compensator).

The sensors do not care: the ZWFS bench + battery at 30 deg (dev
resolution, `zwfs_dm96/runs/bs30_dev`) reproduces the 7 deg rows
(vector 0.9943 / 5 pm, 0.9958 / 1 vs 0.9942 / 5, 0.9961 / 2; G1 / G3
pass; the arm's channel phase difference 1.69 mrad vs 1.63).  The
plate's diattenuation rises from 0.5% to 10% but UNIFORMLY over the
pupil -- a non-term with the laser on the eigenaxis (V3's rule).  The
interferometer's polarization-snapshot form is the one thing that feels
the angle (item 4).

**Dave ruled 22.5 deg (2026-09-15).**  At 22.5 the plate's uniform diattenuation is 5.5% (7 deg: 0.5%; 30: 10%), the vector sensor's channel phase difference 1.66 mrad (1.63 / 1.69), the transmitted beam's shift in the 2.6 mm plate 0.36 mm (`runs/bs22_dev`).

## 2. The builder and the Stage-A solve (yours)

1. **Extend Stage A to the node parts.**  Each of L1 (at `D_L1_BS`),
   the input polarizer, the compensator (`D_BS_CMP`), the output QWP and
   analyzer, and L2 gets the constraint `d sin(2 AOI) >= beam_r +
   part_r + MARGIN`; the solve returns the smallest angle that clears
   ALL of them (near-normal preferred, as now) and the report prints the
   clearance table from `dmg_bench_clearance` -- gate: every part
   >= +25 mm.  With the current distances that solve gives about 27
   deg; Dave picks the round number.
2. **Forward `D_RECOMB` / `D_RC_L2`** (or add `D_OUT`, the output optics'
   distance behind the splitter) so the output QWP and analyzer sit at
   the clear distance just ahead of L2, in collimated space.
3. **The arm QWPs at the retro end for BOTH passes** (the "In" record at
   the same position as the "Out", `D_QWP` from the retro), so the
   layouts show one plate where it is.
4. **Re-solve the OAP folds** for the new angle (the Stage-A fold solve).

## 3. Substrates and thicknesses (Dave): real parts, drawn and traced

Dave, from the drawings (which are to scale): "Is the compensator really
< 1 mm thick?  Are the lenses?"  The model's parts are the 56 mm rig's
scaled by 1.714: the splitter and compensator are 2.6 mm plates, L1 is
5.1 mm at the center with a 1.8 mm edge, L2 7.8 mm with 1.6 mm, the
field lens 4.2 mm; the polarizers, plates and analyzer have no thickness.
Make them real: a 103 mm plate splitter and compensator at 10 mm (the
flatness a 4-inch plate needs; the transmitted beam shifts 1.4 mm at
22.5 deg, which the builder's chief-ray path carries), lenses with a
real edge (3-5 mm for 103 mm singlets), and:

Every thin element is a plate of real thickness: the input polarizer,
the two arm QWPs, the output QWP, the analyzer, the vector sensor's QWP,
and the mask (the etched dimple plate, the pinhole plate, the
metasurface).  Model them as what they are -- two refracting faces
(fused silica or quartz, 2-3 mm) around the ideal element -- so the
layouts show the glass and the trace carries it:
- in the collimated legs the plates add path only: one QWP per arm
  (balanced), the input and output optics common to both arms;
- the mask plate sits in the F/4.2 converging beam: spherical
  aberration W040 = t (n^2 - 1) NA^4 / (8 n^3) = 19 nm (0.03 wave) for
  2 mm at NA 0.12, and a focus shift t (1 - 1/n) = 0.7 mm that the tail
  retune absorbs (the mask is etched on the plate's exit face: put the
  faces before the sandwich's entrance sphere, the mask plane in air
  behind them);
- the vector QWP in the F/3.6 pupil-image leg: 34 nm for 2 mm, on an
  image whose resolution requirement is one actuator.
Builder: a `'substrate', [n t]` option on `add_polarizer`,
`add_waveplate` and the mask seat.  Then the tail retune (tg96_tail /
the ZWFS S1 rounds) with the plates in.

## 3b. The analyzer (Dave): how is it implemented?

The snapshot form's four analyzer orientations must be simultaneous.
A rotating stage would make it a sequential scan (the PZT form's drift
class without its step error).  Implement it as a polarization camera
(a micro-polarizer array at 0 / 45 / 90 / 135 deg on the pixels: at the
9.4 mm pupil image a 3.45 um polarization sCMOS gives 2700 pixels across,
680 per orientation); the model's rotating ideal analyzer between four
frames is that measurement.  The deck says so; the parts list names the
camera.  The vector sensor's two cameras stay (its two channels are
circular states, split by the cube).

## 4. The camera (Dave): a real pixel pitch

The record's pupil image is 9.4 mm across 385 modeled pixels = 24.4 um
per pixel.  State the camera: at that image a 6.5 um sCMOS (2048 x
2048, 13.3 mm) puts 1450 pixels across the pupil, binned 4 to the
modeled 360; a 13.5 um camera 700; a 24 um CCD 390.  Shrinking the image
to 385 pixels of a 3-5 um camera (1.3-1.9 mm) would need an F/0.7 field
lens -- not realistic; keep the image and bin.  Consequence for the
light: smaller pixels HELP the well-depth problem (1e14 photons per
measurement over 1.6e6 pixels is 6e7 electrons each: 600 frames at a
1e5-electron well, 6 s at 100 fps, against 100 s for 1e5 pixels).  The
model's 385 pixels per pupil is the sampling floor, not the camera.
Layouts: draw the camera body at its real size; the parts list names the
camera, the pitch, the binning.

## 5. After the round: what re-runs

- Layouts and parts, all rigs (lens, OAP, P/SRI with TO, vector,
  pinhole), the recipe as before; the clearance table in every report.
- The snapshot form's polarization at the new angle (`tg_aoi_ladder`:
  the plate's diattenuation is 0.149 beta^2 only at small angles; give
  the number at the chosen angle, the analyzer sweep's correction, and
  the cube's).
- CCL: the arm maps V3 and the analyzer V4 on the new bench (a few
  runs), and the mask sensors' bench + battery + one loop as the gate
  that the rows hold (they did at 30 deg without substrates).
- TO: the P/SRI bench through `dmg_bench_clearance` (its pickoff and
  recombiner are at 45 deg; check the reference arm's node).

Report: one file, numbers first, the clearance table, run tags; the
deck's front-end, interferometer, vector and pinhole layout slides are
rebuilt from your figures.

## 6. The order (each item whole, committed, before the next)

1. `D_RECOMB` / `D_RC_L2` forwarded; the arm QWPs at the retro end; the
   clearance table printed by the report (`dmg_bench_clearance`); the
   lens rig re-emitted at 22.5 deg; layouts redrawn -- gate: every part
   >= +25 mm.  Commit.
2. Stage A extended to the node parts (the solve reproduces 22.5 as the
   ruled angle's neighbor; the OAP folds re-solved).  Commit.
3. Thicknesses: the 10 mm splitter and compensator, real lens edges; the
   tail retune; the mask-sensor gate run (bench + battery).  Commit.
4. Substrates on the polarizers / QWPs / analyzer / masks; the tail
   retune again; the gate run.  Commit.
5. The camera: pitch and binning in the parts lists, the sensor body in
   the layouts.  Commit.
6. The snapshot form's polarization at 22.5 deg (`tg_aoi_ladder`).
   Commit.
7. The OAP and P/SRI rigs through the tool (with TO).  Commit.
Then the report's status table says "all done" and Dave pushes.
