# BRIEF for TO: closing the gauge deck's lanes

CCL for TO, 2026-09-15 (Dave relays).  Branch `dev-candidate`, both repos,
shared tree on this 32 GB box.  Follows `BRIEF_to_reflective.md` (items
1-3, 5-7 done; 4 read on the seed tail) and takes over the realism round
(`BRIEF_ccmac_bench_realism.md` items 3-8, CCMac's until its budget ran
out).  Seven items, cheapest-complete-first; each lands whole and is
committed before the next.  Item 7 is the new arc Dave opened today and is
last only because it is largest; he may pull it forward.

## 0. Rules (unchanged, restated because a session can die)

- Section 0 of the realism brief in full: order cheapest-complete-first,
  commit after every item into the shared tree (`git add` your files only,
  foreign-hunk check on shared files), every run longer than minutes is a
  detached batch job chained in a `runs/<name>seq.sh`, the report written as
  you go with a status table at the TOP, small context, clean stop when the
  budget ends.
- **This box has 32 GB: ONE model-1024 MATLAB at a time.**  Every wrapper
  waits (`tg96_batch.sh` / `zwfs_batch.sh` / `pdi_batch.sh`, the lock in
  zwfs runs); never set the NOWAIT bypass here.  Kill by PID only, never by
  name (memory feedback_shell_self_kill).  Sentinels grep the wrapper's own
  `] exit` marker and are bounded.
- No `--amend`, no push (Dave pushes), no re-render of a deck figure (deck
  figures are the tools' own PNGs).  American English in every report.
- **Item 0 is real work: ~70 of your files are UNTRACKED** (`git status
  --short mmacos/templates/40_benches`): `REPORT_bench_realism.md`; the
  `runs/oap22d/*` re-render (modified, uncommitted); fourteen `*_tail.mat`
  in tg_psi_dm96_oap (tailA, oapfix, oapfixb, oapifo, oapifol, oapdesc,
  objseed2/3, objwin2/3, polA/B); run dirs tailA, tailB, oapfixb, oapifol,
  oapdraw/2, nodesolve, node22s/v, clear22, foldsmoke, lens22h; the sketch
  PNGs in lens22/lens22g/oap22/oap22d/oapdraw3/oapifo/oapifo2; oapifo2's
  vlayout + clearance PNGs; seq scripts; pdi_dm96's psriclear/psriclear2
  mat and pdi_layout.in; zwfs_dm96's oapsens22/oapsens22n .mat.  Rule:
  every tag the reports cite gets its run dir committed (grep the three
  reports for the tag); tail mats cited in section 4 (tailA/B, objseed3,
  objwin3, oapifo/oapifol) committed as evidence, the rest deleted;
  sketch PNGs deleted unless a report shows them; the oap22d re-render
  committed only if it is the auto-placed one (else `git checkout` it).
  An uncommitted artifact is an unverifiable claim.  Everything untracked
  under 40_benches is yours EXCEPT a few zwfs_dm96 leftovers from earlier
  lanes (`.batch.lock`, `mask385/macos_param.txt`, `pcam193r_perframe/`,
  `pdi_dev2/3`, `pseq6.sh`, `v3dev_oap0`, `vafter*.sh`): leave those.
  CCL's work is all committed.  **If you are a fresh session:** start with
  that `git status`, then the status table at the top of
  `tg_psi_dm96_oap/REPORT_reflective.md` (sections 4.5-4.6 and 5 are the
  live findings; 4.1-4.4 are RETRACTED attributions, do not re-derive
  them), then `runs/*.log` tails for exit codes; nothing else needs to be
  re-read.

## 0b. URGENT, AHEAD OF ITEM 1: the CTB DM model's beam diameter (Dave 2026-09-15 17:30)

Dave needs the CTB story out soon and an interim report on CTB + e2e6m +
the DM gauge by Friday 2026-09-18 (CCL writes it).  Your finding is in it,
so the FIX and the RE-SCORE come first, by **Thursday 09-17 noon**, so the
CTB deck and the report can absorb the numbers.  Your `ctb_beam_probe.m`
(untracked, 16:32) is the arbiter: commit it with its printed table as the
first line of the record.  Then:

**DAVE'S RULING (2026-09-15 18:30): the beam stays 42.75 mm; the DM
becomes a REAL device -- 64 x 64 at 1 mm pitch (the HCIT-class 64 mm DM)
-- and the controlled set is every actuator with some influence inside
the beam.**  Reach: 42.75 / (2 x 1) = 21 cycles across the beam (21 lam/D;
the 3-15 lam/D annulus sits inside it with margin).  Not chosen: enlarging
the beam to 64 mm to keep 32 cycles (a bench change), or a 0.67 mm pitch
(not a device).

1. `ctb_dm.m`, `ctb_dm_jacobian.m`, `ctb_efc_physics.m`: `nact` default
   64, `pitch_mm` default 1.0 (FIXED, no longer beam/nact), `beam_d_mm`
   default = the PROBED diameter (42.75 if the probe agrees with the deck
   header: Aperture 8.485e-3 x 2519.1 mm = 21.4 mm radius).  **The
   controlled set (Dave 2026-09-15): every actuator whose influence
   reaches inside the beam, however weakly** -- actuators beyond the
   footprint still push the edge of the pupil.  So the active rule widens
   from "centre within beam_r + 1 pitch" (influence 0.12 at the edge) to
   "centre within beam_r + WIN pitches" (WIN = 3, the stamp's own support;
   influence 0.12^9 ~ 5e-9 at the edge for the outermost ring), and the
   measurement decides their weight: poke them ALL (expect ~pi x 24.4^2 ~
   1870 per DM, ~3700 pokes, ~25 min at 512) and print the Jacobian
   column-norm ladder by distance of the centre from the beam edge (-3 ..
   +3 pitches).  EFC's Tikhonov solve gives a weak column little command
   by itself; the ladder shows where "some influence" ends on this
   influence function (Gaussian, 12% at one pitch), and that number, not
   a rule, goes in the report.
   Doc strings corrected ("gate1b probe" was the RADIUS).  Callers that
   pass their own `beam_d_mm`/`nact`/`pitch_mm` (`ctb_dst_2c`,
   `ctb_dst_s1*`, `ctb_vvc`, `ctb_study`): audit each; none may keep 21.3
   or 32.  Gate in `tests/tCtbDm.m`: the active set covers the traced
   footprint, the pitch is 1.0, the lattice spans 64 mm.
2. Regenerate: `ctb_dm_jacobian` at N=512 (~25 min) and the N=1024
   configurations the deck cites (`ctb_dm_jacobian_N1024_*` fingerprints:
   hp, ann, 2c_mono_hp, nb615-655 -- `ctb_study` derives them; read its
   config list before starting and chain the whole set in one detached
   sequence, ONE MATLAB at a time on this box, after lensuw2 exits).
   Every `.fp.json` is re-committed; the `.mat` stay gitignored.
3. Re-score deck_ctb slides 9-13 (EFC hard 2.9e-7 -> 8.1e-9; vortex loop
   1.7e-8 -> 6.8e-15; polarization residual 1.1e-15; bandwidth 8e-13 ->
   5.4e-11; vector vortex) through `ctb_study` -- the same configs, the
   new DM.  Expect the aberration-free floors to MOVE (reach 21 instead
   of 32 cycles, the 3-15 lam/D annulus still inside it; shallower or
   deeper is the measurement; the deck's slide-9 text "0.67 mm pitch =
   beam/32; 880 actuators" is rewritten from the new count).  Table in
   `CTB_PROP_STATUS.md`: old / new per slide, and the README line
   ("880 active actuators ... beam radius") corrected.  Commit per stage.
4. Hand CCL the old/new table + regenerated figures by Thursday noon;
   CCL updates deck_ctb.md and the interim report.  deck_ctb's slides
   carry a DRAFT banner until then (CCL adds it now).

Item 7's step 0 is this item; steps 1-3 of item 7 wait for it.  Item 1
(the vector rows) runs AFTER the Jacobian sequence unless the box is idle.

## 1. The vector pair on the redesigned reflective rig: rows, overcoat, verdict

Your item-5 G4 (bare Al 20/25 deg 634 pm; no coating 208; the record's
5/9 deg 19.6) is the UNCALIBRATED absolute reading, and it sits on a scale
the zwfs record already has (README V3, runs/v3s_p*): G4 absolute error
59 / 178 / 597 / 1877 pm at 0.01 / 0.03 / 0.1 / 0.3 rad rms of channel
phase difference, linear at 6 nm per rad.  So 634 pm is ~0.1 rad and 208
is ~0.035; **through the matrix on the working surface the rows HOLD to
0.1 rad** (single 0.9948 / 4 pm, grid 1.0013 / 4, dense 441 pm, ladder
0.9980 / 0.9889) and bend at 0.3 (grid 1.0068 / 8, dense 962, loop gain
-10%); the polarimetrically calibrated solver (`v_cal 'map'`) reproduces
the ideal record at 0.3 rad.  The redesign is AT the edge, and the
coating's 3x is what puts it there.  Runs, all zwfs_run on the OAP rig at
model 1024 / 193 (the oapsens22 settings):

1. **The rows:** `'stages',{'bench','battery'}`, readings V only,
   `mask.v_arm 'engine'`, `mask.v_cal 'ideal'` -> the V rows on the 30 nm
   surface (single / grid / dense / ladder).  Then `mask.v_cal 'map'`.
   Print `dmg_arm_maps`' channel phase difference rms (the bench stage
   prints it under V3) so the rung is read directly, not inferred.
2. **The overcoat at a quarter wave of 632.8 nm** on both OAPs
   (`bench.coat_oap`, MgF2 n 1.38, physical thickness 632.8/(4*1.38) =
   114.6 nm) vs bare: G4 and the rows.  The engine's measured overcoat
   rule (macos_f90/CLAUDE.md, "overcoat quarter-wave reversal"): at the
   TRUE quarter wave of the working wavelength the coating's
   cross-polarization is 0.05x of bare; the record's "protected Al" was
   not at a quarter wave of 632.8.  Expected: the coating's 3x -> ~1.1x,
   total ~230 pm ~ 0.04 rad.
3. **The loop** (`vloop`-class, same seeds as vloop193) only if the
   uncalibrated rows in 1 bend.

Deliverable: one line for the deck's redesigned-rig slide -- "vector pair
on the redesigned rig: rows X / Y pm uncalibrated, Z with the calibrated
bench, overcoat W" -- and the item-6 verdict restated: the fold lever is a
POLARIZATION lever for the vector reading only (halving OAP2's fold
roughly quarters the channel phase), pulled only if 1 bends and 2 does not
fix it.  Report: REPORT_reflective.md section 5, the status table updated.

## 2. Item 4's remainder on the seed tail: loop, descent, and the 120 nm wrap

`oapifo2` (seed tail, record resolution) reads: rows 0.9915 / 2.4 pm and
0.9905 / 161 pm, cross-talk 0.006-0.02 (the record's 0.22-0.48 was the
misfed collimator), null 28.9 pm.  What remains of item 4:

1. `oapifol` (bench + loop) and `oapdesc` (the descent, `loop.start_rms`
   ladder) on the seed tail, record resolution, tags of their own, each
   copying the seed tail under its tag (`bench.tail_from_mat`) so nothing
   falls back to another bench's tail.  The servo's photons-per-cycle for
   3 pm under noise and the 2 pm walk; the descent from 100 / 200 nm rms.
2. **The break ladder wraps at 120 nm** where the record's reflective rig
   held to 120 and the lens rig never flags.  Two things before this goes
   on a slide: (a) re-run the ladder with the unwrapper on
   (`battery.unwrap true`, the IFO captured from 150 nm with unwrap alone
   in CCMac's lens_deck runs -- compare like with like); (b) say WHY the
   seed tail wraps earlier: its magnification is 10.437 DM-mm per
   detector-mm against the lens rig's 9.879 (fewer pixels per actuator ->
   a larger phase step per pixel at the same surface), or the four-step's
   roll-off -- one number decides (the wrapped fraction vs pixels per
   actuator).
3. Then the interferometer's rows of record on the reflective rig are
   `oapifo2` + `oapifol` + `oapdesc`, and they REPLACE the record's on the
   deck's "Lenses or off-axis mirrors: the record" slide (CCL does the
   deck; you write the three lines).

## 3. The tail tuner: gate its winner, park the why

You disproved the objective fix cleanly (objwin3 reads 0.0338 and scores
better than objseed3 at 0.9809).  Do not chase the mechanism now; make the
tuner unable to hand a bench a tail that does not read:

- `tg96_tail` runs ONE single-actuator battery row through the affine on
  its winner (the quantity the battery measures: recovered gain in
  actuator space) and refuses the winner when that gain is below 0.95,
  falling back to the geometric seed with a printed reason.  Cost: one
  extra trace per tune.  Gate: a test that hands it the objwin3 tail and
  sees the refusal; the lens rig's tuned tail passes (0.9968).
- README of tg_psi_dm96_oap: the reflective rig's tail of record is the
  geometric seed (why); the tuner's objective is known to prefer
  non-reading tails on the OAP rig (the open why, with the two tags).
- The lens rig's retunes in item 4 below go through this gate.

## 4. Realism items 3-5: substrates, thicknesses, the camera

The realism brief's items 3, 3b, 4 verbatim (read them there); in this
order, each with the mask-sensor gate run (bench + battery, S / V / P on
the lens rig at 22.5, the gate22_193 settings) after it:

1. The 10 mm splitter and compensator (the transmitted beam's 1.4 mm
   shift carried by the builder's chief path), real lens edges (3-5 mm on
   the 103 mm singlets); the lens rig's tail retune THROUGH item 3's gate;
   gate run.  Commit.
2. `'substrate', [n t]` on `add_polarizer`, `add_waveplate` and the mask
   seat: two refracting faces (fused silica, 2-3 mm) around the ideal
   element; the mask plate's faces BEFORE the sandwich's entrance sphere
   (W040 19 nm + 0.7 mm focus shift for 2 mm at NA 0.12 -- the retune
   absorbs the shift; print the W040 the trace carries); the vector QWP in
   the F/3.6 pupil-image leg (34 nm for 2 mm).  Retune through the gate;
   gate run; `tBench` green.  Commit.
3. The camera: pitch and binning in every parts list (the 6.5 um sCMOS
   binned 4 at the 9.4 mm pupil image = the modeled 385; the polarization
   camera's 3.45 um / 680 per orientation for the snapshot form), the
   sensor body drawn at size in every layout (lens, OAP, vector, pinhole,
   P/SRI).  Commit.

The analyzer (3b) is decided: a polarization camera; the deck says so.
Report: `REPORT_bench_realism.md` (commit it first -- item 0), status
table on top.

## 5. Realism item 6: the snapshot form's polarization at the built angles

`tg_aoi_ladder` at 22.5 deg on the lens rig (the plate's diattenuation is
0.149 beta^2 only at small angles: give the number at 22.5, the analyzer
sweep's correction, the cube's) and on the reflective rig at ITS angles
(20 / 25 deg OAP folds + the 22.5 plate).  One table each: the
uncorrected arm rotation, the corrected, the residual through the matrix.
Report section + commit.

## 6. Realism item 8: the interferometer's station-by-station figure

The sensors have theirs (`zwfs_run`'s `stations_fig_`, in the deck: S / V
/ P at 22.5, and yours on the OAP rig in oapsens22).  Make the same for
the interferometer in `tg96_run`, both rigs: two rows (the flat, the 30 nm
working surface) x seven stations -- the mirror command; the test-arm
pupil field at the detector (intensity); the reference-arm field; two of
the four phase-stepped frames (0 and pi/2); the recovered phase map as
surface (nm); the raw map minus the engine's field (pm).  One figure
`<tag>_stations.png`, 1800 px wide, 10-13 pt, the pupil cropped to its
box, Interpreter none.  Tags: the lens rig's `lens_deck` settings and the
reflective `oapifo2`.  Commit; CCL puts them on the deck.

## 7. The coronagraph field servo (Dave 2026-09-15; NOTES_gauge_in_coronagraph.md)

Read the note first (macos root): a pupil dichroic at the apodizer pulls
out-of-band light to a vector Zernike or pinhole reading of the
coronagraph's INPUT field, amplitude and phase, and the two DMs servo that
field to the one recorded when the dark hole was dug -- holding the hole
against everything upstream of the pickoff (telescope misalignment,
segment and mirror figure drift), not only the DMs.  The note's section
C+ has the four bounds (what the pickoff cannot see; DM2's Fresnel
amplitude authority, 530x stroke at 2 cycles across the beam, 8x at the
actuator Nyquist on the CTB; the out-of-band transfer; stellar photons:
10 pm per 200 s on a V=5 star from the deck's servo cost).  Model it on
the CTB deck (`30_instruments/bench_ctb`, `ctb_dcr.in`, the DM model
`ctb_dm.m`, the EFC chain `ctb_efc.m`) in three steps, each its own run
and report section, new runner `ctb_field_servo.m` + `REPORT_field_servo.md`
in bench_ctb, sheet-driven like every runner:

1. **The reading at the apodizer conjugate.**  Both DMs as grid surfaces
   (ctb_dm), a dichroic pickoff after OAP3 (a Reference plane for now; the
   plate comes with item 4's substrate option), the vector pair with its
   clear frame (`mask.v_clear`) at a focus + pupil reimage behind the
   pickoff (the zwfs_dm96 sandwich, `dmg_zwfs_gauge`), out of band
   (632.8 nm against the CTB's science band).  Gate: the reading
   reproduces the engine's complex field at the apodizer to the V5 record
   (0.016 / 0.128 pm through 5 / 20% amplitude dips).
0. **First, the DM model's beam diameter** (your finding, note section
   6): `beam_d_mm` from the traced footprint (42.75 mm), pitch 1.336 mm,
   the N=512 Jacobian fingerprint regenerated, hard-occulter + vortex EFC
   re-scored, README corrected.  Commit before step 1.
2. **Separability of the two DMs.**  The multiplexed matrix
   (`battery.calib_mode 'matrix'`) over BOTH DMs' actuators as columns;
   the DM1-vs-DM2 column cross-talk vs spatial frequency.  On the real
   1.34 mm pitch the note predicts 12% amplitude conversion at the
   actuator Nyquist (16 cycles) and 50% only at 33, so expect the two DMs
   NOT to separate inside the band from one conjugate; the measurement is
   how much cross-talk the matrix carries per frequency, and the
   second-conjugate reading (DM2 reimaged behind the same pickoff) is the
   baseline for step 3, not a variant.  Rows for DM1 pokes, DM2 pokes, and
   both, at one and at two conjugates.
3. **The servo.**  `dmg_loop` with the input field as the state and both
   DMs as actuators; a drift injected UPSTREAM (a 10 urad tilt of OAP1,
   a 1 nm figure walk on OAP2 -- `run_simulator`'s pattern), the hold at
   gain 0.5, photons per cycle from the note's stellar table; the
   dark-hole contrast scored through `ctb_efc`'s chain before the drift,
   after it, and after the hold.  The regularization's placement follows
   the stroke table.  One figure: contrast vs cycle for drift alone and
   drift + hold.

Deliverable: the three numbers Dave will ask for -- the reading's
accuracy at the apodizer, the frequency above which the DMs separate,
and the contrast the hold recovers per cycle at the stellar photon rate.

## Order and the status table

| item | what | report |
|---|---|---|
| 0 | commit the untracked record files | -- |
| 1 | vector pair on the redesign: rows, overcoat, verdict | REPORT_reflective 5 |
| 2 | item 4's loop + descent + the 120 nm wrap explained | REPORT_reflective 4 |
| 3 | tail tuner gated by a battery row; README | REPORT_reflective 4.7 |
| 4 | realism 3-5: thicknesses, substrates, camera | REPORT_bench_realism |
| 5 | realism 6: snapshot polarization at the built angles | REPORT_bench_realism |
| 6 | realism 8: the interferometer's station figure, both rigs | REPORT_bench_realism |
| 7 | the coronagraph field servo, three steps | bench_ctb/REPORT_field_servo |

Put this table at the top of REPORT_reflective.md as the live status (the
realism and field-servo reports link back to it), and mark each item done
/ running (tag, started, expected end) / not started as you go.  CCL folds
each landed item into the deck; Dave pushes.

## CCL reply on the OAP ladder (2026-09-15 evening)

Agreed on the statistic and the mechanism: max|h|/(lambda/4) saturates the
moment one pixel wraps; the beyond-fold FRACTION of the lit pupil is the
meter (the ZWFS battery's fold0/fold), and differencing two wrapped absolute
maps is wrong by lambda/2 at exactly the pixels a poke pushes across the
boundary.  That is the same finding that produced `dmg_unwrap` and
`battery.unwrap` on 2026-09-13, so the fix is not new machinery: the tg96
ladder reports the beyond-fold fraction, and the CAPTURE numbers of record
for both rigs are the unwrapped ones -- the deck's interferometer capture row
("100-600 nm of wavefront to 2 pm with unwrapping alone", lens_deck) already
is.  Land the meter edit in tg96_run when the in-flight runs release it.

One thing the control has to answer before "smaller capture range" comes
off the deck: the record's lens ladder (raw, no unwrap) HELD at 120 nm
(1.0258 / 19.3 pm) and broke at 240 (1.93), while the OAP rig broke at 120.
If the erfc of the commanded surface were the whole story, both rigs would
break at the same rung; they did not.  So read lensuw2 by its MEASURED
beyond-fold fraction at 120 nm against the OAP rig's, not against the
19% prediction.  If the lens rig's measured fraction at 120 is also ~19%
and it still reads, the arithmetic is not the whole mechanism; if its
fraction is lower (a tilt or piston term in the OAP rig's difference map
would add to the wrapped fraction at the same commanded rms, and so would
the 10.44 vs 9.88 magnification through the pixel gradient), the number
that differs is the finding.  Either way: two lines in the report -- the
measured fraction per rig at 60 / 120 / 240, and the unwrapped ladder for
both -- and the deck's slide 3 and redesigned-rig lines change to whatever
those say (CCL edits; you write the two lines).
