# BRIEF for TO (Opus): the reflective front end, designed from scratch on the 22.5-deg bench

From CCL for Dave, 2026-09-15.  The deck's "lenses or off-axis mirrors"
slide is a table with no credible drawing: the reflective rig's layout
figure (CCMac's `oap_vlayout.png`) shows the two off-axis parabolas away
from the traced rays, and its numbers were taken on the 7-deg bench that
is not buildable.  Dave: describe the reflective version from scratch.
Your lane (the PDI) is complete; this is your next.  Standing rules as
before (dev-candidate; one model-1024 MATLAB at a time on the Linux box,
two on a 64 GB box; every number in a committed report; American
English).  **Work in increments that survive a dead session** -- section
0 of `BRIEF_ccmac_bench_realism.md`, verbatim: cheapest-complete-first,
commit per item, detached runs, the report written as you go with a
status table on top, a clean stop when the budget ends.  CCMac's session
died at its first item yesterday; nothing was lost because nothing had
been started -- keep it that way by landing whole items.

## 1. The requirement

The same bench as the lens rig, with the two lenses replaced by two
off-axis parabola sections, everything else unchanged:
- source: the filtered HeNe, 632.8 nm; a 103 mm collimated beam (the DM
  aperture 96 mm plus margin: `R_TO_AP` 51.4 mm);
- collimator OAP1 (replaces L1, f 857 mm) and focuser OAP2 (replaces L2,
  f 429 mm): the internal focus at F/4.2 at the mask seat, the tail
  (mask seat, field lens, camera at the pupil image) re-tuned for it;
- the splitter node at **22.5 deg** (Dave's ruling; `zwfs_params` /
  `tg96_params` carry it: `BS_AOI` 22.5, `D_RECOMB` 150, `D_RC_L2` 55,
  compensator at 200 mm), the 700 mm DM leg, the reference arm on its
  piezo -- the builder's `'oap'` mode keeps the splitter, both arms and
  the recomb plane geometrically as the lens rig and folds only the
  source -> OAP1 and OAP2 -> detector legs, in the splitter's plane;
- **buildable**: every part clears every beam it is not in by >= 25 mm
  (`dm_gauge_lib/dmg_bench_clearance`, the parts' aperture + an 8 mm
  mount), INCLUDING the two OAP bodies, the source and the camera on
  their folded legs; the OAPs' off-axis distance and fold angle come out
  of that, not the other way round;
- coatings: bare and protected aluminum (`bench.coat_oap`; CCMac's
  stacks in `zwfs_params`).

## 2. What is known (CCMac's 7-deg work, `tg_psi_dm96_oap/REPORT_oap.md`
##    and `REPORT_gauge_ifo.md` section 4; reproduce or supersede)

- Fold angles from the Stage-A fold solve: OAP1 5 deg (off-axis 149 mm),
  OAP2 9 deg (132 mm); astigmatism scales as the fold angle squared,
  halving OAP2's angle lengthens its leg 1.33x.
- The mask seat trim solved on the trace: 6.14 mm from the lens default
  (the collimator's residual defocus refocused by OAP2, which is fed
  on-axis); best-focus blur 1.1 lamF/D, fold coma 0.82 lamF/D at 9 deg
  (0.17 / 0.31 / 0.47 / 0.65 / 0.82 at 1 / 3 / 5 / 7 / 9 deg).
- The interferometer on it: single 0.996 / 2 pm; 52-site 0.958; dense
  0.950 / 2.1 nm; cross-talk 0.18 (lens < 0.06); servo 3.4e13 noise-only
  and NEVER under the 2 pm walk (4.1 pm floor); noiseless step 27.6%;
  does not capture (stalls at 5.9 / 24 / 53 nm from 60 / 150 / 300).
- The mask sensors on it (`zoap`): the stepped dimple survives (1.000 /
  4 pm; 47 sites 0.967; capture 38 nm), the vector pair fails its gate
  (19.6 pm on 100 nm pokes), the pinhole fails (94 pm): the more
  focus-critical the mask feature, the worse the fold coma.
- Alignment: 10 urad of OAP tilt = 95-103 nm of null shift; 10 um of
  decenter = 15-16 nm.

## 3. Why the earlier drawing is wrong, and the first thing to fix

In `oap_vlayout.png` the OAP bodies sit off the traced rays.  Hypothesis
to test first: an off-axis parabola element's vertex (`VptElt`, what
`view_rx` draws the body at) is the PARENT parabola's vertex, 132-149 mm
from the section the beam actually uses, so the viewer places the body
at the vertex and the rays through the section.  If so, the fix is in
`view_rx` (`elt_geom_`: for a conic whose ray footprint is far from its
vertex, draw the body at the footprint's centroid with the local sag) --
gate it on the lens rig being unchanged.  If the rays are the problem
(the OAP not folding as built), that is a builder bug and it comes
first.  Do not put a drawing in the deck until the bodies sit on the
beam.

## 4. Deliverables, in this order (each whole, committed)

1. **The diagnosis of section 3** and its fix (viewer or builder), with
   the lens rig's figure bit-unchanged as the gate.  Commit.
2. **The design**: fold angles and off-axis distances from a clearance
   solve that includes the node parts and the OAP bodies (extend Stage
   A in `tg96_run`, or solve with `dmg_bench_clearance` in a loop:
   smallest fold angles that clear everything by 25 mm); the seat trim
   re-solved on the trace; the tail re-tuned (`tg96_tail` for 'oap').
   The clearance table printed in the report.  Commit.
3. **The layout**, in the recipe (`zwfs_vlayout.m`): the train from
   above with the OAP sections drawn where the beam hits them, the node
   panel, the tail panel; the parts list (OAP1 / OAP2: parent f,
   off-axis distance, fold angle, section size, coating; the folds'
   mounts).  Commit.
4. **The interferometer on it**: the rows on the 30 nm surface (matrix
   on it), the servo, the descent -- `tg96_run('bench.optics','oap',
   ...)`, the same stages as the lens record.  Commit.
5. **The mask sensors on it**: `zwfs_run('bench.optics','oap',
   'bench.coat_oap','bareAl', 'readings',{'S','V','P'}, 'stages',
   {'bench','battery'})` with `mask.v_arm 'engine'`; if a gate fails
   (G4 / G5), report the number and stop, as before.  The stations
   figure (`stations_fig_`) for each reading on the OAP rig.  Commit.
6. **The fold-angle lever**: the same at half the OAP2 angle (the longer
   leg): does the pinhole recover?  One table.  Commit.
7. **Your own bench through the clearance tool**: the P/SRI's pickoff
   and recombiner at 45 deg, its reference arm's node, on the 22.5-deg
   front end.  Commit.

## 5. Report

`tg_psi_dm96_oap/REPORT_reflective.md` (a new file; CCMac's REPORT_oap.md
stays as the 7-deg history): numbers first, the status table on top,
the clearance table, run tags, figures.  The deck's "lenses or off-axis
mirrors" slide and a new "reflective rig: layout" slide are rebuilt from
it; CCL assembles.

## 6. After this brief: the rest of the bench-realism round is yours

CCMac will not return before this work concludes (Dave 2026-09-15).
`BRIEF_ccmac_bench_realism.md` items 3-8 -- real thicknesses (10 mm
splitter and compensator, lens edges), substrates on the plates and
masks, the camera's pitch and body, the snapshot form's polarization at
22.5 deg, the interferometer's station-by-station figure -- are yours
after the reflective front end, in that order, each whole and committed.
Items 1 and 2 are done (resources 9b2181c, 1338920).

## CCL reply to item 5's first result (2026-09-15 15:30) -- the vector pair is not lost yet

TO's G4 on the redesigned rig: bare Al at 20/25 deg 634 pm, no coating 208 pm,
the record's 5/9 deg 19.6 pm; capture and the scalar readings untouched, so
a polarization-only effect.  Right, and it has a place on a scale that
already exists: **G4 is the UNCALIBRATED absolute reading** (`v_cal 'ideal'`),
and the V3 design scan in `zwfs_dm96/README.md` (runs/v3s_p*) priced exactly
this term -- G4 absolute error 59 / 178 / 597 / 1877 pm at 0.01 / 0.03 /
0.1 / 0.3 rad rms of channel phase difference (linear, 6 nm per rad).  So:

| rig | G4 | = channel phase difference | where that sits on the V3 scan |
|---|---|---|---|
| record 5/9 deg, bare Al | 19.6 pm | ~3 mrad | below the scan's first rung |
| redesign 20/25 deg, no coating | 208 pm | ~0.035 rad | between the 0.03 and 0.1 rungs |
| redesign 20/25 deg, bare Al | 634 pm | ~0.1 rad | AT the 0.1 rung |

and the scan's verdict at those rungs (all measured, README V3): **through
the matrix on the working surface the rows HOLD to 0.1 rad** (single 0.9948
/ 4 pm, grid 1.0013 / 4, dense 441 pm, ladder 0.9980 / 0.9889) and bend at
0.3 (grid 1.0068 / 8, dense 962, loop gain -10%); the polarimetrically
calibrated solver (`v_cal 'map'`) reproduces the ideal record to the digit at
0.3 rad (0.057 pm).  The bare-Al redesign therefore sits at the edge of
"nothing reaches the matrix", and the coating's 3x is the part that puts it
there.  Three runs decide it, in this order:

1. **The rows, not the gate:** the V reading's battery on the 30 nm surface
   (`'stages',{'bench','battery'}`, `mask.v_arm 'engine'`, `mask.v_cal 'ideal'`)
   on the redesign with bare Al -- if single/grid/dense hold within the
   scan's 0.1-rad numbers, the vector pair survives the buildable folds and
   the deck's recommendation stands unchanged.  Then the same with
   `mask.v_cal 'map'` (the calibrated bench): expected ideal.
2. **The overcoat at a quarter wave of 632.8 nm** on both OAPs
   (`bench.coat_oap`; physical thickness ~114.6 nm MgF2 at 632.8 -- see the
   engine's overcoat-reversal note: at the true quarter wave the coating's
   cross-polarization is 0.05x of bare; the record's "protected Al" was NOT
   at a quarter wave of the working wavelength).  That should take the
   coating's 3x to ~1.1x and the total to ~230 pm ~ 0.04 rad, comfortably
   inside.  Print `dmg_arm_maps`' channel phase difference rms so the rung
   is read directly instead of inferred from G4.
3. **The loop** (`vloop` on the redesign) only if 1 bends.

Item 6's lever is then a polarization lever, not a blur one: if 1 bends
at bare Al and 2 does not fix it, halving OAP2's fold (25 -> 12.5 deg)
roughly quarters the channel phase (diattenuation ~ AOI^2) at the price of
the clearance solve -- the trade TO's item 6 reframed away for blur
comes back for the vector reading alone.  The pinhole and the stepped
scalar dimple are unaffected either way (TO's own numbers).

Process notes acknowledged: the sentinel must grep the batch wrapper's
own exit marker (memory feedback_shell_self_kill: bound the loop AND
match the specific marker); the runner's second battery (48x48) is in
every zwfs run's timing.
