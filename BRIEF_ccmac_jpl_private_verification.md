# BRIEF for CCMac: JPL-private verification pass on IRIS and OPTIIX (cold start)

From CCL for Dave, 2026-09-08.  Supersedes `BRIEF_ccmac_sens_core.md`.
Written for a fresh session: everything you need is here or pointed to;
nothing depends on remembering an earlier conversation.

## 0. Purpose, in one paragraph

Three things changed today that no public deck exercises the way the
JPL decks do, and Dave wants evidence on IRIS and OPTIIX (he provides the
prescriptions; they never leave the JPL side) before he merges
`sens-core` to `dev-candidate` and files the MR to `dev`:
(1) the four sensitivity supervisors now share one core (`sens-core`);
(2) the engine's exit-pupil finder changed definition (FEX is
frame-independent: the medial pupil, not the tangential one), and an
element STOP no longer mirrors the source frame on left-handed decks;
(3) the supervisors REFUSE to read a wavefront anywhere but a pupil.
Your job is to run each JPL deck through three tests, classify every
difference as expected-by-design or STOP-AND-REPORT, and send Dave one
table.  You adjudicate nothing; a difference outside the expectation
table is a finding, not physics to explain.

## 1. What changed since your last look (SHAs you will pin against)

macos `dev-candidate` (engine tip `44fc362`; later commits are records):
- `eb84095` FEX's Rx-order warning removed (cosmetic); the CLI STOP
  accepts multi-value answers on one line and accepts Segment elements;
  `stop_info_set` (so `macos.stop`) accepts Segments too.
- `82d8148` FEX/SXP frame-independent: four probes (+/-5e-6 rad about two
  axes perpendicular to the chief), mean crossing = the MEDIAL pupil.
  The legacy single probe about xGrid was the TANGENTIAL pupil.  On any
  off-axis deck the EP radius moves by half the tangential/sagittal
  split; symmetric on-axis decks are bit-identical.  Every FEX run now
  prints `axis1 .. axis2 .. medial ..` and `T/S split .. +/- asymmetry`.
- `44fc362` element STOP preserves the source frame's handedness.
  Before: `macos.stop(<elt>)` on a left-handed deck (`xGrid = -1 0 0`
  with ChfRayDir +z -- e5hex1, 6MST, every Telescope emission) flipped
  xGrid, and on a SEGMENTED source that mirrored the ray grid against the
  segment map: 732 of 985 rays obscured on e2e6m s3 after `stop(1)`, 962
  after a Reflector stop, 2 after an object-space stop.  Now: no flip.
resources `dev-candidate` `3034fff`: your `wf_elt_auto` error variant with
the element-TYPE criterion (Return/Reference at nElt-1 or refuse), your
`tEpDomeGate`, the five FEX-definition pins re-pinned to medial values.
resources `sens-core` `35274ac` = all of the above + the shared core
`private/dw_multi_core.m` + its preflight (`macos:<front>_multi:noPupil`
before the field loop; `flatPupil` warning on a flat reference) + the
three tDwDx tests rebuilt on committed fixtures.  tDwDx 25/25 +
tEpDomeGate 5/5 + the six sibling classes green there.

Records if you need the why: `macos/REPORT_fex_probe_frame_independent.md`
(45-deck pre/post table), `macos/REPORT_ep_dome_review.md` (the pupil
ruling, measured), `macos/BRIEF_ccmac_ep_dome_round2.md` (rounds 2-3),
`macos/macos_f90/CLAUDE.md` sections "FEX probe is FRAME-INDEPENDENT",
"Element STOP preserves the source frame's HANDEDNESS", "STOP on a
Segment element".

## 2. Setup on the JPL Mac (two engines, three resources trees)

Engines (gfortran release; `source ./makems.sh release gfortran` in each
checkout; each gets its own `build_release_gfortran/`):
- OLD = macos `82ced2b` in a worktree: `git worktree add ~/dev/macos_old
  82ced2b`.  This already contains the NS flow-of-light fix (`6bab7af`),
  so ONLY today's changes sit on the engine axis.
- NEW = macos `dev-candidate` at `44fc362` or later (the tip is fine).

Resources trees (one mex each, relinked against the engine named):
- PRE = `git worktree add ~/dev/res_pre eda6c5a` -- the last commit
  before the supervisor rearchitecture.  Mex: NEW engine.
- POST = `sens-core` at `35274ac`.  Mex: NEW engine.
- POST_OLDENG = `git worktree add ~/dev/res_post_old sens-core`.  Mex:
  OLD engine.  (Only for Test 2.)
Mex relink: in each tree's `mmacos/`, `make FC=gfortran
MACOS_BUILD_DIR=<that engine's build_release_gfortran>` (unset any
inherited FC first; the Makefile rejects f77).  Confirm with
`ls -la src/mmacos.mexmaci64` per tree before running anything.

Discipline: ONE `matlab -batch` per invocation, one MATLAB at a time
(state leaks across model sizes inside a process).  Model size 256 by
default; raise to 512 for big-grid decks (mGridMat caps at 256 regardless
of model size).  If a deck loads slowly, that is normal for IRIS.

## 3. Deck intake card (per deck, BEFORE any A/B; NEW engine, POST tree)

Run once and record -- this card decides which expectations apply:
```matlab
macos.init(256); nE = macos.load_rx(rx);
i1 = macos.get_elt_info(nE-1); i0 = macos.get_elt_info(nE);
c  = macos.get_src_csys();  lh = dot(c.xDir, cross(c.yDir, c.zDir)) < 0;
t  = macos.trace(); rs = macos.get_ray_status(t.nRays);
fprintf('nElt %d | nElt-1 %s (kr %.4g) | nElt %s | left-handed %d | rays %d obscured %d lost %d\n', ...
  nE, i1.type, macos.get_elt_kr(nE-1), i0.type, lh, t.nRays, nnz(rs.status==1), nnz(rs.status>=2));
f = macos.fex(1);   % prints axis1 / axis2 / medial and the T/S split
```
Card fields: nElt; type at nElt-1 (**Return/Reference = pupil deck;
anything else = bare-focal deck**); stop declaration in the header
(`ApStop= x y z` object-space, or an element stop); left-handed frame
(yes/no); segmented source (`nSeg` in the header); NS elements present;
nominal ray survival; the FEX print (axis1, axis2, medial, T/S split).
Then pick the field half-widths: start at `fx = fy = 1e-5` rad and
confirm on the PRE tree that no field is fully vignetted (its
`dw_dx_multi` still hard-errors on that; sens-core only warns).

## 4. Test 1 -- the supervisor axis (the merge gate)

PRE vs POST, SAME NEW-engine mex on both sides, per deck:
```bash
matlab -batch "addpath('<POST>/mmacos/tools/sens_core_ab'); sens_core_ab('<PRE>/mmacos/mmacos_setup.m',  '<deck>', 'out_pre/<deck>',  1e-5, 1e-5)"
matlab -batch "addpath('<POST>/mmacos/tools/sens_core_ab'); sens_core_ab('<POST>/mmacos/mmacos_setup.m', '<deck>', 'out_post/<deck>', 1e-5, 1e-5)"
matlab -batch "addpath('<POST>/mmacos/tools/sens_core_ab'); sens_core_ab_compare('out_pre/<deck>', 'out_post/<deck>')"
```
Expectation table:

| deck card | PRE | POST | verdict |
|---|---|---|---|
| pupil deck (Return/Reference at nElt-1) | runs | runs | must be BYTE-IDENTICAL or error-parity on every set |
| bare-focal deck | runs (reads the powered optic: the dome) | `macos:<front>_multi:noPupil` | EXPECTED -- record "bare-focal"; go to Test 3 |
| flat Return/Reference at nElt-1 | runs | runs + `flatPupil` warning on stdout | byte-identical structs expected; note the warning |
| field fully vignetted | hard error `emptyOPD` | warning, completes with 0 rows | EXPECTED, but you chose fx/fy too wide -- shrink and rerun |

Anything else that differs: STOP-AND-REPORT with deck, set, field,
max|diff|.  Do not adjudicate.  (Dave 2026-09-07: on segment-class NS
decks the supervisor axis is byte-identical by expectation.)

## 5. Test 2 -- the engine axis (OLD vs NEW engine, same POST supervisor)

Per deck, two cheap probes on each engine (POST_OLDENG vs POST), one
MATLAB each, then one sens_core_ab pair if the deck is a pupil deck:

a. Plain trace + ray status at the nominal field AFTER the stop the
   supervisors would use (object-space `stop_obj` if the header says so,
   else `macos.stop(<stop elt>)`): rays, obscured, lost.
   Expectation: identical, EXCEPT a left-handed deck harvested with an
   ELEMENT stop, where OLD obscures many rays and NEW does not (the
   handedness fix).  Report both counts.  Object-space stops: identical.
b. `macos.fex(1)` radius and vertex on both engines.  Expectation:
   NEW - OLD = half the T/S split the NEW print shows, i.e.
   |dR| ~ |axis1 - axis2| / 2, along the chief ray; on-axis symmetric
   decks 0.  A |dR| far from that number is STOP-AND-REPORT.
c. Pupil decks only: `sens_core_ab` on POST_OLDENG vs POST and compare.
   Expectation: the Jacobians DIFFER (the reference sphere moved along
   the chief by dR); report max relative column change per family.
   That number is the FEX-definition effect on this deck -- Dave wants
   it, it is not a finding.  Only a difference on a deck whose FEX did
   not move (b. gave 0) is STOP-AND-REPORT.

## 6. Test 3 -- the pupil-read audit (Dave's ruling: read at the pupil)

For every deck, from the intake card:
- Pupil deck: nothing to do beyond Test 1.  If the pupil is FLAT, check
  it sits in collimated space (marginal-ray directions parallel to the
  chief at nElt-1); if not, say so -- a flat reference in converging
  space is tilt-blind like the FocalPlane.
- Bare-focal deck: the supervisors refuse by design.  Place a pupil and
  rerun Test 1 on the pupiled copy:
  `python3 <POST>/mmacos/tools/ep_dome_probe/make_pupil_deck.py <deck> <deck>_pupil.in`
  then load, set the deck's stop, trace, `macos.fex(1)`, `macos.save_rx`
  -- see `<POST>/mmacos/tools/ep_dome_probe/dome_probe.m` for the exact
  call sequence.  If the deck carries a Reference at a collimated pupil
  plane (a "shared pupil" element), passing `'exit_pupil_elt', <that
  element>` is the alternative; record which you used.
- One shape check per pupiled deck (any single segment/optic tilt DOF):
  the column must be a bipolar ramp over that element's footprint, not a
  flat piston.  `dome_probe.m` shows the plane-fit metric; a min/max of
  opposite sign over the footprint is the two-line version.

## 7. What to ask Luis for

Decks with strongly curved NS optics or gratings, and any historically
odd deck -- our corpus is segment-class NS (aperture-partitioned segments
on a shallow common base, one crossing per ray), which is exactly where
"byte-identical" is expected and therefore proves the least.

## 8. Deliverable

One report to Dave: per deck, the intake card, Test 1 verdict (identical
/ error-parity / bare-focal-refused / STOP-AND-REPORT), Test 2 numbers
(ray counts old/new, FEX dR vs half-split, max relative Jacobian change),
Test 3 read surface used, notes.  Attach the `sens_core_ab_compare`
logs.  Deck names and numbers only -- no prescription content leaves the
JPL side.  On green (all differences inside the tables above), Dave
merges `sens-core` -> `dev-candidate`, deletes the branch, and files the
MR to `dev`.  Budget: do ONE IRIS deck end-to-end first (intake, Test 1,
Test 2, Test 3) and send that before batching the rest, so the protocol
gets corrected on the cheapest deck.

## 9. Pointers

`mmacos/tools/sens_core_ab/README.md` (harness rules),
`mmacos/doc/SENSITIVITY_TOOLS.md` (the stack map),
`macos/BRIEF_luis_round3.md` (why the supervisor arc exists),
`macos/REPORT_fex_probe_frame_independent.md`,
`macos/REPORT_ep_dome_review.md`, `macos/BRIEF_ccmac_ep_dome_round2.md`.
