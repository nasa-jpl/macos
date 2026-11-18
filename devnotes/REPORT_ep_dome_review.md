# REPORT: critical look at the EP-dome sensitivity fix (CCMac 4d915f8) -- 2026-09-08

**CCMac's fix (resources dev-candidate 4d915f8, local on his Mac, not
reachable from here):** when no exit-pupil element is supplied, the four
single-DOF supervisors read the OPD at the FocalPlane (nElt) instead of
at nElt-1 whenever nElt-1 is a powered optic and the deck ends in a
FocalPlane, with a warning.  Rationale: the "FP-relation doctrine" --
the ray-grid OPD at the focal plane IS the exit-pupil-referenced
wavefront (PLAN_DESIGN_LAYER Sprint-2B note; add_pupil's docstring).
Validation offered: the segment-tilt map at the FP is "bipolar" (pos/neg
extent 1.33) where the nElt-1 read gave a one-signed dome (21.5).

## Verdict

**The diagnosis of the dome is right; the fix is wrong, and the doctrine it
rests on is wrong for exactly the DOFs a sensitivity Jacobian is built
from.**  The OPD at a FocalPlane is the optical path to each ray's
LANDING POINT on that plane.  A perfect image that merely moves (a tilt,
of the whole beam or of one segment) leaves every path equal at its own
displaced focus, so the focal-plane read is BLIND to tilt: a rigid tilt
reads as ~0, a segment tilt reads as a segment PISTON (the path
difference between the segment's displaced focus and the others'), not
as the bipolar ramp across the segment that the exit-pupil wavefront
carries.  "Bipolar" is not the test: a mean-referenced piston on one
segment is bipolar too (that segment up, the other eighteen down by a
nineteenth).  The right reference is a sphere at the exit pupil, and the
only correct behaviour on a bare-focal deck is to refuse (or to place
one), not to substitute the focal plane.

## Measurement (macos dev-candidate + handedness fix; mmacos; model 256)

Deck `templates/80_end_to_end/e2e6m/s3_imager_full.in` (26 elements,
segmented TMA + Pickoff + OAPim + Imager).  A copy with the add_pupil
pair inserted before the Imager (flat Return at the image, spherical
Return seeded at OAPim, then `stop(1)` + off-axis FEX as add_pupil does)
gives the exit-pupil read.  Segment Rx tilt alpha = 1e-6 rad; maps in
metres; the segment footprint is the piston-perturbation flat-top.

| DOF | read surface | rms | max | min | shape in the footprint |
|---|---|---|---|---|---|
| Seg1 Rx | OAPim (nElt-1, powered) | 8.0e-6 | +9.7e-6 | +6.2e-6 | one-signed dome (CCMac's symptom) |
| Seg1 Rx | FocalPlane (CCMac's fix) | 4.1e-8 | +4.14e-8 | +3.93e-8 | FLAT piston, 5% variation |
| Seg1 Rx | exit-pupil sphere | 5.9e-7 | +9.3e-7 | -1.21e-6 | bipolar RAMP (= 2 alpha x half-segment) |
| Seg8 Rx | OAPim | 7.7e-6 | +9.4e-6 | +5.9e-6 | one-signed dome |
| Seg8 Rx | FocalPlane | 2.3e-6 | +2.32e-6 | +2.29e-6 | FLAT piston, 1% variation, 4x the true rms |
| Seg8 Rx | exit-pupil sphere | 6.0e-7 | +7.6e-7 | -1.26e-6 | bipolar RAMP |
| global field tilt 1e-6 | FocalPlane | 4.4e-10 | | | blind (3000x below the pupil read) |
| global field tilt 1e-6 | exit-pupil sphere | 1.4e-6 | | | tilt (uniform-disk estimate 1.7e-6) |

So the focal-plane "tilt sensitivity" of a segment is a segment piston
whose size depends on the segment's distance from the axis (Seg1 4e-8,
Seg8 2.3e-6), with no ramp -- a column that says "piston" where the DOF
is "tilt".  For control design that is not a mis-scaled column; it is the
wrong basis vector.  The exit-pupil read is the textbook answer:
2 x 1e-6 rad x ~0.5 m half-segment = +/-1e-6 m.

## What the frozen-EP framing got right, and what CCMac's falsification missed

CCMac falsified "frozen-EP vs reset_xp" on `s2_segmented` (which HAS a
placed pupil and reads bipolar even frozen).  That shows the reset is
not the mechanism on a deck that has a pupil; it says nothing about
`s3_imager_full`, which has NO pupil element -- there is nothing for
reset_xp to reset, the guard's noPupilElt verdict fires, and the read
lands on OAPim.  The missing pupil IS the problem; the read surface is
its symptom.

## A second, pre-existing defect that made the correct path look broken

Placing a pupil on this deck with add_pupil's own recipe (`macos.stop(1)`
then FEX) obscured 732 of 985 rays -- the reason "the hand-built EP
sphere kills corner fields" and the reason the exit-pupil read looked
unusable.  Isolated: any ELEMENT stop on this deck mirrored the source
frame (xGrid -1 -> +1: UpdSrcGrid rebuilt it right-handed) while the
segment->element map did not move, so rays hit segments they were not
mapped to.  Segment stop 732 obscured, Reflector stop at OAPim 962, at
M2 814; object-space stop 2; the 253 survivors were the five x=0-column
segments.  Fixed in the engine (sourcsub.F UpdSrcGrid preserves the
deck triad's handedness): Segment stop now 2 obscured, frame intact,
e5hex1 SAVE after `stop elt 1` keeps `xGrid= -1 0 0`.  Every supervisor
harvest that passed an explicit `stop_elt` on a left-handed segmented
deck was silently losing rays to this; results on such decks are worth
a re-run.

## The handle he needed was already in the deck
`s3_imager_full` carries `SharedPupil` (element 23, a flat Reference in
the collimated space after OAP1) -- a REAL pupil plane.  A flat reference
normal to the chief ray in collimated space is a legitimate wavefront
reference, and reading the same Seg8 tilt there gives the exit-pupil
map: correlation 0.9975 with the EP-sphere read, rms ratio 1.03; the
global 1e-6 rad tilt reads 1.385e-6 there vs 1.40e-6 at the sphere.  For
a TELESCOPE-only harvest (the stated purpose of harvest_tel_sens) the
read surface is `exit_pupil_elt = 23`, no pupil placement needed; the
imaging leg's own exit pupil needs the add_pupil pair.

## Recommendation for the fix

1. Do NOT merge the FP-read branch.  Replace it: a powered nElt-1 with no
   placed pupil is an ERROR (`macos:dw_dx:noPupil`) naming
   `Telescope.add_pupil` / the FP-Return-before-ExitPupil recipe -- or,
   if the harvest is allowed to modify the deck, place the pupil
   automatically (the recipe is three engine calls) and read at it.
   Silently reading a tilt-blind surface and returning it as a Jacobian
   column is the failure the fix set out to prevent.
2. Retract the doctrine sentence in PLAN_DESIGN_LAYER / add_pupil's
   docstring: the FP OPD equals the exit-pupil wavefront only modulo the
   terms a displaced perfect image absorbs (tilt and distortion).  That
   is harmless for the optimiser's WFE objective and fatal for a
   rigid-body Jacobian.
3. The gate CCMac proposed (warning fires + map bipolar) would pass the
   wrong map.  Gate on SHAPE: within the perturbed segment's footprint
   the tilt column must be a ramp (plane fit captures >90% with a
   near-zero mean), and a global field tilt must produce a tilt of the
   expected rms at the read surface.  Both are in the committed tool
   `mmacos/tools/ep_dome_probe/` (to be promoted into tDwDx).
4. The `tPupilFindMethod` 3/7 -> 9/1 flip on his Mac is not confirmable
   from here; his mex was not built against a current engine.  Result on
   this tree is in the gate log of this session (see the final message).
5. `test_object_space_apstop_deck_needs_no_stop_elt` failing on both his
   trees is the one thing to walk with Dave next, as he asked.
