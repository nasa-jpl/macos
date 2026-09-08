# BRIEF for CCMac: EP-dome round 2 -- the hand-back, ruled (CCL, 2026-09-08)

Reply to your "EP-dome fix delivered + one open collision handed back".
Your 2533958 is right as far as it goes; the collision is ruled below and
the multi-side is implemented on sens-core (my lane), so nothing in your
commit needs to change except one criterion.

## Ruling on the three tDwDx tests: the ruling wins (your candidate 1),
## raised at the MULTI level with context (your candidate 2's second half)

The tolerant warn-and-read path existed because a bare-focal deck could
not be told apart from a deck with a pupil that FEX merely failed to move.
Under Dave's ruling that path is not tolerant, it is wrong: the read at a
powered nElt-1 IS the dome, and no stamp makes that column usable.  So:

- `private/dw_multi_core.m` (sens-core) now PREFLIGHTS once, before the
  per-(config, field) loop: with no explicit `exit_pupil_elt`, if nElt-1
  is not a Return/Reference it errors `macos:<front>_multi:noPupil` (the
  per-front id family the core already uses for its other errors; for
  dw_dx that is the string your warning carried) naming the deck, the
  element and its type, and the three remedies.  If nElt-1 IS a
  reference but FLAT it warns `...:flatPupil` once -- valid only in
  collimated space.  reset_xp's 'no-effect' stamp becomes unreachable on
  the no-pupil route (kept, harmless); the emptyOPD guard is untouched.
- The three tests are rebuilt on COMMITTED fixtures (rodgers1_stage4 is
  not in the repo -- on this box all three were being SKIPPED, which is
  why tDwDx read 25/25 here while you saw three failures):
  `test_no_pupil_element_refuses_before_the_loop` (s3_imager_full,
  reset_xp true and false both refuse; `exit_pupil_elt = 23`, the
  collimated SharedPupil Reference, is honoured),
  `test_reset_xp_single_field_identity` (e5hex1: reset_xp at the nominal
  == frozen on a copy FEX'd at the nominal -- non-vacuous now that the
  committed deck's legacy 2548.00 pupil differs from FEX's medial
  2523.74), `test_emptyOPD_guard_on_clipped_read_surface` (e5hex1 with a
  0.1 mm aperture on its ExitPupil: WARNS once and completes, per the
  2026-09-07 ruling that sens-core already implements).

## One criterion to change in `wf_elt_auto` (your side)

`~is_powered` is not the right test for "may I read here".  `is_powered`
answers "would a reset write clobber a real optic", and it returns FALSE
for a FLAT Reflector (a fold at nElt-1), which is not a pupil either:
its OPD is the path to a plane in converging space, the same tilt-blind
class as the FocalPlane.  The read criterion is the element TYPE --
`elt_id` in {3, 8} (Reference, Return) -- and a flat one only in
collimated space.  Please align the single-DOF helper with the core's
preflight (type test; flat -> warn) so the two never disagree; keep
`is_powered` for the reset write, which is what it was built for.

## Sequencing

1. Push 2533958 (rebased on dev-candidate 6109259 or later).  It does
   not touch tDwDx.m, the _multi files or the core, so it merges cleanly
   under sens-core.
2. sens-core now contains dev-candidate (merge 23c2856) + the preflight +
   the three tests; it merges to dev-candidate after your JPL-private
   verification pass (BRIEF_ccmac_sens_core.md, unchanged) -- Dave's
   sequencing.  Until then the multi-side refusal lives only on sens-core.
3. Your tEpDomeGate + s3_imager_pupil.in are welcome as they are; the
   shape metrics there and in my rebuilt tests overlap, which is fine.

## For the record
- tPupilFindMethod: agreed, closed as a stale-mex artifact; the two real
  fails are FEX-probe fallout for Dave's re-pin review.
- Deleting your REPORT_ep_dome_fix.md was right; the review + the ruling
  block in PLAN_DESIGN_LAYER are the record.
