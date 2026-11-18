# BRIEF for CCMac: EP-dome fix -- Dave's ruling and the replacement (2026-09-08)

Read first: `macos/REPORT_ep_dome_review.md` (the measurement),
`mmacos/tools/ep_dome_probe/` (the probe, runnable), the resolution block in
`macos/PLAN_DESIGN_LAYER.md` under "Optimization target = exit-pupil WFE".

## Ruling (Dave, 2026-09-08)

**The wavefront is read at the PUPIL by default.  That is what sets the OPD
for a PSF calculation, or any other diffraction calculation.**  Reading it at
the FocalPlane is not an alternative reference: the focal-plane OPD is the
path to each ray's landing point and is blind to tilt (global 1e-6 rad
field tilt: 4.4e-10 m rms at the FP vs 1.4e-6 at the exit-pupil sphere;
a 1e-6 rad segment tilt: a flat 2.3e-6 piston on the segment at the FP vs
the +/-1e-6 bipolar ramp at the sphere).  The "FP-relation doctrine" is
retracted; it held only modulo the terms a displaced perfect image absorbs.

Consequence for 4d915f8: **do not merge the FP-read branch of
`wf_elt_auto`.**  The diagnosis (dome = read at a powered nElt-1) stands;
the powered test via `reset_xp_guard('is_powered')` stands; the
delegation argument for the _multi drivers stands.  Only the fallback
changes.

## Replacement decision tree for `wf_elt_auto`

1. `exit_pupil_elt >= 0`: honoured verbatim (unchanged).
2. Auto-select, in this order:
   a. nElt-1 is a Return/Reference (a placed pupil: add_pupil's pair, or a
      bench pupil) -> read at nElt-1 (unchanged).
   b. Else, if the deck names a real-plane pupil the caller passed (a flat
      Reference normal to the chief in COLLIMATED space -- s3's
      `SharedPupil`, elt 23, reproduces the sphere read to corr 0.9975)
      -> read there.  Telescope-only harvests on e2e6m want exactly this:
      `exit_pupil_elt = 23`.
   c. Else -> ERROR `macos:dw_dx:noPupil` naming the two remedies
      (`Telescope.add_pupil` / the FP-Return-before-ExitPupil recipe, or
      `exit_pupil_elt` at a collimated pupil Reference).  Never the FP.
      If Dave prefers auto-placement over an error, the recipe is three
      engine calls (see `tools/ep_dome_probe/make_pupil_deck.py` +
      FEX); note it MODIFIES the deck in the session, so it must be loud
      and the placed element must be reported in `out`.

## Gate (replaces "warning fires + map bipolar")

Bipolar is not a valid test: a mean-referenced segment piston is bipolar.
Gate on SHAPE, both non-vacuous against the FP read:
- segment tilt column, within the segment footprint (the piston-poke
  flat-top): a plane fit captures > 90 % of the variance with a near-zero
  mean, and the column rms is within 20 % of 2*alpha*half_segment;
- global field tilt column: rms at the read surface within 30 % of the
  expected tilt (D*theta/sqrt(12) for a filled disk; for s3 the sphere
  read gave 1.40e-6 per 1e-6 rad).
Both are computed in `tools/ep_dome_probe/dome_probe.m`; lift them into
tDwDx on a committed bare-focal deck with the add_pupil pair added by
`make_pupil_deck.py`.

## Engine facts you need (both on macos dev-candidate now)

- `44fc362`: element STOPs preserve the source frame's handedness.  Before
  it, `macos.stop(<any element>)` on a left-handed segmented deck (e2e6m,
  e5hex1, every Telescope emission) mirrored the source grid and obscured
  most rays (732/985 on s3 after `stop(1)`) -- the "add_pupil kills
  rays / corner fields" symptom.  Rebuild the engine before re-measuring.
- `82d8148`: FEX's probe is frame-independent (four probes, medial
  pupil).  Off-axis decks' EP radii moved (j18 3e-4, e5hex1 1 %);
  tPupilFindMethod loses two pins to it (cross-config vertex invariance;
  the e5hex1 FEX-vs-cone-fit gap is now 1.3 mm, i.e. the finders agree).
  Your 3/7 -> 9/1 flip is not reproducible here; on this tree the class
  is 8/10 with both fails attributed to the FEX change, and
  `test_object_space_apstop_deck_needs_no_stop_elt` is one of them, not
  a separate pupil_find defect.
- `eb84095`: STOP accepts Segment elements (add_pupil's default
  `stop_elt = 1` now takes the element path on segmented decks).

## Merge note

sens-core (resources branch, `f144b3a`) hoisted the four _multi
supervisors into `private/dw_multi_core.m`; the single-DOF files your
helper edits are untouched there, so the collision risk is only the
helper's own call sites.  Rebase 4d915f8 onto dev-candidate `38e4808`+.
