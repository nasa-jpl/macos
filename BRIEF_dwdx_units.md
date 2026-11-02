# BRIEF for TO: dw_dx translation columns per BASE UNIT (Dave's ruling 2026-10-06)

From Dave via CC.  Fresh session; run on **Opus 5.5**.  Read first: root
`CLAUDE.md`, `CURRENT_SLICE.md` (top entries), the memory index's
`feedback_units_change_consumer_audit` and `feedback_meaningful_tests`
entries (ask CC if you cannot see the memory), and the header of
`mmacos/src/+macos/dw_dx.m` lines 14-60.

## The ruling

`macos.dw_dx` (and the family that shares its convention: `dw_dgrid`,
`dw_dsurf`, `dw_dz_zernike` where they emit rigid-body columns, and
`dwdx_for_current_source`) emits translation columns in **OPD BaseUnits
per SI metre** (the 2026-08-25 convention).  Luis compared against GMI's
dwdx, which is per BaseUnit (per mm on his deck), and read ours as
"1000x too large".  Dave's ruling: add `'trans_output'` with values
`'base'` (OPD-BaseUnits per BaseUnit of translation) and `'si'`
(per metre, today's behaviour), **and make `'base'` the default** --
Luis and every GMI-era consumer expect per-BaseUnit.  Rotations stay
per radian.  The OPD side stays in BaseUnits (that is what 08-25 fixed;
do not touch the numerator).

## What this is NOT

Not a one-line default flip.  The 2026-07 units-change lesson (memory
`feedback_units_change_consumer_audit`): a convention change judged on
a metre deck (CBM = 1) is vacuous -- every consumer looks fine because
the factor is 1.  The audit and the gate must run on an **mm deck**
(`e5hex1.in` is mm: `tDwDx/test_per_dof_delta...` asserts
`out_bu.cbm == 1e-3`).

## The work

1. **Implement** `'trans_output'` in `dw_dx` (default `'base'`): after the
   Jacobian is formed, scale the translation columns (dof_idx 3..5) by
   `cbm` (metres per BaseUnit) so a column reads OPD-BaseUnits per
   BaseUnit of translation; record `out.trans_output`.  Mirror the knob
   where the sibling tools emit rigid-body columns (`dwdx_for_current_source`
   is the shared core -- put the scaling in ONE place).  `dcdx`
   (the centroid columns, fixed this morning to be about the element)
   follows the same rule: per BaseUnit of translation.
2. **Audit every consumer of translation columns** -- the list below is
   `grep -rln 'dwdx\b' --include=*.m src design templates challenges`
   minus `dw_dx.m` itself.  For each: does it multiply a translation
   column by a state in metres (then it needs `cbm` once, or
   `'trans_output','si'`), or by a state in BaseUnits (then the new
   default is what it wanted), or does it only report?  Classify in a
   table in the report: file, use, state units, action, done.  The
   ones that feed CONTROL LOOPS are the ones that bite:
   `design/runners/run_simulator.m` (Tikhonov WFC, the MET loop),
   `run_compare.m`, `run_met.m`, `run_sensitivities.m`, `System.m`
   (`sensitivities` -> `vary`/`optimize`), `met_layout_opt.m`,
   `jacobian_check.m`, `dw_multi_core.m`, `apply_opd_convention.m`, the
   e2e chains (`templates/80_end_to_end/e2e*/s*_*.m`, `r*_*.m`), the
   segmentation examples (`e5_seg*.m`), `run_dwdx_multi.m`,
   `e5hex1/run_dwdx.m`, `e5hex2_refzern/verifyall.m`, the design-layer
   API examples.  Where a consumer writes a record or a figure whose
   numbers are in the deck or a report (e2e6m, rodgers), say so: the
   number changes by 1/CBM and the record's units line must change
   with it -- do not re-run those campaigns, flag them.
3. **Gate, factor != 1:** in `tDwDx`, on `e5hex1.in` (mm): (a) the
   default `'trans_output','base'` translation column equals the
   `'si'` column times `cbm` to 1e-12 relative; (b) a translation poke's
   column has the physical magnitude -- a 1 mm Tz of the focal plane
   element (or a flat) changes the OPD by a known amount in mm; (c)
   `'si'` reproduces today's numbers bit-for-bit (pin them from the
   current build BEFORE you change anything); (d) rotations untouched.
   Plus one consumer gate on an mm deck: `run_sensitivities` /
   `jacobian_check` (or the e2e s4 stage on its smallest fixture)
   producing the same physical wall `dwdx*x + w0` as before for a state
   given in BaseUnits.  The must-fail leg: the pre-change tool fails (a).
4. **Documentation:** the `dw_dx.m` header (replace the 08-25 paragraph's
   translation sentence; keep the history line), the README of
   `templates/50_sensitivities`, a cheatsheet line in
   `mmacos/CLAUDE.md` if one exists there for dw_dx units, and
   `reference_rmswfe_units`-style one-liner for CC's memory in your
   report.  Luis's reply (`macos/DRAFT_email_luis_spot2.md`) promised
   "per-BaseUnit columns directly" -- when this lands, tell CC so the
   note to Luis says the default changed and how to get the old one.

## Rules

`ps -C MATLAB` before any run (CC's gauge runs may be live; one
model-1024 MATLAB on this box).  Do not relink the shared mex (pure
MATLAB change).  Full fast suite before the commit.  American spellings.
Commit locally by path on resources `dev-candidate`; message CC at the
end of step 2 (the audit table, before you change any consumer) and at
the end.  Push only on Dave's word.
