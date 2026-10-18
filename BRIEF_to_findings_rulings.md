# BRIEF to TO — the journal findings, ruled (Dave, 2026-10-18)

Source of truth: PLAN_CONSOLIDATION.md 1.2, each item now carries its RULED line.
Rule 1 applies: every fix ships with its gate, and the gate has a must-fail leg
(red on the pre-fix binary); the corpus A/B before the commit for anything that
can move a trace; both compilers.

## Fix now (the freeze admits these: correctness defects on ordinary commands)

1. **(p) ATMOS.** (i) Put the screen's phase generation under SEED: `AddAtmos`
   passes `NoiseSeed(2)` to `AtmosPhase` -- find out why three runs differed
   (was SEED issued? does `AtmosPhase` honour it?), then gate: SEED n, ATM,
   OPD twice -> identical.  (ii) Placement: Dave asks whether the screen should
   be placed only on a Reference surface; the code's own comment says "this
   should be at a pupil", and the loss at a Reflector is the gradient tilt
   applied in a non-pupil frame on a curved element -- measure, then refuse the
   placements that cannot work with a NAMED message (rule 5: say why).
   (iii) Then write a one-page reconsideration of the model for Dave: the screen
   already tilts each ray by its finite-difference gradient and multiplies the
   field by the phase, so "phase without tilt" is not the current state -- say
   what IS, and what a physical model at a non-pupil plane would need.
2. **(o) Filters.** SFILT's file must read back through RFILT with the star
   types it was written with.  Gate: NFILT -> SFILT -> RFILT -> MULTISPEC on
   the Cassegrain (`FILTERS_MULTISPEC.jou` pins today's failure; flip its
   EXPECT).  Which side is wrong, writer or reader, is yours to measure against
   the `.filt` format the manual describes.
3. **(s) Coherent IMG.** Normalize the coherent branch of `SysImging`
   (sunsub.F) so the coherent image of a unit-amplitude object carries the
   PSF's flux.  Gate: `IMG.jou` already pins the incoherent sum rule (sum = 64
   x the PSF sum); add the coherent sum against the incoherent on the same
   object and re-pin.
4. **(d) SYSPROP.** OPD must not set the pixel scale.  SYSPROP's pixel lines
   require a propagation from the source or from the end of the previous
   propagation (INT / PIX) and say so when none has run.  This lands in the
   engine once: SYSPROP is the api's `sysprop` setter as well as the CLI
   command.  Gate: `SYSPROP_after_OPD.jou` (flip to the refusal) and
   `EFL_SYSPROP.jou` (unchanged after INT).  **Also from Dave: SAVE should be
   able to write the current propagation type** -- CC reads this as the
   trace-state single that names the last propagation (the SAVE-policy item of
   1.1b); confirm the reading with Dave before writing it.

## Later, in the plan (phase 3 engine lane)

5. (e) TPR fills `RayPosHist` itself (no api entry today; CLI-only).
6. (k) DOPD after ZABERR / ZRM: fix the previous-OPD bookkeeping.
7. (n) METCALC despace: probe WHY first (does the target point ride with the
   element frame, or is the perturbation unseen?), then decide.
8. (m) LSPOT under NONE and (t) giza's PNG multiplication sign: fix together.
9. (f) e5hex1's `flat.txt`: repurpose an existing FLAT grid file as the fixture.
10. The older recorded items of 1.2: SETP cannot select SPOT; SAVE garbage for
    a numeric `PropType= 0`; CALIB's empty "Zernike modes" headings; GRIN's
    `GradCoef=` branch; the FEX radius doubling as the far-field distance; the
    Refractor coated-branch radiometry audit; the multi-value READ keywords.

CC gates each landing from the clean worktree as before (load, cmd, ab, both
compilers) and folds the cmdref / manual words.
