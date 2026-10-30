# BRIEF to TO — the journal findings, ruled (Dave, 2026-10-18)

Source of truth: PLAN_CONSOLIDATION.md 1.2, each item now carries its RULED line.

> **STATUS 2026-10-18 (CC took the brief):** items 1-4 DONE, LOCAL on dev-candidate -- (p) ATMOS (root cause
> `RFOURN`'s COMPLEX*8 declaration; stencil 3x; placement refused by name off the source grid frame; record +
> model memo `cli_tests/records/atmos_fix_2026-10-18.md`), (o) filters (the writer), (s) coherent IMG, (d) SYSPROP
> (`dxEltPix`; resources twins re-pinned).  Each item's DONE line in PLAN 1.2 carries its gate and red leg.
> **RULED and DONE 2026-10-18 afternoon (Dave's three answers; CC, LOCAL):** (1) 1.2 (u) -- a deck with no
> output-coordinate block loads: the loader proposes the BEAM frame at the final surface, the CLI asks, the api takes
> it (gate `OUTCORD_missing.jou`); (2) the ATMOS memo's proposals are PLAN §4.6; (3) SAVE writes the CALIB settings
> and LOAD announces them (gate `CALIB_save_roundtrip.jou`); and item 4's SAVE ruling -- a numeric `PropType=`
> resolves to its name, 0 = Geometric (gate `PROPTYPE_numeric.jou`).  The CCMac pre-release brief Dave asked for is
> `BRIEF_ccmac_prerelease.md`.  Items 5-10 stay in the plan.
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
   `EFL_SYSPROP.jou` (unchanged after INT).  **Also from Dave (ruled 2026-10-18): SAVE resolves a numeric
   `PropType=` into its named form** -- `PropType= 0` IS Geometric, and SAVE
   writes `Geometric` (3 -> FarField, and so on through the name table), never
   the garbage the writer prints today for a legacy numeric deck (the older
   1.2 item).  No trace-state keyword; the SAVE policy is untouched.  Gate: a
   deck with every numeric PropType value loads, SAVEs to the names, and the
   saved deck re-loads to the same propagation chain, both compilers.

### ATMOS, after CC's check of the math (2026-10-18, `cli_tests/records/atmos_kolmogorov_check_2026-10-18.md`)

The generator `AtmosPhase` and the RNG are RIGHT (replicated exactly; Kolmogorov structure function
within 20 %, Noll's tilt variance to the digit).  The defect is in `AddAtmos`'s coupling of the screen to
the rays: the screen the rays receive has a pixel-to-pixel checkerboard (lag-1 structure function 70-160x
Kolmogorov), an amplitude that varies by 1e5 with the seed, and on ifx a run-to-run non-determinism
under a fixed seed (uninitialized read -- the ifx debug build shows three different outcomes).  Start
there: the ray-index -> R1 mapping, `rsamp = dxElt(iElt)`, the gradient stencil's 3x constant, and
`ztdir = 0`.  Gates: (a) `seed n; atm; opd 1` deterministic on both compilers; (b) the OPD-1 map's
structure function at lags 1-8 within a factor 1.5 of 6.88 (r/r0)^(5/3) for three seeds; (c) the rms
over the disk at 3-12 rad for r0 0.25 m on the 4 m Cassegrain; (d) spot rms at the focal plane = EFL x
the gradient rms within 10 % (the 3x stencil factor gone).  Then the placement rule and the model memo.

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
