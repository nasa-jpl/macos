# MACOS engine notes — The OPD reference (chief ray vs aperture mean), opdValid, SPOT centring

> **Read this file before touching:** `SUBROUTINE OPD`, `opd_ref_*`, `opd_mask`, `spcOption`, any consumer of `OPDMat`.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## OPD reference: the CHIEF RAY by default, the aperture mean the option (2026-10-10)

*(Documented: Manual §5.2.2 (the OPD definition and the reference paragraph: both references, `UseChfRay4OPD`, OPDREF, the printed line, the lost-chief fallback, `opd_mask`), §4 "Editing .In-Files" + SAVE, §3 "Session options and LOAD"; cmdref OPD / OPDREF, bindings `opd_ref` / `opd_mask`.  Below: the ENGINE record -- where the flags live, what moved, the measurements.  Trimmed 2026-10-16.)*

Dave's ruling (`BRIEF_to_opd_reference.md`): the engine matches the manual --
a ray's OPD is its optical path minus the CHIEF ray's.  `SUBROUTINE OPD`
(tracesub.F; the `propsub.F` twin fills OPDMat during propagation) takes the
chief branch when `LUseChfRayIfOK .AND. LRayOK(1)`, else subtracts `DAvgl`,
the mean over every valid ray.  **Two levels, set where they belong (rule 6):**
`traceutil_mod`'s `LOpdRefChief` is the SESSION option, set `.TRUE.` ONCE at
allocation (`opd_ref_init_opts`, called from `initialise_all` BEFORE
`ray_mod_init`; the CLI's `OPDREF CHIEF|MEAN` changes it and it survives
loads); `LUseChfRayIfOK` is the value in force for the loaded deck,
restored to the session option by every load (`ray_mod_init_vars`), then
overridden per deck by `UseChfRay4OPD= Y|N` and for the loaded deck by
`opd_ref_set`.  `iOpdRefSrc` (0 default / 1 keyword / 2 OPDREF / 3 api)
feeds the line EVERY OPD run prints: `OPD reference: CHIEF RAY (default)`,
`... APERTURE MEAN (UseChfRay4OPD= N)`, `... (OPDREF MEAN)`, `...
(opd_ref_set)`, `... OPDRefRayLen= (Rx)`, and the fallback `OPD reference:
APERTURE MEAN (chief ray lost)` when the chief failed GEOMETRICALLY (no
error -- guards warn).  SAVE writes `UseChfRay4OPD= N` when the mean is in
force; the default writes nothing.  `OPDREF` sits BEFORE the OPD branch in
`macos_cmd_loop.inc` because `LCMP(command,'OPD',3)` is a PREFIX match --
the pre-change CLI reads `opdref` as `opd`.  `design_optim.F` still forces
`.FALSE.` inside optimisation loops and restores the saved value (now the
chief).  **History:** until 2026-10-10 every path ran the mean -- the LOAD
handler set `.TRUE.` ~40 lines before `MBFile6`, whose first act is
`reinitialise_variables()` -> `ray_mod_init_vars()` -> `.FALSE.` (Luis
Marchen's 2026-08-19 diagnosis; that dead line is deleted).
**What the CLI prints does NOT move:** RMS is mean-removed in both branches,
P-V is max-min, and `Average OPD is` prints `DAvgl` = mean(OPL - chief OPL),
the same under either reference -- it IS the constant by which the
chief-referenced map exceeds the mean-referenced one.  Only the per-ray map
moves.  Corpus A/B (406 decks, both compilers,
`cli_tests/records/ab_opdref_moved_51a5624.txt`): status / nPass / lost /
RMS / P-V / round trip identical everywhere; 397 (ifx) / 396 (gfortran)
maps moved by their `Average OPD`; the 6 unmoved are flat collimated trains
(`DAvgl == 0`).  The mmacos `dw_*` drivers pass `opd_ref='mean'` explicitly
and re-applied it after every load; per CC's amendment (a71b43f) their
default is now `'chief'` too, and they call `opd_ref` only when the request
differs from the one in force (`opd_ref_set` dirties the cached trace).
MEASURED (dw_dx, elt 1, 6 DOFs, chief minus mean): every column moves by
ONE CONSTANT on EVERY deck -- segmented (e5hex1, e5pie, e2e6m s2: up to
4.3x the column rms, spread <= 7e-12 of the shift) AND monolithic (e5mono,
e2e6m s1, Cass: 0.3-0.7x the column rms, spread <= 7e-8).  NOT round-off on
monolithic decks, as first expected: a rigid-body poke moves the chief's
own path, so the chief-referenced column carries minus the CHIEF's response
and the mean-referenced one minus the BUNDLE-MEAN response.  Anything that
removes piston is unaffected; only absolute piston moves.
**Zero is not "no ray" any more.**  Under the chief reference a valid ray
whose path equals the chief's reads EXACTLY 0 -- the value `OPDMat` holds
where there is no ray -- so every consumer that found the pupil as `W ~= 0`
silently dropped it (measured: 4 of 276 rays on the tDesignTelescope TMA;
`macos.m2v` built every `dw_*` harvest's rows that way, and `remove_ptt`
then fit planes on the wrong row set).  `traceutil_mod%opdValid` is set at
EXACTLY the statements that write `OPDMat` (tracesub `OPD` x3, propsub x2;
it follows `iObsOpt`, unlike `opdRayMask`, which is always LRayOK.AND.
LRayPass) and exposed as `opd_mask_get` (mmacos `macos.opd_mask`, pymacos
`opd_mask`).  The harvest kernels, the multi-field core, the Telescope WFE
paths and `seg_footprints` now take their rows from it; a map with no
exact-zero valid ray gets the identical rows (same column-major walk).
**Never test `W ~= 0` for validity in new code.**
**The engine's own Zernike fit followed (2026-10-12):** `zern_fit` (zern_wf.F)
took its rows, and ZCO/ZRM/ZAB/`zern_modes_coef` their pupil masks, where
`OPDin /= 0`; now `zern_valid` = `opdValid` (the map's shape must be the OPD
grid's, else the value test).  Blast radius, ifx, 406 decks: 118 of the 380
producing a map carry an exact-zero valid ray -- 88 exactly ONE (the grid ray
coincident with the chief), 30 many (mostly round-off-flat maps; the flat
window `tst_nsref_seq.in` has 724 of 740).  CALIB's LM evaluations run under
the MEAN (`funcs_app` forces it), so WFE_ZMODE is unmoved (`tCalibTargets`).  An all-zero map on
valid rays now fits to zero coefficients (it printed nothing: `inz==0`).
Gate `cli_tests/cmd/ZCOEF_zero_rays.jou` (chief fit == mean fit but piston).

**Why it matters (segmented pupils).**  `DAvgl` is one scalar over the
whole aperture, so perturbing ONE segment shifts it by
`(N_k/N_total)*(mean local response)` and that constant is subtracted
from every ray: unperturbed segments report a spurious piston and the
perturbed one is biased by the same amount.  Measured on `e5hex1`
(7 hex segments, model 128, OPD at the exit-pupil Return, `Tz=1e-8 m` on
one segment): unpoked segments piston by `+2.849e-06` (16.7% of peak);
under the chief reference, exactly `0`, and the poked segment's peak
recovers by that same constant (`1.711e-05 → 1.996e-05`).  **The chief
ray is ALIVE at all five `dw_dx_multi` fields on that deck**
(`LRayOK(1)=1`) — the gate was the flag, not the chief ray, so do not
reach for the "chief dies on segmented decks" explanation without
measuring `ray_info_get`'s `ok_trace(1)` first.

**An OBSCURED chief ray still serves.**  The branch gates on `LRayOK(1)`
— GEOMETRIC — not `LRayPass(1)`.  `LRayOK` is cleared only by a surface
miss / bracket failure (`CTRACE`'s `GO TO 98`) or a non-sequential
dispatch failure; obscuration sets `L1`/`LRayPass` and leaves the
intersection valid (same principle as the SPOT `LocalCoord` fix above).
Measured at the exit pupil: `CassWithExitPupil` and `Rx_Cass_FarField`
have `LRayOK(1)=1, LRayPass(1)=0, RayStatus=Obscured` and the chief
reference IS available (map shifts by an exact constant, std 4.06e-28);
`e5pie` / `e5pie_polyap` / `e5hex1` have a fully unobscured chief.  No
deck checked has a geometrically dead chief.  A pre-fix note claimed
`e5pie` "loses exactly one ray — the chief"; that was the STRUCTURAL
`nPassRays = nRay − 1` (OPD loops `DO iRay=2,nRay`, so the chief is never
written into `OPDMat`) misread as a lost ray.  **Consequence for any RMS compared at 1e-3 (TO, 2026-10-30, IODC 2010 evaluator):** the
`opd()` / `OPDMat` RMS is over every valid ray EXCEPT the chief, so it sits below the RMS of the full
launched set by the chief's share -- 0.1 % on a 741-ray grid, 1e-4 at 7081 (the pupil-centre sample
reads 0 under the chief reference while the mean does not).  `ray_info_get`'s per-ray OPL at the
exit-pupil Return includes ray 1; `mmacos/challenges/iodc2010/iodc2010_eval.m` takes its RMS from that
set and keeps the map's RMS beside it, and `iodc2010_gate` pins the identity (map RMS = per-ray set
without ray 1) to 2e-16.  Not a defect; a convention to state.  Gate:
`tOpdRef/test_an_obscured_chief_ray_still_serves`.

`api`: `opd_ref_set(OK,use_chief)` / `opd_ref_get`; mmacos `macos.opd_ref`,
pymacos `opd_ref` (2026-10-10).  Acts on the LOADED deck (a load restores the
default), and dirties the cached trace via `modified_rx` because OPDMat is
filled DURING the trace.  Gates: `mmacos/tests/tOpdRef.m` (12),
`pymacos/tests/test_opd_ref.py` (6), `cli_tests/cmd/OPD_reference.jou`
(twin deck `ZGD_test_files/tst_opdref_mean.in`).

Related `init` fix: `macos_api_mod`'s `init` rebuilt the engine
(`macos_init_all`, which zeroes `nElt`) without clearing the package
flag `rxLoaded`, so `SystemCheck()` — the guard every wrapper opens with
— kept passing on a wiped model.  Now reset in the rebuild branch only
(an unchanged model size is a no-op and the Rx genuinely survives).
`param_mod_init` answers an unsupported model size with `stop`, which
kills the HOST process when the engine is a mex — `macos.init` now
screens against `macos.model_sizes()` first.

## SPOT / dcdx were CHIEF-RAY-centred: the CLI after every LOAD, mmacos's dw_dx by request (2026-10-06, Luis)
*(Documented: Manual §5 SPOT (SPCenter, shape vs line of sight), §3 "Session options and LOAD" (the five modes, the LOAD notice); cmdref SPCenter / CENTRoid.  Below: the mechanism, the rule and the audit.)*
`spcOption` (the `SPCENTER` option: 1 = spot about the element vertex,
2 = about the chief ray) is a SESSION option.  The 2024 re-init pass
(a9b6f9e) put `spcOption = 0` into `elt_mod_init_vars`, which every Rx
LOAD runs, and the SPOT branch (`macos_cmd_loop.inc`) tests `IF
(spcOption.EQ.1)` with an ELSE that centres on the chief -- so after any
load the CLI's SPOT (its plot and the `TEXT`/`BINARY`/`MAT` spot file)
had the chief ray's position SUBTRACTED from every ray.  A rigid-body tilt
of a mirror therefore moved the printed "Chief ray location" line and
left the spot where it was: Luis's OPTIIX FSM 2e-7 rad test (old MACOS
spot displaced, new MACOS spot unmoved), and the dcdx (centroid
sensitivity) mismatch downstream of the primary -- a PM perturbation
changes the spot's SHAPE, a fold or FSM only its position.  MRESET and
the start-up path set 1 but the load reset ran after them (same trap as
`LUseChfRayIfOK`: the `.TRUE.` set before `MBFile6` never survived
`reinitialise_variables`).  Fixed: `spcOption = 1` at allocation in
`elt_mod_init`, removed from `elt_mod_init_vars`; `SPCENTER` still sets
it.  **The bindings were never affected**: `spot_cmd` issues `SPC` with
the caller's `ref_pos` on every call.  Gate: `scratchpad/spot_cli_gate.py`
pattern (pty CLI; `TEXT`; `SPOT 6 TOUT` on `Rx_Cass_FarField`; `PERTURB
3` by 2e-7 rad; the written spot's mean must move with the chief; the
`--chfray` leg puts the same binary in the pre-fix state and must fail).
**The same hole on the mmacos side, deliberately asked for:** `dw_dx`,
`dw_dgrid`, `dw_dsurf` and `dw_dz_zernike` formed their centroid channel
(dcdx, line of sight) from `macos.spot(..., 'at','chief')`, so a fold or
FSM tilt gave dcdx ~ 0 there too.  Now `'at','elt'`, and `macos.spot`'s
default is `'elt'` (heritage, = pymacos `vpt_center=True`).  Gate
`tDwDx/test_dcdx_of_a_rigid_tilt_is_the_chief_displacement`: the
Cassegrain secondary tilted about x gives dcdx = -6.9 m/rad == the traced
chief's displacement (pre-fix ~0).  A centroid measured about the chief
ray is a SHAPE measure, never a line-of-sight one.
Rule, third time: a value a command or the start-up sets as a USER
OPTION does not belong in `elt_mod_init_vars` -- list of such options
moved so far: the glass catalog, `spcOption`, and (2026-10-10) the OPD
reference's session default `LOpdRefChief` -- the per-load value
`LUseChfRayIfOK` is restored FROM it, never to a constant.  (2026-10-10)
`zcoType` was not moved but its per-load reset now goes to the RIGHT
value, ZernCircular (its command SZCO needs a loaded deck, so per-load
is correct).  **The full audit (2026-10-11, PLAN_CONSOLIDATION 3.3):**
509 per-load resets, 197 with another writer; the test is not 'who else
writes it' but 'can a user set it with NO deck loaded and expect it to
stay'.  A command gated on `ifLoad` (SZCO, SRT, BWK, GBS, LOS, BEAM) and
a per-ELEMENT flag (REGRID) are deck state -- per-load is right.  Five
modes failed the test: POLARIZATION, VECTOR/SCALAR, LNEG, CENTROID, CIR.
**MOVED 2026-10-11 (Dave's ruling (a) + a LOAD notice):** `ifPol`, `ifVecDif3`,
`ifCentroid`, `ifLNegOK`, `ifPGColor` are set in `macos_mod_init` (allocation),
no longer in `macos_mod_init_vars` nor the LOAD handler; every LOAD calls
`ModesInForceNote` -> `LOAD: session modes in force: ...` (nothing when all are
default).  The bindings get the SAME persistence (CC 2026-10-11, on Dave's
one-semantics rule -- flagged in PLAN 4a for him to override): `ModesGet` is the ONE
source, read by the CLI notice and by api `modes_get`; mmacos `macos.load_rx` warns
`mmacos:modesInForce`, pymacos `load` a RuntimeWarning, with the same list
(`macos.modes` / `pymacos.modes`).  The two latent reliances in the suites (628 mmacos
fast + 45 pymacos pol tests): `tTgPol`'s and `tTgPol2`'s scalar detector geometry
(`detGeom`), after a test that left polarization on with crossed polarizers (a 0/0
centroid) -- each now turns polarization off itself.  `STATUS` prints all five;
`cli_tests/cmd/MODE_*.jou` (6, incl. `MODE_load_notice`) pin it.  The CLI `RESET`
still returns four of them to default; it had set `ifCentroid=.TRUE.` (a sixth
site the 2026-08-27 chief-ray ruling missed) -- now `.FALSE.`.  IMGMODE already survives
(`ifAstroMode` lives in plot_mode_mod, outside every `*_init_vars`).

