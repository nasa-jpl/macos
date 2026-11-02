# NOTE to CCMac — (ff), the arm64 COORD frame: source-level bisection in LocalCoord (CC, 2026-10-20)

No new pair: the engine has not changed since 3ead7473.  Dave's word: use the idle cycle on (ff), the one engine
item that reproduces only on your platform.  TO owns the fix; you own the measurement that names the line.  Record
it in `cli_tests/records/mac_ff_localcoord.md` (and a candidate patch as a diff inside that record if (b) below
lands it); push on Dave's word as usual.  No engine change from your lane.

## What is known (PLAN_CONSOLIDATION.md 1.2 (ff))
COORD_ACOOR.jou on CassWithExitPupil.in: zLocal at element 6 reads (4.6e-12, 0.9724, 0.2335) on arm64 at every
gfortran -O1+ build, (0, 0, 1) at -O0 and on every Linux build.  You ruled out, one build each: uninitialized
REAL (-finit-real=snan -ffpe-trap=invalid), uninitialized INTEGER/LOGICAL/derived, automatics vs SAVE,
vectorization, FP contraction.  Bit-identical from -O1 up, so it is one transformation the optimizer is entitled
to by the language -- which usually means the SOURCE breaks a rule the optimizer relies on.

## The suspect (read `macos_f90/tracesub.F` 2925-3142, SUBROUTINE LocalCoord)
`ChfRayDir` / `ChfRayPos` are MODULE arrays (use-associated), not arguments.  LocalCoord saves them
(`oldChfRayDir(1:3) = ChfRayDir(1:3)`, ~2995), then builds the frame by differential chief rays: it perturbs the
module array IN PLACE through calls that receive it as an actual argument --
`CALL DADD(ChfRayDir, oldChfRayDir, th, 3)` and `CALL DUNITIZE(ChfRayDir)` (~3054-3055, again ~3107-3108), traces,
reads `RayPos` back, restores (`ChfRayDir(:) = oldChfRayDir(:)`, several sites), and forms `xLocal` from
`DSUB(xLocal, RayPos, crLocalPos, 3)` + `DUNITVEC` (~3048-3049, ~3071-3079), `yLocal` the same way (~3101-3126),
then zLocal.  Two things an -O1 optimizer may legitimately do differently from -O0 here: (1) a module variable
passed as an actual argument while the callee (or the trace it triggers) also reaches it by USE is ALIASING the
standard forbids -- the compiler may keep `ChfRayDir` in registers across the call or reorder the restore;
(2) `RayPos` is a module array whose base is passed (`DEQUATE(crLocalPos, RayPos, 3)`, ~3025) -- the same class.

## The measurement (one Mac session, both at -O0 and -O1, same deck, `coord 6`)
1. Instrument a COPY of tracesub.F (your gate clone, not the shared tree): after each stage print, in '(A,3ES24.16)',
   `oldChfRayDir`, `ChfRayDir` (after DADD and after DUNITIZE), `crLocalPos`, `RayPos(1:3,1)` after each trace,
   `v1`, `v2`, `xLocal`, `yLocal`, `zLocal`.  Build at -O0 and -O1 (only tracesub.F needs the flag), run
   COORD_ACOOR.jou's `coord 6`, diff the two prints.  The FIRST divergent line is the finding; the earlier lines
   prove the inputs are identical.
2. Experiment (a): the whole tree at -O1 except tracesub.F at -O0 -- if the frame is right, the fault is inside
   LocalCoord (or a routine it calls that tracesub.F defines); if it is still wrong, it is in a callee compiled
   elsewhere (the DADD/DUNITIZE/DEQUATE family in mathsub.F, or the trace in CTRACE).
3. Experiment (b): in the copy, perturb a LOCAL instead of the module array -- `tDir(:) = oldChfRayDir(:) + th(:)`,
   unitize `tDir`, then `ChfRayDir(:) = tDir(:)` -- at both sites, and read `RayPos` into a local before the
   DSUB.  If (b) makes -O1 agree with -O0, the aliasing is the defect and the patch is legal Fortran TO can take as
   is; put the diff in the record.  If (b) does not change it, the record's first divergent line is what TO needs.
4. Keep the legs honest: the -O0 build of the same copy must still print (0, 0, 1) with the instrumentation in
   (otherwise the prints moved the fault -- that is a finding too, say so).

Deliverable: the record with the diff of the two print streams, the (a)/(b) outcomes, and the gfortran version
line.  Linux (CC) will mirror steps 1-3 on x86_64 at -O1 to show the prints agree there, once you have the
first divergent line, so TO sees both platforms in one place.
