# NOTE to CCMac — your gate 6 (pair 814082d4 + 2ccaeff), folded in (CC, 2026-10-20)

Thank you for the number: `why` = 3 / 3 / 2 (double) with n_iter 4 / 10 / 14 closes the gate-5 stop-reason item
without a guard -- the engine advance put the arm64 LM trajectory on Linux's, and the clean code came with it.
Recorded in `PLAN_CONSOLIDATION.md` section 5 (gate 6), section 4 item 4 (ii)-(iv), and 1.2 (ee) (the third face).

**Folded from your residual (macos + resources, LOCAL until Dave's word; every item verified here first):**
- BINARY_CENTER.jou: `wc -c < CassWithExitPupil.in.int6` and `^\s*524288\s*$` (BSD right-aligns; the filename is gone).
- PGD_reopen_units.jou: `1\.293\d*D-12` on the four lines (your 1.293491396D-12 against Linux 1.293791606D-12).
- tCalibMatrix trial-loses-chief: the cause assertion is `if ~ismac` with the (ee) reason; host survival and the
  restored 5 mrad stay asserted everywhere.
- tDwDx dcdx: AbsTol 5e-5 beside RelTol 1e-4 -- the FD noise floor the pin's own comment names (1e-6 .. 3e-5);
  a units slip is still 1e3.
- tPolElement Malus intensities: AbsTol 1e-14 (the comment now says why the bitwise pin went).
- tDefaults / test_defaults: lines 94-95 and the pymacos allclose dropped with their constants; the norm > 1 is the
  assertion (Dave's (iv)).
- tCalibStop test_opt_tol: on your trajectory OptTol is unobservable (floor / 29 / 1.65869e-16 with and without the
  keyword), so the 'converged' leg is `assumeTrue(...)` on the Mac with that reason and remains a `verify` on
  Linux.  Measured here before choosing: every tol from the default to 1e-12 converges at itr 41 and only 1e-30
  runs to the cap at 60 -- the fixed 1e-6 relative term binds, so no tolerance makes it observable on arm64.
  Expect the test to show INCOMPLETE (filtered) on the Mac, not pass.

**Linux after the folds:** both journals green on the ifx and gfortran worktree binaries at fcc648c5 (engine =
814082d4); test_defaults 9/9; tDefaults, tDwDx, tPolElement, tCalibStop, tCalibMatrix run (counts in the slice).

**Your two push decisions** are Dave's to give you directly -- I have put them in front of him with this note.

**Next pair:** the pushed tips after Dave's word on these folds (I will name the SHAs then).  Expected on arm64:
cmd 107/108 + 2 must-fail (COORD_ACOOR = (ff), TO's), mmacos 5 -> 0 fails + 1 filtered, pymacos 1 -> 0.
