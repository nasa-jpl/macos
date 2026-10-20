# Mac release gate (6) — gfortran, Apple Silicon

Pair: macos `814082d4` + MACOS_resources `2ccaeff`. The engine is unchanged in
LocalCoord from gate 5 (`5aff70dd`) but carries the day's dialog fixes (BINARY
RECL via IOLENGTH + CENTER count + debug-print removal (jj)/(kk)/(ll), giza
`/null` reopen (pp), the units reads (qq), the glass model (rr), the NEW
defaults (ss), EXPORT after a refused BUILD (mm)). Resources `2ccaeff` folds
CCMac's gate 3–5 re-pins (CALIB_SPOT_mixed z, GPERTURB/WLENS_SWL/OBS_POLYGON/
ZABERR bounds, mixed_SPOT 1e-8, tDwDx 1e-5, tPolElement 1e-14, spot_centre
allclose removed, trial-loses-chief either-outcome, **tOpdRef bound 1e-13**) and
drops the `rx_fixture_path` `$HOME/dev` fallback. Host: macOS 26.6, gfortran
16.1.0 + Apple clang 21, MATLAB R2024a arm64, XQuartz. Run 2026-10-20, clean
clone. Baseline: Linux `fcc648c5` (load) / 103-journal cmd.

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release gfortran — OK, arm64, 189 api syms — OK — ifx n/a on arm64
- 3a load (163) — 163/163 ok, rt 156/7 n/a, 1 ulp settle, **0 movers** vs fcc648c5 — 160/160 — 2 edge-ray only
- 3b cmd — 108 journals, **105 pass**, 2 must-fail red by design, 3 fail — 103/103 — 3 (below; 1 engine, 2 arm64/BSD pins)
- 4 mmacos fast — **652 pass / 5 fail** (gate 5: 650/7) — all pass — 5 (below)
- 5 pymacos — **6940 pass / 1 fail / 2 skipped(darwin)** (gate 5: 6939/2) — pass — 1 (spot_centre twin)
- 6 /png + /xw — renders clean (legible labels, white bg, no .notdef/inversion); giza /null reopen OK — pass — none
- 7 self-containment — engine CLI + pymacos run with HOME away, gate-clone env only — — — none (fallback dropped)

## What this pair fixed on Mac (vs gate 5)
- **The five folded journal pins are green**: CALIB_SPOT_mixed (z accepts `-9.9…E-01`),
  GPERTURB, WLENS_SWL, OBS_POLYGON, ZABERR all pass now. cmd 97→105.
- **tOpdRef green** (bound 1e-13, folded in resources): the TMA fixture's four chief-path
  rays at 3.5527e-15 fall below 1e-13, the real signal (5.0449e-13) above. Not in the fail set.
- **tCalibStop cap/floor marshalling RESOLVED**: the raw `why` is in-range on arm64 now —
  T=4 → why=3/'cap'/n_iter 4; T=10 → why=3/'cap'/n_iter 10; T=20 → why=2/'floor'/n_iter 14.
  `class(why)=double`; `macos.calib_stats` maps 3→'cap', 2→'floor' correctly. n_iter 4/10/14
  **matches Linux** (gate 5's empty string / n_iter=29 does NOT reproduce). The engine advance
  `5aff70dd`→`814082d4` shifted the arm64 LM trajectory to track Linux (see the two residuals below).
- The three bindings re-pins (mixed_SPOT 1e-8, tDwDx-rms 1e-5, tPolElement-array 1e-14) and the
  spot_centre allclose removal and trial-loses-chief either-outcome are in force; mmacos 650→652.

## Still failing — cmd (3)
- **COORD_ACOOR** — the (ff) gfortran-16 arm64 `-O1+` LocalCoord miscompile; zLocal still
  `(4.62e-12, 0.9724, 0.2335)` (76.5° about x). LocalCoord untouched this pair. Pins the Linux frame. TO's.
- **BINARY_CENTER** — NEW, **coreutils-portability (not engine)**: the journal probes the file size with
  GNU `stat -c '%s %n'`; BSD `stat` rejects `-c` (`illegal option -- c`), so EXPECT `^524288 …` never
  matches. The engine (jj)/(kk) fixes WORK on arm64: the DIRECT ACCESS file is created with no
  "Write exceeds length" RECL crash (pre-fix gfortran crashed here), no DBINOUT/debug prints, CENTER
  reports "Beam centered in 4 iterations". Fix the journal (`wc -c < FILE` + `^\s*524288`, or a
  BSD/GNU-agnostic probe) — same class as the EXPORT `wc -l` and the `sed -i` sweep.
- **PGD_reopen_units** — NEW, **near-zero arm64 pin (not engine)**: the EXPECT `(RMS OPD error is\s+1\.2937…D-12){4}`
  reads `1.293491396D-12` on arm64 (×4), the recurring Cassegrain near-zero (Linux `1.293791606D-12`),
  the same 1e-12 class as WLENS_SWL/ZABERR. The (pp) giza `/null` reopen and (qq) units reads
  themselves pass (0 "Unknown device"). Widen the pattern to `1\.293\d*D-12` or a magnitude bound.

## Still failing — mmacos 5 / pymacos 1 (all arm64 pins or the (ee) mechanism; no engine defect)
- **tCalibStop/test_opt_tol_keyword_reaches_the_stop** — the 'converged' sub-test, an **arm64 LM-trajectory
  difference**. The default-tol TMA conic solve (`design.Telescope` TMA, dof `[…1]`, fields [0.5 1.0]′)
  stops **'floor' at n_iter=29, merit 1.65869e-16** on arm64 — x86_64 takes two flat accepted steps and
  stops 'converged'. tol=1e-30 gives the IDENTICAL floor/29/1.65869e-16, so `st{2}.n_iter > st{1}.n_iter`
  also fails (lines 77 and 81). Round-off steers the quadratic tail to the ten-rejection floor before a
  flat accepted step. Not marshalling (`why` is clean). CC's call: accept 'floor' as a valid arm64
  outcome for this leg and relax the "runs longer" assert, or make 'converged' detection arch-robust.
- **tCalibMatrix/test_a_trial_that_loses_the_chief_ray_does_not_kill_the_host** — the **(ee) unit-6-across-FFI
  mechanism**. On this pair arm64 now DOES lose the chief at 5 mrad (like Linux — the same trajectory
  shift), taking the CALIB_FAILED branch; the host survives (st=0) and the optics restore (TILT 5.0e-3).
  But `e.message` is the generic **"mmacos: calib_run failed"** (the rtn_flag path) — the engine's
  descriptive cause "no trial step was accepted -- N of K trials lost the chief ray" is Fortran unit-6
  (stdout), lost across the mex FFI on macOS — so the `verifyNotEmpty` cause assertion fails. Same fix as
  (ee): api-route the cause, or Darwin-guard the specific-cause assertion (keep host-survival + restored-tilt).
- **tDwDx/test_trans_output_si_reproduces_the_pre_change_numbers** — near-zero translation dwdx entries at
  AbsTol 1e-6: (3,2) 1.421e-6 vs -2.487e-6 (sign flip), (8,2) 2.274e-5 vs 3.126e-5, (9,2) -5.68e-6 vs
  1.03e-5, (10,2) 1.066e-6 vs 9.24e-6. Same arm64 round-off class as the tDwDx-rms test CC re-pinned at
  AbsTol 1e-5, but this is the separate array-reproduction test (still 1e-6). Widen to ~1e-5 / magnitude-bound.
- **tPolElement/test_normal_incidence_unchanged_by_the_material_flip** — `isequaln` (bit-exact) on
  0.75 and 0.869234279364796, arm64 off by 1 ulp (±2.2e-16). The pre-existing material-flip pin; CC
  re-pinned the OTHER tPolElement test (the 12453-ray compare) to AbsTol 1e-14, this one is still bit-exact.
- **tDefaults/test_spot_centre_default_is_the_element** (+ pymacos twin) — lines 94–95 still
  `verifyEqual(s.shift(3:4)/s.shift(1:2), …, RelTol 1e-9)`. The shift's 2nd component is the
  sign-unstable near-zero-denominator quantity the test's own comment flags (arm64 **6242.98 vs -5877.05**,
  sign flip); the 1st component is also 7.2e-7 relative off (> the 1e-9 pin). The `norm(…)>1` assertion CC
  kept passes; the leftover exact-pin lines 94–95 are the fail. Drop/loosen them (keep norm as the quantity).

## (ee) / (ff) — the two engine-side items
- (ff) persists unchanged (LocalCoord untouched). TO: per-file `-O0` on `tracesub.F` in an `-O1` build;
  argument aliasing in LocalCoord fits the signature as well as codegen (CC's reading, recorded for TO).
- (ee) now has a THIRD face (trial-loses-chief cause string), same root as the subprocess/capfd class:
  the embedded engine's unit-6 does not cross the mex/.so FFI on macOS. The Darwin guards hold (3 mmacos
  subprocess + 2 pymacos capfd filtered/skipped); the api-route lift would also fix the cause string.

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate6_*.log; cmd transcripts cli_tests/records/cmd_gfortran_814082d4_transcripts/
  (in the gate clone ~/mac_gate_814082d4); the two numbers' probes /tmp/mac_gate6_scripts/.
- No engine *defect* except (ff). The cmd BINARY_CENTER/PGD fails are a BSD `stat` and a near-zero pin;
  the mmacos/pymacos fails are arm64 LM-trajectory (tCalibStop converged, trial-loses-chief), the (ee)
  cause string, and three near-zero/bit-exact pins (tDwDx-array, tPolElement-flip, spot_centre leftover).
