# BRIEF_to_bindings §E — the Phase B engine items, 2026-10-17 (TO)

Every number below was read from a run in this session; logs in the session scratchpad.  Engine
edits were made in side worktrees (`~/dev/macos_to_e1`, `~/dev/macos_to_e3`) and committed before
each step ended (E.6); binding gates ran on side mexes / `.so`s linked to those builds (the shared
mex and `.so` were not relinked).  Repos: macos `dev-candidate` from 3ee3081, MACOS_resources
`dev-candidate` from 087714b.

## E.1 — CALIB field weights were 0 after every load (macos f41422e, 410aa9d; resources 9c8abf8)

- Reproduced on the shared 220b68f build (`tst_calib_rigid.in`): FOVS `Weight = 0.0`, chisq 7.981324523439985E-043;
  with `OptFOVWt= 1.0` 7.981324523439981E-033.
- Fix: `fov_weight(:) = 1d0` in `dopt_init_vars` (the LOAD handler's 1.0 kept).  Gate
  `cmd/CALIB_fov_weight.jou` (Weight 1.0, chisq 7.98e-33; `OptFOVWt= 0.5` gives 3.990662261719991E-033,
  exactly half): green ifx + gfortran, red on the shared tree's 220b68f build (= 3ee3081's engine; four
  expectations missed).
- Effect, every CALIB in tCalib / tCalibStop / tCalibMatrix / tCalibTargets / tCalibBlindFields /
  tAsphCalib / tBeamRows / tOptFex / tDefaults, the pymacos calib tests and the CLI CALIB journals,
  pre vs post: 75 of 84 solves the same iterations and stop; every tCalibMatrix cell the same path
  with chisq x1e10 except **mixed x SPOT**; the CLI WFE/ZWF journals unchanged (SETP already set
  weight 1); the Rx_AsphCalib SPOT floor stops 19 -> 17 and 15 -> 18 iterations (still 'floor').
- **mixed x SPOT** (one spot-radius residual, two DOFs, a rank-1 Jacobian): LM is invariant to a
  uniform row scale in exact arithmetic, but this path is round-off chaotic -- `OptFOVWt=` 0.5 and
  2.0 end bit-identical, 1.0 / 1e-3 / 1e-6 each elsewhere (tilt 1.456e-4 / 1.171e-4 / 1.549e-4);
  only the old 1e-10 met T=10 (tilt 1.76e-9).  At weight 1, ifx = gfortran to every digit: T=10
  1.456e-4 (cap), T=15 3.27e-7, T=18..22 'gaussj: singular matrix (2)' in the post-solve covariance
  call -> stop 'failed' with tilt 1.19e-8, T=25 1.6e-8, T=30 1.22e-9, T=35 6.0e-10 (h^4 9.7e-10),
  all cap.  Re-pinned to T=35 in tCalibMatrix, test_calib_matrix.py and CALIB_SPOT_mixed.jou (the
  deck's T=10 stays for its WFE / ZWF cells).  **Open, recorded (notes/calib.md):** the undamped
  post-solve covariance call can hit gaussj's exact-zero pivot on a rank-1 alpha and report a
  finished solve as 'failed'.
- After: tCalibMatrix 18/18, test_calib_matrix.py 18/18 (the host-survival subprocess now honours
  `PYMACOS_SRC`); the eight other CALIB classes as before.
- CLI gates at 09cd0b0 (= f41422e + a comment line), from the worktree: cmd 73 journals each
  compiler, 69 pass + the four TERM-sensitive journals (below) which pass under `TERM=dumb` on both;
  2 must-fail red by design.  load ifx: PASS, 0 decks with a changed OPD / nPass / lost vs 8d04261.
- Docs: Manual Table 13 `OptFOVWt=` = "1 each"; notes/calib.md "Field weights default to 1".

## E.2 — the MRESET hazard's binding face: NOT shared (measurement only)

- Static: every allocate in the bindings' size-change path (`elt_mod_init`, `macos_mod_init`,
  `src_mod_init`, `ray_mod_init`, `dopt_init`, SMACOS's scratch set, api `init`) deallocates or
  `allocated()`-guards every name it allocates (24 statements; a scripted audit, name level).  The
  CLI's `CRWalk` / `CRWalkProjLoc` are used only by BWK (CLI-only) and SMACOS never allocates them.
  The audit's other hits are routine locals, the CLI start-up allocate, and `setbeam.inc` (NPSOL
  path only, not built).
- Runtime (pymacos, one process, the E.1 side `.so`): 128 -> 512 -> 128 -> 256 -> 1024 -> 128 -> 512 on
  Rx_Cass_FarField with trace + OPD + spot + intensity per leg: no allocate / deallocate failure;
  intensity grids 128 / 512 / 256 / 1024 as the size; the repeat legs at 128 (twice) and 512
  bit-identical (rms, every valid OPD pixel, every spot point, every intensity pixel).  e5hex1 the
  same (128 -> 512 -> 128 -> 256 -> 1024 -> 128).  No fix, no gate added.

## E.3 — SYSPROP's pixel lines after an OPD (macos b81021f; resources 256b980, ddb35bb)

- The consumer: `DeriveSysProps` reads `dxElt(iElt)`; DXCALC writes it from a diffraction
  propagation (INT / PIX) OR, after an OPD, the ray-grid spacing (LOS reads that on purpose).
  Measured on the manual's Cassegrain (pre-fix): OPD alone -> pitch 4.4094850E-12, lambda/D
  1.2741640E+06 px; INT alone -> 5.2290937 px, 9.8614033E-03 arcsec/px, pitch 1.0744514E-06;
  OPD then INT -> the INT numbers; **INT then OPD -> the 4.4e-12 numbers again**.  Binding face
  (Rx_Cass_FarField, model 128): lamD_px 3.3130544e6 after opd(), 1.0177264 after intensity(),
  3.3130544e6 after intensity() then opd().
- Fix: `dxEltDiffr(iElt)` (elt_mod, own deallocate/allocate pair, reset at load) set by DXCALC's
  three diffraction branches, cleared otherwise; the pixel lines need it; the hint names INT / PIX.
  dxElt itself unchanged.  Gate `cmd/SYSPROP_after_OPD.jou`: green ifx + gfortran, red on the E.1 worktree's
  09cd0b0 build (= 410aa9d's engine; three missed); EFL_SYSPROP.jou green.  Twins: mmacos `tSysProp/test_px_not_from_an_OPD` 4/4 on the
  fix's mex, the new leg red on the pre-fix mex; pymacos `test_sysprop.py` red pre-fix, green with
  `test_api_coverage` (161 passed).

## E.4 — TPR (macos 42c7ade)

- Pre-fix, CassWithExitPupil: every TPR position 0 without `SaveRayPosH= Y`; and after RRAY the
  four empty slots were traced as "ray 0" (`iRay = 0` rows; NaN with the key: RayPos(:,0)).
- Fix: TPR prints RayPos at the end element (what RayPosHist copies), formatted, and skips slots
  <= 0 (a failed chief no longer stops the loop before the marginal rays).  Gate
  `RRAY_PRAY_TPR_FSR.jou` pins the numbers the old TPR printed WITH the key (ray 63
  -9.075851181705730E-11 8.250726490910409E-12 1.5; ray 108 x -5.768541200268373E-11): green ifx +
  gfortran, red on the 09cd0b0 build (pre-E.4).  TPR is CLI-only (no binding twin).

## E.5 — UDSrfParam six vs ten (macos 72ea712, docs only)

Six is the width everywhere (parser `ReadRealsWarn(...,6)`, MODIFY, the element-entry prompt,
SAVE, every kernel dummy `param(6)`); the models read at most three; `mUDSrfParam = 10` is only the
allocation.  The manual's UserDefined section now documents UDSrfParam and lists the 11 types the
code dispatches (its old list called type 1 a 349-actuator mirror; the code's is 1313).

## E.6 / found on the way

- **The four "cwd" journals are TERM-sensitive** (EXPORT, RX_FID_LS, SESSION_INFO, SHELL_LONG): under
  `TERM=xterm-256color` readline writes bracketed-paste escapes (`ESC[?2004l`) at the start of output
  lines and their `^`-anchored EXPECTs miss; under `TERM=dumb` all four pass on both compilers.  My
  2026-10-16 record blamed the cwd -- wrong.  One-line harness fix for CC: set `TERM=dumb` in the
  child before `execv` (cli_load_gate.py `Cli`).

## Suites at the end (the combined build: worktree at 42c7ade = E.1 + E.3 + E.4 on 5bb836f)

- CLI, `TERM=dumb`, from the worktree: cmd 77 journals each compiler, 76 pass + `METCALC.jou`, whose
  worktree copy was CC's pre-79802d7 version (CC has since re-pinned it as measured; the current
  journal passes on both compilers); 2 must-fail red by design.  load 114/114 both, round trip same
  114, 0 decks with a changed OPD / nPass / lost vs 8d04261; A/B 8d04261 -> 42c7ade 0 movers both.
- mmacos fast (size 128, side mex on the combined engine): 641 pass, 0 fail; 33 filtered by
  assumption because the worktree has no SegMirMaker binary (tSegMirMaker, tSegmentRx, tEdgeSensors,
  tMet/test_add_met_stewart_truss, three tRunCompare) -- not run, not affected by these changes.
- pymacos (tests/, excl. proper_compare + sensitivities, side .so): 6895 pass, 5 fail -- all five
  subprocess tests that import the worktree's src/ (no .so there); with the side .so placed there
  the six files (test_calib_targets, test_rx_line_endings, test_rx_refuse, test_rx_short_coef,
  test_zrngr_roundtrip, test_calib_matrix) pass 27/27.  So 6900 / 0.
- The shared trees' mex and .so are NOT relinked: they do not carry E.1 / E.3 until CC's next
  shared rebuild.
