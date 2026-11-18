# Mac release gate — gfortran, Apple Silicon

Pair: macos `22750d4b` + MACOS_resources `e211fe4b` (dev-candidate, nasa-jpl).
Host: Apple M-series, macOS 26.6, gfortran 16.1.0 + Apple clang 21.0.0, cmake
4.4.0 / ninja 1.13.2, Python 3.13.14, MATLAB R2024a (arm64), XQuartz.
Run: 2026-10-19, clean clone outside `~/dev`. Comparison baseline: CC's Linux
gfortran records at `dfc2b476` (= `22750d4` engine content).

## Counts (check — Mac — Linux — differs?)
- 1 clean clone, both repos — both at SHA, clean — — — no
- 2 engine+CLI build (gfortran) — OK, bin/macos arm64, libsmacos.a 189 api symbols — OK (ifx+gfortran) — build-file guard only (makems.sh); ifx n/a on arm64
- 3a CLI load (160 core) — 160/160 ok, rt 153 same / 7 n/a, 8 ulp — 160/160, 153/7 — 2 edge-ray decks
- 3b CLI cmd journals — 70/102 pass (short TMPDIR), 2 must-fail red — 102/102, 2 red — 32 (6 pty-wrap, 21 numeric pins, 2 harness)
- 4 mmacos fast — 645 pass / 15 fail — all pass (tRunCompare excl.) — 15 arm64 tolerance pins
- 5 pymacos (PROPER excl.) — 6935 pass / 8 fail — ~6924 pass — 8 arm64 pins (twins of #4)
- 6 graphics /png + /xw — both render; PNG white-bg + labels — pass — finding 1.2(t) not reproduced on these plots
- 7 self-containment — engine+pymacos run with ~/dev renamed away — — — 4 hidden-dep notes

## Named Mac-only differences (one line each)
- load tst_FF_refr.in — npass 51 vs 60 (9 edge rays lost); arm64-vs-x86_64 FreeForm-refractive edge, not a miscompute
- load tst_cmds_refr.in — npass 56 vs 60 (4 edge rays lost); same edge class
- cmd (6 pty-wrap, pass with TMPDIR=/tmp/cg): EMACS_VI, FILTERS_MULTISPEC, NEW, WINDOW_GBS, WLENS_SWL, ZABERR_ZRM_DOPD
- cmd (2 harness): GPERTURB.jou → untracked fixture pymacos/tests/Rx/Rx_e5hex1.in (tree has e5hex1.in; Rx_e5hex1.in only under GMI/regression); OBS_POLYGON.jou → GNU `sed -i -e` fails on BSD sed
- cmd (21 numeric EXPECT pins on near-zero/precision values, trailing-digit arm64): ATMOS_DMBUILD, ATMOS_seed, BWK, CALIB_SPOT_mixed, CALIB_failed_trials, CALIB_fov_weight, CENTER_FFP_PFP_FDP, COORD_ACOOR, EXPORT, FOV_SETC_SAOPT, IMG, METCALC, MRESET, OPD_reference, OUTCORD_missing, PLOCATE_NOPLOC, PREAD, PROPTYPE_numeric, RRAY_PRAY_TPR_FSR, TRACE_keeps_rx, UDSINIT
- mmacos 15 (tolerance/isequaln/string pins; RelErr ~1e-5..1e-8 vs tol 1e-9..1e-15): tCalibMatrix[mixed_SPOT], tCalibMatrix[trial-loses-chief], tCalibStop, tDefaults(opd_reference), tDefaults(spot_centre), tDwDx, tOpdRef, tPolElement(material-flip, pre-existing Mac pin), tRxLineEndings, tRxRefuse, tRxShortCoef, tSaveFixedPoint, tSrt, tTmaLongslit(pole-frame), tTmaLongslit(record-rescores)
- pymacos 8 (binding twins): test_calib_matrix[mixed_SPOT], test_defaults(opd_reference), test_defaults(spot_centre), test_rx_parse_warnings(look_alike_header), test_save_fixed_point, test_srt, test_zrngr_roundtrip(finite_deck_rejects_nothing), test_zrngr_roundtrip(nscount)

## Build-file changes made (Darwin-guarded; macos repo only)
- makems.sh — portable NJOBS (`sysctl hw.logicalcpu`; nproc is Linux-only); skip the bundled-readline bootstrap on Darwin (CMakeLists links Homebrew readline on macOS)

## Self-containment hidden-dep notes
- cli_tests/corpus_core.txt / corpus_committed.txt carry absolute paths → run make_corpus.sh on any non-author machine (clone both repos as siblings)
- run_mmacos_tests.sh defaults MACOS_HOME → $HOME/dev/macos/macos_f90; pymacos/mmacos default MACOS_BUILD_DIR → ~/dev/macos — a clone elsewhere must set them (unset MACOS_HOME → MEX SIGSEGV)
- GPERTURB.jou references an untracked fixture (above)
- OBS_POLYGON.jou GNU-sed (above)

## Not run / notes
- ifx: no Apple-Silicon arm — gfortran only (expected).
- mmacos tRunCompare + several sensitivity/segmentation tests: filtered from the fast suite (as on Linux); tRunCompare also needs SegMirMaker under build_release_ifx (impossible on arm64).
- No engine-source defect found; all deltas are arm64 round-off through Linux-minted pins, or harness/environment portability.
