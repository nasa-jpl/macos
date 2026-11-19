# Mac release gate (2) — gfortran, Apple Silicon

Pair: macos `3ff1d55` + MACOS_resources `cd006b2e` (dev-candidate, nasa-jpl).
This pair carries the glyph fix (`1632450d`), the harness fixes (`7acb2911`),
the cross-architecture re-pins (`cd006b2e`), and the lifted 128-char path limit.
Host: Apple M-series, macOS 26.6, gfortran 16.1.0 + Apple clang 21.0.0, MATLAB
R2024a (arm64), XQuartz. Run 2026-10-19, clean clone outside `~/dev`.
Engine == CC's Linux gate `4d7b1aa1` (only records differ 4d7b1aa1→3ff1d55);
comparison baselines `load_gfortran_4d7b1aa1.csv` / `cmd_gfortran_4d7b1aa1.csv`.

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA, clean — — — no
- 2 build (gfortran) — OK, bin/macos arm64, libsmacos.a 189 api symbols — OK — makems already Darwin-patched (985913a); ifx n/a on arm64
- 3a CLI load (161 core) — 161/161 ok, rt 154 same / 7 n/a, 1 ulp settle — 160/160 — 2 edge-ray decks; +Rx_e5hex1.in (new fixture); ulp-settle 8→1 (psiElt dead band)
- 3b CLI cmd — 87/102 pass (was 70), 2 must-fail red — 102/102 — 15 fail: ~9 GNU-sed-in-journal, ~5 near-zero pins, 1 CALIB margin, 1 fixture
- 4 mmacos fast — 650 pass / 10 fail (was 15) — all pass — 10: 3 subprocess-empty, re-pins under-widened, 2 case items
- 5 pymacos (PROPER excl.) — 6939 pass / 4 fail (was 8) — ~pass — 2 subprocess-empty, 1 margin, 1 unstable-quantity
- 6 graphics /png + /xw — render; glyph FIXED (1x10-10, was box) — pass — none (1.2(t) confirmed fixed)
- 7 self-containment — engine + pymacos run with ~/dev renamed away — — — same 4 hidden-dep notes (all now mitigated except the sed sweep)

Progress vs pair 1 (22750d4): cmd 70→87 pass, mmacos fail 15→10, pymacos fail 8→4; glyph fixed; ulp-settle 8→1.

## Still failing — classified for round 2

### GNU-sed-in-journal (NOT pins; 14 cmd journals use bare `sed -i`)
Bare `sed -i '...'` is BSD-incompatible (BSD reads the next token as the backup
suffix); GNU `a\`/`i\` append-insert and `{}` address blocks also differ. The
`!`-escape silently no-ops the deck edit on macOS, so the test setup never
happens. CC fixed 4 journals (→ `-i.bak`) but 14 still use bare `sed -i`.
Failing from this class: BWK (`a\ CalBWk= Y`), WLENS_SWL (`a\ OptWaveLen=`),
PROPTYPE_numeric, FOV_SETC_SAOPT, METCALC, IMG, OBS_POLYGON (its `-i.bak` is
set but a `-e '1,/addr/{/x/d}'` block still breaks BSD sed: `sed: 1: ...`).
Fix: make journal deck-edits portable (perl -i -pe, or stage edits in the
Python cmd-gate), not just `-i.bak`.

### Near-zero / precision EXPECT pins not re-pinned or under-widened (arm64)
- WINDOW_GBS: elt-5 beam diameter mac `2.1431740984723464E-010` vs EXPECT `2.14310100…D-10` (5th digit; ~1e-10).
- EXPORT: Zernike coef mac `0.285169123616263D-10` vs `0.28515551…D-10` (5th digit; ~1e-10/1e-11).
- ZABERR_ZRM_DOPD: `1.293491396D-12` vs `1.293791606D-12` (Cassegrain near-zero).
- COORD_ACOOR: CC's re-pin pattern assumes zLocal≈(0,0,1) but it is (4.6e-12, 0.9724, 0.2335); xLocal (-1,0,1.98e-11), yLocal (1.93e-11,-0.2335,0.9724) — pattern shape wrong.
- CALIB_fov_weight: merit chisq mac `2.28e-29` vs EXPECT `~8e-33` — pin the exponent to a magnitude bound, not E-033.

### mmacos re-pins applied to the wrong comparison / too tight
- tDwDx: still > reltol 1e-4 on at least one dw/dx entry (give the entry a wider reltol or 5-fig pin).
- tPolElement: dead band set on the OPD twin, but the failing compare is a 12453-index array at ~2.2e-16 — widen that array's compare.
- tOpdRef: `|W|<1e-18` bound set, but the assertion that still fails is a chief-path ray COUNT `> min` reading 0 on arm64.
- tCalibStop: stop-reason char array empty / iteration count differs on arm64.

### Subprocess-empty — CC's macos_home_.m fix did NOT resolve it on Mac
tRxShortCoef, tRxRefuse, tRxLineEndings (mmacos) and test_rx_parse_warnings
(::look_alike, ::optasph_above_varelt) (pymacos) still capture **empty** output.
Passing MACOS_HOME to the nested `matlab -batch` (mmacos) / the Python
subprocess (pymacos) was not the root cause: a `system()` launched from inside
MATLAB, and the pymacos subprocess, still cannot run the engine on macOS
(DYLD/SIP env-stripping of the child, or the nested MATLAB cannot load the mex).
These need a different approach (or skip-on-Darwin with a reason).

### Case items (Dave's rulings — still failing on Mac)
- tCalibMatrix/test_a_well_posed_cell_removes_its_spoil[mixed_SPOT] + pymacos twin: near-singular SPOT cell, arm64 LM lands tilt `1.5997e-05` (vs `< 1e-8`); the margin widening wasn't applied here (or 1.6e-5 exceeds it).
- tDefaults/test_spot_centre + pymacos twin: SPOT_SHIFT 2nd component mac `6305.26` vs lin `-5938.92` — opposite sign, unstable near-zero-denominator quantity; needs the "other quantity" rethink, not a pin.
- tCalibMatrix/test_a_trial_that_loses_the_chief_ray: expected non-empty "N of K lost the chief" but got empty — on arm64 the trial step did not lose the chief.

## Confirmed fixed on Mac this pair
- Finding 1.2(t): `spot 5` on /png now reads `1x10-10` (ASCII x), was the .notdef box.
- Re-pins that landed: tSrt/test_srt, tDefaults opd_reference/test_opd_reference, tSaveFixedPoint/test_save_fixed_point, tTmaLongslit ×2, test_zrngr_roundtrip ×2 (now /tmp-staged), and the 11 cmd journals (ATMOS*, MRESET, OPD_reference, OUTCORD, PLOCATE, PREAD, RRAY, TRACE_keeps_rx, UDSINIT, CENTER_FFP_PFP_FDP).
- psiElt dead band: Mac load ulp-settle 8 → 1.
- GPERTURB fixture Rx_e5hex1.in now present (but the deck needs GridFile `flat.txt` staged too — still fails on `File flat.txt does not exist`).

## Not run / notes
- ifx: no Apple-Silicon arm — gfortran only.
- Raw logs: /tmp/mac_gate2_{makems,load,cmd,mmacos,pymacos_test}.log; cmd transcripts in cli_tests/records/cmd_gfortran_3ff1d55_transcripts/.
- No engine-source defect; all deltas are arm64 round-off, GNU-sed-in-journal, subprocess-env, or the two case items.
