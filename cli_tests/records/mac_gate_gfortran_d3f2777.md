# Mac release gate (3) — gfortran, Apple Silicon

Pair: macos `d3f2777` (engine == `725e0389`, ImageDemo repaired; `d3f2777` is the
NOTE commit on top, empty engine diff) + MACOS_resources `629a2d2`. This pair
carries CC's round-2 fixes: perl `-i -pe` journals, the flat.txt fixture, the
re-pins, the Darwin guards on the subprocess/capfd tests, the case decisions,
the ORS/Gaussian engine fixes. Host: macOS 26.6, gfortran 16.1.0 + Apple clang
21, MATLAB R2024a arm64, XQuartz. Run 2026-10-19, clean clone outside `~/dev`.
Baselines: Linux `4d7b1aa1` (load/cmd); ImageDemo moves legitimately (repair).

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA, clean — — — no
- 2 build release+debug gfortran — OK, arm64, 189 api syms — OK — makems already patched; ifx n/a
- 3a CLI load (161) — 161/161 ok, rt 154 same / 7 n/a, 1 ulp — 160/160 — 2 edge-ray; ImageDemo 578/163 == repaired Linux
- 3b CLI cmd — 96/102 pass (70→87→96), 2 must-fail red — 102/102 — 7 fail (below)
- 4 mmacos fast — 650 pass / 7 fail (15→10→7), 3 subprocess filtered-by-assumption — all pass — 7 (below)
- 5 pymacos (PROPER excl.) — 6939 pass / 2 fail (8→4→2), 2 capfd skipped(darwin) — pass — 2 (below)
- 6 /png + /xw — render; glyph still 1x10-10 — pass — none
- 7 self-containment — engine CLI runs with ~/dev renamed away — — — none new

Progress: cmd 87→96, mmacos fail 10→7, pymacos fail 4→2. The (ee) guards work: 3 mmacos
subprocess tests "Filtered by assumption", 2 pymacos capfd tests skipped(darwin), both WITH reason.

## Two investigations CC asked for

### (ff) COORD_ACOOR — an -O2 arm64 miscompile in LocalCoord, NOT a pin
`coord 6` on CassWithExitPupil.in, three ways:
- (a) RELEASE, right after `old`: xLocal (-1,0,1.98e-11), yLocal (1.93e-11,-0.2335,0.9724), zLocal (4.6e-12,0.9724,0.2335) — a frame rotated 76.5° about x (Linux: (-1,0,0)/(0,-1,0)/(0,0,1)).
- (c) DEBUG build (`makems.sh debug gfortran`, -O0 -fcheck=all): **correct** — (-1,0,0)/(0,-1,0)/(0,0,1), **no bounds violation reported**.
- (b) after `opd 6`: moot — `opd 6` opens a device prompt that eats the next line; (a) already shows the rotation pre-trace, so it is not trace state.
Conclusion: the rotation is **release-only (-O2) on arm64**, correct at -O0, with no array-bounds error → an uninitialized local exposed at -O2, or -O2 codegen, in `LocalCoord` (tracesub.F from ~2924; locals del/th/v1/v2, xLocal/yLocal/zLocal INTENT(OUT)). TO's lane. The journal keeps the Linux frame as its pin until fixed.

### (ee) subprocess/capfd empty — unified root cause: embedded engine's unit-6 stdout
The standalone CLI prints `MBFile6: Unidentified header string = OptFooKey` on Mac (its journal passes). The SAME engine loaded as the pymacos `.so` emits that warning to **nothing** — terminal, Python fd-redirect, subprocess PIPE, and `GFORTRAN_UNBUFFERED_ALL=y` all empty. So it is not MACOS_HOME and not buffering: **gfortran unit-6 (stdout) reaches the process from the standalone executable but not from a dynamically-loaded `.so`/mex on macOS** (it does on Linux — why capfd works there). The mmacos nested-`matlab -batch` tests fail for the same reason (the mex's unit-6 output never reaches the captured stdout). Portable fix: route warning text through the api (a `load_warnings_get`), not Fortran stdout across the FFI. CC's Darwin guards (skip/assume) are the right stopgap; this is the mechanism for lifting them.

## Still failing — round 3 (CC's round-2 widenings tripping the next bound on arm64)

### cmd (7)
- COORD_ACOOR — the -O2 bug above (pins the Linux frame).
- EXPORT — **BSD `wc -l` leading space**: journal `! wc -l` → `    3514 /tmp/...`, EXPECT `^3514 \S+` (anchored) fails. New coreutils-portability item (same class as the sed sweep): use `^\s*3514` or `wc -l < file`.
- CALIB_SPOT_mixed — tilt `1.5997e-05` (E-05); the binding was widened to tilt<1e-4 but the journal's psiElt EXPECT still allows only E-06..E-09. Widen the journal pattern to match the binding ruling.
- GPERTURB — fixture flat.txt now staged + read (fix worked); residual is a numeric pin: RMS `3.622346865D-05` vs `3.622346888D-05` (8th sig fig; rule-2, cut to 6 figs).
- WLENS_SWL — perl fix works (2 wavelengths now); residual EXPECT `RMS OPD error is 1.293791606D-12` (Cassegrain near-zero) not re-pinned here; Mac `1.293491396D-12`.
- ZABERR_ZRM_DOPD — first RMS line re-pinned OK; secondary `3.566438007D-11` vs `…435264D-11` and the dOPD/WFE ~1e-11 lines not re-pinned.
- OBS_POLYGON — perl fix works (polygon obscuration applied, 145 rays lost); residual near-zero RMS pin (`1.277068625D-12`).

### mmacos (7) — exact Mac values for the next cut
- tCalibMatrix mixed_SPOT — `tilt<1e-4` passes, but the later `rms < 1e-10` fails: Mac **2.7255e-09**. (pymacos twin same.)
- tCalibMatrix trial-loses-chief — expected non-empty "N of K lost the chief"; got EMPTY — on arm64 the trial step did not lose the chief even with the explicit 5 mrad tilt.
- tCalibStop — stop-reason char array not equal + n_iter below the pinned minimum (range floor|cap, 11..40); arm64 stop reason / iteration count differ.
- tDwDx — AbsTol 1e-6 + RelTol 1e-4 still exceeded at index 15 (3,2): Mac **1.4210854715202e-06** vs **-2.48689957516035e-06** (both ~1e-6 near-zero, sign-flip; Error 3.9e-6 > AbsTol 1e-6) — a near-zero dw/dx element; widen abstol to ~1e-5 or magnitude-bound it.
- tOpdRef — `|W| < 1e-15` set, but the surviving assertion is a chief-path ray COUNT `> min` reading 0 on arm64.
- tPolElement — AbsTol 1e-15 still exceeded by some ray in the 12453-index Ex/Ey/Ez compare (first rows are ~2.2e-16, but a worse one is beyond 1e-15); widen to ~1e-14 or a few-ulp relative.
- tDefaults spot_centre — the stable `norm>1` assert was added but the old brittle `allclose(s[2][:2], SPOT_SHIFT_CHIEF, rtol=1e-9)` is still present: Mac `6305.26` vs `-5938.92` (sign-unstable near-zero-denominator 2nd component); remove/replace that assert.

### pymacos (2)
- test_calib_matrix[mixed_SPOT] — `rms < 1e-10` fails at Mac **2.7255e-09** (twin of mmacos).
- test_spot_centre — the brittle `allclose(s[2][:2], SPOT_SHIFT_CHIEF)` still present (6305 vs -5938; twin of mmacos).

## Confirmed working this pair
- Glyph (1.2 t): `spot 5` /png reads `1x10-10` (ASCII x).
- perl `-i -pe` journal rewrite: WLENS_SWL/OBS_POLYGON/BWK/PROPTYPE_numeric/FOV_SETC_SAOPT/METCALC/IMG now apply their deck edits on macOS.
- GPERTURB flat.txt fixture staged + read.
- ImageDemo repaired: Mac 578 pass / 163 lost == repaired Linux.
- The subprocess/capfd Darwin guards: filtered/skipped with reason (ee pass condition).

## Not run / notes
- ifx: no Apple-Silicon arm — gfortran only.
- Raw logs /tmp/mac_gate3_*.log; cmd transcripts cli_tests/records/cmd_gfortran_d3f2777_transcripts/.
- No engine-source *defect* except (ff) the -O2 LocalCoord miscompile; everything else is arm64 round-off, coreutils-portability, the (ee) unit-6 mechanism, or re-pins/sub-asserts to finish.
