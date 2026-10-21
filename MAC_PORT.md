# MACOS on macOS (Apple Silicon) — build guide + release-gate record

This file is in two parts:
1. **How to build and test on a Mac** — matches the build files as committed
   (gfortran-only; the CMake / Makefile / `makems.sh` Darwin arms are all in
   the tree).
2. **Validation record** — the last clean-clone release-gate run on real
   hardware, its counts, and the Mac-only deltas from Linux.

---

## Part 1 — How to build and test on a Mac

### The one big fact: gfortran end-to-end
`ifx` has no Apple-Silicon (arm64) target, so **gfortran is the compiler** on
macOS (`gfortran-16` from Homebrew GCC; it is also the MathWorks-supported
Linux Fortran compiler, and the project's gfortran build is clean). The C
compiler is **Apple clang** (the CMake `if(APPLE)` arm downgrades the legacy
C's implicit-declaration errors). The mmacos mex is a **direct gfortran link**
against `libmx/libmex/libmat` (not MATLAB's `mex` driver), so there is no
`mex -setup FORTRAN` gate to satisfy — the Mac mex is a link-recipe port.

### Toolchain (Homebrew + a venv)
```
brew install gcc cmake ninja cairo readline      # GCC 16 → gfortran-16, bare `gfortran` symlink
brew install python@3.13                          # pymacos needs ≥3.13
python3.13 -m venv ~/.venv/pymacos && ~/.venv/pymacos/bin/pip install numpy pytest
```
- MATLAB must be the **native arm64** build (R2023b+); the mex arch is `maca64`.
- **XQuartz** must be installed for the interactive `macos` CLI (`/xw`); the
  headless library + bindings do not need it.
- Debugger is **lldb** (gdb is painful on Apple Silicon).
- A **case-sensitive** clone volume is safest (this tree keeps `.f`/`.F`
  distinct); the default case-insensitive APFS happens to work today because no
  such pair currently collides.

### Build files — all Darwin-guarded, Linux/Windows unchanged
These arms are **already committed**; nothing needs editing on a fresh clone:
- **`CMakeLists.txt`** — `if(APPLE)` adds `-Wno-implicit-function-declaration
  -Wno-implicit-int` to the C flags; the GNU-ld-only `-Wl,--disable-new-dtags`
  is wrapped in `if(NOT APPLE)`; the readline block uses **Homebrew's**
  `libreadline.a` on macOS (`brew --prefix readline`), warning-but-continuing
  if absent; MATLAB auto-detected; vendored **pixman** arch-gated (x86 SSE off
  on arm64, generic C paths).
- **`makems.sh`** — Darwin-guarded: a portable job count (`sysctl
  hw.logicalcpu`, since `nproc` is Linux-only) and the Linux-only bundled-
  readline bootstrap is **skipped on Darwin** (CMake links Homebrew readline).
- **`mmacos/Makefile` + `mmacos/CMakeLists.txt`** — Darwin arm: `maca64` /
  `.mexmaca64`, bundle + `-Wl,-undefined,dynamic_lookup`, plain `-lmx -lmex
  -lmat`, `-mcpu=native`, `-lc++`, libgfortran on rpath; `FC` default via
  `ifeq ($(origin FC),default)`.
- **`pymacos/src/cmake/CMakeLists.txt`** — `elseif(APPLE)` → MACOS64, gfortran
  flags, `-E -cpp -P` preprocess.

### Build commands (verified 2026-10-19)
```sh
# 1. engine + CLI (gfortran) — libsmacos.a + bin/macos
cd <clone>/macos && source ./makems.sh release gfortran
#    outputs: build_release_gfortran/{lib/libsmacos.a, bin/macos}

# 2. mmacos mex (shell FC may be gfortran-16; pass FC=gfortran explicitly)
cd <clone>/MACOS_resources/mmacos
rm -f src/mmacos.mexmaca64
MACOS_BUILD_DIR=<clone>/macos/build_release_gfortran FC=gfortran make
MACOS_HOME=<clone>/macos/macos_f90 \
MACOS_BUILD_DIR=<clone>/macos/build_release_gfortran ./run_mmacos_tests.sh fast

# 3. pymacos (f2py; venv python 3.13)
cd <clone>/MACOS_resources/pymacos/src/cmake
PATH="$HOME/.venv/pymacos/bin:$PATH" cmake -B build_mac -S . \
  -DCMAKE_Fortran_COMPILER=gfortran-16 \
  -DMACOS_BUILD_DIR=<clone>/macos/build_release_gfortran \
  -DPython_EXECUTABLE=$HOME/.venv/pymacos/bin/python
PATH="$HOME/.venv/pymacos/bin:$PATH" cmake --build build_mac
cd <clone>/MACOS_resources/pymacos/tests
MACOS_HOME=<clone>/macos/macos_f90 \
PYTHONPATH=<clone>/MACOS_resources/pymacos/src \
~/.venv/pymacos/bin/python -m pytest --ignore=proper_compare -q
```

### The load-bearing gotchas (still bite)
- **`macos_param.txt` + MEX = SIGSEGV.** `init` loads `macos_param.txt` from
  cwd, exe dir, `$MACOS_HOME`, then a dead compile-time `BUILD_LOC`. `matlab
  -batch` runs where there is no param file, and inside the MEX a *missing*
  param file turns the engine's clean abort into a hard SIGSEGV on the MCR
  thread (empty fault stack — looks like codegen, is not). **Fix: set
  `MACOS_HOME` at the clone's `macos_f90/`.** A standalone Fortran driver shows
  a clean "could not be found, bye!" for the same condition.
- **`run_mmacos_tests.sh` defaults `MACOS_HOME` to `$HOME/dev/macos/macos_f90`**
  (and `MACOS_BUILD_DIR` to `build_release_gfortran`). A clone elsewhere must
  set these, or the MEX SIGSEGVs per the gotcha above.
- **CLI corpus lists carry absolute paths.** `cli_tests/corpus_core.txt` /
  `corpus_committed.txt` are committed with the author's absolute paths; on any
  other machine **run `cli_tests/make_corpus.sh` first** (it rebuilds the lists
  from the local tree and resolves `MACOS_RESOURCES` as the sibling clone, so
  clone both repos side by side).
- **BSD vs GNU in test journals.** A `cli_tests/cmd/*.jou` that shells out with
  `sed -i -e` (GNU form) fails on BSD/macOS sed (`sed -i '' -e` is the Mac
  form). One journal (`OBS_POLYGON.jou`) hits this.
- **macOS `$TMPDIR` is long** (`/var/folders/…/T/…`); the cmd-gate pty is 80
  columns and sets no winsize, so echoed long paths wrap and break `# EXPECT`
  regexes. Run the cmd gate with a short `TMPDIR` (e.g. `/tmp/cg`) on Mac.
- **`timeout` is not a macOS command** (use `gtimeout` from coreutils, or a
  per-process cap).
- The `-reentrancy=none` ifx thread-at-exit fix is irrelevant here (no Intel
  runtime). Bundled `libreadline.a` arrow-key history is replaced by Homebrew
  readline (above).

---

## Part 2 — Validation record

### 2026-10-20 (8) — re-gate (pair macos `3ead7473` + resources `00697be`)
Full record: `cli_tests/records/mac_gate_gfortran_3ead7473.md`. **Engine changed**
since gate 7 — the Zernike-repair line and the CLI startup model-size prompt
(Dave's rulings 1.2 (tt)–(yy)). On Mac: load **164/164** (0 movers vs fcc648c5);
cmd **111/112** (only COORD_ACOOR = (ff); 2 must-fail red; the four new journals
ZERN_refused / ZERN_normhex / ZERN_mode_count / STARTUP_size all pass); mmacos fast
**656/0**; pymacos **6961/0** (+2 skip; the two Zernike test files 20/20);
/png renders clean; self-contained. **The two Mac asks both agree:** NormHex mode-10
SRT sag `5.00004121928100D+00` matches Mahajan & Dai to every printed digit (no
re-pin), and STARTUP_size's shell leg exits with status **1** on macOS (`five-bad
rc=1`, `eof rc=1`, no forrtl). The (tt)–(yy) engine changes land clean on Apple
Silicon; no cross-architecture `ZerntoMon*`/hzdef difference. (ff) persists unchanged (TO).

### 2026-10-20 (7) — re-gate (pair macos `9a971d25` + resources `bcf351c`)
Full record: `cli_tests/records/mac_gate_gfortran_9a971d25.md`. **Engine
byte-identical to gate 6** (`macos_f90/` + CMakeLists untouched since `814082d4`),
so this gate confirms the manual/Zernike/binding-veneer additions land clean on
Apple Silicon. On Mac: load **163/163** (0 movers vs fcc648c5); cmd **107/108**
(only COORD_ACOOR = (ff); 2 must-fail red); mmacos fast **656/0**; pymacos
**6955/0** (+2 skip) — up 6941→6955 as TO's **14 new Zernike-ordering tests all
pass on arm64** (`test_zernike_orderings.py`: the `ZerntoMon*` sag/sign/norm fits
reproduce Linux — no cross-architecture surfsub difference); /png renders clean;
self-contained. The four Zernike findings (tt)–(ww) in front of Dave show the same
on arm64 as Linux (nothing to add from the Mac). (ff) persists unchanged (TO).

### 2026-10-20 (6) — re-gate (pair macos `814082d4` + resources `2ccaeff`)
Full record: `cli_tests/records/mac_gate_gfortran_814082d4.md`. Engine unchanged
in LocalCoord from gate 5; carries the day's dialog fixes and the folded gate 3–5
re-pins + **tOpdRef bound 1e-13**. On Mac: load **163/163** (0 movers vs
fcc648c5); cmd **105/108** (2 must-fail red; 3 fail = COORD_ACOOR (ff), plus two
NEW non-engine pins — BINARY_CENTER a BSD-`stat -c` portability item (the (jj)/(kk)
engine fixes themselves work: file written, no RECL crash, CENTER 4 iters) and
PGD_reopen_units a near-zero `1.293491396D-12` arm64 pin); mmacos **652/5**,
pymacos **6940/1** (+2 skip). **tCalibStop cap/floor marshalling RESOLVED** — the
raw `why` is in-range on arm64 (T=4→3/cap/4, T=10→3/cap/10, T=20→2/floor/14,
`class double`), matching Linux; the gate-5 empty string does not reproduce (the
engine advance shifted the arm64 LM trajectory to track Linux). **tOpdRef green.**
Residual (all arm64 pins or the (ee) FFI mechanism, no engine defect): tCalibStop's
'converged' sub-test (arm64 reaches 'floor' at 29, not 'converged'), trial-loses-chief's
cause string (arm64 now fails at 5 mrad like Linux + restores the tilt, but the
descriptive cause is Fortran unit-6 lost across the mex — (ee); the generic
"calib_run failed" survives), tDwDx array-reproduction (near-zero translation at
AbsTol 1e-6), tPolElement material-flip (isequaln 1-ulp), spot_centre leftover
exact-pins (sign-unstable 2nd component). (ff) persists unchanged (TO).

### 2026-10-19 (5) — re-gate (pair macos `5aff70d` + resources `629a2d2`)
Full record: `cli_tests/records/mac_gate_gfortran_5aff70d.md`. Adds the SAVE DM
fix (AOExample frame restored, UDSINIT re-pinned) + CC's two gate-4 journal
widenings. On Mac: **EXPORT and UDSINIT now pass** (BSD-`wc` fix + the DM fix
work on arm64); cmd 96→**97/103**; load clean 161/161; mmacos 650/7 and pymacos
6939/2(+2 skip) unchanged (resources 629a2d2 identical). CALIB_SPOT_mixed is the
one short of CC's expected 98 — its psiElt **3rd** component (normal ≈ −1) prints
`-9.99…E-01` on arm64 (near-singular tilt), not `E+00`; the pattern's 3rd group
needs that form. The (ff) COORD −O2 bug was bounded by elimination (snan,
int/logical/derived init, -fno-automatic, -fno-tree-vectorize, -ffp-contract=off
all still rotated; correct only at −O0) → a deterministic gfortran-16 arm64 −O1+
codegen miscompile in `LocalCoord`, TO's lane.

### 2026-10-19 (4) — re-gate (pair macos `36699d3` + resources `629a2d2`)
Full record: `cli_tests/records/mac_gate_gfortran_36699d3.md`. The macos tip
advanced past pair 3 (SRS gg, ORS+RSTRACE ii, manual rebuild); resources
`629a2d2` is unchanged, so the binding residual is identical to pair 3 by
construction. Verified: the **gg/ii engine fixes work on arm64** — the new
ORS/RSTRACE/SRS + BEAM journals (cmd 102→103) all pass, load clean 161/161
(ImageDemo now matches Linux). COORD (ff) persists (`LocalCoord` untouched —
still the −O2 rotation, correct at −O0); the (ee) subprocess/capfd Darwin
guards hold (filtered/skipped with reason). Counts: load 161/161, cmd 96/103,
mmacos 650/7, pymacos 6939/2 (+2 skipped). Residual = the round-3 re-pins/
sub-asserts + the two engine items (−O2 LocalCoord, unit-6-across-FFI).

### 2026-10-19 (3) — re-gate (pair macos `d3f2777`/`725e0389` + resources `629a2d2`)
Full record: `cli_tests/records/mac_gate_gfortran_d3f2777.md`. Deltas: cmd **87→96**
pass, mmacos fail **10→7**, pymacos fail **4→2**; glyph still `1x10-10`; ImageDemo
repaired (Mac 578/163 == Linux). The (ee) subprocess/capfd Darwin guards work
(filtered/skipped with reason). Two investigations landed: **(ff) COORD** is an
**`-O2`-only arm64 miscompile in `LocalCoord`** — the frame is correct under the
`-O0` debug build and shows no bounds violation, so an uninitialized local or
`-O2` codegen (TO's lane); **(ee)** the subprocess/capfd empties are **one root
cause** — the gfortran engine's unit-6 (stdout) reaches the process from the
standalone executable but not from a loaded `.so`/mex on macOS, so warning text
should go through the api, not Fortran stdout across the FFI. Residual is round-3
re-pins/sub-asserts and one coreutils item (BSD `wc` leading space in EXPORT.jou).

### 2026-10-19 (2) — re-gate after the fixes (pair macos `3ff1d55` + resources `cd006b2e`)
Second clean-clone gate, on the pair carrying the glyph fix, harness fixes,
cross-architecture re-pins and the lifted 128-char path limit. Full record:
`cli_tests/records/mac_gate_gfortran_3ff1d55.md`. Headline deltas from pair 1:
cmd **70→87** pass, mmacos fail **15→10**, pymacos fail **8→4**; finding 1.2(t)
**fixed** (`spot 5` on /png reads `1x10-10`, was the box); psiElt dead band
dropped the Mac load ulp-settle **8→1**. Residual, for round 2: 14 cmd journals
still use bare `sed -i` (BSD-incompatible — the deck edit silently no-ops);
five bindings subprocess tests still capture empty (the `macos_home_.m`
MACOS_HOME fix did not resolve the nested-`matlab -batch` / Python-subprocess
launch on macOS); a few re-pins applied to the wrong comparison or too tight
(tDwDx, tPolElement, tOpdRef, COORD_ACOOR, CALIB_fov_weight); and two Dave
case-items (mixed_SPOT near-singular margin, spot_centre unstable quantity).

### 2026-10-19 — pre-release gate, clean clone, Apple Silicon
- **Hardware / OS:** Apple M-series, macOS 26.6 (Darwin 25.6).
- **Pair gated:** macos `22750d4b` + MACOS_resources `e211fe4b`
  (`dev-candidate`, nasa-jpl remotes; the gate-green pair CC named — macos's
  `22750d4` is "load 160/160, cmd 102/102" on Linux ifx+gfortran).
- **Toolchain:** gfortran 16.1.0 (Homebrew GCC 16), Apple clang 21.0.0,
  cmake 4.4.0, ninja 1.13.2, Python 3.13.14 (venv + numpy), MATLAB R2024a
  (arm64), XQuartz on `/opt/X11`.
- **Clean clone:** both repos cloned fresh from nasa-jpl into a scratch dir
  (not `~/dev`, not the sandbox); each `git checkout <sha>`.

| # | Check | Mac (gfortran) | Linux (gfortran) | Mac-only delta |
|---|---|---|---|---|
| 1 | Clean clone, both repos | both at SHA, clean | — | none |
| 2 | Engine + CLI build | OK, `bin/macos` arm64, `libsmacos.a` 189 `__macos_api_mod_MOD_` | OK (ifx+gfortran) | one build-file guard (`makems.sh`); no ifx on arm64 |
| 3a | CLI load (160 core decks) | **160/160 ok**, rt 153 same / 7 n/a, 8 ulp | 160/160, 153/7 | 2 edge-ray decks (below) |
| 3b | CLI cmd journals | **70/102 pass**, 2 must-fail red | 102/102, 2 red | 32 fails, all explained (below) |
| 4 | mmacos fast | **645 pass / 15 fail** | all pass (tRunCompare excl.) | 15 arm64 tolerance pins |
| 5 | pymacos (PROPER excl.) | **6935 pass / 8 fail** | 6924 pass | 8 arm64 tolerance pins (binding twins of #4) |
| 6 | Graphics `/png` + `/xw` | both render; PNG white-bg + labels | pass | finding 1.2(t) not reproduced on these plots |
| 7 | Self-containment | engine + pymacos run with `~/dev` renamed away | — | 4 hidden-dep notes (below) |

**No engine-source defect was found.** Every Mac-only difference is either
arm64-vs-x86_64 floating-point round-off surfacing through a Linux-minted tight
tolerance / exact-string pin, or a test-harness / environment portability
point. The engine's optical output matches Linux gfortran to ~8–9 significant
digits on every real-signal deck; large relative diffs are confined to columns
whose value is numerical zero (collimated-system RMS ~1e-12, OPD averages
~1e-15).

**#3a — the 2 edge-ray decks:** `tst_FF_refr.in` (npass 51 vs 60) and
`tst_cmds_refr.in` (npass 56 vs 60) are FreeForm-refractive decks whose rays on
an aperture/TIR edge flip pass/fail on arm64 vs x86_64 — the same class as the
known compiler-sensitive `Grating_example_001`. Not a miscompute (the RMS shift
tracks the lost rays).

**#3b — the 32 cmd fails:** 6 are pty line-wrap artifacts of the long macOS
`$TMPDIR` (they pass with `TMPDIR=/tmp/cg`); 21 are over-pinned `# EXPECT`
regexes on near-zero / precision-sensitive numeric prints differing in the
trailing digits (e.g. the recurring Cassegrain `1.293791606D-12`); 2 are
harness portability — `GPERTURB.jou` references an **untracked** fixture
`pymacos/tests/Rx/Rx_e5hex1.in` (the committed tree has `e5hex1.in`;
`Rx_e5hex1.in` lives only under `GMI/regression/`), and `OBS_POLYGON.jou` uses
GNU `sed -i -e` in a `!`-escape.

**#4 / #5 — the pin fails** are tight-tolerance / `isequaln` / exact-string
assertions minted on Linux (RelErr ~1e-5…1e-8 vs tolerances 1e-9…1e-15, mostly
on near-zero or precision-sensitive quantities): `tDefaults`, `tSrt`,
`tSaveFixedPoint`, `tCalibMatrix[mixed_SPOT]`, `tCalibStop`, `tDwDx`, `tOpdRef`,
`tTmaLongslit`, `tRxLineEndings`, `tRxRefuse`, `tRxShortCoef`, `tzrngr_roundtrip`
and the pre-existing `tPolElement` material-flip 1-ULP pin. The pymacos 8 are
the binding twins of the same tests. `tRunCompare` and several
sensitivity/segmentation tests are filtered out of the fast suite (as on
Linux); `tRunCompare` also cannot run on arm64 (it needs `SegMirMaker` under
`build_release_ifx/`).

**#6 — graphics:** `/png` renders a white-background figure with legible axis
labels, title, and wedge (ray layout + GRAY OPD map inspected; the giza
inverted-black-background bug is absent). `/xw` opens under XQuartz, draws, and
quits with no X11 error. Finding 1.2(t) (multiplication sign on PNG labels) was
**not reproduced** on the layout or OPD plots — they use plain / `E-02`
notation, not a `×10ⁿ` axis multiplier, so they do not exercise it; the
specific Linux plot that shows 1.2(t) is needed to compare.

**#7 — self-containment:** with both `~/dev` trees renamed away, the engine CLI
loaded + traced a deck from the clone (using the clone's `macos_param.txt`) and
pymacos imported — the clone + explicit clone env is self-contained. The hidden
dependencies to fix for a non-JPL clone are the four notes in Part 1's gotchas
(corpus-list absolute paths; `MACOS_HOME`/`MACOS_BUILD_DIR` defaults pointing at
`~/dev`; the `GPERTURB` untracked fixture; the `OBS_POLYGON` GNU-sed).

### 2026-07-22 — first hardware pass (superseded, kept for the deltas)
GCC 16 / gfortran-16, MATLAB R2024a arm64. libsmacos.a built (then 155
`macos_api_mod` symbols; now 189), mmacos **226 pass / 0 fail / 3 skipped**
(suite has since grown to ~660 selected), pymacos import+trace only (no full
suite then), the full `macos` CLI built after the pixman arch-gate fix and
reached the `/xw` prompt. The 2026-07-22 pass established the Darwin build arms;
the 2026-10-19 run is the first full three-suite + graphics + self-containment
gate.
