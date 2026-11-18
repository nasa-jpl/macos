# Mac release gate (5) — gfortran, Apple Silicon

Pair: macos `5aff70d` + MACOS_resources `629a2d2`. This pair adds (over gate 4):
the SAVE DM fix `ee7d8a1c` (AOExample's frame restored, UDSINIT re-pinned), the
manual's generated transcripts, and CC's two journal widenings from gate 4
(EXPORT `^\s*3514`, CALIB_SPOT_mixed tilt 1e-4). Resources `629a2d2` unchanged
from gates 3–4. Host: macOS 26.6, gfortran 16.1.0 + Apple clang 21, MATLAB
R2024a arm64. Run 2026-10-19, clean clone. Baseline: Linux `ee7d8a1c`.

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release+debug — OK, arm64, 189 syms — OK — ifx n/a
- 3a load (161) — 161/161 ok, rt 154/7, 1 ulp — 160/160 — 2 edge-ray only (ImageDemo matches; AOExample frame-restore adds no mover)
- 3b cmd — 97/103 pass (gate4: 96), 2 must-fail red — 103/103 — 6 fail (below)
- 4 mmacos fast — 650 pass / 7 fail, subprocess filtered — all pass — 7 (same as gate 4)
- 5 pymacos — 6939 pass / 2 fail / 2 skipped(darwin) — pass — 2 (same as gate 4)
- 6 /png + /xw — render, glyph 1x10-10 — pass — none
- 7 self-containment — engine CLI + pymacos import, ~/dev renamed — — — none

## What this pair fixed on Mac
- **EXPORT.jou now PASSES** — CC's `^\s*3514` accepts BSD `wc -l`'s right-aligned count.
- **UDSINIT.jou PASSES** — the SAVE DM fix / AOExample frame restore works on arm64 (CC's `opd 16` = 1.2675e-6 is the UDSINIT DM-restore sequence; a bare `old AOExample.in; opd 16` is the nominal 3.163831900D-07, which matches Linux's 3.163831991D-07 to arm64 round-off).
- Load clean, ImageDemo matches Linux, no new movers from the SAVE-DM/manual changes.

## Still failing — cmd (6)
- **CALIB_SPOT_mixed** — CC's tilt widening to 1e-4 fixed the psiElt 2nd component (now matches E-05), but the **3rd component (the mirror normal ≈ −1) prints `-9.9999999987204913E-01` (E-01) on Mac** where the pattern expects `-\d\.\d+E\+00`. The near-singular tilt 1.6e-5 makes z = −√(1−tilt²) ≈ −0.99999999987, just under 1.0, so it formats `9.99E-01`, not `1.0E+00`. Fix: the pattern's 3rd group must accept `(-\d\.\d+E\+00|-9\.9+\d*E-01)`. (Same near-singular story as the binding's `rms<1e-10`→2.7e-9.)
- **COORD_ACOOR** — the −O2 LocalCoord bug (see below); pins the Linux frame.
- **GPERTURB, WLENS_SWL, ZABERR, OBS_POLYGON** — the 4 near-zero/precision pins; CC was waiting on the Mac values (now in `mac_gate_gfortran_d3f2777.md`, pushed): GPERTURB RMS `3.622346865D-05`, WLENS_SWL `1.293491396D-12`, ZABERR secondary `3.566438007D-11`, OBS_POLYGON `1.277068625D-12`.

## Still failing — mmacos 7 / pymacos 2 (resources 629a2d2 unchanged → identical to gates 3–4)
mixed_SPOT rms 2.7255e-9 (>1e-10), tDwDx near-zero element 1.42e-6 vs -2.49e-6 (>AbsTol 1e-6), tPolElement a ray >AbsTol 1e-15, tOpdRef chief-path count 0, tCalibStop stop-reason/n_iter, spot_centre brittle allclose (6305 vs -5938) still present, trial-loses-chief empty on arm64.

## (ff) COORD — the -O1+ miscompile, bounded by elimination (TO's lane)
`coord 6` on CassWithExitPupil.in. Correct ONLY at -O0; the 76.5°-about-x rotation is bit-identical at every build ≥ -O1. Eliminated as cause:
- uninitialized REAL: `-O0 -finit-real=snan -ffpe-trap=invalid` — correct, NO trap.
- uninitialized INTEGER/LOGICAL/derived: `-O1 -finit-integer=-99999 -finit-logical=false -finit-derived` — bit-identical rotated (not it).
- automatic/SAVE: `-O2 -fno-automatic` — rotated (not it).
- vectorization: `-O2 -fno-tree-vectorize` — rotated (not it).
- FP reassociation / FMA: `-O2 -ffp-contract=off` — rotated (not it).
- opt level: wrong at -O1 and -O2, correct at -O0.
So it is a deterministic gfortran-16 arm64 code-generation difference in `LocalCoord` (tracesub.F from ~2924) that first appears at -O1, not explained by uninitialized memory, aliasing-of-automatics, vectorization, or FP contraction. Source-level: compare the -O0 vs -O1 intermediate values through LocalCoord's frame construction (the y/z Gram-Schmidt off xLocal). The journal keeps the Linux frame until fixed.

## (ee) subprocess/capfd — root cause unchanged (gate 3): embedded engine's unit-6 stdout does not reach fd 1 from a loaded .so/mex on macOS; api-route the warnings. Darwin guards hold (3 mmacos filtered, 2 pymacos skipped, with reason).

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate5_*.log + the COORD probe logs /tmp/mac_{coord,noauto,initint,nofma}_*.log; cmd transcripts cli_tests/records/cmd_gfortran_5aff70d_transcripts/.
- No engine *defect* except (ff) the -O1+ LocalCoord miscompile.
