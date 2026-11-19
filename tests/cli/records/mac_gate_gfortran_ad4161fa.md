# Mac release gate (9) — gfortran, Apple Silicon — NO ENGINE MISS LEFT

Pair: macos `ad4161fa` + MACOS_resources `00697be`. **Engine changed by ONE line since gate 8
(`3ead7473`): the (ff) `LocalCoord` fix (`0efc2dec`), from CCMac's source-level bisection
(`cli_tests/records/mac_ff_localcoord.md`).** `git diff 3ead7473..ad4161fa -- macos_f90` = `tracesub.F`
(the fix) + `notes/pupil_fex_stop.md` (the method). Host: macOS 26.6, gfortran 16.1.0 + Apple clang 21,
MATLAB R2024a arm64, XQuartz. Run 2026-10-21, clean clone. Baseline: Linux `fcc648c5` (load) / `0efc2dec` (cmd).

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release gfortran — OK, arm64, 189 api syms — OK — ifx n/a
- 3a load (164) — 164/164 ok, 1 ulp settle, **0 movers** vs fcc648c5 — 160/160 — 2 edge-ray only
- 3b cmd — 115 journals, **113 pass, 0 genuine fail**, 2 must-fail red by design — 113/113 — **NONE** (COORD cleared)
- 4 mmacos fast — **656 pass / 0 fail** (94 classes; the Darwin/env filters unchanged) — all pass — none
- 5 pymacos — **6961 pass / 0 fail / 2 skipped(darwin)** — pass — none
- 6 /png + /xw — renders clean — pass — none
- 7 self-containment — engine CLI (Cassegrain RMS OPD 1.293491396D-12) + pymacos import, HOME away — — — none

## (ff) CLOSED — COORD_ACOOR passes on the Mac for the first time
The one-line fix (my candidate, landed verbatim at `0efc2dec` with a comment crediting the bisection):
the y-rotate branch now recomputes `yLocal` from the rotated ray's `RayPos` before `DUNITVEC`, as the
x-branch does. On arm64 at -O2:
- **COORD_ACOOR.jou PASSES** (coaxial Cassegrain): zLocal `(1.98e-11, 1.98e-11, 1.0)` ≈ `(0,0,1)` — the
  76.5° rotation is gone. First time green on Apple Silicon across all nine gates.
- **COORD_offaxis_frame.jou PASSES** (the new gate the fix exposed): opt_example zLocal
  `(2.83e-12, -0.94994, 0.31244)`, e5hex1 zLocal `(2.08e-12, 0.04398, -0.99903)` — both the element normal
  to ≥5 figures, matching Linux's five-figure pins (no last-digit divergence to report).

**The bisection found a real cross-platform bug, not just an arm64 quirk.** CC's landing note records it:
pre-fix Linux (x86_64, every build) read opt_example.in's focal-plane frame **1.4° off** its normal and
e5hex1.in's **12° off** — the stale translate-stage residual is real on non-coaxial decks; the coaxial
Cassegrain of COORD_ACOOR had hidden it with an exactly-zero residual. One line fixes both platforms.

## The Mac line: no engine defect remains
Across gates 1→9 the arm64 residual was progressively classified to: round-off pins (now keyed to shared
digits), coreutils-portability (BSD sed/wc/stat in journals), long-$TMPDIR pty wrap, the (ee)
unit-6-across-FFI mechanism (Darwin-guarded with reason), and the two engine items (ff) and (ee-as-defect).
**(ff) was the last engine miss, and it is now fixed** — on both architectures. The only remaining Mac-only
deltas are harness/environment (the Darwin-guarded subprocess/capfd tests, ifx-absent) and arm64 round-off
within the shared-digit pin policy; every real optical signal agrees with Linux to 8–9 figures and every
engine-command frame (COORD included) now matches.

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate9_*.log; cmd transcripts in the gate clone
  ~/mac_gate_ad4161fa/macos/cli_tests/records/cmd_gfortran_ad4161fa_transcripts/.
- No engine defect. This gate confirms (ff) fixed on arm64 and the whole 7-check suite clean; the method
  (instrument, diff -O0 vs -O1, first divergent intermediate) is now `notes/pupil_fex_stop.md`.
