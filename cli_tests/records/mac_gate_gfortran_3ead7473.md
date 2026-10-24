# Mac release gate (8) — gfortran, Apple Silicon

Pair: macos `3ead7473` + MACOS_resources `00697be`. **Engine CHANGED** since gate 7
(`9a971d25`): the Zernike-repair line and the CLI startup model-size prompt, all
ruled by Dave 2026-10-20 (PLAN_CONSOLIDATION 1.2 (tt)–(yy)). Engine files touched:
elt_mod.F, hzdef.inc, iosub.inc, macos_api_mod.F90, macos_cmd_loop.inc, macos_mod.F,
macosio.F, msmacosio.inc. Host: macOS 26.6, gfortran 16.1.0 + Apple clang 21, MATLAB
R2024a arm64, XQuartz. Run 2026-10-20, clean clone. Baseline: Linux `fcc648c5` (load).

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release gfortran — OK, arm64, 189 api syms (engine rebuilt; mex 1800704b) — OK — ifx n/a
- 3a load (164) — 164/164 ok, rt 157/7 n/a, 1 ulp, **0 movers** vs fcc648c5 — 160/160 — 2 edge-ray only
- 3b cmd — 114 journals, **111 pass**, 2 must-fail red by design, 1 fail (COORD_ACOOR = (ff)) — 112/112 — 1 ((ff), TO's)
- 4 mmacos fast — **656 pass / 0 fail** (94 classes; the Darwin/env filters unchanged) — all pass — none
- 5 pymacos — **6961 pass / 0 fail / 2 skipped(darwin)** (gate-7 base 6955 + 6 new Zernike tests) — pass — none
- 6 /png + /xw — renders clean (spot diagram legible, white bg) — pass — none
- 7 self-containment — engine CLI (Cassegrain RMS OPD 1.293491396D-12) + pymacos import, HOME away — — — none

## The Zernike-repair line + startup prompt — all green on arm64
- **(tt)/(vv) NormAnnularNoll / NormAnnular refused** (load, MOD, api, and the ZCO output fit) + **ExtFringe modes
  47–49/58–64 refused by mode**: journal `ZERN_refused.jou` (6 legs) PASS.
- **(uu) NormHex OPENED** (ZernType_MaxMode 45, hzdef.inc now double): journal `ZERN_normhex.jou` PASS. **The Mac ask —
  the SRT sag of mode 10 — matches to every printed digit**: arm64 Segment Length `5.00004121928100D+00` (sag
  4.1219281000687e-05 m = 1e-4·4√(70/103)·0.5³·cos0) == Mahajan & Dai's closed form `5.0000412192810…`; the (0,0.5)
  cos(3A)=0 ray meets z=0 (`5.00000000000000D+00`). No last-digit divergence → no 5-figure re-pin needed.
- **(xx) coefficient past a type's mode count refused by mode** at every face: journal `ZERN_mode_count.jou` PASS.
- **(yy) the CLI startup model-size prompt** refuses an unrecognized answer and re-asks (≤5), end of input / fifth
  refusal stop with exit status 1: journal `STARTUP_size.jou` PASS. **The Mac ask — the shell leg's exit status —
  reads 1 on macOS**: `five-bad rc=1` and `eof rc=1`, no `forrtl` / Fortran runtime error (the pre-fix CLI crashed here).
- **(ww) pymacos veneer limits read the generated constants**: `test_zern_veneer_limits.py` 4/4; the updated
  `test_zernike_orderings.py` 16/16 — the two Zernike files 20/20 on arm64. No cross-architecture `ZerntoMon*`/hzdef
  difference (the sag/sign/normalization fits reproduce Linux).

## cmd — the single fail (unchanged, TO's)
- **COORD_ACOOR** — the (ff) gfortran-16 arm64 `-O1+` LocalCoord miscompile; zLocal still the 76.5°-about-x rotation
  (LocalCoord untouched by this engine change). Pins the Linux frame. All of gate 6's folds remain green; nothing regressed.

## Count note
CC's NOTE predicted pymacos 6951 (6941 + 6 + 4), but the gate-7 base was 6955 (not 6941), so the actual is 6955 + 6 =
**6961** (test_zernike_orderings 14→16 = +2, test_zern_veneer_limits +4 → the "6 new" are counted against the gate-7
base, which already carried the first 14 Zernike tests). 0 failures either way.

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate8_*.log; cmd transcripts in the gate clone
  ~/mac_gate_3ead7473/macos/cli_tests/records/cmd_gfortran_3ead7473_transcripts/.
- No engine *defect* except (ff). The (tt)–(yy) engine changes all land clean on Apple Silicon; the two numbers CC
  asked for (NormHex mode-10 sag digits; STARTUP_size exit status) both agree with Linux / read 1 on the Mac.
