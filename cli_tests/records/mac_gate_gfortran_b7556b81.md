# Mac release gate (10) — gfortran, Apple Silicon

Pair: macos `b7556b81` + MACOS_resources `d436b9c`. Engine changed (TO's cleanup round 1 + three of
Dave's 2026-10-22 rulings): 31 macos commits above the gate-9 base `e474a137`. The Mac line was declared
free of engine defects at gate 9; this pair tests that it stays so. Host: macOS 26.6, gfortran 16.1.0 +
Apple clang 21, MATLAB R2024a arm64, XQuartz. Run 2026-10-22, clean clone. Baseline: gate-9 `ad4161fa` / Linux `0f9983eb`.

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release gfortran — OK, arm64, 189 api syms — OK — ifx n/a
- 3a load (165) — 165/165 ok, 1 ulp, **the 15 GridFile movers present with Linux's numbers** — 160/160 — the GridFile movers (expected)
- 3b cmd — 118 journals, **117 pass**, 2 must-fail red by design, **1 fail = COAT_legacy near-zero pin** — 117/117 — 1 (a round-off re-pin, below)
- 4 mmacos fast — **656 pass / 0 fail** (94 classes) — all pass — none
- 5 pymacos — **6962 pass / 0 fail / 3 skipped(darwin)** — pass — none (the new side-tree guard passes; save-quiet skips (ee))
- 6 /png + /xw — renders clean — pass — none
- 7 self-containment — engine CLI (Cassegrain RMS OPD 1.293491396D-12) + pymacos import, HOME away — — — none

## (f) GridFile-beside-the-deck — the 15 movers confirmed on arm64, Linux's numbers
The load A/B gate-9 `ad4161fa` → this pair moves exactly the 15 grid decks of
`ab_gridfile_moved_0f9983eb.txt` (they were silently loading without their grids; now the grid beside the
deck is read). arm64 reads the SAME numbers as Linux:
- `SegDemo3.in` npass 49318→**49286, lost 32, RMS 6.657848853D+00** (6.66 mm) == Linux.
- `SegDemo3data.in` →**49284, lost 34, RMS 1.425715627D-01** (0.143 mm) == Linux.
- `tst_ZGD2.in` RMS →1.861175691D-04, `tst_eol_grid_lf.in`/`_crlf.in` →3.793426102D-05,
  `tst_ZGDdata2.in`/`tst_short_zerncoef.in`/`FFSegDemoData.in` →their grid-found RMS, the four NSRef decks
  →4.989388559D-09 — all == Linux to every printed figure.
- `SegDemo.in` the one exception, a near-zero: arm64 `6.179096053D-12` vs Linux `6.179211807D-12` (6 figures
  shared, a numerical-zero `D-12` — round-off within the shared-digit policy, not a mover discrepancy).
So (f) lands identically on arm64: no number differs beyond the one near-zero. **Nothing to send.**

## (aaa) LOH_tokens, and the other cleanup journals — all green on arm64
- **LOH_tokens.jou PASS** — the 11- and 100-token `!` lines come back whole; the 142-token line is refused
  by name (`** Command line refused (!): 142 tokens, at most 128 on a line`) and the next line runs. The
  refusal text and the next-line-read both read correctly on macOS (CC's check).
- **GRIDFILE_beside.jou, PROPTYPE_numeric.jou, DOPD_identical.jou** PASS; ZABERR_ZRM_DOPD green.

## cmd — the one fail: COAT_legacy.jou (a near-zero re-pin, NOT an engine miss)
The new (zz) coating journal pins the Cassegrain near-zero RMS `RMS OPD error is\s+1\.2937\d*D-12`; arm64
reads **`1.293491396D-12`** (Linux `1.293791606D-12`) — the recurring `1e-12`-class Cassegrain OPD where the
two architectures share ~4 figures. Identical class to PGD_reopen_units / WLENS_SWL / ZABERR, which CC already
widened. **Fix: widen the pin to `1\.293\d*D-12` (or a magnitude bound), exactly as the others.** This is the
one thing for CC. No (ff)-class engine miss remains — the Mac line stays engine-defect-free.

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate10_*.log; cmd transcripts in ~/mac_gate_b7556b81/macos/cli_tests/records/cmd_gfortran_b7556b81_transcripts/.
- The new pymacos `test_side_tree_is_loaded.py` passes; `test_save_new_file_quiet.py` skips on darwin (the (ee)
  unit-6-across-FFI class — the api save_rx stray-command text is Fortran stdout, lost across the FFI on macOS).
- No engine defect. The one cmd red is a Linux-minted near-zero pin missing the arm64-aware widening.
