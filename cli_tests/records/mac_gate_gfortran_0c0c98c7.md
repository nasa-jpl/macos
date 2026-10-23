# Mac release gate (11) — gfortran, Apple Silicon — the (ee) message buffer

Pair: macos `0c0c98c7` + MACOS_resources `0b39b0a`. Engine changed: TO's round 2, the **api message buffer**
(`emsg_mod`) that resolves PLAN_CONSOLIDATION 1.2 (ee) — the root cause found at gate 4: Fortran unit 6 does
not reach fd 1 across the FFI on macOS, so the embedded engine's warnings/refusals/cause-lines were invisible
to the bindings. 228 print sites now also write to a buffer the api exposes (`messages_count/get/clear`).
Host: macOS 26.6, gfortran 16.1.0 + Apple clang 21, MATLAB R2024a arm64, XQuartz. Run 2026-10-23, clean clone.
Baseline: gate-10 `b7556b81` / Linux `492998db`.

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release gfortran — OK, arm64, **192 api syms** (the 3 new messages_* entries) — OK — ifx n/a
- 3a load (165) — 165/165 ok, 1 ulp, **0 movers** vs gate-10 — 160/160 — none
- 3b cmd — 118 journals, **118 pass, 0 genuine fail**, 2 must-fail red by design — 118/118 — **none** (COAT_legacy bounded now)
- 4 mmacos fast — **659 pass / 0 fail** (94 classes; the 3 ex-guarded Rx classes now run+pass) — 695/0 — count differs by class selection only
- 5 pymacos — **6971 pass / 0 fail / 1 skipped** (only the capfd parity leg, by design) — 7014/0/1 — count by collection; 0 fail both
- 6 /png + /xw — renders clean — pass — none
- 7 self-containment — engine CLI (RMS OPD 1.293491396D-12) + pymacos import, HOME away — — — none

## (ee) RESOLVED — the headline: 6 MATLAB + 2 pytest formerly-guarded tests now PASS on the Mac
The message buffer lifts every Darwin guard that existed only because unit-6 text did not cross the FFI.
All now PASS on arm64 (verified each directly):
- **MATLAB (6):** `tRxRefuse` 1/0, `tRxShortCoef` 5/0, `tRxLineEndings` 1/0 (the nested-`matlab -batch`
  subprocess tests — a child echoes the buffer tagged `MSG|`), `tSaveNewFileQuiet` 1/0 (the (v) save_rx
  stray-command), `tMessages` 3/0 (the new buffer test), and `tCalibMatrix/trial-loses-chief` (18/0 — the
  cause line "N of K trials lost the chief ray" now reaches the test through the buffer; the `if ~ismac` gone).
- **pytest (2):** `test_rx_parse_warnings` (both legs: look-alike header key, OptAsph-above-VarElt) and
  `test_save_new_file_quiet` — `skipif darwin` gone, both PASS.
- mmacos fast rose 656→659 (the 3 ex-filtered Rx classes now run); the Rx classes are no longer in the
  Filtered-by-assumption list.

Only one darwin skip remains, **by design** (CC said not to report it): `test_messages.py`'s capfd parity
leg `test_the_buffer_holds_exactly_what_was_printed` — it compares the buffer against captured stdout, which
macOS cannot supply. And `tCalibStop/test_opt_tol_keyword_reaches_the_stop` stays `assumeTrue` — that is the
arm64 LM *trajectory* from gate 6 (floor at 29, not a buffer issue), as CC noted.

## Pins — all green on arm64
- **COAT_legacy.jou PASSES** — the Cassegrain collimated RMS `1.293491396D-12` is bounded `< 1e-11` now, in
  all seven journals that carried that numerical zero. The gate-10 red is closed.
- The post-tilt `3.566438007D-11` in ZABERR_ZRM_DOPD / DOPD_identical stays at five figures (a real response,
  not a zero) and reads correctly on arm64 — it did NOT move, so no finding there.

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate11_*.log; cmd transcripts in ~/mac_gate_0c0c98c7/macos/cli_tests/records/cmd_gfortran_0c0c98c7_transcripts/.
- No engine defect. (ff) was closed at gate 9; (ee) is closed here — **the Mac line now has no engine miss AND
  no Darwin-guarded test left except the one by-design capfd parity skip and the tCalibStop trajectory leg.**
  The root cause I reported at gate 4 (unit-6-across-FFI) is fixed at its source.
