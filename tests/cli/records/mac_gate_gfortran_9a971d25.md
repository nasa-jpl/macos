# Mac release gate (7) — gfortran, Apple Silicon

Pair: macos `9a971d25` + MACOS_resources `bcf351c`. **The engine is byte-identical
to gate 6** — `git diff 814082d4..9a971d25 -- macos_f90 CMakeLists.txt` is empty,
so the binary, load and ab legs are a formality (confirmed below). What is new
since gate 6: the pair-6 folds (already confirmed on arm64 at gate 6), the manual
work (Appendix E Zernike orderings, figures 23/30/31, section command blocks), and
on the resources side **TO's `pymacos/tests/test_zernike_orderings.py` — 14 new
pymacos tests** that probe the engine's sag one Zernike mode at a time (291
(type,mode) pairs) and fit (n,m)/sign/normalization. The mmacos change is help
text only (`get/set_elt_zrn_type`, `get_elt_zrn` ids 1–11). Host: macOS 26.6,
gfortran 16.1.0 + Apple clang 21, MATLAB R2024a arm64, XQuartz. Run 2026-10-20,
clean clone. Baseline: Linux `fcc648c5` (load).

## Counts (check — Mac — Linux — differs?)
- 1 clean clone — both at SHA — — — no
- 2 build release gfortran — OK, arm64, 189 api syms — OK — ifx n/a; engine == gate 6 (no macos_f90 diff)
- 3a load (163) — 163/163 ok, rt 156/7 n/a, 1 ulp, **0 movers** vs fcc648c5 — 160/160 — 2 edge-ray only
- 3b cmd — 108 journals, **107 pass**, 2 must-fail red by design, 1 fail (COORD_ACOOR = (ff)) — 103/103 — 1 ((ff), TO's)
- 4 mmacos fast — **656 pass / 0 fail** (94 classes; the Darwin/env filters unchanged from gate 6) — all pass — none
- 5 pymacos — **6955 pass / 0 fail / 2 skipped(darwin)** (6941 + the 14 Zernike tests) — pass — none
- 6 /png + /xw — renders clean (spot diagram legible, white bg) — pass — none
- 7 self-containment — engine CLI (Cassegrain RMS OPD 1.293491396D-12) + pymacos import, HOME away — — — none

## The one new thing worth a line — the 14 Zernike-ordering tests PASS on arm64
`test_zernike_orderings.py` 14/14. The `ZerntoMon*` sag / sign / normalization fits
(ANSI, Malacara, BornWolf, Noll, Fringe, ExtFringe and the Norm* variants) that TO
measured on Linux reproduce on arm64 — **no cross-architecture `ZerntoMon*` difference**
(CC flagged that a failure here would be a surfsub.F finding, not a re-pin; there is none).
The full pymacos suite rose 6941 → 6955 with no regression.

## cmd — the single fail (unchanged from gate 6)
- **COORD_ACOOR** — the (ff) gfortran-16 arm64 `-O1+` LocalCoord miscompile; zLocal still
  the 76.5°-about-x rotation. LocalCoord untouched (engine == gate 6). Pins the Linux frame. TO's.
All of gate 6's folds remain green here (BINARY_CENTER `wc -c`, PGD_reopen_units `1\.293\d*D-12`,
and the five bindings re-pins); nothing regressed.

## Findings TO's Zernike probe raised (PLAN_CONSOLIDATION 1.2 (tt)–(ww)) — not seen differently on arm64
(tt) ExtFringe orders 11–14 inert, (uu) NormHex unloadable, (vv) NormAnnularNoll ignores its
ratio, (ww) pymacos's own 9/45 veneer limit vs the api's 11/66. All await Dave's rulings; the
14 passing tests exercise what the manual's Appendix E measured, and arm64 shows them the same
as Linux (nothing to add on these from the Mac).

## Not run / notes
- ifx: gfortran only on arm64.
- Raw logs /tmp/mac_gate7_*.log; cmd transcripts in the gate clone
  ~/mac_gate_9a971d25/macos/cli_tests/records/cmd_gfortran_9a971d25_transcripts/.
- No engine *defect* except (ff). The engine is byte-identical to gate 6, so this gate is really a
  confirmation that the manual/Zernike/binding-veneer additions land clean on Apple Silicon.
