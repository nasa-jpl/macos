# Mac release gate (14) — gfortran, Apple Silicon — the pre-consolidation release pair (part A)

Pair: macos `0984cd34` + MACOS_resources `ef2d7f9a` (`dev-candidate`, nasa-jpl remotes). This is part A of
the release gate (BRIEF from CCL 2026-11-02): the full clean-clone gate on the pair the release is cut from,
before Andy's repo consolidation. Part B (re-clone the rewritten repo, gate from `main`+`dev`) follows the
rewrite and is the Mac's release sign-off.

Host: macOS 26.7.2 (25G332), gfortran 16.1.0 (Homebrew GCC 16), Apple clang 21.0.0, MATLAB R2024a arm64,
Python 3.13.14 (venv + numpy 2.5.1), XQuartz on /opt/X11. Run 2026-11-02, clean clone into
`~/mac_gate_0984cd34/{macos,MACOS_resources}` (case-sensitive volume). Baseline: gate-11 `0c0c98c7` / the
load record `fcc648c5`.

**Gate numbering:** this is gate **(14)**, not (12). The intervening pairs CC queued — pair 12 (macos
`08b92085` + resources `446bd45`) and pair 13 (macos `c6c83d55`) — were **never gated on this Mac**: no
clone dir, no logs, no record, no commit exists for either (every gate through 11 has a `~/mac_gate_<sha>`
dir; these two do not). Their `NOTE_to_ccmac_pair12/13.md` describe anticipated runs. Per Dave's ruling
(2026-11-02) the unrun pairs consume their numbers; this pair is (14). No records were fabricated for 12/13.

## Counts (check — Mac — Linux (brief) — differs?)
- R1 clean clone — both at pair SHA, 0 dirty, case-sensitive volume — — — no
- R2 engine+CLI build — OK, arm64 `bin/macos` + `libsmacos.a`, **194 api syms** — OK (ifx+gfortran) — ifx n/a on arm64
- R2 mmacos mex — OK, `mmacos.mexmaca64` (direct gfortran link) — OK — none
- R2 self-containment — PASS: no external addpath in src/design; `~/matlab` OFF path (restoredefaultpath) → tDesignTelescope **73/0**; pymacos clone PYTHONPATH only — — none
- R3 /png — **renders clean**: 800×600, white-bg + ink (92.9% white / 0.1% ink / 967 colors), NOT giza-inverted — pass — none
- R3 /xw — not re-exercised (zero giza/plot source change since gate 11; gate 11 confirmed /xw opens under XQuartz) — pass — n/a
- R4 pymacos (PROPER excl.) — **6985 pass / 0 fail / 1 skip** — 7023 (incl. proper_compare) — count by collection; **0 fail both**
- R5 load (165 core) — **165/165 ok, 0 movers** vs fcc648c5, 3 must-fail red; roundtrip 158 same/7 n/a; settle 127 same/30 differ/1 ulp — 160/160 (Mac corpus = 160 + 5 Mac-only) — none
- R5 cmd — **122 pass, 0 genuine fail**, 3 must-fail red (calib_matrix_pin, glasselt_prefix, selftest) — 122 + 3 must-fail — **none**
- R5 mmacos fast — **757 pass / 1 fail / 46 incomplete** (804 enumerated, 4 shards) — 803/0 (TO -j3) — 1 arm64 near-zero pin; filter set differs
- R6 record — this file + MAC_PORT.md Part 2 + 1 Part 1 edit (make_corpus) — — —

## The one mmacos fail — a cross-architecture near-zero, NOT a defect, NOT re-pinned
`tSessionModes/test_fex_follows_the_session_axis`: `verifyEqual(s.xp, s0.xp, 'AbsTol', 0)` comparing the
omitted-mode FEX (`fex()`) against explicit CENTROID (`fex(0)`). All components agree but index 2:
**Actual `1.89100164151503e-24`, Expected `9.50724582198163e-25`** — both numerical dust ~1e-24 against
~1e3-scale vector components, failing an absolute tolerance of exactly 0. Per Dave's rule and
`cli_tests/README.md` ("Cross-architecture pins") this is a **magnitude bound** (both `< 1e-23`, i.e. zero
to ~27 figures below the signal scale), reported as a finding, not re-pinned. Same signature seen in the
2026-10-29 full-suite run. Test-side over-strictness (`AbsTol 0` on a near-zero component); the fix is
test-side (bound index 2), a judgment call for Dave/CCL, not an engine change.

## The 46 incomplete (filtered on arm64, as on Linux) — not failures
By class: tRunCompare 15 (needs ifx `SegMirMaker`, impossible on arm64), tSegmentRx 12, tRunMet 5,
tSegMirMaker 3 (`SegMirMaker` binary not built on this clean clone — assumeTrue skip), tRunSensitivities 3,
tEdgeSensors 3, tBudgetKernel 2, tRunSegmentation / tMet / tCalibStop (the OptTol LM-trajectory leg) 1 each.
The segmentation/MET classes and tRunCompare are filtered out of the fast suite on Linux too.

## What the pair added since gate 11, confirmed landing clean on arm64
GridSrfOrder=3 C1 interp (tGridInterp in CTB_512, not in fast); exact glass-name match + N-SF5/N-F2
(205-glass table); session modes_set CENTROID/LNEG + OBS-across-load (tSessionModes — the one near-zero
above is the only snag); header-key rewrites anchored `(?m)^` (tRxHeaderRewrite, the OptRayGrid=/yGrid=
sibling-clobber I found on 2026-10-29 — green here); real lenses macos.design.Lens + singlet/doublet/triplet
(tLens/tLensSinglet/tLensDoublet/tLensTriplet — pass); the camera family (tCassCamera/tSmallCamera — pass);
e2e stages 9–11 and the sharded runner (`-j N`, see finding 2). api syms 192→194.

## Findings for CCL (tooling portability — neither an engine defect; Part 1 items)
1. **`cli_tests/make_corpus.sh` aborts before writing `corpus_core.txt`** on macOS. `set -euo pipefail` +
   the `is_rx()` trailing `grep -q` returning nonzero (last sorted file has no `nElt`/`Element`) kills the
   script before the line-32 core-list grep, leaving `corpus_core.txt` the stale committed Linux copy (160
   decks, `/home/dcr/...`). Worked around by regenerating core by hand from the resolved
   `corpus_committed.txt` (→ the correct 165 Mac decks). Real fix: `|| true` on the is_rx grep or drop
   pipefail around it; ordering-dependent so it can bite Linux. **Part 1 edit made** (documented in the
   gotchas), diff to CCL for the Linux leg.
2. **Sharded `run_mmacos_tests.sh -j N` clamps to 1 shard on macOS.** The MemAvailable probe is Linux
   `/proc`-style and returns 0 GB on Darwin, so `-j 4` printed `0 GB available -> 1 of 4 shards at once` and
   ran serially (correct results, ~2419 s wall not ~832 s). Fix: a Darwin branch using `sysctl hw.memsize`
   / `vm_stat` (same family as makems.sh's `hw.logicalcpu`). The runner echoes "(re)building
   src/mmacos.mexa64" but builds `.mexmaca64` correctly.

## Already resolved (brief asked; prior Mac work covered it)
`cli_tests/cmd/OBS_POLYGON.jou` was made portable at gate 5 (`54481537`): `perl -i -pe`, not `sed -i`;
passes here with its EXPECT-NOT guard. No Part 1 edit needed.

## Not run / notes
- ifx: gfortran only on arm64 (no Intel target).
- /xw not re-exercised (no graphics source change; X server launches on demand under XQuartz).
- SegMirMaker binary not built on this clone (engine+mex+pymacos only) → tSegMirMaker/tRunCompare skip.
- Raw logs: /tmp/g14_*.log; shard results /tmp/mm_shards.5899/.
- **No engine-source defect found.** Every Mac-only delta is arm64 round-off through a Linux-minted
  zero-tolerance pin (the one tSessionModes near-zero), a compiler/environment filter (ifx-only classes,
  SegMirMaker), or gate-tooling portability (the two findings above).
