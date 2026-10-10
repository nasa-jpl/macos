# PLAN_CONSOLIDATION — close the open items, test the CLI, finish the documentation, build in nothing new

Dave, 2026-10-08: "step back and pick up the PLAN items, including delayed
fixes and features (section 0) and a full CLI test suite, along with
documentation completion and updates (manual and command reference) — and
be sure not to build in any more problems."

This file sequences the work; it does not replace `PLAN.md` (the record
of every item, 123 open checkboxes) or `PLAN_DESIGN_LAYER.md`.  Every item
below points at its PLAN section.  Check items off HERE as they land and
in `PLAN.md` where they live.

## 0. The rules that keep this from making new problems

Each one is a lesson paid for this year (the memory and `macos_f90/CLAUDE.md`
carry the incidents).  They are the definition of done for every item below.

1. **A fix ships with a gate, and the gate has a must-fail leg** run against
   the pre-fix binary or a deliberately broken input.  A gate that cannot
   fail is not a gate (the grating-aperture, s/p-sign and dcdx cases).
2. **Shared parsers get gated on every surface they serve:** the CLI, mmacos
   and pymacos share `msmacosio.inc` / `iosub.inc` / `macos_cmd_loop.inc`;
   a CLI-only fix is a half fix.  The CLI suite (§2) is the missing surface.
3. **A change to a root pick, a surface routine or a reference convention
   runs the corpus A/B first** (`fs_fix/scan_engine` pattern: every deck
   loads on both binaries; RMS / P-V / nPass / lost identical to 10 digits
   except the decks the change is meant to move, each explained).  Element
   indices come from the ENGINE, never from text-parsing `.in` files.
4. **Units changes get a consumer audit on a factor ≠ 1 fixture** (mm deck),
   never judged on a metre deck.
5. **Guards warn once per run and count; they do not error** on behaviour a
   long-standing tool tolerated.  A hard error needs Dave's sign-off.
6. **User options never live in `elt_mod_init_vars`** (the glass catalog,
   `spcOption`); per-load resets are prescription state only.
7. **Both build trees** (ifx and gfortran) rebuild for every engine change;
   gfortran's stricter checks have found real bugs; the mex relinks only
   when no other lane's MATLAB is running.
8. **Fixed-form lines stay under 72 columns**; >132 truncates silently.
9. **No `--amend`, no force, no edit of a running batch sequence; push only
   on Dave's word**; every number quoted in a report or deck comes from a
   committed run record.
10. **Freeze on new features on `dev-candidate` while this plan runs**, except
    what an item below needs.  New capability requests go into `PLAN.md`
    §4 with a date and wait.

## 1. The open items, classed (pointers into PLAN.md and CLAUDE.md)

### 1.1 Engine defects, section 0 (PLAN §0) — do first
- [ ] IRIS `save_rx` → reload SIGSEGV on the real ZrnGrData grids (§0 line 29).
- [ ] lensarr trace-time overrun (§0 line 42).
- [ ] `trace(26)` then `trace(27)` wrong first OPD on the jwst zoom deck (§0 line 51) — likely the same family as the restarted-trace `PrevNonSeg` fix (CLAUDE 2026-10-03); check that first.
- [ ] `ApStop` and the other unsaved variables in SAVE (§0 line 55).
- [ ] `ZernCoef=` / `MonCoef=` / `FFZernCoef=` short-line READ guards via `ReadRealsPad` (§0 line 270; the `AsphCoef=` fix is the template).
- [ ] `ifLNsrf` root metric in RefSrf / ObsSrf / PolElt / IntSrf (§0 line 271; `LNsrfRoot` exists; the three-case rule in CLAUDE.md applies).
- [ ] `LUseChfRayIfOK` global default (§0.x) — Dave's ruling needed (segmented decks want the chief).

### 1.1a Found by the CLI load gate on day one (2026-10-09, core corpus, ifx and gfortran)
- [ ] **SAVE → load → SAVE is a 1-ulp 2-cycle on tilted unit vectors** (`psiElt`, `xGrid`/`yGrid`): 9 of 91 core decks (CornerCube, HOEExample, opt_example ×3, Rx_Mask_Parabolas_glb, Rx_CornerCube, pymacos e5hex1 …) never reach a fixed point — the unitise-at-load of a printed 17-digit vector alternates between two neighbours, the same disease `OrthoSrcFrame` cured for the source frame with the dead band (CLAUDE.md "Re-traces are IDEMPOTENT").  The CLAUDE.md comments entry already calls the psiElt wobble "pre-existing".  Engine item for TO after the five in the brief: the same `DeadBandUlp` rule at the psiElt/xGrid unitise on load; gate = the suite's `roundtrip` column going `ulp` → `same` on those nine.
- [x] **`docs/macos-manual/examples/SegDemo.in` is refused by the validator** — FIXED 2026-10-09 (the stray heading removed; the extractor scripts no longer exist in the tree, so the file is the fix; the other 22 examples end on a keyword line and load).  Original:: the Appendix-A extractor left the next section's heading (`A.6 Near-Field Propagation Example`) after the deck, and a non-keyword line after a blank following a `Tout=` block reads as "blank line inside multi-row block".  Docs lane (CCMac): the extractor strips trailing text; check the other 22 manual examples the same way (they load today).  Validator question for later: should a line with no `=` and no digit end the block instead?
- [ ] **Two fixtures trace to NaN OPD**: `tst_save_keys.in` (the lensarr item, brief item 3) and `tst_block_comment.in` (the same deck with comment blocks added, so the same NaN — one item, not two).
- [ ] **`pymacos/tests/Rx/Grating_example_001.in` passes 7303 rays on ifx and 7295 on gfortran** (RMS OPD 177.168 vs 177.155): eight rays' pass/fail depends on the compiler's arithmetic — rays on a knife edge of the concave grating's aperture (the deck CLAUDE.md notes was "never ray-compared").  Not a defect by itself (a ray exactly on an edge is a coin toss) but a fixture that cannot anchor a cross-compiler record; either find the edge and move it, or accept it as the one known compiler-sensitive deck.  (`Rx_Coro_FPM_Zern_vortex_oversized_noLyot.in` differs at 5e-6 of a 1e-12 OPD: noise.)
- [ ] `ZGD_test_files/SegDemo3broken.in` refuses by design (a broken-deck fixture) — move it to `cli_tests/must_fail/` as a third permanent negative control.

### 1.2 Deferred engine items recorded only in `macos_f90/CLAUDE.md`
- [ ] `NSRefractor` passes a null grid frame to `GridSrf` (piston-only figure on refractive NS grids).
- [ ] `ZernTypeL = 11` (ExtFringe) has no `ZerntoMon` converter: add one or an error.
- [ ] FEX radius doubles as the far-field propagation distance: split the surface leg from the plane leg if a physical-optics leg ever needs it (decide: document or implement).
- [ ] `Refractor`'s coated branch oblique radiometry vs ray density (open audit question).
- [ ] `RayFailElt` stamps `nElt+1` for obscuration instead of the clipping element.
- [ ] Multi-value `READ(VALUE,*)` keywords other than AsphCoef/AnaCoef (same as 1.1).

### 1.3 Optimizer (PLAN §3) — the urgent two, then the matrix
- [ ] `WFE_ZMODE_TARGET` and `OPL_TARGET`: implement or remove (§3.1, marked urgent).
- [ ] `SPOT_TARGET` reporting 0 iterations (§3.1).
- [ ] Keyword fixes §3.2 (`OptTol`, `OptMxItrs` prefix, `OptAsph` naming, the MBFile6 warning).
- [ ] Control flow §3.3 (convergence check every iteration; one `nls_optim_dvr` call site; `quit` after CALIB; `dopt_init_vars` audit).
- [ ] The target × DOF regression matrix §3.4 (12 cells) — each cell is a CLI-suite case (§2).
- [ ] Asphere + Zernike on one mirror as a surface option; EPFIX (§3.1, Dave 2026-10-05) — features; they wait for the freeze to lift unless Dave rules otherwise.

### 1.4 Tools and runners (outside the engine)
- [ ] The gauge runner's shared `ctx.msk` includes the reference-only annulus (tg96_run `arm_setup_`); its pixel-rms metrics are diluted. Fix in the runner after the talk's record is frozen; actuator-space rows are unaffected.
- [ ] `dyson_ladder` meniscus seeding (fixed by TO 2026-10-08; confirm the gate).
- [ ] Gauge-deck migration into `templates/40_benches/gauge_deck` (memory: planned 2026-09-15).

### 1.5 Documentation (PLAN §6 + the command reference)
- [ ] Manual chapters `docs/macos-manual/src/00..09`: math re-entry and the three figures still pending (memory `project_manual_markdown_migration`); §6.2 new chapters (polarization, physical-optics medium awareness, the prescription validator, comments, SAVE round-trip, STOP on segments, FEX conventions) — each written FROM the CLAUDE.md entry it documents, then the entry trimmed to a pointer.
- [ ] Command reference `docs/macos-manual/cmdref/`: Phase A catalog exists; Phase B = every command with syntax, defaults, conventions table, one example journal and a pointer to its CLI-suite case (§2.3 makes the example and the test the same file).
- [ ] `templates/00_INDEX.md` and each template README current (`tma_longslit`, `gauge_deck`, dyson5).
- [ ] §6.3 worked examples: the sensitivity workflow (§1.3), one telescope design, one spectrometer, one gauge bench — each a runner that exists today.

## 2. The CLI test suite (PLAN §1.1, built here)

The CLI is the one surface with no regression suite; the bindings have 570
fast tests and pymacos 6601.  Three layers, one script
`macos/cli_tests/run_cli_tests.sh [load|cmd|ab|full]`, both compilers, the
results written as records so a diff is the report.

### 2.1 Load gate (`load`) — every deck in the corpus
`ZGD_test_files/` (48), `MACOS_resources/GMI/test_ff/` (21), the mmacos and
pymacos `tests/Rx/` fixtures, the templates' emitted decks, the manual's
examples.  For each: LOAD, `opd nElt`, SAVE, reload, SAVE again; record
(RMS, P-V, nPass, lost, byte-identity of the second round trip).  The
`fs_fix/scan_engine` harness is the seed; it already drove 406 decks.
Must-fail: a deck with a known short `ZernCoef=` line, until 1.1 lands.

### 2.2 Command gates (`cmd`) — one journal per cmdref entry
pty-driven (`scratchpad/spot_cli_gate.py` / `fex_cli_gate.py` pattern:
readline needs a tty; the model-size prompt comes first; sub-prompts are
answered by the journal).  Each `cli_tests/cmd/<COMMAND>.jou` carries its
asserted lines (`# EXPECT: <regex>`) and its must-fail twin where the
command has a known failure mode.  Priority order: the commands the
bindings cannot reach (MOD, VALIDATE, JOURNAL, LOG, the interactive
prompts), then STOP/FEX/SXP/SPOT/OPD/PERTURB/CALIB, then the rest of the
catalog.  The §3.4 optimizer matrix lives here as 12 journals.

### 2.3 A/B gate (`ab`) — the previous binary against this one
Builds the previous tagged commit into `build_ab/` and diffs the `load`
records; any changed deck must be named in the commit message.  Runs
before any engine commit touches `surfsub.F`, `elemsub.F`, `tracesub.F`,
`propsub.F`, `msmacosio.inc`, `iosub.inc` or `macos_cmd_loop.inc`.

### 2.4 Non-vacuity
Every gate carries one case that fails on a pre-fix binary (kept in
`cli_tests/prefix/` as a rebuilt commit hash) or on a broken input.  A
suite with no red leg on record is not accepted.

## 3. Sequencing (four phases; dates assume a start 2026-10-09)

| phase | weeks | engine lane (TO) | docs lane (CCMac) | CC | gate to exit |
|---|---|---|---|---|---|
| **1 Harness + §0** | 10-09 → 10-16 | §2.1 load gate on the corpus, both compilers; §1.1 defects in the order listed, each with its must-fail leg | cmdref Phase B skeleton (one template page; the journal-as-example convention) | reviews each fix against §0's rules; corpus A/B on every surface change | load gate green on 100 % of the corpus; §1.1 closed or Dave-ruled |
| **2 Commands + optimizer** | 10-16 → 10-30 | §2.2 command gates by priority; §1.3 urgent two; §3.2–3.3; §1.2 deferred items | cmdref pages written as the journals land (same file); manual chapters from CLAUDE.md entries | A/B gate §2.3 wired; weekly full-suite run (fast mmacos + pymacos + CLI) | every cmdref command has a journal; §3.4 matrix green |
| **3 Documentation** | 10-30 → 11-13 | the §1.2 items that remain; the gauge runner mask; the design-layer PLAN's own §0 | manual math and figures; §6.2 chapters; §6.3 worked examples run from the templates; CLAUDE.md entries trimmed to pointers once documented | reads every page against the running tool (doc = test) | manual builds; each example reproduces its record |
| **4 Release gate** | 11-13 → 11-20 | both compilers clean; Windows configure smoke (memory: Luis's CMake incidents); self-containment sweep (no `addpath` outside the repos; suite on a clean checkout) | release notes from the gates' records | NPSOL-history question to Andy (memory: `project_npsol_history_exposure`); promotion `dev-candidate` → `dev`/`opt-dev` per PLAN §9 | fast + CLI + pymacos suites green on a clean clone; Dave's word to push |

Lanes are the current ones (TO engine; CCMac docs; CC review, gates and
this plan); Dave rules the open questions below and signs off each phase.

## 4. Rulings needed from Dave

1. Feature freeze on `dev-candidate` for the four phases (rule 10) — yes/no, and the exceptions (asphere+Zernike surface; EPFIX; the 1.5k template run in flight).
2. `LUseChfRayIfOK` global default (PLAN §0.x): chief for segmented decks?
3. IRIS real-grid reload: fix the engine, or document the grids as unsupported for SAVE round-trip?
4. The CLI suite's home: `macos/cli_tests/` (engine repo, beside `ZGD_test_files/`) is proposed.
5. Mac port (memory `project_mac_port`): in phase 4 or after.

## 4a. Rulings taken so far
- 2026-10-09: Dave started phase 1 ("take it slow and careful"); TO on §1.1 per `BRIEF_to_consolidation_p1.md`, CC on §2.1.  The suite's home is `macos/cli_tests/` (ruling 4, by use).  Rulings 1–3 and 5 still open; IRIS and `LUseChfRayIfOK` not started.
- 2026-10-09 (rule 5): **Dave approved TO's proposal for the no-natural-zero reads** — the ~113 multi-value geometry / integer-list keywords (psiElt, VptElt, RptElt, ChfRayDir, xGrid, TElt rows, ApVec/ObsVec groups, ZernModes, EltGrp …) and ~90 scalars that today abort the Fortran runtime on a short line or a non-numeric token will instead FAIL THE LOAD cleanly: one message naming the key and the element, `LOAD_SUCCESS` false, the host alive.  Zero-padding is wrong there (a zero normal, a zero-radius aperture, element 0).  Item 1b in TO's lane, after item 2.
- 2026-10-09: **Dave approved line-ending normalization at LOAD** (item 5's fix): the reader shared by the CLI, mmacos and pymacos accepts LF, CRLF and CR-only decks and strips the trailing CR from every record — CRLF (Windows, `autocrlf` checkouts) is the primary case and carries the GridFile-tab class of bug on filename/name values; CR-only (classic Mac; 13 corpus decks, all failing on ifx) is the crash.  SAVE stays native (LF here, CRLF on Windows) — it is the normalization by construction.  Gate: one deck in three endings → identical OPD on both compilers and a byte-identical second SAVE; the judge names an endings-only difference `eol`.  The same gate runs in the phase-4 Windows smoke.

## 5. Bookkeeping

- `CURRENT_SLICE.md` carries the phase in flight; this file carries the checklist; `PLAN.md` stays the record.
- Every closed item gets its line in the nested `CLAUDE.md` turned into a one-line pointer to the manual page that now documents it (the gotcha file shrinks as the manual grows).
- Weekly, CC runs the three suites and writes one line per suite here with the date and counts.
- 2026-10-09 night, CLI load gate, FULL committed corpus (1258 decks), binaries at 52b1409: ifx and gfortran both 1255 ok / 3 timeout_load (three `pupilsim_redo_lens_*_deck.in` run artifacts under tg_psi_dm96_oap/runs — model-512 decks at model 256? — to classify); round trips 1117 same / **129 at 1 ulp** / **9 differ, all `dyson5_t5f_*_e2e.in` joined decks** (a real SAVE round-trip difference, TO's domain: what the join emits that SAVE rewrites); NaN OPD 3 (`tst_save_keys`, `tst_block_comment`, `dyson5_t3_t600_off10_s3`); 5 decks differ between compilers (3 bench_ctb at 1e-9 relative, the grating edge rays, a 1e-12 OPD).  Records `cli_tests/records/corpus_{ifx,gfortran}_52b1409.csv`.  The driver's pty-fd leak (stopped the first attempt at deck 511) fixed en route.
- 2026-10-09 CLI load gate, core corpus (91), binaries at 8cf28c3 (clean worktree): ifx 89 ok / 2 load_fail (SegDemo3broken by design; SegDemo.in extractor artifact), round trip 80 same / 9 ulp; gfortran identical counts; 2 decks differ between compilers (Grating_example_001 8 edge rays; a 1e-12 OPD at noise).  Must-fail decks: both refused on both compilers.  Records `cli_tests/records/load_{ifx,gfortran}_8cf28c3.csv`.
