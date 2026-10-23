# mmacos manual, Part I: report for Dave (TO, 2026-10-23)

The brief's form (BRIEF_to_mmacos_manual.md, Report): the SHAs, which chapters have generated
artefacts, the page count, the open questions in one list, and the two defaults stated as defaults.
All commits are on `dev-candidate` in `~/dev/macos`; CCL pushed each one after reviewing it, on your
standing word. Nothing in MACOS_resources was changed: every template ran unchanged from a copy, or
from a copy with a single block cut where noted, and the library's records were checksummed before
and after.

## Part I: seven chapters, all generated, 45 pages

`make all pdf` builds `docs/mmacos-manual/build/mmacosMan.pdf` (45 pages) with the MACOS manual's
pipeline, reused by reference. Every listing, figure and quoted number comes from a chapter driver's
run (`tools/transcripts/*.m`, `make transcripts`). Each run is gated: pins hold to 5 significant
figures on a signal or a bound on a numerical zero, and listings may drift only within the printed
precision. A failed gate writes nothing.

| Ch. | Title | What runs | Chapter SHA |
|---|---|---|---|
| 1 | A Telescope | e2e stage 1, 4 m f/18 folded Korsch TMA (14 min) | b0e3c18a |
| 2 | Instruments | e2e stage 2 Offner imager (49 min); `coro_planet_demo` surfaces page; `bench_layout`; ImageDemo counter-example | ad35e34c |
| 3 | Segmentation | e2e stage 3, pie (7) and hex2 (19) (87 s) | 77ee8bcd |
| 4 | Sensitivities | e2e stage 4, dwdx / dwdz / dwdgrid (11 min) | f79c2580 |
| 5 | Laser Metrology | e2e stage 5 truss and live layout search (52 min) | 44a51c35 |
| 6 | Closing the Loop | e2e stages 6 and 7, plus the chapter's own OPD-reference check | 81be4f34 |
| 7 | An Error Budget | the e5hex1 budget template (off the thread, see question 1) | 00c45a83 |

Pipeline and gate commits along the way: b0e3c18a (the pipeline), 06497da7 (bounds by the
finite-difference noise floor), dbc2aea2 (a bound pin fails closed; its test is red before the fix
and green after). Plan findings: 9b6de6c0, fecf9ecf, 8f373f32, 7f51397a, 1b9f8ff4. The thread's
long stages ran twice where it mattered (chapters 1, 2, 4 and 5); each pair agreed exactly.

## The two defaults (overturn either in one line)

1. **The thread is e2e** (`templates/80_end_to_end/e2e`, the 4 m f/18 folded Korsch, all seven
   stages), not e2e2, the 3 m f/20 design-only flow, which chapter 1 points to.
2. **The forty uncalled helpers go in one table in Part II**, not deleted. Part II is not written
   yet; it is next.

## Open questions, in one list

1. **The thread's error budget.** e2e has no budget stage, and its group table (a sigma per error
   line, the control set, the noise) is your choice. Chapter 7 therefore runs the library's e5hex1
   budget, the same seven-segment class. Should an e2e budget spec be written, and by whom?
2. **(kkk) The OPD reference for sensitivity harvests.** The engine default became 'chief' on
   2026-10-10; `run_sensitivities` harvests at 'mean'; `run_compare` takes the default. So the e2e
   engine-against-model gate reads w ≈ 1 against the 2.6e-3 it read in August. Measured: segment 1
   piston +100 nm gives 188.9 nm rms under 'chief' (178.4 nm of it piston), 61.98 under 'mean', and
   61.98 in the linear column. Which reference should the harvest use? (CCL's reading of the
   one-reference rule: the engine's, 'chief'.)
3. **(bbb) Regenerate e2e's committed outputs?** The committed stage-1 products are stale against a
   live run (conics, ±1′ WFE, the hole 0.095 → 0.232 m), and the staleness reaches every later stage.
   (fff) is the live pie clipping 0 gap rays against the committed 376. (ggg) is the committed
   stage-5 report, made on the `opt_preset` path, which the committed script does not take; the live
   search finds a different center-hexagon pattern (MET-only WEM 3.877 against 3.707).
4. **(ccc), (ddd) The two coronagraph templates.** `coro_walkthrough` shows no dark hole, and
   00_INDEX names it as the T3 coronagraph entry. `coro_planet_demo`'s planet page shows no planet
   PSF. Nothing either one reads has changed in git since June, so the engine or the dispatcher moved
   under them. Neither is diagnosed.
5. **(eee) ImageDemo.in has no `BaseUnits`.** The camera template should declare `BaseUnits= mm`.
6. **(hhh), (iii), (jjj) Three figure defects in the library** (the figure-producers-rot class):
   - `run_met`'s metric chart drops its last bar under the log axis;
   - `met_view`'s hub inset is labeled mm whatever the deck's units;
   - `run_budget`'s chart leaves out the rows that make its total.
   The manual works around each one and says so in the caption. Fix them?
7. **The design_report "EFL" label.** It is the local plate scale at the field point: 72.32 m
   against SYSPROP's 72.14 m in chapter 1. A one-word label fix in my lane, if you agree.
8. **The authorship / acknowledgment block** in the front matter is the MACOS manual's, as a
   placeholder.
9. **Carried from the bindings round, still yours:** the `elt_csys` update flag, and the perturb
   frame (PLAN_CONSOLIDATION 1.4a).

The full measurements behind each finding are in `docs/mmacos-manual/OPEN_QUESTIONS.md` and in
PLAN_CONSOLIDATION 1.2, (bbb) to (kkk).

## The findings by letter (PLAN_CONSOLIDATION 1.2; measured, none fixed)

- **(bbb)** e2e's committed stage-1 outputs are stale against a live run, and the staleness reaches every later stage.
- **(ccc)** `coro_walkthrough` no longer shows a coronagraph; 00_INDEX names it as the T3 entry.
- **(ddd)** `coro_planet_demo`'s planet page shows no planet PSF; its surfaces page is sound.
- **(eee)** ImageDemo.in declares no `BaseUnits` (CBM = 1, so `dpup_m` reads 380, the deck's mm).
- **(fff)** A live e2e stage 3 clips 0 gap rays on the pie, against the committed 376; source rays 12128 against 12520.
- **(ggg)** The committed `s5_report.txt` came from `opt_preset`, not from the committed script's search.
- **(hhh)** `run_met`'s metric chart drops the "edge+MET opt" bar under its log-axis floor.
- **(iii)** `met_view`'s hub-plane inset is labeled mm whatever the deck's units (met_view.m:199-200).
- **(jjj)** `run_budget`'s per-group chart leaves out the actuator and sensing rows that make its total.
- **(kkk)** `run_compare` mixes OPD references: the engine's 'chief' against the harvest's 'mean' (188.9 against 61.98 nm on one poke).

## One thing I got wrong, caught before commit

The pipeline's bound pins printed their bound without checking it. In chapter 6's first run that
produced "w < 0.01" for a measured 0.945. Nothing false was committed, and every committed bound pin
holds its rule under the fixed checker. The helpers now fail closed, and the regen script checks
every bound on the run that records it (dbc2aea2).
