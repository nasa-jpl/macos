# NOTE to CCMac — fixes from the manual regeneration and the Linux gate, branch `mr_merge_regen_fixes` (CC, 2026-11-19)

This branch sits on `mr_merge` at `b5b37c34`. It changes MATLAB code, tests, two build defaults,
page specs and records. It changes **no engine source and no generated binding file**, so nothing
needs rebuilding: your `mr_merge` build of the engine, the mex and pymacos is the build to test
with. Read this note from `devnotes/` in your clone. Dave's word is to push the commits as they
land; when the Mac checks below pass, open the pull request with base `mr_merge`.

## Why it exists

Two things were run on Linux on 2026-11-18: every page of the mmacos manual was regenerated, and
the full gate was run on `mr_merge`. Both were run under a stand-in home directory that has no
old `~/dev/MACOS_resources` checkout, which is the state of anyone who clones today. That exposed
what the old checkout had been hiding on Dave's machine.

## The commits

1. **`4fa663a1` Restore the two tuned tail workspaces the DM-gauge bench pages load.**
   `lens_tail.mat` and `redo96_lens_tail.mat` in `templates/40_benches/tg_psi_dm96_oap/` were
   untracked (`*.mat` is ignored) and the fresh clones dropped them. Without `lens_tail.mat`,
   `tg96_run` falls back to a scaled seed and the page's five pins move.
2. **`aaa454fe` Re-pin the pupil_id manual page.** Its deck `tma_onaxis.in` was regenerated on
   10-29, five days after the page's pins were recorded. Three pins move in the fourth digit.
3. **`8d8292fc` Size field-scan apertures from the fan's global positions.**
   `Telescope.realize_apertures` read `macos.draw_rays` output as global coordinates. Those are
   projections on the source axes, which rotate with the field angle, so an off-axis field's
   footprint was displaced by −z·tan(field) and the realized aperture clipped the edge field.
   It now reads `macos.draw_rays3d`. New test
   `tApertureFrame/test_realized_apertures_pass_the_scan_edge_field`. The `rc_onaxis` page is
   re-pinned: its field balance had been silently leaving out the 4.0 arcmin field.
4. **`f7a1ce57` Record the command-line gate on the merged tree.** Linux records only.
5. **`b9e123ef` Run the two CTB DM test classes in a scratch folder.** `tCtbDm` and `tCtbProbe`
   rewrote the tracked `bench_ctb/ctb_dm.in` on every run, which fails the full suite's
   clean-tree guard.
6. **`e04a33e3` Replace paths into the old MACOS_resources checkout.** Code that
   found files through `~/dev/MACOS_resources` now builds the path from its own location: the
   eight `tests/proper_compare` files, the tools `fex_one`, `vh_cass_probe`, `oc_ladder` and
   `run_pol_validation`, 44 template and challenge scripts (their `addpath` lines), the
   coronagraph experiment and walkthrough scripts, 18 run scripts under `40_benches/*/runs`,
   and the mmacos Makefile and CMakeLists smoke-test deck default. Four page specs (`e2e2`,
   `e5_pie`, `coro_walkthrough`, `bench_ifo_dm`) link the repository folders their templates
   climb to. The commit message says 58 scripts and 19 run scripts; 44 and 18 are the measured
   counts. 100 files, 160 lines in, 149 out.

7. **`f6a5ea6d` Commit the DM deck as the current generator writes it.** `bench_ctb/ctb_dm.in`
   gains `GridSrfOrder= 1` on each DM (the engine default, written by `ctb_dm_rx` since 10-31).
   The bench's own functions regenerate this deck in place, so it showed as modified after any
   bench run. Two lines; nothing to check beyond `git status` staying empty after a
   `ctb_dm_jacobian` call without a deck.

8. **`177fe171` Check the CTB deck generator on the source deck the hand decks were made from.**
   `tCtbProp/test_generator_reproduces_hand_decks` failed on `mr_merge` itself. The engine had
   not moved: commit `d3e4770f` (2026-10-29) regenerated the generator's source deck
   `ctb_planar_stageF.in` a day after the pins were set, and that moved the optics (105 of 614
   lines of the generated compact deck), not only the four radii. The test now runs the
   generator on the earlier source, kept as the fixture `tests/Rx/ctb_planar_stageF_2026-08.in`.
   No deck or pin changes. On the Mac: `./run_mmacos_tests.sh tCtbProp`, 8 pass expected.

## What to do on the Mac

No rebuild. `git fetch origin && git switch mr_merge_regen_fixes`, then, with **no**
`~/dev/MACOS_resources` reachable (rename it for the run if you still have one):

1. `./run_mmacos_tests.sh tApertureFrame` — 6 pass expected.
2. `./run_mmacos_tests.sh tCtbDm` and `./run_mmacos_tests.sh tCtbProbe` — 12 and 6 pass, and
   `git status` must be empty after each.
3. The PROPER classes and the overcoat test, one class per process:
   `tProperCompareCassFF` (4), `tProperCompareCassFFAberrations` (6), `tPolContrast` (15),
   `tProperCompareCoroApodizer` (1), `tProperCompareCoroNFprop` (3),
   `tProperCompareCoroDMPhase` (3), `tProperCompareCoroPhase3` (6). On Linux all 23 PROPER
   tests and the overcoat test failed with "Rx file not found" before commit 6, on any home
   directory without the old checkout.
4. `./run_mmacos_tests.sh -j 4 fast`.
5. The manual pages `p2_rc_onaxis`, `p2_pupil_id`, `p2_tg_psi_dm96_oap`, `p2_cophi_dm96`, and
   the four that could not run before: `p2_e2e2`, `p2_e5_pie`, `p2_coro_walkthrough`,
   `p2_bench_ifo_dm` (the last needs the pymacos Python environment; commit 6 looks for
   `app/pymacos/.venv`, then `~/.venv/pymacos`, then `python3`).

Send me any pin that differs, with the digits MATLAB prints.

## What I expect on the Mac

- Steps 1 to 3 pass with the counts given. They do not depend on anything the Mac lacks.
- Step 4: Linux gives 817 pass, 0 fail, 0 incomplete.
- Step 5: on Linux all eight pages pass with the committed pins. `p2_rc_onaxis` prints
  `0.036 0.034 0.028 0.019 0.016 0.035 0.064` waves over the seven fields and is
  diffraction-limited out to 4.0 arcmin. `p2_coro_walkthrough` runs at model size 1024 and
  needed more than 8 GB here.
- The template scripts no longer add `~/dev/MACOS_resources/mmacos/src` to the MATLAB path. If a
  Mac result changes because of that, it was reading the old checkout's source before; tell me
  which script.

## Known, not fixed here

- **The full suite's model-1024 PROPER group cannot run four shards at once in 20 GB.** The
  runner's clamp assumes 4 GB per shard. Run that group at `-j 1`, or set `MM_SHARD_GB=8`.
- **Two `tE2eSensing` grid-DM tests are skipped** on a fresh checkout (they need local stage-11
  products). The sharded summary prints a skipped test under "FAIL"; read the count line.
- **Left alone on purpose:** run logs and `.json` sidecars that record the old absolute paths,
  `devnotes/`, the history sections of the top-level documents, the generator's
  "AUTO-GENERATED by MACOS_resources/..." header (changing it regenerates `mmacos_gen.F` and
  forces a mex relink), and two Fortran comments (an engine rebuild for a comment).
