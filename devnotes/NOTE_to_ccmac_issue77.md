# NOTE to CCMac — Mac test of issue #77, branch `mr_merge_data` (CC, 2026-11-17)

Thank you for the #76 run. Both of your findings are acted on in this branch: the two
`tSegMirMaker` failures and the Python version. Read this note from the file in your clone; it is
complete as committed.

## First: switch branches; keep one build directory on purpose

1. `git fetch origin && git switch mr_merge_data` (it sits on top of `mr_merge_external`).
2. **Do not delete `build_release_gfortran/` this time.** Rebuild on top of it. That tree has the
   data directory cached as `.../macos_f90` from #75 and #76, and this branch must heal it: the
   configure step should print a warning that `MACOS_DATA_DIR=.../macos_f90` does not hold
   `macos_param.txt`, and `build_release_gfortran/generated/build_loc.inc` should then end in
   `/data`. Please report what you see. If it does not heal, delete the tree and say so.
3. Rebuild the mex, pymacos and SegMirMaker after the engine, as always.
4. Keep `MACOS_HOME` unset and `MACOS_resources` out of the way.

## What the branch changes

1. **The two start-up files moved**: `macos_f90/macos_param.txt` and `macos_f90/macos_glass_list.txt`
   are now `data/macos_param.txt` and `data/macos_glass_list.txt`. The engine's built-in directory
   is `data/`.
2. **`MACOS_DATA_DIR` is no longer cached with a default.** It is an optional override, empty by
   default; the default is worked out on every configure. A cached value that no longer holds
   `macos_param.txt` is dropped with a warning (item 2 above tests exactly this).
3. **`tSegMirMaker` compares numbers, not bytes** (your two failures). Same text, every number to
   1e-12, integers exactly, values below 1e-9 treated as numerical zeros; the helper is
   `mmacos/tests/private/text_numbers_match.m`. A value printed to 8 or more significant digits may
   also differ by one unit in its last printed digit. If either test still fails on arm64, the
   message names the first number that differs: send that line, with both values, and do not widen
   the tolerance yourself.
4. **The generators run on Python older than 3.12** (your finding): the backslash inside an
   f-string is gone from `mmacos/src/gen_mex_wrappers.py` and from
   `docs/mmacos-manual/tools/gen_part2.py`. Please run the generator once with Apple's `python3`
   (3.9) to confirm.
5. **Code that names the files** follows the move: `macos.design.glass`, the nested-MATLAB helper
   `macos_home_`, the Zemax converter's default catalog, `tools/gen_glass_builtin.py`.

## What to run

The same checks as #76: engine build, the probe, CLI corpus regeneration, CLI load and command
gates, mmacos fast suite, pymacos suite, SegMirMaker build, the bindings-generator diff. In the
probe, the parameter file and the glass table should now be under `<clone>/data/`, and
`MACOS_HOME=<clone>/macos_f90` should give a two-line warning for each file and then `data/`.

## Mac-specific things I could not check

- **The stale-cache healing** on your build tree (item 2 of the first section).
- **`tSegMirMaker` on arm64** with the new comparison. On Linux: 4 pass.
- **Python 3.9** on the two generator files. There is no Python older than 3.12 on the Linux box; I
  checked by scanning for the two forms 3.12 added (a backslash, or the outer quote reused, inside
  an f-string expression) across the repository's tools, and by the generator's output being
  unchanged.
- **Anything of yours that names `macos_f90/macos_param.txt`**, such as a probe driver or a
  `startup.m`.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Engine, `smacos_dvr`, GMI and mmacos mex build, on the existing build trees | all succeed; the stale cached data directory is dropped with a warning and `data/` is recorded |
| `macos` started in an empty directory, gfortran and ifx | reads both files from `data/` |
| `MACOS_HOME` pointing at `macos_f90/` | warns for each file, then uses `data/` |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the #76 records |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| mmacos fast suite | 815 pass, 0 fail (run before the `tSegMirMaker` change; that class was then rerun alone: 4 pass, one of them new) |
| `smacos_dvr`, gfortran and ifx | output identical to the #75 run apart from the parameter-file line, which now names `data/` |
| GMI regression, gfortran and ifx | 6 of 6 pass with each |
| Zemax converter tests | 34 pass, 10 skipped (they need decks that are not in the repository) |
| Generators | `gen_mex_wrappers.py` and `gen_glass_builtin.py` reproduce their committed files exactly |

Not run for this change: the mmacos model-256, 512 and 1024 groups, the manual builds and the manual-transcript journals.

## Known, and not from this branch

- The gfortran `macos` spins at 100% CPU at the prompt when its input is a pipe.
- The manuals' transcripts still show `macos_f90/macos_param.txt` in the start-up banner, and the
  three manual PDFs still carry the two-repository text; both are redone once, after #82.

## Reporting

As before: your counts beside the Linux counts, every difference named. If it is green, open the
PR from `mr_merge_data` into `mr_merge` after the #76 PR is merged, with the text below and your
Mac line added. If it is not green, report to Dave before opening it.

## PR text

```markdown
# Move the runtime data files to data/

Closes #77. Branch `mr_merge_data`, on top of `mr_merge_external` (#76).

## What changes for a user

Nothing has to be done. The two files the engine reads at start-up, `macos_param.txt` and `macos_glass_list.txt`, move from `macos_f90/` to a new top-level `data/` directory, and a build finds them there by itself.

A `MACOS_HOME` that still points at `.../macos_f90` no longer holds the files. The engine prints a two-line warning for each file and then uses `data/`. Remove `MACOS_HOME` from the shell start-up file to silence it.

An existing build tree keeps working: the next build re-records the directory.

## What is in the commit

- **The move.** `git mv` for the two files; content unchanged, tracked-file total unchanged.
- **CMake.** The default data directory is `data/`. `MACOS_DATA_DIR` is now an optional override that is empty by default: the default is worked out on every configure and is not stored in the CMake cache. Issue #75 had stored it, and a build tree made then kept pointing at `macos_f90/` after the files moved; the first test run of this branch caught that. A cached value that no longer holds `macos_param.txt` is now dropped with a warning.
- **Code that names the files.** `macos.design.glass`, the nested-MATLAB test helper, the Zemax converter's default catalog, and `tools/gen_glass_builtin.py` (which reproduces `glass_builtin.f90` exactly from the moved file).
- **Two more sibling-clone paths**, missed in #76: `mmacos/src/gen_mex_wrappers.py` (`binding_classes`) and `mmacos/src/api_binding_gap.py` read the engine source from this repository.
- **CLI test runner.** Its "built from an edited tree" check also watches `data/`.
- **`tSegMirMaker` compares numbers, not bytes.** On the Mac two of its tests failed: they compared SegMirMaker's output byte for byte with references made on x86_64, and the arm64 build prints a sixteenth significant digit that differs by one. The test now requires the same text and every number equal to 1e-12 (integers exactly; values below 1e-9 are numerical zeros), through a new helper `text_numbers_match`. A new test shows what the comparison still refuses: a change in the seventh digit, an integer, the text.
- **The generators run on Python older than 3.12.** `gen_mex_wrappers.py` and the manual tool `gen_part2.py` each had a backslash inside an f-string, which Apple's Python 3.9 rejects.
- **Documents.** `BUILD.md`, `README.md`, the mmacos README and the engine notes say `data/`; nine comments and READMEs that named `macos_f90/macos_glass_list.txt` are corrected.

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone.

| Check | Result |
|---|---|
| Engine, `smacos_dvr`, GMI and mmacos mex build, on the existing build trees | all succeed; the stale cached data directory is dropped with a warning and `data/` is recorded |
| `macos` started in an empty directory, gfortran and ifx | reads both files from `data/` |
| `MACOS_HOME` pointing at `macos_f90/` | warns for each file, then uses `data/` |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the #76 records |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| mmacos fast suite | 815 pass, 0 fail (run before the `tSegMirMaker` change; that class was then rerun alone: 4 pass, one of them new) |
| `smacos_dvr`, gfortran and ifx | output identical to the #75 run apart from the parameter-file line, which now names `data/` |
| GMI regression, gfortran and ifx | 6 of 6 pass with each |
| Zemax converter tests | 34 pass, 10 skipped (they need decks that are not in the repository) |
| Generators | `gen_mex_wrappers.py` and `gen_glass_builtin.py` reproduce their committed files exactly |

Not run for this change: the mmacos model-256, 512 and 1024 groups, the manual builds and the manual-transcript journals.

## Not covered

- **Mac and Windows.** Not run here.
- **The manuals' transcripts** still show `<macos>/macos_f90/macos_param.txt` in the start-up banner; they are regenerated with the manual rebuild after #82.
- **`macos_f90/ChfRayPosHist.txt`** is the output of a past command, committed by accident; it is not a data file and was left where it is.
```
