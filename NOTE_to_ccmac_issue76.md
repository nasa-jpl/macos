# NOTE to CCMac — Mac test of issue #76, branch `mr_merge_external` (CC, 2026-11-17)

Thank you for the #75 run: every count beside the Linux count, and your `startup.m` finding is
acted on below. Read this note from the file in your clone; it is complete as committed, with the
Linux results and the PR text in it.

## First: switch branches and delete the build directories

Your clone from the #75 run is good; there is no need to clone again.

1. `git fetch origin && git switch mr_merge_external` (it sits on top of `mr_merge_macos_home`).
2. **Delete every build directory in the clone** (`build_release_gfortran/` and any other `build_*`,
   `pymacos/src/cmake/build*`, `segmirmaker/build_*`), by full path, and the built mex and pymacos
   library. The third-party objects move inside the build tree and readline is found in a new
   place; an old build tree hides a wrong path.
3. Keep `MACOS_HOME` unset and `MACOS_resources` out of the way, as for #75.

## What the branch changes

1. **Third-party code moved to `external/`** by `git mv`, 5,571 files, content unchanged:
   `macos_f90/giza` -> `external/giza`, `macos_f90/readline-8.2` -> `external/readline-8.2`,
   `macos_f90/fits_build` -> `external/fits_build`, `macos_f90/slsqp` -> `external/slsqp`.
   `macos_f90/` now holds the 112 engine files plus `notes/` and `Archive/`.
2. **Build files follow:** `CMakeLists.txt`, the four `make*.sh` scripts, `segmirmaker/CMakeLists.txt`.
   The nested cheatsheets are now `external/giza/CLAUDE.md` and `external/slsqp/CLAUDE.md`.
3. **`tDataDir` no longer depends on `startup.m`.** Your finding: `matlab -batch` runs `startup.m`
   before the leg, so a `startup.m` that sets `MACOS_HOME` undid the test's clean environment. The
   leg now clears the variable inside MATLAB (`setenv('MACOS_HOME', '')`; an empty value is unset to
   the engine). Checked here with a `startup.m` that sets it.
4. **The documents describe one repository**: the top-level `README.md` is rewritten;
   `HOW_TO_COMPILE.md`, `BRANCHES.md`, `MAC_PORT.md`, the package READMEs and cheatsheets and both
   manuals' sources are corrected. `MAC_PORT.md` now gives `<clone>/macos/mmacos` and
   `<clone>/macos/pymacos`; please read its Part 1 commands once and say if one is wrong.
5. **Generators and tools read this repository**, not a sibling clone: `mmacos/src/gen_mex_wrappers.py`,
   `docs/macos-manual/tools/gen_cmdref.py`, the manual `Makefile` and tools, the polarization report
   tool, and the Zemax converter's default glass catalog.
6. **GMI regression references re-made** (chief-ray OPD reference; 6 of 6 on Linux with both
   compilers). Not a Mac item unless you now build GMI.
7. **`.gitignore`** lists readline's eleven generated files (a Linux build leaves them; the Mac
   uses Homebrew's readline).

## What to run

The same checks as #75, in the same order: engine build, the probe, CLI corpus regeneration, CLI
load and command gates, mmacos fast suite, pymacos suite. Then SegMirMaker, and
`python3 mmacos/src/gen_mex_wrappers.py --out-dir <scratch>` followed by a `diff` of the four
files it writes against `mmacos/src/mmacos_gen.F`, `mmacos/src/mmacos_gen_cmds.txt`,
`pymacos/src/pymacos/api_gen.py` and `pymacos/src/cmake/source/pymacos_f2py.f90` (identical here).

## Mac-specific things I could not check

- **Readline.** On Linux the scripts build the bundled `external/readline-8.2` once. On the Mac the
  engine and SegMirMaker use Homebrew's readline and the bundled tree is not built; confirm nothing
  on the Mac path still looks under `macos_f90/readline-8.2`.
- **giza under `external/`** with Apple clang: the include paths in `CMakeLists.txt` changed from
  `macos_f90/giza/src` to `external/giza/src`, for both the source and the build tree.
- **SegMirMaker.** Its Linux readline default changed from `~/dev/macos/macos_f90/readline-8.2` to a
  path relative to the repository. The Mac branch of that file was not touched.
- **`tDataDir` with a `startup.m` that sets `MACOS_HOME`.** If you still have the backup of yours,
  put it back for one run of `./run_mmacos_tests.sh tDataDir`; it should pass now. Then remove it again.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Engine and `smacos_dvr` build, gfortran and ifx | all succeed |
| Clean copy of the tracked files, gfortran | readline builds once in `external/readline-8.2`, the engine builds and runs |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the #75 records |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| mmacos fast suite | 815 pass, 0 fail |
| `smacos_dvr`, gfortran and ifx | output identical to the #75 run |
| GMI regression, gfortran and ifx | 6 of 6 pass with each (it was 2 of 6 before the references were re-made) |
| SegMirMaker build, gfortran | succeeds |
| Bindings generator | reproduces the committed files exactly |

The engine, `smacos_dvr`, GMI, the mmacos mex, pymacos and SegMirMaker were all rebuilt from scratch, with every earlier build tree deleted, before the final run of these checks.

Not run for this change: the mmacos model-256, 512 and 1024 groups, the PROPER comparisons outside the pymacos suite, the manual builds and the manual-transcript journals.

## Known, and not from this branch

- The gfortran `macos` spins at 100% CPU at the prompt when its input is a pipe.
- The three manual PDFs still carry the two-repository text; they are rebuilt once, after #82.
- Some templates, tests and Makefiles still default to `~/dev/macos/...` (for example
  `MACOS_BUILD_DIR`); a clone elsewhere sets them, as before.

## Reporting

As before: your counts beside the Linux counts, every difference named (for pymacos, say whether
the PROPER comparisons were collected: Linux collects 7,031 with them and 6,988 without). If it is
green, open the PR from `mr_merge_external` into `mr_merge` after the #75 PR is merged, with the
text below and your Mac line added. If it is not green, put your report in a note on this branch.

## PR text

```markdown
# Move the bundled third-party code to external/

Closes #76. Branch `mr_merge_external`, on top of `mr_merge_macos_home` (#75).

## What changes for a user

Nothing at run time. The engine directory `macos_f90/` now holds only MACOS's own files: 112 source and data files, plus `notes/` and `Archive/`. The 5,571 files of bundled third-party code are under `external/`:

| Was | Is | Files |
|---|---|---|
| `macos_f90/giza/` | `external/giza/` | 5,394 |
| `macos_f90/readline-8.2/` | `external/readline-8.2/` | 160 |
| `macos_f90/fits_build/` | `external/fits_build/` | 13 |
| `macos_f90/slsqp/` | `external/slsqp/` | 4 |

A build tree made before this change still works, but delete it and rebuild: the giza objects move from `build_*/macos_f90/giza/` to `build_*/external/giza/`, and readline is built once in its new place.

## What is in the commit

- **The move.** `git mv` for the four directories. Every file arrived with the same content, mode and relative path (compared by blob hash before and after), and the tracked-file total is unchanged at 12,277.
- **Build.** `CMakeLists.txt`, `makems.sh`, `makesd.sh`, `makeall.sh`, `makegfortran.sh` and `segmirmaker/CMakeLists.txt` point at `external/`. SegMirMaker's readline path was `~/dev/macos/macos_f90/readline-8.2`; it is now relative to the repository.
- **CLI test runner.** Its "built from an edited tree" check now also watches `external/` and `cmake/`.
- **Pointers and documents.** The root and nested `CLAUDE.md` files, `README.md`, `HOW_TO_COMPILE.md`, `DEV_FILES.md`, the giza build notes, `platform_requirements.txt`, the SegMirMaker README, one manual appendix, and the two legacy builders in `macos_f90/`.
- **GMI regression references re-made.** Four of the six GMI tests failed, on `dev` as well. The references dated from before 2026-10-10, when the engine's default OPD reference changed from the aperture mean to the chief ray, so every map had moved by one constant (e5hex1: -8.2589e-06 mm, the engine's printed "Average OPD"). With the mean forced back the e2e_pie maps matched the old files to 3e-14; the e5hex1 maps matched to 8.7e-11 mm, six units in the last place of the 66,573 mm path, which is round-off from engine changes since. The new references are chief-referenced; the gfortran and ifx mex files reproduce them bit for bit. `GMI/regression/README.md` records this.
- **`tDataDir` made independent of the user's `startup.m`.** On the Mac the test failed where MATLAB's `startup.m` sets `MACOS_HOME`: `matlab -batch` runs it before the test leg. The leg now clears the variable inside MATLAB.
- **Readline's build output is ignored.** The bundled readline is configured and built in its own source directory, which left eleven generated files showing as untracked after every Linux build. They are now listed in `.gitignore` by exact path; no tracked file is hidden by the new lines (the count of tracked files that match an ignore rule is unchanged at 90).
- **The documents describe one repository.** About 120 path references in 48 Markdown files were rewritten by fixed rules (`MACOS_resources/mmacos/...` becomes `mmacos/...`, and so on), and the passages that describe the set-up in prose were rewritten by hand: `HOW_TO_COMPILE.md` (clone, branch and pull steps), `BRANCHES.md` (a new first section; the older sections are kept as the record), the root `CLAUDE.md`, `MAC_PORT.md`, the package READMEs and cheatsheets, and the installation and introduction text of both manuals' sources. Dated records (briefs, notes, plans, logs, status files) were left as written.
- **Generators and tools find this repository.** The bindings generator `mmacos/src/gen_mex_wrappers.py` looked for the engine in a sibling clone named `macos`; it now reads this repository's `macos_f90/`, and it reproduces the four committed binding files and `25_defaults.md` exactly. The same correction is made in `docs/macos-manual/tools/gen_cmdref.py`, the manual `Makefile`, four manual tools, the polarization report tool and the Zemax converter's default glass catalog.
- **Top-level README rewritten for one repository.** A section on what the repository holds, a directory table that matches the tree, one clone and nothing to set, every `MACOS_resources/...` path corrected, and the `macos_param.txt` passages corrected (they sent the reader to `mmacos/`, which does not hold the file).

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone.

| Check | Result |
|---|---|
| Engine and `smacos_dvr` build, gfortran and ifx | all succeed |
| Clean copy of the tracked files, gfortran | readline builds once in `external/readline-8.2`, the engine builds and runs |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the #75 records |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| mmacos fast suite | 815 pass, 0 fail |
| `smacos_dvr`, gfortran and ifx | output identical to the #75 run |
| GMI regression, gfortran and ifx | 6 of 6 pass with each (it was 2 of 6 before the references were re-made) |
| SegMirMaker build, gfortran | succeeds |
| Bindings generator | reproduces the committed files exactly |

The engine, `smacos_dvr`, GMI, the mmacos mex, pymacos and SegMirMaker were all rebuilt from scratch, with every earlier build tree deleted, before the final run of these checks.

Not run for this change: the mmacos model-256, 512 and 1024 groups, the PROPER comparisons outside the pymacos suite, the manual builds and the manual-transcript journals.

## Not covered

- **Mac and Windows.** Not run here. The Windows commands in the README were corrected by reading, not by running.
- **Other builders.** `cmake/CMakeLists.txt` (the Windows draft) names none of the moved directories and was not touched. The legacy `macos_f90/CMakeLists.txt` and `Make-Linux-x86_64.T` have the new paths but were not run.
- **The manual PDFs are not rebuilt.** `macosMan.pdf`, `macosCmdRef.pdf` and `mmacosMan.pdf` still carry the two-repository text; their sources are corrected. One rebuild after the last reorganization step (#82) covers every path change.
- **The generated command reference is behind its sources**, for reasons that predate this branch: regenerating adds one journal to `10_engine_commands.md` and the `design.Lens` entries to `30_higher_level.md`. Those two files are left as committed.
- **Other absolute paths.** Some templates, tests and Makefiles still default to `~/dev/macos/...` (for example `MACOS_BUILD_DIR`); a clone elsewhere sets them, as before.
```
