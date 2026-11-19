# NOTE to CCMac — Mac test of issue #81, branch `mr_merge_src` (CC, 2026-11-18)

The engine directory is renamed: `macos_f90/` is `src/`. Read this note from `devnotes/` in your
clone.

## First: switch, then rebuild everything in place

1. `git fetch origin && git switch mr_merge_src` (it sits on top of `mr_merge_tests`).
2. **Keep `build_release_gfortran/` and rebuild in it**: `source ./makems.sh release gfortran`.
   On Linux a tree made before the rename reconfigures and recompiles from `src/` with no
   complaint. Say if yours does not; deleting the tree and building again is the fallback.
3. **Rebuild what links the engine library**: the mmacos mex, pymacos, GMI, SegMirMaker.
4. **Delete a leftover `<clone>/macos_f90/` directory by full path**, if the switch leaves one
   (untracked files of yours would be all it holds).
5. Keep `MACOS_HOME` unset. If any start-up file of yours still sets it to `.../macos_f90`, the
   engine now warns twice at start-up and carries on; remove the setting.

## What the branch changes

1. **`git mv macos_f90 src`** (127 files).
2. **Six comment lines in the engine source, no code**: three pointers to `src/CLAUDE.md`, three
   to `tests/cli/cmd/...` (held over from #80). This is why the engine is rebuilt.
3. **Build**: `CMakeLists.txt` (`smacos_src_dir`), `makegmi.sh` and `makeall.sh`
   (`macossrc_dir`), and the Windows draft `cmake/CMakeLists.txt`.
4. **Code that reads the engine source**: `gen_mex_wrappers.py`, `api_binding_gap.py`,
   `repo_paths.py` (`ENGINE_DIR`), `tApiCoverage.m`, the manual generators, the two tools in
   `tools/`, `tests/cli/tools/hzdef_check.py`.
5. **`tests/cli/run_cli_tests.sh`**: the "+dirty" check watches `src`.
6. **Documents**: both `CLAUDE.md` files, `src/notes/`, README, build guides, manual sources.

## What to run

The full gate: engine build in the existing tree, the empty-directory probe, CLI load and command
gates, the mmacos fast suite with `tCalib` and `tViewRx`, the pymacos build and suite, the GMI
regression, SegMirMaker build, and the bindings-generator diff
(`python3 app/mmacos/src/gen_mex_wrappers.py --targets mex,f2py,veneer --out-dir <scratch>`).

## Mac-specific things I could not check

- **An existing Mac build tree after the rename** (item 2 above).
- **The GMI build on the Mac**: `makegmi.sh` passes `macossrc_dir=<clone>/src` to GMI's Makefile.
- **Editor and debugger settings** in your clone that name `macos_f90`.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Engine and `smacos_dvr` rebuild in the existing trees, gfortran and ifx; GMI and mmacos mex | all succeed |
| `macos` started in an empty directory | finds the parameter file and the glass table in `data/` |
| `MACOS_HOME` naming the old `macos_f90` directory | warns, then uses `data/` |
| CLI load gate, gfortran and ifx | 166 of 166 decks load with each; all 12 result columns of every deck equal the record taken before the rename |
| CLI command journals, gfortran and ifx | 122 of 122 pass with each; the 3 must-fail journals still fail |
| mmacos fast suite | 816 pass, 0 fail, 0 incomplete |
| `tCalib` and `tViewRx` (outside the fast suite) | 9 pass and 7 pass |
| pymacos build and suite | 7,030 pass, 0 fail, 1 skipped |
| GMI regression, gfortran | 6 of 6 pass |
| `smacos_dvr`, gfortran and ifx | output identical to the #79 run |
| Generators | `gen_mex_wrappers.py` reproduces the four binding files and `25_defaults.md` exactly from `src/`; `gen_cmdref.py` reproduces `20_bindings.md` and `40_suite_coverage.md` |
| Clean checkout, gfortran | a fresh clone builds the engine and SegMirMaker from `src/`, finds its data files, and leaves no untracked file behind |

## Known, and not from this branch

- About 60 templates and a few tests still default to `~/dev/MACOS_resources/...`.
- Two generated command-reference files are behind their sources (one journal, the
  `design.Lens` entries); they wait for the manual rebuild after #82.
- The gfortran `macos` spins at 100% CPU at the prompt when its input is a pipe.

## Reporting

As before. Open the PR from `mr_merge_src` with base `mr_merge_tests` (see
`devnotes/NOTE_to_ccmac_pr_bases.md`), with the text below and your Mac line added.

## PR text

```markdown
# Rename macos_f90/ to src/

Closes #81. Branch `mr_merge_src`, on top of `mr_merge_tests` (#80).

## What changes for a user

The engine source directory is renamed:

| Was | Is | Tracked files |
|---|---|---|
| `macos_f90/` | `src/` | 127 |

After pulling this change:

- **Rebuild.** An existing build tree works: run the same build script again and it recompiles from `src/`. Then rebuild the interfaces you use (the mmacos mex, pymacos, GMI, SegMirMaker), because they link the engine library.
- **Nothing has to be set.** The engine finds `macos_param.txt` and the glass table in `data/`, as it has since #75 and #77. A `MACOS_HOME` that still points at the old `macos_f90` directory prints a warning and the engine uses its own data directory.
- **Personal scripts and editor settings** that name `macos_f90/...` need `src/` in its place.

This is the rename that #75 was done for: no user setting depends on the directory's name.

## What is in the commit

- **The rename.** `git mv macos_f90 src`. Every file arrived with the same mode and relative path, and the same content apart from the six comment lines below.
- **Build.** `CMakeLists.txt` compiles the engine from `src/`; `makegmi.sh` and `makeall.sh` pass `src/` to the GMI build; the Windows draft `cmake/CMakeLists.txt` names `../src`.
- **Six comment lines in the engine source, no code.** Three pointers to the engine's `CLAUDE.md` say `src/CLAUDE.md`. Three pointers to command journals say `tests/cli/cmd/...`; they were held over from #80 so that #80 needed no engine rebuild.
- **Code that reads the engine source** finds it in `src/`: the bindings generator and its gap check, the pymacos and mmacos API-coverage tests (`repo_paths.py`, `tApiCoverage.m`), the command-reference and appendix generators, the glass-table generator, the error-message tool, and the Hex-Zernike check in `tests/cli/tools/`.
- **CLI test runner.** Its "built from an edited tree" check watches `src/`.
- **Manual tools.** The Zernike-table target no longer sets `MACOS_HOME` to the engine directory, which has held no data file since #77. The transcript cleaner hides the clone path of the data files under `data/`, where the engine now prints them.
- **`release-exclude.txt`.** One path re-anchored; the release strip removes the same 575 files, compared name by name.
- **Documents.** `CLAUDE.md` at the top and in `src/`, the nine engine notes in `src/notes/`, the README, the build guide, `BUILD.md`, `MAC_PORT.md`, the package READMEs and cheatsheets and the manual sources name `src/`. Dated records were left as written.

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone. Every binary was rebuilt for this change, in build trees made before the rename.

| Check | Result |
|---|---|
| Engine and `smacos_dvr` rebuild in the existing trees, gfortran and ifx; GMI and mmacos mex | all succeed |
| `macos` started in an empty directory | finds the parameter file and the glass table in `data/` |
| `MACOS_HOME` naming the old `macos_f90` directory | warns, then uses `data/` |
| CLI load gate, gfortran and ifx | 166 of 166 decks load with each; all 12 result columns of every deck equal the record taken before the rename |
| CLI command journals, gfortran and ifx | 122 of 122 pass with each; the 3 must-fail journals still fail |
| mmacos fast suite | 816 pass, 0 fail, 0 incomplete |
| `tCalib` and `tViewRx` (outside the fast suite) | 9 pass and 7 pass |
| pymacos build and suite | 7,030 pass, 0 fail, 1 skipped |
| GMI regression, gfortran | 6 of 6 pass |
| `smacos_dvr`, gfortran and ifx | output identical to the #79 run |
| Generators | `gen_mex_wrappers.py` reproduces the four binding files and `25_defaults.md` exactly from `src/`; `gen_cmdref.py` reproduces `20_bindings.md` and `40_suite_coverage.md` |
| Clean checkout, gfortran | a fresh clone builds the engine and SegMirMaker from `src/`, finds its data files, and leaves no untracked file behind |

Not run for this change: the ifx GMI mex, the other mmacos groups (model 256, 512, 1024, PROPER), the manual builds.

## Not covered

- **Mac and Windows.** Not run here. The Windows draft `cmake/CMakeLists.txt` and the Windows commands in the README were corrected by reading.
- **The legacy builders inside `src/`** (`CMakeLists.txt`, `Makefile`, `Make-Linux-x86_64.T`) reach their neighbours by relative path and needed no edit; they were not run.
- **The manual PDFs** are rebuilt once, after #82. Two generated command-reference files are behind their sources for reasons that predate this series (one journal, the `design.Lens` entries); they were corrected by rule, not regenerated.
```
