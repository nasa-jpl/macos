# NOTE to CCMac — Mac test of issue #75, branch `mr_merge_macos_home` (CC, 2026-11-17)

## First: delete your macos clone and clone it again

Do this before anything else.

1. **Check the old clone for work that exists nowhere else**: `git status`, `git stash list`, and
   `git log --branches --not --remotes --oneline`. Copy anything you need to keep to a directory
   outside the clone. Do not push from the old clone.
2. **Delete the old macos clone** (the whole directory, by its full path).
3. **Clone again**: `git clone git@github.com:nasa-jpl/macos.git` (not `--recursive`, see below),
   then `git switch mr_merge_macos_home`.
4. **Move your MACOS_resources clone out of the way** (rename it or delete it). Its content is now
   inside macos, and a test that still reads the old directory must fail, not pass.

Why: the history of `nasa-jpl/macos` was rewritten on 2026-11-15/16 and `mr_merge` was force-updated
on 2026-11-17, so an older clone holds commits that no longer exist on the remote. A pull into it
gives a tangle, and a push from it would put the old history back.

## Where we are

The group met on 2026-11-17 and agreed the plan. `mr_merge` (at `0a526602`) is now only the merge
of MACOS_resources into macos: **one repository**, with `mmacos/`, `pymacos/`, `GMI/`,
`segmirmaker/`, `rx_converter/` and `optical_design/` at the top level beside `macos_f90/`. The
reorganization follows as GitHub issues #75 to #82, one branch and one PR per issue, in order,
each branch starting from the one before. Luis takes #83 (CMake on every platform) after #82.

This note is for the first one: **#75, "Fix MACOS_HOME error paths"**, on branch
`mr_merge_macos_home`. CC has run it on Linux. Dave asks you to run it on the Mac and, when it is
green, open the pull request from `mr_merge_macos_home` into `mr_merge`. The PR text is at the end
of this note.

## Set-up that differs from every earlier pair

- **One clone.** There is no MACOS_resources clone any more, and nothing in this branch may read
  one. The fresh clone from the first section is the only tree on the machine.
- **Do not clone with `--recursive`.** The merge left a submodule entry (`pymacos/src/macos`)
  whose commit no longer exists; a recursive clone fails. It is a known leftover, not part of #75.
- **`MAC_PORT.md` still gives the two-repository paths.** Read `<clone>/MACOS_resources/mmacos` as
  `<clone>/mmacos`, and the same for `pymacos`.
- **Unset `MACOS_HOME`** in the shell you test from (`unset MACOS_HOME`), and do not pass it on any
  command line. That is the point of the test.

## What the branch changes

1. **The engine finds its own data files.** CMake writes the data directory (`macos_f90/`) into
   `build_*/generated/build_loc.inc`, and `param_mod.F` includes it. The search order is unchanged:
   working directory, executable's directory, `$MACOS_HOME`, then this built-in directory. The
   fourth was dead until now (`'BUILD_LOC'`, a macro name inside quotes). `BUILD.md` has the
   section "Where the engine finds its data files".
2. **No parameter file anywhere is exit status 1**, not 0 (`macos_io_failure`).
3. **Nothing in the tree sets or requires `MACOS_HOME`.** `run_mmacos_tests.sh` no longer exports
   it; 56 template asserts and 20 batch-script defaults are gone; six pymacos tests and four mmacos
   tests find the checkout from their own location (`pymacos/tests/repo_paths.py`,
   `mmacos/tests/private/repo_root_.m`).
4. **Two new tests:** `pymacos/tests/test_data_dir.py` and `mmacos/tests/tDataDir.m` (in the fast
   suite). Each starts the engine in an empty directory with `MACOS_HOME` unset, then with a stale
   value.
5. **Seven one-line path fixes the merge needs** (CLI corpus generator and command gate, `makegmi.sh`,
   `makeall.sh`, the two CMake defaults for `GMI/` and `mmacos/`).

## What to run

1. **Engine build**, gfortran: `source ./makems.sh release gfortran`. Then
   `cat build_release_gfortran/generated/build_loc.inc` — it must hold your clone's `macos_f90`
   path, split across lines if it is longer than 56 characters.
2. **The probe.** From an empty scratch directory, `MACOS_HOME` unset, start
   `build_release_gfortran/bin/macos` in a terminal, give a model size, and check for
   `Using parameter file "<clone>/macos_f90/macos_param.txt"` and the glass-table line. Repeat with
   `MACOS_HOME=/no/such/dir`: a two-line warning, then the same parameter file.
3. **CLI gates:** regenerate the corpus lists first (`cli_tests/make_corpus.sh`; the committed
   lists still carry Linux absolute paths), then `MACOS_CLI_TREE=. ./run_cli_tests.sh load` and
   `cmd`.
4. **mmacos:** build the mex, then `./run_mmacos_tests.sh fast` (set `MACOS_BUILD_DIR` if the
   clone is not `~/dev/macos`; the runner still defaults that one).
5. **pymacos:** build, then the suite from `pymacos/tests`.

## Mac-specific things I could not check

- **`tDataDir` starts its nested MATLAB with `env -u MACOS_HOME`** inside a `system()` call. BSD
  `env` has `-u`, and a `system()` from MATLAB on macOS drops the parent environment in any case,
  but this is the first test to rely on it. It reads only MATLAB's own `INIT_OK` lines, not the
  Fortran unit-6 output.
- **The four nested-MATLAB tests** (`tRxRefuse`, `tRxShortCoef`, `tRxLineEndings`,
  `tSaveNewFileQuiet`) still pass `MACOS_HOME` to the child through `macos_home_()`. Its fallback
  was `<parent>/macos/macos_f90`, which is right only for a clone named `macos`; it is now this
  checkout's `macos_f90`. Please run one of them from a clone with another name.
- **`INCLUDE 'build_loc.inc'`** under your gfortran and with a `/Users/...` path.
- **The Finder-launched MATLAB case** that `run_pupil_id.m` used to guard (no shell profile, so no
  `MACOS_HOME`). Its guard is removed; `run_pupil_id` should now simply work there.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Engine build, gfortran and ifx | both succeed |
| `macos` started in an empty directory | finds `macos_param.txt` and the glass table in its own tree; before the change it stopped with status 0 |
| `MACOS_HOME` naming a directory that does not exist | warns, then uses its own tree |
| A build configured with a nonexistent data directory | stops with status 1 |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the records taken before the change |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,028 pass, 0 fail, 1 skipped; the new `test_data_dir.py` (2 tests) was added after that run started and passes on its own |
| mmacos fast suite | 815 pass, 0 fail (the earlier 814 plus `tDataDir`) |
| `smacos_dvr`, gfortran and ifx | output identical to the run before the change |
| GMI, gfortran | builds; the regression gives 2 pass, 4 fail, the same four tests and the same numbers as before the change and as `dev` (stale references) |

Not run for this change: the mmacos model-256, 512 and 1024 groups, the PROPER comparisons, the manual builds and the manual-transcript journals.

## Known, and not from this branch

- The gfortran `macos` spins at 100% CPU at the prompt when its input is a pipe and not a
  terminal. `dev` does the same.
- The GMI regression fails four of six tests against stale references, on `dev` as well.
- The top-level `README.md` still describes two repositories.

## Reporting

As before: your counts beside the Linux counts, and every difference named. If it is green, open
the PR with the text below. If it is not, do not open it; put your report in a note on this branch
and Dave will pull it.

## PR text

Title and body, to paste as they are. Add your Mac line to the test results.

```markdown
# Make the engine find its own data files, so MACOS_HOME is optional

Closes #75. Branch `mr_merge_macos_home`, on top of `mr_merge` at `0a526602`.

## What changes for a user

Nothing has to be set. A build now finds `macos_param.txt` and `macos_glass_list.txt` in the tree it was built from. `MACOS_HOME` still works and still wins over the built-in directory, so an existing shell start-up file needs no edit. A `MACOS_HOME` that points at a directory without the file prints a warning and the engine uses its own directory instead. That is what makes a later rename of `macos_f90/` cost users nothing.

## Why it did not work before

The engine has always had a fourth place to look, after the working directory, the executable's directory and `$MACOS_HOME`: the build directory. It never worked. `param_mod.F` set it with `build_loc = 'BUILD_LOC'`, and a macro name inside quotes is not replaced by either preprocessor, so every library carried the literal text `BUILD_LOC`. The top-level CMake build did not define the macro in any case.

## What is in the commit

One commit, in two parts.

**Part 1: point the build scripts and the CLI tests at one repository.** Seven one-line path fixes that the merge itself needs; nothing here is about `MACOS_HOME`. Without them the CLI test list cannot be generated, the command journals stop after 46 of 122, and `makegmi.sh` and `makegfortran.sh` fail.

- `cli_tests/make_corpus.sh`, `cli_tests/cli_cmd_gate.py`: the resources tree is this repository.
- `cli_tests/make_corpus.sh`: a missing `~/dev/MACOS_sandbox` no longer aborts the script (an old bug the merge exposed).
- `makegmi.sh`, `makeall.sh`, `CMakeLists.txt`, `macos_f90/CMakeLists.txt`: `GMI/` and `mmacos/` are in this repository.

**Part 2: make the engine find its own data files.**

- **Engine.** CMake writes the data directory into `build_*/generated/build_loc.inc` and `param_mod.F` includes it (`cmake/macos_build_loc.cmake`). The path is written in short pieces, so a long clone path cannot run past the fixed-form line limit. `-DMACOS_DATA_DIR=<dir>` overrides it at configure time.
- **Exit status.** With no parameter file anywhere, the engine now stops with status 1. It was status 0, which let a pymacos run that died part-way read as passed.
- **Tests.** Six pymacos tests and four mmacos tests used `MACOS_HOME` to find the engine *source* and the manual's example decks, with `~/dev/macos/macos_f90` as the fallback. They now find the checkout from their own location (`pymacos/tests/repo_paths.py`, `mmacos/tests/private/repo_root_.m`).
- **Test runner.** `run_mmacos_tests.sh` no longer sets `MACOS_HOME`. If the caller has it set, the runner says so in one line.
- **Templates and batch scripts.** 56 template scripts asserted that `MACOS_HOME` was set, and 20 batch scripts defaulted it to `~/dev/macos/macos_f90`. Both are removed. `macos.design.glass` and the nested-MATLAB tests fall back to this checkout, not to a sibling clone named `macos`.
- **New tests.** `pymacos/tests/test_data_dir.py` and `mmacos/tests/tDataDir.m` start the engine in an empty directory with `MACOS_HOME` unset, and again with a stale value. Both fail on a build made before this change.
- **Documents.** `BUILD.md` has a new section, "Where the engine finds its data files". The package READMEs, `MAC_PORT.md` and the template READMEs no longer tell the reader to set `MACOS_HOME`.

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone, so nothing could read an old `macos` or `MACOS_resources` tree.

| Check | Result |
|---|---|
| Engine build, gfortran and ifx | both succeed |
| `macos` started in an empty directory | finds `macos_param.txt` and the glass table in its own tree; before the change it stopped with status 0 |
| `MACOS_HOME` naming a directory that does not exist | warns, then uses its own tree |
| A build configured with a nonexistent data directory | stops with status 1 |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the records taken before the change |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,028 pass, 0 fail, 1 skipped; the new `test_data_dir.py` (2 tests) was added after that run started and passes on its own |
| mmacos fast suite | 815 pass, 0 fail (the earlier 814 plus `tDataDir`) |
| `smacos_dvr`, gfortran and ifx | output identical to the run before the change |
| GMI, gfortran | builds; the regression gives 2 pass, 4 fail, the same four tests and the same numbers as before the change and as `dev` (stale references) |

Not run for this change: the mmacos model-256, 512 and 1024 groups, the PROPER comparisons, the manual builds and the manual-transcript journals.

## Not covered

- **Mac.** Not run here; CCMac runs it from this branch (`NOTE_to_ccmac_issue75.md`).
- **Windows.** Not run. The Windows build skips the executable-directory and `MACOS_HOME` lookups, so it should now find the files through the built-in directory; the top-level README still tells Windows users to `cd` to a directory that holds `macos_param.txt`, which stays correct.
- **Other builders.** `cmake/CMakeLists.txt`, `macos_f90/CMakeLists.txt` and `macos_f90/Makefile` do not write the include. They compile as before and have no built-in directory. To give one the directory, call `macos_write_build_loc`, define `MACOS_BUILD_LOC_INC`, and add the include path.
- **A binary copied away from its build tree** still needs `MACOS_HOME`, and inside MATLAB a missing parameter file still ends MATLAB. A clean error there needs an engine change of its own.
```
