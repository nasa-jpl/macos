# NOTE to CCMac — Mac test of issue #79, branch `mr_merge_app` (CC, 2026-11-17)

This is the move with the widest reach so far: the three interfaces go into `app/`. Read this note
from `devnotes/` in your clone.

## First: switch, then rebuild the interfaces in their new place

1. `git fetch origin && git switch mr_merge_app` (it sits on top of `mr_merge_devnotes`).
2. **Keep `build_release_gfortran/`.** The engine did not move. If your tree was configured with
   GMI or mmacos built through CMake, the configure step should warn that the cached
   `GMI_SRC_DIR` or `MMACOS_SRC_DIR` no longer holds the sources and switch to `app/...`; say if
   you see that.
3. **Delete the old interface build products by full path**: any `pymacos/src/cmake/build*`
   directory that came along under `app/pymacos/src/cmake/` (its CMake cache names the old source
   directory and the configure fails on it), the built pymacos library, and the mmacos mex. Then
   rebuild the mex and pymacos from `app/mmacos` and `app/pymacos/src/cmake`.
4. **MATLAB and Python paths** are now `<clone>/app/mmacos/src`, `<clone>/app/GMI` and
   `<clone>/app/pymacos/src`. Check your `startup.m`, your venv's path file and any probe driver.
5. Keep `MACOS_HOME` unset and `MACOS_resources` out of the way.

## What the branch changes

1. **`GMI/`, `mmacos/`, `pymacos/` are `app/GMI/`, `app/mmacos/`, `app/pymacos/`** (`git mv`, 4,282
   files, content unchanged). `segmirmaker/`, `rx_converter/` and `optical_design/` stay at the top.
2. **The leftover submodule is gone**: the `pymacos/src/macos` entry, `.gitmodules`, and the check
   for it in pymacos's CMake file. Please try `git clone --recursive` once into a scratch
   directory; it should now succeed.
3. **`.gitignore`**: eleven rules re-anchored under `app/`.
4. **Code that climbs to the repository root from inside a package** climbs one level more:
   `repo_root_.m`, `repo_paths.py`, twelve mmacos test classes, `macos.design.glass`,
   `macos.design.segmirmaker_run`, the bindings generator, a few tools and templates.
5. **CLI tests**: `$RES` expands to `<clone>/app`; no journal changed. Regenerate the corpus lists
   as always; the core list should again have 166 decks, 70 of them under `app/`.
6. **`MAC_PORT.md`** gives `<clone>/macos/app/mmacos` and `<clone>/macos/app/pymacos`.

## What to run

The full gate, as for #76: engine build, the probe, CLI corpus regeneration, CLI load and command
gates, the mmacos fast suite, the pymacos suite, SegMirMaker build, the bindings-generator diff
(now `python3 app/mmacos/src/gen_mex_wrappers.py --out-dir <scratch>`; the four files are
`app/mmacos/src/mmacos_gen.F`, `app/mmacos/src/mmacos_gen_cmds.txt`,
`app/pymacos/src/pymacos/api_gen.py`, `app/pymacos/src/cmake/source/pymacos_f2py.f90`).

## Mac-specific things I could not check

- **The classes that were filtered on the Mac before SegMirMaker was built there** (`tRunCompare`,
  `tSegmentRx`, `tRunMet`, `tEdgeSensors`, `tMet`, `tRunSegmentation`, `tRunSensitivities`): they
  find `segmirmaker/test_in` and the SegMirMaker binary through the changed path code.
- **`segmirmaker_bin.m` and `macos.design.segmirmaker_run`**: they look for the binary in
  `<clone>/segmirmaker/build_release_gfortran/` first on the Mac.
- **The mmacos `Makefile` and `run_mmacos_tests.sh` defaults** still say
  `$HOME/dev/macos/build_release_gfortran`; unchanged, but your clone path decides whether you must
  set `MACOS_BUILD_DIR`.
- **Recursive clone**, item 2 above.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Engine and `smacos_dvr` build, gfortran and ifx; GMI and mmacos mex build | all succeed |
| pymacos build and suite, from `app/pymacos` | 7,030 pass, 0 fail, 1 skipped |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; matched by name across the new paths, every deck gives the same result as before the move |
| CLI command journals, gfortran and ifx | 122 of 122 pass, with no journal edited |
| mmacos fast suite, from `app/mmacos` | 816 pass, 0 fail |
| `tCalib` and `tViewRx` (outside the fast suite; their paths changed) | 9 pass and 7 pass |
| GMI regression, gfortran, from `app/GMI` | 6 of 6 pass |
| `smacos_dvr`, gfortran and ifx | output identical to the #77 run |
| Zemax converter tests | 34 pass, 10 skipped (they need decks that are not in the repository) |
| Generators | `gen_mex_wrappers.py` reproduces the four binding files exactly from `app/mmacos/src`; `gen_cmdref.py` runs |

## Known, and not from this branch

- About 60 templates and a few tests still default to `~/dev/MACOS_resources/...`; broken by the
  merge, to be fixed as a change of its own.
- 91 Markdown links inside the templates and challenges already pointed at files that are not in
  the repository.
- The gfortran `macos` spins at 100% CPU at the prompt when its input is a pipe.

## Reporting

As before. Open the PR from `mr_merge_app` with base `mr_merge_devnotes` (see
`devnotes/NOTE_to_ccmac_pr_bases.md`), with the text below and your Mac line added.

## PR text

```markdown
# Move GMI, mmacos and pymacos into app/

Closes #79. Branch `mr_merge_app`, on top of `mr_merge_devnotes` (#78).

## What changes for a user

The three interfaces move one level down:

| Was | Is | Files |
|---|---|---|
| `GMI/` | `app/GMI/` | 36 |
| `mmacos/` | `app/mmacos/` | 4,034 |
| `pymacos/` | `app/pymacos/` | 212 |

After pulling this change:

- **MATLAB path:** `<macos>/app/mmacos/src` and `<macos>/app/GMI` (it was `<macos>/mmacos/src` and `<macos>/GMI`).
- **Python path:** `<macos>/app/pymacos/src`.
- **Rebuild the interfaces.** The engine build tree can stay; the mex files, the pymacos library and its CMake build directory must be made again in the new place. Delete `pymacos/src/cmake/build` if it came along: its cache names the old source directory.
- **Personal scripts** that name `mmacos/...`, `pymacos/...` or `GMI/...` under the repository need `app/` in front.

`segmirmaker/`, `rx_converter/` and `optical_design/` stay at the top; the issue names only the three interfaces.

## What is in the commit

- **The move.** `git mv` for the three directories. Every file arrived with the same content, mode and relative path (blob hashes compared before and after).
- **`.gitignore`.** Eleven rules that named `GMI/`, `mmacos/` or `pymacos/` paths are re-anchored under `app/`. The set of tracked files that match an ignore rule is the same 90 files before and after, compared name by name. Without this, eleven rules would have stopped applying in silence.
- **The leftover submodule is removed.** The merge had left a submodule entry, `pymacos/src/macos`, pointing at a commit that no longer exists, with `.gitmodules`; a recursive clone failed on it. The entry, `.gitmodules` and the check for that directory in pymacos's CMake file are gone. A recursive clone now works.
- **Build.** `CMakeLists.txt`, `makegmi.sh`, `makeall.sh` and the legacy `macos_f90/CMakeLists.txt` find GMI and mmacos under `app/`. A build tree made before the move has the old GMI or mmacos directory in its CMake cache; that value is replaced, with a warning.
- **Code that climbs to the repository root from inside a package** now climbs one level more: the test helpers `repo_root_.m` and `repo_paths.py`, twelve mmacos test classes that find `segmirmaker/test_in`, `macos.design.glass`, `macos.design.segmirmaker_run`, the bindings generator, and four tools and templates. Five of these also still pointed at a sibling clone named `macos`, which the merge had already broken; they now read this repository (`tCalib`, `tViewRx`, `tDrawFrameGuard`, `view_rx_demo.m`, `coro_planet_demo.m`).
- **Tools outside `app/`** that read the packages: the command-reference generator, the manual `Makefile` and tools, the Zemax and CODE V converter tests.
- **CLI tests.** `$RES` in the command journals now expands to `app/`; the journals themselves are unchanged. The deck list was regenerated and names the same 166 decks.
- **`release-exclude.txt`.** Fourteen paths gain `app/`; the release strip removes the same 572 files.
- **Documents.** Paths to the three packages are rewritten in 68 Markdown files, each rewritten path checked to exist under `app/`; the layout diagrams in `README.md` and `HOW_TO_COMPILE.md` are redrawn.

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone.

| Check | Result |
|---|---|
| Engine and `smacos_dvr` build, gfortran and ifx; GMI and mmacos mex build | all succeed |
| pymacos build and suite, from `app/pymacos` | 7,030 pass, 0 fail, 1 skipped |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; matched by name across the new paths, every deck gives the same result as before the move |
| CLI command journals, gfortran and ifx | 122 of 122 pass, with no journal edited |
| mmacos fast suite, from `app/mmacos` | 816 pass, 0 fail |
| `tCalib` and `tViewRx` (outside the fast suite; their paths changed) | 9 pass and 7 pass |
| GMI regression, gfortran, from `app/GMI` | 6 of 6 pass |
| `smacos_dvr`, gfortran and ifx | output identical to the #77 run |
| Zemax converter tests | 34 pass, 10 skipped (they need decks that are not in the repository) |
| Generators | `gen_mex_wrappers.py` reproduces the four binding files exactly from `app/mmacos/src`; `gen_cmdref.py` runs |

Not run for this change: the ifx GMI mex, the other mmacos groups (model 256, 512, 1024, PROPER), SegMirMaker rebuild (its files did not move), the manual builds.

## Not covered

- **Mac and Windows.** Not run here. The Windows commands in the README are rewritten by rule, not run.
- **About 60 templates and a few tests still name `~/dev/MACOS_resources/...`** in a default path. They were broken by the merge, before this branch, and are a change of their own.
- **91 Markdown links inside the templates and challenges** already pointed at files that are not in the repository; this branch did not add to them.
- **The manual PDFs** are rebuilt once, after #82.
```
