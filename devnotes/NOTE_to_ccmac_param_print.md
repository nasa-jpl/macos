# NOTE to CCMac — the print fix, branch `mr_merge_param_print` (CC, 2026-11-18)

This branch sits on top of `mr_merge_manuals`. It changes engine source, so it needs a Mac gate,
but a short one. Read this note from `devnotes/` in your clone.

## What it is

With the ifx build in a clone whose path is long, the engine split two messages across lines:
`Using parameter file "..."` at start-up and `Input file ... being loaded.` when a deck is loaded.
Both were list-directed writes to the terminal, which ifx wraps at 80 columns. Neither wraps now.
For a path that fitted before, the text is unchanged, character for character.

## What changed

1. **`src/macosio.F`**: the CLI's deck-load line goes through the message routine (`emsgRec` and
   `emsg_put`), like its twin in `smacosio.F`. Made with `tools/emsg_convert.py`.
2. **`src/param_mod.F`**: the parameter-file line, the stale-`MACOS_HOME` warning and the
   "build's data directory" line use explicit formats.
3. **`tests/cli/cmd/STARTUP_size.jou`**: finds its binary with `ps -ww`. Plain `ps` cut a long
   path at 80 columns, and with gfortran in a long path the binary was then not found.
4. **`src/notes/api_bindings.md`**: the record of 1 and 2.

## What to do on the Mac

1. `git fetch origin && git switch mr_merge_param_print`.
2. Rebuild the engine in the existing tree: `source ./makems.sh release gfortran`. Then rebuild
   the mex and pymacos as usual: `param_mod.F` is part of the engine library.
3. `tests/cli/run_cli_tests.sh load` and `cmd`. Expected: 166 of 166 with no deck changed, and
   122 of 122.
4. The mmacos fast suite and the pymacos suite. `tDataDir` and `test_data_dir.py` are the two
   that read the start-up messages.

## What I expect on the Mac

No visible change. gfortran never wrapped these lines, so the Mac build should print what it
printed before. Two things are worth a look:

- **`STARTUP_size.jou`**: `ps -ww` is the same option on macOS, but it is the only line here that
  depends on the platform's `ps`.
- **The CLI's deck-load line now passes through the message routine.** If any journal on the Mac
  shows that line twice, or not at all, send its transcript.

## Linux results (CC; a clone path of 51 characters, no short link, `MACOS_HOME` unset)

| Check | Result |
|---|---|
| Engine and `smacos_dvr` rebuild, gfortran and ifx; GMI and mmacos mex | all succeed |
| `macos` started in an empty directory, ifx and gfortran | the parameter-file line is one line, identical between the compilers |
| `MACOS_HOME` naming a directory that does not exist | the warning is one line per message, then the engine carries on |
| CLI command journals from the long path, ifx | 122 of 122 pass. Before this change 11 of 122 failed; with only the `param_mod.F` half, 9. |
| CLI command journals from the long path, gfortran | 122 of 122 pass (121 before the `ps -ww` change) |
| CLI load gate, gfortran and ifx | 166 of 166; no deck changed |
| MACOS manual dialogs regenerated from the long path, ifx | every listing and dialog reproduces with no change, so no printed text moved for the paths the manual shows |
| mmacos fast suite | 816 pass, 0 fail, 0 incomplete |
| pymacos build and suite | 7,030 pass, 0 fail, 1 skipped |

The mmacos and pymacos runs were made on the `param_mod.F` change, which is the only part the libraries contain. The command-line program's line was converted afterwards, and the command-line checks in the table were all repeated on that final build.

## Reporting

As before. Open the PR from `mr_merge_param_print` with base `mr_merge_manuals`, with the text
below and your Mac line added. No issue number.

## PR text

```markdown
# Print long paths on one line

Branch `mr_merge_param_print`, on top of `mr_merge_manuals`. No issue number.

## What changes for a user

With the ifx build in a clone whose path is long, the engine printed two start-up and load messages across two or three lines:

      Using parameter file "
     /home/.../data/macos_param.txt"

      Input file
     /home/.../tests/ZGD_test_files/opt_ex
     ample.in being loaded.

They now print on one line each, as the gfortran build always printed them. For a path that fitted in 80 columns nothing changes: the text is the same, character for character.

## Why

Both messages were written with list-directed output (`print *`, `WRITE(*,*)`) straight to the terminal. The Intel compiler starts a new line there for an item that does not fit in 80 columns, and splits an item longer than a line.

This is also why the command-line journals had to be run from a short path with the ifx build: nine of the 122 could not find a deck's "being loaded" line when the clone path was long.

## What is in the commit

- **`src/macosio.F`**: the command-line program's "Input file ... being loaded." line goes through the engine's message routine, as the library's copy of the same line already did. That routine exists for this: it writes the line to an internal record, where the 80-column rule does not apply, and prints the record. The change was made with the project's conversion tool, `tools/emsg_convert.py`, which gains one line in its table.
- **`src/param_mod.F`**: the parameter-file line, the two lines of the stale-`MACOS_HOME` warning and the "build's data directory" line use explicit formats, which are not wrapped. They are not sent through the message routine, because they print at start-up in every MATLAB and Python session and would add a line to each session's message buffer.
- **`src/notes/api_bindings.md`** records both, and why they differ.
- **`tests/cli/cmd/STARTUP_size.jou`**: the journal found its own binary with `ps`, which cuts the command line at the terminal width, 80 columns under the test runner. With gfortran in a long path the binary was then not found. It now uses `ps -ww`.

## Test results (Linux)

All runs from a clone whose path is 51 characters long, with no short link, `MACOS_HOME` unset.

| Check | Result |
|---|---|
| Engine and `smacos_dvr` rebuild, gfortran and ifx; GMI and mmacos mex | all succeed |
| `macos` started in an empty directory, ifx and gfortran | the parameter-file line is one line, identical between the compilers |
| `MACOS_HOME` naming a directory that does not exist | the warning is one line per message, then the engine carries on |
| CLI command journals from the long path, ifx | 122 of 122 pass. Before this change 11 of 122 failed; with only the `param_mod.F` half, 9. |
| CLI command journals from the long path, gfortran | 122 of 122 pass (121 before the `ps -ww` change) |
| CLI load gate, gfortran and ifx | 166 of 166; no deck changed |
| MACOS manual dialogs regenerated from the long path, ifx | every listing and dialog reproduces with no change, so no printed text moved for the paths the manual shows |
| mmacos fast suite | 816 pass, 0 fail, 0 incomplete |
| pymacos build and suite | 7,030 pass, 0 fail, 1 skipped |

The mmacos and pymacos runs were made on the `param_mod.F` change, which is the only part the libraries contain. The command-line program's line was converted afterwards, and the command-line checks in the table were all repeated on that final build.

Not run for this change: the GMI regression (the library's numerical code is untouched), the large-model mmacos groups.

## Not covered

- **Mac and Windows.** Not run here. gfortran did not wrap these lines, so the Mac build should print what it printed before.
- **Other list-directed messages that carry a path.** These two are the ones the journals and the manual's dialogs read. A search for others was not made.
```
