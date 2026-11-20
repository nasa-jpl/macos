# NOTE to CCMac — the Windows build branch with the console and end-of-input fixes, branch `mr_merge_windows_console` (CCL, 2026-11-20)

This branch is Luis's `mr_merge_windows_build` (four commits, `a6255619`), merged with `mr_merge`
at `9f77de6f`, plus three commits of mine. **It changes engine source**, so it needs a rebuild of
the engine, the mex and pymacos before any check. It is not a pull request yet: Luis tests the two
fixes on Windows first. What is asked of the Mac is below; the main question is the last commit,
which exists because gfortran behaves differently from ifx and I had tested only ifx.

## What the branch changes, against `mr_merge`

| Commit | What | Files |
|---|---|---|
| `b588adab`, `0e45bcdc`, `6fd47108`, `a6255619` (Luis) | Windows: pymacos from the top-level CMake, `--lower` for f2py **on every OS**, `/heap-arrays`, two tests made independent of the path separator and of `/tmp` | `CMakeLists.txt`, `app/pymacos/src/cmake/CMakeLists.txt`, `test_sint.py`, `test_zrngr_roundtrip.py` |
| `9d47df59` | A build **without readline** stops at end of input instead of prompting forever | `src/macosio.F`, `tests/cli/tools/eof_gate.sh` |
| `871683a8` | `MACOS_CONSOLE`: where the engine's console lines go inside a host | `src/emsg_mod.F90`, `src/macos_api_mod.F90`, `CMakeLists.txt` (an ifx option), `test_console.py` |
| `2eda4ec7` | The same, made safe with gfortran | `src/emsg_mod.F90`, `test_console.py`, `tests/cli/tools/console_gate.sh` |

No generated binding file changes: `MACOS_CONSOLE` is an environment setting read at the api's
`init`, not an api routine. The notes are `src/notes/api_bindings.md` (the `MACOS_CONSOLE` entry)
and `src/notes/cli.md`.

## Why the last commit, and why the Mac matters

`MACOS_CONSOLE` = unset or `on` (the console), `off` (the null device), or a file name. It works by
reconnecting Fortran unit 6. On Linux with gfortran 13.3 I found:

- gfortran's runtime **ignores** a failed console write; ifx ends the process. So with gfortran
  the setting only quiets or logs the engine.
- Reconnecting unit 6 when the host had **already closed** its output faulted (SIGSEGV), until
  unit 6 was closed first.
- A file that could not be opened left `WRITE(*)` going to a new file `fort.6`, until the file was
  opened first on a unit of its own.

These are runtime behaviors, and the Mac's gfortran is a different version on a different
system. The gate below runs them directly.

## Checks to run on the Mac

1. **Build** the engine (gfortran), the mex and pymacos from this branch. Luis's CMake changes
   touch the top-level file and the pymacos one; `--lower` is now passed to f2py on every OS.
2. **`tests/cli/tools/console_gate.sh gfortran`** from the repository root. It links
   `src/emsg_mod.F90` into a ten-line driver and runs ten legs. On Linux all ten pass. If a leg
   fails, send its line as printed; "dead (want live)" on the two legs with the console closed and
   the setting unset or `on` would mean your gfortran ends the process on a failed write, which
   is a fact worth having and only needs the gate's expectation changed.
3. **`python -m pytest app/pymacos/tests/test_console.py`**: 9 tests. Linux (ifx library): 9 pass.
   - The negative control reads which runtime the library is linked to and expects survival from
     a gfortran library.
   - Three legs read the engine's lines from the child's stdout. If those lines reach no readable
     stream on the Mac (as with `capfd` in `test_messages`), those legs skip themselves; say how
     many skipped.
4. **pymacos batch A** as usual. Linux: 398 pass (389 before the branch, and 9 new).
5. **mmacos**: `tMessages` and `tMacosPkg`, and `tMacosPkg` once more with `MACOS_CONSOLE=off` in
   the environment. Linux: 3, 25 and 25 pass.
6. **The CLI gate** (`tests/cli/run_cli_tests.sh full`). Linux, before the last commit, which the
   CLI does not reach: load 167 of 167 on both compilers with every number equal to the record at
   `b5b37c34`; command journals 122 of 122 with gfortran; transcripts 30 of 30.
7. **Optional, `tests/cli/tools/eof_gate.sh`**: it needs a macos built with
   `-DMACOS_FORCE_NO_READLINE=ON` (the build line is in the script's header). Linux: three legs
   pass. The Mac's normal build has readline and is not changed by that commit.

## Known, not from this branch

From a clone with a long path (91 characters to the deck, on Linux) the ifx command journal
`VALIDATE.jou` fails: ifx breaks the line that prints the deck's path at 80 columns. An ifx build from before
this branch does the same. It does not affect gfortran.
