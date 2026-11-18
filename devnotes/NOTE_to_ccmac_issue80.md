# NOTE to CCMac — Mac test of issue #80, branch `mr_merge_tests` (CC, 2026-11-18)

The engine's own tests move under `tests/`: `cli_tests/` becomes `tests/cli/` and `ZGD_test_files/`
becomes `tests/ZGD_test_files/`. The package tests stay where they are. Read this note from
`devnotes/` in your clone.

**Order of the stack**: #79 (`mr_merge_app`), then #87 (`mr_merge_tools`, the `tools/` move, note
`devnotes/NOTE_to_ccmac_issue87.md`), then this branch. If you have not done #87, do it first.

## First: switch; nothing is rebuilt

1. `git fetch origin && git switch mr_merge_tests` (it sits on top of `mr_merge_tools`).
2. **Keep every build**: the engine, the mex, pymacos, SegMirMaker. No source or build file changed.
3. **Delete the leftover `<clone>/cli_tests/` directory by full path**, if the switch leaves one.
   Git moves the tracked files only; your untracked records and logs stay behind in it. Move any
   record you want to keep into `tests/cli/records/` first.
4. Keep `MACOS_HOME` unset and `MACOS_resources` out of the way.

## What the branch changes

1. **The move** (`git mv`, content unchanged). The 94 prescriptions are byte for byte the same.
2. **The deck lists are no longer tracked, and you no longer regenerate them.**
   `tests/cli/run_cli_tests.sh` calls `make_corpus.sh` at the start of every `load`, `corpus` and
   `extended` run. Drop the "run `make_corpus.sh` first" step from your gate. `make_corpus.sh`
   writes the extended list only when called with the argument `extended`.
3. **`prev_record.py`** reads a record's history by name in any directory, with renames followed,
   so the pick of the previous record survives the move.
4. **`judge_load.py`** prints how many decks are in both records, and a note when some are not.
   Expect that note on your first run: your earlier records name the decks under
   `<clone>/ZGD_test_files/`, so 75 of the 166 have no partner. It is information, not a failure.
5. **`.gitattributes`**: the three line-ending rules are re-anchored
   (`tests/ZGD_test_files/tst_eol_*.in`, `tests/ZGD_test_files/tst_short_coef_cr.in`,
   `app/pymacos/tests/Rx/Rx_Eol_*.in`). The pymacos line had not followed pymacos into `app/`.
6. **Journals**: 57 command journals and 3 manual-transcript journals name `$ROOT/tests/...`.
7. **Scripts in `tests/cli/`** climb one level more to the repository root.

## What to run

1. `tests/cli/run_cli_tests.sh load` and `cmd`, with no list regeneration beforehand.
2. The mmacos fast suite and `tCalib`.
3. The pymacos suite.
4. The line-ending check below.

## Mac-specific things I could not check

- **Line endings on your checkout.** Please run
  `git ls-files --eol -- 'tests/ZGD_test_files/tst_eol_*' tests/ZGD_test_files/tst_short_coef_cr.in 'app/pymacos/tests/Rx/Rx_Eol_*'`
  and say what it prints. Every one of the nine lines should end `attr/-text`, and the `crlf`
  files should show `i/crlf w/crlf`.
- **`make_corpus.sh` on macOS on every run**: it is the script that once aborted there under
  `set -euo pipefail`. It gained an optional argument and two lines at the end; say if the
  runner stops with "make_corpus.sh failed".
- **The command journals with a longer path.** Each deck path is six characters longer
  (`tests/`). On Linux the ifx journals need a short clone path because the engine wraps long
  lines; if a journal fails on the Mac only at a wrapped path, that is the cause.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Deck lists written by the runner | the same 166 core decks and 1,521 in all as before the move |
| CLI load gate, gfortran and ifx | 166 of 166 decks load with each; matched across the moved path, all 12 result columns of every deck equal the record taken before the move |
| CLI command journals, gfortran and ifx | 122 of 122 pass with each; the 3 must-fail journals still fail |
| mmacos fast suite | 816 pass, 0 fail, 0 incomplete |
| `tCalib` (outside the fast suite; it reads a `ZGD_test_files` deck) | 9 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| `prev_record.py` across a move, scratch clone | same pick as before the move for all 10 kinds |
| Fresh clone, load gate | a fresh clone has no deck lists; the runner writes them (166 core) and the gfortran load gate gives 166 of 166; the previous record it picks is the latest committed one; the three line-ending rules apply |

## Known, and not from this branch

- Three comment lines in the engine source still say `cli_tests/cmd/...`; they are corrected
  with #81.
- Records are matched by absolute path, so a comparison across clones is partial; the gate now
  says so.
- About 60 templates and a few tests still default to `~/dev/MACOS_resources/...`.
- The gfortran `macos` spins at 100% CPU at the prompt when its input is a pipe.

## Reporting

As before. Open the PR from `mr_merge_tests` with base `mr_merge_tools` (see
`devnotes/NOTE_to_ccmac_pr_bases.md`), with the text below and your Mac line added.

## PR text

```markdown
# Gather the engine's tests under tests/

Closes #80. Branch `mr_merge_tests`, on top of `mr_merge_tools` (#87).

## What changes for a user

The two top-level test directories move under `tests/`:

| Was | Is | Tracked files |
|---|---|---|
| `cli_tests/` | `tests/cli/` | 737 |
| `ZGD_test_files/` | `tests/ZGD_test_files/` | 94 |

After pulling this change:

- **Run the command-line tests from `tests/cli/`**: `tests/cli/run_cli_tests.sh load`, and so on.
- **There is no "regenerate the deck lists" step any more.** The runner writes them itself on every run (see below).
- **Personal scripts and journals** that name `ZGD_test_files/...` or `cli_tests/...` need `tests/` in front, and `cli_tests` becomes `tests/cli`.
- Nothing is rebuilt. No engine source, interface or prescription changed.

## Why only these two move

`cli_tests/` and `ZGD_test_files/` are the engine's own tests. They belong to no package, and they sat at the top of the repository.

The other tests stay where they are, for now:

| Tests | Tracked files | Why they stay |
|---|---|---|
| `app/mmacos/tests/` | 249 | Run by `app/mmacos/run_mmacos_tests.sh`, which expects them beside the package; the test helpers find the repository and the fixtures from their own position. |
| `app/pymacos/tests/` | 196 | Run by `pytest` from the package; `conftest.py` and `repo_paths.py` find the package and the engine from their own position. |
| `app/GMI/regression/` | 19 | Its script and reference files sit beside the GMI sources they test. |
| `tools/rx_converter/*/tests/` | 21 | Each converter carries its own tests and fixtures. |

Moving them would separate each interface from its tests, change the command every user of those runners types, and need a second round of path edits like #79, with no gain to a user. A package that carries its own tests is also what a Python or MATLAB user expects to find. This can be revisited; nothing here prevents it.

## What is in the commit

- **The move.** `git mv` for the two directories. Every moved file arrived with the same content, mode and relative path (blob hashes compared before and after). The 94 prescriptions in `tests/ZGD_test_files/` are byte for byte unchanged.
- **The deck lists are no longer tracked.** `corpus_core.txt` and `corpus_committed.txt` hold absolute paths. The committed copies named the decks of the two old clones on one machine. On any other machine they named files that do not exist, and on that one machine the load gate tested the old clones and reported green without touching the branch. `run_cli_tests.sh` now calls `make_corpus.sh` at the start of every `load`, `corpus` and `extended` run, and the two lists are ignored by git. The sandbox scan for the extended list runs only when it is asked for.
- **The load gate picks the right previous record after a move.** `prev_record.py` dates a record by the commit that added it. A move makes git report every record as added by the moving commit, so all dates tied and the pick fell to the alphabetically last name. It now reads the history of the record's name in any directory, with renames followed. Tested in a scratch clone: after the move the old script picks a different record for 8 of 10 kinds; the new one picks the same record as before the move for all 10.
- **The load gate says how much it compared.** Records match decks by full path, so a deck under a moved directory, or in another clone, has no partner in the previous record and nothing is compared for it. `judge_load.py` now prints how many decks are in both records, and a note when that is fewer than all of them. It was silent before.
- **Line-ending fixtures protected again.** `.gitattributes` exempts the CR and CRLF test prescriptions from line-ending conversion by path. The two `ZGD_test_files` lines are re-anchored under `tests/`. The `pymacos/tests/Rx/Rx_Eol_*.in` line had not followed pymacos into `app/` in #79, which was a miss in #79; it is corrected here. The stored bytes of all nine fixtures were checked and are intact.
- **Scripts that find the repository from `tests/cli/`** climb one level more: `run_cli_tests.sh`, `make_corpus.sh`, `cli_cmd_gate.py`, `tools/hzdef_check.py`.
- **Journals.** 57 command journals and 3 manual-transcript journals name `$ROOT/tests/ZGD_test_files/...` or `$ROOT/tests/cli/...`. Nothing else in them changed.
- **Interface tests that read these directories**: three mmacos test classes and two pymacos test files build the new path; the others changed only in comments.
- **Manual tools and documents.** The command-reference generator and five other manual tools read `tests/cli`; the generated command reference, the engine notes (`macos_f90/notes/`), the README, the Mac guide and the manual sources name the new paths. Dated records were left as written.

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone, on the engine binaries of #79 (the engine source is unchanged).

| Check | Result |
|---|---|
| Deck lists written by the runner | the same 166 core decks and 1,521 in all as before the move |
| CLI load gate, gfortran and ifx | 166 of 166 decks load with each; matched across the moved path, all 12 result columns of every deck equal the record taken before the move |
| CLI command journals, gfortran and ifx | 122 of 122 pass with each; the 3 must-fail journals still fail |
| mmacos fast suite | 816 pass, 0 fail, 0 incomplete |
| `tCalib` (outside the fast suite; it reads a `ZGD_test_files` deck) | 9 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| `prev_record.py` across a move, scratch clone | same pick as before the move for all 10 kinds |
| Fresh clone, load gate | a fresh clone has no deck lists; the runner writes them (166 core) and the gfortran load gate gives 166 of 166; the previous record it picks is the latest committed one; the three line-ending rules apply |

Not run for this change: an engine rebuild (no source or build file changed), the GMI regression, the manual builds and the manual-transcript journals, the `corpus` and `extended` runs of the load gate.

## Not covered

- **Mac and Windows.** Not run here.
- **Three comment lines in the engine source** (`macos_cmd_loop.inc`, `tracesub.F`) still say `cli_tests/cmd/...`. Changing them would mean rebuilding and re-gating the engine for a comment; they are corrected with #81, which renames the engine directory.
- **Comments inside 27 prescriptions** still say `cli_tests/...`. The prescriptions are test fixtures and were kept byte for byte.
- **The manual PDFs** are rebuilt once, after #82. The generated command reference was corrected by rule, not regenerated.
- **Records are still matched by absolute path**, so a comparison across clones or across a moved directory is partial. The gate now says so; making records independent of the clone path is a change of its own.
```
