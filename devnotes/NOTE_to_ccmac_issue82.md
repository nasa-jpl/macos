# NOTE to CCMac — Mac check of issue #82, branch `mr_merge_remove_scratch` (CC, 2026-11-18)

The last step of the reorganization: `scratchpad/` and `demo_session/` are removed. Read this note
from `devnotes/` in your clone.

## First: switch; nothing is rebuilt

1. `git fetch origin && git switch mr_merge_remove_scratch` (it sits on top of `mr_merge_src`).
2. **Keep every build.** No source, build file, test or prescription changed.
3. If the switch leaves a `<clone>/demo_session/` or `<clone>/scratchpad/` directory behind, it
   holds only untracked files of yours. Delete it by full path, or keep what you want from it first.

## What the branch changes

1. **`git rm -r scratchpad demo_session`**: 342 files, 137 MB of working tree. They stay in the
   history.
2. **`release-exclude.txt`**: the `demo_session` line is dropped.
3. **Three pointers**: `src/notes/opd_reference.md` names `tests/cli/cmd/SPCENTER.jou` as the
   standing gate for the SPOT centring (the scratch script was its first form); one comment in a
   `dyson5` script and one line in the `pdi_dm96` README no longer point into `demo_session/`.

## What to run

A light check is enough, since nothing that is built or tested changed:

1. `git status` is clean after the switch, and `git ls-files | wc -l` gives 11,941.
2. `tests/cli/run_cli_tests.sh load`: 166 of 166, no deck changed.

Run more if you want the full gate on the final tree; I would, once, on the branch that ends the
series.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Tracked files | 11,941 with the note for this step, down from 12,282 |
| Release strip | the same files as before, less `demo_session/` |
| Deck lists written by the runner | the same 166 core decks and 1,521 in all |
| CLI load gate, gfortran and ifx | 166 of 166 decks load with each; no deck changed against the #81 records, and all 166 were compared |

## What follows this branch

- **GMI on the Mac**: a branch `mr_merge_gmi_mac` on top of this one, with its own note. Dave has
  asked for it; it answers your two blockers (`uname -i`, `MATLAB_DIR`) and the rest of the
  Linux-only Makefile.
- **The manuals**: the three PDFs are rebuilt once on the final tree, on a branch of their own.

## Reporting

As before. Open the PR from `mr_merge_remove_scratch` with base `mr_merge_src` (see
`devnotes/NOTE_to_ccmac_pr_bases.md`), with the text below and your Mac line added.

## PR text

```markdown
# Remove scratchpad/ and demo_session/

Closes #82. Branch `mr_merge_remove_scratch`, on top of `mr_merge_src` (#81).

## What changes for a user

Nothing. Neither directory was part of a release, and no code, test, build file or manual reads either of them.

| Removed | Tracked files | Size | What it held |
|---|---|---|---|
| `demo_session/` | 341 | 137 MB | the working area for presentations: slide sources, 28 PowerPoint files, 240 figures, the slide builder |
| `scratchpad/` | 1 | 4 kB | one test script, `spot_cli_gate.py` |

A clone is 137 MB smaller to check out. The files stay in the history of the repository.

## What was checked before removing

- **Nothing reads them.** A search of every tracked file outside the two directories found no code, test, build file or tool that opens a path in either. Three comments and two dated notes name a file in them.
- **The one test script is covered.** `scratchpad/spot_cli_gate.py` was the first check that the command-line SPOT is centred about the element after a LOAD. The command journal `tests/cli/cmd/SPCENTER.jou` asserts the same thing in the standing suite, on both compilers.
- **The presentation material is kept elsewhere** by its authors; it is not something a user of the code needs.

## What is in the commit

- **The removal.** `git rm` of the two directories, 342 files.
- **`release-exclude.txt`.** The `demo_session` line is dropped. The release strip removes the same files as before, less the 341 that are gone (compared name by name).
- **Pointers.** The engine note that named the scratch script now names `SPCENTER.jou` and says where the script went. One comment in a `dyson5` script and one line in the `pdi_dm96` README no longer point into `demo_session/`. Dated records were left as written.

## Test results (Linux)

No engine source, interface, test or prescription changed, so the binaries of #81 were used as they were.

| Check | Result |
|---|---|
| Tracked files | 11,941 with the note for this step, down from 12,282 |
| Release strip | the same files as before, less `demo_session/` |
| Deck lists written by the runner | the same 166 core decks and 1,521 in all |
| CLI load gate, gfortran and ifx | 166 of 166 decks load with each; no deck changed against the #81 records, and all 166 were compared |

Not run for this change: the command journals, the mmacos and pymacos suites, the GMI regression, the manual builds. None of them reads a removed file.

## Not covered

- **Mac and Windows.** Not run here.
- **A third-party paper** (`demo_session/responses/Jorge_VVC_SPIE_2024.pdf`) leaves the tree with this change but stays in the history, as everything removed here does.
- **The manual PDFs.** With this step the reorganization of the tree is complete, and the three manuals can be rebuilt once for all of it.
```
