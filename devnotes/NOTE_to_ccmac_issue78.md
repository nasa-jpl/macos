# NOTE to CCMac — Mac test of issue #78, branch `mr_merge_devnotes` (CC, 2026-11-17)

Thank you for the #77 run: the stale build tree healed on the Mac as hoped, and Apple's Python 3.9
runs the generator. The one failure you reported is fixed on the #77 branch itself (`c36c542d`, the
line-scale comparison), which this branch sits on.

## First: where the notes are now

**This branch moves the developer notes into `devnotes/`**, this note among them. From here on a
note to you is `devnotes/NOTE_to_ccmac_issue<N>.md`, not a file at the top of the repository. The
three earlier ones (`issue75`, `issue76`, `issue77`) moved there too. `devnotes/NOTE_to_ccmac_pr_bases.md`
is the note on how to open the PRs stacked.

Then, as before: `git fetch origin && git switch mr_merge_devnotes`. Keep your build directories;
nothing compiled changed. Keep `MACOS_HOME` unset and `MACOS_resources` out of the way.

## What the branch changes

1. **140 top-level entries (145 files) moved to a flat `devnotes/`** by `git mv`: every `BRIEF_*`,
   `NOTE_*`, `PLAN*`, `RUNBOOK_*`, `REVIEW_*`, `REPORT_*`, `DRAFT_*` file, `CURRENT_SLICE.md`,
   `DEV_FILES.md`, `ENGINE_ISSUES.md` and the `NOTES/` directory. The top of the repository goes
   from 176 entries to 36.
2. **`release-exclude.txt`** has one line, `devnotes`, in place of 21; the release strip removes
   the same 570 files as before.
3. **Stays at the top:** `CLAUDE.md` (it ships), `release-exclude.txt` (the strip tool reads it
   there) and `demo_session/` (issue #82 removes it).
4. **Pointers:** the root `CLAUDE.md` now says `devnotes/PLAN.md`, `devnotes/PLAN_CONSOLIDATION.md`,
   `devnotes/PLAN_DESIGN_LAYER.md`, `devnotes/CURRENT_SLICE.md`; `BRANCHES.md`, the mmacos and pymacos
   cheatsheets and the README table follow. A note cited elsewhere by its bare name is unchanged;
   look for it in `devnotes/`.

## What to run

Nothing here is compiled or executed, so a short run is enough:

1. `python3 tools/release_strip.py | head -1` should say 572 files to strip (the 570 from before the
   move and the two notes this commit adds).
2. `./run_mmacos_tests.sh tSegMirMaker` (the #77 fix, which this branch carries): 4 pass expected.
3. The CLI load gate and `./run_mmacos_tests.sh quick`, as a check that nothing read a moved file.
4. Anything of yours that opens a note by path (`PLAN.md`, `CURRENT_SLICE.md`, a `NOTE_*`): it is in
   `devnotes/` now. Say which, if any, you had to change.

A full fast suite is welcome but not asked for.

## Linux results (CC, `MACOS_HOME` unset, a home directory holding only this clone)

| Check | Result |
|---|---|
| Release strip (`tools/release_strip.py`, dry run) | removes the same 570 files as before the move, plus the two notes this commit adds (572) |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the #76 records |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| mmacos fast suite | 816 pass, 0 fail (the 815 plus the comparator test added in #77) |
| GMI regression, gfortran | 6 of 6 pass |

Not run for this change: the ifx GMI mex, the mmacos model-256, 512 and 1024 groups, the manual builds.

## Reporting

As before. Open the PR from `mr_merge_devnotes` with base `mr_merge_data` (see the PR-bases note),
with the text below and your Mac line added.

## PR text

```markdown
# Move the developer notes into devnotes/

Closes #78. Branch `mr_merge_devnotes`, on top of `mr_merge_data` (#77).

## What changes for a user

Nothing at build time or run time. The top of the repository goes from 176 entries to 36: the 140 developer-only entries (145 files: plans, briefs, notes, reviews, runbooks, drafts, reports, and the `NOTES/` directory) are now in one directory, `devnotes/`. None of them ships in a release, before or after.

## What is in the commit

- **The move.** `git mv` for the 140 entries, into a flat `devnotes/`. Every file arrived with the same content, mode and name (blob hashes compared before and after). The notes cite one another by bare file name, so a flat directory keeps those citations pointing at a neighbour.
- **`release-exclude.txt`.** One new line, `devnotes`, replaces 21 lines that named top-level files. The patterns that also match files deeper in the tree (`BRIEF_*.md`, `STATUS_*.md` and so on) stay. The release strip removes exactly the same 570 files as before the move, compared name by name.
- **What stays at the top:** `release-exclude.txt` itself, because the strip tool and the release steps read it there; `CLAUDE.md`, which ships and must be at the top to be read; and `demo_session/`, which issue #82 removes.
- **Pointers.** The root `CLAUDE.md`, `BRANCHES.md`, the mmacos and pymacos cheatsheets and the README's directory table say `devnotes/`. Six links inside `PLAN.md` to engine files gain `../`. A new note goes in `devnotes/`.
- **Citations elsewhere are unchanged.** About 60 source files and 70 documents cite a note by its bare name in a comment (`PLAN.md §0`, `REVIEW_POL_SP_SIGN_2026-07-27.md`). None of them reads the file; I checked every script, Makefile, tool and test for an open or a path built from one of the 140 names and found none. The names still identify the files, now in `devnotes/`.

## Test results (Linux)

All runs had `MACOS_HOME` unset, under a home directory that holds only this clone. Nothing in this change is compiled, so the binaries are those of #77.

| Check | Result |
|---|---|
| Release strip (`tools/release_strip.py`, dry run) | removes the same 570 files as before the move, plus the two notes this commit adds (572) |
| CLI load gate, gfortran and ifx | 166 of 166 decks load; no deck changed against the #76 records |
| CLI command journals, gfortran and ifx | 122 of 122 pass |
| pymacos suite | 7,030 pass, 0 fail, 1 skipped |
| mmacos fast suite | 816 pass, 0 fail (the 815 plus the comparator test added in #77) |
| GMI regression, gfortran | 6 of 6 pass |

Not run for this change: the ifx GMI mex, the mmacos model-256, 512 and 1024 groups, the manual builds.

## Not covered

- **Mac and Windows.** Not run here. Nothing in this change is compiled or executed.
- **Developer notes below the top level** (for example `mmacos/PLAN_EXAMPLES_REORG.md`, the template status files) stay where they are; the issue is about the top-level files, and the strip list still removes them by pattern.
```
