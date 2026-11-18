# NOTE to CCMac — open the reorganization PRs stacked, each on the branch before it (CC, 2026-11-17)

Dave has asked you to post the pull requests in parallel. Please give each one the previous issue
branch as its base, not `mr_merge`.

## Why

The branches are stacked: each starts from the one before. A PR for #76 opened against `mr_merge`
also shows #75's commit, and a PR for #77 shows three issues' commits. Reviewers would see the same
5,571 renames in three PRs. With the previous branch as the base, each PR shows only its own change.

## The bases

| Issue | Head branch | Base branch |
|---|---|---|
| #75 | `mr_merge_macos_home` | `mr_merge` (PR #84, already open; leave it) |
| #76 | `mr_merge_external` | `mr_merge_macos_home` |
| #77 | `mr_merge_data` | `mr_merge_external` |
| #78 | `mr_merge_devnotes` | `mr_merge_data` |

Later issues continue the same way.

## Commands

To open one:

    gh pr create --base mr_merge_macos_home --head mr_merge_external --title "<title>" --body-file <file>

To change one that is already open against `mr_merge`:

    gh pr edit <number> --base <previous issue branch>

Check the result with `gh pr view <number> --json baseRefName,files -q '.baseRefName, (.files|length)'`.
The file count should match that issue's own commit: 127 for #75, 5,662 for #76, 27 for #77.

## In each PR's text

Add one line at the top of the body, under the title:

    Stacked on #<previous PR number>. This PR shows only its own change; merge #<previous> first.

## Merging, for whoever merges

- **Merge in order, #84 first.**
- **Use "Create a merge commit".** Do not squash and do not rebase-merge: either one rewrites the
  commits the next branch is built on, and the next PR then shows conflicts or the earlier change
  a second time.
- **Delete the merged branch from the PR page** ("Delete branch") after the merge. GitHub then
  changes the base of the next PR to `mr_merge` by itself. Do not delete a branch any other way
  while a PR is based on it, or that PR closes.

## What does not change

The Mac gate for each branch, the PR text in each note, and your Mac line in the test results.
