# Runbook: the October 2026 release moves (for Andy; drafted by CC for Dave, 2026-10-21)

Six steps, in this order, nothing run before the Friday meeting approves it.  Every count below was re-measured
on 2026-10-21 (`git log --all --diff-filter=A --name-only`): macos history holds 423 NPSOL paths and 1564 pgplot
paths; MACOS_resources history holds 21 OPTIIX/IRIS prescription paths (0 at HEAD); `.git` is 550 MB (macos) and
2.2 GB (resources).  Who runs what: steps 1, 2 and 6 need admin on nasa-jpl (Andy); steps 3-5 are ordinary pushes
Dave can make after Andy's go.  The branch rules are `BRANCHES.md`; the strip lists are `release-exclude.txt` in
each repo (amended 2026-10-21, reviewed in `RELEASE_EXCLUSIONS_REVIEW.md`).

> **STATUS after the group meeting of 2026-10-26 (Dave, Andy and the group).**  The room set the path.  It keeps
> the technical content below as the checklist, but changes WHO runs the moves and HOW the two repositories join:
>
> 1. **We consolidate** macos and MACOS_resources first: orphaned files removed, documentation finished to release
>    standard.  `PLAN_CONSOLIDATION.md` section 3b lists what "done" means for this announcement.
> 2. **Andy writes his proposed steps on the GitHub wiki**: promote `dev-candidate` to `dev`, then `dev` to `main`
>    with the redactions we specified (`release-exclude.txt` in both repos, `RELEASE_EXCLUSIONS_REVIEW.md`); history
>    preserved as far as possible while the NPSOL, OPTIIX and IRIS sources are removed from it; the macos and
>    MACOS_resources histories MERGED (not step 5's snapshot).  **We review the wiki page** against steps 1-6 below.
> 3. **Andy executes** after Dave announces that 1 and 2 are done.
> 4. **We test the new repository**: every clone deleted and re-cloned first (the re-clone block at the end of
>    step 2), then the gate list of step 4 run on the new clone.
> 5. When it passes, **Andy makes a new release.**
>
> So steps 1, 2, 3, the merge-and-strip of 4, 5 and 6 are Andy's; ours are the re-clone and the gates.  What to look
> for when we review his wiki page: (a) the step-1 archive bundle exists and is verified BEFORE any rewrite; (b) the
> removal lists are drawn from history (step 2's `git log --all --diff-filter=A` recipe) and checked on EVERY ref
> of BOTH repos before the histories are joined -- joining first would carry a removed path across; (c) a merged
> history means the resources tree moves under a subdirectory of macos (`git filter-repo --to-subdirectory-filter`
> on the rewritten resources mirror, then `git merge --allow-unrelated-histories`), and the 153 path fix-ups of
> step 5 are still needed afterwards, as one gated commit, whoever makes it; (d) the strip on `main` uses the two
> `release-exclude.txt` lists (folded into one once the trees are joined), by `tools/release_strip.py` or an
> equivalent that is idempotent; (e) protection on `main`/`dev` off for the force-push and back on after, and the
> GitHub-support purge requested; (f) the announcement that no one pushes from a pre-rewrite clone.  Still
> undecided by the room, to ask: whether the 2267 run artifacts leave the resources history (step 2's optional
> block).
>
> **Rules in force after the release** (the room, 2026-10-26; recorded in `BRANCHES.md`): a temporary branch per
> bug fix or feature; human-readable commits (clear, concise, plain language, readable without deep MACOS
> experience); PRs reviewed, merged and their branch deleted within a week; otherwise the order of events as
> written here.

## Step 1 -- archive the full history (before anything is rewritten)
```bash
A=~/macos-archive-$(date +%Y%m%d); mkdir -p $A
for R in macos MACOS_resources; do
  git clone --mirror git@github.com:nasa-jpl/$R.git $A/$R.mirror
  git -C $A/$R.mirror bundle create $A/$R-all.bundle --all
  git bundle verify $A/$R-all.bundle          # must say "is okay"
done
du -sh $A                                      # ~2.8 GB; copy it OFF the laptop (JPL share or external drive)
```
The bundles are the only copy of the removed history after step 2.  They stay internal.

## Step 2 -- rewrite history on mirror clones, then everyone re-clones
Needs `git filter-repo` (`pip install git-filter-repo`).  Work on a SECOND mirror, never on the archive.
```bash
W=~/rewrite-$(date +%Y%m%d); mkdir -p $W; cd $W
git clone --mirror git@github.com:nasa-jpl/macos.git macos.git
git clone --mirror git@github.com:nasa-jpl/MACOS_resources.git MACOS_resources.git

# the exact paths to remove, from history itself (case-insensitive), one file per repo
cd $W/macos.git
git log --all --diff-filter=A --name-only --format= | grep -i 'npsol\|pgplot' | sort -u > ../macos-remove.txt
wc -l ../macos-remove.txt                      # 2026-10-21: 1987 (423 + 1564)
cd $W/MACOS_resources.git
git log --all --diff-filter=A --name-only --format= | grep -i 'optiix\|iris' | sort -u > ../resources-remove.txt
wc -l ../resources-remove.txt                  # 2026-10-21: 21 -- READ the list: 'iris' must not catch an innocent name
```
Review both lists by eye (the resources one is short; the macos one is `npsol/`, `macos_f90/npsol/`, the netlib
blas/lapack that came with it, and every `pgplot` path).  If the room also votes to drop the 2267 run artifacts
from resources history, add them to `resources-remove.txt` first:
```bash
cd $W/MACOS_resources.git
git ls-tree -r --name-only dev-candidate | grep -E '\.(mat|png|log|pptx)$' | grep -v '^pymacos/tests/\|^mmacos/tests/' >> ../resources-remove.txt
sort -u -o ../resources-remove.txt ../resources-remove.txt
```
(exclude the test fixtures' own .mat/.png by that filter; check the list before going on -- the runners regenerate
every artifact, the fixtures they do not).  Then:
```bash
for R in macos MACOS_resources; do
  cd $W/$R.git
  git filter-repo --invert-paths --paths-from-file ../${R/MACOS_resources/resources}-remove.txt --force
  # verify: every ref reads 0 for the removed names
  for ref in $(git for-each-ref --format='%(refname)'); do
    n=$(git log $ref --diff-filter=A --name-only --format= | grep -ci 'npsol\|pgplot\|optiix\|iris'); echo "$ref $n"
  done | awk '$2!=0{bad=1; print} END{exit bad}' && echo "$R: clean on every ref"
  git reflog expire --expire=now --all && git gc --prune=now --aggressive && du -sh .
done
```
Force-push the rewritten mirrors.  On GitHub first: Settings -> Branches -> disable protection on `main`/`dev`
(re-enable after); note any open PRs (their refs/pull/* keep old objects until GitHub support purges them -- ask
support to run the sensitive-data purge for both repos, and to clear forks' caches).
```bash
cd $W/macos.git           && git push --mirror --force
cd $W/MACOS_resources.git && git push --mirror --force
```
Then the announcement: **everyone deletes their clones and re-clones; no push from a pre-rewrite clone, ever**
(a push from an old clone re-introduces the removed history -- that is what happened to the July rewrite).  Dave's
machine today holds three independent clones (`~/dev/macos`, the shared working tree; `~/dev/macos_cli_base`, the
clean gate tree; `~/dev/macos_to_hex`, TO's one-off gate tree): delete all three and re-clone ONCE, as one clone
plus a worktree that shares its objects -- the gate tree then cannot drift from the repository it gates, and the
third tree is not recreated (a second gate tree, if two lanes ever gate at the same hour again, is a one-line
`git worktree add` at the time):
```bash
# KEEPERS FIRST (Dave 2026-10-27, "we'll have to rethink before recloning"): the clones hold files that
# the new repository will not -- every untracked file (the BRIEF_/DRAFT_/REPORT_ drafts at the macos root,
# scratch runs, Dave's edit copies) and, on the release surface, everything the two release-exclude.txt
# lists strip (the six template REPORT_*.md among them, Dave 2026-10-27: not shipped, kept locally).
# Copy them out before any rm, then decide file by file what comes back into the re-cloned tree.
K=~/dev/keep_$(date +%Y%m%d); mkdir -p $K/macos $K/MACOS_resources
(cd ~/dev/macos && git ls-files --others --exclude-standard -z | rsync -a0 --files-from=- . $K/macos/)
(cd ~/dev/MACOS_resources && git ls-files --others --exclude-standard -z | rsync -a0 --files-from=- . $K/MACOS_resources/)
(cd ~/dev/MACOS_resources && git ls-files -z -- 'mmacos/templates/*/*/REPORT_*.md' | rsync -a0 --files-from=- . $K/MACOS_resources/)
du -sh $K                                      # read the listing before the rm below
rm -rf ~/dev/macos ~/dev/macos_cli_base ~/dev/macos_to_hex ~/dev/MACOS_resources     # after step 1's bundle is safe AND $K is read
git clone git@github.com:nasa-jpl/macos.git ~/dev/macos && cd ~/dev/macos && git checkout dev-candidate
git worktree add ../macos_cli_base dev-candidate       # the clean gate tree; run_cli_tests.sh's default
git clone git@github.com:nasa-jpl/MACOS_resources.git ~/dev/MACOS_resources   # until step 5 folds it in
```

## Step 3 -- dev-candidate -> dev, both repos (fast-forward), ENGINE FIRST
**Precondition: step 2's re-clone is complete** -- the clone these commands run in was made AFTER the force-push
(`git log -1 --format=%H origin/dev-candidate` must match the rewritten mirror's tip, and `git log --all --diff-filter=A
--name-only --format= | grep -ci npsol` must read 0 in it).  Steps 3-6 are ordinary pushes from that clone; none of
them needs another re-clone.
```bash
for R in macos MACOS_resources; do
  cd ~/dev/$R && git fetch origin
  git checkout dev && git merge --ff-only origin/dev-candidate && git push origin dev
done
```
If `--ff-only` refuses, dev has a commit dev-candidate lacks: stop and look (`git log origin/dev-candidate..dev`).
Engine first because a resources veneer promoted ahead of its engine fails at prescription load, not at build.

## Step 4 -- dev -> main with the exclusions applied, macos first, then resources
`main` has its own history (last touched 2026-07-24 / 07-22), so it is a merge, then the strip, then one commit:
```bash
cd ~/dev/macos && git fetch origin && git checkout main
git merge --no-ff origin/dev -m "Release 2026-10: dev into main"
python3 tools/release_strip.py            # dry run: 481 files on 2026-10-21; read the list once
python3 tools/release_strip.py --apply && git commit -m "Release strip (release-exclude.txt)"
git grep -il 'npsol' -- . | wc -l         # 0 (sources); BRANCHES.md/notes may MENTION the word -- read any hit
```
The strip tool is `tools/release_strip.py` (dry run by default; a pattern is an exact path, a directory, or a glob;
a glob without `/` matches basenames anywhere, so `BRIEF_*.md` catches briefs in subdirectories).  Idempotent: the
next release merges dev again and strips again.  Resources: the same four commands in `~/dev/MACOS_resources`
with `python3 ../macos/tools/release_strip.py` (19 files on 2026-10-21).

Release gates on main BEFORE pushing it, from a CLEAN clone of the local main (`git clone ~/dev/macos /tmp/rel`):
- engine builds with ifx and gfortran (`source ./makems.sh release`, `... release gfortran`);
- `cli_tests/run_cli_tests.sh load` and `cmd` with `MACOS_CLI_TREE=.` (160/160, every journal but the two must-fail);
- mmacos `./run_mmacos_tests.sh fast` and the pymacos suite on the clean clone;
- the self-containment sweep of BRANCHES.md (no `addpath` reaching outside the repos; `HOME` moved away);
- `docs/macosMan.pdf` and `docs/macosCmdRef.pdf` present and the manual builds (`make all pdf`);
- `git ls-files | grep -c CLAUDE.md` = 5 on macos main (root, macos_f90, giza, slsqp, manual) -- the agent set ships;
  `git ls-files | grep -c 'PLAN\|CURRENT_SLICE\|BRIEF_'` = 0.
Then `git push origin main` for macos, and only then the same for resources.

## Step 5 -- import MACOS_resources into macos as a snapshot (no history)
> 2026-10-26: the room chose a MERGED history instead (Andy's wiki steps); this snapshot recipe stays as the
> fallback and for its list of directories and path fix-ups, which a history merge needs just the same.

```bash
cd ~/dev/macos && git fetch origin && git checkout -b import-resources origin/dev
S=~/dev/MACOS_resources; SHA=$(git -C $S rev-parse --short origin/dev)
for d in GMI mmacos pymacos segmirmaker optical_design rx_converter doc; do rsync -a --exclude .git $S/$d/ $d/; done
cp $S/release-exclude.txt resources-release-exclude.txt   # fold its 19 patterns into release-exclude.txt, prefixed by their directory
git add GMI mmacos pymacos segmirmaker optical_design rx_converter doc && git commit -m "Import MACOS_resources dev @$SHA as a snapshot (no history: the archive bundle holds it)"
```
Then the path fix-ups, one commit, gated: 153 tracked files in macos name `../MACOS_resources` or `MACOS_resources`
(2026-10-21: `git grep -c MACOS_resources -- . ':!*.md' ':!docs/'`).  The ones that matter at run time:
`cli_tests/cli_cmd_gate.py` (the `MACOS_RESOURCES` default), `cli_tests/make_corpus.sh`, `corpus_core.txt` /
`corpus_committed.txt` (absolute deck paths), `CMakeLists.txt` and `.gitattributes`; the `cli_tests/records/*` hits
are historical records and stay.  On the resources side after the copy: `mmacos/mmacos_setup.m`, the pymacos
`CMakeLists.txt` (`MACOS_BUILD_DIR`), `rx_fixture_path.m`, `tests/context.py` -- each assumes `../macos`; grep
`'\.\./macos'` in the imported tree.  Gate the import branch exactly as step 4's list, then PR -> dev, then dev -> main
by step 4 again (the strip list now carries the resources patterns).

## Step 6 -- archive MACOS_resources
After step 5 is on macos `main`: GitHub Settings -> Archive this repository (read-only, not deleted; keep for a few
months, then decide).  The side branches `develop`, `pol-ifo`, `tg96-oap` (4 commits not on dev-candidate) are in
the step-1 bundle; delete them before archiving or leave them -- the room's call, no code depends on them.

## Review checklist for Andy's wiki page (CC 2026-10-27; tick each against his draft)

Archive and rewrite
- [ ] A verified `git bundle --all` of BOTH repos exists before any rewrite, copied off the laptop (step 1).
- [ ] The removal lists are drawn from history (`git log --all --diff-filter=A --name-only`), reviewed by eye, and the
      `iris` match checked for innocent names (step 2).
- [ ] After the rewrite, the verification loop reads 0 for npsol / pgplot / optiix / iris on EVERY ref of EACH repo,
      BEFORE the two histories are joined.
- [ ] The 2267-class run artifacts: the room's decision (stay / leave history) is stated, and if they leave, the list
      excludes the test fixtures (`pymacos/tests/`, `mmacos/tests/`).
- [ ] Branch protection off for the force-push and ON again after; GitHub support asked for the purge; open PRs noted.
- [ ] The announcement says: delete every clone, re-clone once, never push from a pre-rewrite clone.

Joining the repositories
- [ ] Resources moves under a subdirectory of macos with history (`filter-repo --to-subdirectory-filter`, then
      `merge --allow-unrelated-histories`), and the subdirectory name is stated.
- [ ] The path fix-ups of step 5 are scheduled as one gated commit after the join (the `../MACOS_resources` and
      `../macos` assumptions in `cli_tests/`, `CMakeLists.txt`, `.gitattributes`, `mmacos_setup.m`, the pymacos
      CMake, `rx_fixture_path.m`, `tests/context.py`).
- [ ] The two `release-exclude.txt` lists are folded into one, each pattern prefixed by its directory.

Promotion and release
- [ ] `dev-candidate` -> `dev` is a fast-forward in both repos (or the divergence is named).
- [ ] `dev` -> `main` is merge, strip (`tools/release_strip.py --apply` or an equivalent that is idempotent), one commit.
- [ ] The strip is checked on the result: 5 `CLAUDE.md` on macos main (root, macos_f90, giza, slsqp, manual) + 4 in
      the resources subtree; zero `PLAN_* / CURRENT_SLICE / BRIEF_* / NOTE_* / DRAFT_* / RUNBOOK_*`; `git grep -il
      npsol` reads 0 in sources.
- [ ] The step-4 gate list runs on a CLEAN clone of the new `main` before the release is cut, on Linux and on the Mac.
- [ ] The release notes name the version (4.00) and the re-clone requirement.

Dry-run state on 2026-10-27 (CC): macos strip 494 of 7802 files, 5 CLAUDE.md kept, no process file kept; resources
strip 39 of 6318 after adding the `BRIEF_*.md / NOTE_*.md / PLAN_*.md / STATUS_*.md / SESSION_STATE_*.md /
CURRENT_SLICE.md` patterns (the named-file list had missed 17 dyson5 briefs, 2 notes and `design/PLAN_CONFIGURATIONS.md`);
the self-containment sweep finds no `addpath` reaching outside the repos, and the two `matlab_dir` options
(`segment_grid_basis`, `gs_zernike_segment_basis`) are deprecated-and-ignored defaults, not dependencies.

## What this runbook does NOT do
Nothing here is reversible after step 2's force-push except from the step-1 bundle; nothing here is run by CC.
