# Runbook: the October 2026 release moves (for Andy; drafted by CC for Dave, 2026-10-21)

Six steps, in this order, nothing run before the Friday meeting approves it.  Every count below was re-measured
on 2026-10-21 (`git log --all --diff-filter=A --name-only`): macos history holds 423 NPSOL paths and 1564 pgplot
paths; MACOS_resources history holds 21 OPTIIX/IRIS prescription paths (0 at HEAD); `.git` is 550 MB (macos) and
2.2 GB (resources).  Who runs what: steps 1, 2 and 6 need admin on nasa-jpl (Andy); steps 3-5 are ordinary pushes
Dave can make after Andy's go.  The branch rules are `BRANCHES.md`; the strip lists are `release-exclude.txt` in
each repo (amended 2026-10-21, reviewed in `RELEASE_EXCLUSIONS_REVIEW.md`).

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
three Claude worktrees (`~/dev/macos`, `~/dev/macos_cli_base`, `~/dev/macos_to_hex`) included.

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

## What this runbook does NOT do
Nothing here is reversible after step 2's force-push except from the step-1 bundle; nothing here is run by CC.
