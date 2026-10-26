# MACOS branch model and release strategy

> Moved verbatim out of the root `CLAUDE.md` on 2026-10-15 (the root is the every-session
> rules + the index; this file is read when you merge, push, promote or prepare a release).  Add new branch rules here.

## Branch model (UPDATED 2026-07-28 — `sls-dev` retired, `dev` is integration)
- **`dev`** — **integration + public developer branch for new work**
  (the role `sls-dev` used to hold).  Multi-person WIP lands here; all
  developer-facing files (CLAUDE.md, PLAN*.md, `.claude/`, internal
  fixtures) live here and are stripped from `main`.  **Both `dev` and
  `main` are public** per the public-release strategy below.
- **`sls-dev`** — **RETIRED 2026-07-28, deleted in BOTH repos.**  Andy
  deleted it on `nasa-jpl/macos`; matched on `MACOS_resources` the same
  day (`dev` fully contained it — `dev` = old `sls-dev` + the
  expose-beam beam-API commits).  Do not recreate it; new work → `dev`.
- **`opt-dev`** — release target.  Accepts bug fixes only.  Promotes
  to `main` for the public release (with NPSOL source tree removed
  at promotion time).
- **`release-candidate`** — **frozen** at `19bfbf8` (the mZern slice-
  overrun cherry-pick).  No new commits.  Pre-existing references to
  it elsewhere (PLAN.md, scripts) should be retargeted to `opt-dev`
  over time.
- **`main`** — public release surface; promoted from `opt-dev`/`dev` at
  release time per the public-release strategy below.

Day-to-day: push new features to `dev`; cherry-pick bug fixes to
`opt-dev`; let `dev` accumulate until a promotion gate.  The two repos
(`~/dev/macos` + `~/dev/MACOS_resources`) are NOT branch-identical —
`nasa-jpl/macos` is slimmed (as of 2026-07-28 remote heads: `dev`,
`fixREADME`, `main`, `pol-core`) while `MACOS_resources` still carries
the wider working set (`dev`, `main`, `opt-dev`, `pol-core`, `pol-ifo`,
`release-candidate`, `develop`, `dr-dev`, `dr-dev2`, `ifo-l2`).  The
shared, deliberate facts are that **neither repo has `sls-dev` or
`bench-builder`** — both retired 2026-07-28.  Always
`git ls-remote --heads` to see the live set before assuming.

### Merge ordering: engine first, resources second (Dave 2026-07-28)
`MACOS_resources` may be updated **without review**, but as a best
practice **do not merge a resources branch to `dev` until the engine is
updated to match.**  Engine (`macos`) lands first or in the same beat;
resources follows.

**Why it matters:** a veneer promoted ahead of its engine support fails
at **prescription-load time, not at build time** — so the resources test
suite can be fully green while the paired `dev` engine cannot parse the
Rx keywords that veneer emits.  Nothing on the resources side catches it.

**How to apply.** Before merging a resources branch, ask what engine
features its new emitters/wrappers depend on and confirm those are
already on macos `dev`.  Purely additive, engine-independent work
(examples, plotting, docs, pure-MATLAB layers) is exempt.

Worked example — the 2026-07-28 promotions, which followed the rule:
- **PR #10** (`bench-builder` → resources `dev`) was safe: its only
  engine dependency, the BEAM API (`beam_set`/`beam_get`), was already
  on macos `dev` via **PR #64**.
- The pol emitters (`Bench.add_polarizer` / `add_waveplate`, which emit
  `PolAxis=` / `Retardance=` and need Phase-3 `TrPolarizerElt(15)` /
  `WavePlateElt(18)`) were correctly **left on `pol-ifo`**, because that
  engine support is on macos `pol-core`, not macos `dev`.
- Consequence for the in-flight arc: resources `pol-ifo` → `dev` is
  **gated on** macos `pol-core` → `dev`.  Promote the engine first; the
  two are sequenced, not independent.

## Public-release strategy (UPDATED — Dave 2026-07-22; lands Friday ~2026-07-24)
- **Single repo, history rewritten.**  `nasa-jpl/macos` goes **public
  again** with its **history REWRITTEN to erase all presence of
  non-public code — especially NPSOL, but also pgplot, etc.** (scrubbed
  from history, not just deleted at HEAD).  This SUPERSEDES the earlier
  two-repo snapshot model.
- **Branches:** `main` is updated to **sls-dev functionality**; a public
  **`dev`** branch is retained.  All developer-facing files (CLAUDE.md,
  PLAN*.md, `.claude/`, the BRIEF/REPORT/REVIEW/NOTES process files,
  `docs/Archive/dev_optimization_surfsub` -- the list is `DEV_FILES.md` +
  `release-exclude.txt`) **stay on `dev` but are stripped from `main`**.
  `ZGD_test_files/` is engine test fixtures and stays on BOTH.  **Both
  branches are public.**
- **Users re-clone.**  Andy wants everyone to delete their local clones
  and re-clone (the history rewrite makes old clones diverge).  Do NOT
  push from a stale pre-rewrite clone afterward — it would reintroduce
  scrubbed history.
- **Release gate — self-containment sweep (Dave 2026-09-03, from
  Luis's zernike_mode catch):** before any public cut, grep BOTH
  trees for `addpath`/external-path references reaching outside the
  repos, and run the mmacos suite on a clean checkout with no
  `~/matlab` on the path.  A dependency every JPLer happens to have
  is invisible until a non-JPLer runs it — the pattern:
  `segment_grid_basis`'s Noll path needed `~/matlab/zernike_mode.m`;
  fixed self-contained as `macos.noll_mode` (98d3320, A/B 1e-10,
  gate `tNollMode`).
- Pre-rewrite full-history safeguard: `~/macos-archive-YYYYMMDD/`
  (`git bundle --all` per repo + worktree snapshots; e.g. the 2026-07-22
  archive taken before this rewrite).  See the agent MEMORY branch-model
  entry for the current tips.

## Repo consolidation and the rules after it (group meeting 2026-10-26, Dave + Andy et al.)
- **Order of events:** (1) this program consolidates macos and MACOS_resources -- orphaned files
  removed, documentation to release standard (`PLAN_CONSOLIDATION.md` 3b); (2) Andy writes the
  execution steps on the GitHub wiki -- `dev-candidate` -> `dev`, `dev` -> `main` with the
  redactions in `release-exclude.txt`, history preserved as far as possible with NPSOL / OPTIIX /
  IRIS sources removed, the two repos' histories merged -- and we review them against
  `RUNBOOK_release_moves_2026-10.md`; (3) Andy executes after Dave announces 1 and 2 done;
  (4) we delete every clone, re-clone, and test the new repo; (5) Andy makes the release.
- **Rules in force from that release on:**
  - a. a **temporary branch per bug fix or new feature** (no more work directly on the branch
    of record; the 2026-09-07 interim policy ends with the release);
  - b. **human-readable commits** -- clear, concise, plain language, readable by someone without
    deep MACOS experience (subject line = what changed and why, in words; the SHAs, gate counts
    and plan citations go in the body or the records, not the subject);
  - c. **PRs reviewed, merged and their branch deleted within a week**;
  - d. otherwise the order of events as specified in the runbook.

