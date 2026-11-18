# Developer-only files — strip from public `main`, keep on `dev`

Manifest of developer/agent process files in the **macos** repo that should
be **removed from the public `main` branch** but **preserved on `dev`** for the
2026 public-release split.  These hold agent instructions, planning/working
state, and internal audits — *not* external user documentation (the manual,
command reference, README, and build guide stay on `main`).

This file is itself dev-only — it appears in its own strip list below.
`release-exclude.txt` is the machine-readable mirror of the strip list
(one path per line), consumed by the strip step:
`git rm -r -q --pathspec-from-file=release-exclude.txt`.
Keep the two in sync; the .txt is also dev-only.

## Strip list (paths from repo root)

```
CLAUDE.md
CURRENT_SLICE.md
PLAN.md
PLAN_DESIGN_LAYER.md
PLAN_POLARIZATION.md
POLARIZATION_PHASE0_AUDIT.md
ENGINE_ISSUES.md
MAC_PORT.md
DEV_FILES.md
release-exclude.txt
NOTES
REVIEW_POL_2026-07-26.md
REVIEW_POL_2C_2026-07-27.md
REVIEW_POL_ELEMENTS_2026-07-27.md
REVIEW_POL_EXTERNAL_2026-07-28.md
REVIEW_POL_IFO_SLICE1_2026-07-27.md
REVIEW_POL_IFO_SLICE2_2026-07-27.md
REVIEW_POL_IFO_SLICE3_2026-07-28.md
REVIEW_POL_OVERCOAT_CHROMATIC_2026-07-28.md
REVIEW_POL_RADIOMETRIC_2026-07-28.md
REVIEW_POL_SP_SIGN_2026-07-27.md
macos_f90/CLAUDE.md
external/giza/CLAUDE.md
external/slsqp/CLAUDE.md
macos_f90/SAVE_KEYWORD_AUDIT.md
docs/macos-manual/CLAUDE.md
docs/macos-manual/FIGURE_RESCUE_LOG.md
docs/macos-manual/RECONCILE_4_01.md
docs/macos-manual/src/_dropped_legacy_index.md
docs/Archive/dev_optimization_surfsub
BRIEF_*.md
DRAFT_email_*.md
REPORT_*.md
RUNBOOK_*.md
NOTES_*.md
MERGE_HANDOFF_*.md
demo_session
NOTE_*.md
MR_*.md
docs/macos-manual/audit_4.1beta.txt
GMI/CLAUDE.md
mmacos/CLAUDE.md
pymacos/CLAUDE.md
segmirmaker/CLAUDE.md
segmirmaker/.vscode/launch.json
segmirmaker/.vscode/tasks.json
doc/mmacos_pymacos_parity.md
doc/STYLE_REPORTS.md
mmacos/design/PLAN_AFOCAL4.md
mmacos/design/PLAN_TMA_E2E2.md
mmacos/challenges/afocal4/STATUS_S4B.md
mmacos/challenges/afocal4/STATUS_S4C.md
mmacos/challenges/rodgers1/OPTFEX_REDO_LIST.md
mmacos/challenges/rodgers1/SESSION_STATE_2026-07-30.md
mmacos/templates/40_benches/vsg_wip/PLAN_VSG2_MODELS.md
mmacos/templates/40_benches/vsg_wip/VSG2_review.md
mmacos/templates/30_instruments/bench_ctb/CTB_PROP_STATUS.md
mmacos/templates/40_benches/bench_ifo_dm/l2_trade/PLAN_IFO_PUPIL_RELAY.md
mmacos/templates/90_polarization/bench_ifo_pol/CURRENT_SLICE.md
optical_design/CORONAGRAPH_DESIGN_AGENT_GUIDE.md
optical_design/OPTICAL_DESIGN_AGENT_GUIDE.md
optical_design/CORONAGRAPH_DESIGN_RULES.md
optical_design/TELESCOPE_DESIGN_REFERENCE.md
optical_design/SPECTROMETER_DESIGN_REFERENCE.md
optical_design/AGENT_NOTES.md
optical_design/coronagraph_layout.md
optical_design/fixtures/telescope_design_fixtures.md
optical_design/fixtures/tma_fixture.md
optical_design/fixtures/afocal_tma_fixture.md
mmacos/BRIEF_ep_dome_ccmac_handback.md
mmacos/PLAN_EXAMPLES_REORG.md
mmacos/templates/30_instruments/bench_ctb/BRIEF_dst_defects.md
```

## Gate: dry-run the strip before EVERY promotion (dev -> main)

The list is a snapshot of the tree it was audited against; new process
files arrive with every arc.  Before each `dev -> main` promotion run the
strip as a dry run and read the SURVIVORS at the root:

```
git rm -r -n --cached --ignore-unmatch --pathspec-from-file=release-exclude.txt \
  | sed "s/^rm '//; s/'$//" | sort > /tmp/stripped
git ls-files | sort | comm -23 - /tmp/stripped | awk -F/ 'NF==1'     # root survivors
```

Every root survivor must be a user file (README, HOW_TO_COMPILE, LICENSE,
CMakeLists, make*.sh, platform_requirements, giza_build_notes).  Anything
else is a new family: add it (as a glob if it is a family) and re-run.
Two facts about `--pathspec-from-file`, both measured 2026-09-29:
- it takes NO comment or blank lines -- `#` is a literal path (skipped by
  `--ignore-unmatch`) and an EMPTY line is `fatal: empty string is not a
  valid pathspec`, which aborts the whole strip.  The .txt is a bare list;
  the reasons live HERE.
- a root glob such as `BRIEF_*.md` is ROOT-ONLY (it does not reach
  `demo_session/BRIEF_x.md`); nested families need their directory
  (`demo_session` is listed whole for that reason).

Re-run the gate again after the MACOS_resources import lands in this repo
(its list merges into this one; its paths gain no prefix if the
directories land at top level).

2026-09-29 re-audit against `dev-candidate` (the branch being promoted to
`dev`), and why:
- `NOTE_*.md` (2: the CCL close-out and CCMac welcome-back notes),
  `MR_*.md` (the 2026-09-09 merge-request text) and
  `docs/macos-manual/audit_4.1beta.txt` (the manual's structural audit
  worksheet) -- found by the gate above as root/process survivors of the
  first 2026-09-29 pass.  `NOTE_` (singular) is a different family from
  `NOTES_*.md`.
- `BRIEF_*.md` (44), `RUNBOOK_*.md` (5), `REPORT_*.md` (7), `DRAFT_email_*.md`
  (3), `NOTES_*.md`, `MERGE_HANDOFF_*.md` -- the agent<->agent briefs, Mac
  runbooks, internal review/finding reports and email drafts of the
  Aug-Sep 2026 gauge, sensitivity and afocal4 arcs.  Same class as
  `REVIEW_POL_*`.  Listed as GLOBS: git pathspecs glob, so the strip step
  covers new members without another edit here.
- `demo_session/` (216 files) -- deck builds (`deck_*.pptx`, edit copies,
  baselines, sidecars), their figures and the pptx tooling.  Working decks
  for JPL-internal discussions, not user documentation.
- NOT listed, deliberately: `ZGD_test_files/` -- the engine's own test
  fixtures (tst_save_keys.in, tst_block_comment.in, the FreeForm decks);
  stripping them would break the gates on `main`.  The root CLAUDE.md's
  "internal ZGD fixtures" wording predates the proprietary-Rx purge
  (OPTIIX/IRIS are gone) and should be read as this manifest.
  `README.md` and `HOW_TO_COMPILE.md` are user docs and stay.

2026-08-06 additions and why (all arrive via pol-core → dev, PR #67):
- `NOTES` (whole directory) — the agent↔agent hand-off channel
  (TO/CCL/CCMac notes).
- `PLAN_POLARIZATION.md` — sprint plan (working state).
- `POLARIZATION_PHASE0_AUDIT.md` — internal audit (same class as
  SAVE_KEYWORD_AUDIT).
- `REVIEW_POL_*.md` (10 files) — internal review packets.

## Explicitly KEPT on `main` (external documentation — do NOT strip)

- `README.md`, `HOW_TO_COMPILE.md`
- `docs/macos-manual/src/*.md` — the user manual (00–09, 90, 91, 93)
- `docs/macos-manual/cmdref/{00_orientation,01_cli_syntax,02_smacos_syntax,03_mmacos_syntax,04_pymacos_syntax,10_engine_commands,20_bindings,30_higher_level}.md`
- `docs/macos-manual/README.md`, `docs/macos-manual/cmdref/README.md` — doc-build guides (ship with the manual)
- `docs/macos-manual/polval/*.md` (2026-08-06, via pol-core) — the
  polarization validation report chapters: external validation
  documentation, part of the manual tree.

## Not tracked in this repo (won't be on `main` *or* `dev` unless committed to `dev` first)

- `.claude/`, `macos_f90/.claude/` — agent config dirs (currently untracked)
- `macos_f90/.fortls` — Fortran language-server config (untracked)
- The agent memory (`MEMORY.md` + `memory/*.md`) lives OUTSIDE the repo at
  `~/.claude/projects/-home-dcr-dev-macos/memory/` — not in the tree at all.

## Judgment calls (default shown; flip if desired)

- `docs/macos-manual/src/_dropped_legacy_index.md` — listed as STRIP (a dropped
  legacy index artifact, not part of the live manual). Keep only if it still
  serves a purpose.
- `docs/macos-manual/polval/` — listed as KEEP (validation report = user
  evidence of engine correctness). Flip to STRIP if it reads as internal.

## Carryover from MACOS_resources

2026-09-29 re-audit against `dev-candidate` (dry run of the strip, then
the process-looking survivors), and why:
- `mmacos/BRIEF_ep_dome_ccmac_handback.md`, `bench_ctb/BRIEF_dst_defects.md`
  -- agent<->agent briefs (same class as the macos repo's `BRIEF_*.md`).
- `mmacos/PLAN_EXAMPLES_REORG.md` -- the 2026-08-18 reorg plan (working
  state, like the `design/PLAN_*.md` entries).
- KEPT, deliberately: `mmacos/design/PLAN_CONFIGURATIONS.md` (linked from
  the zoom_5x5 README -- a plan, but a user doc points at it; stripping
  it breaks that link).
- The `REPORT_*.md` files inside templates (`bench_ctb/REPORT_field_servo`,
  `pdi_dm96/REPORT_gauge_pdi`,
  `tg_psi_dm96_oap/REPORT_{bench_realism,gauge_ifo,oap,reflective}`) were
  kept through 2026-10-27 as the benches' records (each linked from 1-4
  READMEs).  **Dave 2026-10-27: they do not ship** -- `REPORT_*.md` is on
  the exclude list by pattern; the files stay tracked on dev-candidate and
  in the local clones, and the README links to them break on the release
  surface.  What becomes of them is to be rethought before the re-clone
  (the runbook's keepers step).
- Gate: the dry-run recipe in the macos repo's `DEV_FILES.md` ("Gate:
  dry-run the strip") applies here verbatim; the .txt is a bare list
  (no comments/blank lines -- an empty line aborts `git rm`).

2026-08-06 additions and why:
- `doc/STYLE_REPORTS.md` — internal report/deck style gate (agent process).
- `mmacos/design/PLAN_*.md` — sprint plans (working state).
- `afocal4/STATUS_S4B.md`, `STATUS_S4C.md` — session status blocks.
- `rodgers1/OPTFEX_REDO_LIST.md`, `SESSION_STATE_2026-07-30.md` — working
  lists / session state.
- `vsg/PLAN_VSG2_MODELS.md`, `VSG2_review.md` — plan + agent-voiced review
  transcript.
- `bench_ctb/CTB_PROP_STATUS.md` — work-in-progress hand-off record
  (SESSION blocks; the user-facing content lives in the bench_ctb README
  and Coro_propagation_summary, both KEPT).
- `bench_ifo_dm/l2_trade/PLAN_IFO_PUPIL_RELAY.md`,
  `bench_ifo_pol/CURRENT_SLICE.md` — plan / in-flight state.
- `optical_design/fixtures/afocal_tma_fixture.md` — matches its two
  stripped sibling fixture notes.

## `optical_design/` — KEEP the code + fixtures (FUNCTIONAL, do NOT strip)

The shipped MATLAB design layer (`mmacos/+macos/+design/Telescope.m`,
`seidel_seed.m`) and the test suite (`tDesignTelescope` via
`mmacos/tests/private/design_fixture_path.m`) are ported/validated against
these — stripping them breaks tests:

- `optical_design/fixtures/{telescope_design_fixtures,tma_fixture,jwst_anchor}.json`
- `optical_design/{seidel,coronagraph_layout,make_fixtures,make_tma_fixture}.py`  (provenance tooling — KEPT per Dave)
- `optical_design/coronagraph_layout.json`

Only the agent-voiced `.md` docs in `optical_design/` are stripped (above).
`optical_design/README.md` is a user-facing README documenting the kept
fixtures/tools and **stays on `main`**; the previous agent-voiced README was
preserved as `optical_design/AGENT_NOTES.md` (dev-only, in the strip list).

## Explicitly KEPT on `main` (external documentation — do NOT strip)

- `README.md`, `GMI/README.md`, `mmacos/README.md`, `pymacos/README.md`
- `mmacos/doc/README.md`, `mmacos/doc/elt_grid_add.md` (interface doc)
- `mmacos/design/README.md` + `mmacos/design/runners/README.md`
- `mmacos/templates/00_INDEX.md` + `mmacos/templates/*/*/README.md`
- `mmacos/challenges/README.md`
- `mmacos/sensitivities/README.md`
- `mmacos/tools/fex_sweep/README.md`
- `GMI/regression/README.md`, `mmacos/tests/README.md`, `pymacos/tests/proper_compare/README.md` (test-run guides — user-usable)

2026-08-06 additions (KEPT — study records and technical documentation):
- `mmacos/challenges/rodgers1/{PACKET,README}.md`,
  `mmacos/challenges/rodgers2/{PACKET,README}.md` — benchmark-study records,
  retractions in place (the study precedent).
- `mmacos/challenges/afocal4/{FORM_STUDY,RESULTS,README}.md` — study
  records, same family as the packets.
- `mmacos/templates/80_end_to_end/e2e2/E2E2_REPORT.md` + e2e2/relay_followon READMEs
  — worked-example report.
- `mmacos/templates/30_instruments/bench_ctb/{README,Coro_propagation_summary}.md`
  + `bench_ctb/dst2/README.md` — house rules + the diffraction recipe of
  record (technical documentation, not process).
- `mmacos/templates/40_benches/bench_ifo_dm/l2_trade/TRADE_NOTE.md` — trade-study
  record (result-first).
- `mmacos/tools/pol_{external_anchor,sp_sign_probe,validation_report}/README.md`
  — tool run guides.
- `*.fp.json` fingerprint sidecars anywhere — functional (large-file
  policy), never stripped.

## Not tracked here (won't be on `main` *or* `dev` unless committed to `dev`)

- Agent memory lives OUTSIDE both repos at
  `~/.claude/projects/-home-dcr-dev-macos/memory/` (`MEMORY.md` + `memory/*.md`).
