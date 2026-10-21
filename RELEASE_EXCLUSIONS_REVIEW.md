# Release exclusions, amended 2026-10-21 (for review before the dev-candidate -> dev -> main moves)

**Principle (Dave 2026-10-21):** what an AI agent needs to WORK ON the code ships on `main`; what records HOW WE
WORKED stays on `dev`.  The strip lists (`release-exclude.txt` in both repos) are amended accordingly in this
commit; `BRANCHES.md` describes how they are applied.  Simulated on the trees at HEAD: macos keeps
6903 of 7383 tracked files (480 stripped); resources keeps 6266 of 6285 (19 stripped) -- the
2267 run artifacts are a separate decision for the room, not in these lists.

## Now SHIPS (moved off the strip lists)

| Repo | File | Why an agent (or a Mac user) needs it |
|---|---|---|
| macos | `CLAUDE.md` (root) | the tree-wide rules, build summary, branch model, the subsystem index; its post-compaction directive now says "if present" for the plans, slice and memory |
| macos | `macos_f90/CLAUDE.md` | the engine index: which `notes/` file holds which cluster of gotchas |
| macos | `macos_f90/giza/CLAUDE.md` | the plotting layer's gotchas (raster inversion, glyphs, device reopen) |
| macos | `macos_f90/slsqp/CLAUDE.md` | the constrained optimizer's dispatch and scaling rules |
| macos | `macos_f90/notes/*.md` (9 files; never excluded, stated here for completeness) | the engine record: surfaces, rx_io, calib, cli, pupil/FEX/STOP, OPD reference, polarization, api bindings, hygiene |
| macos | `docs/macos-manual/CLAUDE.md` | the manual's build rules: generated appendix, transcripts, figures, the product-name convention |
| macos | `MAC_PORT.md` | the Apple-Silicon build guide and the clean-clone validation record (was on the strip list although it is the Mac user's entry point) |
| resources | `mmacos/CLAUDE.md`, `pymacos/CLAUDE.md`, `GMI/CLAUDE.md`, `segmirmaker/CLAUDE.md` | the per-package gotchas and build/test wiring; the two binding files' directives now say "if present" for the plans and memory |
| resources | `optical_design/OPTICAL_DESIGN_AGENT_GUIDE.md`, `CORONAGRAPH_DESIGN_AGENT_GUIDE.md`, `AGENT_NOTES.md` | agent-facing by name: how to drive the design layer |
| resources | `optical_design/CORONAGRAPH_DESIGN_RULES.md`, `TELESCOPE_DESIGN_REFERENCE.md`, `SPECTROMETER_DESIGN_REFERENCE.md`, `coronagraph_layout.md` | the design rules and references the templates and challenges cite |
| resources | `optical_design/fixtures/*.md` (3) | the golden design fixtures the tests compare against |
| resources | `doc/mmacos_pymacos_parity.md` | which binding has which function |

Agent MEMORY files never ship (Dave's ruling 2026-10-19); Appendix D of the manual points at the CLAUDE.md set.

## Stays EXCLUDED (unchanged unless noted)

| Repo | Pattern(s) | What it is |
|---|---|---|
| macos | `PLAN_*.md` (NEW pattern: it now catches PLAN_CONSOLIDATION, PLAN_ATMOSPHERE, PLAN_ERROR_BUDGET, which the old list of three named plans missed), `PLAN.md`, `CURRENT_SLICE.md` | the program's sequencing and in-flight state |
| macos | `ENGINE_ISSUES.md`, `DEV_FILES.md`, `POLARIZATION_PHASE0_AUDIT.md`, `macos_f90/SAVE_KEYWORD_AUDIT.md`, `NOTES/` | working audits and the old notes directory |
| macos | `REVIEW_POL_*` (10 files) | the polarization review packets |
| macos | `BRIEF_*.md`, `DRAFT_email_*.md`, `REPORT_*.md`, `RUNBOOK_*.md`, `NOTES_*.md`, `NOTE_*.md`, `MR_*.md`, `MERGE_HANDOFF_*.md`, `demo_session/` | the lane-to-lane process files and the demo session |
| macos | `docs/macos-manual/FIGURE_RESCUE_LOG.md`, `RECONCILE_4_01.md`, `src/_dropped_legacy_index.md`, `audit_4.1beta.txt`; `docs/Archive/dev_optimization_surfsub/` | the manual's working files and the archived patch set |
| macos | `release-exclude.txt`, `RELEASE_EXCLUSIONS_REVIEW.md` (this file) | the strip machinery itself |
| resources | `DEV_FILES.md`, `release-exclude.txt`, `segmirmaker/.vscode/*` | working files and editor config |
| resources | `doc/STYLE_REPORTS.md` | the internal deck/report style guide |
| resources | `mmacos/design/PLAN_AFOCAL4.md`, `PLAN_TMA_E2E2.md`, `mmacos/PLAN_EXAMPLES_REORG.md`, `templates/.../PLAN_VSG2_MODELS.md`, `PLAN_IFO_PUPIL_RELAY.md`, `templates/90_polarization/bench_ifo_pol/CURRENT_SLICE.md` | plans and a slice inside the resources tree |
| resources | `challenges/afocal4/STATUS_S4B.md`, `STATUS_S4C.md`, `challenges/rodgers1/OPTFEX_REDO_LIST.md`, `SESSION_STATE_2026-07-30.md`, `templates/.../CTB_PROP_STATUS.md`, `VSG2_review.md` | status and session records |
| resources | `mmacos/BRIEF_ep_dome_ccmac_handback.md`, `templates/30_instruments/bench_ctb/BRIEF_dst_defects.md` | briefs |

Test fixtures (`ZGD_test_files/`, `cli_tests/`, the pymacos and mmacos `tests/`) stay on both branches as before.

## Dangling references checked
Every file that ships was grepped for the names of files that do not.  The three post-compaction directives
(macos root, mmacos, pymacos) now read "if present"; the manual's `RECONCILE_4_01.md` line and the polarization
note's audit pointer are marked "dev branch only".  The remaining mentions are citations of the form
"PLAN_CONSOLIDATION 1.2 (tt)" or "BRIEF_to_bindings A.2" in the notes -- provenance of a rule, declared as such
in the root `CLAUDE.md`, not files an agent is told to open.
