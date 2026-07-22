# Developer-only files — strip from public `main`, keep on `dev`

Manifest of developer/agent process files in the **macos** repo that should
be **removed from the public `main` branch** but **preserved on `dev`** for the
2026 public-release split.  These hold agent instructions, planning/working
state, and internal audits — *not* external user documentation (the manual,
command reference, README, and build guide stay on `main`).

This file is itself dev-only — it appears in its own strip list below.

## Strip list (paths from repo root)

```
CLAUDE.md
CURRENT_SLICE.md
PLAN.md
PLAN_DESIGN_LAYER.md
ENGINE_ISSUES.md
MAC_PORT.md
DEV_FILES.md
macos_f90/CLAUDE.md
macos_f90/giza/CLAUDE.md
macos_f90/slsqp/CLAUDE.md
macos_f90/SAVE_KEYWORD_AUDIT.md
docs/macos-manual/CLAUDE.md
docs/macos-manual/FIGURE_RESCUE_LOG.md
docs/macos-manual/RECONCILE_4_01.md
docs/macos-manual/src/_dropped_legacy_index.md
docs/Archive/dev_optimization_surfsub/README.md
```

## Explicitly KEPT on `main` (external documentation — do NOT strip)

- `README.md`, `HOW_TO_COMPILE.md`
- `docs/macos-manual/src/*.md` — the user manual (00–09, 90, 91, 93)
- `docs/macos-manual/cmdref/{00_orientation,01_cli_syntax,02_smacos_syntax,03_mmacos_syntax,04_pymacos_syntax,10_engine_commands,20_bindings,30_higher_level}.md`
- `docs/macos-manual/README.md`, `docs/macos-manual/cmdref/README.md` — doc-build guides (ship with the manual)

## Not tracked in this repo (won't be on `main` *or* `dev` unless committed to `dev` first)

- `.claude/`, `macos_f90/.claude/` — agent config dirs (currently untracked)
- `macos_f90/.fortls` — Fortran language-server config (untracked)
- The agent memory (`MEMORY.md` + `memory/*.md`) lives OUTSIDE the repo at
  `~/.claude/projects/-home-dcr-dev-macos/memory/` — not in the tree at all.

## Judgment calls (default shown; flip if desired)

- `docs/macos-manual/src/_dropped_legacy_index.md` — listed as STRIP (a dropped
  legacy index artifact, not part of the live manual). Keep only if it still
  serves a purpose.
