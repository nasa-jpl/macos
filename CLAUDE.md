# MACOS Source Tree

NASA/JPL optical ray tracing code. Legacy Fortran, some files date to the 1980s.
Fixed-form source: .F files use the C preprocessor, .f files do not.

> **Post-compaction / post-upgrade — re-read the docs first.** After a
> context compaction or a tooling upgrade, before resuming build or
> engine work, re-read: **this file**, `PLAN.md`, `PLAN_CONSOLIDATION.md`
> (the active sequencing, 2026-10-08), `PLAN_DESIGN_LAYER.md`,
> `CURRENT_SLICE.md` (in-flight state), and the agent `MEMORY.md`
> (build/test workflow entries).
>
> **THEN re-read the nested `CLAUDE.md` for whatever subsystem you are
> resuming — and, for the engine core, the `macos_f90/notes/<topic>.md` it points to** — see the index below. Nested CLAUDE.md files load on demand
> when CC reads a file in their directory and are **NOT auto-re-injected
> after compaction**; if you resume a giza / optimizer / engine-core task
> without first opening a file there, that subsystem's gotchas are not in
> context. Re-read the matching nested file explicitly. This directive
> heads each working-folder CLAUDE.md.

> **Why this file is short.** Subsystem gotchas live next to their source
> (load-on-demand, lower per-session context cost, higher adherence than
> one long file). Only genuinely cross-cutting, every-session rules stay
> here. When you add a gotcha, put it in the nested file for its subtree,
> not here — unless it applies tree-wide.

---

## Subsystem cheatsheets — load on demand (and re-read after compaction)

| Working in… | Read | Covers |
|---|---|---|
| `macos_f90/giza/`, `pgplotsub.F`, any plot/display path | `macos_f90/giza/CLAUDE.md` | raster inversion + PGCTAB CI-0/1 clobber, non-raster black-on-white, labels/titles/wedge, annotation lines, supersampling/idle-window, GRAY CR/CG/CB 0-based bounds, GRAINI caching, PANEL_ENV viewport, smacos_dvr graphics init, draw_rays data-only |
| `macos_f90/slsqp/`, `design_slsqp_optim.F`, constrained optim | `macos_f90/slsqp/CLAUDE.md` | SLSQP-default / NPSOL-opt-in dispatch, variable pre-scaling underflow fix (do **not** bump `dtt`), A/B reference, deferred setbeam port |
| engine core: `macos_api_mod.F90`, IO/parse (`msmacosio.inc`/`iosub.inc`/`macosio.F`), `elt_mod.F`, `macos_cmd_loop.inc`, `funcsub.F`, `surfsub.F`, `mathsub.F` | `macos_f90/CLAUDE.md` (a 54-line INDEX + the every-session rules) → then the topic file it names in `macos_f90/notes/` (`hygiene`, `surfaces`, `pupil_fex_stop`, `opd_reference`, `rx_io`, `calib`, `cli`, `api_bindings`, `polarization`) | the engine record, split 2026-10-15: each notes file holds its cluster verbatim; **nothing in `notes/` is auto-injected — read the one you touch** |
| pymacos wrappers / regression | `../MACOS_resources/pymacos/CLAUDE.md` | (existing) intensity/complex_field/dx_at/apodize, PROPER-compare; **cross-ref**, do not duplicate here |
| coronagraph test prescriptions | `../MACOS_resources/optical_design/CORONAGRAPH_DESIGN_RULES.md` + the `Element= Obscuring` mask rule | **cross-ref** |

> **CC:** If a gotcha spans two subtrees (e.g. a wrapper that touches both
> `macos_api_mod` and pymacos), put the engine-side detail in the nested
> file nearest the Fortran and leave a one-line pointer in the other.

---

## Build (summary — the full notes are `BUILD.md`)
CMake-based; scripts in `macos/` take `[debug|release] [gfortran]` (defaults release, ifx);
each combination gets its own `build_{release|debug}[_gfortran]` directory.  Giza is the only
PGPLOT provider; NPSOL is gone from `main`; SLSQP is the constrained path.

| Script | Targets |
|---|---|
| `source ./makeall.sh` | macos + smacos + smacos_dvr + GMI |
| `source ./makems.sh` | macos + libsmacos.a (engine-only; the fast path) |
| `source ./makesd.sh` | smacos_dvr |
| `source ./makegmi.sh` | GMI.mexa64 (gfortran default; `ifx` opt-in) |
| `source ./makegfortran.sh` | all four via gfortran |

Rules that bite: C compiler = gcc; never leave a lowercase `surfsub.f` beside `surfsub.F`;
smacos_dvr has its own module dir; `-reentrancy=none` is load-bearing on the ifx GMI link;
`source ./makems.sh debug gfortran` is the only build that bounds-checks every line.  Intel
RPATH, the gfortran-portability list (the ifx-only constructs gfortran rejects), the GMI
build choice and the debug-build notes are in `BUILD.md`.

## Branches (summary — the model and the public-release strategy are `BRANCHES.md`)
`dev-candidate` is where this program's work lands (both repos); `dev` is the public
integration branch, `opt-dev` the bug-fix release target, `main` the public release surface;
`sls-dev` and `bench-builder` are RETIRED (do not recreate).  **Engine first, resources
second**: never merge a resources branch whose emitters need engine keywords that are not on
the matching engine branch — the failure is at prescription-load time, not build time.
Always `git ls-remote --heads` before assuming the live set.  Push only on Dave's word.

## Key files for current work
- macos_f90/elt_mod.F        : per-element data arrays and SrfType constants
- macos_f90/surfsub.F        : all surface intersection routines
- macos_f90/elemsub.F        : calls surface routines during ray trace (FindSrf, Reflector, Refractor)
- macos_f90/param_mod.F      : array dimension parameters (mElt, mGridMat, etc.)
- macos_f90/iosub.inc        : ChkDf2 defaults, PrtSingleEltInfo output (included by macosio.F)
- macos_f90/msmacosio.inc    : prescription file reader (included by macosio.F and smacosio.F)
- macos_f90/macosio.F        : interactive UI dialog for element entry
- macos_f90/tracesub.F       : ray trace loop (Reflector + Refractor + FindSrf call sites)
- macos_f90/propsub.F        : propagation (Reflector + Refractor + FindSrf call sites)
- macos_f90/srtrace.F        : single-ray trace (Reflector call sites)
- macos_f90/funcsub.F        : CPERTURB, CPRead, CPERTURB_GRP (perturbation routines)
- macos_f90/macos_ops.F      : CPERTURB_2 (perturbation, macos_ops_mod)
- macos_f90/lnk_pert.inc     : LnkEltCPERTURB (linked-element perturbation)
- macos_f90/macos_api_mod.F90: language-neutral SMACOS-call wrapper layer (`MODULE macos_api_mod`) — backbone for pymacos and mmacos; compiled into libsmacos.a

## Conventions (new code)
- IMPLICIT NONE throughout
- Sequential DO loops only — no DO CONCURRENT
- DBLE() not FLOAT() for real conversion
- WRITE(*,*) + STOP instead of PAUSE
- Fixed-form column layout: code in cols 7-72, continuations in col 6
- Compiler flag -132 is set; lines may extend to col 132, but keep to 72 for readability
- `!-->` / `!<--` markers bracket added code blocks — follow the same pattern
- BLOCK construct is fine for local variable scoping inside ELSE IF branches

## Whitespace / editing notes
- Source files use mixed tabs and spaces for indentation (legacy; do not reformat).
  Fixed-form is column-sensitive so bulk retabbing risks introducing subtle bugs and
  destroys git blame. Leave indentation as-is.
- When the Edit tool fails due to whitespace ambiguity (tab vs spaces mismatch),
  use a Python script with exact string replacement instead.
