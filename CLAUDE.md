# MACOS Source Tree

NASA/JPL optical ray tracing code. Legacy Fortran, some files date to the 1980s.
Fixed-form source: .F files use the C preprocessor, .f files do not.

> **Post-compaction / post-upgrade — re-read the docs first.** After a
> context compaction or a tooling upgrade, before resuming build or
> engine work, re-read: **this file**, `PLAN.md`, `PLAN_DESIGN_LAYER.md`,
> `CURRENT_SLICE.md` (in-flight state), and the agent `MEMORY.md`
> (build/test workflow entries).
>
> **THEN re-read the nested `CLAUDE.md` for whatever subsystem you are
> resuming** — see the index below. Nested CLAUDE.md files load on demand
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
| engine core: `macos_api_mod.F90`, IO/parse (`msmacosio.inc`/`iosub.inc`/`macosio.F`), `elt_mod.F`, `macos_cmd_loop.inc`, `funcsub.F`, `surfsub.F`, `mathsub.F` | `macos_f90/CLAUDE.md` | §0 closed-hygiene bugs, CALIB wrappers, FreeForm=14 + user input, per-ray status, PERTURB, SXP, ZernType parse + apply-dispatch, deferred PolyObsVec, MOD crash guard, prescription validator, CLI sub-prompts, LOG cleanup, name-collision traps, macos_api_mod wrapper layer + buffer-wrapper template |
| pymacos wrappers / regression | `../MACOS_resources/pymacos/CLAUDE.md` | (existing) intensity/complex_field/dx_at/apodize, PROPER-compare; **cross-ref**, do not duplicate here |
| coronagraph test prescriptions | `../MACOS_resources/optical_design/CORONAGRAPH_DESIGN_RULES.md` + the `Element= Obscuring` mask rule | **cross-ref** |

> **CC:** If a gotcha spans two subtrees (e.g. a wrapper that touches both
> `macos_api_mod` and pymacos), put the engine-side detail in the nested
> file nearest the Fortran and leave a one-line pointer in the other.

---

## Build
CMake-based. All scripts live in `macos/` and accept `[debug|release] [gfortran]`
in any order (defaults: release, ifx). Each combination gets its own build directory.

PGPLOT was dropped on release-candidate: Giza is the only PGPLOT-API
provider now, and the legacy `[giza|pgplot]` script argument is gone.
NPSOL was removed on `main` (PR #44); the bound-constrained path is
now provided by SLSQP on `opt-dev` (default) and `main`'s unconstrained
LM (`nls_optim_dvr` in `design_optim_mod`) remains available for
non-constrained problems.

| Script | Targets |
|---|---|
| `source ./makeall.sh` | macos + smacos + smacos_dvr + GMI (all four) |
| `source ./makems.sh` | macos + libsmacos.a |
| `source ./makesd.sh` | smacos_dvr (builds macos/smacos too if needed) |
| `source ./makegmi.sh` | GMI.mexa64 (requires macos+smacos built first; defaults to gfortran — see GMI build choice section) |
| `source ./makegfortran.sh` | macos + smacos + smacos_dvr + GMI via gfortran (all four) |

All build scripts accept `npsol` to enable `-DUSE_NPSOL=ON` (this
branch only) — adds the NPSOL-backed constrained-optimization path
and links `libnpsol.a` + `liblapacklib.a` + `libblaslib.a`.  Default
OFF.

Build directory naming: `build_{release|debug}[_gfortran][_npsol]`

- CMakeLists.txt: top-level `CMakeLists.txt` (macos_f90 sources).
  (The legacy `macos_f90/npsol/CMakeLists.txt` is gone with rem_npsol.)
- CMakePresets.json: debug and release presets for VS Code CMake Tools integration.
- C compiler must be gcc (not icx) — legacy C files use implicit function declarations.
- smacos_dvr re-compiles macos_mod.F with -DCMACOS (smacos_lib's copy lacks CMACOS-only symbols like ifPGColor).
- **smacos_dvr MUST have its own `Fortran_MODULE_DIRECTORY` (`mod_smacos_dvr`).** It and the macos exe both compile macos_mod.F (+macos_vars_mod/macosio/pgplotsub) with -DCMACOS; emitting both into the same `mod_macos/` makes two **Ninja** rules generate `mod_macos/macos_mod.mod` → the Fortran dyndep step aborts ("multiple rules generate …"). Unix Makefiles tolerates the duplicate — why it stayed latent (`makems.sh` leaves BUILD_SMACOS_DVR OFF; a build dir that fell back to Makefiles never hits it). Fixed for issue #56 (sls-dev `84f0b9c`, opt-dev `cf01617`). The 4 legacy `makefile_*.sh` (pgplot/npsol/intel-14, dead) are deleted — **cmake is the one build path; ninja is OPTIONAL** (cmake falls back to plain `make`).
- GMI Makefile is **FC-aware** (`.fc_stamp`, gitignored): GMI.o/GMIG.o live in a flat dir shared by makegmi (gfortran default) + makeall (ifx); switching compilers used to reuse stale .o → `undefined reference to _gfortran_*`. The stamp forces an object rebuild when FC changes (sls-dev `5c729f0`/opt-dev `c2f3c19`).
- GMI.mexa64 has two build paths: (1) the standalone `MACOS_resources/GMI/Makefile` — more robust across MATLAB versions — which `makeall.sh`/`makegmi.sh` invoke with `MACOS_BUILD_DIR` pointing at the cmake build tree; and (2) a cross-platform cmake path (`add_subdirectory(MACOS_resources/GMI/CMakeLists.txt)`, Luis 2026-06) gated on `-DBUILD_GMI=ON`, used by `makegfortran.sh`, Windows, and opt-in Linux. **`makeall.sh` (both branches) keeps `-DBUILD_GMI=OFF`** so the Makefile is its sole GMI builder (no double build) — opt-dev's makeall kept `-DBUILD_GMI=ON` until `796dc51`, a double build whose cmake leg failed on a box lacking the hard-coded MATLAB (Scott's `GMIG.F fintrf.h` error). New `-DBUILD_MMACOS=ON` builds the mmacos mex the same way (`MMACOS_FC` selects its compiler). **Path gotcha:** the cmake `GMI_SRC_DIR`/`MMACOS_SRC_DIR` defaults must be `${CMAKE_CURRENT_SOURCE_DIR}/../MACOS_resources/{GMI,mmacos}`, NOT `../../../{GMI,mmacos}` — the latter is correct only for pymacos's `src/macos/CMakeLists.txt`; in the top-level macos `CMakeLists.txt` it resolves to `/home/GMI` and FATAL-breaks `makeall`/`makegfortran` configure on Linux. **MATLAB hard-code gotcha (2026-07-01):** the top-level `CMakeLists.txt` hard-coded `MATLAB_ROOT=/usr/local/MATLAB/R2025b`, so any cmake GMI/mmacos build (`-DBUILD_GMI=ON`; makegfortran/Windows) failed `GMIG.F`'s `#include <fintrf.h>` on a box without that exact release; now it auto-detects the latest installed MATLAB (mirrors the Makefile's `sort -V | tail -1`; override `-DMATLAB_ROOT=`), sls-dev `2f4948e`/opt-dev `796dc51`.

## Intel RPATH (self-contained ifx binary + smacos_dvr + GMI mex)
- Top-level CMakeLists.txt embeds RPATH for ifx builds so `/opt/intel/oneapi/...` libs
  are found without sourcing setvars.sh. Applies to macos executable, smacos_dvr, and
  GMI.mexa64: `BUILD_RPATH/INSTALL_RPATH` = `${INTEL_LIB_DIR};/opt/intel/oneapi/mkl/latest/lib`,
  `-Wl,--disable-new-dtags` forces DT_RPATH (transitive) over DT_RUNPATH (direct-only).
  libimf → libintlc transitive dep requires DT_RPATH.
- For gfortran builds: `CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES` (auto-set by CMake) used
  as RPATH so libgfortran/libquadmath are found without a module file.

## gfortran portability + GMI build choice
- `source ./makegfortran.sh release` now builds macos + smacos + smacos_dvr
  + GMI clean. The build was broken on release-candidate after the FreeForm /
  ZernType / EltType-table work landed; restored 2026-05-24 by fixing a batch
  of ifx-only constructs the gfortran compiler rejects:
  - `CHARACTER,DIMENSION(:),ALLOCATABLE :: EltName*N` — VAX `name*length`
    syntax. Use `CHARACTER(LEN=N), DIMENSION(:), ALLOCATABLE :: EltName`.
  - `(/'foo','bar quux'/)` char array constructors with mixed-length literals.
    Wrap with `[CHARACTER(LEN=Nmax) :: ...]` (Fortran 2003).
  - 2D INTEGER PARAMETER built from `[[a,b],[c,d],...]` — ifx flattens, gfortran
    sees rank-1. Use `RESHAPE([flat...], [rows,cols], ORDER=[2,1])` to preserve
    visual row-major layout.
  - `DO CONCURRENT` (forbidden by project convention anyway) — convert to
    sequential `DO` (with `IF` guard for the mask form).
  - Preprocessor directives indented past column 1 — `#define`, `#if`, `#else`,
    `#endif` must start in column 1 for gfortran's cpp. ifx's fpp is lenient.
    Bulk-fix with a `re.sub` to strip leading whitespace from these lines.
  - `Function FOO` with no `()` for a parameterless function — add empty parens.
  - `IF (LOGICAL == .TRUE.)` — use `.eqv.`.
  - `derived.component` access — use `%`.
  - REAL → LOGICAL implicit coercion (`logical_var = real_val`) — wrap with
    `(INT(real_val).EQ.1)`.
  - Slice assignment to an unallocated module-scope allocatable
    (`arr(:) = 0` when `allocated(arr)` is .FALSE.) — ifx silently no-ops,
    gfortran SIGSEGVs. Guard with `if (allocated(arr)) arr(:) = 0`. Caught in
    `src_mod_init_vars` (ds1/ds2 are allocated lazily by sourcsub).
- Build choice for GMI (`makegmi.sh` defaults to gfortran):
  - **gfortran (default):** `source ./makegmi.sh` (or `makegmi.sh release`).
    Mex lands at `~/dev/MACOS_resources/GMI/GMI.mexa64`. Requires
    `build_release_gfortran/` populated by `makegfortran.sh release` first.
  - **ifx (opt-in):** `source ./makegmi.sh ifx`. Requires `build_release/`
    populated by `makems.sh release` first.
- Both compilers run the GMI regression suite green (6/6) with bit-identical
  numeric results AND exit MATLAB cleanly (exit 0). They use the SAME source
  tree; the standalone `MACOS_resources/GMI/Makefile` and the top-level
  `CMakeLists.txt` have per-compiler conditionals for the compile/link flags
  (ifx: `-fpp -132 -gen-interfaces -fp-model strict -xHOST -shared-intel
  -reentrancy=none`;
  gfortran: `-cpp -ffixed-line-length-132 -march=native
  -fallow-argument-mismatch -std=legacy`, no Intel runtime).
- `-reentrancy=none` on the ifx link line is **load-bearing for ifx**: it
  switches the Intel Fortran runtime from the default multi-threaded variant
  (`libifcoremt.so.5`) to the single-threaded one (`libifcore.so.5`).
  `libifcoremt` keeps worker threads parked in the host process across mex
  calls; once MATLAB's `clear mex` unloads the DSO that owned their wake-up
  callback, those threads die at process exit by jumping to a now-unmapped
  function pointer. Symptom: `SIGSEGV` with `RIP=...e2c0` inside a
  `clone3/start_thread` stack — pure thread-spawn-to-freed-memory, no module
  of our own on the stack. GMI doesn't use the Intel Fortran runtime's
  internal threading anyway, so `-reentrancy=none` is a strict win. Without
  this flag, the ifx mex still produces correct numerical results but
  SIGSEGVs at MATLAB process exit, after the regression summary has printed.
  Diagnosed via stack-frame symbolization
  (`libc+641168 = start_thread`, `libc+1219692 = clone3`).
- Makefile gotcha: gfortran's `ld` is stricter about lib-after-object
  ordering than ifx's bundled linker. The Makefile splits compile flags
  (`LDFLAGS`, before objects) from trailing `POST_LIBS` (`-l:libmx.so` etc.,
  after objects/archives) so `--no-undefined` resolves under both. Also
  links `libfitslib.a` when present — gfortran needs it for FITS-related
  symbols; ifx is happy without.
- Use gfortran when investigating "weird ifx-only" behavior — its stricter
  rejection of latent bugs (uninitialized allocatables, etc.) surfaces
  problems that ifx silently hides. The Zernike-apply zero-response bug
  (see ZernType section) was caught exactly this way.

## Debug builds
- `source ./makems.sh debug`    — macos + smacos, -O0 -check all (ifx)
- `source ./makeall.sh debug` — all four targets, debug
- CMake debug uses -check all,noarg_temp_created (suppresses harmless array temporary warnings).
- VS Code: no `.vscode/` config exists currently (the old launch.json/tasks.json
  pointing at the retired `build_debug_giza` tree is gone; that tree and the other
  legacy build dirs — giza/pgplot/joint/smoke variants — were deleted 2026-07-02).
  If recreated, point launch targets at `build_debug[_gfortran]/` and set
  `debug.allowBreakpointsEverywhere: true` for breakpoints in .F files
  (requires ms-vscode.cpptools).
- smacos_dvr: compiles smacos_dvr.F with -DCMACOS, links against smacos_lib.a
  plus macosio.o/pgplotsub.o/macos_vars_mod.o/macos_mod.o from MACOS_OBJS.

## Build notes
- Never leave a surfsub.f (lowercase .f) in macos_f90/ alongside surfsub.F — cmake
  may pick it up and try to build it as a standalone executable.
- Use `source ./makems.sh` for macos_f90-only changes (faster than makeall.sh).
  cmake recompiles only changed files; delete the build directory for a clean rebuild.
- macos.F line 59: `use propsub_mod128` was a stale name; fixed to `use propsub_mod`.
- GMI (Matlab mex): built via standalone MACOS_resources/GMI/Makefile (invoked by
  makeall.sh — its sole GMI builder — and makegmi.sh), OR via the cross-platform
  cmake add_subdirectory path (`-DBUILD_GMI=ON`; makegfortran.sh / Windows). MATLAB
  auto-detected under /usr/local/MATLAB/R*. ifx 2025.3 quirk: -C flag implies
  -fsanitize=memory; use -check all instead.
- jGridSrf mapping: tracesub.F, propsub.F, srtrace.F use nGridMat(iElt).GT.0
  (not SrfType checks) so all grid-using surfaces get the correct GridMat slot.

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
(`~/dev/macos` + `~/dev/MACOS_resources`) are kept branch-symmetric —
as of 2026-07-28 both carry the same set (`dev`, `main`, `opt-dev`,
`bench-builder`, `pol-core`, `pol-ifo`, `release-candidate`, …) and
neither has `sls-dev`.

## Public-release strategy (UPDATED — Dave 2026-07-22; lands Friday ~2026-07-24)
- **Single repo, history rewritten.**  `nasa-jpl/macos` goes **public
  again** with its **history REWRITTEN to erase all presence of
  non-public code — especially NPSOL, but also pgplot, etc.** (scrubbed
  from history, not just deleted at HEAD).  This SUPERSEDES the earlier
  two-repo snapshot model.
- **Branches:** `main` is updated to **sls-dev functionality**; a public
  **`dev`** branch is retained.  All developer-facing files (CLAUDE.md,
  PLAN*.md, `.claude/`, internal ZGD fixtures, `docs/Archive/`) **stay on
  `dev` but are stripped from `main`**.  **Both branches are public.**
- **Users re-clone.**  Andy wants everyone to delete their local clones
  and re-clone (the history rewrite makes old clones diverge).  Do NOT
  push from a stale pre-rewrite clone afterward — it would reintroduce
  scrubbed history.
- Pre-rewrite full-history safeguard: `~/macos-archive-YYYYMMDD/`
  (`git bundle --all` per repo + worktree snapshots; e.g. the 2026-07-22
  archive taken before this rewrite).  See the agent MEMORY branch-model
  entry for the current tips.

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

---

> **Migration note (delete after the split lands):** sections not shown
> above moved to nested files per the index. Move text verbatim — do not
> paraphrase engine gotchas. After moving, each nested file gets the same
> post-compaction directive header pointing back here. Target lengths:
> root ≤ ~180 lines; each nested file holds its cluster. Verify with the
> doc test: give CC three representative prompts per subsystem against
> old vs. split, confirm equal-or-better rule-following.
