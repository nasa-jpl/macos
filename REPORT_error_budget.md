# REPORT — optical error-budget tool (`macos.budget` + `run_budget`)

**For CCL — review request.**  CCMac, 2026-10-15 (Dave travelling, out of
direct contact).  Resources-side, engine-independent; **pushed to
`origin/dev-candidate` (MACOS_resources)** as `3e85c0db..a0961dcb` (5
commits).  Nothing on the engine (`macos`) changed.

## What Dave asked

Plan + build a generic tool to turn the WFE (or SPOT) **sensitivity
matrices** (`dw_dx`, `dw_dz_zernike`, `dw_dsurf`, `dw_dgrid`) into **optical
error-budget coefficients**, following Redding 2014 (SPIE), *"Optical error
budgeting using linearized ray-trace models"*
(`macos_sandbox/WFE_budgets/Redding_LinRayTrace.pdf`): a properly-formed
RSS budget IS a covariance analysis,
`alpha = sqrt(trace(Pu·D_g·Sigma_g·D_g'·Pu')/n)`, roll-up
`WFE = sqrt(sum sigma_i^2 alpha_i^2)`; Monte-Carlo is the identical answer
when assumptions match.  The consuming spreadsheet is
`macos_sandbox/WFE_budgets/HabExB_WFE_budget.xlsx` (and the ATLAST one) —
a "WFE" sheet (σ inputs, Pre/Post-WFC + Metrology columns) over a
"Coefficients" sheet keyed `coef.zg0/zs0/u0/l0/w0/Xi0` for two control
sequences (iSeq=1 WFC, iSeq=2 metrology).

Dave's scoping decisions (via AskUserQuestion): **WFE/OPD first** (kernel
kept metric-agnostic so SPOT drops in later); **uncontrolled first**
(`Pu=I`) to validate the covariance≡MC core before layering controls;
**both** outputs (a `coef.*` struct/.mat/CSV that drops into the
HabEx/ATLAST sheets AND a fresh self-contained `.xlsx`); grouping
**auto-proposed then user-edited**, the group table authored from **Excel,
CSV, or a `.m`** interchangeably.  Full plan:
`~/.claude/plans/linked-sauteeing-oasis.md`.

## What landed (v1 — uncontrolled, WFE)

All in MACOS_resources.

- **Kernel** `mmacos/src/+macos/+budget/`: `assemble` (load a `_sens.mat`
  nested `ox/oz/og/os`, a flat `save_dw_flat`, or a live `dw_d*_multi`
  struct; unify per-column metadata; Gram `G=Jk'*Jk`), `cov_budget`,
  `mc_budget`, `project_control` (the `Pu` hook), `covariance_model`
  (per-group Σ in native units; normal/uniform/worstcase),
  `select_columns`, `spec_from_table`, `suggest_groups`, `read_spec`,
  `rollup`, `to_coef`, `write_csv`, `write_xlsx`, `zernike_decomp`.
- **Runner** `mmacos/design/runners/run_budget.m` — the pipeline's budget
  stage: assemble → group table → covariance + Monte-Carlo → report +
  self-contained `.xlsx` (WFE + Coefficients) + CSV + `coef.*` `.mat` +
  per-group figure.
- **Worked example** `mmacos/templates/50_sensitivities/e5hex1_budget/`
  (`run_budget_e5hex1.m`, `budget_spec_e5hex1.m`, `README.md`) on the
  public 7-segment hex deck.
- **Tests** `mmacos/tests/tBudgetKernel.m` (7/7) + `tRunBudget.m` (2/2),
  both added to `SUITE_FAST`.
- **Docs**: `mmacos/CLAUDE.md` "Error budgets" section (gotchas),
  `mmacos/design/runners/README.md` budget stage + runner row.

**Design heart:** the covariance/MMSE/Monte-Carlo kernel already existed in
`run_met.m` (`rmsw = sqrt(trace(P·G)/nw)`, `G=Dk'*Dk`; `post = X -
X·H'·((H·X·H'+R)\(H·X))`; the MC draws).  I **lifted that Gram-trick
kernel out, metric-agnostic, across all four Jacobian channels** rather
than reinventing it.  Per-group it is `trace(Sigma_g·G_gg)/nw` with `G_gg`
the group's diagonal Gram block; controls survive as `(Pu·Jk)'(Pu·Jk)` at
`Jk`'s size — the dense `Nw×Nw` is never formed (the same discipline as
[trace(D·P·D') OOM] — always `trace(P·G)`).

## Verification

- **`tBudgetKernel` 7/7** (pure algebra, no engine): cov ≡ MC per group &
  total + RSS identity + `WFE = σ·α` (the paper's Fig. 6 gate); control
  projection reduces WFE and cov ≡ MC under control; `cov_budget`
  reproduces `run_met`'s `sqrt(trace(X·G)/nw)` to 1e-12 (kernel
  refactor-safety); uniform/worstcase cov ≡ MC; unit reconciliation
  (nrad/nm via `cbm`) + orphan/overlap warnings; Zernike-mode recovery.
- **`tRunBudget` 2/2** (real `dw_dx_multi` harvest on `e5hex1.in`, model
  128): end-to-end budget **104.2 nm (cov) vs 104.6 nm (MC, 800 draws) =
  0.4%** — Fig. 6 agreement on real engine data; all artifacts written;
  both workbook sheets present.

## For CCL to review / gate

1. **Run the Linux gate.**  `./run_mmacos_tests.sh fast` (or at least
   `tBudgetKernel` + `tRunBudget`) on a Linux engine build.  `tBudgetKernel`
   is platform-independent algebra; `tRunBudget` depends on the engine
   trace of `e5hex1.in` — the cov≡MC agreement should hold identically.
2. **No A/B needed:** I touched **no existing source** — only new
   `+budget/` files, a new runner, two new test classes, two additive
   `SUITE_FAST` entries, and additive doc sections.  Regression risk is
   nil; still, confirm the full fast suite stays green on Linux.
3. **Build note (Mac):** the local `build_release_gfortran` was stale
   (missing newer `macos_api_mod` symbols — `opd_mask_get`, `modes_get`,
   `calib_stats_get`, `grid_dropped_get`); a clean `rm -rf
   build_release_gfortran && makems.sh release gfortran` was needed before
   the mmacos mex would build.  (`makegfortran` had left `BUILD_GMI=ON` in
   the CMakeCache, which FATALs on the Mac's `/Applications` MATLAB path.)
   Not a code issue — flagging in case Linux sees the same stale-lib
   symptom after today's `dev-candidate` engine pull.

## Deferred (next phases, in order)

1. **Controls phase** — the `Pu` WFC pseudo-inverse stage
   (`G=pinv(Ju)`, already prototyped in `project_control`) + the second
   **metrology/drift MMSE stage** (reuse `run_met`'s `post`/`K`), filling
   the `coef.*` Post + iSeq=2 slots and the HabEx two-stage columns.
2. **SPOT (`dcdx`) metric** — kernel already takes `(J, w0)`; feed the
   `compute_los` spot Jacobians + a SPOT adapter.
3. **pymacos twin**; **`mode_budget`** multi-field per-Zernike-mode rows
   (`zernike_decomp` utility is in); **byte-match `coef.*`** to the HabEx
   legacy keys before freezing `to_coef`/`write_xlsx`.

Awaiting Dave's word on whether to start the controls phase.  Memory:
`project_error_budget_tool`.
