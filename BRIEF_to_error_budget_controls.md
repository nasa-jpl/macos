# BRIEF — the error-budget tool, controls phase then the workbook (TO)

CC for TO, 2026-10-15, on Dave's word ("brief TO on the error budget work and have them
add the controls to the coefficient generation; the next step is the spreadsheet").  You
start from a cleared session; this file and `macos/PLAN_ERROR_BUDGET.md` are the whole
hand-off.  Read, in order: `PLAN_ERROR_BUDGET.md` (§0 the model in the paper's equations,
§2 your phase, §3 the next), `REPORT_error_budget.md` (CCMac's v1 — what exists, verified
how), `MACOS_resources/mmacos/CLAUDE.md` "Error budgets" (the gotchas), the paper
`MACOS_sandbox/error_budgets/Redding_LinRayTrace.pdf` (10 pages; Eqs. 19–33 are the
phase), the HabEx workbook `MACOS_sandbox/error_budgets/HabExB_WFE_budget.xlsx` (the shape
you are generating toward; `openpyxl` is NOT in the system python — read it with
zipfile+xml or the pymacos venv), then the code: `mmacos/src/+macos/+budget/*.m`,
`mmacos/design/runners/run_budget.m`, `run_met.m` (the MMSE kernel you will reuse),
`mmacos/templates/50_sensitivities/e5hex1_budget/`, tests `tBudgetKernel`, `tRunBudget`,
`tRunMet`.  Both repos `dev-candidate`, pushed through macos `6c62e6c` / resources
`8a1434b`; state SHAs when you commit.  Resources-only work — the engine does not change.

## A. Controls phase (plan §2)

A.1 **WFC stage (iSeq = 1).**  `run_budget('control', ...)` with the group table's `Control`
    column selecting the actuator DOFs — **Dave's ruling: 6 DOF per PM segment + the SM's
    6 DOF**, no DM/FCA unless the deck has them; **WFS σ 3 nm rms** (phase retrieval on image
    data); **actuator σ 1 nm** per DOF (rotations: 1 nm of edge motion — write the conversion
    you use); per error group the post-control coefficient via
    `project_control` (Eq. 23/30 — exists, prototype), plus the two rows HabEx carries under
    "WFC": the sensing-error row (Eq. 32, `α = sqrt(trace(D_u G G' D_u')/n)` with the WFS
    noise σ) and one actuator-error row per control group (Eq. 31).  Gram discipline
    throughout (`(Pu Jk)'(Pu Jk)` at Jk's size; never `Nw × Nw`).  Gates: `Pu Pu = Pu` and
    `Pu Ju = 0`; cov ≡ MC under control on the e5hex1 harvest (a Monte-Carlo with the
    pseudo-inverse controller in the loop, Eq. 22, including the sensing and actuator
    noise draws) — extend `tRunBudget`; the uncontrolled numbers unchanged to the digit.
A.2 **Metrology stage (iSeq = 2).**  `+budget/project_met`: the drift group table (the
    "Drift Inputs" σ), `dmet_dx` (from `run_met`'s artifact or a live harvest), the gauge
    noise R → the MMSE posterior `post = X − X H'(H X H' + R)⁻¹ H X` (lift `run_met.m`'s
    kernel; do not rewrite it), then the WFC projection on top (post-metrology = metrology
    AND control).  Outputs `alpha_met`, `alpha_met_post`.  Gates: `cov_budget` with the met
    stage reproduces `run_met`'s post-metrology number on its own example to 1e-10; cov ≡
    MC with the estimator in the loop; a drift on a DOF the metrology does not see passes
    through at `alpha_met_post = alpha_post` (the must-fail leg: remove the gauge rows and
    the post-met column must rise).
A.3 **`to_coef` / `rollup` filled** (schema unchanged): the NaN slots become numbers, the
    sensing / actuator rows appear with `.kind`, `rollup` carries Post WFC / Pre-MET /
    Post-MET / Total (= `sqrt(PostWFC² + PostMET²)`) and the summary rows (per section,
    Everything Else, Design WFE input, Payload, Requirement, Margin).  `write_xlsx` (v1's
    two-table export) and `write_csv` follow.

## B. Workbook phase (plan §3) — after A

`+budget/write_workbook`: the HabEx-form calculator with LIVE formulas (σ in white cells,
`E = C*Coefficients!D` …, `K = SQRT(SUMSQ(F,J))`, the summary rows, the Coefficients sheet
with the raw `coef.<key>0 / <key>f` blocks the D..H cells point at, the block-diagram sheet's
`SQRT(SUMSQ(WFE!…))` tree, stubbed PSD / StrTherm / Log).  First measure whether MATLAB's
`writecell` writes `=`-strings as formulas on this box; if not, template workbook +
`openpyxl` in the pymacos venv, driven from MATLAB — the plan §3 "Writer" paragraph.  Gate:
headless LibreOffice recalculation gives `Payload WFE == coef.total` to 1e-9, and a σ edited
in the file re-rolls.  Example: `e5hex1_budget` grows the workbook (the 7 segments as the
control set; say in the README that the DM is absent).

## C. How we work

One item in flight; gate first, red on the pre-change code where a red is possible (a side
checkout of resources at `8a1434b` serves as the pre-change tree — no engine rebuild
needed); fast suite after each landing (`./run_mmacos_tests.sh fast`, or the budget classes
+ `tRunMet` between); commits local with SHA + branch to CC; CC gates; Dave says push.
Numbers only from output you have read; `ps -C MATLAB` before a model-1024 run (none of this
needs one; e5hex1 is model 128).  Do not touch CCMac's kernel signatures without a note in
the commit — they are the plan's contract; add, don't rename.

## D. Report

One message to CC per landed item (A.1, A.2, A.3, B): SHAs, the gate numbers (cov vs MC
per stage, the `run_met` reproduction), anything measured that changes the plan.  Things
that need Dave (the workbook's section layout for a given deck; anything beyond the
control set and σ's he has ruled) go to him through CC with the measurement attached.
