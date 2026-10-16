# PLAN_ERROR_BUDGET — from sensitivity Jacobians to the error-budget workbook

Owner: Dave.  Lanes: CCMac (v1, 2026-10-14/15), TO (controls + workbook phases, from
2026-10-15), CC (gates, this plan).  Method: Redding 2014 (SPIE), "Optical error budgeting
using linearized ray-trace models" (`MACOS_sandbox/error_budgets/Redding_LinRayTrace.pdf`);
the consuming artifact: `MACOS_sandbox/error_budgets/HabExB_WFE_budget.xlsx` (the HabEx-B
WFE calculator — Fig. 5 of the paper is its ATLAST ancestor).  Everything is resources-side
(`MACOS_resources/mmacos`); the engine is unchanged by this plan.

## 0. The model, in the paper's terms (what every phase must agree with)

- Linear model `w = (dw/dx) x + (dw/dz) z + (dw/du) u + w_nom` (Eq. 13); the Jacobians are
  the `dw_d*` harvests (`run_sensitivities`), columns in native units (rad, m, BaseUnit).
- **WFC (iSeq = 1):** pseudo-inverse controller `G = pinv(dw/du)` (Eq. 20), residual
  `w1 = Pu (dw/dx) x0 + Pu (dw/dz) z0 + (dw/du) G dw_est + (dw/du) du` (Eq. 22) with the
  control projection `Pu = I − (dw/du) G` (Eq. 23).  Three post-control terms, each a
  coefficient: the UNCORRECTABLE error (fitting error, `α = sqrt(trace(Pu D Σ D' Pu')/n)`,
  Eq. 30), the SENSING error that leaks through the controller (`α = sqrt(trace(D_u G G' D_u')/n)`,
  Eq. 32), and the ACTUATOR error (`α = sqrt(trace(D_u D_u')/n)`, Eq. 31).
- **Metrology / drift (iSeq = 2):** the second control sequence — the drifts after the last
  WFC update, held by a metrology-driven estimator: `run_met.m`'s MMSE update
  `post = X − X H' (H X H' + R)⁻¹ H X` (X = the drift covariance in state space, H = the
  metrology Jacobian `dmet_dx`, R = gauge noise), then the same `Pu D post D' Pu'` roll-up.
  Pre-metrology = the drift σ through the uncontrolled coefficient; post-metrology = through
  the estimator + controller.
- **Covariance ≡ Monte-Carlo** (Eq. 27 vs Eqs. 22+24; Fig. 6) is the correctness gate of
  every phase; a phase that cannot show it on real engine data is not done.
- Roll-up `WFE = sqrt(Σ σ_i² α_i²)` (Eq. 33), with the design residual and any
  non-ray-trace terms (PSD figure, thermal/strain sheets) RSS'd in as rows the system
  engineer owns.

## 1. v1 — uncontrolled, WFE (CCMac, DONE, pushed resources 3e85c0d..a0961dc)

`+macos/+budget/` kernel (`assemble`, `cov_budget`, `mc_budget`, `project_control`,
`covariance_model`, `select_columns`, `spec_from_table`, `suggest_groups`, `read_spec`,
`rollup`, `to_coef`, `write_csv`, `write_xlsx`, `zernike_decomp`), runner `run_budget`,
example `templates/50_sensitivities/e5hex1_budget/`, gates `tBudgetKernel` 7 +
`tRunBudget` 2 (e5hex1: 104.2 nm cov vs 104.6 nm MC).  CCL's Linux gate 2026-10-15: green.
Report `REPORT_error_budget.md`.  Carried into phase 2 as-is: the group table schema
(`Name Elements DOFs Modes Param Channel Dist Sigma Units Min Max Cov Control`), the
unit-homogeneous-group rule, the Gram discipline (`trace(Σ_g G_gg)/nw`, never `Nw × Nw`),
and the `coef.*` shape with its NaN Post / iSeq=2 slots.

## 2. Controls phase — fill the Post and iSeq = 2 columns (TO, `BRIEF_to_error_budget_controls.md`)

2a. **WFC stage, iSeq = 1.**  Inputs: the `Control` column of the group table, the WFS
    noise σ and the actuator noise σ per control group.  **Dave's defaults (2026-10-15):
    control set = 6 DOF per PM segment + the SM's 6 DOF** (rigid-body actuators; no DM / FCA
    unless the deck carries them); **WFS σ ≈ 3 nm rms** (phase retrieval on image data, the
    Eq. 32 sensing row); **actuator σ ≈ 1 nm** per control DOF (the Eq. 31 rows; a rotation
    DOF's 1 nm is the surface motion at the segment edge — state the conversion).  These
    are the sheet's starting values, editable in the workbook's white cells.  Outputs per error group: `alpha` (pre, exists), `alpha_post` (Eq. 30 with
    `project_control`, exists as a prototype), plus two new ROWS the HabEx sheet carries under
    "WFC": the sensing-error row (Eq. 32) and one actuator-error row per control group
    (Eq. 31).  Gate: cov ≡ MC under control on e5hex1 (extend `tRunBudget`), the Eq. 31/32
    rows reproduced by a Monte-Carlo with the controller in the loop, and
    `project_control`'s `Pu` checked idempotent (`Pu Pu = Pu`) and null on the control span.
2b. **Metrology stage, iSeq = 2.**  Reuse `run_met.m`'s kernel verbatim (its `post`, its
    `K`, its base-vs-SI equivalence test) through a `+budget/project_met` that takes the
    drift group table (the "Drift Inputs" σ column), `dmet_dx` (`run_met`'s artifact or a
    live `dmet_dx` harvest) and the gauge noise R, and returns the post-metrology Gram per
    group.  Outputs: `alpha_met` (= `alpha` by construction when the drift acts on the same
    DOF), `alpha_met_post`; the WFC stage's `Pu` applied after the estimator (the HabEx
    "Post-metrology" column is metrology AND control).  Gate: cov ≡ MC with the MMSE
    estimator in the loop; `cov_budget` with the met stage reproduces `run_met`'s own
    post-metrology number on its e5 example to 1e-10 (the refactor-safety gate v1 already
    uses for the uncontrolled kernel).
2c. **`to_coef` filled, schema unchanged**: `.alpha_post`, `.alpha_met`, `.alpha_met_post`
    no longer NaN; new entries for the sensing / actuator rows (`.kind = 'wfs' | 'act'`).
    `rollup` gains the Post WFC / Pre-MET / Post-MET / Total columns with the HabEx
    arithmetic (Total = `sqrt(PostWFC² + PostMET²)` per row).

## 3. Workbook phase — generate the HabEx-form spreadsheet (TO, after 2; spec from the sandbox example)

The deliverable is a workbook the system engineer EDITS: σ inputs in white cells, every
other cell a live formula or a tool-written coefficient, so a changed σ re-rolls without
MATLAB.  Generated from `coef.*` + the group table by `+budget/write_workbook` (v1's
`write_xlsx` stays as the plain two-table export).  Layout, from `HabExB_WFE_budget.xlsx`:

- **Sheet `WFE`** (the calculator).  Row 1: date serial (A1), title (B1), section heads
  "Initial Alignment" (C1), "Metrology Control" (G1), "Total WFE" (K1).  Rows 2–3: column
  heads — `B Degree of Freedom | C Values D Units (Inputs) | E Pre WFC F Post WFC (RMS WFE nm)
  | G Values H Units (Drift Inputs) | I Pre-metrology J Post-metrology | K Post Metrology RMS
  WFE (nm)`.  One row per budget line, grouped by SECTION (A = section label on its first
  row: Optical Telescope Assembly, then one per optic, WFSC Unit, Metrology, WFC, WFE
  Summary).  Per row: `E = C*Coefficients!D`, `F = C*Coefficients!E`, `I = G*Coefficients!G`,
  `J = G*Coefficients!H`, `K = SQRT(SUMSQ(F, J))`.  Non-ray-trace rows (PSD figure, thermal,
  strain) point at their own sheets (`PSD!H2`, `StrTherm!B26`) — generated as pass-through
  rows with the formula hooks present and the sheets stubbed, the numbers the engineer's.
  **WFE Summary** rows: per-section `SQRT(SUMSQ(E/F/I/J/K rows))`, "Everything Else", the
  "Design WFE" input row (C68 → F68), "Payload WFE" `= SQRT(SUMSQ(K66:K68))`, "WFE
  Requirement" (input), "WFE Margin" `= 1 − K69/K70`.  Row numbers are whatever the group
  table produces; the formulas reference the generated rows (never hard-coded 66–71).
- **Sheet `Coefficients`.**  Same rows as `WFE`: `B name | C units | D Pre WFC | E Post WFC |
  F units | G Pre-MET | H Post-MET`, and to the right the raw `coef.*` blocks the tool
  writes: `J/K/L` = `coef.<key>0:` value (+ its MC check) for iSeq = 1, `N/O/P` = the same
  key for iSeq = 2; `D..H` cells are formulas pointing at those blocks (`D4 = K12`), as the
  HabEx sheet does, so the provenance of every coefficient is visible.  Key naming follows
  HabEx (`<dof><group>0` pre-control, `<dof><group>f` post-control: `zg0`, `zs0`, `u0`,
  `l0`, `w0`, `Xi0`, `X00`, `Q00`, `zsmid0`, `zshi0`, …, and the `f` twins); `to_coef`
  already sanitizes group names — the byte-match to the HabEx keys is item 5.
- **Sheet `WFE BD V4` (block diagram).**  The roll-up tree: Total / Static (post-WFC) /
  Drift (post-metrology) / Reserve, with `/400` wavelength rows and per-branch
  `SQRT(SUMSQ(WFE!…))` formulas over the section rows.  Generated as formulas over the
  section ranges; the box layout (merged cells, borders) is template furniture — copy it
  from a template workbook rather than draw it in code (see the writer note).
- **Sheets `PSD`, `StrTherm`, `Log`.**  Stubbed with their headers and the formula hooks
  the `WFE` rows expect; `Log` gets the generation line (date, sens artifact, SHAs).
- **Writer.**  MATLAB's `writecell` writes a string beginning with `=` as TEXT on Linux
  (verify on this box first — if it writes a formula, use it); otherwise: ship a
  TEMPLATE workbook (`+budget/templates/wfe_budget_template.xlsx`, the HabEx furniture with
  empty rows) and fill it with the pymacos venv's `openpyxl` (formula strings are formulas
  there) from a MATLAB-written JSON/CSV of rows — one tool, two steps, both in the repo.
  Gate: open the product headless (`soffice --headless --convert-to csv` recalculates on
  load) and check `Payload WFE` equals MATLAB's `coef.total` to 1e-9 relative; change one σ
  in the file, re-open, check the re-roll — a workbook whose formulas do not recompute is
  not the deliverable.
- **Example.**  `e5hex1_budget` grows the workbook; a second example on a HabEx-class
  deck (the e2e TMA or the design layer's `tma_longslit`) with PM segments as the control
  set, so the WFC rows are not vacuous (e5hex1's 7 segments suffice if the DM is absent —
  say so in its README).

## 4. After the workbook (order Dave's)

- **SPOT (`dcdx`) metric**: the kernel takes `(J, w0)`; a SPOT adapter over
  `compute_los`'s Jacobians, the same two stages, a `LOS` sheet.
- **pymacos twin** of `run_budget` (reads the same `_sens.mat`, writes the same workbook).
- **`mode_budget`**: per-Zernike-mode rows across fields (`zernike_decomp` is in).
- **Byte-match `coef.*`** to the HabEx keys on the ATLAST/HabEx sheets before freezing
  `to_coef` / `write_workbook`.

## 5. Rules carried from the rest of the program

One item in flight; its gate first, red before green; cov ≡ MC on real engine data is the
definition of done for every stage; no engine change; the shared mex is not relinked while
another MATLAB runs; commits local with SHA + branch to CC; Dave says push.
