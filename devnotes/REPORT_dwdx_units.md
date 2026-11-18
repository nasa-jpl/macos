# REPORT: dw_dx translation columns per BaseUnit ('trans_output') -- TO, 2026-10-06

Brief: `BRIEF_dwdx_units.md`.  Status: **implemented, gated, committed locally: resources dev-candidate 199342c (NOT pushed); fast suite 558/0; tRunCompare 5/5, tRunMet 4/4 (model 512).**  CC approved the helper approach after the audit (sec. 2).

## 0. Baseline pinned before any change

e5hex1.in (mm, cbm = 1e-3), model 128, `dw_dx(..., 'dofs', 0:5, 'compute_los', true)`:
10245 x 66 Jacobian + 66 x 2 dcdx, saved (scratchpad `pin_si_pre.mat`) for the
bit-for-bit `'si'` gate.  Sample: Elt 1 Tz column rms 6.947352e+02 (OPD-mm per
metre) -> 0.6947 OPD-mm per mm under the new default.

## 1. Where the scaling goes (ONE place)

All four channel kinds that carry translations (RigidBody, FocalPlane,
GroupedRigidBody, Source) take the poke in SI metres and use dof_idx 3..5
for translations.  `dw_dx` already hands `dwdx_for_current_source` an
`output_scale_fn(ch)` that multiplies BOTH the dwdx column and the dcdx row;
it is the identity today.  `trans_output='base'` makes it `cbm` for
dof_idx 3..5 (any channel kind).  `dw_dx_multi` forwards the knob through
`F.single` and records it.  `dw_dgrid` / `dw_dsurf` / `dw_dz_zernike` emit
NO rigid-body columns (they only mention dw_dx in comments) -> nothing there.
`dwdx_for_current_source` (the FD core) stays untouched.

## 2. Consumer audit (classified by READING the code, not by running it)

State units: what the consumer multiplies a translation column by.
"adapter-needed" = breaks by 1/cbm on an mm deck under the new default.

| file | use | state units | action | done |
|---|---|---|---|---|
| design/runners/run_sensitivities.m | harvests dw_dx_multi -> jac .mat (ox); prints per-segment column norms | n/a (producer) | none; ox carries `trans_output` (printed translation norms scale by cbm on mm decks -- the intended change) | - |
| design/runners/run_compare.m | `Bc*cbm` -> OPD-m per (rad\|m); pokes `poke_trans` in SI m via ch.apply; exports `dwdu` in SI | SI m | **adapter-needed**: convert ox translation cols to per-m when `ox.trans_output=='base'` (missing field = legacy .mat = 'si') | yes |
| design/runners/run_simulator.m | `Bc*cbm`; `Bc*X`, `Bc*Uh` (Tikhonov WFC, MET loop); X history in SI m/rad | SI m | **adapter-needed** (same helper) -- CONTROL LOOP | yes |
| design/runners/run_met.m | `D = ox.dwdxall` cols; merit trace(X*Gdt), X covariance in SI (rad\|m); feeds met_layout_opt | SI m | **adapter-needed** (same helper) -- the rot/trans balance of the MET merit moves by cbm^2 otherwise | yes |
| src/+macos/+design/met_layout_opt.m | takes D, header says "SI" | SI m (contract) | none if run_met converts; header already says SI | - |
| design/src/jacobian_check.m | `md = A0(:,q)*a`, `a = d_trans` in SI m | SI m | **adapter-needed**: scale the translation poke to the column's denominator (a/cbm when 'base') | yes |
| src/+macos/+design/System.m | `sensitivities` returns `out.rigid = dw_dx(...)`; vary/evaluate/optimize use FD traces, NOT the Jacobian | reports only | doc line: rigid.dwdx is per BaseUnit by default (`out.rigid.trans_output`) | yes |
| src/+macos/dw_dx_multi.m / private/dw_multi_core.m | forward + stack | n/a | forward `trans_output`, record it | yes |
| src/+macos/private/apply_opd_convention.m | orient/sign only | n/a | none | - |
| design/src/flag_zero_norm_channels.m | per-element block max-column RMS vs median of live blocks (relative) | n/a | none (relative; block max is a rotation column on every deck checked) | - |
| design/src/drop_channels.m | column removal | n/a | none | - |
| sensitivities/group_exhibit.m | group/member column-norm RATIO, same convention both sides | n/a | none in code; header units sentence -> "per BaseUnit (default) / per metre ('si')"; the driver comments in zoom_5x5 + run_dwdx_multi that say it "divides the group TRANSLATION columns by CBM" are STALE (it does not) -- fix the comments | yes |
| sensitivities/run_dwdx_multi*.m, save_dw_*.m, plot_dw_* | harvest via run_sensitivities, save, plot (autoscaled) | reports | none in code; records change (see sec. 3) | - |
| templates/50_sensitivities/zoom_5x5/run_dwdx_5zoom_5fov.m + README | report; header + README say "per SI METRE" | reports | doc update (README sec. ~104, ~162; header ~101, ~127) | yes |
| templates/50_sensitivities/run_dwdx_multi/run_dwdx_multi.m | report; header "per SI metre" | reports | doc update | yes |
| templates/50_sensitivities/e5hex1/run_dwdx.m, verifyall.m; e5hex2_refzern/verifyall.m | harvest + autoscaled tiles + dwdxall-vs-dwdzall alias check | reports | none | - |
| templates/10_telescopes/design_layer_api/example_*_from_rx.m | column rms print; `[rigid.dwdx, zern.dwdz]` join | reports | none (join is shape-only); printed rms moves x cbm on mm decks | - |
| templates/20_segmentation/e5_seg/e5_seg.m, e5_seg_metopt.m | OWN FD (`(wp-w0)*cbm/h`), not dw_dx | own SI | none (not a dw_dx consumer) | - |
| templates/80_end_to_end/e2e/s4,s6,s7 | run_sensitivities/compare/simulator; s7 `Bc` with comment "cbm=1, e2e" | SI m, **metre decks** | covered by the runner helper; factor 1 -> numbers unchanged | - |
| templates/80_end_to_end/e2e6m/s4,s5; e2e6m_r2/r0*,r3*,r4* | harvest; `A0 = per_field_dwdx` in the ridge/timeseries with states in m/rad | SI m, **metre decks** (every .in BaseUnits= m) | factor 1 -> numbers unchanged; s5/r4 basis_ would need the helper on a non-metre deck -- add it for robustness | yes |
| challenges/rodgers1/* | "rigid" = design-layer M2/M3 decenter/tilt readback, NOT dwdx | n/a | none (not a consumer) | - |
| src/+macos/+design/configs_from_table.m | DOF names only | n/a | none | - |
| tools/sens_core_ab/sens_core_ab.m | A/B of two harvests, same convention both sides | n/a | none | - |

Proposed adapter (one helper, used by the four adapter-needed consumers +
s5/r4): `macos.dwdx_trans_si(ox)` (or a design/src private) returning the
Jacobian blocks with translation columns divided by cbm when
`ox.trans_output=='base'`, identity when 'si' or when the field is absent
(every jac .mat written before this change).  This keeps old artifacts
loadable with no regen.

Alternative (smaller blast): `run_sensitivities` pins `'trans_output','si'`
so the runner pipeline is untouched -- rejected unless CC/Dave prefer it,
because run_sensitivities is also the harvest Luis-style templates
(zoom_5x5, run_dwdx_multi) go through, and they would then still emit per-metre.

## 3. Records whose translation numbers would change (x cbm) on a re-run (do NOT re-run; flag)

- mm decks only.  `templates/50_sensitivities/zoom_5x5/*_sens_report.txt`
  and `run_dwdx_multi/dwdx_multi_e5hex1_sens_report.txt` (translation column
  norms; the README's quoted group/segment RATIOS are unchanged -- same
  convention both sides).  `sensitivities/*.mat`, e5hex1 harvest .mats.
- e2e / e2e6m / e2e6m_r2: metre decks -> bit-unchanged.  rodgers: not a consumer.

## 4. What landed

- `dw_dx`: `'trans_output'` `'base'` (default) | `'si'`; the scale lives in its
  `output_scale_fn` (x cbm on dof_idx >= 3, all channel kinds) -- applied to
  the dwdx column AND the dcdx row.  `out.trans_output` recorded.  Header:
  the 08-25 paragraph's translation sentence replaced, history kept.
- `dw_dx_multi`: forwards + records `trans_output`.
- NEW `macos.dwdx_trans_per_metre(ox)`: identity when `trans_output` is absent
  (every pre-change harvest) or 'si'; on 'base' divides the translation
  columns of dwdx / dwdxall / per_field_dwdx and the dcdx / dcdx_per_field
  rows by cbm, sets 'si'; idempotent.
- ONE line at the ox load in: run_compare, run_simulator, run_met,
  jacobian_check, e2e6m/s5_timeseries, e2e6m_r2/r4_timeseries, e2e/s7_simulate.
- Docs: dw_dx + dw_dx_multi headers, System.sensitivities doc, zoom_5x5 README
  (units paragraph; the PM table is marked as a pre-10-06 per-metre harvest),
  zoom_5x5 + run_dwdx_multi drivers, sensitivities/README + two drivers,
  GroupedRigidBodyChannel, mmacos/CLAUDE.md unit-conventions bullet.
- `group_exhibit`: its EMITTED report units line was hard-coded "per SI METRE";
  it now says per BaseUnit / per SI METRE from `out.trans_output` (absent =
  per metre).  The two STALE driver comments ("the helper divides the group
  TRANSLATION columns by CBM") are fixed -- it never did.

## 5. Gates (all on mm decks, factor 1e-3)

tDwDx (e5hex1):
- (a) `test_trans_output_default_is_base_x_cbm`: base = si x cbm, RelTol 1e-12
  (measured 0.0 max rel), dcdx rows too; (d) rotation columns + rows bit-equal.
- (c) `test_trans_output_si_reproduces_the_pre_change_numbers`: measured
  BIT-FOR-BIT before committing the gate -- the full 66-column + dcdx harvest
  and the 12-column subset, HEAD-622ee52 dw_dx vs new 'si': `isequal` true.
  Committed pin = column rms + dcdx at 1e-12 (robust to an ulp engine rebuild).
- (b) `test_segment_piston_is_two_cos_aoi_per_baseunit`: Seg2 Tz under the
  chief reference: 1/7 of the rays at |c| in [2cos8deg, 2] = [1.9805, 2]
  BaseUnits per BaseUnit (measured 1.9864-1.9958), the rest < 1e-6.
- `test_trans_per_metre_helper`: legacy (no field) bit-identity, base -> x 1/cbm
  at 1e-12, idempotent, single + multi-field.
Consumers:
- tJacobianCheck `test_closes_identically_on_base_and_legacy_si_harvests`:
  n_mod / n_eng equal at 1e-12, translation rel < 0.05.
- tRunCompare `test_base_and_legacy_si_jac_compare_identically` (e5mono mm):
  run_compare w_rms_t / w_rel / dwdu and run_simulator u / rms_wfe_unc /
  rms_wfe_corr identical for a 'base' harvest and its legacy per-metre twin;
  w_rel < 0.05; the loop bites (corrected < 0.2 x uncorrected).
- tRunMet `test_base_and_legacy_si_jac_give_the_same_merits`: every merit +
  exported dwdx equal at 1e-10 for a 'base' twin of the synthetic jac.

MUST-FAIL, measured against the pre-change code / a missing adapter:
- old dw_dx default on gate (b): |c| = 1986.4-1995.8 vs band [1.98, 2] -> FAIL;
  gate (a): old dw_dx rejects 'trans_output' -> FAIL.
- jacobian_check handed per-BaseUnit columns WITHOUT the adapter (field
  stripped): rel = [1.8e-5 1.9e-5 NaN 0.999 0.999 0.999] -- the 1 - cbm
  signature on Tx/Ty/Tz; with the adapter [.. 9.5e-4 9.6e-4 1.3e-5].

## 6. One-liner for CC's memory (reference_rmswfe_units style)

`dw_dx` / `dw_dx_multi` translation columns + dcdx rows are OPD-BaseUnits per
BaseUnit by default since 2026-10-06 (`'trans_output','base'`; `'si'` = per
SI metre, the 08-25..10-06 default, = base / cbm); rotations per rad;
numerator BaseUnits; SI consumers call `macos.dwdx_trans_per_metre(ox)`
(identity on pre-change harvests).

Note for Luis (DRAFT_email_luis_spot2.md): the default CHANGED -- dw_dx columns
are now per BaseUnit (per mm on his deck) out of the box; `'trans_output','si'`
returns the previous per-metre columns.
