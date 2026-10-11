# CALIB LM loop: old vs new, every CALIB in the CALIB classes (PLAN_CONSOLIDATION 3.3, 2026-10-11)

**Old** = engine a8c08bf with the iteration/merit REPORTING added and nothing else (the old
schedule: checks at int(.21T/.5T/.71T), T > 4 only, exit `itr_done .and. itr > 0.7T`, `itr_done`
a static never initialised).  **New** = this commit: flat = `|d rms| <= OptTol` (absolute,
`rms = sqrt(chisq/rows)`) AND `d chisq <= 1e-6 chisq` on two consecutive accepted steps; floor =
ten consecutive rejected trials; cap = T.  Both gfortran mex, one MATLAB process per class
(`run_mmacos_tests.sh <class>`), in this order: tAsphCalib tBeamRows tCalib tCalibTargets
tCalibBlindFields tDesignTelescope tSpectrometerRx (tCalibBlindFields reaches no LM run;
tAsphCalib's failure case raises before the log).  Merit = the final LM chisq.  `d sqrt` =
`sqrt(chisq_new) - sqrt(chisq_old)`, an UPPER BOUND on the change of the rms merit for any row
count; OptTol = 1e-12 on every case (traced).

**Result: every case within OptTol of the old loop or better** (worst: fourth_mirror_psi_parity
/ add_fold, `d sqrt` 1.8e-14, +2.7e-6 relative chisq, ~1e-15 m of WFE).  Iterations 822 -> 615.
Stop reasons: 14 cap, 14 floor, 2 converged.

**'cap' where the old loop exited converged** (all T = 10, merit identical or better):
tAsphCalib/test_asphere_step..., tCalib x5 (runs_baseline, converges_after_perturbation,
returns_expected_schema, set_var_elt_via_name_list, set_var_elt_via_positional_mask), tCalibTargets
zmode legs x2.  Mechanism: each is a QUADRATIC solve already at the round-off floor by itr 3-4; no
accepted step is ever flat, and the ten-rejection floor would land at ~itr 14, after the cap.
The old loop's '8' was its 0.7T exit, not a convergence test it could fail.

| class | case | old itrs | new itrs | new stop | old chisq | new chisq | rel chisq | d sqrt | verdict |
|---|---|---|---|---|---|---|---|---|---|
| tAsphCalib | test_asphere_step_scales_with_the_coefficient_and_calib_converges | 8 | 10 | cap | 1.7333369499e-43 | 1.7333369499e-43 | +0.00e+00 | +0.00e+00 | same |
| tBeamRows | test_direction_rows_ride_on_a_spot_target | 15 | 13 | floor | 4.6331826745e-21 | 4.6331826745e-21 | +0.00e+00 | +0.00e+00 | same |
| tBeamRows | test_per_field_position_targets_score_every_field | 15 | 14 | floor | 6.2499829422e-32 | 6.2499829422e-32 | +0.00e+00 | +0.00e+00 | same |
| tBeamRows | test_centroid_option_scores_the_ray_centroid_not_the_chief | 15 | 10 | floor | 0.0000000000e+00 | 0.0000000000e+00 | +0.00e+00 | +0.00e+00 | same |
| tBeamRows | test_centroid_option_scores_the_ray_centroid_not_the_chief #2 | 15 | 15 | floor | 0.0000000000e+00 | 0.0000000000e+00 | +0.00e+00 | +0.00e+00 | same |
| tBeamRows | test_centroid_option_scores_the_ray_centroid_not_the_chief #3 | 15 | 15 | floor | 0.0000000000e+00 | 0.0000000000e+00 | +0.00e+00 | +0.00e+00 | same |
| tCalib | test_calib_runs_baseline | 8 | 10 | cap | 2.0679515314e-35 | 2.0679515314e-35 | -2.58e-16 | -7.70e-34 | better |
| tCalib | test_calib_converges_after_perturbation | 8 | 10 | cap | 1.9246557550e-50 | 1.9246557550e-50 | +0.00e+00 | +0.00e+00 | same |
| tCalib | test_calib_returns_expected_schema | 8 | 10 | cap | 2.0679515314e-35 | 2.0679515314e-35 | -2.58e-16 | -7.70e-34 | better |
| tCalib | test_calib_set_var_elt_via_name_list | 8 | 10 | cap | 2.0679515314e-35 | 2.0679515314e-35 | +0.00e+00 | +0.00e+00 | same |
| tCalib | test_calib_set_var_elt_via_positional_mask | 8 | 10 | cap | 2.0679515314e-35 | 2.0679515314e-35 | +0.00e+00 | +0.00e+00 | same |
| tCalib | test_calib_set_iter_and_tol | 1 | 1 | cap | 1.6750407422e-32 | 1.6750407422e-32 | +0.00e+00 | +0.00e+00 | same |
| tCalibTargets | zmode_leg | 8 | 10 | cap | 2.4458612778e-43 | 1.5451187973e-43 | -3.68e-01 | -1.01e-22 | better |
| tCalibTargets | zmode_leg #2 | 8 | 10 | cap | 5.7913038096e-53 | 5.7913038096e-53 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_optimize_freeform_corrects_injected | 43 | 13 | floor | 1.6131400552e-28 | 1.6131400552e-28 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_convex_conic_seed_zero_and_optimizes | 106 | 32 | floor | 7.0547176220e-15 | 7.0547176220e-15 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_base_sphere_zernike_correction | 120 | 120 | cap | 1.3376382671e-09 | 1.3376382671e-09 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_tma_multifield_optimize_native | 43 | 20 | floor | 2.3320688402e-16 | 2.3320690245e-16 | +7.90e-08 | +6.03e-16 | within OptTol |
| tDesignTelescope | test_field_bias_optimize_corrects_offaxis | 43 | 24 | floor | 4.4895156831e-16 | 4.4895156831e-16 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_optimize_rigid_body_reemits_moved_design | 30 | 30 | cap | 1.0913139245e-16 | 1.0913139245e-16 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_optimize_accepts_2d_field_cross | 43 | 34 | floor | 3.2201388619e-16 | 3.2201402832e-16 | +4.41e-07 | +3.96e-15 | within OptTol |
| tDesignTelescope | test_optimize_elts_subset | 29 | 12 | floor | 2.3110957175e-16 | 2.3110970567e-16 | +5.79e-07 | +4.40e-15 | within OptTol |
| tDesignTelescope | test_optimize_per_element_dofs | 30 | 30 | cap | 1.6603751923e-16 | 1.6603751923e-16 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_fourth_mirror_psi_parity | 60 | 41 | converged | 1.6586821376e-16 | 1.6586866779e-16 | +2.74e-06 | +1.76e-14 | within OptTol |
| tDesignTelescope | test_add_fold_after_optimize_keeps_conics | 60 | 41 | converged | 1.6586821376e-16 | 1.6586866779e-16 | +2.74e-06 | +1.76e-14 | within OptTol |
| tDesignTelescope | test_optimize_freeform_lmon_survives_readback | 22 | 13 | floor | 2.2735293273e-28 | 2.2735293273e-28 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_optimize_freeform_lmon_survives_readback #2 | 8 | 10 | floor | 2.2735293273e-28 | 2.2735293273e-28 | +0.00e+00 | +0.00e+00 | same |
| tDesignTelescope | test_optimize_area_grid_dedup | 43 | 45 | floor | 5.6838374310e-16 | 5.6838374310e-16 | +0.00e+00 | +0.00e+00 | same |
| tSpectrometerRx | TestRunner.evaluateMethodCore | 1 | 1 | cap | 1.6429499114e-12 | 1.6429499114e-12 | +0.00e+00 | +0.00e+00 | same |
| tSpectrometerRx | TestRunner.evaluateMethodCore #2 | 1 | 1 | cap | 1.1374663585e-09 | 1.1374663585e-09 | +0.00e+00 | +0.00e+00 | same |

**Rules this table rejected on the way** (full-precision full-T traces, scored exactly; the
derivation is in macos_f90/CLAUDE.md "CALIB's LM stopping rule"):
- floor K = 3: base_sphere_zernike_correction stopped at itr 9, chisq 481x the old (6.44e-7 vs
  1.34e-9; WFE 1.16e-4 vs 5.28e-6 m); convex_conic WFE +9.1e-12 m.  The longest mid-solve
  rejection run a real decrease later followed is 6 (convex_conic, itr 9-14).
- flat by the absolute rms test alone: five tDesignTelescope cases 'converged' at itr 3 with
  chisq +4..39% (fourth_mirror 1.6587e-16 -> 2.2984e-16, ~1e-10 m rms vs OptTol 1e-12 m) -- the
  first steps of an LM still ramping lambda are ~1e-3 relative = ~1e-12 m on a nm-class merit.
- a strict relative-chisq reading of the tolerance (1e-12): K >= 15 needed (runs of 11-14
  rejections at a converged design's residual, then 1e-7..6e-7 relative gains), floor saving 8%.
