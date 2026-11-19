# A trace must not rewrite the prescription: the re-run of 3b's measurement (PLAN_CONSOLIDATION 1.1b, 2026-10-12)

Per deck (the 406-deck A/B list `decks_all.txt`, model 256): LOAD -> SAVE a -> `opd nElt` -> SAVE b in ONE CLI
process; the two saves compared key by key per element.  Pre = prefix10 (e8c80ea) ifx; post = this commit, ifx and
gfortran.  (3b's first measurement, gfortran, counted xObs on 105.)

| key rewritten by a trace | pre ifx | post ifx | post gfortran |
|---|---:|---:|---:|
| `PropType` | 3 | 3 | 0 |
| `SegXgrid` | 14 | 14 | 14 |
| `Surface` | 4 | 0 | 0 |
| `xGrid` | 48 | 48 | 48 |
| `xObs` | 106 | 0 | 0 |
| `yGrid` | 62 | 62 | 62 |
| decks loaded | 404 | 404 | 403 |

- **`xObs` and `Surface`: 0 after the fix** -- the item.  Trace numbers (ab_corpus_seed, 406 decks: status /
  nPass / lost / RMS / P-V) pre vs post: 0 movers on ifx AND gfortran.
- `xGrid` / `yGrid` / `SegXgrid`: unchanged, the documented one-time settle of a hand-written source frame
  (PLAN 1.1a 6a: settle-once, not a trace rewrite of deck intent).
- `PropType` (ifx only, pre AND post, so not this item): KeckDemo.in, mars.in, wideAngle.in declare the legacy
  NUMERIC `PropType= 0`; SAVE's name lookup prints out-of-table bytes for it -- garbage that differs between
  the two saves.  A SAVE-writer defect, recorded in PLAN_CONSOLIDATION 1.2.
- gfortran loads one deck fewer: `old_Rx/ape.in` (`KcElt= -1.003981000D`, a D with no exponent) is refused
  by the gfortran build pre-fix too -- a compiler difference, not this change.

Pre-fix decks rewritten (`Surface`):

- ~/dev/MACOS_sandbox/old_Rx/interpSurf.in
- ~/dev/MACOS_sandbox/old_Rx/Int_Surf.in
- ~/dev/MACOS_sandbox/old_Rx/MillsCross.in
- ~/dev/MACOS_sandbox/old_Rx/MillsFilled.in

Pre-fix decks whose `xObs` a trace rewrote (106):

- ~/dev/MACOS_res_dev/GMI/e2e_pie_met.in
- ~/dev/MACOS_res_dev/GMI/regression/Rx/e2e_pie_met.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_epd4060_stage1.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_epd4060_stage2.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_epd4060_stage3.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_epd4060_stage4.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_stage1.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_stage2.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_stage3.in
- ~/dev/MACOS_res_dev/mmacos/challenges/rodgers1/rodgers1_stage4.in
- ~/dev/MACOS_res_dev/mmacos/templates/10_telescopes/sz_tma/sz_tma.in
- ~/dev/MACOS_res_dev/mmacos/templates/10_telescopes/tma_offaxis/tma_offaxis.in
- ~/dev/MACOS_res_dev/mmacos/templates/10_telescopes/tma_onaxis/tma_onaxis.in
- ~/dev/MACOS_res_dev/mmacos/templates/20_segmentation/e5_pie/e5pie_polyap.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_pil150.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_pil75.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/zoom_5x5/dwdgrid_5zoom_5fov_jwst_ote_designc_grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/zoom_5x5/jwst_ote_designc.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_mono_back.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_mono_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_mono_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_back.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d040_back.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d040_dm.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d040_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d040_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d070_back.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d070_dm.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d070_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d070_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d110_back.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d110_dm.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d110_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d110_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_dm.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r3_grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r3_met.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r3_metopt.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s1_telescope.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s2_segmented.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_imager_ep.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_imager_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_mono_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_mono_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_seg_full.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_seg_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s4_grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/e2e_hex2.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/e2e_pie.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/e2e_pie_met.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/e2e_pie_metopt.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/s1_telescope.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/s2_instrument.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/s2_variants/dy-0.70/s2_instrument.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/s2_variants/dy-1.05/s2_instrument.in
- ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e/s4_grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/s2_test_aoi15.in
- ~/dev/MACOS_res_dev/pymacos/tests/Rx/e2e_pie.in
- ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Mask_Parabolas_glb.in
- ~/dev/MACOS_res_dev/pymacos/tests/sensitivities/e2e_pie/e2e_pie.in
- ~/dev/MACOS_sandbox/coro/hcit/d8Yb.in
- ~/dev/MACOS_sandbox/ctb/ctb_pil150.in
- ~/dev/MACOS_sandbox/ctb/ctb_pil75.in
- ~/dev/MACOS_sandbox/Design/RC_ex/rc_6m_f20.in
- ~/dev/MACOS_sandbox/Design/RC_ex/RC_ex2.in
- ~/dev/MACOS_sandbox/Design/RC_ex/RC_ex.in
- ~/dev/MACOS_sandbox/e2e_m2_12/s1_telescope.in
- ~/dev/MACOS_sandbox/ns_griddata_variants/iris_dp_conic.in
- ~/dev/MACOS_sandbox/ns_griddata_variants/iris_dp_GD.in
- ~/dev/MACOS_sandbox/ns_griddata_variants/optiix_dp_conic.in
- ~/dev/MACOS_sandbox/ns_griddata_variants/optiix_dp_gd.in
- ~/dev/MACOS_sandbox/ns_griddata_variants/optiix_dp_zgd.in
- ~/dev/MACOS_sandbox/old_Rx/6MST.in
- ~/dev/MACOS_sandbox/old_Rx/6MST_segV3.in
- ~/dev/MACOS_sandbox/old_Rx/6MST_wfs_segV3.in
- ~/dev/MACOS_sandbox/old_Rx/ape.in
- ~/dev/MACOS_sandbox/old_Rx/btc3.in
- ~/dev/MACOS_sandbox/old_Rx/btcNonSeg.in
- ~/dev/MACOS_sandbox/old_Rx/dmt6mono.in
- ~/dev/MACOS_sandbox/old_Rx/dmt6seg1313dm_centered.in
- ~/dev/MACOS_sandbox/old_Rx/eac1_opt_met_vcs.in
- ~/dev/MACOS_sandbox/old_Rx/eac2_7seg.in
- ~/dev/MACOS_sandbox/old_Rx/iris_dp_ZGD.in
- ~/dev/MACOS_sandbox/old_Rx/j18dcWithStop.in
- ~/dev/MACOS_sandbox/old_Rx/j18sa.in
- ~/dev/MACOS_sandbox/old_Rx/j18sc.in
- ~/dev/MACOS_sandbox/old_Rx/keckFF.in
- ~/dev/MACOS_sandbox/old_Rx/mars.in
- ~/dev/MACOS_sandbox/old_Rx/wideAngle.in
- ~/dev/MACOS_sandbox/slides/rodgers1_epd4060_stage1.in
- ~/dev/MACOS_sandbox/slides/rodgers1_epd4060_stage4.in
- ~/dev/MACOS_sandbox/slides/rodgers1_stage1.in
- ~/dev/MACOS_sandbox/slides/rodgers1_stage4.in
- ~/dev/MACOS_sandbox/VSG/handoff_CCL_2026-07-24/opt_work/vsg_base_clean.in
- ~/dev/MACOS_sandbox/VSG/handoff_CCL_2026-07-24/opt_work/vsg_optimized.in
- ~/dev/MACOS_sandbox/VSG/handoff_CCL_2026-07-24/vsg_rough.in
- ~/dev/MACOS_sandbox/xp_tst/ab_cli/cli_C.in
- ~/dev/MACOS_sandbox/xp_tst/ab_cli/cli_LL.in
- ~/dev/MACOS_sandbox/xp_tst/ab_cli/cli_LR.in
- ~/dev/MACOS_sandbox/xp_tst/ab_cli/cli_UL.in
- ~/dev/MACOS_sandbox/xp_tst/ab_cli/cli_UR.in
- ~/dev/MACOS_sandbox/xp_tst/j18sa.in
- ~/dev/MACOS_sandbox/xp_tst/jwst_ote_designc.in
- ~/dev/macos/ZGD_test_files/e2e_pie.in
