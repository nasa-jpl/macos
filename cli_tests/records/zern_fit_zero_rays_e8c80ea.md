# zern_fit: exact-zero VALID rays in the corpus (PLAN_CONSOLIDATION 1.2, 2026-10-12)

Measured BEFORE the fix, on an instrumented side build of e8c80ea (ifx, model 256; a
WRITE of COUNT(opdValid==1 .AND. OPDMat==0) after SUBROUTINE OPD -- not committed), over
the 406-deck A/B corpus (`decks_all.txt`, the list the phase-1/2 A/B records used; its
resources paths point at the `~/dev/MACOS_res_dev` worktree).  Each deck: LOAD, `opd nElt`
under the default CHIEF-ray reference.  These are the rays `zern_fit` dropped from every
ZCO / ZRM / ZAB fit (its rows were `OPDin /= 0`).

| | decks |
|---|---|
| loaded | 404 of 406 (2 load_fail by design) |
| produced an OPD map (the count line printed; the other 24 printed neither it nor an nPassRays line) | 380 |
| >= 1 exact-zero valid ray | **118** |
| exactly 1 (not located per deck; on odd on-axis grids it is plausibly the ray at the chief's pupil point) | 88 |
| more than 1 | 30 (6 of them every ray: flat maps) |

## More than one

| zero rays | valid rays | RMS OPD | deck |
|---:|---:|---|---|
| 86 | 1257 | 9.905059341D-16 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_imager_leg.in |
| 224 | 1257 | 1.498120362D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_mono_back.in |
| 224 | 1257 | 1.498120362D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m/s3_seg_back.in |
| 308 | 1257 | 2.718555296D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_mono_back.in |
| 308 | 1257 | 2.718555296D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_back.in |
| 122 | 1257 | 2.797240158D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d040_back.in |
| 294 | 1257 | 2.978039662D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d070_back.in |
| 204 | 1257 | 2.960857553D-15 | ~/dev/MACOS_res_dev/mmacos/templates/80_end_to_end/e2e6m_r2/r1_seg_d110_back.in |
| 12453 | 12453 | 0.000000000D+00 | ~/dev/MACOS_res_dev/mmacos/tests/Rx/Rx_PolElt.in |
| 12453 | 12453 | 0.000000000D+00 | ~/dev/MACOS_res_dev/mmacos/tests/Rx/Rx_PolElt_Ref.in |
| 9899 | 12453 | 1.910074152D-16 | ~/dev/MACOS_res_dev/mmacos/tests/Rx/Rx_PolElt_Tilt.in |
| 12453 | 12453 | 0.000000000D+00 | ~/dev/MACOS_res_dev/mmacos/tests/Rx/Rx_Refract.in |
| 1126 | 12453 | 6.368106544D-17 | ~/dev/MACOS_res_dev/mmacos/tests/Rx/Rx_Refract45.in |
| 48592 | 48592 | 0.000000000D+00 | ~/dev/MACOS_res_dev/mmacos/tests/Rx/Rx_VecChain.in |
| 388 | 2828 | 2.948851589D-12 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Coro.in |
| 412 | 3096 | 2.931599397D-12 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Coro_DM.in |
| 142 | 1020 | 3.016202550D-12 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Coro_FPM_Zern_vortex.in |
| 92 | 708 | 3.103888975D-12 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Coro_FPM_Zern_vortex_oversized.in |
| 10905 | 65536 | 2.795428840D-12 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Coro_FPM_Zern_vortex_oversized_noLyot.in |
| 412 | 3096 | 2.931599397D-12 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Coro_noLyot.in |
| 31981 | 50615 | 1.974830222D-14 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Mask_Parabolas.in |
| 16424 | 50615 | 2.934418196D-14 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_Mask_Parabolas_glb.in |
| 48592 | 48592 | 0.000000000D+00 | ~/dev/MACOS_res_dev/pymacos/tests/Rx/Rx_VecChain.in |
| 19 | 1184 | 8.441579404D-11 | ~/dev/MACOS_sandbox/Design/RC_ex/RC_ex.in |
| 19 | 1184 | 8.441579404D-11 | ~/dev/MACOS_sandbox/Design/RC_ex/RC_ex2.in |
| 10 | 1184 | 8.441641823D-11 | ~/dev/MACOS_sandbox/Design/RC_ex/rc_6m_f20.in |
| 388 | 2828 | 2.948851589D-12 | ~/dev/MACOS_sandbox/coro/hcit/coro.in |
| 4 | 6286 | 6.621617200D-08 | ~/dev/MACOS_sandbox/coro/hcit/d8Yb.in |
| 50617 | 50617 | 0.000000000D+00 | ~/dev/macos/ZGD_test_files/tst_ZGDdata2.in |
| 2173 | 3001 | 6.903706658D-16 | ~/dev/macos/docs/macos-manual/examples/HOEExample.in |

## Exactly one

- ~/dev/MACOS_res_dev/GMI/ff_pie.in
- ~/dev/MACOS_res_dev/GMI/regression/Rx/Rx_e5hex1.in
- ~/dev/MACOS_res_dev/mmacos/templates/20_segmentation/e5_pie/e5pie_polyap.in
- ~/dev/MACOS_res_dev/mmacos/templates/20_segmentation/e5_seg/e5_seg_met.in
- ~/dev/MACOS_res_dev/mmacos/templates/20_segmentation/e5_seg/e5_seg_metopt.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_dcr.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_dcr_gen.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_dm.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_pil150.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_pil75.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_opt.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_prop.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_prop_dcr.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_seed.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageA.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageB2.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageB3.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageC.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageD.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageE.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_planar_stageF.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/ctb_s2s_dcr_gen.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/dst2/dst2_F01_norm.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/dst2/dst2_F02_norm.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/dst2/dst2_F03_norm.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/dst2/dst2_F04_norm.in
- ~/dev/MACOS_res_dev/mmacos/templates/30_instruments/bench_ctb/dst2/dst2_F05_norm.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_ifo/ifo_ref_arm.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_ifo/ifo_test_arm.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_ifo_dm/ifo_dm_base_arm.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_ifo_dm/ifo_dm_ref_arm.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_ifo_dm/ifo_dm_test_arm.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_layout/bench_layout_opt.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_layout/bench_layout_seed.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_layout/bench_layout_stageA.in
- ~/dev/MACOS_res_dev/mmacos/templates/40_benches/bench_layout/bench_layout_stageB.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/e5hex1/e5hex1.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/e5hex1/e5hex1_grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/e5hex2_refzern/e5hex2.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/e5hex2_refzern/e5hex2grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/run_dwdgrid_multi/e5hex1_grid.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/run_dwdsurf_multi/e5hex1.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/run_dwdx_multi/e5hex1.in
- ~/dev/MACOS_res_dev/mmacos/templates/50_sensitivities/run_dwdz_multi/e5hex1.in
- ~/dev/MACOS_res_dev/mmacos/templates/60_visualization/view_rx_demo/e5hex1.in
- ~/dev/MACOS_res_dev/mmacos/templates/60_visualization/view_rx_demo/e5mono_met.in
- ~/dev/MACOS_res_dev/mmacos/templates/60_visualization/view_rx_demo/e5pie.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/ifo_ref_pol.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/ifo_test_pol.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/s2_ref_aoi15.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/s2_ref_aoi45.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/s2_test_aoi15.in
- ~/dev/MACOS_res_dev/mmacos/templates/90_polarization/bench_ifo_pol/s2_test_aoi45.in
- ~/dev/MACOS_res_dev/pymacos/tests/Rx/Grating_example_001.in
- ~/dev/MACOS_res_dev/pymacos/tests/Rx/e5hex1.in
- ~/dev/MACOS_res_dev/pymacos/tests/sensitivities/e5hex1/e5hex1.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/e5hex2.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/e5mono.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/e5mx.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/e5pie.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/e5seg1.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/ff_hex2.in
- ~/dev/MACOS_res_dev/segmirmaker/test_in/xx.in
- ~/dev/MACOS_sandbox/VSG/handoff_CCL_2026-07-24/opt_work/vsg_base_clean.in
- ~/dev/MACOS_sandbox/ctb/ctb_pil150.in
- ~/dev/MACOS_sandbox/ctb/ctb_pil75.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_opt.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_seed.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageA.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageB2.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageB3.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageC.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageD.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageE.in
- ~/dev/MACOS_sandbox/ctb/ctb_planar_stageF.in
- ~/dev/MACOS_sandbox/old_Rx/Int_Surf.in
- ~/dev/MACOS_sandbox/old_Rx/MillsCross.in
- ~/dev/MACOS_sandbox/old_Rx/MillsFilled.in
- ~/dev/MACOS_sandbox/old_Rx/interpSurf.in
- ~/dev/macos/ZGD_test_files/e5hex1.in
- ~/dev/macos/ZGD_test_files/e5hex2.in
- ~/dev/macos/ZGD_test_files/e5pie.in
- ~/dev/macos/ZGD_test_files/ff_hex2.in
- ~/dev/macos/ZGD_test_files/opt_example_asph.in
- ~/dev/macos/ZGD_test_files/tst_SegMirMaker.in
- ~/dev/macos/ZGD_test_files/tst_ZGD2.in
- ~/dev/macos/docs/macos-manual/FreeFormPieAperture.in
- ~/dev/macos/docs/macos-manual/examples/Luneberg.in
