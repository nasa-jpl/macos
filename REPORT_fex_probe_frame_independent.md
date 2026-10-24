# REPORT: FEX probe made frame-independent (2026-09-08)

**Ask (Dave):** "Make FEX's probe frame-independent."  Trigger: an element
STOP rebuilds the source frame right-handed (`define_local_csys` in
UpdSrcGrid) while an object-space STOP keeps the deck frame; on e5hex1
(deck `xGrid = -1 0 0`) that flipped xGrid, and FEX's single differential
probe `th = 5d-6*xGrid` flipped with it: EP radius 2548.0019 vs 2549.5813
(1.58 mm along the chief ray) for the SAME chief ray.

## Change
`FEXProbeCross` (tracesub_mod; used by FEX and SXP): four differential
chief rays, +/-5e-6 rad about two orthonormal axes perpendicular to the
source chief ray (azimuth seeded by xGrid, yGrid fallback), each crossed
with the chief ray (`FindCrossPt`, point ON the chief); the EP vertex is
`cr1pos + mean(v)*cr1dir`.  The +/- pair is a central difference, so the
term linear in the probe angle -- the frame-SIGN dependence -- cancels
exactly; the two axes average the tangential and sagittal pupils (the
medial pupil), which is rotation-invariant at first order, so the frame
AZIMUTH no longer matters.  Telecentric test = largest chief/probe sine;
all-probes-LOST is reported separately from PARALLEL; 1-3 lost probes
warn and average the survivors.  Every run prints the two axis pupils,
the medial, the T/S split and the +/- asymmetry, e.g. e5hex1:
```
 ***** FEX: EP crossing from4 probes (+/-axis1, +/-axis2), mm along the chief ray:
 *****   axis1  2548.79157      axis2  2498.67947      medial  2523.73552
 *****   T/S split 5.011E+01  +/- asymmetry 7.897E-01
```

## Verification
- e5hex1: object-space stop and Segment stop now give the SAME radius
  (2523.7355168 both; CLI and pymacos), where the legacy probe gave
  2548.0019 vs 2549.5813.
- Independent check with the LEGACY engine on 90-degree-rotated source
  frames (the legacy probe then measures the OTHER axis): e5hex1
  tangential 2548.00 / sagittal 2498.68 (split 49.3 mm; four-probe
  medial 2523.74); j18sc 2997.89 / 3000.29 (2.4 mm); 6MST -13730.9 /
  -12690.6 (1040 mm; four-probe spread 1043 mm).  The splits are the
  engine's own pupil astigmatism, not a probe artefact.
- Symmetric on-axis decks are unchanged to round-off (table).

## Corpus blast radius (45 decks: old_Rx, manual examples, pymacos/mmacos
fixtures, sensitivity templates; ifx, model 512; stop = deck ApStop or
`stop elt <first Reflector> 0,0`; pre = single-probe engine, post = this)

```
deck                                            f_pre         f_post     rel_df  |dVpt| mm    |dpsi|    spread
6MST                                       616.985771   13209.967810   2.04e+01  5.209e+02   0.0e+00       nan
6MST_segV3                                 235.162199     273.138154   1.61e-01  5.083e+02   0.0e+00       nan
6MST_wfs_segV3                              20.816411      20.772033   2.13e-03  4.438e-02   0.0e+00       nan
MillsCross                                      NOFEX          NOFEX
ape                                       1934.796953    1934.796953   0.00e+00  0.000e+00   0.0e+00       nan
btc3                                            CRASH          CRASH
btcNonSeg                                       NOFEX          NOFEX
dmt6mono                                   559.526546     559.871451   6.16e-04  3.449e-01   0.0e+00       nan
dmt6seg1313dm_centered                     561.734116     559.574596   3.84e-03  2.160e+00   0.0e+00       nan
eac1_opt_met_vcs                                NOFEX          NOFEX
eac2_7seg                                52597.595970   23314.171480   5.57e-01  4.255e+03   0.0e+00       nan
interpSurf                                      CRASH          CRASH
iris_dp_ZGD                              13796.387870   13805.760820   6.79e-04  9.373e+00   0.0e+00       nan
j18dcWithStop                             3066.575321    3068.384497   5.90e-04  1.809e+00   0.0e+00       nan
j18dc_double_pass                         3066.575321    3068.384497   5.90e-04  1.809e+00   0.0e+00       nan
j18mono                                   3037.064406    3037.968348   2.98e-04  9.039e-01   0.0e+00       nan
j18sa                                     3035.338204    3037.622431   7.53e-04  2.284e+00   0.0e+00       nan
j18sc                                     2997.889981    2999.081593   3.97e-04  1.192e+00   0.0e+00       nan
keckFF                                          CRASH          CRASH
lst3zern                                  1355.758453    1353.267585   1.84e-03  2.491e+00   0.0e+00       nan
nngx                                            NOFEX          NOFEX
wfpc3                                           NOFEX          NOFEX
wideAngle                                       CRASH          CRASH
CassWithExitPupil                            6.790668       6.790668   0.00e+00  0.000e+00   0.0e+00       nan
CoroExample                                  2.379164       2.379164   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Cass_FarField                             6.790668       6.790668   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro                                   7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_DM                                7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_FPM                               7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_FPM_Zern                          7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_FPM_Zern_vortex                   7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_FPM_Zern_vortex_oversized         7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_FPM_Zern_vortex_oversized_noLyot    7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_FPM_noLyot                        7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Coro_noLyot                            7096.674622    7096.674622   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Mask_Parabolas                            0.000000       0.000000   0.00e+00  0.000e+00   0.0e+00       nan
Rx_Mask_Parabolas_glb                        0.000000       0.000000   3.58e-01  2.066e-10   0.0e+00       nan
Rx_e5hex1                                 2548.001852    2523.735517   9.52e-03  2.427e+01   0.0e+00       nan
e2e_pie                                      1.976054       1.858813   5.93e-02  1.172e-01   0.0e+00       nan
e5hex1                                    2548.001852    2523.735517   9.52e-03  2.427e+01   0.0e+00       nan
Rx_Cass_NS                                   5.561146       5.560146   1.80e-04  1.000e-03   0.0e+00       nan
Rx_Cass_NSseq                                6.790668       6.790668   0.00e+00  0.000e+00   0.0e+00       nan
SegDemo3conic                                2.764995       2.764995   0.00e+00  0.000e+00   0.0e+00       nan
dwdgrid_5zoom_5fov_jwst_ote_designc_grid    2997.889981    2999.081593   3.97e-04  1.192e+00   0.0e+00       nan
jwst_ote_designc                          2997.889981    2999.081593   3.97e-04  1.192e+00   0.0e+00       nan
```

Reading: `rel_df` = relative change of the written EP radius; `|dVpt|` =
EP vertex move (mm, along the chief ray); `|dpsi|` = 0 everywhere (the
chief ray is untouched).  NOFEX = the deck did not load or FEX aborted on
BOTH binaries (pre-existing: MillsCross/btcNonSeg/nngx/wfpc3 fail to
load, eac1 aborts); CRASH = loader `forrtl severe (59)` on BOTH binaries
(btc3, interpSurf, keckFF, wideAngle) -- pre-existing, not chased here.

**Groups.**  (a) Bit-identical: every rotationally symmetric on-axis deck
(Cassegrains, the Rx_Coro family, SegDemo3conic, ape, Rx_Mask_Parabolas).
(b) Off-axis decks move by half their T/S pupil split: j18 family
3e-4..7.5e-4 (0.9-2.3 mm), dmt6mono 6e-4, iris_dp 7e-4, lst3zern 1.8e-3,
dmt6seg 3.8e-3, e5hex1 9.5e-3 (24 mm), e2e_pie 5.9e-2 (its EP sits at
x=-0.48 m after the Offner relay; split 0.235 m on a 1.98 m leg).
(c) Exotic: 6MST T/S split ~1 m on a 13.7 m radius -- the medial
crossing lands 96 mm from the next element, the beam-footprint guard
autoswitches to the legacy iEm1->EP leg (13210) where the single probe
had written 617 -- both fragile; 6MST_segV3 16%.  eac2_7seg: ALL FOUR
probe chief rays are LOST before element 45 (the rebuilt engine says so
explicitly), so FEX takes the station fallback (23314); the legacy 52597
-- the value the 2026-07-03 compat pass already flagged as a material
shift for review -- was computed from a lost probe ray that the legacy
code never checked.  Neither number is an exit pupil; that deck needs a
chief ray that survives to its EP return before FEX means anything.

## Pinned tests that move (need a REVIEWED re-pin, not done)
`tFocalSurface` on the jwst zoom fixture: `test_null_radii_are_pinned`
(1e-7: 3017.58 -> 3018.29..3018.74), `test_null_radii_match_the_ab_report_
stop_orders_now_agree` (5e-3 abs vs REPORT_wnom_cli_ab V4), `test_emitted_
deck_loads_and_its_fex_radius_follows_the_fit` (3017.5444 -> 3018.4638).
All three fail by exactly the medial shift (+0.7..1.2 mm).
`tPupilFindMethod` (found later the same day, 8/10): `test_field_scope_
places_per_combo_tilt_absorbing_spheres` pinned the written vertex
INVARIANT across the two zoom configurations to 1e-6 mm -- true for the
tangential probe (the FSM tilts about y, deflecting the beam in x, and
the y-plane crossing is symmetric under that), not for the medial, whose
sagittal half sees the deflection: measured 5.5e-4 mm.  `test_object_
space_apstop_deck_needs_no_stop_elt` pinned a 10-40 mm gap between
FEX and pupil_find's cone-convergence station on e5hex1 ("this deck's
smear is ~23 mm"); the medial FEX moved 24 mm toward that station and
the gap is now 1.3 mm -- the two independent finders AGREE, and the
"smear" was mostly the tangential-vs-medial offset.  **All five were
re-pinned to the medial values on 2026-09-08 (Dave's ruling), each with the
legacy value and the mechanism recorded in the test.**  Green:
tStopReload 2/2 (now asserts obj == Segment stop at 1e-9), tOptFex 3/3,
tPupilMap 13/13, tVeneerXP 5/5, tSpot 7/7, tDwDx 25/25 (the dw_dx
sensitivity channels reset_xp through FEX on e5hex1 -- their gates pin
Jacobian STRUCTURE and identities, not the EP radius, so they hold).
pymacos: module relinked, e5hex1 stop smoke agrees obj/Segment to 1e-16.
GMI regression not run: its fixtures use object-space ApStop with
ifFEX=0, so neither change reaches it.

## Not changed (scoped out, recorded in macos_f90/CLAUDE.md)
XPS keeps its single per-ray probe (its vertex can now differ from FEX by
the probe-sign term on asymmetric decks); the STOP-ELT entrance-pupil
crossing, FPP/PFP and the WINDOW/PLOCATE beam-frame probes define a frame
from xGrid by design.  The ELT-stop handedness flip itself remains.
