<!--
deck_gauges_material.md -- RAW MATERIAL for "DM Surface Gauge Comparison"
(JPL HWO WFS&C discussion group).  Numbers only, one section per outline slide.
No slide text.  Sources, read in full 2026-09-14:
  R-IFO  = MACOS_resources/mmacos/templates/40_benches/tg_psi_dm96_oap/REPORT_gauge_ifo.md
  R-OAP  = .../tg_psi_dm96_oap/REPORT_oap.md      RM-IFO = .../tg_psi_dm96_oap/README.md
  R-PDI  = .../40_benches/pdi_dm96/REPORT_gauge_pdi.md   RM-PDI = .../pdi_dm96/README.md
  RM-ZW  = .../40_benches/zwfs_dm96/README.md     D-ZW  = macos/demo_session/deck_zwfs.md
  PLAN   = macos/BRIEF_gauge_deck.md (7, 7.1, 7.1b, 7.2, 7.3, 10, 11)
Figure paths ls-checked; pixel sizes from `file`.  American English.
Readings: lens IFO / OAP IFO (Twyman-Green four-step); L ZWFS linear 1 frame;
I+ ZWFS exact 1 frame w/ sign prior; S ZWFS stepped 4 frames; V vector pair
2 frames; P stepped pinhole (common path); PF P/SRI (waveguide reference arm).
-->
# Material for "DM Surface Gauge Comparison"

## 1. Title
Nothing to extract.  PLAN 5.1: title **DM Surface Gauge Comparison**; audience the JPL HWO WFS&C discussion group.

## 2. The requirement + the four scores
- DM: **96x96 actuators, 1 mm pitch, 96 mm aperture** (48x48 @ 2 mm alongside). [RM-ZW]
- Working surface of record: **30 nm rms** random, seed 7; ladder rungs 30/40/50/60/80/100/120/160/240/480 nm.
- On-orbit spec (Dave, S11): "the DM surface needs to remain constant to **<< 10 pm**, with frequent remeasurement and closed-loop DM actuator servo control". [RM-ZW S11]
- Ultimate measurement-error target **~1 pm** (Dave 2026-09-04). [RM-ZW Rulings]
- Initial figure to capture: **100-200 nm WFE = 50-100 nm of surface**. [PLAN 5.5, 7]
- Wrap limit of every phase reading: **+-158 nm of surface** (+-pi at 632.8 nm double pass). [R-PDI 8]
- Photon scale: 1e14 photons at 633 nm = **31 uJ**; a 1 mW laser at the arm's 25% throughput = **0.13 s** per measurement. [D-ZW]

The four scores, verbatim from D-ZW's "How the numbers are scored" slide:
- *Gain*: "recovered change divided by the true change.  1.0 is perfect.  Gain is quoted after calibration."
- *Floor*: "the spread (rms) of the actuators that were not changed, in picometers."
- *SNR*: "recovered change divided by the floor.  Above 5 counts as detected."
- *Working surface*: "a random DM shape, 30 nm rms unless stated, present during the test."
- *Measurement / photons*: "a measurement is one DM shape measured once ... Photon budgets are quoted per measurement, all its frames summed."  Measuring a *change* costs two measurements.
- *Capture range to 10%* [R-PDI 2]: "the largest working-surface rms at which a 10 nm change on the 47 grid sites reads within 10% of its size (gain in 0.9...1.1), log-interpolated between ladder rungs."
- *Closed-loop one number* [RM-ZW S11]: photons per cycle to hold **3 pm rms**, gain 0.5, 60 cycles.

## 3. Shared front end
Bench parameters [R-IFO 6; RM-IFO; RM-ZW]:
- source: spatially filtered HeNe **632.8 nm**, beam radius **51 mm** collimated.
- collimator L1 **f 857 mm, dia 103 mm**, AR.  focuser L2 **f 429 mm, dia 103 mm**, AR.  field lens FL **f 43 mm, dia 21 mm**, AR.
- beamsplitter plate: AOI **7 deg**, **2.6 mm** thick, 50/50 (polarizing).  Compensator plate matched, **171 mm** from BS, AR.
- DM (test object): **96x96, 96 mm aperture, 1 mm pitch, 700 mm leg**, protected Al.  Reference flat + PZT: dia >= 103 mm, **~564 mm** leg, protected Al.
- camera: **385 px per pupil** (1 Mpix class); 193 px per pupil = development sampling.  Detector px per actuator **2.5** at 193 rays, **5.0** at 385.
- mask seat at the internal focus of the detector leg, **F/4.17**; lens-rig MASK_TRIM **-5.58 mm** off the thin-lens focus.
- DM -> detector magnification **10.158 DM-mm per detector-mm** (ZWFS/PDI test arm); 10.690 on the P/SRI bench.  Pupil image **9.4 mm** dia, **32.4 mm** behind FL.
- scale off the 56 mm v1 rig: s = 96/56 = **1.714**.

Figures: `/home/dcr/dev/macos/demo_session/figs/zwfs_render_rig.png` **2657 x 875** (D-ZW slide 3; identical copy at `.../zwfs_dm96/zwfs_render_rig.png`).  `.../tg_psi_dm96_oap/lens_vlayout.png` **3958 x 3250** (recipe, 3 panels, 2026-09-14).
**MISSING:** a redone *front-end-only* figure in the recipe.  PLAN 4 lists "shared front end | zwfs_render_rig.png -- redo in the recipe" as owed by CCL; it does not exist.  `lens_vlayout.png` is the nearest (whole front end plus reference arm).
Parts = the bullets above. [R-IFO 6; RM-IFO]

## 4. Interferometer, best configuration (lens rig)
**Phase-shift form picked: the HYBRID** -- polarization snapshot for the change measurements, PZT four-step as the absolute calibrator. [R-IFO 3]  Numbers that pick it:

| form | error | number | tag |
|---|---|---|---|
| PZT, 2% step error | single-actuator gain | **0.9743** (-2.6%); grid 0.9905; dense 0.9885; floor 2->5 pm | `lens_deck_se2` |
| PZT, 5% step error | single-actuator gain | **0.9543** (-4.6%); grid 0.9962; dense 0.9915; floor ->9 pm | `lens_deck_se5` |
| PZT 2% in hold | 3 pm noise-only / walk | **5.4e12 / 1.7e13** (record 5.5e12 / 2.0e13) -- common mode | `loop_lens_se2` |
| PZT + camera 1/f walk 1e-3 of signal, 25% intra | 3 pm noise-only | **6.1e12** vs 5.4e12 (+13% light); spectrum [<4, 4-12, >12 cyc/ap] = [0.02, 0.07, 0.93] pm; exactly immune at cam_intra 0 | `loop_lens_cam` |
| PZT + DM within-scan walk 2 pm/cyc, 25% intra | 3 pm under walk | **1.6e13** (intra 0: 1.7-2.0e13); ss 3.11 vs 3.21 pm at 1e13 | `loop_lens_intra` |
| snapshot v1 (plate) | BS diattenuation | test arm rotated **+7.479 deg**, PSI scale gain **1.11661 (+11.7% high)**; nulled to **1.00000** by analyzer sweep + **+3.768 deg** waveplate clock, fringe visibility 0.996; residual 4-theta term 8.9e-4 of the fringe, **1.7e-14 nm** in the differential | `tg_psi_dm` |
| snapshot v2 (MacNeille cube) | diattenuation structurally absent | arms **5.3e-6 deg** from orthogonal; gain **0.999999**; visibility **1.000000**; **2.27x** delivered power; naive odd stack H(LH)^4 = **2.11e-2 (2.1%) R_p**, symmetric (1/2H L 1/2H)^4 = **R_p 0**, T_p/T_s **2382:1** | `tg_psi_dm_v2` |

Rows on the 30 nm surface, matrix measured ON it, tag `lens_deck` [R-IFO 1]: single 10 nm at the hold-out site **0.9911 / 2 pm / SNR 6076**; 1 nm on the **52** grid sites (4:8:96) **0.9895 / 1 pm / SNR 700**; dense random 10 nm **0.9895, err 168 pm**.
> NOTE (site count): R-IFO 1 -- "the lit set puts **52** actuators on the every-8th grid ... the count is the tg96 pupil's, **not 47**."  Every ZWFS/PDI row uses **47**.  Do not merge silently.

Loop (`loop_lens`, D7) [R-IFO 4; R-OAP D7]: 3 pm noise-only **5.5e12** photons/cycle; 3 pm under the 2 pm walk **2.0e13**; thermal floor **13.1 pm** (bias 12.9 = lag; >12 cyc/ap band 8.7 pm); noiseless 1 nm step floors at **87 pm ~ 8.6%**, rising 66 -> 87 over the last 30 cycles (10 nm step 864 pm, linear); rho 0.511, tau 1.5; sig_n **11.97 pm at 1e12** -> 1 pm at **~1.4e14**.  Noise vs theory sig_n*sqrt(g/(2-g)): 7.01/6.91, 2.22/2.19, 0.70/0.69, 0.22/0.22 pm at 1e12...1e15.
Layout figures: `.../tg_psi_dm96_oap/lens_vlayout.png` **3958 x 3250** (whole train; BS/compensator/polarization node; tail focuser->mask seat->field lens->camera).  `.../runs/loop_lens/loop_lens_loop.png` **2026 x 844**.  Regenerate with `tg96_run('stages',{'bench','figs'})`.
Parts, lens rig [R-IFO 6]: the 9 shared-front-end items in section 3, plus -- **snapshot form only:** input polarizer, arm QWPs (test/reference), output QWP, analyzer -- quarter-wave, azimuths **45 / 0 / 45 / 0 / 0 deg**, count 5; **v2 snapshot:** MacNeille cube **12.7 mm**, ZnS/cryolite on n_g **1.655**, symmetric stack, replacing plate BS + compensator + ideal polarizers.

## 5. Zernike sensor, best configuration (stepped reading S, matrix on the surface)
| row | S @ 193 rays (`matbase`, `pdi193fbase`) | S @ 385 rays (`matbase385`) |
|---|---|---|
| single 10 nm | **0.989 / 5 pm / SNR 2160** (`pdi193fbase` 0.9885 / 5 / 2120) | **0.99 / 5 pm** |
| 47 sites @ 1 nm | **0.999 / 4 pm / SNR 284** (`pdi193fbase` 0.9993 / 4 / 284) | **1.00 / 3 pm** |
| dense random 10 nm | **0.98 / 0.68 nm** | **0.98 / 0.67 nm** |
| ladder 30/40/50/60 nm (30 nm matrix) | **1.08 / 0.98 / 0.92 / 0.70** (floors 243 / 228 / 399 / 664 pm) | -- |
| capture range aging from 30 nm | **42 nm** (`cap385`) | -- |
| N(1 pm) | flat matrix **5.4e13** (`pdi193f`) / **5.6e13** (`v193noise`) / **1.0e14** (`rec193full`); on-surface **8.8e13** (`noise193_b30`) or **9.3e13** as quoted in R-PDI | -- |
| loop 3 pm noise-only / 2 pm walk | **2.6e12 / 7.5e12**; thermal 10.05 pm; step **0.000 pm**; rho 0.61-0.63 | `loop385`: 2.6e12 / 7.7e12, thermal 10.1 pm |

Figures: `demo_session/figs/zwfs_render_rig.png` **2657 x 875**; `demo_session/figs/zwfs_mask385_mask.png` **1230 x 515** (= `.../zwfs_dm96/runs/mask385/mask385_mask.png`); older `.../zwfs_dm96/zwfs_layout.png` **1467 x 539**, `.../zwfs_mask_fig.png` **1387 x 604**.
Parts -- mask substrate + camera [RM-ZW; D-ZW]: fused-silica plate, 3x3 array of etched dimples (VSG2 nine-spot), one in the beam at a time, class Thorlabs W4101FT1; **etch depth 346.2 nm** in fused silica = 2*pi*(n-1)d/lambda = **1.571 rad** (quarter wave) at 632.8 nm; **dimple of record 2.0 lamF/D = 5.3 um** at F/4.2 (**7.92 px** at the mask plane at both 193/1024 and 385/2048), enclosing **70.7%** of the focused light; other spots 1.06 lamF/D = **2.79 um**, 1.22 = **3.22 um**, 3.0 = **7.9 um**; sampling rule **dimple >= 6 px**; camera **385 px per pupil**, 5.0 detector px per actuator; **4 frames** (three etch depths + a clear frame).

## 6. Vector Zernike sensor (reading V)
Rows, 30 nm surface, matrix on it (`v193base`, 193 rays) [RM-ZW V1; D-ZW]:

| row | V | S (contrast) |
|---|---|---|
| single 10 nm | **0.9935 / 4 pm / SNR 2830** (`pdi193fbase` 0.9935 / 4 / **2835**) | 0.989 / 5 / 2160 |
| 47 sites @ 1 nm | **0.9992 / 3 pm / SNR 345** | 0.9993 / 4 / 284 |
| dense random 10 nm | **0.9999 / 0.33 nm** | 0.984 / 0.68 nm |
| ladder 30/40/50/60 nm, raw | **0.997 / 0.996 / 0.991 / 0.981** (floors 20 / 41 / 74 / 113 pm) | 0.999 / 0.908 / 0.859 / 0.657 |
| capture range aging | **70 nm** (68 nm with the flat matrix) | 42 nm |
| N(1 pm) | **4.7e13** (`v193noise`); on-surface @30 nm **6.1e13** (`noise193_b30`) | 5.6e13 / 8.8e13 |
| loop 3 pm noise-only / 2 pm walk | **1.5e12 / 5.3e12**; thermal 9.9 pm; step **0.000 pm**; rho 0.509 | 2.6e12 / 7.5e12 |
| G4 fold gate, 100 nm pokes (12 nm rms) | **0.053 pm** (single +phi frame errs by **9.0 nm**) | -- |
| local sensitivity move with the surface | 0.15 rel / 0.90 amplitude | S 0.29 / 0.73; L 0.60 / 0.49 |

Figures: `demo_session/figs/zwfs_vlayout.png` **2438 x 1368** (both channel decks traced; = `.../zwfs_dm96/zwfs_vlayout.png`); `demo_session/figs/crop_zwfs_vlayout_tail.png` **1808 x 1020** (focus -> two cameras from above; purple = channel A transmitted, orange = channel B reflected).
Parts [RM-ZW]: geometric-phase (half-wave) **metasurface** in the etched plate's seat, +pi/2 on one circular state and -pi/2 on the other; **quarter-wave plate** behind the field lens, fast axis **45 deg** between the cube's s and p; **12.7 mm cemented MacNeille cube** (`Bench.pbs_cube` + `pbs_macneille`); **two cameras**, both **20.7 mm** behind the cube at the pupil image (the cube's glass lengthens the 32.4 mm image distance by **5.0 mm** on both ports); two cameras rather than a Wollaston because the 9.4 mm pupil image 32 mm behind FL would need a **17 deg** split, past a calcite prism's limit; decks `zwfs_v_camA.in` / `zwfs_v_camB.in` (the engine does not split rays).

## 7. Point-diffraction, best configuration
**Picked: the stepped pinhole P, with a pinhole-only (shutter) frame per state and the five-frame Schwider-Hariharan scan** -- six frames a measurement, all common path, one plate in the mask seat. [R-PDI 0]
Why P and not PF [R-PDI 0]: PF's reference is surface-independent, gain inside **0.7%** over a **16x** range of working surface (`pfdeck`); but P + a shutter frame buys the same range in the common path (**1.02 / 1.06 / 1.13 at 120 / 240 / 480 nm**, `pdi193state`) for one extra frame; PF costs **~2x** the light traced (**5.1e12 / 1.5e13** vs P **2.3e12 / 7.0e12** photons per cycle), a second arm to build and balance, and its own reference-arm drift.  PF's role is **capture**, not hold.

| row (30 nm surface, matrix on it) | P (`pdi193fbase`) | PF traced (`pfdeck`) | PF synthesized (`pdi193fbase`) |
|---|---|---|---|
| single 10 nm | **0.9935 / 4 pm / SNR 2790** | 0.9935 / 2 pm / 4895 | 0.9935 / 4 / 2842 |
| 47 sites @ 1 nm | 0.9992 / 3 pm / 345 | 0.9924 / 1 pm / 878 | 0.9992 / 3 / 346 |
| dense random 10 nm | 0.9985 / 338 pm | 0.9926 / 144 pm | 1.0002 / 330 pm |
| N(1 pm), flat matrix | **3.3e13** | 2.00e14 | 1.9e14 |
| N(1 pm), on-surface @30 nm | **9.8e13** | 3.72e14 | -- |
| loop 3 pm noise-only / 2 pm walk | **2.3e12 / 7.0e12** | 5.1e12 / 1.5e13 | 7.0e12 / 2.5e13 |
| 2% step error, 12 nm figure, 5-frame SH | **4.9 pm** (4-step LS 421 pm) | **2.2 pm** (4-step LS 251 pm) | -- |

Pinhole diameter of record **2.0 lam/D** [R-PDI 3a-3c]:

| | 2.0 lam/D, model 1024 / 193 rays | 1.0 lam/D, model 2048 / 385 rays |
|---|---|---|
| physical diameter | **5.27 um** (F/4.17) | **2.64 um** |
| surround transmission t (amplitude) | **0.719** (0.52 in power) | **0.282 / 0.283** |
| pinhole coupling eta / throughput | 0.682 / **0.821** | 0.240 / 0.294 |
| px at the mask plane | 7.92 (PASS, >=6) | 3.96 (NOT MET) |
| single 10 nm | 0.9935 / 4 pm / SNR 2790 | 0.9937 / 4 pm / SNR 2967 |
| capture range | 62 nm | 120 nm+ |
| N(1 pm) on-surface / 3 pm under the walk | **9.8e13** / **< 1e13** | 2.1e14 / 3.6e13 |
| tags | `pin20_1024`, `pin20_loop` | `pin10_2048`, `pin10_loop` |

P/SRI bench balance (`psri_layout_fig`) [R-PDI 1a]: chief optical paths **4453.1061 mm each, difference 0.00e+00 mm**; compensator **21.857 mm** of n=1.5 glass; exit chiefs coincide to **0.00e+00 mm**, camera planes to **2.3e-13 mm**; 3210 of 3210 rays per arm at 65 rays.  Lr1 f **300 mm**, **F/2.92** on the 102.9 mm beam, conic **-0.5784**; pinhole seat **+1.096 mm** past the thin-lens focus; pinhole **3.69 um** (2.0 lam/D at F/2.92), **7.84 px** across.  Coupling eta **0.677**, throughput **0.806**, visibility **0.873**, flat reads 3.3e-16 rad rms.
Layout figures, CURRENT (recipe, 2026-09-13): `.../pdi_dm96/pdi_layout.png` **2438 x 1368** (common-path form); `.../pdi_dm96/psri_layout.png` **2438 x 1368** (panel 1 both arms, panel 2 the reference-arm node cropped); `.../pdi_dm96/psri_render.png` **2438 x 1030**.
**STALE copies in `demo_session/figs/`** (pre-2026-09-13 `Bench.sketch`, E-number labels at 9-11 pt; md5 differs -- do not use): `pdi_layout.png` 3662 x 1310, `psri_layout.png` 3232 x 1840, `psri_render.png` 3125 x 1094; and `pdi_layout_tail.png` 3068 x 1624, which **no longer exists** in `pdi_dm96/` (the recipe puts the crowded node in the same figure).
Parts -- common-path pinhole, reading P [RM-PDI]: pinhole substrate = fused-silica plate in the mask seat, pinhole **5.27 um** (2.0 lam/D at F/4.17, 632.8 nm); surround attenuated to **t = 0.719** in amplitude (**0.52** in power), pinhole clear; phase stepping through 0, pi/2, pi, 3pi/2 (4 frames) **or** the 5-frame Schwider-Hariharan scan -pi...pi; camera at the reimaged pupil, ~2.5 detector px per actuator at 193 rays (~20 um px, 1 Mpix class, at 385).  Fallback at 1.0 lam/D: **2.64 um** pinhole, **t = 0.282**, throughput 0.29.
Parts -- P/SRI, reading PF [RM-PDI]; all Mach-Zehnder incidences **45 deg**: pickoff **BS2** plate 2.57 mm, n 1.5, CA >= 103 mm, **R 0.60 / T 0.40**; **Lr1** f **300 mm**, CA 102.9 mm -> **F/2.92**, plano-convex n 1.5, conic **-0.5784**, AR; **pinhole / waveguide seat 3.69 um**, opaque surround, at the true focus **+1.096 mm** beyond the thin-lens one; **photonic phase shifter**, thermo-optic, on the waveguide chip, reference arm only; **Lr2** = Lr1 mirrored about the pinhole, same f, same conic; folds **M1** (test) and **M3** (reference), flat 150 mm, AOI 45 deg, protected metal, M3 solved for coincident exit chiefs; **compensator** (test arm) **21.857 mm** of n = 1.5 glass at normal incidence, CA >= 103 mm, AR; **BS3** plate 2.57 mm, n 1.5, 50/50; **one camera** at the pupil image, chiefs **2e-12 mm** apart, camera planes **2e-13 mm** apart.

## 8. Performance side by side -- 30 nm surface, matrix measured ON it
| reading | single 10 nm: gain / floor / SNR | 47 (52) sites @ 1 nm: gain / floor | dense 10 nm: gain / err | tag |
|---|---|---|---|---|
| **lens IFO** | **0.9911 / 2 pm / 6076** | **0.9895 / 1 pm** (52 sites, SNR 700) | **0.9895 / 168 pm** | `lens_deck` |
| **OAP IFO (bare Al)** | **0.9956 / 2 pm / 5301** | **0.9576 / 1 pm** (52 sites, SNR 683) | **0.9497 / 2099 pm** | `oap_deck` |
| **L** | **1.04 / 21 pm / 492** (193 r); 1.05 / 23 pm (385 r) | 1.05 / 22 pm / SNR 47 | 1.04 / 4.7 nm | `matbase`, `matbase385` |
| **I+** | **0.78 / 74 pm** (385 r); 0.69 (193 r) | 0.90 / 35 pm | 0.86 / 5.2 nm | `matbase385` |
| **S** | **0.989 / 5 pm / 2160** (193 r); 0.99 / 5 pm (385 r) | **0.9993 / 4 pm / 284** | **0.984 / 0.68 nm** | `matbase*`, `pdi193fbase` |
| **V** | **0.9935 / 4 pm / 2830** (2835 in `pdi193fbase`) | **0.9992 / 3 pm / 345** | **0.9999 / 0.33 nm** | `v193base`, `pdi193fbase` |
| **P** | **0.9935 / 4 pm / 2790** | **0.9992 / 3 pm / 345** | **0.9985 / 338 pm** | `pdi193fbase` |
| **PF** traced | **0.9935 / 2 pm / 4895** | **0.9924 / 1 pm / 878** | **0.9926 / 144 pm** | `pfdeck` |
| **PF** synthesized | 0.9935 / 4 pm / 2842 | 0.9992 / 3 pm / 346 | 1.0002 / 330 pm | `pdi193fbase` |

Caveats for the slide: the IFO rows use **52** grid sites, the rest **47**.  `pfdeck` runs on the P/SRI's own two decks (magnification 10.690) and `pdi193fbase` on the ZWFS test arm (10.158), so the floor and SNR columns carry a bench term as well as a reference term [R-PDI 1b].  IFO rows come from `tg96_run` stage `deck`; ZWFS/PDI at 193 rays unless a 385 value is given.
> DISAGREEMENT with PLAN 3 (draft table, superseded by the reports): PLAN 3 lists lens IFO **0.99 / 4.8 pm** single and **0.988 / 686 pm** dense from `lens_base` on a **16 nm** base, OAP **0.996 / 2.8** and **0.95**.  R-IFO 1 is the re-run on the 30 nm surface (table above).  Use R-IFO.
Figures: `demo_session/figs/zwfs_m2048_battery_n96.png` **1588 x 699**; `demo_session/figs/pdi193fbase_battery_n96.png` **1588 x 699**.

## 9. Capture range
(a) Calibration AGING from 30 nm, 10 nm change on the 47 (52) grid sites -- capture range to 10%: **I+ 36 nm**, **S 42**, **L 44** (`cap385`); **P 62** with the flat |b|^2 (`cap385p`); **V 70** (`cap385`); **P with a shutter frame 480 nm+** (`pdi193state`); **PF 480 nm+**, within 1% at 100 nm (`cap385p`); **lens IFO 322 nm single site / 480 nm+ on the 52-site grid**, gain 0.99-1.03 to 480 (`lens_deck`); **OAP IFO 480 nm+ / 480 nm+** (`oap_deck`).
Ladder gains (aging, 385 rays) at 30 / 40 / 50 / 60 / 80 / 100 nm [RM-ZW; R-PDI 2a]: L 1.027 / 0.956 / 0.828 / 0.675 / 0.389 / 0.19; I+ 0.962 / 0.866 / 0.655 / 0.540 / 0.261 / 0.10; S 0.999 / 0.909 / 0.862 / 0.658 / 0.222 / 0.06; V 0.997 / 0.996 / 0.991 / 0.982 / 0.827 / 0.06; P 0.9970 / 0.9858 / 0.9678 / 0.9425 / 0.4674 / 0.0618 then **-0.0147 (folded)** at 120; PF 0.9967 / 0.9989 / 1.0012 / 1.0034 / 1.0078 / 1.0122, 1.0166 at 120, 1.0431 at 240, **1.0962** at 480.
> DISAGREEMENT -- lens IFO aging range, three published values: (i) **322 nm** single / **480 nm+** grid, R-IFO 2, **the current record**, after the 2026-09-13 differential bug fix (the deck differentials were `measf(base+dev) - measf(base)`, a difference of two separately-wrapped absolute maps; now the wrapped phase difference); (ii) **120 to 158 nm** "its four-step wrap", D-ZW capture slide from `lens_base` single site, whose ladder is 0.992 / 1.001 / 1.017 at 30 / 60 / 120 nm and 1.93 at 240 [RM-ZW]; (iii) **120-158 nm** single-site and **60-120** for the OAP, PLAN 3 (pre-fix).  RM-ZW's pre-fix OAP ladder: 0.996 / 0.998 at 30 / 60, 0.959 at 120 (map correlation 0.32), 0.475 at 240.

(b) Matrix RE-MEASURED on the surface, 1 nm on 47 sites -- gain / floor:

| reading | 60 nm | 90 nm | 120 nm | 160 nm | tag |
|---|---|---|---|---|---|
| L | 0.962 / 25 pm | 0.954 / 21 | 0.963 / 17 | 0.968 / 13 | `cap385_b*` |
| I+ | 0.827 / 150 | 0.950 / 20 | 0.962 / 17 | 0.968 / 13 | `cap385_b*` |
| S | 0.934 / 59 | 0.960 / 29 | 0.970 / 22 | 0.968 / 20 | `cap385_b*` |
| V | 1.000 / 3 | 0.957 / 47 | 0.969 / 23 | 0.967 / 21 | `cap385_b*` |
| P | 1.0003 / 3 / SNR 297 | 0.9581 / 41 / 23 | 0.9687 / 19 / 51 | 0.9670 / 17 / 56 | `cap385p_b*` |
| PF | 0.9999 / 3 / 378 | 1.0007 / 2 / 415 | 1.0014 / 2 / 458 | 1.0021 / 2 / 516 | `cap385p_b*` |
| lens IFO | gain **0.99**, floor 1 pm, corr **0.9999**, flat to 160 nm | | | | `lens_deck` |
| OAP IFO | gain **0.958**, corr **0.980** to 160 nm | | | | `oap_deck` |

(c) Photon cost of operating off null -- N(1 pm) at 30 / 60 / 120 / 160 nm, matrix measured on each surface, 193 rays, 6 realizations per point:

| reading | 30 nm | 60 nm | 120 nm | 160 nm | ratio to 30 nm | tag |
|---|---|---|---|---|---|---|
| L | 9.2e13 | 1.8e14 | 6.2e14 | 4.6e14 | ~5x | `noise193_b*` |
| S | 8.8e13 | 6.1e15 (fold outlier) | 1.3e15 | 3.4e15 | -- | `noise193_b*` |
| V | 6.1e13 | 1.1e14 | 2.2e15 | 2.4e15 | ~40x | `noise193_b*` |
| P | 9.77e13 | 1.73e14 | 2.95e15 (**x30**) | 1.75e15 (**x18**) | 18-30x | `noise193p_b*` |
| PF | 3.72e14 | 4.13e14 | 6.58e14 | 3.67e14 (**x0.99**) | within 1% | `noise193p_b*` |
| lens IFO | **2.5e14** (S5 form), stable ~1.4-2.5e14 to 160 nm | | | | -- | `lens_deck` |
| OAP IFO | **6.2e14** (~2.5x the lens), stable ~6-8e14 | | | | -- | `oap_deck` |

R-PDI 2c flags P's 120 nm point as non-monotonic -- read it as "tens of times worse", not a number (the floor column does the same: 41 pm at 90 nm, 19 at 120, 17 at 160).  D-ZW: at 160 nm V needs 2e15 photons per measurement = **0.8 mJ = 3 s of a 1 mW laser** at 25% throughput.
Figure: **none specific to capture range**; D-ZW slides 11-12 are tables only.

## 10. Photons for 1 pm and the closed-loop hold
N(1 pm) per measurement, matrix on the FLAT [RM-ZW S8/V1; R-PDI 1c, Concl. 2]: **P 3.3e13** (`pdi193f`); F 3.7e13 (`rec193full`); **V 4.7e13** (`v193noise`, `pdi193f`); L 5.4e13 (`rec193full`); **S 5.4e13** (`pdi193f`) / **5.6e13** (`v193noise`) / **1.0e14** (`rec193full`); I 6.4e13; I+ 8.8e13 (8.4e13 with a noise-free prior); PF synthesized 1.9e14 (`pdi193f`); **PF traced 2.00e14** (`pfdeck`); **lens IFO 2.5e14** S5 form / **~1.4e14** loop sig_n form / **8e14** the legacy S5 estimate quoted in D-ZW; OAP IFO 6.2e14 S5, ~4.9e14 loop form.
> NOTE (S's three values): `pdi193f` 5.4e13 and `v193noise` 5.6e13 are the S5 single-actuator scenario on the flat matrix; `rec193full` 1.0e14 is the same stage at development sampling in the ZWFS record run.  On-surface matrix at 30 nm: **8.8e13** (`noise193_b30`, RM-ZW) vs **9.3e13** as quoted in R-PDI Conclusion 2.  Quote the tag with the number.

Closed loop -- photons per cycle to hold 3 pm (gain 0.5, 60 cycles, matrix on the 30 nm surface, drift seed 77 shared):

| reading | noise only | 2 pm/actuator walk | thermal floor (5 pm/cyc ramp) | noiseless 1 nm step at cycle 60 | sig_n at 1e12 | tag |
|---|---|---|---|---|---|---|
| **V** | **1.5e12** | **5.3e12** | 9.9 pm (9.87) | **0.000 pm** | 6.4 pm | `vloop193` |
| **L** | 2.1e12 | 7.3e12 | 27.6 pm (27.64) | 1.2 pm, still falling | 7.4 pm | `loop193` |
| **P** | **2.3e12** | **7.0e12** | 9.86 pm | **0.000 pm** | 7.8 pm | `ploop193` |
| **S** | 2.6e12 | 7.5e12 | 10.0 pm (10.05) | **0.000 pm** | 8.2 pm | `loop193` |
| **PF** traced | 5.1e12 | 1.5e13 | 9.88 pm | **0.000 / 0.000 pm** | 11.7 pm | `pfdeck_loop` |
| **lens IFO** | **5.5e12** | **2.0e13** | **13.1 pm** | **87 pm (~8.6%), rising 66->87** | 11.97 pm | `loop_lens` |
| **PF** synthesized | 7.0e12 | 2.5e13 | 9.9 pm | 0.000 pm | 13.5 pm | `ploop193` |
| **OAP IFO** | **3.4e13** | **never** (floors 4.1 pm) | **39.0 pm** | **276 pm (~27.6%), rising 210->276** | 22.2 pm | `loop_oap` |
| **I+** | diverges | diverges | diverges | 99 nm (1 nm step); 178 nm (10 nm step) | -- | `loop193` |

Other loop numbers: contraction rho -- V/P/PF **0.509-0.511**, S 0.61-0.63, L 0.82-0.83, lens IFO 0.511 (tau 1.5).  Held-residual spectrum under the walk [<4, 4-12, >12 cyc/ap]: V/PF **0.25 / 0.72 / 2.20 pm**; lens IFO none 0.22, walk 2.29, thermal 8.69 pm; OAP thermal >12 band alone 33 pm.  Thermal lag theory rate/(gG) ~ 10.4 pm at G 0.98.  385-ray cross-check (`loop385`): L 2.1e12 / 7.2e12, S 2.6e12 / 7.7e12, thermal 24.7 / 10.1 pm.  7.5e12 photons/cycle = **2.4 uJ** = 9 ms of a 1 mW laser at 25% throughput [D-ZW].
Loop figures: `demo_session/figs/crop_zwfs_loop193_left.png` **796 x 677** (S residual per cycle, four light levels, + noiseless 1 nm step); `demo_session/figs/crop_zwfs_vloop193_right.png` **808 x 681** (hold error vs photons/cycle, S green vs V purple, none/walk/ramp); `.../zwfs_dm96/runs/loop193/loop193_loop.png` **1609 x 699**; `.../runs/vloop193/vloop193_loop.png` **1609 x 699**; `.../tg_psi_dm96_oap/runs/loop_lens/loop_lens_loop.png` **2026 x 844**; `.../runs/loop_oap/loop_oap_loop.png` **2026 x 844**; `demo_session/figs/zwfs_rec193full_noise.png` **749 x 593** (noise vs photons per measurement, per reading, with the 1 pm line).

## 11. Capturing the initial figure (the descent)
Largest start reaching 3 pm in 40 cycles, matrix measured at the start, identical at 1e13 and 1e15 photons/cycle, unwrap OFF and ON (`cap_nouw`, `cap_uw`; S/V/P rows bit-identical between the arms) [R-PDI 9a]: **L 30 nm, S 30, V 60, P 60, PF 60** -- unchanged by unwrapping.
Opening-differential diagnostic (`cap_nouw`, 193 rays; wrapped-map rms, 2-pi residues inside the mask, largest wrapped gradient rad/px):

| start (truth) | S | V | P | PF |
|---|---|---|---|---|
| 60 nm (30) | 16.4 nm, **8 res**, 3.14 | 25.6, 0, 1.97 | 23.3, 0, 1.77 | 27.0, 0, 2.08 |
| 100 (70) | 22.3, **0 res**, 1.83 | 26.5, 0, 2.10 | 26.3, 0, 2.09 | 59.8, **70 res**, 3.14 |
| 150 (120) | 24.1, 0, 1.58 | 28.5, 0, 1.79 | 28.4, 0, 1.78 | 75.4, **928 res**, 3.14 |
| 200 (170) | 24.2, 0, 1.62 | 28.7, 0, 1.88 | 28.5, 0, 1.88 | 77.2, **2580 res**, 3.14 |
| 300 (270) | 24.0, 0, 1.48 | 28.4, 0, 1.66 | 28.3, 0, 1.65 | 77.8, **4928 res**, 3.14 |

S, V and P past ~60 nm return **24-28 nm whatever the truth is**, zero residues: blind (focal reference collapsed), not folded.  Only PF folds.
PF from a 100 nm start (200 nm WFE), 1e15 photons/cycle [R-PDI 9b]: neither unwrap nor recal **64 095 pm** (`cap_nouw`); unwrap only **5 383 pm**, rho 0.615, 10 nm at cycle 6 (`cap_uw`); recalibrate only **63 520 pm** (`cap_nouw_recal`); **both 0.245 pm**, k(10 nm) 6, k(3 pm) **26** (28 at 1e13) (`cap_uw_recal`, `descent193`).
The recommended configuration (P + shutter frame, unwrapped, recalibrated every 10 cycles) from a 100 nm surface [R-PDI 9c]: **0.207 pm at 1e15**, **2.066 pm at 1e13**, 3 pm at cycle **26-27**, contraction **0.716**, **3** recalibrations (`cap_state_uw_recal`, `descent193s`).  At 150 nm: 10 nm at cycle 22, ends at **9 324 pm** -- ceiling between 100 and 150 nm.  Without recal, P-shutter tracks PF to four digits (5 383.2 vs 5 383.4 pm with unwrapping; 64 092 vs 64 095 without).
Recalibration cadence (`descent193f`, K 20, 1e15): every 2 cycles **6.5 pm** for **9** recals; every 5 **11.2 pm** for **3**; every 10 is enough -- the binding constraint is the cycle count.
The interferometer captures (lens rig, PZT four-step, unwrap ON, no recal; `descent_lens`) [R-IFO 5]:

| start rms (WFE) | r(1) | k to 10 nm | k to 3 pm | r(K) at 1e13 / 1e15 | rho |
|---|---|---|---|---|---|
| 60 nm (120) | 29.7 nm | 3 | 16 | 2.15 / 0.22 pm | 0.51 |
| 100 nm (200) | 69.7 nm | 4 | 18 | 2.14 / 0.21 pm | 0.53 |
| 150 nm (300) | 119.7 nm | 6 | 21 | 2.12 / 0.21 pm | 0.57 |
| 200 nm (400) | 169.7 nm | 7 | 25 | 2.11 / 0.21 pm | 0.62 |
| 300 nm (600) | 269.7 nm | 25 | 43 | 2.09 / 0.21 pm | 0.84 |

Recal-every-10 vs never are identical (300 nm: k(3 pm) 45 vs 43; r(K) 2.25 vs 2.09 pm).  **The OAP does NOT capture** (`descent_oap`, bare Al): from 60 nm it reaches 10 nm at cycle 3 then stalls at **~5.9 nm**; from 150 / 300 nm it never reaches 10 nm (stalls at **24 / 53 nm**).
Within-scan DM drift (`intra193_0` vs `intra193`), hold error pm rms over lit at 1e15 [R-PDI 4b; PLAN 7.1] -- 2 pm walk, DM still: L 2.38, S 2.32, V 2.32, P 2.32, PF 2.33; drift across the scan: L 2.38, **S 1.57**, V 2.32, **P 1.58**, **PF 1.71**.  5 pm thermal ramp, DM still: L 27.64, S 10.05, V 9.87, P 9.86, PF 9.88; across the scan: L 27.64, **S 7.16**, V 9.87, **P 6.70**, **PF 7.30**.  Stepped readings improve **26-32%**; L and V unchanged to the digit (the plumbing gate).  IFO equivalent (`loop_lens_intra`): 3 pm under the walk from **1.6e13** (intra 0: 1.7-2.0e13).  Camera within-scan drift has the opposite sign: S **5.4**, P **5.3**, PF **10.7 pm** at 1e15 (`pcam193ri`) against L **10.8 nm** and V **89 pm** (`pcam193r`).
P/SRI reference-arm walk (`rw193_1e3 / 1e2 / 1e1`; P the common-path control, bit-identical across all three) [R-PDI 5; PLAN 7.1b]: hold noise-only at 1e13 -- P 1.45 pm, PF 2.50 / 2.51 / 2.53; at 1e15 -- P 0.14, PF 0.25 / 0.25 / 0.25; under the 2 pm DM walk at 1e13 -- P 2.72, PF 3.42 / -- / 3.44; **photons per cycle to hold 3 pm under the walk: P < 1e13, PF 4.8e13 / 4.9e13 / 5.0e13** -- **4% for a hundredfold range of walk**.
Unwrapper (`dm_gauge_lib/dmg_unwrap`, Ghiglia & Romero JOSA A 11, 107 (1994)); gates tDmgLoop G13, **15/15**: wrapped ramp 1.5 waves exact to **8.5e-14 rad**, 0 residues; band-limited 1.5 waves PV on a disc **6.7e-13 rad**, 0 residues, max gradient 0.90 rad/px, 16 PCG iterations; 3.0 waves PV **1.3e-12 rad**, 1.81 rad/px; never-wrapped map passes through at 5.4e-14; 40 waves PV reports **644 residues**.  ~30 ms per differential on a ~200x200 box.  mmacos fast suite **469 pass / 0 fail**.  Non-disturbance: `uwoff_ref` bit-identical to `pfsmoke_ref`, v3dev G4 = 0.296 pm.
Figure: `demo_session/figs/gauge_modes_flow.png` **1800 x 1000** (tool `demo_session/gauge_modes_flow.m`).

## 12. Lenses vs OAPs
| | lens | OAP (bare Al) |
|---|---|---|
| single 10 nm differential (flat) | 0.9916 / 2.2 pm | 0.9950 / 2.1 pm |
| single 10 nm (30 nm surface) | 0.9911 / 2 pm / SNR 6076 | 0.9956 / 2 pm / SNR 5301 |
| 52-site 1 nm (30 nm surface) | 0.9895 / 1 pm / 700 | 0.9576 / 1 pm / 683 |
| dense random 10 nm (flat) | 0.9894 (169 pm) | 0.9497 (2099 pm) |
| modal transfer gain / cross-talk | 0.96-0.99 / **< 0.06** | 0.63-0.95 / **~0.18** (0.42 uncoated, RETIRED) |
| flat-DM null | **0.134 nm** | **13.09 nm** (cancels in the differential) |
| loop 3 pm noise-only / 2 pm walk | 5.5e12 / 2.0e13 | 3.4e13 / **never** (4.1 pm floor) |
| thermal floor / noiseless step | 13.1 pm / 87 pm (8.6%) | 39 pm / 276 pm (27.6%) |
| N(1 pm) | 2.5e14 (S5 form) | 6.2e14 (~2.5x) |
| descent | captures 60-300 nm start to ~2 pm in <= 43 cycles | does not capture (stalls 5.9 / 24 / 53 nm) |
| capture range aging | 322 nm single / 480 nm+ grid | 480 nm+ / 480 nm+ |

Alignment sensitivity, D4, 10 um decenter / 10 urad tilt one at a time [R-OAP D4]: OAP1 decenter null shift 16.1 nm, **1.61 nm/um**, single-diff gain 0.9929, resid 2.4 pm; OAP1 tilt 94.6 nm, **9.46 nm/urad**, 0.9929, 5.8 pm; OAP2 decenter 15.3 nm, **1.53 nm/um**, 0.9932, 2.3 pm; OAP2 tilt 102.7 nm, **10.27 nm/urad**, 0.9915, 2.5 pm.
Coatings [R-OAP D5, item B] -- uncoated (ideal reflector) -> bare Al -> protected Al: flat/single 0.9948/2.2 -> 0.9950/2.1 -> 0.9950/2.1 pm; flat/dense 0.7486/4848 -> **0.9497/2099** -> 0.9488/2118 pm; dark-25% columns **0.05 -> 0.82 -> 0.82**; flat-DM null 12.893 -> 13.089 -> 13.092 nm.  L1 retardance mean pi (3141.6 mrad) ideal / 3139.7 mrad real; **retardance VARIATION 0.00 mrad ideal, 0.55 mrad bare Al, 0.53 mrad protected Al**; fringe visibility 0.999 where lit; central-band lit fraction 0.00 / 0.04 / 0.03.  Fold-angle lever: astigmatism scales as fold angle squared; OAP2 9 deg -> 6.4 deg halves it but costs a **~1.33x longer optical leg** (lateral offset drops below the 126 mm the source/camera bodies need).
Seat trim and fold coma [PLAN 7.2; R-IFO 4; `oap_focus_probe.m`, model 512, flat DM]: solved seat trim **MASK_TRIM 6.14 mm** (lens default -5.58 mm), root = OAP1's residual collimation defocus (0.056 deg convergence, focus ~52 m) refocused by OAP2 (F2^2/52 m ~ 6-10 mm); OAP2 fed **exactly on-axis** (0.000 deg); best-focus ray blur vs OAP2 AOI **0.17 / 0.31 / 0.47 / 0.65 / 0.82 lamF/D** at **1 / 3 / 5 / 7 / 9 deg**; with the trim solved at model 1024 **G1 = 2.1e-15, G3 = 4.2e-16, both PASS**, best-focus blur **2.94 um (1.1 lamF/D)**, peak/sum **0.0104** (gate 0.01); the earlier "8 um / 3 lamF/D" was pure defocus (at NA 0.069 a 0.11 mm trim error adds ~7.6 um; peak/sum falls 0.0104 -> ~0.002 within +-0.1 mm).  D1 window placement: **lens 100.00% within 2 px** (median 0.07 px); **OAP 72.77%**, window containment **98.9%** (median 0.16 px).
The other gauges on the OAP front end (`zoap`, bare Al, matrix ON the 30 nm surface) [R-IFO 4; PLAN 7.2]: **L** single 10 nm 1.044 / 22 pm, 47-site grid 0.998 / SNR 5.9, dense 0.985, capture 42 nm; **S** single **1.000 / 4 pm**, grid **0.967 / SNR 221**, dense 0.952, capture **38 nm** (lens-rig S: 0.989 / 5 pm, grid 0.9993, dense 0.984, capture 42 nm); **V G4 FAIL -- 19.6 pm** against the 12 pm gate on the 100 nm pokes; **P G5 FAIL -- 94 pm** against 12 pm.
Figures: `.../tg_psi_dm96_oap/oap_vlayout.png` **3958 x 3250**; `.../lens_vlayout.png` **3958 x 3250**; `.../runs/loop_oap/loop_oap_loop.png` **2026 x 844**.  `.../tg_psi_dm96_oap/runs/oap/d1_picture.png` **2236x1478** (D1 error + column-norm picture, both rigs; R-OAP item 4).
Parts, OAP rig (replaces L1, L2) [R-IFO 6; RM-IFO]: **OAP1** collimator f **857 mm**, AOI **5 deg**, off-axis **149 mm** (+22 mm clearance margin), bare/protected Al; **OAP2** focuser f **429 mm**, AOI **9 deg**, off-axis **132 mm** (+6 mm margin), bare/protected Al.  Off-axis distance = f*|sin(180-2*AOI)|.  BS, DM, reference flat, field lens, camera and polarization optics as the lens rig; the tail is re-tuned for the OAP focuser.

## 13. Systematics priced
**Interferometer** [R-IFO 3; R-OAP item B]: PZT step error 2% / 5% -- single-actuator gain **0.9743 / 0.9543**, floor 2 -> 5 / 9 pm, **negligible in hold** (3 pm at 5.4e12 vs the record 5.5e12, common mode) (`lens_deck_se2/_se5`, `loop_lens_se2`).  Camera 1/f walk (1e-3 of signal, 25% within scan) -- 3 pm at **6.1e12** vs 5.4e12, **+13% light**; exactly immune at cam_intra 0 (`loop_lens_cam`).  BS diattenuation (v1 plate) -- arm rotated **+7.479 deg**, gain **+11.7%**, corrected to **1.00000**; residual 4-theta term **8.9e-4** of the fringe, **1.7e-14 nm** differentially.  MacNeille cube R_p (v2) -- naive odd stack **2.11e-2 (2.1%)**, symmetric period **0**, T_p/T_s **2382:1**, arms 5.3e-6 deg from orthogonal, gain 0.999999.  Coating retardance (OAP fold, L1) -- **VARIATION 0.55 mrad** bare Al against **0.00 mrad** for the ideal reflector (the singular idealization); a floor, not a gain error; no analogue on the lens rig.
**ZWFS scalar** [RM-ZW S7-S10; D-ZW]: model defocus -- the exit sphere carried **23.86 mm** against the entrance sphere's **352.7 mm**, a **4.86 m** defocus of the pupil image; round trip 0.159 -> **1.8e-15**; pupil-brightness modulation under a 30 nm state **29% rms -> 4e-16**; poke kernel peak 0.27 / ring -0.16 -> 0.75 / no ring; the Talbot null near 30 cyc/ap (predicted 34) gone; the correction alone took the single-actuator floor **744 -> 67 pm** (SNR 9 -> 80), the grid case SNR **1.46 -> 14**, dense 42 -> 9.6 nm.  Sampling -- mask px per lam/D = **0.74 * MODEL / NGRID**, detector px per actuator ~ NGRID; 193/1024 gives 7.92 px and 2.5 px/act, 385/1024 gives 3.96 and 5.0, **385/2048 gives 7.92 and 5.0 (the compliant run)**; test-actuator gain 0.946 / 0.996 / 0.996; the stencil-site fix took own-site 0.9576 -> **0.9915**, the record 0.900 -> 0.946 and the 385-ray value 0.935 -> **0.9963**.  Color -- five colors (480 / 532 / 632.8 / 700 / 780 nm) through the one 346.2 nm etch (dimple 2.10 / 1.88 / 1.57 / 1.41 / 1.27 rad; 2.64 / 2.38 / 2.00 / 1.81 / 1.62 lam/D); on the corrected model the combination's minimum transfer is **0.991** against the best single **0.962** -- **3%, not 3x**; rows neutral for I+/S (single-on-base SNR 269 -> 183, grid 31.8 -> 32.4) and **worse for L** (480 nm reads negative, |c| = 1.73); range stays chromatic, **780 nm the best single color**, best pair 700+780 (`rec193full`).
**Vector ZWFS** [RM-ZW V2/V3/V4; PLAN 7.3]: **V2 metasurface retardance error** -- uncalibrated G4 bias on a 12 nm figure, leak in phase **153 / 380 / 750 / 1461 pm** at err 0.02 / 0.05 / 0.10 / 0.20 rad and in quadrature **0.7 / 4.5 / 18 / 75 pm**; a three-number calibration on the flat's two images recovers kappa and eta to five digits and the gate reads **0.048 pm**; rows through the on-surface matrix uncalibrated at err 0.10 and 0.20 are the ideal numbers (**0.9934 / 4 pm**, grid 0.9992 / 3, dense 0.9999 / 331 pm), loop identical (`v2g_*`, `v2loop`).  **V3 arm polarization aberration** -- lens-rig diattenuation **5.1e-3 mean, 1.1e-3 rms, 7.8e-3 max**, retardance **0.9 mrad mean, 0.7 rms, 2.8 max**, channel phase difference **1.63 mrad rms (PV 8.0)** with the laser at 45 deg, 2.31 at 0 deg, **2.1e-5 at 90 deg**; uncalibrated absolute G4 **9.0 pm at 45 deg** (17.2 at 0, 11.2 at 90; ideal 0.053; AR-coated 5.0); `'amp'` calibration 0.11 pm at 90 deg, 9.5 at 45; design scan G4 **59 / 178 / 597 / 1877 pm** at 0.01 / 0.03 / 0.1 / 0.3 rad of channel phase (**6 nm per rad**) and **226 / 677 / 2268 / 7106 pm** at the same amplitude ratios (**3.8x per unit**); through the on-surface matrix rows hold to 0.1 and bend at 0.3; loop at 0.3 rad contraction 0.554 against the ideal 0.509 (**10% gain loss**), steps -> 0.000 pm, 3 pm held below 1e13 (`v3arm*`, `v3s_*`, `v3loop`).  **V4 analyzer (QWP + cube)** -- cube extinction leak l_A **4.2e-4** uniform, l_B **6.3e-4** varying to **2.6e-3** with the cone's AOI; coherent term |c| = delta/2 for a retardance error delta and = theta for an azimuth error theta; an ideal plate in the converging beam leaves **1.5e-3 of zero mean**:

| analyzer | uncalibrated G4 (100 nm pokes) | single 10 nm on the 30 nm surface | grid 1 nm |
|---|---|---|---|
| ideal | **0.30 pm** | 0.9942 / 5 pm | 0.9961 / 2 pm |
| cube alone | **3.6 pm** | -- | -- |
| ideal-plate term, bound | **23 pm** | 0.9938 / 5 | 0.9963 / 2 |
| plate lambda/300 | **162 pm** (1.4% of the figure) | 0.9919 / 6 | 0.9975 / 2 |
| plate 1 deg azimuth | **264 pm** | -- | -- |
| plate lambda/100 | **487 pm** | -- | -- |

Linear at **15 pm per unit of |c|**; |c| = 1.05e-2 at lambda/300, 3.1e-2 at lambda/100, 1.7e-2 at 1 deg.  Spec: zero-order plate at **lambda/300**, azimuth to **1 deg**.  Dev resolution 512 / 65; record 1024 / 193 = `an193_*` (`v4seq.sh`).  `'none'` reproduces the record bit-for-bit (an_ref = uwoff_ref).
**PDI** [R-PDI 4b, 5, Concl. 6/7; RM-PDI]: step error 2% -- four-step least squares turns a 12 nm figure into **421 / 251 pm** (P / PF) and the flat reads 6.9e-2 / 5.7e-2 rad rms instead of 1e-15, while the **five-frame Schwider-Hariharan** scan reads **4.9 / 2.2 pm** with the flat at 1.3e-5 rad rms and the differential rows error-free to the digit (`pdi193se_ls`, `pdi193se_sh5`).  Camera drift -- the paper's 0.13 e/px/cycle is **invisible to every reading** (`pcam193`: L 1.39 / 0.14 pm at 1e13 / 1e15, S 1.52 / 0.15, V 1.17 / 0.12, P 1.45 / 0.14, PF 2.50 / 0.25, with and without), because at 1e13 photons per measurement a lit pixel collects **~3e8 photons per frame** so an electron is **1e-4** of its shot noise; in the relative form L imprints **10.8 nm** and V **89 pm** while S/P/PF read their noise-only values to the digit (`pcam193r`); with the whole step inside each scan S **5.4**, P **5.3**, PF **10.7 pm** at 1e15 (`pcam193ri`).  Reference-arm walk (PF only) -- **4.8e13 / 4.9e13 / 5.0e13** photons/cycle for 1e-3 / 1e-2 / 1e-1 rad/cycle, **4% over a hundredfold range**, P bit-identical (no such arm); the mechanism is that a path-length change is a **piston**, the one mode the actuator estimator nulls (a rank-one term since S10); NOT benign for absolute E-field reconstruction, and a reference arm that **tilts** is outside the model (`rw193_*`).  Reference MOTION -- traced arm moving with the state **5.889 pm = 0.045%** of a 12 983 pm figure, traced once and frozen **0.000 pm**, synthesized LP01 **0.000 pm**; shape change under the 30 nm surface traced **0.0103 (1.03%)**, FFT pinhole surrogate at 2 lam/D 0.0059, synthesized **0 by construction**; differentially the gain is identical to 3 digits with **SNR -12%** and dense floor **+12%**, range unchanged; in loop **no fixed error**, both noiseless steps 0.000 pm (`pfdeck`, `pfdeck_frz`, `pfdeck_loop`).
Figure: `demo_session/figs/pdi193fbase_pdi.png` **1681 x 980** (focal spot with pinhole / dimple / mode; reference amplitudes; reference motion by diameter; visibility maps).

## 14. Recommendation -- each report's opening
**R-IFO, "The interferometer's role in the deck: CAPTURE":** "The interferometer's edge is CAPTURE, not hold" -- every start **60-300 nm surface (120-600 nm WFE)** converges to the **~2 pm** hold floor in tens of cycles with no recalibration, while the focal-plane sensors' capture dies by **~160 nm** of surface.  "In HOLD the sensors beat it": S and V hold 3 pm at **1.5-2.6e12** photons per measurement against the interferometer's **5.5e12**, and with **no fixed noiseless-step error** against the four-step's **~9%**.  "So the deck's line is the **HYBRID BENCH**: capture with the interferometer, hold with the sensor."  "**Lens, not OAP**" -- the reflective fold leaves **~0.18** residual cross-talk against the lens's **< 0.06**, which walls both the hold and the capture.  Best *internal* IFO form: the **polarization snapshot calibrated by a PZT four-step**.
**R-PDI, section 0:** "The point-diffraction approach's best configuration is the **STEPPED PINHOLE P**, with a pinhole-only (shutter) frame per state and the **five-frame Schwider-Hariharan scan**" -- six frames, all common path, one plate in the mask seat.  The numbers: **0.9935 / 4 pm / SNR 2790**; N(1 pm) **3.3e13** flat matrix, **9.8e13** on-surface (level with S at 9.3e13) and **3.8x** cheaper than the P/SRI's 3.7e14; 3 pm held from **2.3e12** noise-only and **7.0e12** under the walk; steps to **0.000 pm**; capture **1.02 / 1.06 / 1.13 at 120 / 240 / 480 nm**; 2% step error **4.9 pm**; capture to **100 nm (200 nm WFE)** with unwrapping and re-calibration.  The P/SRI's "place in the deck is as the **CAPTURE** instrument, not the hold instrument".
**D-ZW conclusions:** "**The polarized dimple is the best reading on every line**: two simultaneous frames, no fold, 3 pm held from **5.3e12** photons per cycle, and a calibration that stays within 2% out to 60 nm rms."  "Calibrate with a measured response matrix, on the working surface": single actuator on the flat **0.994 / 4 pm**; on the working surface the stepped reading **0.99 / 5 pm** and the linear one-frame reading usable (**1.05 / 23 pm**).  "In closed loop the stepped reading has no fixed error."
**PLAN 7.1 consequence:** a self-referenced reading (Zernike dimple scalar or vector; stepped pinhole) captures to **~60 nm** of surface (120 nm WFE); an externally referenced one (interferometer flat, P/SRI waveguide) captures the whole **100-200 nm WFE** range with unwrapping.  So the bench needs the hybrid, the P/SRI alone, or the sensor with a second color.

## 15. The modes (PLAN 11.1)
- **Mode 0, ground flat**: the ground-calibrated voltage map (~0 WFE on the ground); on orbit the residual is launch, gravity release and thermal change, **100-200 nm WFE**.
- **Mode 1, image-based phase retrieval**: the WFS&C loop's own focal-plane retrieval, no gauge.  Wraps at the **half-wave** level with high-spatial-frequency WFE -- the same failure class as the gauges' wrap, **one wave earlier**.  Its reach sets what mode 2 must capture.
- **Mode 2, capture**: an externally referenced reading (interferometer or P/SRI) with unwrapping, re-calibrated on the surface as it moves; or the sensor with a second color.  Ends at the hold regime (**~30 nm** of surface, then the matrix measured there).
- **Mode 3, closed-loop hold**: the sensor at picometers (stepped or vector Zernike, or the pinhole), **gain 0.5**, the matrix re-measured on the held surface when the calibration ages.
- **Recalibration events** between modes: the matrix on the current surface (photon cost **5-40x** of a measurement); the flat re-taken.
- Capture limits on the arrows: **30 / 60 / 100 / 150 nm** of surface (PLAN 7.1).
- Figure: `demo_session/figs/gauge_modes_flow.png` **1800 x 1000** (MATLAB, 16-19 pt; graphviz is not installed on the build box).  Tool `demo_session/gauge_modes_flow.m`.

## 16. The complex amplitude (PLAN 11.2)
Which readings give it: ZWFS linear L and exact I **no** (one frame; the solve assumes the flat's amplitude); ZWFS stepped S **yes** (its clear frame IS the pupil intensity, the three depths give E conj(b)); vector V **yes with the clear frame** (see below); pinhole P and P/SRI PF **yes** (the phase-stepped solve is the complex field against a known reference, Dube 2024); interferometer four-step **yes** (fringe modulation = test amplitude x the known reference's).
V5 measured (`an193_clear`, `v5seq.sh`) [PLAN 11.2; RM-ZW V5]: the pair ALONE is ambiguous -- its two images are two circles whose intersections are mirror images about the reference wave, so it measures **A sin(phi)** and **|A cos(phi) - b|**; 100 nm pokes reach **1.9 rad** and cross it; the iterated solve diverges (**10% in amplitude after five passes**) and reads **5.8 nm** on the G4 pokes; exact on the flat (**2e-15**), and a 0.3 rad blob with the true reference wave reads **3e-4 rad / 8e-4 amplitude**.  The pair **plus the state's CLEAR FRAME** reads both exactly: through a **5% / 20%** pupil amplitude dip **0.35 / 0.50 pm**, against **241 / 967 pm** for the phase-only solve of record, with the 30 nm-surface rows unchanged to the digit (**0.9942 / 5 pm**, **0.9961 / 2 pm**).  Deck statement: the vector sensor gives the complex pupil field from **three frames** (the two images + the clear frame), the stepped sensor from its four steps plus the clear frame, the pinhole and the interferometer from their steps -- **no separate camera**.  Runner `mask.v_clear`, gate G9 `mask.v_dip`.  Not wired: the noise and loop stages' frames (they read the pair only).

## 17. One bench, all modes (PLAN 11.3)
| switch | part | modes |
|---|---|---|
| the mask seat translates | one substrate with the etched dimples, the pinholes (with their attenuated surrounds), the metasurface and a clear window (the VSG2 nine-spot idea) | Zernike / vector Zernike / pinhole / clear (interferometer, phase retrieval) |
| reference-arm shutter | the interferometer's reference flat on its PZT stays built | interferometer on / sensors (arm shuttered) |
| quarter-wave plate in or out | the MacNeille cube and camera B stay behind the field lens; with the laser p-polarized the cube transmits (**98%**) to camera A in every scalar mode; the plate in makes the vector split | vector Zernike / all others |
| flip-in pickoff plate | the P/SRI's second arm (Lr1, pinhole, Lr2, folds, compensator, waveguide, BS3) on its own breadboard behind a flip-in plate | P/SRI / all others |
| in-arm quarter-wave plates in or out; analyzer | the polarization-snapshot form on the plate rig; the PZT form needs neither; the hybrid uses both | interferometer forms |

The **v2 cemented-cube** interferometer is the one form that does not switch in (it replaces the plate splitter) -- backup slide.
Figure: `demo_session/figs/gauge_one_bench.png` **1800 x 1000** (tool `demo_session/gauge_one_bench.m`, 2026-09-14 19:00).  PLAN 11.3 assigns the universal-bench drawing to TO (`psri_bench` + `zwfs_vlayout` recipes combined); the existing file is the CCL-side version -- QA before use.

## 18. Future work -- a second color (PLAN 10.1)
- Unwrapping resolves a wrapped differential only while neighboring pixels differ by **less than pi**: at **4 px per actuator** a 100-200 nm rms actuator pattern is **0.7-1.4 rad per pixel** and unwraps; a figure of **more than a wave** with actuator-scale structure (quilting, print-through, a stuck actuator) does not.
- Synthetic wavelength lambda1*lambda2/|lambda1-lambda2|: **632.8 + 700 nm = 6.6 um**, unambiguous surface range **1.6 um** (double pass); **632.8 + 780 nm = 3.35 um**, range **0.84 um**.
- Noise cost **lambda_syn / lambda = 5-10x** -- the coarse color pair captures, the single color holds.
- Machinery exists: the color stage runs one physical mask at five wavelengths (phase **1.57 rad at 632.8 nm**, **1.42 at 700**; stepped depths scale as (n-1)/lambda); a two-color IFO is the classic form.
- Measurement to add: the start-rms ladder with a two-color coarse solve feeding the single-color loop, per reading; the number is the largest capturable start and its light.

## 19. Future work -- other approaches (PLAN 10.2; one line each, none modeled)
- **Shack-Hartmann or modulated pyramid** as the capture stage: range of many waves, no wrap, sensitivity far from picometers; hands off to the gauge once inside its range.
- **Phase diversity** (two defocused pupil images, no mask): wide range, iterative, common path; a candidate for capture with the existing camera and a translation stage.
- **White-light (low-coherence) scanning** in the interferometer: absolute surface with no wrap, slow; a one-time capture tool.
- **Model-based large-figure solve**: the DM's influence functions as the basis of a nonlinear (iterative) fit to the sensor's frames; the actuator-space estimator already carries the basis, the nonlinearity is the addition.
- **Vector Zernike with a polarization camera** (micro-polarizer array): one camera instead of the cube and two; costs a quarter of the pixels per image.
- **Heterodyne or lock-in detection**: immunity to slow drift and 1/f electronics, at the cost of a frequency-shifted reference (the interferometer and the P/SRI can carry it; the common-path sensors cannot).
- **Direct actuator metrology** (capacitive, optical) as the coarse reference the optical gauge is calibrated against.

## 20. Future work -- as-built (PLAN 10.3)
Model today: ideal optics, a perfect camera, photon noise, DM and camera drift.  As-built adds, in the order they are likely to matter:
1. **The camera's throughput sets the measurement time, not the laser.**  **1e14** photons per measurement over **~1e5** pupil pixels = **1e9 electrons per pixel**; a **1e5**-electron well means **1e4 co-added frames** -- at **100 frames per second, 100 s per measurement**, against the **0.13 s** the laser needs.  Price: well depth, frame rate, read noise per frame (which enters as sqrt(frames) x read noise), bit depth (quantization at **1e-4** of the signal), gain nonlinearity and pixel-response nonuniformity, persistence between frames.
2. **Optical surface errors:** each lens, OAP, plate and the mask substrate with a typical figure (**lambda/10 to lambda/20 PV, a 1/f^2 spectrum**) as GridData on the element; the sensor's reference core sees the low orders (the S7 lesson: **0.16 of defocus moved every number**); the mask substrate's flatness and the dimple's etch-depth error (the metasurface's is V2).
3. **Alignment and stability:** the sensitivities exist (D4: **10 um, 10 urad**); add mask centering drift on the focal spot, thermal expansion of the bench (a leg length per degree), source pointing and wavelength drift (mask phase goes as 1/lambda), laser polarization angle (V3), and for the two-arm instruments the non-common-path air and mount drift.
4. **More drift terms in the loop:** actuator hysteresis and creep, command quantization (**14-16 bit: 0.1-0.5 nm steps** as a floor), influence-function error (absorbed by a matrix measured through the sensor), vibration within a stepped scan (`loop.intra`).
5. **The error budget in the JPL form:** per reading, fixed terms (calibratable, the residual after calibration), drift terms weighted by the servo bandwidth, and noise terms (photon, read, quantization), summed in quadrature to a held error at a stated light and time.  The deck's comparison table becomes its first column.

## 21. Run it yourself
Interferometer [RM-IFO; R-IFO "Reproduce"] -- sheet/runner `tg96_params.m` / `tg96_run.m`, model 1024 ~11 GB, one at a time:

    tg96_run                                             % lens rig
    tg96_run('bench.optics','oap','tag','oap')           % reflective rig
    ./tg96_batch.sh lens_deck "'stages',{'bench','deck'},'battery.noise',true"
    ./tg96_batch.sh oap_deck  "'bench.optics','oap','bench.coat_oap','bareAl','stages',{'bench','deck'},'battery.noise',true"
    ./tg96_batch.sh loop_lens "'stages',{'bench','loop','figs'}"
    matlab -batch "tg96_tail('tag','oap','bench.optics','oap')"   % tail retune, OAP first

Zernike sensor [RM-ZW; D-ZW] -- sheet/runner `zwfs_params.m` / `zwfs_run.m` (+ `zwfs_run_figs.m`, `zwfs_run_batch.m`, `zwfs_batch.sh`, `macos_param_2048.txt`); equivalence gate `runs/rec193` vs `zwfs_s7iter_report.txt`, 64 row/ladder lines, **8 differ in the last printed digit** (pcg tolerance), every gate identical:

    out = zwfs_run;                                      % bench + battery + figs
    out = zwfs_run('tag','ng385','NGRID',385);           % 1 Mpix-class camera
    zwfs_run('MODEL',2048,'NGRID',385,'param_file','macos_param_2048.txt')
    zwfs_run('stages',{'battery','color','noise','figs'})
    ./zwfs_batch.sh loop193 "'stages',{'bench','loop','figs'}"    % ~80 min

Point-diffraction [RM-PDI] -- sheet/runner `pdi_params.m` / `pdi_run.m`; code shared with `zwfs_run` + `../dm_gauge_lib` (nothing copied); `pdi_batch.sh` is serialized with `zwfs_batch.sh` on the same lock, **one engine MATLAB at a time on the box**:

    P = pdi_params;  out = pdi_run(P);                   % the record sheet
    pdi_run('pdi.DIA_LAMD', 1.0, 'stages', {'bench','battery','figs'})
    ./pdi_batch.sh TAG "pdi_params, 'stages',{'bench','loop','figs'}"

The other gauges on the OAP front end [R-IFO "Reproduce"]:

    zwfs_run('bench.optics','oap','bench.coat_oap','bareAl','readings',{'L','S','V','P','PF'}, ...
             'stages',{'bench','battery','noise','loop'},'mask.v_arm','engine')

## Backup
**B1. The three IFO phase-shift forms** [R-IFO 3]: *PZT four-step* -- reference-arm phase stepped in time (sequential); gets wrong phase-step miscalibration and within-scan camera / DM drift; priced by `pzt.step_err` 2% / 5% and `dmg_loop` cam.  *Polarization snapshot* -- four analyzer channels at once (simultaneous); gets wrong the polarization systematics, no within-scan drift; v1 plate **+11.7%** (correctable to 1.00000), v2 cube R_p, coating retardance.  *Hybrid* -- snapshot for the change, PZT for the absolute step; each half removes the other's error; the recommended form.  The snapshot's weakness is that its four analyzer azimuths are fixed design constants, so a systematic in their realization is a fixed gain it cannot self-calibrate; the PZT's is the sequential scan.  R-IFO: "the four-step is robust to all three sequential-form systematics in hold mode ... the hybrid's value is mainly the absolute calibration, not drift immunity."
**B2. OAP rig loop** -- slide 12, plus [R-OAP D7]: OAP sig_n **22.2 pm at 1e12** -> 1 pm at **~4.9e14** (3.5x the lens); the walk's cross-talk bias is **3.25 pm at 1e15**, adding in quadrature to the 2.4 pm walk floor to give the 4.1 pm floor.  The lens rig is ~**2.6x** the ZWFS photons per measurement (5.5e12 vs 2.1e12 for L noise-only; 2.0e13 vs 7.5e12 for S under the walk).
**B3. The scalar Zernike readings that lost, and the fold** [RM-ZW S7/S10/S11]: **L** holds walk and noise like S (2.1e12 / 7.3e12) but thermal **27.6 pm** and still creeping at cycle 60 (9 pm lag + 26 pm of >12 cyc/ap error the loop imprints; 15 pm with the worst 100 actuators removed; two 4-actuator clusters at 0.5-0.8 nm); its noiseless steps decay with rho 0.82-0.83 (1.2 pm left at cycle 60 from 1 nm, 9.8 from 10 nm).  **I+ DIVERGES** in the loop -- a 1 nm step grows to **99 nm** in 60 cycles, 10 nm to **178 nm**, and the noise-only loop at 1e15 wanders to 0.9 nm; actuators whose footprint sits beyond the quarter-wave fold read with the wrong sign, and `P.loop.rmax` now stops such runs at 1 um.  **The fold**: on the 30 nm surface **7.8%** of pupil pixels sit past it; the true beyond-fold fraction is **7.7 / 12.6 / 16.9 / 20.6%** at 30 / 40 / 50 / 60 nm rms, which the refined prior finds to 4 digits (the plain one finds 4.7 / 7.2 / 8.6 / 8.7%) and gets **99.99%** of pixels right against the plain stepped map's 3% miss; with the oracle b and the true branch the inversion is exact to **2.8e-14**.  **Fold crossings under a change** (`fold_diag`): 3 of 2258 beyond-fold pixels move under a single 10 nm change (0.01% of the mask), 10 under the 47-site 1 nm grid, **946 (3.3%)** under a dense 10 nm random change.
**B4. Pinhole-diameter and P/SRI-vs-pinhole trades** -- slide 7, plus [R-PDI 3b]: the range is the DIAMETER's, not the sampling's -- `pdi193d1` runs 1.0 lam/D at **1024 / 193** (the same sampling as the 2.0 lam/D leg) and also does not fold (**1.02 / 1.06 / 1.13 at 120 / 240 / 480 nm**); that run's t_auto 0.28, eta_pin 0.24, throughput 0.29, visibility 0.94, G5 **0.002 pm**, N(1 pm) at the camera **2.9e14** against S 1.0e14.  Run-record note: both model-2048 runs ended **exit 137 -- NOT a failure**, the documented model-2048 crash AT EXIT after everything is written; both reports end "run complete" (123.6 and 69.5 min).
**B5. Drift -- camera, DM, within-measurement**: the paper's camera number (`pcam193`, 0.13 e/px/cycle = ~1 e over the run) is invisible to every reading; the relative form (`pcam193r`, 1e-3 of the mean photons per lit pixel per frame, one scale per scan) gives L **10.8 nm** and V **89 pm** with S/P/PF at their noise-only values; the whole step within each scan (`pcam193ri`) gives S **5.4**, P **5.3**, PF **10.7 pm** at 1e15.  DM drift within a scan has the opposite sign (slide 11).  I+ on this base floors at **885 pm** regardless (its fold-flipped sites).
**B6. V2 / V3 numbers** -- slide 13.
**B7. The model correction and sampling story** -- slide 13, plus [RM-ZW S7 Block A]: legacy unmasked round trip **0.159** on a flat DM; SPH2PL focal factor S = **8.60e-5 rad/px^2**; z_eff = Z1(Z1-Z2)/Z2 = **4860 mm**; Talbot null predicted at **34.1 cyc/ap** against the S3 record's 28-32; on the 30 nm base |E|/|E_flat| std **0.290** (range 0.017..2.38).  Corrected: round trip **1.8e-15**, b surrogate against the engine's Eb **2.0e-15**, amplitude modulation **4e-16**.
**B8. The coating story (item B)** -- slide 12, "Coatings".
**B9. Provenance -- run tags by directory.**  `.../tg_psi_dm96_oap/runs/`: descent_lens, descent_oap, lens, lens_base, lens_deck, lens_deck_se2, lens_deck_se5, lens_render, lens_sketch, loop_lens, loop_lens_cam, loop_lens_intra, loop_lens_se2, loop_oap, oap, oap_bareAl, oap_base, oap_coat, oap_deck, oap_jones, oap_render, oap_sketch, oap_vor, zoap.  `.../zwfs_dm96/runs/` (ZWFS plus the pre-2026-09-13 PDI records): an193_*, cal_base, cap385, cap385_flat, cap385_b60/b90/b120/b160, fold_diag, ks_hold*, loop193, loop385, m2048, m2048_lat, mapdiag, mask193, mask385, mat193, mat193b, mat193c, mat385, matbase, matbase385, ng385, ng385g, ng385_lat, ng385s3, noise193_b30/b60/b120/b160, pcam193, pcam193i, pcam193r, pcam193r_perframe, pcam193ri, pdi193, pdi193d1, pdi193f, pdi193fbase, pdi193se_ls, pdi193se_sh5, pdi193state, ploop193, rec193, rec193full, rec193_lat, v193base, v193flat, v193noise, v2g_*, v2e10a, v2e10afit, v2e20a, v2loop, v3arm*, v3s_*, v3loop, vlad193, vloop193.  `.../pdi_dm96/runs/` (everything from 2026-09-13): cap385p, cap385p_b60/b90/b120/b160, cap_nouw, cap_nouw_recal, cap_uw, cap_uw_recal, cap_state_nouw, cap_state_uw, cap_state_uw_recal, descent193, descent193f, descent193s, intra193, intra193_0, noise193p_b30/b60/b120/b160, pfdeck, pfdeck_frz, pfdeck_loop, pfdeck_smoke, pfdeck_smoke2, pfsmoke_ref, pin10_2048, pin10_loop, pin20_1024, pin20_loop, rw193_1e1, rw193_1e2, rw193_1e3, sm_cap, sm_cap2, sm_knobs, sm_psri_nl, uwoff_bat, uwoff_ref.  Sequence drivers gsmoke.sh -> gseq1.sh -> gcap.sh -> gseq3.sh -> gseq4.sh -> gseq2.sh (gmaster.sh, gmaster3.sh, gmaster4.sh; gmaster2.sh superseded).  Gates: tDmgLoop **15/15**, mmacos fast suite **469 pass / 0 fail**, tBench 9/9 (`'lens'` byte-identical, no engine work).
**B10. Figure inventory -- all ls-checked.**  In `/home/dcr/dev/macos/demo_session/figs/`: zwfs_render_rig.png 2657x875 (shared front end); zwfs_mask385_mask.png 1230x515 (mask + focal spot); zwfs_response.png 1443x604 (first light, legacy model); zwfs_s7iter.png 1593x699 (model correction); zwfs_m2048_battery_n96.png 1588x699 (fully sampled battery); crop_zwfs_rec193full_color_L.png 618x606 (color); zwfs_rec193full_noise.png 749x593 (photons); crop_zwfs_loop193_left.png 796x677 (loop, S); crop_zwfs_vloop193_right.png 808x681 (loop, S vs V); zwfs_vlayout.png 2438x1368 (vector layout, both channels); crop_zwfs_vlayout_tail.png 1808x1020 (vector tail crop); zwfs_poke_triptych.png 1858x554 and zwfs_defocus_triptych.png 1864x554 (backup, legacy model); pdi193fbase_battery_n96.png 1588x699; pdi193fbase_pdi.png 1681x980; gauge_modes_flow.png 1800x1000 (slide 15); gauge_one_bench.png 1800x1000 (slide 17, QA before use).  In the template tree: tg_psi_dm96_oap/lens_vlayout.png 3958x3250 and oap_vlayout.png 3958x3250; tg_psi_dm96_oap/runs/loop_lens/loop_lens_loop.png 2026x844 and runs/loop_oap/loop_oap_loop.png 2026x844; pdi_dm96/pdi_layout.png 2438x1368, psri_layout.png 2438x1368, psri_render.png 2438x1030 (all **current**); zwfs_dm96/runs/loop193/loop193_loop.png 1609x699 and runs/vloop193/vloop193_loop.png 1609x699; zwfs_dm96/zwfs_layout.png 1467x539 and zwfs_mask_fig.png 1387x604 (older); tg_psi_dm96_oap/runs/oap/d1_picture.png 2236x1478.
**MISSING / STALE:** no redone **front-end-only** figure in the recipe (PLAN 4 lists it as owed by CCL).  `demo_session/figs/pdi_layout.png` (3662x1310), `psri_layout.png` (3232x1840) and `psri_render.png` (3125x1094) are the **pre-2026-09-13** `Bench.sketch` versions with E-number labels at 9-11 pt and differ by md5 from the current `pdi_dm96/` files -- use the `pdi_dm96/` copies.  `demo_session/figs/pdi_layout_tail.png` (3068x1624) exists but the file was **retired** in `pdi_dm96/`.
