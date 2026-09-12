<!--
deck_pdi.md — point-diffraction interferometer readings on the DM gauge
bench: the stepped pinhole (common path) and the phase-shifting
self-referenced interferometer of Dube et al. 2024 (waveguide reference,
photonic phase shifter), beside the Zernike sensor's stepped and
polarized-dimple readings on the same DM truth.  DRAFT — pending Dave's
sign-off; the builder suppresses the export mark on DRAFT decks.
Build: python3 make_brief_slides.py deck_pdi.md
Template: deck_zwfs.md / deck_tg_fang.md (register per
MACOS_resources/doc/STYLE_REPORTS.md; figures = the runner's own PNGs,
unmodified).
2026-09-12 v1: literature and model, dev gates, record 1 (flat matrix,
runs/pdi193), record 2 (fiber reference, runs/pdi193f/fbase/state/d1/se),
the camera-drift loop (runs/pcam193, pcam193i).  Numbers marked [R2] /
[CAM] are filled from those runs' reports as they land.
Sources: templates/40_benches/zwfs_dm96 (README "P / PF", runs/pdi*),
dm_gauge_lib/dmg_pdi_gauge.m, macos/BRIEF_pdi_campaign.md,
~/dev/MACOS_sandbox/pSRI/Dube_pSRI.pdf + photonic-sri-simulation/.
-->

# A point-diffraction interferometer for the DM gauge
The same 96 mm bench and the same deformable mirror: a stepped pinhole at the focus, and the phase-shifting self-referenced interferometer with a waveguide reference (Dube et al. 2024), read exactly and compared with the Zernike sensor's stepped and polarized-dimple readings
D. C. Redding, with Claude Code.
September 2026.
DRAFT — pending review.  Status: both readings built into the sensor runner and gated; the flat-matrix record measured; the fiber-reference record, the 1 λ/D pinhole, the step-scheme trade and the closed-loop rows with a camera-drift knob in progress.

## The idea: a reference made from the beam's own focal core, kept clean | A pinhole passes only the core of the spot; what diffracts from it is a smooth wave whose shape barely depends on the aberration that made the spot
::: left
- **Common-path form (P):** a pinhole in a partly transmitting plate at the focus.  The core passes the pinhole and spreads over the pupil image as the reference; the rest of the beam passes the plate, attenuated so the two amplitudes match.  Stepping the pinhole's phase gives four frames and a textbook four-step solve: exact, linear, amplitude and phase, no small-phase limit, no sign fold.
- **The Zernike sensor is this instrument with a phase dimple and a clear plate.**  At unit plate transmission and the dimple's diameter, the stepped pinhole reading reproduces the sensor's stepped reading to 5×10⁻¹⁵ (check G7).  What the pinhole changes is the reference: smaller and attenuated, so it moves less with the surface.
- **The price is light:** the plate transmission squared, 0.72² here; the sensor keeps all of it.
::: right
![The focal spot on the flat DM with the 2 λ/D pinhole (red), the sensor's dimple (green, the same size) and the waveguide mode's 1/e contour (brown); the reference amplitudes across the pupil image; how much each reference moves under the 30 nm working surface, by pinhole size; the fringe visibility of both readings on the flat.  Dev-resolution run.](figs/pdi_dev3_pdi.png){h=3.6}
~ Smartt & Steel 1975 (PDI); Medecki, Tejnil, Goldberg, Bokor 1996 (phase-shifting PDI); Naulleau 1999 (EUV PS/PDI, reference accurate to λ/330).

## The P/SRI: the reference from a single-mode waveguide, stepped in the chip | Dube et al. 2024 (SPIE 13092-178): a non-common-path interferometer that a coronagraph layout can carry, measuring both phase and amplitude and high-pass filtering the camera's slow drift
::: left
- **Layout:** a plate beamsplitter sends 60% of the light to a reference arm that focuses onto a single-mode waveguide in a photonic chip; only the waveguide's mode comes back out — the point-diffraction reference — phase-shifted thermo-optically, recollimated by two OAPs and recombined with the 40% test beam.  Out of band, with a laser and a notch filter in the science channel.
- **Reconstruction:** the five-frame Schwider–Hariharan scan; sums with weights that add to zero give the field's cosine and sine per pixel, hence amplitude and phase.  The change between two measurements is taken by complex division, so a 20 pm change on top of a wrapped absolute phase comes back exact (their example: 10⁻¹⁴ nm).
- **Why phase shifting:** any camera pattern constant within one scan subtracts; Roman's Zernike sensor is limited by ~1 electron per pixel of camera drift over 12 hours.
- **Modeled here as their code has it:** the exact step-index mode (V 2.3, b 0.5, core radius 0.5 λ/D), the coupling of each state's focal field into it, the 60/40 split, the four-step and the five-frame schemes, an optional shutter frame of the reference alone.
::: right
- **What moves with the surface, and what does not.**  The waveguide fixes the reference's shape; the surface changes only how much light couples in — a scalar.  Under the 30 nm working surface the coupling drops to 0.86 with 0.025 rad of phase (record run); the differential phase divides the scalar out.
- **The pinhole reference by size (record 1, 30 nm surface):** shape change 0.10% at 0.5 λ/D, 0.26% at 1, 0.58% at 1.5, 1.06% at 2 (the dimple), 2.6% at 3; the amplitude scale 0.86 at every size — the Strehl ratio.
- **Coupling and light:** 59% of the focal light couples into the mode on the flat DM; with the 60/40 split the reference is budget-limited to 73% of the visibility-1 amplitude; visibility 0.86, throughput 0.75.
~ Their next steps, stated in the paper: noise, non-common-path drift, detector systematics.  The model here is ideal in the arm; the camera-drift knob (slide 7) is the first of those.

## The setup and the scoring: nothing changes but the mask and the solve | The sensor's test arm, the same tuned tail, the same DM truth, the same actuator-space score and the same closed loop
::: full
- **Bench:** the TG96 test arm alone with the mask seat inside the reference-sphere bracket at the internal focus; the camera at the DM image.  The pinhole plate replaces the dimple plate; the P/SRI adds nothing to the bench — its reference is synthesized from the traced pupil field, exactly as the paper's model does.
- **Checks that run every time (all pass at 10⁻¹⁵ to 10⁻¹⁰):** the reference model against the engine's own field; the flat DM reads zero; the 100 nm sparse pokes beyond the sensor's fold come back to 0.3 pm (pinhole) and 0.000 pm (waveguide) against the engine's phase; the pinhole reading at unit transmission equals the stepped Zernike reading.
- **Scored as the sensor is:** a change of known actuators on a working surface, before and after, fitted through the measured response matrix (every actuator poked once in sparse grids) in actuator units; gain, floor (pm), SNR; photons per measurement counted at the camera, the throughput printed beside them.
| reading | frames per measurement | reference | light kept |
| Zernike stepped (S) | 4 | the dimple's core, moves with the surface | all |
| Zernike polarized pair (V) | 2 at once | the same, re-computed | all |
| stepped pinhole (P) | 4 | pinhole-diffracted, iterated | 0.82 |
| P/SRI waveguide (PF) | 4 (or 5) | the waveguide mode, fixed | 0.75 |
~ Photons "per measurement": one DM shape measured once, all frames summed.  Sampling: the pinhole gets the dimple's rule (≥ 6 pixels across at the focus; 3.96 per λ/D at 1024/193, so 1 λ/D needs the 2048 grid).

## Record 1: on the working surface the three exact readings agree, and the stepped sensor does not | A 10 nm actuator change on a 30 nm surface reads at 0.93–0.94 for the polarized pair, the pinhole and the waveguide; 0.75 for the stepped dimple (flat-DM calibration)
::: full
| 96×96, calibration on the flat DM (runs/pdi193) | S | V | P | PF (pinhole-shaped reference) |
| 20 nm on the flat: gain / floor pm | 0.996 / 4 | 0.994 / 4 | 0.994 / 4 | 0.994 / 4 |
| 10 nm on the 30 nm surface: gain / error pm / SNR | 0.752 / 47 / 444 | 0.939 / 25 / 419 | 0.933 / 25 / 420 | 0.940 / 25 / 418 |
| 1 nm on 47 actuators, same surface: gain / SNR | 0.829 / 82 | 0.997 / 83 | 0.992 / 84 | 0.998 / 83 |
| 10 nm random on the same surface: gain / error pm | 0.818 / 2119 | 0.999 / 1072 | 0.993 / 1069 | 1.001 / 1073 |
| working surface at which the gain drops below 0.8 | 40 nm | 120 nm | 120 nm | > 120 nm (0.77 at 120, 0.55 at 240) |
| photons per measurement for 1 pm (camera) | 5.4×10¹³ | 4.7×10¹³ | 3.3×10¹³ | 2.1×10¹⁴ |
- **The pinhole and the waveguide read as the polarized pair does** — exact solves — at four sequential frames instead of two simultaneous ones.
- **Per photon the stepped pinhole is the cheapest of the four** (1.6× fewer than the stepped dimple): its reference carries more of the modulation.  The apodized first reference of the waveguide reading cost 4×; record 2 re-prices it with the mode's own shape.
- **Range:** the pinhole's reference amplitude collapses with the Strehl ratio like the dimple's, so both fold at 120 nm; a reference that does not depend on the surface keeps reading.
~ Flat-DM calibration (the record's first form); the on-surface calibration that took the stepped sensor to 0.99 / 5 pm is record 2's.  Errors in pm over the lit actuators; SNR = recovered change / floor.

## Record 2: calibrated on the working surface, the three exact readings are the same instrument | A 10 nm change on the 30 nm surface reads 0.9935 with a 4 pm floor for the polarized pair, the pinhole and the waveguide; the waveguide reading has no range limit
::: left
| 96×96, response matrix measured on the 30 nm surface (runs/pdi193fbase) | S | V | P | PF |
| 10 nm on the surface: gain / floor pm / SNR | 0.989 / 5 / 2120 | 0.994 / 4 / 2835 | 0.994 / 4 / 2790 | 0.994 / 4 / 2842 |
| 1 nm on 47 actuators: gain / floor pm | 0.999 / 4 | 0.999 / 3 | 0.999 / 3 | 0.999 / 3 |
| 10 nm random: gain / error pm | 0.984 / 681 | 1.000 / 331 | 0.999 / 338 | 1.000 / 330 |
| gain at 60 / 120 / 240 / 480 nm rms | 0.66 / fold / fold / fold | 0.98 / fold / fold / fold | 0.94 / fold / fold / fold | 1.01 / 1.02 / 1.06 / 1.13 |
| photons per measurement for 1 pm, at the camera | 5.4×10¹³ | 4.7×10¹³ | 3.3×10¹³ | 1.9×10¹⁴ |
| light kept (divide the row above by it for incident photons) | 1 | 1 | 0.82 | 0.75 |
- **Shutter frame** (the reference alone, once per state, 5 frames): [R2].
- **A 1 λ/D pinhole:** [R2].
- **A 2% step-size error,** four-step least squares against the five-frame Schwider–Hariharan scan: [R2].
::: right
![The record run's own figure: the focal spot with the 2 λ/D pinhole (red), the dimple (green) and the waveguide mode's 1/e contour (brown); the reference amplitudes across the pupil image; the reference's motion under the 30 nm surface by pinhole size (shape, red; amplitude, black; the waveguide's coupling, brown); the fringe visibility of both readings on the flat DM.](figs/pdi193fbase_pdi.png){h=3.3}
- **Range without a fold:** the self-referenced readings lose their reference amplitude with the Strehl ratio and fold at 120 nm rms; the waveguide reference is fixed, so the reading holds gain 1.0–1.1 to 480 nm (the floor grows to 1.6 nm because the response matrix was measured at 30 nm).
- **Light:** with the paper's 60/40 split the reference arm returns 35% of the light and the test beam keeps 40%; the modulation is a smaller share of the detected flux than the stepped pinhole's, whose reference rides on the same beam — 4× more camera photons per picometer, before the arm's own loss.
~ Same seeds, same actuator-space scoring, same code as the sensor rows (dm_gauge_lib).  "fold" = gain −0.01: the reading returns nothing.  P at unit transmission and the dimple's size equals S to 5×10⁻¹⁵ (check G7).

## In closed loop: the hold metric with a camera that drifts | [CAM] filled from runs/ploop193, pcam193, pcam193i
::: left
- **The loop:** the DM held at the working surface at gain 0.5 for 60 cycles through one reading; the differential to the set point's frames fitted through the matrix; scored as the steady-state hold error.  The ZWFS rows on the same seeds: stepped 2.6×10¹² photons per cycle to hold 3 pm against noise, 7.5×10¹² against a 2 pm walk; polarized pair 1.5 / 5.3×10¹².
- **The camera knob:** every pixel's offset random-walks 0.13 electrons per cycle (about 1 electron over the run, the paper's 12-hour figure compressed), constant within a scan.  Readings whose step weights sum to zero (S, P, PF) subtract it exactly; the single-frame readings and the simultaneous pair see the walk and put it on the DM.
- **Within-scan drift** (the whole step across the frames): what the zero-sum readings pay: [CAM].
::: right
| reading | 3 pm held, noise only | 3 pm held, 2 pm walk | thermal floor | camera walk, 10¹⁵ photons: hold error |
| stepped dimple S | 2.6×10¹² | 7.5×10¹² | 10.0 pm | [CAM] |
| polarized pair V | 1.5×10¹² | 5.3×10¹² | 9.9 pm | [CAM] |
| linear L (1 frame) | 2.1×10¹² | 7.3×10¹² | 27.6 pm | [CAM] |
| stepped pinhole P | [CAM] | [CAM] | [CAM] | [CAM] |
| P/SRI waveguide PF | [CAM] | [CAM] | [CAM] | [CAM] |
~ dm_gauge_lib/dmg_loop, shared with the interferometer; drift seed 77 on the full actuator grid.

## Run it yourself | One parameter sheet and one script; the sensor runner carries the two point-diffraction readings as options
::: full
- **Two files:** pdi_params.m returns zwfs_params with the readings S, V, P, PF and the camera-drift loop kind set; pdi_run.m runs zwfs_run on it.  Every pinhole and waveguide setting is a documented knob (P.pdi): pinhole size, plate transmission, phase steps and scheme, shutter frame, reference iterations, split ratio, the waveguide's mode, a step-size error; the loop's camera walk and its within-scan fraction (P.loop).
| call | what it does |
| pdi_run | the record: bench checks, the full test set, figures |
| pdi_run('battery.calib_surface','base', 'battery.ladder_sites','grid') | calibrated on the working surface |
| pdi_run('pdi.DIA_LAMD',1.0, 'MODEL',2048, 'NGRID',385, 'param_file','macos_param_2048.txt') | a 1 λ/D pinhole, fully sampled |
| pdi_run('pdi.scheme','sh5', 'pdi.step_err',0.02) | the five-frame scan with a 2% step error |
| pdi_run('stages',{'bench','loop','figs'}, 'loop.drifts',{'cam'}) | the camera-drift hold rows |
| ./zwfs_batch.sh TAG "pdi_params, the same arguments" | unattended, memory-capped, logged; one engine run at a time |
~ Records: templates/40_benches/zwfs_dm96/runs/pdi193 (record 1), pdi193f / pdi193fbase / pdi193state / pdi193d1 / pdi193se_* (record 2), ploop193 / pcam193 / pcam193i (loop); README "P / PF"; plan and literature in macos/BRIEF_pdi_campaign.md.

## Conclusions | [R2/CAM] to be stated once record 2 and the loop rows land
::: full
- **The stepped pinhole is the stepped Zernike sensor with a smaller, attenuated reference:** the same frames and the same algebra (equal to 5×10⁻¹⁵ at unit transmission); it reads exactly where the dimple's four-step does not (0.93 vs 0.75 on the 30 nm surface), for 1.6× fewer camera photons, and pays 18% of the light at the plate.
- **A reference that does not depend on the surface is what extends the range:** the waveguide reading keeps reading past the surfaces at which every self-referenced reading folds.
- **[R2]** the on-surface calibration, the 1 λ/D pinhole, the scheme trade.
- **[CAM]** which readings the camera walk reaches, and what within-scan drift costs the zero-sum ones.
- **What the model does not yet carry:** the reference arm's own drift (non-common path), detector nonlinearity and persistence — the paper's own list.

# Backup

## Provenance and records | Every number re-derives from a committed script
::: full
- Code: templates/40_benches/dm_gauge_lib/dmg_pdi_gauge.m (both readings, the reference models, the step schemes), dmg_loop.m (the camera drift), zwfs_run.m (readings P / PF, checks G5–G7, the loop's 'cam' kind), zwfs_run_figs.m (the point-diffraction figure), pdi_params.m / pdi_run.m; tests/tDmgLoop.m (10 checks, including the camera immunity).
- Runs: zwfs_dm96/runs/pdi_dev, pdi_dev2, pdi_dev3 (dev-resolution checks); pdi193 (record 1); pdi193f, pdi193fbase, pdi193state, pdi193d1, pdi193se_ls, pdi193se_sh5 (record 2); ploop193, pcam193, pcam193i (loop).
- Paper and model: ~/dev/MACOS_sandbox/pSRI/Dube_pSRI.pdf; photonic-sri-simulation (coupleThroughPhotonicChipRecollimate.m, singleModeField.m, PSIAccumulate.m, differentialPSIReconstructAmpPhs.m).
- Plan and literature: macos/BRIEF_pdi_campaign.md.
