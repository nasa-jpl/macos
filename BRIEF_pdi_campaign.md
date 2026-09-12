# BRIEF: point-diffraction interferometer readings on the DM gauge bench

Dave 2026-09-12: "While CCL and CCMac work on tg96 and ZWFS, let's look at
another sensor, using a point-diffraction IFO approach -- see papers by
Brandon Dube for examples."  This brief is the plan and the running
record; numbers land here and in `zwfs_dm96/README.md` as they are
measured.  The pattern is the ZWFS campaign's (`BRIEF_zwfs_campaign.md`):
same DM truth, same battery, same loop, one more instrument.

## The literature, and the paper meant (`~/dev/MACOS_sandbox/pSRI/Dube_pSRI.pdf`)

- **Dube, Nejadriahi, Sidick, Jewell, Redding, Lou, Basinger, "Absolute
  and differential complex E field reconstruction by phase shifting
  interferometry", Proc. SPIE 13092, 130926F (2024)** -- the paper the ask
  points at, with its MATLAB model beside it
  (`pSRI/photonic-sri-simulation/`).  What it says, and what the model
  here takes from it:
  - HWO's coronagraph needs a few pm rms in phase and ~1% in amplitude
    over a ~12 h cycle; Roman's LOWFS error budget is dominated by camera
    drift (~200 pm-equivalent, ~1 electron per pixel over 12 h).  Temporal
    PSI high-pass filters such 1/f noise at the measurement repetition
    rate: with an integer number of modulation cycles the weights sum to
    zero, so any pattern constant within one scan (bias, fixed pattern)
    subtracts.  Continuous operation: a triangle-wave phase shifter,
    scans processed in alternating order.
  - Reconstruction: de Groot's DFT at the known modulation frequency,
    C = sum c_n I_n, S = sum s_n I_n (five-frame Schwider-Hariharan,
    alpha pi/2, s = {0 2 0 -2 0}, c = {-1 0 2 0 -1}); the complex field
    C + iS gives amplitude AND phase (eq. 5-11); the DIFFERENTIAL by
    complex division, cos/sin of (theta1 - theta0) from C1 C0 + S1 S0 and
    S1 C0 - C1 S0 (eq. 15-19) -- no unwrapping, integer arithmetic until
    the last step.  Both DMs of a series pair are disambiguated only
    because amplitude is sensed too.
  - The instrument, the P/SRI: a NON-common-path interferometer; a plate
    beamsplitter splits test and reference arms, the reference arm
    focuses onto a single-mode waveguide in a photonic chip (a clean
    reference "in a manner similar to a point-diffraction
    interferometer"), thermo-optic / strain tuning in the chip is the
    phase shifter, OAPs recollimate, a second beamsplitter recombines; run
    out-of-band with a monochromatic laser, a notch filter in the science
    channel.  Named after Medecki's PS/PDI.  Their code: the reference
    pupil field is the inverse transform of the exact step-index LP01
    mode (V 2.3, b 0.5, core radius 0.5 lam/D at the focus, the Thorlabs
    UV-fiber set), scaled by the overlap coupling of the state's focal
    field into the mode (a scalar; the shape is the mode's), beamsplitter
    R 0.6 to the photonic channel, T 0.4 to the test arm.  The reference
    intensity is to be measured by a shutter that blocks the test arm.
  - Calibration: a sparse, windowed interaction matrix built in situ by
    multiplexed actuator pokes (their Fig. 4; 2048^2 camera, 19 px per
    actuator, 5x5-actuator windows) -- Dave's measured response matrix
    (S10) is the same object, already the campaign's default.
  - Simulated example: EAC1 pupil + Roman HLC DM solutions uptiled to
    96x96, sensed at 399 nm; absolute field reconstructed (a systematic
    in the intensity from the Gaussian reference, calibratable by the
    shutter frame); a 5 mK DM thermal soak (2.6% per K, ~20 pm)
    reconstructed differentially to 1e-14 nm on top of a wrapped absolute
    OPD.  No noise, no non-common-path drift, no detector systematics
    yet; "a later paper".
  - Companions: Redding et al., "Wavefront sensing and control for a
    future Habitable Worlds Observatory", SPIE 13092, 130921P (2024); Lou
    et al. SPIE 13129-22 (2024); the FY24 SRTD poster.
- **Dubost, Bharmal, Dubbeldam, Myers, "Concept validation of a high
  dynamic range point-diffraction interferometer for wavefront sensing in
  adaptive optics", Appl. Opt. (2022), arXiv 2204.04940** -- the
  pupil-modulated PDI (m-PDI / CAWS): a Ronchi grating splits the beam,
  order 0 through a 2.5 lam/D pinhole (reference, b0 = 0.42 of the
  amplitude), order +1 through a larger window (test), interfered with a
  spatial carrier and demodulated from ONE frame by a Fourier sideband.
  Visibility 0.96, throughput 0.28 (0.56 for the mask alone), +-pi range
  (unwrappable, 5.1 pi PV measured), closed loop to 7.7 nm rms (central
  72% of the pupil, 47.7 nm over all), 77 nm FWHM polychromatic to the
  same residual.  Their eq. 5 gives the slope limit A < lam/(6 pi T kappa)
  set by the sideband window; the second-harmonic filtering explains a
  20% high-frequency transfer drop at 1 rad.  The single-frame carrier
  form is a later reading here.
- Classical PDI: Smartt & Steel (1975); phase-shifting PDI Medecki,
  Tejnil, Goldberg, Bokor, Opt. Lett. 21, 1526 (1996); EUV PS/PDI
  Naulleau et al., Appl. Opt. 38, 7252 (1999), reference-wave accuracy
  lam_EUV/330; PSI algorithms de Groot (1995), Surrel (1996), Servin
  (2009), Hariharan-Oreb-Eiju (1987); self-referencing interferometers
  for AO (Notaras & Paterson 2007, closed loop through optical vortices).
- The ZWFS is a PDI: Dubost say so in their introduction, and the model
  below makes it literal -- the stepped pinhole reading at unit surround
  transmission and the dimple's diameter IS the stepped Zernike reading S
  (gate G7).  The scan for the ZWFS side is `REPORT_zwfs_lit_scan.md`.

## What the PDI changes, physically

The Zernike sensor's reference is the core of the spot passed through a
phase dimple -- it MOVES with the state (14% under the 30 nm working
surface, the S7 "moving core"), which is why the exact readings iterate
b, why the one-frame readings fold at the quarter wave, and why the
campaign needed the stepped (S) and vector (V) readings to get an exact
solve.  A point-diffraction interferometer makes its reference from a
PINHOLE: only the core passes, the diffracted wave is a clean sphere
whose shape depends weakly on the aberration, and the surround is
attenuated so the two amplitudes match (visibility ~ 1).  With phase
steps the solve is the classical four-step: exact, linear, +-pi range,
amplitude AND phase, no fold, no branch prior.  What it costs is light:
the surround transmission t^2 (a few percent to tens of percent), or in
the fiber form the pickoff to the reference arm plus the pinhole's own
coupling loss.  The comparison the campaign can make, on the same DM
truth with the same actuator-space scoring:

| reading | frames | reference | solve | light |
|---|---|---|---|---|
| ZWFS S (stepped dimple) | 4 | moves with the state (|b|^2 calibrated on the flat) | rank-2 stepped retrieval | all of it |
| ZWFS V (polarized dimple) | 2 | moves; b iterated | exact pair | all of it |
| **P** (stepped pinhole, common path) | K (4) | pinhole-diffracted, nearly fixed; iterated | exact four-step | t^2 of it |
| **PF** (fiber reference, photonic steps) | K (4) | fixed by construction | exact four-step | 1 - pickoff, + the reference |

## The model (`dm_gauge_lib/dmg_pdi_gauge.m`, readings P / PF in `zwfs_run`)

Same bench as the ZWFS (the TG96 test arm, the FocalMask inside the
symmetric NF sandwich, the detector at the DM-conjugate pupil), so the
only change is the mask and the solve.

- **P, pinhole**: mask V_k = t + (e^{i theta_k} - t) D at the FocalMask,
  D the gray-edged pinhole disk (`zwfs_mask` at phase 0), t the surround
  amplitude transmission (`pdi.t_surr` 'auto' = the reference's rms
  amplitude / the beam's), theta_k the steps (`pdi.thetas`, default 0
  pi/2 pi 3pi/2).  Detector field t E + c_k b, c_k = e^{i theta_k} - t,
  b = T(D Ti(E)) the state's own pinhole reference (the FFT surrogate,
  gated against the engine at 1e-10).  Per pixel I_k = A + B cos theta_k
  + C sin theta_k with B = 2t(Re X - |b|^2), C = 2t Im X, X = E conj(b),
  A = t^2|E|^2 + (1+t^2)|b|^2 - 2t^2 Re X: a 3-parameter fit over the K
  steps; |b|^2 from the flat's pinhole-only frame (`pdi.b2` 'flat',
  iterated with the phase) or a pinhole-only frame per state ('state',
  K+1 frames, exact); the reference's phase iterated from the estimate
  `pdi.NITER` times as the ZWFS exact readings do.  Phase relative to the
  flat = angle(X) + angle(b) - angle(E0).
- **PF, the P/SRI (fiber reference)**: detector field s E + a kappa
  e^{i theta_k} R, s = sqrt(1 - f) (`pdi.pickoff` f = beam power to the
  reference arm, 0.6 as in the paper), R = the recollimated LP01 mode of
  the waveguide at unit rms over the pupil (`pdi.ref_shape` 'fiber':
  V 2.3, b 0.5, core radius 0.5 lam/D -- their set; 'pinhole' = the
  first idealization, the pinhole-diffracted flat field), kappa = the
  state's coupling into the mode relative to the flat's (a complex
  SCALAR: the shape is fixed by construction), a the reference amplitude
  (`pdi.a_ref` 'auto' = min of the visibility-1 match and the pickoff
  budget a^2 sum|R|^2 <= f |c0|^2).  One trace per state, the K frames
  synthesized; X = E conj(R) exact in one pass, no |b|^2 degeneracy;
  the solver takes kappa = 1 (`pdi.b2` 'flat', the paper's differential
  mode) or |kappa| from a shutter frame ('state').  Schemes: `pdi.scheme`
  'ls' (least squares over `pdi.thetas`) or 'sh5' (the paper's five-frame
  Schwider-Hariharan, de Groot weights); `pdi.step_err` miscalibrates the
  steps in the frames only.  IDEAL: no arm drift, no detector drift (what
  the photonic modulation buys on hardware); both belong to the loop
  stage as knobs.
- Differential = the WRAPPED phase difference (V1's lesson); height
  S_CONV phi lam/(4 pi); photons counted at the DETECTOR (the campaign's
  currency) with `throughput` = detected/incident printed so any number
  can be restated in incident photons.

## Stages and gates (each prints; the report is the record)

- **Bench (G5-G7, in `stage_bench_`)**: G5 exact beyond the one-frame
  fold on the same 100 nm sparse pokes as V's G4 (< 0.1% of the figure);
  G6 the reference's motion under the 30 nm working state
  |b_state - b_flat|/|b_flat| for pinhole diameters 0.5 .. 3 lam/D
  against the dimple's -- the PDI's argument, MEASURED; G7 the pinhole
  reading at t = 1 and the dimple's diameter reproduces the stepped
  Zernike reading S (< 1e-2): same instrument, same frames.  Sampling
  budget: the pinhole gets the dimple's rule (>= 6 px at the mask plane;
  3.96 px per lam/D at 1024/193, so 1.5 lam/D warns and 1 lam/D needs
  MODEL 2048).
- **Battery** (matrix calibration, the five rows, the ladder): P and PF
  beside S and V.  Expected: exact readings, so the flat rows track V;
  the ladder should hold to the lam/2 wrap like V, with no core collapse
  for a small pinhole -- if the 60 nm whole-pupil cliff is gone for P,
  that is the reference-stability argument made in actuator space.
- **Noise**: N(1 pm) per reading at the detector, then / throughput.  The
  PSI form spends its photons on K frames and on a reference that carries
  no signal; expect P and PF to sit right of L/V per DETECTED photon and
  further right per INCIDENT photon.  The number to report is both.
- **Loop**: the hold metric rows beside S and V (same seeds, `dmg_loop`).
  Then the PF-specific knob: a reference-arm drift (piston is invisible;
  tilt and focus are not) per cycle, priced in the hold error -- the
  non-common-path cost the common-path P does not pay.
- **Sweep** (steer-dependent): pinhole diameter (0.5 .. 2 lam/D at MODEL
  2048), surround t about 'auto', K = 3/4/5 steps, pickoff f; the
  single-frame carrier form (Dubost) if a tilted reference is wanted.

## Dev-resolution gates (RUN 2026-09-12, `runs/pdi_dev`: model 512, NGRID 65, 48x48, 2 lam/D)

- P (pinhole 2 lam/D): surround t_auto 0.7215 (the 2 lam/D disk passes
  eta_pin 0.69 of the light -- at the dimple's size this is barely a
  pinhole; the classical regime needs 1 lam/D or less), throughput 0.84,
  visibility 0.55 on the flat; the flat reads 2e-15 rad; G5 1.87 pm of an
  11.8 nm sparse-poke figure (1.6e-4, PASS; the residual is the reference
  iteration at NITER 3, b2 'flat').
- PF (fiber reference, pickoff 0.5): a budget-limited (0.0090 vs the
  visibility-1 match 0.0125), visibility 0.92, throughput 0.85; G5 0.000
  pm -- exact, as the algebra says.
- G7: 5.35e-15 -- the pinhole reading at t = 1 and the dimple's diameter
  IS the stepped Zernike reading S.
- G6 as first written (total relative change of the pinhole reference
  under the 30 nm state) read 0.142 at 0.5 lam/D .. 0.167 at 3 lam/D,
  dimple 0.154: nearly flat in diameter, because the metric is dominated
  by a Strehl-class AMPLITUDE drop of the reference (exp(-sigma^2/2) ~
  0.84 at 0.6 rad rms) that the |b|^2 calibration / iteration absorbs.
  Split into scale and SHAPE before the record run; the shape number is
  the PDI's argument.

## The P/SRI re-pinned to the paper (dev gate `runs/pdi_dev2`, 2026-09-12)

PF with the LP01 reference: coupling eta_c 0.587 on the flat; pickoff 0.6,
a budget-limited (0.0081 vs the visibility-1 match 0.0112), visibility
0.863, throughput 0.752; G5 0.000 pm (exact); under the 30 nm working
state the coupling |kappa| = 0.849, arg 0.005 rad -- the reference's only
motion, an amplitude scale + a piston; the differential phase does not
see the scale at all (it divides out of angle(X)), so PF should hold on
any base until the fringes vanish -- the record ladder tests that.

## Record 1 (`runs/pdi193`, model 1024, flat matrix; PF here = the pinhole-shaped first idealization)

- G6 reference SHAPE change under the 30 nm state: 0.10% (0.5 lam/D),
  0.26% (1.0), 0.58% (1.5), 1.06% (2.0 = the dimple), 2.56% (3.0); the
  complex scale 0.86 at every diameter (the Strehl).  The PDI argument in
  one table: a 1 lam/D pinhole's reference moves four times less in shape
  than the dimple's.
- Rows (flat matrix; g / e pm): on the 30 nm surface V, P, PF agree --
  single 0.939 / 0.933 / 0.940 (25 pm), grid 0.997 / 0.992 / 0.998 (12
  pm), dense 0.999 / 0.993 / 1.001 -- where S reads 0.752 / 0.829 / 0.818.
  On the flat all four 0.994-0.996.
- Ladder (hold-out, flat matrix): P tracks V to 60 nm (0.78 vs 0.85) and
  folds with it at 120 nm (both -0.02: the reference's amplitude collapses
  with the Strehl, and P's |b|^2 'flat' assumption breaks); the
  fixed-reference PF holds 0.77 at 120 nm and 0.55 at 240.
- Noise, N(1 pm) photons per MEASUREMENT at the detector: S 5.4e13, V
  4.7e13, **P 3.3e13**, PF(pinhole shape) 2.1e14; throughput P 0.82,
  PF 0.84 (divide for incident photons).  The stepped pinhole beats the
  stepped dimple by 1.6x per photon; the first fiber idealization paid
  4x for its reference apodization -- record 2 re-prices it with the
  LP01 reference.

## Record 2 (`runs/pdi193f` flat matrix, `runs/pdi193fbase` matrix on the 30 nm surface; PF = the P/SRI with the LP01 reference)

- PF, the P/SRI: coupling eta_c 0.587 on the flat; a budget-limited
  (0.73 of the visibility-1 match at pickoff 0.6); visibility 0.863,
  throughput 0.752; G5 0.000 pm; |kappa| 0.860 with 0.025 rad under the
  30 nm surface.
- Flat matrix: PF reads as V and P on every row (single 0.940 / 25 pm,
  grid 0.998, dense 1.001); the ladder holds 0.77 at 120 nm, 0.55 at 240
  where S / V / P fold; N(1 pm) at the camera 1.9e14 (S 5.4e13, V
  4.7e13, P 3.3e13) -- per detected photon the waveguide form is the
  most expensive of the four, because with the paper's 60/40 split 60%
  of the light goes to an arm that returns 59% of it as reference and the
  test beam keeps 40%: the modulation is a smaller fraction of the
  detected flux than the stepped pinhole's, whose reference rides on the
  same beam.  (In incident photons: divide by 0.75 / 0.82.)
- **Matrix on the 30 nm surface (the campaign's operating point):**
  single 10 nm -- S 0.9885 / 5 pm / SNR 2120, V 0.9935 / 4 / 2835,
  **P 0.9935 / 4 / 2790, PF 0.9935 / 4 / 2842**; grid 1 nm -- all four
  0.999 / 3-4 pm; dense random 10 nm -- S 0.984 / 681 pm, V 0.9999 /
  331, P 0.9985 / 338, PF 1.0002 / 330.  The three exact readings are
  indistinguishable at the operating point; S carries twice the dense
  error.
- **Ladder (47 grid sites, matrix from 30 nm):** at 60 nm S 0.66, V
  0.98, P 0.94, PF 1.006; at 120 nm S / V / P fold (-0.01), **PF 1.02**;
  240 nm PF 1.06; 480 nm PF 1.13 (floor 1.6 nm: the 30 nm matrix
  applied to a 16x larger surface).  A reference that does not depend
  on the surface has no fold: the P/SRI's range is the wrap of the
  DIFFERENCE, not of the surface.

## Record 2, continued (`runs/pdi193state`, `runs/pdi193d1`; matrix on the 30 nm surface)

- **Shutter frame per state (`pdi.b2 'state'`, 5 frames):** the pinhole
  reading P becomes PF's twin on every row AND the ladder -- single
  0.9935 / 4 pm, grid 0.9992 / 3, dense 1.0002 / 330, and **1.02 / 1.06 /
  1.13 at 120 / 240 / 480 nm**: P's fold at 120 nm in record 2 was the
  flat |b|^2 assumption, not the pinhole.  Cost: a fifth frame; noise
  priced with the on-surface matrix in this run: S 9.3e13, P(5) 1.3e14,
  PF 2.7e14 photons per measurement at the camera.
- **A 1 lam/D pinhole (`pdi.DIA_LAMD 1.0`; 3.96 px at the mask plane,
  the budget line warns):** t_auto 0.28, eta_pin 0.24, throughput 0.29,
  visibility 0.94; G5 0.002 pm (the reference iteration converges 100x
  better than at 2 lam/D); rows identical to PF's (0.9935 / 4; 0.9992 /
  3; 1.0002 / 330); **no fold with the flat |b|^2: 1.02 / 1.06 / 1.13 at
  120 / 240 / 480 nm** -- the classical PDI regime buys the P/SRI's
  range in the common path, at 29% of the light (N(1 pm) at the camera
  2.9e14 vs S 1.0e14 in this run; 1e15 incident).
- The 2% step-error runs (`pdi193se_ls`, `se_sh5`) aborted at gate G5
  (0.1% exactness), which a deliberate step error is meant to violate;
  the gate is now informational when `pdi.step_err` is set, and the two
  runs are re-queued.

## The loop rows (`runs/ploop193`, S11 seeds, matrix on the 30 nm surface)

| reading | 3 pm held, noise only | 3 pm held, 2 pm walk | thermal floor | noiseless step at cycle 60 |
|---|---|---|---|---|
| ZWFS linear L (1 frame) | 2.1e12 | 7.3e12 | 27.6 pm | 1.2 pm, still falling |
| ZWFS stepped S (4) | 2.6e12 | 7.5e12 | 10.0 pm | 0.000 |
| ZWFS polarized pair V (2) | 1.5e12 | 5.3e12 | 9.9 pm | 0.000 |
| **stepped pinhole P (4)** | **2.3e12** | **7.0e12** | 9.9 pm | 0.000 |
| **P/SRI waveguide PF (4)** | **7.0e12** | **2.5e13** | 9.9 pm | 0.000 |

Both PDI readings contract at 0.509 per cycle (gain 0.98, as V) with no
fixed error; the stepped pinhole costs the same light as the stepped
dimple in the loop; the P/SRI costs 3x (its single-shot noise 1.7x P's
at every level: the 60/40 split).  The held residual's spectrum under
the walk is V's (0.25 / 0.72 / 2.2 pm in the < 4 / 4-12 / > 12 cycles
per aperture bands).  Loop figure: runs/ploop193/ploop193_loop.png.

## Camera drift in the loop (`runs/pcam193`; the relative form pending)

The paper's number -- an offset random-walking 0.13 electrons per pixel
per cycle, ~1 e over the 60-cycle run -- is INVISIBLE to every reading at
1e13 and 1e15 photons per cycle: hold error, bias and spectrum equal the
noise-only rows to the printed digit (L 1.39 / 0.14 pm, S 1.52 / 0.15, V
1.17 / 0.12, P 1.45 / 0.14, PF 2.50 / 0.25 at 1e13 / 1e15, with and
without).  A lit pixel collects ~3e8 photons per frame at 1e13 per
measurement, so an electron is 1e-4 of its shot noise; the PSI immunity
argument lives in Roman's photon-starved LOWFS regime (1e2-1e3 per pixel
per frame, integrated over 12 h).  The knob therefore has a relative
unit now (`loop.cam_unit 'rel'`: a fraction of the mean photons per lit
pixel per frame, per cycle -- a bias / gain drift scaled to the signal);
`runs/pcam193r` (1e-3 per cycle, constant within a scan) and
`pcam193ri` (the whole step within each scan) are queued behind the
peer session's ZWFS V3 sequence; `pcam193i` (the electron form,
within-scan) is invisible as well.  **First relative-form run
(`runs/pcam193r_perframe`, 1e-3 per cycle, 14:32) -- the single-frame
reading L imprints the walk at 10.8 nm (floor 1.4 pm without it), the
pair V 89 pm, and the zero-sum readings S / P / PF 32 / 34 / 130 pm --
NOT the exact immunity the unit test guarantees.  Cause, mine: the
'rel' scale was each FRAME's own mean, so the frames of one scan got
different offsets; a camera bias is the same electrons on every frame.
Fixed to one scale per scan (the mean over the reading's frames);
`pcam193r` re-queued; the numbers above are superseded.**  The unit-test gate (tDmgLoop
G8) already shows the mechanism: a zero-sum reading is exactly immune, a
single-frame reading imprints the walk.

## Decision points for Dave

1. Which paper / layout is meant: if the 2024 paper's reference is a
   coronagraph-internal pickoff, PF is the model; if it is a pinhole at
   an intermediate focus, P is.  Both run; the record can carry either.
2. Pinhole diameter of record (2.0 lam/D = the dimple, clean comparison,
   budget PASS; or 1.0 at MODEL 2048, the classical PDI regime).
3. Whether the PDI gets its own campaign directory or stays as readings
   of `zwfs_dm96` (as built: readings P / PF, opt-in, the record runs
   untouched by default).

## Records

`zwfs_dm96/runs/pdi_dev` (dev-res gates), `runs/pdi193` (record),
README section "P / PF", this brief's numbers, memory
`project_tg96_gauge`.  Push only on Dave's review.
