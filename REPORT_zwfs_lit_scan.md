# ZWFS literature scan -- modeling improvements for the DM-gauge model (2026-09-09)

Dave's ask: scan the ZWFS literature (Wallace, Jewell, Serabyn, ...) for
modeling-improvement ideas, and: would a physical (PZT) phase shift help
the interferometer?  Read in full: Ruane/Wallace/Steeves JATIS 2020
(arXiv 2010.10541), N'Diaye ZELDA II (1606.01895), Darcis multi-lambda
2025 (2508.15458), Doelman vector-Zernike 2019 (1811.08805); abstracts /
HTML: Wallace 2011 SPIE 8126, Chambouleyron phase-shifted 2024
(2409.04547), Haffert nonlinearity 2024 (2401.08090), HiCAT mid-order
2024 (2409.03411), Keck vZWFS 2024 (2404.08728), Wallace Keck piston
2022 (2205.02241), Moore & Redding NLZWFS for LUVOIR 2018 (SPIE 10698),
PIAA-ZWFS 2026 (2606.28136), Steeves Optica 2020 (OPN summary).

## The one import that targets our weak results

**Model-based reconstruction with an ITERATED reference wave.**  Ruane
eq. 36-37 / N'Diaye eq. 3-4 give the per-pixel exact (second-order)
solution from the masked image, A = sqrt(I0) (dimple-offset frame) and
b; Doelman (eq. 9-10), Chambouleyron 2024 and Haffert 2024 iterate:
estimate phi, re-propagate the estimated field through the mask model
to update b, repeat 3-5x; Darcis 2025 and Haffert generalize to gradient
descent through the forward model.  Our frozen-reference linear reading
IS the un-iterated case, and the "reference core moves with the working
state" effect we measured in S2b (defocus-on-base 0.743; linear gain
0.66 on the 30 nm state; the grid-on-base non-detection) is the frozen-b
error.  Our advantage: the engine gives exact E0 and Eb, so an FFT
surrogate for b (Ruane eq. 27-28) can be validated to round-off before
it is trusted.  Expected: grid-on-base recovers WITHOUT the depth ladder;
the 1.16-1.18 over-correction disappears.  In-house precedent: Moore &
Redding 2018 (polychromatic nonlinear ZWFS for LUVOIR).

## The rest, ranked

2. JPL sensitivity + systematic budget in their form: beta_p (Ruane eq. 4)
   from our measured E0/Eb predicts the S5 photon numbers analytically;
   the four systematic scale factors (pupil calibration 1/2, b error 1,
   dimple depth gamma/chi', initial phase chi/chi'); the asymmetric
   dynamic range -(theta)..(pi-theta) and the sensitivity collapse toward
   the negative bound (their Fig 2b) = our base-dependent gain, modelable.
3. Differential protocol as practiced: alternate DM states every 1000
   frames, median of differences, 4.4e5 frames -> 1 pm (Ruane);
   flat-vs-waffle alternation (Steeves, 1.6 pm in 4.3 s).  Add a drift
   term + the alternating-median estimator to our S5 noise model.
4. Calibrate the DM through the sensor's OWN image model (Ruane 4.4/4.6:
   match ZWFS images of poke grids to a simulation incl. propagation
   distances; out-of-conjugate actuators appear as donuts).  Our IFO's
   detector is NOT at the DM conjugate after the null-tuned tail -- this
   is the model-matched version of our registration doctrine.
5. Spot size revisited with the better reconstructor: Ruane 1.06 lam/D
   (uniform illumination, b accuracy), HiCAT 2.5, Keck 2.4, ours 2.0;
   "larger spot = more sensitivity but worse b accuracy" -- an iterated b
   removes the penalty; rerun our spot sweep on it.
6. Broadband/multi-color: Ruane saw no degradation to 25% bandwidth and
   notes broad band averages sub-actuator diffraction; flags the missing
   broadband reconstructor -- Darcis supplies it (joint per-lambda fit).
   Our S6 adds what they only imply: the fine-scale nulls migrate as
   1/lambda.  Merge = joint multi-lambda model fit instead of the
   equal-weight Wiener (also removes the stepped-reading regression).
7. Vector ZWFS before committing to the metasurface: Doelman's exact
   two-image phase+amplitude solution (eq. 7-8) and leakage terms
   (HWP/QWP retardance offsets, PBS rotation = calibratable cos^2 gain
   loss); Keck/SEAL limits were differential defocus between the two
   pupils and polarization crosstalk; fabricated shifts 0.30pi/0.68pi vs
   +-0.5pi.  Geometric phase is achromatic -> only the lam/D size scales
   with color -> multi-color gains shrink.  Model it with our pol machinery.
8. Segmented mirrors: HiCAT/Keck per-segment PTT interaction matrices;
   14-bit DM quantization steps; Keck underestimated piston 2-4x below
   50% Strehl.  The sensor sees segment piston, not global piston.
9. Physical phase stepping belongs on the ZERNIKE side: Wallace 2011
   all-reflective sensor with a dynamic, arbitrary core phase shift read
   by four consecutive measurements (phase AND amplitude);
   Chambouleyron's +-pi/2 pairs -> cos/sin of phi -> arcsine, iterative,
   ~1 rad rms range.  Our S2b depth ladder is the static three-mask
   version.
10. FALCO carries Ruane's reconstructors as open source -- cross-check.

## What we do that the literature does not (say it to Fang)
Exact engine-measured E0 and Eb (they compute b analytically or from a
flat pupil); the ray-affine registration doctrine; actuator-space
scoring through the DM's own influence model; the multi-depth structural
result (|c|^2 = -2 Re c: two observables per pixel); the spatial-
frequency-resolved null-migration measurement across colors (S6).

## PZT phase shift for the interferometer -- verdict
In the MODEL it changes nothing: the four-step is exact for ideal
elements and the floors (46 pm single actuator, 4.2 nm dense random)
are geometric imaging transfer (S6 proved them color-independent).  On
HARDWARE it changes the error budget: removes the polarization
systematics (invisible scale from diattenuation / waveplate azimuth:
11.7% at 45 deg, 0.15% at 7 deg) and gives an ABSOLUTE phase scale from a
known displacement -- the one systematic the polarization gauge cannot
see; costs 4 sequential frames (drift/vibration) vs one polarization-
camera snapshot; 5/7-step algorithms (Hariharan/Schwider) cancel linear
step error; no polarization optics (~2x light).  Recommendation: HYBRID
-- polarization snapshot for the differential measurements, a PZT on the
reference flat as scale calibrator + cross-check.  Model experiment
(~1 h): reference flat translated lam/8 per frame on the non-polarizing
two-arm bench, same battery, then +-2% step errors, 4-step vs 5-step.

## Suggested first ZWFS task
Implement the iterated-b reconstructor (item 1) as `dmg_zwfs_gauge`
reading #3 ("iterated"): per-pixel exact solve (Ruane eq. 36-37) with
A from the dimple-offset frame, b from an FFT surrogate of the mask
model, 3-5 iterations; gate the surrogate against the engine's Eb at
1e-10; then rerun the S3/S4 base rows and the grid-on-base scenario.
Success = grid-on-base SNR >= 5 from ONE frame + the hold-out gain
within 3% of 1 without the Wiener over-correction.
