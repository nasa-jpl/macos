# ATMOS against Kolmogorov theory -- the check of 2026-10-18 (CC, for Dave's ruling on 1.2 (p))

## The generator (`AtmosPhase`, sourcsub.F) is RIGHT

- Filter: sqrt(PSD) = 0.151 r0^(-5/6) f^(-11/6) -- Noll's 0.023 r0^(-5/3) f^(-11/3) (rad^2 per (cycle/m)^2),
  sqrt(0.023) = 0.1517.  Times the frequency measure 1/(N dx) so each mode carries PSD df^2.  The DC bin
  takes the first bin's value (piston is removed afterwards).
- Noise: Hermitian-symmetric unit-variance complex Gaussians, the four self-conjugate bins real with
  sqrt(2); a pair (k, -k) contributes 2 PSD df^2 as it must.  Inverse transform unnormalized, as the
  measure in the filter requires.  Phase (rad at the r0 wavelength) -> OPD by lambda/2pi: achromatic.
- Tilt: the least-squares piston and tilt over the disk are removed, then a random tilt of edge-phase
  variance pTilt^2 * 1.8 (D/r0)^(5/3) is re-added -- Noll's tilt term: Delta_1 - Delta_2 = 1.030 - 0.582 =
  0.448 per axis in a_2 units, x4 for the phase at the edge = 1.79.  Exactly the code's 1.8.  This is the
  standard remedy for the FFT screen's missing low-frequency power (no subharmonics).
- Measured by an exact Python replica with the engine's own RNG (`cli_tests/tools/atmos_screen_emulation.py`),
  seeds 12345 / 4321 / 777, 256-grid, dx 0.08 m, r0 0.25 m: tilt-removed rms over the 4 m disk 5.4 / 11.3 /
  10.3 rad, D_phi(1 px) 0.86 / 0.80 / 1.06 rad^2 against Kolmogorov 1.03, residual imaginary part 0.

## The RNG (`RNG` / `GAUSS` / `RNGIN`, mathsub.F) is RIGHT

Two LCGs (moduli 32717, 32713) XORed; replicated exactly: unit-variance Gaussians, no short period,
no degeneracy when `SEED n` loads both states with the same n (the moduli differ, the states diverge at
once).  The start-up seed is time-based (`macos_init.inc`: 29 * (StartTime mod 1000)), which is why an
unseeded run is not reproducible.

## What the ENGINE's rays receive is WRONG (gfortran release, `seed n` before `atm 1 0.25 1e-6 0 ray`, Cassegrain)

| seed | screen at element 1 (OPD 1 after ATM) | lag-1 D_phi / Kolmogorov | adjacent-difference sign | OPD at the focal plane |
|---|---|---|---|---|
| 12345 | 8.0 rad (1.28 um) | 160 | -1.7 (checkerboard) | 24.7 um |
| 4321 | 2 409 rad (0.38 mm) | 9e6 | -0.9 | 16.8 mm |
| 777 | 306 481 rad (49 mm) | 2e11 | -1.1 | 1.17 m |
| unseeded | 5.7 rad | 71 | -0.9 | -- |

Deterministic per seed on gfortran (release and debug give 2.469520091e-05 at the focal plane every
time for seed 12345); on ifx the SAME seeded sequence gives different results run to run (release: 4 of
4 lost every ray; debug: 1.03 m, lost, 2.47e-5, lost) -- an uninitialized read that only ifx's stack
exposes.  Lag 4-16 D_phi sits at 0.5-0.8 of Kolmogorov (the FFT screen's known low-frequency deficit);
the lag-1 excess with an alternating sign is a pixel-to-pixel artifact the generator does not have, so it
is made where the rays READ the screen: `AddAtmos` (ray index -> R1 index, the 6-term gradient stencil,
`rsamp = dxElt(iElt)`), or the OPD's own path after ATM.  The "lost every ray at element 2 (obscured)"
outcome is the same defect at a larger amplitude: the rays are thrown off the primary.

## Three more facts for the fix

1. `C1 = 0.5/dx` was right for the one-row central difference it was written for; the stencil now sums
   THREE rows, so the ray tilt is 3x the gradient.  The spot-vs-gradient ratio measured 0.7-0.9, not 3,
   but the measured gradient is dominated by the checkerboard, so the ratio is not yet a test.
2. `ztdir = 0` (jzlou): the direct field multiply `WFElt *= exp(i k R)` is disabled; the screen reaches a
   propagation only through `CumRayL`.  Decide which path is the model (both would double-count).
3. The screen's pitch is `dxElt(iElt)` of the square ray grid and the rays are mapped by grid index, so
   the model is valid only on a flat element at a pupil with the full square grid (the code's own
   comment: "this should be at a pupil!") -- the placement rule Dave asked about.
