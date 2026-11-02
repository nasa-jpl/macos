# NOTE — Kent Wallace's ZWFS reconstruction against the gauge's readings

CC, 2026-10-13.  Source: J. K. Wallace, "Zernike Reconstructor: Iterative ArcSin" (25 Oct
2022; `MACOS_sandbox/DM_gauge/Zernike WFS Reconstruction.pdf`) + his 2026-10-12 mail (the
vector form).  Compared with `dm_gauge_lib/dmg_zwfs_gauge.m` (readings I, I+, V) and the
record `zwfs_dm96/zwfs_s7iter_report.txt`.

## Kent's scalar recipe (the PDF)

Measure I_p (off dimple), I_z (on dimple), darks D; I'_p = I_p − D, I'_z = I_z − D.  Iterate:
A_p = √I'_p; E_p = A_p e^{iφ_i}; e = FT[E_p]; e' = e × circular amplitude focal-plane mask
(the dimple's footprint); E_b = FT⁻¹[e']; I_b = |E_b|², β = arg E_b;
φ_{i+1} = π/4 + arcsin[(I'_z − I'_p − 2 I_b) / (2 √(I'_p · 2 I_b))] + β; stop when
Δφ = φ_{i+1} − φ_i < threshold.  Surface = λ φ / 2π (his λ_nm φ_rad / 2π; ours carries the
double pass in S_CONV).

## Kent's vector recipe (the mail)

Two pupil images through dimple phases +π/2 and −π/2:
Φ = atan2[(I_{z,−π/2} − I_{z,+π/2}), (I_{z,−π/2} + I_{z,+π/2} − 4 I_b − 2 I_p)], then ordinary
phase unwrapping for range; I_b and β iterated per image as above.

## What the gauge already runs

**Scalar: Kent's recipe IS reading I** (`solveI_`, 2026-09-09).  Same model
I = |E + c b|², same per-pixel solve — ours is written as cos(φ − Θ) = (I − A² − |c|²|b|²) /
(2A|c||b|), Θ = arg c + arg b − θ₀, principal branch (φ − Θ) ∈ [−π, 0]; for the quarter-wave
dimple c = e^{iπ/2} − 1 = √2 e^{i3π/4} it is exactly his π/4 + arcsin[...] + β with the same
branch (φ − β ∈ [−π/4, 3π/4]).  Same b: re-propagated from the current estimate through the
dimple's amplitude mask (`C.bsur`, the FFT surrogate of the mask model, gated against the
engine's own Eb).  Same iteration (NITER, default 5; his: until Δφ < threshold).

**Vector: Kent's atan2 is reading V in the ideal-metasurface case.**  With I± = P² + 2b² +
2Pb[∓ … ] the pair gives −4Pb sin(φ − β) and −4Pb cos(φ − β), so his atan2 = (φ − β) + π:
exact per pixel, no branch choice.  Ours solves the same two circles' intersection per pixel
(`solveV_`), and in addition carries the real metasurface's retardance error and leak
(κ±, η) and the arm's polarization maps q_L, q_R (V2 / V3), reducing to his formula at
κ = 1, η = 1, q = 1.

## The three practical differences, and what the record says about each

1. **Amplitude.**  Kent: A_p = √I'_p from the state's own off-dimple image.  Ours: A = |E₀|,
   the flat's amplitude, by default; `solveI_(Ia, C, I0)` takes the state's own unmasked frame
   and uses √I0 — his form.  On a DM-conjugate pupil they are the same thing: record G3,
   |E|/|E_flat| on the 30 nm state, std 4.1e−16.  Where the pupil image is NOT conjugate
   (defocused DM, a real bench), his per-state amplitude is the right one and ours is the
   option to switch on.
2. **Stopping rule.**  Kent: Δφ < threshold.  Ours: fixed NITER, with Δφ recorded per
   iteration (`info.dphi`).  Record: near flat (20 nm hold-out) Δφ = 6e−3, 4e−5, 2e−5, 7e−6,
   3e−6, 1e−6 rad — converged in 3–4 iterations, map error 7.5e−7 rad; the fixed 5 is
   equivalent to a threshold ~1e−6.  A threshold stop is a two-line change; worth adding so
   the bench and the model stop the same way.
3. **Darks.**  Bench practice; the model has no dark frame (read noise is in the camera
   model, no offset).  Nothing to compare.

## The thing the record adds that the recipe does not say

The scalar iteration CONVERGES only where the phase is within the arcsin's fold.  On the
30 nm working surface (φ rms 0.54 rad, 7.8 % of lit pixels beyond the fold) the plain
iterated reading stalls: Δφ 0.41, 0.047, 0.025, 0.025, 0.025 … rad, map error 0.25 rad —
the iteration cannot choose the branch (record G5).  With a branch prior per pixel from the
stepped frames of the same state (reading I+), Δφ 0.49, 0.073, 0.015, 0.004, 0.0013,
0.00045 and the error is 2.8e−4 rad piston-removed.  The residual is a pure PISTON, the
sensor's null (the intensity is invariant under a common phase on E and b), which an
iterated b cannot recover and the engine's oracle b can (G4: 5e−15 rad).  Kent's vector form
has no fold — his "common phase unwrapping" is the range extension on the vector reading,
and it is what reading V does; its limit is the metasurface (V2/V3), not the algorithm.

## What to do

- Cite Kent on the gauge deck's I / V rows (his recipe, 2022; the stepped-prior branch
  choice and the metasurface terms are the gauge's additions) and close deck_cophi's
  "Kent's ZWFS reconstruction when it arrives".
- Add the threshold stop (`opt.DPHI_TOL`) and expose the per-state amplitude as a named
  option (`opt.AMP = 'state'`) in `dmg_zwfs_gauge` — no record moves (both are the default's
  equal on the record bench, by G3 and G6); gate = the s7iter numbers reproduced.
- For Kent: the fold number (7.8 % of a 30 nm surface beyond it) and the piston null are
  the two things his diagram does not carry; the stepped prior and the oracle b are our
  answers.  Draft reply: `DRAFT_email_kent_zwfs.md`.
