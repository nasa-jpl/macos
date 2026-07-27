<!-- GENERATED FILE -- do not edit.
     Source: 90_coverage.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# Gate index

Every claim on this page, its measured value, the truth it is compared
against, and the test that pins it in continuous integration.

| Gate | Claim | Measured | Truth | Pinned by |
|---|---|---|---|---|
| 1.1 | geometry invariant under `ifPol` | 0.000e+00 waves, bitwise 1 | 0 / true | `tPolarization`, this driver |
| 1.2 | coating round-trip identity | 0.000e+00 | 0 | `tPolarization/test_coat_roundtrip_identity` |
| 1.3 | polarization-off is a no-op | 6/6, bit-identical (vs-ref = 0.000e+00) | bit-identical | GMI regression *(external)* |
| 2.1 | lossless Jones pupil is unitary | D 2.61e-15, δ 5.79e-16 rad | 0 | `tJonesPupil/test_unitarity_gate` |
| 2.2 | per-ray coefficients match Fresnel | 1.20e-14 (mag), 3.16e-14 rad (phase) | closed form | `tJonesPupil/test_fold_fresnel_analytic` |
| 2.3 | 2θ symmetry on a symmetric system | azimuth lock 2.18e-13 rad | 0 | `tJonesPupil/test_2theta_symmetry` |
| 2.4 | D basis-invariant; s/p retardance is artifact | 6.70e-16; 247.1× inflation | 0; ≫1 | `tJonesPupil/test_basis_invariance_and_sp_artifact` |
| 2.5 | two-mirror reduces to polarization astigmatism | other terms 7.18e-15, circular 8.86e-16 | 0 | `tJonesPupil/test_pol_zernike_two_mirror_form` |
| 2.5 | on-axis diattenuation extrapolates to zero | 3.21e-05 | 0 | `tJonesPupil/test_pol_zernike_two_mirror_form` |
| 2.6 | decomposition algebra exact | synthetic recovery to 1e-12 | exact | `tJonesPupil/test_pol_maps_synthetic_identity` |
| 3.1 | polarized-scalar ≡ scalar | bitwise 1 | true | `tVecChain/test_polarized_scalar_is_bit_identical` |
| 3.2 | vector ≡ scalar, every leg and state | 6.42e-16 worst case | 0 | `tVecChain/test_vector_equals_scalar_every_state` |
| 3.3 | energy conserved per leg | 0.00e+00, 1.23e-16 | 0 | `tVecChain/test_energy_conserved_per_leg` |
| 3.4 | masks act on the vector path | 1.14e-16 | 0 | `tVecChain/test_mask_throughput_identical_on_vector_path` |
| 3.5 | far-field normalization unified | 1.67e-15 | 0 | `tVecChain/test_far_field_vector_matches_scalar_normalization` |
| 3.6 | scalar physics undisturbed | 4.836e-13 reproduced exactly | committed residual | pymacos PROPER suite *(external)* |
| 4.1 | one mirror matches the perfect-conductor closed form | 6.39e-14 (transverse), 5.03e-14 (longitudinal), retardance 1.07e-16 | closed form | `tPolarization/test_odd_mirror_crosspol_pec_analytic` |
| 4.2 | odd-mirror cross-pol is slope-driven and bounded | slope 1.871, 1.034 × the bound | ρ²; O(sin² AOI) | `tPolarization/test_odd_mirror_crosspol_rho2_law` |
| 5.1 | the decomposition invents no floor | cross 0.000e+00; contrast curve 2.525e-15 | 0 | `tPolContrast/test_reduction_to_scalar_contrast_curve` |
| 5.2 | co/cross split is a unitary basis change | 1.343e-24; closure 1.802e-16 | 0 | `tPolContrast/test_parseval_on_the_split`, `test_energy_bookkeeping` |
| 5.3 | uncoated two-mirror floor | 7.0612e-07 of co-polarized power | §4 ray-level value | `tPolContrast/test_floor_reported_by_component` |
| 5.4 | analyzer is derived, circular state included | 1.110e-16; conjugate analyzer 7.130e+05 | 0; ≫1 | `tPolContrast/test_analyzer_tracks_input_state`, `test_analyzer_would_be_wrong_if_conjugated` |
| 5.5 | floor moves with coating choice | 27.9× (Al), 151.3× (MgF₂/Al) | monotone, ≫1 | `tPolContrast/test_coating_sensitivity` |
| 5.6 | Tranche-1 shortfall is detected, not hidden | carried 0.8412 bare, 0.5653 coated | 1 (would be) | `tPolContrastCoro/test_tranche1_shortfall_is_detected` |
| 5.7 | coronagraph floor, lower bound | 5.787e-13 mean contrast | — | `tPolContrastCoro/test_floor_reported_by_component` |
| 6.1 | polarization-off ≡ Reference surfaces | bitwise 1 | true | `tPolElement/test_unpolarized_bit_identical_to_reference_twin` |
| 6.2 | Malus's law; exact crossed extinction | 1.079e-13 of I₀; crossed exactly zero 1 | cos²θ; 0 | `tPolElement/test_malus_law`, `test_crossed_polarizer_extinction` |
| 6.3 | retarder Stokes, with the SIGN of S₃ | 0.000e+00 (signed S₃/S₀ = −1), linear residual 1.741e-16 | closed form | `tPolElement/test_qwp_linear_to_circular` |
| 6.4 | half-wave 2θ law; two QWPs ≡ one HWP | slope residual 8.882e-16; cascade 1.464e-16 | 2; 0 | `tPolElement/test_hwp_rotates_by_2theta`, `test_two_qwp_equal_one_hwp` |
| 6.5 | retarder is unitary, linear and circular in | 0.000e+00; JᴴJ−I 0.000e+00 | 0 | `tPolElement/test_waveplate_is_unitary` |
| 6.6 | the diffraction grid carries the train | 1.833e-15 | cos²θ | `tPolElement/test_grid_carries_the_polarizing_train` |

# Coverage and gaps

Stating what is *not* here is part of the report. A validation document that
only lists passing gates invites the reader to assume the uncovered ground
was checked.

## Open items in the work this document covers

* ~~The out-of-plane attribution is unverified.~~ **CLOSED** (§3.5). The
  plane-selectable complex-field getter now exists, the per-plane
  contributions are measurable, and the difference decomposes into two
  mechanisms rather than the one originally assumed. What remains
  unexplained is the 2.8983e-04 residual of that decomposition — a shape
  difference between the scalar field and Ex, consistent with their
  different seeds, not further verified and not relied upon.
* **Tranche 1's validity condition is a real restriction, not a formality.**
  It holds for chains whose elements between physical legs are non-polarizing.
  A chain with a coated or reflective surface between legs — an off-axis
  paraboloid between coronagraph masks, or an interferometer recombiner
  followed by folds under near-field propagation — needs Tranche 2 and is
  **not** validated by anything on this page.  §5.6 now quantifies the cost
  on such a chain rather than leaving it qualitative: the diffraction grid
  carries 0.8412 of the ray-level cross-polarized fraction on
  `Rx_Coro`, 0.5653 once its mirrors are coated, and the
  coating sensitivity it reports there has the wrong **sign**.
* **The coronagraph floor in §5.7 is a lower bound on that chain**, and the
  chain is a coaxial model with uncoated normal-incidence mirrors — so the
  number is the numerical-aperture term alone, not the floor of a built
  instrument.  The coating sensitivity that *is* trustworthy is §5.5's, on
  the single-leg two-mirror fixture.
* **The phase relation between the per-ray field and the diffraction grid is
  train-dependent.** The assembly applies no convention bridge at the seed,
  which is correct, but the reason is structural rather than a fixed
  convention: return legs are subtracted by the cumulative-path bookkeeping
  and added by the per-surface field phases, so the relation differs between a
  return-terminated train and a plain trace-to-detector flow. Tranche 2 must
  carry phases explicitly against the path bookkeeping rather than assume a
  sign. The behavioural gates above, not a convention argument, are what carry
  the correctness claim.
* **The coated `Refractor` branch has no analytic gate**, before or after the
  §4 normalization, and it carries a normalization discrepancy against its
  own uncoated branch: the uncoated path multiplies transmission by the
  radiometric factor √(n₂cos θ₂ / n₁cos θ₁) and the coated path omits it.
  Measured with an index-matched single layer (optically a bare interface):
  coated/uncoated |Eₓ| = 0.816442 at normal incidence, exactly 1/√1.5 for
  that substrate. A coated lens therefore under-transmits by ~18% in
  amplitude relative to the same surface uncoated. Recorded rather than
  fixed — it changes results and wants its own decision and gate.
* **Vector mode repurposes the three wavefront planes as Ex/Ey/Ez**, so it
  handles one wavefront only — no multi-wavefront composition concurrently.
  This is a documented constraint, not a defect to work around.
* **The reflective polarizer `RfPolarizer` is not implemented.** The
  off-normal axis convention no longer blocks it (§6.7 settles it and gates
  it at 20° of incidence); what remains is the physics a reflective wire grid
  carries beyond the axis rule — grid reflection efficiency and the
  substrate's own s/p response.
* **The polarizer and waveplate are thin idealizations** — no ray splitting,
  walk-off, face reflections or substrate — so a polarizing beamsplitter is
  modelled as two traces, or better as a coated `Reflector` at its working
  angle. Their outputs are purely transverse; a longitudinal component at
  the surface is discarded rather than tracked. The waveplate's retardance
  is also independent of incidence angle, where a real crystal plate's is
  not — the field-of-view effect behind compound and Pancharatnam designs,
  which would need a birefringent-plate model with o/e indices and a
  thickness.

## Phases not yet represented

Each appends its own evidence section here as it lands; the list is the
outstanding work, not a disclaimer.

| Phase | Deliverable | Evidence it will add |
|---|---|---|
| 2d | interferometer polarization metrology | visibility budget predicted vs simulated; phase-shifting systematic vs input state; the angle-of-incidence trade |
| 3 | vector vortex; reflective polarizer | vortex null depth vs retardance error against the analytic curve. The polarizer and waveplate landed — section 6. `RfPolarizer` is held for grid-efficiency and substrate s/p physics, not for the axis convention, which §6.7 settles and gates; the vortex is a separate design question (retardance and handedness conventions, broadband leakage) |
| 3a Tranche 2 | Jones-through-chain | chains with coated surfaces between legs |
| 4 | spatially variable coatings | per-segment and radial-grade Jones pupils; grid-resolution speckle-floor sweep |

## Ladder rungs not yet climbed

The validation ladder for polarization ground truth runs, in increasing
effort and increasing independence:

1. analytic single-surface Fresnel — **done** (§2.2)
2. unitarity of a lossless Jones pupil — **done** (§2.1)
3. rotational-symmetry 2θ invariant — **done** (§2.3)
4. interferometer self-measurement cross-check — **not done**; needs Phase 2d.
   Puts a known coating in the test arm of a simulated Twyman-Green and
   *measures* its retardance with the phase-shifting pipeline, then compares
   against both the Fresnel analytic and the Jones-pupil decomposition —
   three independent paths to one number, entirely in-house.
5. published two-mirror polarization-aberration results — **partly done**
   (§2.5). The low-order expansion exists and the measured system matches
   the *analytic form* the literature predicts: 2θ azimuthal structure,
   quadratic radial law, no circular component, everything else at
   round-off. What is still missing is a numeric regression against a
   specific published system, which needs that system set up here rather
   than more machinery.
6. cross-check against a commercial polarization-ray-trace code — **not
   done**. This is the only rung that satisfies an outside reviewer on its
   own, and it depends on license access rather than on effort here.

Rungs 1–3 are exact analytic or invariant checks and are the strongest
evidence available without a second implementation. Rung 3 in the *scalar*
sense is already met by the PROPER comparison (§3.6), which is why that gate
is worth as much as it is: it is the one number on this page produced by code
that shares nothing with MACOS.
