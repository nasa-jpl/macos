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
| 7.1 | uncoated transmission IS the power transmittance | 2.220e-16 normal; 2.220e-16 p, 0.000e+00 s at 45° | characteristic-matrix T | `tPolRadiometric/test_uncoated_transmission_is_the_power_transmittance`, `test_uncoated_transmission_oblique_s_and_p` |
| 7.2 | index-matched layer ≡ bare interface | 0.000e+00 normal; 2.220e-16 p, 2.220e-16 s at 45° | 0 | `tPolRadiometric/test_index_matched_layer_equals_bare_interface_normal`, `..._oblique` |
| 7.2 | and on the diffraction grid | 0.000e+00 | 0 | `tPolRadiometric/test_index_matched_layer_at_the_detector_plane` |
| 7.3 | MgF₂ quarter-wave vs the textbook multilayer | 2.220e-16 normal; 1.332e-15 p, 1.110e-16 s at 45° | characteristic-matrix T | `tPolRadiometric/test_mgf2_quarterwave_normal_incidence`, `..._45deg_s_and_p` |
| 7.4 | air-to-air power closure, mixed plate | 2.220e-16 | T₁·T₂ | `tPolRadiometric/test_air_to_air_power_closure_mixed_plate` |
| 7.4 | radiometric factors telescope across a plate | 2.220e-16 | 0 | `tPolRadiometric/test_air_to_air_factors_telescope` |
| 7.5 | polarization state untouched (common real scalar) | 4.441e-16 | factor-free t_p/t_s | `tPolRadiometric/test_scalar_factor_leaves_the_polarization_state_alone` |
| 7.5 | the factor lives inside `ifPol` | 1 | bit-identical | `tPolRadiometric/test_pol_off_is_untouched_by_the_coating` |
| 7.6 | quarter-wave structure survives the scalar | 5.551e-16; ratio flat to 9.992e-16 | characteristic-matrix T(λ) | `tPolRadiometric/test_quarterwave_structure_survives_the_scalar_factor` |
| 8.1 | dielectric-on-metal matches a PUBLISHED model | D 2.828e-14, retardance 4.937e-14 | G. van Harten, F. Snik and C. U. Keller, "Polarization properties of real aluminum mirrors I. Influence of the aluminum oxide layer," PASP 121(878), 377-383 (2009), ±0.01 *(external)* | `tPolExternal/test_vanharten_machinery` |
| 8.1 | p̂ bridge measured, not assumed | PEC r_s/r_p = +1 exactly | +1 | `tPolExternal/test_phat_convention_bridge_is_zero` |
| 8.1 | the ratio estimator is frame-free | M11/M22 vs −M12/M21 | 0 | `tPolExternal/test_ratio_estimator_self_consistency` |
| 8.2 | the 4 nm oxide is distinguishable | exceeds 0.01 | > publication accuracy | `tPolExternal/test_omitting_the_oxide_exceeds_published_accuracy` |
| 8.2 | 4.12 nm vs the rejected 50 nm | exceeds 0.01 | > publication accuracy | `tPolExternal/test_rejected_50nm_oxide_is_excluded` |
| 8.2 | admittance assignment pinned by [1,2] sign | D > 0 and < 0.15 at 70° | published Fig. 1a axis | `tPolExternal/test_admittance_assignment_sign` |
| 8.3 | overcoat trade reverses across the quarter wave *(analytic)* | 0.18 vs 6.35 vs 0.0157 | sign reversal | `tPolExternal/test_overcoat_trade_reverses_with_optical_thickness` |
| 8.3.1 | the same reversal, measured by the **engine** | 0.2035× at 632.8 nm vs 5.2707× at 1 µm | sign reversal, 25.90× | `tPolContrast/test_overcoat_trade_reverses_across_the_quarter_wave_condition` |
| 8.3.1 | the reversal is the film's optical thickness, nothing else | achromatic control 5.2707×, lands on the 1 µm answer to 2.1e-08 | no reversal | `tPolContrast/test_overcoat_trade_reverses_across_the_quarter_wave_condition` |
| 8.3.1 | the quarter-wave condition is wavelength-invariant | 1.9e-08; metal-only leg 1.3e-08 | 0 | `tPolContrast/test_overcoat_trade_reverses_across_the_quarter_wave_condition` |

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
* ~~The coated `Refractor` branch has no analytic gate, and carries a
  normalization discrepancy against its own uncoated branch.~~ **CLOSED**
  (§7). The convention was decided
  (`REVIEW_POL_RADIOMETRIC_2026-07-28.md`) in favour of the incumbent
  power-amplitude rule, the coated branch was brought to it with a single
  factor applied once after the recursion, and the branch now has thirteen
  gates against the Abeles characteristic matrix — seven of which fail
  against the pre-fix engine. What remains open is narrower and stated
  there: transmission into an **absorbing substrate** has no meaningful
  power convention and is deliberately ungated, and the `Reflector`
  transmittance blocks are still `if(.false.)` dead code, so a coated mirror
  reports no transmitted leakage.
* **The cos θ term in the radiometric factor overlaps conceptually with
  beam-area bookkeeping that ray *density* also carries** when the grid is
  seeded for diffraction. Whether that is a double-count in oblique
  refractive seeding is an open audit question that **predates** §7 — the
  uncoated branch has always carried the factor, and the PROPER comparisons
  are mirror trains that never exercised it hard. Recorded here, deliberately
  uncoupled from the §7 landing, which changed nothing about it.
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
