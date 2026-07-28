<!-- GENERATED FILE -- do not edit.
     Source: 60_contrast_floor.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# The polarization contrast floor

The coronagraph question is not what shape the polarization aberration has
but how much light it puts where the dark zone should be. `pol_contrast_floor`
answers it by splitting the detector field into three channels —
**co-polarized**, **cross-polarized** and **longitudinal** — and reporting
the cross-polarized channel, the part no scalar deformable-mirror control can
remove, in peak-normalized contrast.

Two design choices carry most of the risk, and both are gated below.

**The split is taken at the detector, on the engine's own Ex/Ey/Ez component
planes, not through a pupil multiplier.** The chain is linear in the input
Jones state and Tranche 1 propagates all three planes with the identical
scalar kernel, so a spatially uniform analyzer commutes with propagation:
projecting after the propagation is the same operation as projecting before
it. This also avoids the obstacle that stopped two earlier designs — the
Jones pupil is assembled from `RayE` and carries the accumulated optical path
length, so it can never be used as a `WFElt` multiplier, and there is no
fixed conjugation to divide out because the `RayE`↔`WFElt` phase relation is
train-dependent (§3).

**"Co-polarized" is referenced to the mean *output* state, not to the input
state.** A train can rotate polarization geometrically with zero
diattenuation and zero retardance; charging that rotation to the
cross-polarized channel reports an aberration where there is none, since an
observer would simply align the analyzer to the output. The analyzer used is
the dominant eigenvector of the 2×2 pupil coherency matrix
`C_ij = Σ E_i conj(E_j)`, which is phase-insensitive — the common wavefront
cancels inside the product, so unlike a plain pupil mean it does not collapse
on an aberrated pupil — and which by construction minimizes cross-polarized
power. An unpolarized source is two traces with orthogonal input states
summed in intensity, each with its own analyzer; the second state is never
synthesized from the first.

## 5.1 The decomposition invents no floor

The first thing to establish is the negative one. On a train where
polarization does nothing, the machinery must report *nothing* — otherwise
every floor it produces afterwards is partly its own reference-frame
choices.

`Rx_VecChain.in` is that train: collimated, on-axis, flat uncoated planes, so
the ray field direction is a constant unit vector. With an x-polarized input
the cross-polarized channel comes out at 0.000e+00 of the co-polarized
power — **exactly zero**, not merely small — the longitudinal channel with
it, and the co-polarized channel reproduces the scalar (polarization-off) run
to 1.331e-15 of the detector peak. Taken as a contrast curve, the
quantity a coronagraph budget actually reads, the two agree to
2.525e-15 relative across 180 radial bins.

## 5.2 The split is a unitary change of basis

The analyzer and its complement form an orthonormal pair, so the two
transverse channels must account for the transverse intensity pointwise, with
nothing lost between them and nothing double-counted. Measured on
`Rx_Cass_FarField.in`: co + cross equals |Eₓ|² + |E_y|² to
1.343e-24 of the detector peak, and the three channels together
reproduce the engine's own intensity to 1.802e-16. Both are round-off
on quantities that are equal by construction — which is the point: an error
in the projection, the complement, or the channel bookkeeping would show up
here rather than as a plausible-looking floor.

## 5.3 The floor of an uncoated two-mirror train

`Rx_Cass_FarField.in` puts both its mirrors before its single far-field leg,
so the diffraction grid carries the whole train (§5.5) and the number is
physics rather than a bound. With bare perfect-conductor mirrors the floor is
purely geometric — the NA-driven cross-polarization of §4 — and comes to
**7.0612e-07** of the co-polarized power, in a
peak-normalized mean contrast of **9.593e-12** over a 10–40 pixel
annulus.

That number arrives here through a completely different path from §4: from
the propagated grid planes rather than from `RayE`, through a diffraction
leg, with the analyzer derived rather than assumed. It agrees with the
ray-level two-mirror cross-polarized fraction of §4 to the digits printed.

The **longitudinal** channel is not negligible and is reported separately
rather than folded into either transverse one: it carries
2.1065e-03 of the total on this train. Tranche 1's out-of-plane
attribution (§3.5) independently measured the same prescription's
out-of-plane fraction as 1 − f = 2.110e-3, from the pupil rather than the
detector and at a different model size; the two agree to better than a
percent.

![Co-polarized and cross-polarized channels at the detector, and the
cross-polarized contrast curve for three coating choices.](media/polval_2c_channels.png)

## 5.4 The analyzer is derived — and the circular state is what proves it

On a train that does not rotate polarization the mean output state *is* the
input state, so the derived analyzer must track an arbitrary input. It does,
to 1.110e-16 for a circular input.

That gate is not decoration. The coherency matrix is
`C_ij = Σ E_i conj(E_j)`; building it with the conjugation on the wrong
operand yields `conj(C)`, whose dominant eigenvector is the **conjugate**
analyzer. For any *linear* input state the analyzer is real and the two are
identical, so x, y and 45° cases all pass vacuously — while for a circular
state the conjugate analyzer is exactly *orthogonal* to the truth. Projecting
the same detector field on it reports a cross-to-co ratio of
7.130e+05 against the true 1.403e-06. That error was live in
the first draft of this function and this is the gate that found it, which is
why the circular state is in the CI test and not only here.

## 5.5 Coating sensitivity

The §2c deliverable is the floor *and* how it moves with coating choice. On
the two-mirror fixture, coating both mirrors with 200 nm of
aluminium raises the cross-polarized power by **27.9×** over the
uncoated baseline, taking the annulus mean contrast from 9.593e-12 to
**3.345e-10**. Adding a 110 nm MgF₂ overcoat raises it by
**151.3×** of baseline, to **1.992e-09**.

So on this train the coating, not the geometry, sets the floor: the bare
metal costs about a decade and a half over the geometric term, and this
overcoat costs most of another. That is the shape of answer the coronagraph
design rules ask for — the floor is a coating decision before it is an
optical-prescription decision.

> **Corrected 2026-07-28 — read §8.3 before quoting these numbers.** This
> subsection previously described the 110 nm MgF₂ film as "a
> quarter wave, an ordinary protected-aluminium recipe". It is not. This
> fixture runs at `Wavelen = 1.0E-06` m, at which 110 nm of
> MgF₂ is **0.607** quarter-waves; the "632.8 nm" in the coating
> constants' own comments is not this fixture's wavelength. The measured
> **151.3×** is arithmetically right for the film as configured —
> §8.3 reproduces it from the engine's Jones pupil and from an independent
> analytic — but the *generalisation* was backwards: a true quarter-wave
> overcoat **suppresses** the floor to 0.0157× of bare
> aluminium rather than raising it. Treat the number above as a datum for one
> specific sub-quarter-wave film, not as a protected-aluminium design rule.

## 5.6 What this cannot yet measure — the Tranche-1 gap on a real chain

Tranche 1 seeds the three component planes from `RayE` at the **first**
physical-optics leg of a trace and thereafter advances them with a common
scalar phase. A polarizing surface *after* that leg therefore transforms the
rays but not the diffraction grid. `Rx_Cass_FarField` is unaffected — both
its mirrors precede its only leg — but `Rx_Coro.in`, a coaxial chain of
7 reflectors with propagation legs interleaved between them,
is not.

This is measured rather than assumed. The function computes the pupil
cross-polarized fraction twice, once from the grid planes and once from
`RayE` at the same element, and reports the ratio per input state. On
`Rx_Cass_FarField` it is 1.000000 against a ray-level fraction of
7.0612e-07 — full carry, and non-vacuous, since there is real
cross-polarized light to fail to carry. On `Rx_Coro` the grid carries
4.7897e-09 against the rays' 5.6940e-09, a carried fraction of
**0.8412**, and a `macos:pol_contrast_floor:tranche1` warning fires.

Coating the chain makes it worse in the way that matters most. With all
7 mirrors aluminised the carried fraction falls to
**0.5653**, and the reported coating sensitivity comes out at
-0.0322× baseline — while the ray-level cross-polarized fraction
rises from 5.6940e-09 to 9.0336e-09, i.e. +0.587×. The
grid does not merely understate the coating's effect on this chain; it gets
the **sign** wrong, because only the first mirror's coating precedes the seed
leg. Closing this is Phase 3a Tranche 2 (per-ray running Jones applied to the
grid field), and the two tests that pin these numbers are written so that
Tranche 2 has to come back and change them.

## 5.7 The coronagraph number, and how to read it

With that stated, the chain floor itself. `Rx_Coro.in` declares
`nGridpts=511` and is run at model size 512 — below that the engine resets the
grid and then intermittently faults, so the fixture is not usable at 128
despite appearing to run there. With an x-polarized input the focal-plane
cross-polarized channel peaks at **1.2705e-09** peak-normalized
contrast, and over a 20–80 pixel annulus averages **5.787e-13**
(median 1.435e-13) against a co-polarized 4.164e-05 — a
cross-polarized power fraction of 4.7897e-09. Parseval and energy
closure hold at model 512 exactly as at 256 (2.682e-20,
2.197e-16), and the pupil remains fully polarized to
0.99999999, so the analyzer is well conditioned.

**Read that as a lower bound on the uncoated chain, not as the floor of a
built coronagraph.** Six of its seven mirrors sit past the seed leg (§5.6),
and the model's mirrors are all at normal incidence and uncoated, so the
number reflects the numerical aperture alone. The machinery is validated; the
physics it is allowed to see on *this* chain is not yet the whole train.

![Co-polarized and cross-polarized focal-plane channels on the coronagraph
chain, and their contrast curves.](media/polval_2c_coro.png)
