<!-- GENERATED FILE -- do not edit.
     Source: 80_radiometric.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# Coated-refractor transmission radiometry

Every gate up to this point that involved a coating involved a *mirror*.
This section is about the other half of the coating machinery — the
transmitted field through a coated refracting surface — which until this
landing had **no analytic gate at all**, and was wrong.

The defect was not a Fresnel error. It was a disagreement between two code
paths about what the transmitted amplitude *means*. The uncoated branch of
`Refractor` multiplies its Fresnel coefficients by the radiometric factor
√(n₂cos θ₂ / n₁cos θ₁), so that `|TP|²` is the **power** transmittance: the
amplitude carries the ray's power. The coated branch composed plain Fresnel
**field** coefficients through the Airy recursion and applied no such factor.
The two therefore returned different answers for a surface that is
physically the same surface.

The sharpest way to see that is to put a coating on a surface that is not
there. A layer whose index equals the substrate's is optically a bare
interface — the inner boundary has zero reflectance and the layer only adds a
propagation phase — so a coated run must reproduce the uncoated one exactly.
Measured against the pre-fix engine it did not:

| Case | coated / uncoated amplitude, PRE-fix |
|---|---|
| normal incidence | 0.8164965809 = 1/√1.5 |
| 45°, p | 0.7311104457 |
| 45°, s | 0.7311104457 |
| detector-plane intensity | 0.6666666667 = 1/1.5 |

A coated lens under-transmitted by about 18% in amplitude relative to the
identical surface uncoated, and the diffraction grid inherited it. The two
oblique numbers being *equal* is itself informative: a Fresnel error would
separate s from p, whereas a missing common scalar cannot.

## 7.1 The convention, and why it is the incumbent one

The convention kept is the uncoated branch's — `|TP|² = T`, the power
transmittance — and the coated branch was brought to it, rather than the
reverse. Three reasons, in order of weight
(`REVIEW_POL_RADIOMETRIC_2026-07-28.md`):

1. The **internal inconsistency is the bug**, independent of which
   convention is preferable in the abstract. An index-matched layer must
   reproduce a bare interface.
2. The uncoated convention is what **every existing polarized result and
   gate on this page already sits on**. Changing it would churn gated
   behaviour in order to repair a branch that had no gates.
3. It is **self-consistent air-to-air**: across a full transit the per-face
   factors compose to the total power transmittance, which is what the
   detector-side bookkeeping wants (§7.4).

Before comparing anything to the coated branch, the report first pins that
the uncoated branch really does implement that convention. Against the
Abeles characteristic-matrix transmittance — written from Macleod ch. 2,
not from the engine — the uncoated branch agrees to **2.220e-16**
at normal incidence, **2.220e-16** for p and **0.000e+00** for s at
45°. At that angle the two polarizations are genuinely split
(`0.0835` in absolute transmittance), so the oblique agreement is
not two numbers that happen to coincide.

The fix itself is one factor applied **once**, after the recursion, to `TP`
and `TS`:

`TP ← TP · √( Re(n_sub)·cos θ_sub / (n_a·cos θ_inc) )`

with `n_a` the medium the ray is actually travelling in. Never per
interface: the radiometric conversion in the multilayer theorem is
boundary-to-boundary — `T = (n_sub cos θ_sub)/(n_inc cos θ_inc)·|t|²` — and
the interior-layer factors cancel identically. The engine comment carries
the full argument, including why the observation that a plain chain's
per-interface factors happen to telescope is not a licence to write them
that way.

Transmission into an **absorbing substrate** is out of scope and
deliberately ungated: the transmitted wave is evanescent and carries no
propagating power, so no scalar factor is the right one. The use case is
coated lenses, with dielectric substrates.

## 7.2 An index-matched layer is a bare interface

The headline gate, post-fix, on both fixtures and both polarizations:

| Case | coated / uncoated amplitude − 1 |
|---|---|
| normal incidence | 0.000e+00 |
| 45°, p | 2.220e-16 |
| 45°, s | 2.220e-16 |

The oblique pair is not redundant. At normal incidence both cosines are 1
and the factor collapses to √(n_sub/n_a), so a normal-incidence-only gate
exercises the *index* half of the factor and is completely blind to the
`cos θ_sub / cos θ_inc` half. `Rx_Refract45.in` exists for that reason, and
its s/p decomposition is exact by construction — the element is tilted while
the beam stays along the source axis, so the plane of incidence is fixed and
an x-polarized source is pure p, a y-polarized source pure s. The tests
verify that claim from the engine's own ray directions and surface normal
rather than from the prescription's declared numbers.

The same claim holds **on the diffraction grid**: detector-plane integrated
intensity, coated versus uncoated, agrees to **0.000e+00**. That is not a
formality. `CPROPAGATE` re-traces the seed rays through its own element
dispatch chain, so a transmittance can be right at ray level and wrong in
the image — the finding that section 6 was built around. Pre-fix this ratio
was 0.6666666667, so the grid under-reported flux by a third.

## 7.3 A real coating against the textbook

An index-matched layer is a strong self-consistency check but a weak physics
one: it has no interference structure. A quarter-wave MgF₂ layer on glass
does. Against the characteristic-matrix `T`:

| Case | relative residual |
|---|---|
| normal incidence | 2.220e-16 |
| 45°, p | 1.332e-15 |
| 45°, s | 1.110e-16 |

and the stack does what an AR coating is supposed to do — transmittance
**0.9859** against 0.96 for bare glass — so the agreement is not
a stack that silently did nothing.

## 7.4 Air-to-air closure on a coated plate

The air/glass/air parallel plate is where the convention has to compose.
With the **front face coated and the back face bare**, the air-to-air power
transmittance matches the product of the two textbook single-face
transmittances to **2.220e-16**.

The mixed plate is the deliberate choice. With *both* faces coated the two
radiometric factors are √(n_g cos θ_g / cos θ_i) and √(cos θ_o / n_g cos θ_g),
whose product is √(cos θ_o / cos θ_i) = 1 — so a fully coated plate is
air-to-air **invariant under this landing** and cannot detect the defect at
all. It was verified to pass against the pre-fix engine. Mixing a coated
face with an uncoated one breaks the cancellation and the gate becomes
discriminating.

That invariance is still worth measuring, as the composition identity it is:
the air-to-air amplitude through the fully coated plate equals the bare
product of the two composed Fresnel field coefficients, with no radiometric
factor left over anywhere, to **2.220e-16**. This is decision
ground 3 above, stated as a number rather than an argument — and it is
labelled here so that its green is not later mistaken for coverage of the
fix.

## 7.5 What the landing did not move

The factor is a **common real scalar** on `TP` and `TS`, so it cancels
identically in the ratio `t_p/t_s`. The transmitted polarization *state* —
and with it every diattenuation and retardance quantity a Jones-pupil
analysis reads off a coated refractor — is untouched; only the absolute
amplitude changed. Measured against the ratio of the two composed Fresnel
field coefficients, computed with no radiometric factor at all: agreement to
**4.441e-16**, on a ratio that is `0.0213` away from unity.

For the same reason nothing **reflected** can move: `RP`, `RS` and the
`RP1·RP` Airy denominators are untouched by the change, so a coated
surface's reflected field is bit-identical across it. That is what allows
the Phase 2a/2b coated results earlier on this page to stand unamended.

Finally, the whole coated branch is gated on `ifPol`. With polarization off,
a coated and an uncoated run of the same prescription agree **bitwise** in
OPD and in detector intensity (`1`, 1 = identical), and the GMI
regression — which carries no polarization at all — is
6/6, bit-identical (vs-ref = 0.000e+00).

## 7.6 The interference structure is untouched

A scalar multiplier cannot create, shift or flatten a spectral feature, and
this is the cheapest way to confirm that nothing else was disturbed.
Sweeping the wavelength across the quarter-wave design point with the
**physical** stack held fixed — `coat_set` stores physical thickness, so the
sweep is genuinely chromatic — the engine tracks the characteristic-matrix
`T` to **5.551e-16**, and the ratio of the two is constant to
**9.992e-16** across the sweep: whatever the factor is, it is
wavelength-independent, as a real scalar must be.

The structure being preserved is real, not flat: the transmittance varies by
**0.0195** across the sweep and peaks at the design wavelength,
where bare glass would be a featureless 0.96.

## 7.7 Non-vacuity

Rebuilt against the pre-fix engine, `tPolRadiometric` scores **6 pass, 7
fail**. The six that pass are the six that should: §7.1's two
uncoated-convention gates, the 45° fixture's s/p purity check, §7.4's
telescoping identity, §7.5's polarization-state invariance, and the
polarization-off gate. Each of those says so in its own comment, so a future
reader is not left to infer which greens are coverage and which are
controls.
