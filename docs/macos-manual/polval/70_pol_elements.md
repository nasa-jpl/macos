<!-- GENERATED FILE -- do not edit.
     Source: 70_pol_elements.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# Polarizing elements — polarizer and waveplate

Everything measured up to this point describes polarization that a system
*acquires*: the diattenuation and retardance a metal coating imposes at
oblique incidence, and how those propagate to a detector. This section
covers the two elements that impose polarization *deliberately* — an ideal
linear polarizer (`TrPolarizer`) and a linear retarder (`WavePlate`) — which
are what a polarization phase-shifting interferometer is built from, and
what a circular analyzer for a vector vortex coronagraph would be built from.

Both existed only as names before this work: `RfPolarizerElt` and
`TrPolarizerElt` were entries in the element-type table with no trace
dispatch anywhere, so a prescription naming one loaded without complaint and
did nothing at all. `WavePlate` is new. The per-element Jones array
`JmatElt` was allocated, zeroed and deallocated by the engine and never
otherwise touched; it now holds what the element applied.

Geometrically these are `Reference` surfaces — the same conic intersection,
the same path-length accounting, the ray direction unchanged. What they add
is a 2×2 Jones matrix applied to the ray's field in the plane transverse to
the ray, gated on `ifPol`. That equivalence is the first thing measured
below, because it is what makes the elements safe to leave in a prescription
that is being used for scalar work.

## 6.1 With polarization off, they are Reference surfaces

The gate fixture `Rx_PolElt.in` has a twin, `Rx_PolElt_Ref.in`, identical in
every respect except that its four polarizing elements are plain `Reference`
surfaces. Traced with polarization off, the two must agree **bitwise** — not
to a tolerance — in OPD, in detector intensity, and in per-ray status.

Measured: `1` (1 = bitwise identical on all three).

This is deliberately a comparison against a *separate prescription* rather
than an inspection of the source. The claim being tested is that the new code
path does not perturb the geometry, and a twin fixture is the only way to
test it that does not depend on reading `PolElt` and believing it.

## 6.2 Malus's law

With the input polarizer along *x* and the analyzer at θ, the projection onto
the analyzer axis is cos θ and the transmitted intensity is `I₀cos²θ`.
Measured across 25 angles spanning half a turn, the largest deviation from
the closed form is **1.079e-13** of `I₀` (left panel, Figure 7).

The companion number matters as much: the ratio of the smallest to the
largest transmitted intensity across the sweep is **3.749e-33**. An
element that had quietly failed to do anything — the pre-existing behaviour
of this element type — would produce a *flat* curve, which fits `I₀cos²θ`
nowhere but would sail through a test that only checked a mean or a
correlation. The dynamic range is what makes the fit non-vacuous.

At exactly crossed axes the extinction is not merely small, it is
**exact**: with two axes that are orthogonal as floating-point vectors, the
projection of a field lying entirely along one onto the other is identically
zero, and the measured transmitted power is `0` — bitwise
(`1`, 1 = exactly zero). Anything else would mean the
transverse basis is not orthonormal.

## 6.3 The retarder, and the sign of the retardance

The engine propagates a field through a medium as `exp(−i·2πLN/λ)`, i.e.
with `exp(+iωt)` time dependence. The slow axis therefore accumulates the
*more negative* phase, so with the declared axis taken as fast, the element's
Jones matrix in its own eigenbasis is `diag(1, exp(−iδ))`, `δ = 2πR`. That is
read off the same line of `elemsub.F` the conventions table of section 1
already cites — it is not a separate choice.

For linear input at 45° to the fast axis, the algebra gives
`Ex = (1 + e^{−iδ})/2`, `Ey = (1 − e^{−iδ})/2`, hence

```
S1/S0 =  cos 2πR        S2/S0 = 0        S3/S0 = -sin 2πR
```

The middle panel of Figure 7 is the measured sweep against those two curves.
At the quarter-wave point the state is exactly circular:
`|S1|, |S2| ≤ 1.741e-16` and — the number that carries the convention —
the **signed** `S3/S0` sits at −1 to within **0.000e+00**.

The sign is the whole point. A gate written on `|S3|` would be satisfied by
either convention, and would therefore certify a retarder whose fast and slow
axes were swapped. This is the same lesson the coherency matrix taught in
section 5, where a conjugation-order slip was invisible to every linear input
state and only a circular one exposed it.

The A/B is in the same rig: with the retardance set to zero, the emerging
state stays linear and `|S3|/S0` collapses to **0.000e+00**.

## 6.4 The 2θ law and composition

A half-wave plate applies `diag(1,−1)` in its eigenbasis, which is a
*reflection* of the transverse plane about the fast axis. Linear input at 0°
therefore emerges at 2θ when the plate is at θ. Fitting the measured output
orientation against plate angle over a quarter turn gives a slope of 2 to
within **8.882e-16** (right panel, Figure 7). With the retardance
set to zero the same sweep gives slope **9.205e-19** — rotating a
plate that does nothing does nothing.

Orientation is taken from `½·atan2(S2, S1)` rather than from the ratio of
field components, because the hop between elements applies a common complex
propagation phase that a real-part ratio would fold into the answer. The
Stokes pair is quadratic in the field, so the common phase cancels
identically.

Composition is the gate that the elements *cascade* correctly, which no
single-element test can see: two quarter-wave plates on a common axis must
reproduce one half-wave plate, since `diag(1, e^{−iπ/2})² = diag(1, e^{−iπ})`.
Field for field, the difference is **1.464e-16** of peak amplitude. And
the pair is not trivially equal to anything — a *single* quarter-wave plate
at the same axis differs by **3.078e-01**.

## 6.5 Unitarity

A lossless retarder conserves power for every input state. Measured across
linear input at three relative orientations and a circular input (built by
making the first plate a quarter-wave plate at 45°, since an ideal polarizer
cannot deliver a circular state through it), the power ratio across the plate
departs from unity by at most **0.000e+00**, and the element's own
reported Jones satisfies `JᴴJ = I` to **0.000e+00**.

The non-vacuity companion is the polarizer: its Jones is checked to *fail*
the same condition, `‖JᴴJ − I‖ = 1.000e+00`. Without that, a unitarity
test would also be passed by an engine that had stopped applying anything.

## 6.6 The diffraction grid sees the elements — and nearly didn't

The detector-plane intensity obeys Malus as well as the rays do: across four
analyzer angles the plane integral tracks `I₀cos²θ` to **1.833e-15**.

That gate exists because of a defect it caught during development, which is
worth recording as a general hazard rather than an incident. MACOS has **two**
element dispatch chains: `tracesub.F` traces rays for ray-level queries, and
`propsub.F`'s `CPROPAGATE` *re-traces* the rays that seed the diffraction grid
through a separate chain of its own. An element wired into only the first one
satisfies every ray-level gate in this section — Malus, extinction, the
quarter- and half-wave identities, composition, unitarity — and is completely
invisible to `intensity` and `complex_field`. Measured in that state: crossed
polarizers took the ray power to 3.6e-33 while the detector
plane held the full incident 9.69e-01 — the same value it
reports with the analyzer *aligned*, which is the point *(external, captured
2026-07-27; the wiring that produced it no longer exists
in the tree)*.

The failure mode is not a crash and not an obviously wrong number. It reads as
*"polarization has no effect on the image"* — a conclusion a reader of this
report might well have accepted. Both chains are now wired, and this gate is
the tripwire for the next element that touches the field.

Note also the fixture constraint behind the measurement. All four polarizing
elements in `Rx_PolElt.in` sit **before** its single physical-optics leg. That
is the Tranche-1 validity condition of section 3: the grid is seeded from the
ray field at the first physical leg, so a polarizing element placed after that
leg would transform rays and never reach the grid — the same limitation
quantified for the coronagraph chain in section 5.

## 6.7 Scope: normal incidence, and the off-normal convention (now settled)

Every gate above is at strict normal incidence. That is not a convenience.

Away from normal incidence an *ideal* polarizer carries a genuine convention
question. The element's axis is declared in global coordinates and must be
projected into the plane transverse to the ray, and orthographic projection
does not preserve orthogonality — so declaring the **pass** axis and
projecting it is not the same operation as declaring the **block** axis,
projecting that, and transmitting the complement. The two constructions
disagree by **3.56°** of axis orientation at 20° of incidence
(at 45° of azimuth, where the effect is largest; the closed form is
`acos(2cos a/(1+cos²a))`, bounded by `sin²a`).

**The convention is settled (2026-07-27): project the *material* axis.**
The invariant object in a tilted element is the material direction fixed in
it — the absorbing direction (the wires of a wire grid, the dipole chains of
a dichroic sheet) for a polarizer, the crystal fast axis for a waveplate.
This is the *dipole model*, and it is not a modelling taste: exactly this
ambiguity was adjudicated experimentally for a tilted dichroic polarizer by
Korger et al., *Opt. Express* **21**, 27032 (2013), in the dipole model's
favour. For the waveplate the declared (fast) axis *is* the material axis,
so the implementation above is already the settled rule. For the polarizer
the material axis is the block direction — the in-plane complement of the
declared pass axis — and the implementation currently projects the pass axis
instead: **identical at the normal incidence gated here, pending the flip
off normal**, which lands with an off-normal fixture and an engine-driven
gate against the closed form above. Until then the off-normal polarizer
path is unvalidated as well as un-flipped, and this section's numbers are
unaffected either way.

The disagreement vanishes identically at normal incidence, which is why every
number in this section is unaffected by the choice. It **also** vanishes when
the declared axis lies in, or perpendicular to, the plane of incidence:
`0.000e+00°` at the same 20° tilt. That second zero is a trap
worth naming — the most natural way to probe the effect (axis along *x*, tilt
in the *x*–*z* plane) is exactly the degenerate case, and reports a reassuring
zero that means nothing.

`RfPolarizer` — a *reflective* polarizer, which is inherently off-normal —
is still not implemented. The settled convention no longer blocks it, but a
reflective wire grid carries additional physics (grid reflection efficiency,
the substrate's own s/p response) beyond the axis rule, and nothing requires
it yet; it fails loudly at load rather than tracing as absent. For a
beamsplitter at a substantial angle, a coated `Reflector` is in
any case the better model: there the s and p directions are set by the
physical plane of incidence and the thin-film recursion produces real
diattenuation and retardance, gated against the Fresnel closed form in
section 2.

The remaining idealizations are stated rather than measured, since they are
definitions rather than approximations to be bounded: no ray splitting (a
polarizing beamsplitter is two traces), no walk-off, no Fresnel loss at the
plate faces, no substrate thickness. The output field is purely transverse,
so any longitudinal component present at the surface is discarded — which is
what a 2×2 Jones element means, and which is exactly zero for the collimated
normal-incidence fixture used here.

![Figure 7 — Polarizing elements against their closed forms. Left: Malus's
law across 25 analyzer angles. Centre: the Stokes parameters of linear input
through a retarder whose fast axis is at 45°, swept over a full wave of
retardance, against `cos 2πR` and `−sin 2πR`; the signed `S3` curve is what
pins the retardance convention. Right: the half-wave plate's 2θ rotation
law.](media/polelt_gates.png)
