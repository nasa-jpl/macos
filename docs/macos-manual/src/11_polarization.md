## SECTION 11: Polarization

MACOS carries a complete vector-polarization physics: every ray can carry
a complex three-vector electric field, every reflector and refractor
applies its Fresnel or multilayer-coating coefficients to it at the ray's
own angle of incidence, two element types polarize and retard the field,
and the diffraction propagations can carry the three field components
through the whole chain.  Above the engine, the bindings assemble the
Jones pupil and the standard polarization-aberration metrics and the
polarization-limited contrast floor of a coronagraph.  This section is
the user's description of that machinery; every claim in it is pinned by
a named test, and the measurements behind the claims are in the MACOS
Polarization Validation Report (Section 11.8).

### Turning polarization on

Polarization ray tracing is off by default and is a *session mode*: once
turned on it survives every LOAD until it is turned off, and every LOAD
prints a line naming the non-default modes in force.

    MACOS>pol
     Polarization turned on
     Vector diffraction enabled
    Enter Ex0 (real,imag): [0.,0.]: 1,0
    Enter Ey0 (real,imag): [0.,0.]: 0,0

The POLarization command prompts for the complex source field amplitudes
Ex0 and Ey0, in the source grid frame; `POL 1,0 0,0` on one line is the
same.  A prescription can carry the state in the `Ex0Ey0=` keyword
(Section 4) — the state only; turning the trace on remains a command or
API decision.  NOPOLARIZATION turns it off.  On every stock model size the
command also enables vector diffraction; SCALAR and VECTOR (Section 6)
toggle that part separately.  A 0,0 launch is accepted silently and
traces no field, so set the state before reading any polarization result.

From the bindings the same switch is `polarization('on', Ex0, Ey0)` in
mmacos and pymacos (api `pol_set` / `pol_get`), and `vecdif_set` for the
vector-diffraction part.  Changing the polarization state marks the cached
trace dirty, so the next query re-traces with the new state.

With polarization on, the geometry is unchanged: ray positions,
directions and path lengths — and therefore every OPD and spot result —
are bit-identical to the scalar trace.  Polarization adds the field; it
does not move the rays.

### Conventions

Cross-code polarization disagreements are nearly always convention
mismatches.  MACOS's conventions were read off the engine source, are
asserted by the tests, and are the authority for every number in this
section and in the validation report.

| Convention | Value |
|---|---|
| Time dependence | exp(+iωt); a field accumulates exp(−i·2π·L·N/λ) through a medium of complex index N over a path L, so the phase decreases as the path grows |
| Absorbing index | N = n − iκ with κ > 0 meaning loss, as entered in `IndRef=` / `Extinc=` and in `Coating=` layers |
| Jones storage basis | linear (x, y); circular states by a unitary change of basis |
| Jones pupil reference frame | double-pole (Chipman) by default; local s/p and global available as diagnostics |
| Retardance sign | fast axis leads; a retarder's Jones matrix in its own eigenbasis is diag(1, exp(−iδ)) |
| Diattenuation | intensity-based, D = (T_max − T_min)/(T_max + T_min) in [0, 1] |
| Coating thickness and retardance in a prescription | waves at the `Wavelen=` in force when the parser reaches the line; converted to a physical length at load, so a film or a plate is fixed glass and a wavelength sweep is chromatic |
| Field amplitude normalization | the uncoated and the coated transmission both carry the radiometric factor √(n₂ cos θ₂ / n₁ cos θ₁), so \|t\|² is the power transmittance |

Two consequences of the index convention are worth stating.  A mirror
described by the perfect-conductor idiom (`IndRef= 1.0`,
`Extinc= 1.0E+22`) reflects with r_s = r_p = −1 in a fixed transverse
basis and is polarization-neutral: it changes the ray geometry and
introduces neither diattenuation nor retardance, so any stock
perfect-conductor prescription is a unitarity check for free.  A real
metal mirror is modeled with a metal coating stack; an optically thick
metal layer reproduces the bare metal's Fresnel coefficients whatever the
substrate.

### What the trace does with the field

**Launch.**  A collimated source launches every ray with
E = Ex0·xGrid + Ey0·yGrid, uniform across the beam.  A point source
launches each ray in its own frame, y = unit(RayDir × xGrid),
x = y × RayDir, which reduces to the grid frame on the chief ray; the
difference from a uniform frame is of order NA²/2 and is measured, not
assumed.  The flux normalization is a common real scalar that every
Jones quantity carries and no ratio sees.

**At a bare surface.**  The field is decomposed into s and p components
in the local plane of incidence, the surface's complex Fresnel
coefficients for the ray's angle of incidence are applied, and the
field is reassembled on the s and the *reflected* p directions, so the
field stays transverse to the outgoing ray.  The incident medium is the
one the ray is actually traveling in; the other medium is the element's
own complex index.

**At a coated surface.**  An element carrying a `Coating=` stack
(Section 4) is evaluated by the recursive thin-film (Airy) calculation
through every layer at the ray's own angle of incidence, with the stack's
complex indices, the incident medium the ray is in and the element's
index as the substrate.  Reflectors return the reflected field; refractors
the transmitted one.  The transmitted amplitude carries the radiometric
factor once, boundary to boundary, so a coated lens and an uncoated lens
agree on what an amplitude means.  Transmission into a metal substrate
is evanescent and carries no power; coated *lenses* are the use case.
MACOS has a second, older coating description (`nCoatElt` and its
companions, real indices only) that serves the non-sequential refractive
trace and does not enter the polarization path; the `Coating=` stack is
the one that does.

**Two lessons from the engine's history** are kept here because they are
the shape of the errors a reader is most likely to make in another code.
A reflection's p-coefficient sign depends on whether the p direction is
held fixed or follows the outgoing ray; the wrong pairing reflects the
transverse field about the p axis instead of negating it, which is
unitary and an involution — it cancels exactly across any mirror pair,
so every two-mirror and unitarity test passes while one near-normal
mirror turns x-polarized light into equal x and y.  The tell was a
cross-polarization *flat in pupil radius* where an isotropic surface
must give a ρ² law.  And a transmission amplitude's radiometric factor
must be applied once for the whole stack, never per interface: the
interior factors cancel identically, and a per-interface form is right
only while every interface is present, in order, with the right media.

### Polarizing elements

`TrPolarizer`, an ideal linear polarizer, and `WavePlate`, a linear
retarder, are described in Section 4 (the `PolAxis=` and `Retardance=`
keywords, with examples).  Geometrically they are Reference surfaces:
with polarization off they *are* Reference surfaces, and a prescription
with them in it cannot change a scalar result.  With polarization on each
applies a 2×2 Jones matrix in the plane transverse to the ray —
diag(1, 0) for the polarizer, diag(1, exp(−iδ)) for the retarder, in an
eigenbasis built per ray.

The basis is the physical content.  Away from normal incidence the
declared axis has to be brought into the plane transverse to the ray, and
orthographic projection does not preserve orthogonality, so it matters
which axis is projected.  MACOS projects the direction fixed in the
*material* of the element — the absorbing direction of a polarizer (the
wires of a wire grid), the crystal fast axis of a waveplate — and
completes the partner axis by orthonormalization, which is what keeps a
lossless retarder unitary.  The same question was settled experimentally
for a tilted dichroic polarizer in the dipole model's favor (Korger et
al., Opt. Express 21, 27032, 2013).  The effect is 3.6° of transmitted
axis at 20° incidence and 45° azimuth, bounded by sin² of the incidence
angle, and exactly zero at normal incidence and whenever the axis lies in
or perpendicular to the plane of incidence — so the obvious test, axis
along x with the tilt in the x–z plane, reports a reassuring zero that
proves nothing.  `PolAxis=` still declares a polarizer's transmission
axis; the absorbing direction is derived from it.

The idealizations: no ray splitting (a polarizing beamsplitter is two
traces, or better a coated Reflector at 45° where the s/p split is
physical), no walk-off, no face Fresnel loss, no substrate, a purely
transverse output, and a retardance that does not vary with incidence
angle — a real crystal plate's does, which is the field-of-view effect
behind compound and Pancharatnam designs.  A reflective polarizer
(`RfPolarizer`) is a named type that is not implemented.

### Vector diffraction

With polarization on and VECTOR in force, every propagation leg carries
the three Cartesian field components Ex, Ey, Ez as three independent
scalar fields through the same kernel: the far-field sphere-to-plane leg,
every near-field, plane-to-plane, surface-to-focus and DFT leg.  Scalar
mode runs the same code on one plane and is bit-identical to the engine
before vector mode existed.  The grid is seeded from the per-ray field at
the first physical-optics leg and advanced through the following legs;
rays the trace has clipped are extinguished at the seed, so polarization
on with vector diffraction off is bit-identical to polarization off.  A
0/1 mask is a diagonal Jones matrix and is applied to all three planes.

Three things to know when reading a vector-mode field:

- The three wavefront storage planes are repurposed as Ex, Ey, Ez of
  one wavefront, so vector mode and multi-wavefront composition
  (COMPOSE and the broadband commands) cannot run together.
- The components are reachable from the bindings:
  `complex_field(elt, 'plane', k)` returns Ex, Ey or Ez for k = 1, 2, 3
  and refuses the request in scalar mode, where a plane is an unrelated
  wavefront and would read as a component.  The vector intensity is the
  sum of the three; on an off-normal train it differs from the scalar
  intensity by the out-of-plane power and by the separate diffraction of
  Ey and Ez (0.26 % on the Cassegrain example, decomposed and pinned).
- The field's common phase is train-dependent relative to the OPD map.
  MACOS bookkeeps the OPD by subtracting return legs and the per-ray
  field phase by adding every leg, so on a Return-terminated train the
  assembled field matches the scalar convention, while a plain trace to
  a detector measures the exact conjugate of the OPD map.  There is no
  universal bridge; the Jones-pupil layer divides the common scalar out
  and never needs one, and anyone comparing a vector-mode field to a
  scalar one phase for phase must expect the conjugate on such a train.

The one limitation that remains: between two physical-optics legs the
grid field is advanced by a scalar phase, which is exact only when the
intervening elements are non-polarizing (Obscuring, Reference,
FocalPlane).  A coated or reflecting surface *between* two legs — folds
between an interferometer's recombiner and its detector, an off-axis
parabola between two masks — is not yet carried by the grid; the
polarizing elements have the same restriction, being applied to rays
before the seed.  Results on such a chain are lower bounds, and the
contrast-floor function measures the shortfall and says so.

### Polarization analysis

**Polarization analysis.** The per-ray field at any element can be
harvested through the bindings (`ray_field`), and the standard
polarization-aberration analysis is built on top of it: `jones_pupil`
assembles the 2×2 Jones matrix at every pupil point from two traces
with orthogonal source states, and `pol_maps` decomposes it into
diattenuation and retardance maps (with the pupil mean — a state
change — separated from the spatially varying part that drives error
budgets). `pol_zernike` then expands those maps onto a Zernike basis,
giving the standard polarization-aberration terms — piston, tilt,
defocus, astigmatism and up — so that a MACOS result can be compared
with the published literature term by term rather than by map shape.
For an on-axis rotationally symmetric two-mirror system the expansion
reduces, as the theory requires, to pure *polarization astigmatism*:
diattenuation and retardance grow as the square of the pupil radius
with an axis locked to the pupil azimuth, and every other term —
including the entire circular component — sits at round-off. These
functions and their conventions, including the choice of pupil
reference basis (double-pole by default), are documented in the MACOS
Command Reference, Part II. The interactive CLI exposes the underlying
trace and coating machinery; the Jones-pupil analysis layer is
binding-only.

**Polarization contrast floor.** For a coronagraph the question is not
the shape of the aberration but how much light it puts where the dark
zone should be. `pol_contrast_floor` answers it directly: it propagates
with vector diffraction on and splits the detector field into a
co-polarized, a cross-polarized and a longitudinal channel. The
cross-polarized channel is the part no scalar deformable-mirror control
can remove, so its peak-normalized level is the polarization-limited
contrast floor, and the function also reports how that floor moves with
the coating choice.

Two points of method matter. First, the split is taken at the detector
rather than in the pupil. Because the chain is linear in the input
polarization state and all three components propagate with the same
kernel, a spatially uniform analyzer commutes with propagation, so
projecting after the propagation gives the same answer as projecting
before it — and it avoids the Jones pupil, which cannot serve as a
pupil multiplier because it carries the accumulated optical path
length. Second, "co-polarized" is defined against the *mean output*
state, not the input state. A train may rotate the polarization
geometrically while introducing neither diattenuation nor retardance;
charging that rotation to the cross-polarized channel would report an
aberration where there is none, since an observer would simply align
the analyzer to the output. The analyzer used is therefore the dominant
eigenvector of the pupil coherency matrix, which is insensitive to the
common wavefront and by construction minimizes cross-polarized power.

An unpolarized source is modeled as two traces with orthogonal input
states summed in intensity, never as one trace with the second state
inferred from it. The floor reported on a chain that places coated or
reflecting surfaces between two propagation legs is a lower bound, for
the reason given in Section 6; the function measures the shortfall and
warns rather than leaving it to the reader.

The sharpest cheap check of a polarization set-up is the Zernike
expansion on an on-axis rotationally symmetric two-mirror train: the
theory requires pure polarization astigmatism — astigmatism-0 in the
first Pauli component, astigmatism-45 in the second, equal in magnitude,
no circular part, no defocus, a ρ² radial law whose on-axis extrapolation
vanishes.  MACOS puts every forbidden term at 1e-15 of the astigmatic
one.  A broken pupil frame or a wrong s/p sign breaks the mode *pattern*,
not merely the map, which is why this test sees what magnitude tests
miss.

### Design rules the measurements taught

- **Specify a protective overcoat at λ/4 of the working wavelength.**  An
  overcoat is neither a polarization cost nor a benefit in itself; the
  sign of its effect is its optical thickness at the working wavelength,
  and it reverses across the quarter-wave condition.  Measured on the
  Cassegrain with a 110 nm MgF₂ overcoat on aluminium, cross-polarized
  power against the bare metal: 5.27× at 1 µm, 0.20× at 632.8 nm, a true
  quarter wave 0.053× at either.
- **Report two ratios.**  The total cross-polarized power is what a
  designer sees and carries the irreducible geometric term; the
  coating's excess over the uncoated baseline is the coating-only
  quantity and the only one comparable with a pure-Fresnel analytic.
  They diverge most at the quarter wave, where the coating term is
  nearly extinguished and the total sits at 1.5× the uncoated floor.
- **An unpolarized source is two traces** with orthogonal input states
  summed in intensity, never one trace with the second state inferred.
- **Across wavelengths, measure total power**, never a fixed-pixel
  annulus mean, which subtends a different λ/D range at each wavelength
  and mixes the coating effect with the diffraction scale.
- **Write a gate's analytic from the textbook**, never from the engine's
  own expression: a reference transcribed from the code is circular in
  exactly the sign it should check.

### Validation

The MACOS Polarization Validation Report (`docs/macos-manual/polval/`,
built as `macosPolVal`) measures every claim above with a named test
behind each number and the conventions table reproduced byte for byte
from its own conventions section.  Its ladder, from the bottom:

1. polarization on does not move geometry (OPD bitwise unchanged);
2. the two-trace Jones pupil is unitary on a lossless train
   (diattenuation residual 2.6e-15) and matches closed-form Fresnel
   coefficients on a coated fold to 1.2e-14;
3. the two-mirror train reduces to pure polarization astigmatism
   (forbidden terms at 7e-15 of the astigmatic term);
4. after one near-normal perfect-conductor mirror the whole single-
   reflection Jones matrix matches the Born and Wolf closed form at
   2e-15 median, and the cross-polarization grows as ρ² (slope 1.87);
5. coated transmission matches the Abeles characteristic matrix
   (Macleod, ch. 2) at normal incidence and at 45°, the angle that
   exercises the cosine half of the radiometric factor;
6. the polarizing elements obey Malus's law, the retarder's signed
   circular Stokes parameter, the 2θ composition law and unitarity, on
   both element dispatch chains, at 20° incidence;
7. vector diffraction reproduces the scalar intensity to 4.5e-16 on a
   normal-incidence chain for x, 45° and circular input (a test with
   x-polarized input alone passes vacuously);
8. the coated reflector is anchored to a published protected-aluminium
   mirror model (van Harten, Snik and Keller, PASP 121, 377, 2009):
   diattenuation and retardance agree with their curves to 3e-14 and
   5e-14 against their stated ±0.01, driven with *their* indices — index
   tables disagree, and that is not a machinery error.

Every fix the ladder found is pinned by a test that the pre-fix engine
fails; the report names them.

### Idealizations and gaps

No ray splitting at any element; no Mueller calculus and no partially
coherent light (unpolarized light is two traces); no vector vortex
coronagraph element; `RfPolarizer` not implemented; a waveplate's
retardance independent of incidence angle; between-leg polarizing
surfaces not carried by the grid (the Tranche-2 running Jones matrix per
ray is designed, not built); spatially variable coatings not yet
available; vector mode and multi-wavefront composition exclusive.  The
oblique-incidence radiometric factor and the ray-density bookkeeping of a
diffraction seed overlap conceptually; whether that double-counts in
oblique refractive seeding is an open audit question — the mirror trains
the PROPER comparisons used never exercised it.
