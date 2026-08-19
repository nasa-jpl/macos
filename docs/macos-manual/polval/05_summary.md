<!-- GENERATED FILE -- do not edit.
     Source: 05_summary.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# Executive summary

*For a reviewer with no time. Every number here is the same generated token
as in the body — nothing is typed by hand or rounded for the summary — and
nothing is claimed that the body does not measure. Engine
`aaa9ece` on `pol-core`, bindings
`a89b077` on `pol-core`, model sizes
128 / 256 / 512, generated 2026-07-28 07:59:06.*

## Scope

MACOS has carried vector polarization physics — a per-ray complex E-field,
Fresnel and multilayer coatings with complex index, three-component vector
diffraction — for a long time, largely ungated. This package exposes that
physics through a language-neutral API, adds the Jones-pupil and
polarization-aberration layer above it, and puts a named CI test behind
every claim; the eight result paragraphs below are the eight numbered
sections of the report. Four engine defects were found while building it;
all four are fixed and each is named in *Conventions* with the gate that now
pins it. What is **not** covered — the interferometer track, the vector
vortex, the reflective polarizer, per-ray Jones through a chain (Tranche 2),
spatially variable coatings — is listed under *Coverage and gaps*, and where
it limits a number quoted here it is stated inline, not in a footnote.

## Conventions — review priority 1

Cross-code polarization disagreements are almost always convention
mismatches. These were pinned by *deriving them from the engine source*, not
by legislating them, and they are asserted in the tests. No number in this
report is interpretable without them. Reproduced verbatim from *Conventions*,
which is the authority; the two copies are compared byte for byte at build
time.

<!-- CONVENTIONS-TABLE:BEGIN -->
| Convention | Value | Where it comes from |
|---|---|---|
| Time-harmonic phase | `exp(+iωt)`; spatial propagator `exp(−ikz)` | `elemsub.F:387` (`C1 = exp(−i·2π·L·N/λ)`, phase decreases as path grows) and the coating recursion at `elemsub.F:512-516`. Consistent with the independently established interferometer result that field phase advances as optical path shortens, and with the `opd_sign_flip` the pymacos↔PROPER comparison already carries. |
| Absorbing index | `N = n − iκ`, with κ > 0 = loss | As stored in `IndRefArr`/`ExtincArr` and applied as `DCMPLX(n,−κ)`. |
| Jones storage basis | Linear (x, y); circular by unitary change of basis | Phase-2 decision — storing linear and converting is strictly simpler than the reverse. |
| Jones reference frame | Double-pole (Chipman) by default; local s/p and global as diagnostics | The naive local s/p basis carries a coordinate singularity that appears as spurious retardance. Quantified in the basis-artifact gate below. |
| Retardance sign | Fast axis leads | Stated so that vector-vortex handedness (Phase 3) has a fixed reference. |
| Diattenuation | Intensity-based, `D = (T_max − T_min)/(T_max + T_min)` ∈ [0,1] | |
| Coating thickness | Rx `Coating=` layer thickness is **waves at parse-time `Wavelen`**, converted to physical at load; the trace applies phase at the *current* wavelength | `msmacosio.inc:2660`. Broadband sweeps are therefore already correct. `Coating=` must follow `IndRef=` in the file (the parser snapshots the boundary media). The `coat_set` API takes **physical** thickness and sidesteps all of this. |
<!-- CONVENTIONS-TABLE:END -->

## Results

**Foundation — exposure and the Jones pupil (§1, §2).** Turning polarization
on does not move geometry (OPD bitwise unchanged, 0.000e+00 waves).
The two-trace Jones pupil is unitary on a lossless train (residual
diattenuation 2.61e-15), matches closed-form Fresnel coefficients on a
coated fold to 1.20e-14 in magnitude, and on an
on-axis symmetric two-mirror train reduces to the pure **polarization
astigmatism** the literature requires, with every forbidden Zernike term at
**7.18e-15** of the astigmatic term. That last is the sharpest cheap
check here: a broken pupil frame breaks the mode *pattern*, not just the map.

**§4 — a reflection sign convention, corrected.** `Reflector` assembled the
reflected field on a p̂ that follows the outgoing ray while multiplying the
p-component by −r_p, reflecting the transverse field about p̂ instead of
negating it. The operation is an involution (it cancels exactly across a
mirror pair) and unitary (no unitarity gate can see it) — and every fixture
before §4 has an even number of mirrors or none. After one near-normal
mirror an x-polarized beam came out at Py/Px = 1.0163e+00,
where physics allows O(sin²β); post-fix **2.0724e-04**,
growing as ρ² (slope 1.871) where pre-fix it was *flat in pupil
radius*. Two-mirror results are bit-identical across the fix, so prior
two-mirror evidence stands.

**§7 — coated-transmission radiometry, corrected.** Coated and uncoated
`Refractor` transmission disagreed about what the amplitude *means*: the
uncoated branch carries √(n₂cos θ₂/n₁cos θ₁) so |T|² is the power
transmittance; the coated branch omitted it. A coated lens under-transmitted
by ~18% in amplitude and the grid inherited it: on a layer that is optically
a bare interface, the pre-fix detector intensity came out at
**0.6666666667** = 1/1.5 of the uncoated run it had to equal. One
factor applied once after the Airy recursion — never per
interface — restores the incumbent convention; air-to-air power now closes
on a mixed coated/bare plate to 2.220e-16. Nothing reflected moved and
the transmitted polarization *state* is untouched (4.441e-16), which is
what lets the §2 coated results stand unamended.

**§6 — polarizing elements, and the off-normal axis rule.** An ideal linear
polarizer and a linear retarder are implemented (both were name-table-only
stubs) and are bitwise `Reference` surfaces with polarization off; Malus's
law holds to 1.079e-13 of I₀ and the retarder is unitary to
0.000e+00, with the retardance *sign* pinned by the signed circular
Stokes parameter (a |S₃| gate would accept either). Off normal the declared
axis must be projected into the ray's transverse plane, and projection does
not preserve orthogonality, so what gets projected is the direction fixed in
the element's **substance** — the absorbing direction for a polarizer, the
crystal fast axis for a retarder. Not a taste call: the same ambiguity was
adjudicated experimentally for a tilted dichroic polarizer (§6.7 cites it).
Worth **3.56°** of transmitted-axis orientation at 20° AOI, and
identically zero at normal incidence *and* when the axis lies in or normal
to the plane of incidence — so the obvious test tilt is degenerate and its
reassuring zero means nothing.

**§3 — vector propagation across the whole chain.** All three component
planes now propagate through every leg, not only the far-field one, which is
what unblocks the coronagraph and vortex cases. On the polarization-neutral
gate train the vector chain reproduces the scalar intensity after every leg
for x, 45° and circular input — worst case **6.42e-16** against a truth
of exactly zero. Run at x-polarization only, that gate passes vacuously.
Pre-fix the same gates fail at 0.21 .. 0.38 relative error, and the
validated scalar physics is undisturbed (PROPER residual
4.836e-13 reproduced exactly).

**§8 — anchored against publication.** The coating machinery is compared
curve-on-curve with a published protected-aluminium model by reproducing
*their* configuration through `coat_set` — their indices, their 4.12 nm
Al₂O₃ on 220 nm Al, their wavelengths, their 6–70°: **2.828e-14**
(diattenuation) and **4.937e-14** (retardance) against their stated
±0.01 per normalized Mueller element. Using their indices is the
point — aluminium index tables genuinely disagree, and an index-table
difference is not a machinery error.

**§8.3.1 — the overcoat trade reverses across the quarter wave.** The anchor
exposed that the MgF₂ overcoat in §5.5's ladder is 0.6072
quarter-waves at that fixture's own wavelength, not the quarter wave its
source comment claims. Both sides are now engine-measured on the unmoved
fixture: the identical film raises cross-polarized power by
5.2707× at 1 µm and *lowers* it to 0.2035× at
632.8 nm — a reversal of **25.90×**. It is a
measurement, not a unit conversion, because `coat_set` takes *physical*
thickness and the trace divides by the *current* wavelength; the
counterfactual engine that pinned thickness in waves is reachable from the
real one and shows no reversal at all. **Design rule: a protective overcoat
is neither inherently a polarization cost nor a benefit — specify it at λ/4
of the working wavelength.**

**§5 — the contrast floor, and its boundary.** The detector field splits into
co-polarized, cross-polarized and longitudinal channels with the analyzer
*derived* from the pupil coherency matrix rather than assumed to be the
input state, so a purely geometric rotation is not charged to the
cross-polarized channel. On a train where polarization does nothing the
machinery reports exactly zero — which is what makes the floors it reports
afterwards physics rather than its own frame choices. Uncoated two-mirror
floor: 7.0612e-07 of co-polarized power. Coronagraph chain:
**5.787e-13** mean contrast.

> **Read the coronagraph number as a lower bound, not as a built
> coronagraph's floor — this is the package's principal open boundary.**
> Vector mode seeds the component planes at the **first** physical-optics
> leg and thereafter advances them with a common scalar phase, so a
> polarizing surface *after* that leg transforms the rays but not the grid.
> Measured, not assumed: on the two-mirror fixture (both mirrors before the
> only leg) the grid carries 1.000000 of the ray-level
> cross-polarized fraction — full carry. On the 7-mirror
> coronagraph chain it carries **0.8412**, and with all mirrors
> aluminised **0.5653** — at which point the grid does not
> merely understate the coating's effect, it gets the **sign** wrong
> (-0.0322× reported against a ray-level +0.587×),
> because only the first mirror's coating precedes the seed. A run-time
> warning fires. Closing it is Phase 3a Tranche 2, and the tests are written
> so Tranche 2 has to come back and change these numbers.

## Where to go deeper

*Gate index* lists every claim → measured value → truth → CI test; *Coverage
and gaps* lists what is not yet covered and what its absence costs; the
overall plan is `PLAN_POLARIZATION.md`. Numbers marked *(external)* were
measured in another binding, compiler, repository, or against a pre-fix
binary no longer in the tree, and carry their capture date. Review packets,
all in `macos/`:

| Topic | § | Packet |
|---|---|---|
| Conventions; the four engine defects | *Conventions* | `POLARIZATION_PHASE0_AUDIT.md` |
| Jones pupil, aberration maps, vector chain | 1–3 | `REVIEW_POL_2026-07-26.md` |
| Reflected p̂ / Fresnel `r_p` sign | 4 | `REVIEW_POL_SP_SIGN_2026-07-27.md` |
| Contrast floor; Tranche-1 boundary | 5, 5.6 | `REVIEW_POL_2C_2026-07-27.md` |
| Polarizer, waveplate, material-axis rule | 6, 6.7 | `REVIEW_POL_ELEMENTS_2026-07-27.md` |
| Coated-transmission radiometry | 7 | `REVIEW_POL_RADIOMETRIC_2026-07-28.md` |
| External anchor against publication | 8 | `REVIEW_POL_EXTERNAL_2026-07-28.md` |
| Overcoat quarter-wave reversal | 8.3.1 | `REVIEW_POL_OVERCOAT_CHROMATIC_2026-07-28.md` |
