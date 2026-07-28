<!-- GENERATED FILE -- do not edit.
     Source: 85_external_anchor.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# External anchor — protected metal, against publication

Every gate before this one compared the engine against **our own** analytic.
That is the right first move, and it caught real defects — but it cannot
detect a shared misunderstanding, and it leaves the reader with no reason to
believe the numbers outside this document. This section closes that gap for
the configuration class the coronagraph results actually depend on: a
**dielectric film on metal**, i.e. a protected mirror.

The gap was specific. `tJonesPupil`'s Fresnel gate covers an optically
**thick single layer** — algebraically a bare interface. `tPolRadiometric`
covers the Abeles matrix in **transmission**. Neither exercises a real
dielectric-on-metal stack, which is exactly what §5.5's coating ladder is
made of.

## 8.1 The method — reproduce someone else's configuration

The anchor is G. van Harten, F. Snik and C. U. Keller, "Polarization properties of real aluminum mirrors I. Influence of the aluminum oxide layer," PASP 121(878), 377-383 (2009) *(external, captured 2026-07-28)*.

It was chosen because it states **both** halves: a full Mueller matrix
measured by ellipsometry at 14 incidence angles over 6–70° and four
wavelengths, **and** every model input numerically — both indices, the film
thickness, and the film model itself.

The engine is then driven with **their** indices, **their** film thickness,
**their** wavelengths and **their** incidence angles, and compared
curve-on-curve against **their** model. That construction is deliberate.
Aluminium index tables genuinely disagree with one another — the publication
says so itself, and *fits* the imaginary index rather than adopting a table.
A disagreement traceable to index tables is **not** a machinery error, and
this comparison cannot confuse the two.

`macos.coating` takes arbitrary n, κ and physical thickness, so reproducing a
published stack needs no engine change and no fixture surgery.

Two conventions had to be settled rather than assumed:

- **Index sign.** The publication uses `N = n − ik`, k ≥ 0 = loss — the same
  convention the engine stores. This is load-bearing, not cosmetic: the paper
  reports that the opposite sign made their fitted oxide come out at ~50 nm
  instead of ~4 nm.
- **p̂ frame.** Measured, not derived, on the perfect-conductor case, whose
  answer is known exactly: the engine's r_s/r_p ratio comes back at **+1**
  with zero imaginary part at every angle, so the bridge to the publication's
  frame is **zero**. An earlier version of this harness *assumed* a π bridge
  from the ray-following doctrine and was wrong by exactly 180° everywhere.

The analytic is written from Macleod ch. 2 and the publication's stated
equations, **never transcribed from `elemsub.F`** — an "analytic" copied out
of the engine is circular in exactly the coefficient it is meant to check,
which is how the r_p sign defect of §4 survived four years of gates.

## 8.2 Result — the machinery is anchored

Over the publication's four wavelengths and six incidence angles inside
their measured range, comparing **per ray** at each ray's own incidence
cosine:

| quantity | worst deviation | publication's own accuracy |
|---|---|---|
| diattenuation (the [1,2] Mueller element) | 2.828e-14 | 0.01 |
| retardance, in the same Mueller units | 4.937e-14 | 0.01 |

The engine sits inside the published measurement bar by a factor of
**2.0e11**. Both sides implement the same thin-film theory, so
round-off agreement is the expected result; the publication's ±0.01
is what the comparison would have had to beat to mean anything, not the
achieved figure.

**Non-vacuity.** The comparison can fail, and the two ways it can are the
publication's own physics:

- Omitting the 4.12 nm oxide — comparing the engine's oxidized mirror against
  a *bare* analytic — disagrees by **exceeds 0.01**, more than their
  measurement accuracy. That independently reproduces the paper's central
  claim, that a ~4 nm oxide must be modelled at all.
- The historical ~50 nm oxide, which the paper attributes to a sign error on
  the imaginary index, is separable from 4.12 nm by **exceeds 0.01**.

## 8.3 What this changes about §5.5 — the overcoat trade reverses

With the machinery anchored, §5.5's coating ladder becomes checkable
*internally*: the per-mirror coefficients can now be trusted, so the
cross-polarized power they imply can be predicted and compared with what
`macos.pol_contrast_floor` reported.

For an on-axis rotationally symmetric train the mirrors' s/p axes share one
azimuth, so their Jones matrices commute into
`J = c₀I + c₂·[cos2φ, sin2φ; sin2φ, −cos2φ]`, and the cross-polarized
amplitude is set entirely by `ε_tot = c₂/c₀`. The cross **power** ratio
between two coatings is therefore |ε_tot|² between them.

| source | MgF₂/Al ÷ bare Al, cross power |
|---|---|
| independent analytic at the fixture's own wavelength | 6.35 |
| the engine's own Jones pupil | 6.4225 |
| `macos.pol_contrast_floor`, as reported in §5.5 | 5.4238 |

So **the §5.5 number is arithmetically correct.** The three agree; the
residual spread is dark-zone annulus weighting against a pupil mean.

What is *not* correct is the description attached to it. `Rx_Cass_FarField`
runs at `Wavelen = 1.0E-06` m, while the coating constants
are labelled "Al at 632.8 nm" and the MgF₂ thickness "a quarter wave". At the
fixture's actual wavelength, 110 nm of MgF₂ is **0.607**
quarter-waves — not a quarter wave, and not an ordinary protected-aluminium
recipe.

That distinction decides the sign of the design rule, because the overcoat
trade reverses across the quarter-wave condition. The same 110 nm film, the
same aluminium, the same few-degree incidence:

| film | optical thickness | cross power vs bare Al |
|---|---|---|
| 110 nm MgF₂ at 632.8 nm | 0.96 quarter-waves | 0.18 |
| 110 nm MgF₂ at 1 µm *(the fixture)* | 0.607 quarter-waves | 6.35 |
| **true quarter wave** (181.2 nm at 1 µm) | 1.00 quarter-waves | 0.0157 |

A genuine quarter-wave overcoat **suppresses** the polarization floor by
about 1.8 decades. §5.5's sentence — that the protective overcoat "costs most
of another" decade — is therefore withdrawn: it is true of the film as
configured, and false of the recipe that film was described as.

The conclusion does not rest on pinning the fixture's exact incidence angles:
the ratio is flat at **6.34 .. 6.36** across the whole plausible
near-normal range. (A per-element AOI was deliberately *not* quoted — the
two routes available through this binding are both unusable, and the tool
records why.) The crossover where a quarter-wave overcoat stops helping and
starts hurting sits near 35°, so the conclusion is specific to near-normal
trains and must be re-derived for a fold-heavy layout.

**No engine change follows from any of this.** The thin-film machinery
reproduced a published protected-metal measurement to 2.828e-14; what
moved is a fixture's configuration and the prose describing it.

## 8.4 What this section does not claim

- It compares against the publication's **model**, evaluated at the
  publication's stated inputs. Their measured Mueller matrices are published
  as figures, not tables, and digitising a figure would add a transcription
  error larger than anything being measured. The publication's own
  model-to-measurement agreement (±0.01 per element) is the outer
  envelope, and it is quoted as such.
- One wavelength region, one metal, one dielectric. It anchors the
  characteristic-matrix recursion on dielectric-on-metal; it says nothing
  about absorbing overcoats, many-layer stacks, or the transmission side
  (§7 covers that separately).
- The aluminium indices used in §5.5 remain a **model choice**, not an
  anchored one. This section pins the machinery, not the material data.
