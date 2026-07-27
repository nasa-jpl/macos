<!-- GENERATED FILE -- do not edit.
     Source: 50_sp_sign.md.in
     Numbers: generated/numbers.json (MATLAB driver)
     Regenerate: make polval-regen  (docs/macos-manual)
-->
# The reflected-p̂ / Fresnel r_p sign — an engine convention corrected

Every other measurement in this report is taken on a train with an **even**
number of mirrors. That is not a coincidence of fixture choice, and until
2026-07-27 it was load-bearing: the engine's `Reflector` assembled the
reflected field on a p̂ that follows the *outgoing* ray but multiplied the
p-component by **−r_p**, which reflects the transverse field about the local
p̂ instead of negating it. That operation is an **involution** — it cancels
*exactly* across a mirror pair — and it is **unitary**. The stock Cassegrain
fixture has exactly two mirrors, the vector-chain fixture has none, and a
unitarity gate cannot see a unitary error. The defect arrived in the 2022
bulk source import with no recorded rationale; the standard form sat
commented out one line above, annotated as the original.

The correction restores the Born & Wolf `+r_p` in **both** branches: the
uncoated `RP`, and the coated recursion's innermost `RP` plus its per-layer
`RP1`. Substituting `−r` at every interface into the Airy recursion negates
the result while the denominator's double sign cancels, so the coated form
was the same clean sign restoration, and the magnitude-and-ratio fold gate of
§2.2 could never have seen it.

**All two-mirror results in this report are unaffected**, which is itself a
check on the claim: measured on the same fixture and build pair, the
two-mirror cross-polarized power fraction is 7.0612e-07 before
the fix and 7.0612e-07 after it — bit-identical *(external,
captured 2026-07-27; the reproducer runs at model 256, while the
numbers measured below are this driver's at model 128)*. Regenerating this
whole report across the fix moved no physics number: of 67 numeric tokens,
33 are bit-identical and the rest move only within round-off — the largest
is the unitarity gate's peak retardance, which *improves* from 2.46e-15 to
5.79e-16.

## 4.1 One mirror against the perfect-conductor closed form

The stock mirrors carry the perfect-conductor idiom (`IndRef=1`,
`Extinc=1e22`), giving `r_s = −1`, `r_p = +1` in the engine's ray-following
basis. A single reflection of an x-polarized collimated beam is then fixed by
geometry alone. With AOI `a` and pupil azimuth `φ`, exactly — no small-angle
expansion:

```
Ey/Ex = -sin(2φ)·sin²(a) / den        Ez/Ex = -sin(2a)·cos(φ) / den
den   = 1 - 2·sin²(a)·cos²(φ)
```

Both `a` and `φ` are taken from the engine's **ray directions** — `a` from
the stop-to-mirror deflection, `φ` from the outgoing transverse direction —
so the comparison assumes nothing about how the pixel grid maps to the pupil.
The expression is written from the textbook, *not* transcribed from the
engine's own expression: that circularity is exactly why the fold gate of
§2.2 passed both before and after the fix, and it has since been removed
there too.

Over 10188 rays spanning up to 10.49° of incidence, the
engine reproduces the closed form at **median 2.55e-15, max
6.39e-14** relative for the transverse component and
5.03e-14 for the longitudinal one — round-off. The perfect
conductor also introduces no retardance: 1.07e-16.

## 4.2 The half that needs no reference value at all

Cross-polarization from an isotropic, rotationally symmetric mirror is
driven by the local surface slope. It must therefore vanish on axis, grow as
ρ², and stay bounded by O(sin² AOI). This is fixture-independent: it can be
checked without any reference number, and it is the tell that identified the
defect in the first place.

| ρ / ρ_max | 0.25 | 0.50 | 0.75 | 1.00 |
|---|---|---|---|---|
| median \|E_y/E_x\| | 1.812e-03 | 6.262e-03 | 1.407e-02 | 2.408e-02 |

The log-log slope is **1.871** against the required 2; the peak
cross-polarization is 1.034 of the O(sin² AOI) bound, and the
cross-polarized power fraction after one mirror is 2.126e-04.

![One mirror: measured cross-polarization, the perfect-conductor closed form,
and the radial law.](media/polval_spsign.png)

## Non-vacuity, measured

The engine was rebuilt with the sign flipped back and the same gates re-run
*(external, captured 2026-07-27 — the pre-fix binary is not in the
tree)*. Nothing about these gates is satisfiable by a broken engine:

| Quantity | Pre-fix | Post-fix | Gate |
|---|---|---|---|
| closed-form residual (median) | 1.14e+02 | 2.55e-15 | < 1e-11 |
| radial profile, ρ = 0.25 … 1.00 | 0.988 / 1.038 / 1.029 / 1.035 | 1.812e-03 … 2.408e-02 | grows as ρ² |
| radial log-log slope | 0.033 | 1.871 | > 1.7 |
| retardance | 3.9e-10 | 1.07e-16 | < 1e-14 |

and, from the standalone reproducer at model 256 *(external, both cells)*:

| Quantity | Pre-fix | Post-fix |
|---|---|---|
| cross-pol power, ONE mirror | 1.0163e+00 | 2.0724e-04 |
| cross-pol power, TWO mirrors | 7.0612e-07 | 7.0612e-07 |

A single near-normal mirror turning x-polarized light into a roughly 50/50
mixture is not available to physics: the cross-polarization an isotropic
surface can generate is O(sin² β) in amplitude, β the local slope. The
pre-fix profile was also **flat in pupil radius**, which no slope-driven
effect can be.

## Scope of the correction

`Refractor` carries its own copy of the coated recursion and had the same
negated `RP`/`RP1`. There the flip was **internal**, not emitted: the element
assembles its output from the transmission coefficients (p̂ does not flip on
transmission, and `TP1`/`TP` were already standard), and `RP` reaches the
transmitted field only through the products `RP1·RP` in the Airy
denominators, which are invariant under a *consistent* flip. Both were
normalized anyway, so that a future beamsplitter or polarizing-beamsplitter
element does not inherit a flipped reflected branch; the transmitted field
was verified **bit-identical** across the change on a coated singlet, and the
inconsistent flip that a careless edit would produce was built and measured
to move transmitted power by −3.2%, so the invariance is a real result rather
than an untested path.

Two limitations of this section are worth stating plainly. It validates the
**perfect-conductor** single reflection; a coated single mirror is covered
only through the ratio-form Fresnel gate of §2.2, which is on an even-mirror
rig. And the coated `Refractor` branch has **no analytic gate at all** —
neither before this work nor after — while carrying a normalization
discrepancy against its own uncoated branch (see *Coverage and gaps*).
