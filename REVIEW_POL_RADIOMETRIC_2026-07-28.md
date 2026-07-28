# Decision packet — coated-Refractor transmission radiometric convention

**Date:** 2026-07-28 · **Lane:** Fable (convention decision) · **Status:**
DECIDED — implementation + gates are an Opus slice (spec below).

## The finding (session A audit, 2026-07-27)

Coated and uncoated `Refractor` transmission use different amplitude
normalizations. Measured: with an index-matched single layer (optically a
bare interface) the coated/uncoated |Ex| ratio is **0.816442 = 1/√1.5
exactly** at normal incidence — a coated lens under-transmits ~18% in
amplitude relative to the identical uncoated surface.

## Line-read (this review)

* **Uncoated** (`elemsub.F:1146-1149`): `S1 = 2·√(cosθₜ/(cosθᵢ·μ))`,
  `μ = n₁/n₂`, folded into `TP/TS`. Expanding: `TP = t_p^Fresnel ·
  √(n₂cosθₜ/(n₁cosθᵢ))` — the exact radiometric factor, so
  `|TP|² = T_p`, the POWER transmittance. **The incumbent engine
  convention is: the transmitted amplitude carries the ray's power.**
* **Coated** (`elemsub.F:~1200-1262`): innermost seed + per-layer `TP1`
  are plain Fresnel field coefficients composed through the Airy
  recursion; no radiometric factor anywhere.

## The decision

**Keep the incumbent power-amplitude convention; bring the coated branch
to it.** Grounds:

1. The internal inconsistency is the bug, independent of which convention
   is "right" in the abstract: an index-matched layer is physically a
   bare interface and MUST reproduce the uncoated result. Today it is off
   by 1/√n.
2. The uncoated convention is what every existing polarized result and
   gate sat on; changing IT would churn gated behaviour to fix a branch
   that has no gates yet.
3. The power-amplitude convention is self-consistent air-to-air: for a
   full transit chain the per-face factors compose to the total power
   transmittance (the n's and cosines telescope), which is what the
   detector-side bookkeeping wants.

**The fix is ONE factor, applied ONCE, after the recursion:**

```
TP = TP * sqrt( REAL(n_sub)*cos_sub / (na*cos_inc) )     (TS likewise)
```

with `n_sub = nb_arr(nCoat+1)`, `cos_sub = DBLE(ccfb_arr(nCoat+1))`,
`na` the Phase-2-corrected incident medium, `cos_inc = ccfb_arr(0)`.

Do **NOT** add per-interface factors inside the recursion: the standard
multilayer result is `T = (n_sub·cosθ_sub)/(n_inc·cosθ_inc)·|t_composed|²`
— the interior-layer factors cancel identically, so per-layer factors
would double-count. That cancellation is also why the fix cannot change
any REFLECTED quantity (`RP/RS` untouched) and why the r·r Airy
denominators are unaffected.

**Real parts, and a scope line:** the factor is power bookkeeping, so it
uses the real parts of the indices. A metal-substrate `Refractor`
transmission (evanescent) has no meaningful power convention — out of
scope, document at the site. Coated LENSES (dielectric substrate) are the
use case.

## Landing spec (one Opus slice)

1. The one-line factor (×2 for TS), with a comment carrying the
   interior-cancellation argument and the metal-substrate scope line.
2. Gates, textbook-analytic and non-vacuous:
   * index-matched layer ≡ bare interface: coated/uncoated amplitude
     ratio = 1 (the current 1/√1.5 is the built-in A/B);
   * MgF₂-on-glass single layer at normal AND 45°, s and p, vs the
     textbook multilayer `T` (written fresh, not from the engine);
   * air-to-air closure on a coated parallel plate: total transmitted
     power vs analytic;
   * pol-off bit-identity (the factor lives inside `ifPol`);
   * a λ-sweep point pinning that the quarter-wave interference
     structure (which the recursion already had) is unchanged by the
     scalar factor.
3. polval: a short evidence subsection + the coverage-section open item
   closed; the engine `CLAUDE.md` AUDIT-FINDING paragraph flips to
   FIXED with the convention stated.
4. Suites: mmacos full + ifx build + GMI (ifPol-off ⇒ must be
   bit-identical) + pymacos regression.

## Recorded, no action — the deeper audit question

The cosθ part of the radiometric factor overlaps conceptually with
beam-area bookkeeping that ray DENSITY also captures when the grid is
seeded for diffraction. Whether that constitutes a double-count in
oblique refractive seeding is an open audit question that predates this
decision (the uncoated branch has always carried it, and the PROPER
comparisons — mostly mirror trains — never exercised it hard). Left
open deliberately; do not couple it to this fix.
