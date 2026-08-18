# Decision packet — coated-Refractor transmission radiometric convention

**Date:** 2026-07-28 · **Lane:** Fable (convention decision) · **Status:**
**LANDED** 2026-07-28 (Opus). Engine `a5e4288`; gates, fixtures and report
section in the same push. Landing record at the bottom.

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

---

# Landing record (Opus, 2026-07-28)

**Engine:** macos `a5e4288` — `elemsub.F`, one block inside the
`ifPol .and. nCoat/=0` branch, immediately after the `ic` loop. `S2` was a
declared-but-unused `REAL*8` in `Refractor`; nothing else in the routine
moved. Guarded `IF (S2.GT.0d0)` so a pathological argument leaves `TP/TS`
untouched rather than producing a NaN.

**The pre-fix A/B, captured before the engine was touched** (rebuilt
gfortran engine + relinked mex; the numbers are in `external.json` as
`X_RAD_PREFIX_*`):

| | uncoated | coated (index-matched) | coated/uncoated |
|---|---|---|---|
| normal | 0.9797958971 = √0.96 | 0.8000000000 = 2n₁/(n₁+n₂) | **0.8164965809 = 1/√1.5** |
| 45°, p | 0.9957577723 | 0.7280089087 | **0.7311104457** |
| 45°, s | 0.9528833281 | 0.6966629547 | **0.7311104457** |
| detector I | 8.9319673126e-01 | 5.9546448751e-01 | **0.6666666667 = 1/1.5** |

All four are 1.0 post-fix. Two things the capture settled beyond the
headline: the uncoated column matches the textbook power transmittance
exactly (√0.96 at normal incidence, and the 45° pair), so the incumbent
convention was confirmed by measurement and not only by line-read; and the
coated column matches the plain composed Fresnel field coefficient exactly,
so the branch was *internally* correct and missing precisely one scalar.

**Gates:** `tPolRadiometric` (mmacos, 13 tests, added to `SUITE_FAST`), on
two new fixtures `Rx_Refract.in` and `Rx_Refract45.in`. Every analytic is
the Abeles characteristic matrix typed from Macleod ch. 2. Non-vacuity:
**6 pass / 7 fail** against the rebuilt pre-fix engine.

**Two corrections to the spec above, both worth recording:**

1. **"Per-layer factors would double-count" is not quite the right
   argument.** In a plain chain the per-interface factors
   `√(nⱼcθⱼ/nⱼ₋₁cθⱼ₋₁)` *telescope* to exactly the boundary-to-boundary
   factor, so a per-interface implementation would give the same number
   here. The decision stands unchanged — one factor, one place — but the
   engine comment now carries the argument that actually holds: the
   telescoping identity is contingent on every interface being present, in
   order, and paired with the *right* media, and this branch has already
   been bitten there once by the stale `nb_arr(0)` boundary slot. A
   per-interface factor built on that slot would telescope to
   `√(n_sub·c_sub/(nb_arr(0)·c₀))` and be silently wrong for a coated
   element following another refractor.
2. **The air-to-air closure gate, as specified, cannot see this defect.**
   For a parallel plate the two faces' factors multiply to
   `√(cosθ_out/cosθ_in) = 1`, so a fully coated plate is air-to-air
   invariant under the landing — verified: it passes against the pre-fix
   engine. The gate was rebuilt as a **mixed** plate (front face coated,
   back face bare), which breaks the cancellation and does discriminate.
   The both-coated case is retained as the composition identity it is
   (decision ground 3, measured rather than argued) and labelled in the
   test and in the report so its green is not mistaken for coverage.

**Analytic gotcha, recorded because it looked exactly like a failed gate:**
Macleod's `t = 2η₀/(η₀B + C)` is the *tangential* amplitude coefficient. For
p-polarization the tangential component is `E·cosθ`, so it exceeds the
ordinary Fresnel `t_p` by `cosθ_sub/cosθ_inc` — 1.2472 at 45° into n=1.5,
i.e. the size of a plausible radiometric error. `T` is unaffected either
way, which is why every transmittance gate passed while the one
amplitude-ratio gate did not.

**Docs:** `polval/80_radiometric.md.in` (section 7) + gate-index rows; the
coverage open item is closed and replaced by the two narrower ones that
remain (absorbing substrate ungated; `Reflector`'s transmittance blocks
still dead code). The "deeper audit question" above is carried into the
coverage section verbatim as an open item that predates this landing.

**Engine `CLAUDE.md`** flipped from open AUDIT FINDING to FIXED with the SHA.

---

## Fable-lane review of the landing (2026-07-28) — PASS, and both spec corrections ACCEPTED

Verified independently on this box: tPolRadiometric re-run 13/13; the
engine block line-read (factor, S2>0 guard — which also handles the
TIR/evanescent edge sensibly by leaving TP/TS untouched — and the
corrected comment); and three numerical checks:

1. **Correction 1 is right and my spec's argument was wrong.**  The
   per-interface factors DO telescope (∏√(nⱼcⱼ/nⱼ₋₁cⱼ₋₁) =
   √(n_sub·c_sub/(n₀c₀)), verified numerically on a 2-layer stack) — a
   correctly-paired per-interface implementation would give the same
   number, so "double-count" was not the hazard.  The hazard that IS
   real, and that the engine comment now carries, is the media-pairing
   contingency: a per-interface factor built on the parser's stale
   `nb_arr(0)` slot telescopes to the WRONG boundary ratio for a coated
   element following another refractor — the exact trap this branch
   already fell into once.  The decision (one factor, one place, from
   `na`) survives with a better justification than it was issued with.
2. **Correction 2 is right and important.**  The both-coated
   parallel-plate closure I specified is exactly invariant
   (entry × exit factors = √(c_out/c_in) = 1, verified) — it passes the
   PRE-FIX engine, i.e. it was a vacuous gate for this defect, and TO
   measured that rather than arguing it.  The mixed-plate rebuild
   (front coated, back bare) breaks the cancellation and
   discriminates; keeping the both-coated case labelled as a
   composition identity is the right disposition.
3. **The pre-fix capture is independently confirmed**: the 45°
   coated/uncoated ratio must equal the inverse radiometric factor
   1/√(n_g·c_g/c_i) = 0.7311104457, polarization-independent — matching
   TO's measured value (both s and p) to all printed digits, and the
   uncoated column's √0.96 confirms the incumbent power convention by
   measurement.

Also endorsed: the merge-not-rebase call (a5e4288 stamped across six
files; a rebase would have dangled every reference — the f10b234 lesson
applied correctly in the opposite direction), and dropping the 2.6 MB of
metadata-only PNG regenerations after hash-verifying pixel identity.

Reviewer's note for the record: both corrections were to MY spec.  This
is the review loop working in both directions — the implementer
gate-checked the reviewer's arguments, found one wrong and one blind,
and recorded why with measurements.  That standard is now the bar.
