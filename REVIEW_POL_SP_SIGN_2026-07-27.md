# Review packet — reflected-p̂ / Fresnel-r_p sign conflict in `Reflector`

**Date:** 2026-07-27 · **Branch:** `pol-core` (both repos) · **Lane:** Fable
(this is a convention decision in the engine's polarization core, not a
template-follow)

**Status:** diagnosed and verified by experiment; **nothing landed**.  The
engine tree is unmodified — the candidate fix was built, measured, and
reverted (`git checkout elemsub.F`, both compilers back to green,
`run_mmacos_tests.sh fast` 287 pass / 0 fail).

**Why you are being asked:** Phase 2c (contrast floor, Opus worklist item 3)
is blocked on this, and the fix flips a sign in the Fresnel core that every
polarized result rides on.  One decision needed — see *The decision* at the
end.

---

## The claim

`Reflector` (`macos_f90/elemsub.F`) assembles the reflected field in a basis
whose p̂ **follows the outgoing ray**, but multiplies the p-component by a
Fresnel coefficient written in the convention where p̂ does **not** flip.
The two are off by one sign relative to each other.

The three lines, verbatim:

```fortran
430:  CALL DXPROD(pihat,shat,ihat) ! incident P pol vector
431:  CALL DXPROD(prhat,shat,rhat) ! reflected P pol vector      <-- p̂ ∥ OUTGOING ray

454:  !RP=(NbCmplx*ccfa-NaCmplx*ccfb)/(NbCmplx*ccfa+NaCmplx*ccfb) ! dcr's original
455:  RP=(NaCmplx*ccfb-NbCmplx*ccfa)/(NbCmplx*ccfa+NaCmplx*ccfb)  ! fixed by jzlou
456:  RS=(NaCmplx*ccfa-NbCmplx*ccfb)/(NaCmplx*ccfa+NbCmplx*ccfb)

592:  Eout(1)=DCMPLX(prhat(1))*Epr+DCMPLX(shat(1))*Esr            <-- uses prhat
```

Line 455 is `−r_p` in the standard (Born & Wolf / "p̂ follows the ray")
convention; line 454, commented out, is `+r_p`.  At normal incidence line 455
gives `RP = RS`, which is the signature of the *other* convention — the one
where p̂ is held fixed through the surface.  Line 431 says p̂ is not held
fixed.

**Consequence.** Write the transverse field as `E = E_s ŝ + E_p p̂`.  With
p̂_r ≈ −p̂_i near normal incidence and `RP = RS = −1` (the perfect-conductor
idiom), the assembled output is

```
E_out = −E_s ŝ − E_p p̂_r  =  −E_s ŝ + E_p p̂_i
```

which is a **reflection of the transverse field about the local p̂ axis**,
where the physical answer for an isotropic surface is `E_out = −E_in`.  For a
rotationally symmetric mirror p̂ is radial, so the reflection maps a uniform
linear input into a polarization pattern that winds with azimuth — at full
amplitude, with no small parameter in front of it.

---

## The evidence

Reproducer, one command, in the bindings repo:
`MACOS_resources/mmacos/tools/pol_sp_sign_probe/probe_sp_sign.m`.
It runs on **`Rx_Cass_FarField.in`** — the same fixture the Phase-2 gates use
— and reads `RayE` directly (`macos.ray_field`), so nothing in the
diffraction layer is involved.

### 1. One mirror shows it; two cancel it

x-polarized input, perfect-conductor mirrors, model 256.  `elt 2` is after the
Primary (ONE mirror), `elt 3` after the Secondary (TWO):

| | current engine | with `+r_p` restored | change |
|---|---|---|---|
| after 1 mirror, `Py/Px` | **1.0163e+00** | **2.0724e-04** | 4900× |
| after 2 mirrors, `Py/Px` | 7.0612e-07 | 7.0612e-07 | **bit-identical** |

A single mirror at ≤2° AOI turning x-polarized light into a 50/50 x/y mixture
is not available to physics.  The cross-polarization an isotropic surface can
generate is `O(sin²β)` in amplitude, β = local surface slope — ~1e-3 here.
`2.07e-04` in *power* (≈1.4e-2 in amplitude) is in that range; `1.02` is not.

The even-mirror column is the whole story of why this survived: the defect is
an involution (a reflection satisfies `M² = I`), so it cancels **exactly**
across a mirror pair.

### 2. It does not scale with radius — so it is not AOI-driven

After ONE mirror, median `|Ey/Ex|` by pupil radius:

| ρ (px) | 31.9 | 63.8 | 95.7 | 127.6 |
|---|---|---|---|---|
| median `\|Ey/Ex\|` | 1.0142 | 1.0160 | 1.0098 | 1.0051 |

Flat.  Any real polarization effect of an isotropic rotationally symmetric
mirror must vanish on axis and grow as ρ² (it is driven by the local slope).
This is fixture-independent and needs no comparison to a reference value.

### 3. It is not the point-source launch frame

The engine launches point-source rays in a per-ray frame
(`yray = unit(RayDir × xGrid)`, `ssrcray.inc`), which is documented as a
reason per-component maps "aren't globally aligned".  Control experiment:
`Rx_VecChain.in` with `zSource` flipped to finite, reading `ray_field` at
element 1 — an `Obscuring` stop, so the value **is** the launch state,
untouched by any surface:

```
collimated : Py/Px = 0            max|Ey/Ex| = 0
point src  : Py/Px = 1.78e-11     max|Ey/Ex| = 1.02e-05   at NA 0.0045
```

`1.02e-05 ≈ NA²/2` — the launch frame contributes at second order in NA and
is three to five orders below the effect above.  Ruled out.

### 4. Which fix — and why it is not a free choice

Both of these repair normal incidence:

* **A. restore `+r_p`** (line 454), keep `prhat`;
* **B. keep `−r_p`**, assemble with `pihat` instead of `prhat`.

**B is not admissible.**  `pihat = ŝ × î` is perpendicular to the *incident*
direction, not the reflected one, so B emits a field with a component along
the outgoing ray — an unphysical longitudinal term that would then propagate.
**A is the fix.**

### 5. What the candidate fix does to the suite

`run_mmacos_tests.sh fast` against the patched engine: **285 pass, 2 fail.**

Everything carrying physics is unchanged and green — unitarity, the
Fresnel-analytic fold, 2θ symmetry, `pol_zernike` two-mirror form, all 8
`tVecChain`, all 10 `tPolarization`.

The two that move are both **basis-artifact** assertions in `tJonesPupil`:

| test | assertion | current | patched |
|---|---|---|---|
| `test_basis_invariance_and_sp_artifact` | `var_rms.ret` local-sp > 10× double-pole | passes | local-sp 3.346e-3 vs double-pole 3.607e-3 — the two bases now **agree to 8%** |
| `test_pol_zernike_basis_dependence` | local-sp must inflate the retardance expansion | passes | **both bases give exactly 0** |

Read that carefully, because it cuts in favour of the fix rather than against
it.  The documented "local-sp inflates retardance variation 10–250×" was, on
this fixture, substantially *this defect* rather than a pure coordinate-basis
effect — and with `+r_p` the two-mirror perfect-conductor system reports
**zero retardance in both bases**, which is the physically correct answer for
a perfect conductor.  Those two assertions would need to be rewritten (or
re-pointed at a geometry where local-sp really is singular), not preserved.

Note the Fresnel-analytic fold gate **passes unpatched and patched** because
the fold rig is *coated*, so it runs the `nCoat/=0` branch, which the scratch
patch did not touch.  That branch carries the same `−r_p` form (innermost
`RP`, and `RP1` per layer) and would need the same treatment.

---

## Why the current gates are structurally blind to this

Not bad luck — four independent reasons, all of which happen to hold:

1. **`Rx_Cass_FarField` has exactly two mirrors.**  Exact cancellation.
2. **`Rx_VecChain` has no reflectors.**  Nothing to get wrong.
3. **The unitarity gate cannot see a unitary error**, and a reflection about
   p̂ is unitary.
4. **The Fresnel gate is circular in this sign.**  `tJonesPupil.m:166` builds
   its "analytic" `RPa` with the *engine's own* expression, then compares the
   ratio `RS/RP`.  A sign shared by both sides cancels out of the comparison.

Adding a gate is therefore part of the fix: an **odd-mirror-count** fixture,
asserting cross-polarization scales as ρ² and stays below `O(sin²β)`.  I did
not add it now — a test that fails against the shipped engine cannot land.

---

## Scope still to be audited (I did not)

* the **coated** `Reflector` branch (`RP` innermost + `RP1` per layer) — same
  `−r_p` form;
* **`Refractor`** (`elemsub.F:1102-1103`, `1226-1228`) — same `pihat`/`prhat`
  pattern with its own `RP`.  Note transmission is a genuinely different
  question: p̂ does *not* flip on transmission, so `TP`/`TS` (~1569-1591) may
  legitimately need the other sign;
* whether anything downstream **depends** on the current behaviour — the two
  `tJonesPupil` assertions above are the only ones the fast suite surfaces,
  but masks/proper/freeform were not run against the patch.

## Provenance

The flip is not a recoverable decision: `git log -S` puts it in the 2022 bulk
import `e1fa721` ("Commiting all MACOS source files delivered by John Lou"),
already in its current form.  The `! dcr's original` comment line was added
*later* (`50ff00e4`, 2023-11-08) as a note. So there is no commit message,
issue, or test explaining what the flip was meant to fix — worth weighing
before assuming it was load-bearing somewhere unmeasured.

---

## The decision

**Restore the standard `+r_p` (line 454) across `Reflector`'s uncoated and
coated branches, and rewrite the two basis-artifact assertions — yes or no?**

If yes, the follow-on work is mechanical and Opus can take it: apply to both
branches, audit `Refractor`, add the odd-mirror ρ² gate, re-run both
compilers + GMI + PROPER, regenerate `polval` (its unitarity/2θ/Zernike
numbers are all two-mirror, so most should be unchanged — that is itself a
check), and update the conventions table in `macos_f90/CLAUDE.md`.

If no, Phase 2c needs a different justification for trusting the pupil
polarization state on odd-mirror trains, because today it is scrambled at
O(1).

---

## Side findings from the same session (recorded, not decisions)

* **`Rx_Coro.in` must be run at model ≥ 512.**  It declares `nGridpts=511`;
  below that the engine prints `Too many grid points. Resetting npts to N`
  and then intermittently SIGSEGVs in `intensity` (~30–50% of runs at 128,
  deterministically under `trace`).  Dave: grid size must be ≤ model size, or
  use `MREset`.  This **corrects** the note in `PLAN_POLARIZATION` §2c item 3
  saying `Rx_Coro` runs at 128 — it appears to, then doesn't.  At 512 and
  1024 it is fully deterministic (repeat runs bit-identical).
* **`macos.trace(e)` on `Rx_Coro` loses every ray past element 7** (`nOK = 0`,
  RMS OPD at the 9.9999e36 sentinel) while `intensity`/`complex_field`
  propagate fine.  Unrelated to polarization; not chased.
* Two results worth keeping for whenever 2c unblocks, both measured on
  `Rx_Coro` at model 512 and both making the machinery simpler than planned:
  * the chain is **linear in the input Jones state** to 4.2e-16 (45°) and
    5.8e-16 (circular) — `E(a·x̂ + b·ŷ) = a·E(x̂) + b·E(ŷ)` at the detector;
  * therefore a spatially uniform analyzer **commutes with propagation**, so
    the co/cross split can be taken at the detector on the engine's own
    component planes.  No pupil multiplier is needed at all — which retires
    both the Jones-pupil-multiplier design *and* the "divide by the scalar
    run" successor sketched in the plan.  Parseval on that split is exact
    (3.0e-16 … 4.5e-16).
  * incidental: `complex_field(..., 'reset_trace', false)` returns
    bit-identical planes ~100× faster (0.01 s vs 0.83 s at model 512), so
    reading three component planes costs one propagation, not three.

---

## Fable-lane decision (2026-07-27, appended at landing)

**YES — the standard `+r_p` is restored in BOTH `Reflector` branches.
Landed by the Fable lane** (this is the convention-arbitration case the
lane split exists for), with the mechanical tail left for Opus below.

**Independent verification beyond the packet:**

1. **Line-verified the algebra.**  Line 456's `RS` is textbook r_s; the
   commented 454 is textbook Born&Wolf r_p (ray-following p̂); active 455
   is its exact negation.  With `prhat = ŝ×r̂` (431) the assembly (592)
   requires +r_p: PEC at normal incidence then gives `E_out = −E_in`.
2. **The coated branch's flip propagates EXACTLY.**  Substituting −r at
   every interface into the Airy recursion `(r₁+r·C)/(1+r₁·r·C)` negates
   the output while the denominator's double sign cancels — final coated
   `RP` was exactly −(standard multilayer r_p), |RP| unchanged.  So the
   coated fix is the same clean sign restoration (innermost `RP` +
   per-layer `RP1`), and the magnitude/ratio fold gate could never see it.
3. **Internal-inconsistency evidence:** the same file's (dead)
   transmittance code (`TP1`, ~568) is ALREADY in the standard
   convention — the file disagreed with itself; "dcr's original" is the
   self-consistent version.
4. **Probe reproduced bit-for-bit** pre-fix (1.0163 / 7.0612e-07) and
   post-fix (2.0724e-04 / 7.0612e-07 bit-identical).  Post-fix bonus the
   scratch run couldn't show: the radial profile now GROWS as ρ²
   (1.63e-3 → 6.2e-3 → 1.38e-2 → 2.39e-2 at ρ = 32/64/96/128 px; ratios
   3.8/8.5/14.7 vs ρ² = 4/9/16) — slope-driven AOI physics restored.

**One packet conclusion is CORRECTED — the suite prediction §5.**  The
scratch patch touched only the UNCOATED branch, so the Al-coated
secondary stayed flipped: the "two bases now agree to 8%" numbers
describe a HALF-BROKEN system (one fixed mirror + one flipped), not the
fix.  With BOTH branches fixed, the local-sp artifact is real and large
on the correct engine: var_rms ret = 0.8913 (sp) vs 3.6067e-3 (dp) —
ratio **247×**, at the top of the documented 10–250× range, with the sp
mean retardance parking at π/2 (the coordinate artifact in person).
**The two basis-artifact assertions are therefore KEPT UNCHANGED** and
pass non-vacuously; the only test change needed is the fold gate's
analytic `RPa`/`RP2`, transcribed to the textbook form so the phase
comparison stops being circular in this sign (packet §"blind" item 4).

**Refractor is deliberately untouched** (confirming the packet's scope
call): its transmitted output contains only r·r products, which are
invariant under the consistent flip, and its `TP1` is already standard.
Its internal `RP` convention should still be normalized in the audit for
consistency.

**Left for Opus (one session-slice each, per the discipline rules):**
1. Odd-mirror ρ² gate (the probe's radial table is the template; a
   fixture with 1 or 3 mirrors, assert cross-pol grows as ρ² and stays
   under O(sin²β)).
2. Refractor audit (normalize its internal RP/RP1 to standard; verify
   transmitted output bit-identical, which the algebra above predicts).
3. polval regen + re-run PROPER/GMI/masks/freeform (two-mirror numbers
   should be unchanged — that is itself a check; GMI is ifPol-off and
   must be bit-identical).
4. §2c unblocks: proceed with the analyzer-at-detector design (endorsed
   — the linearity measurement at 4e-16 justifies commuting a uniform
   analyzer past the chain, and it retires both pupil-multiplier
   sketches; amend §2c's text accordingly as part of that landing).

---

## Closeout — items 1–3 landed (2026-07-27, Opus lane)

Commits: mmacos `fc2e22e` (odd-mirror gates), macos `25c4386`
(Refractor normalization), plus the polval section and doc updates below.
Item 4 (§2c) is explicitly NOT in this session.

### 1. Odd-mirror gate — stronger than the spec asked for

The spec was "assert cross-pol grows as ρ² and stays under O(sin²β)".
Both are asserted, but the fixture turned out to support something much
sharper: on the perfect-conductor idiom (`IndRef=1`, `Extinc=1e22`) the
ENTIRE single-reflection Jones is fixed by geometry, so the engine can be
compared to a closed form rather than to a bound.  With `a` = AOI and
`φ` = pupil azimuth, exactly (no small-angle expansion):

```
Ey/Ex = -sin(2φ) sin²(a) / den      Ez/Ex = -sin(2a) cos(φ) / den
den   = 1 - 2 sin²(a) cos²(φ)
```

derived from `E_out = (E·p̂_i) p̂_r - (E·ŝ) ŝ` with `r_s = -1, r_p = +1`.
**Measured agreement: median 2.1e-15, max 5.9e-14 relative** over 10188
rays spanning AOI up to 10.5° (longitudinal component 2.7e-15 / 5.0e-14).
Retardance from a perfect conductor: 1.1e-16.

Two properties make this non-circular in exactly the way §"blind" item 4
demands: the expression is written from Born & Wolf, not transcribed from
the engine; and BOTH `a` and `φ` come from the engine's RAY DIRECTIONS
(`a` from the stop→mirror deflection, `φ` from the outgoing transverse
direction), so no pixel-grid-to-pupil mapping is assumed — a transposed
or mirrored grid convention cannot make it pass.

The fixture-free half is asserted too: max cross-pol = 1.03× the
O(sin²AOI) bound, cross-pol power fraction 2.09e-4, radial medians
1.81e-3 / 6.26e-3 / 1.41e-2 / 2.41e-2 at ρ = 0.25/0.50/0.75/1.00,
log-log slope **1.871**.  No new fixture was needed: `Rx_Cass_FarField`
element 2 IS the one-mirror state (RayE/RayDir are the current trace
state, not a per-element history, so `trace(2)` then `ray_field(2)`
gives it — the two-mirror value only appears if you trace to 3 first).

**Non-vacuity, measured rather than argued.**  The engine was rebuilt
with the uncoated `RP` flipped back, at model 128 on the same fixture,
and both tests fail on 7 of their 8 assertions: closed-form residual
median **1.14e+02** / max 1.61e+05 (asserted < 1e-11); radial profile
flat at 0.988/1.038/1.029/1.035 with slope **0.033** (asserted > 1.7);
cross-pol power fraction 3384 (asserted < 1e-3); retardance 3.9e-10
(asserted < 1e-14).  Tree restored and re-verified afterwards: the probe
reproduces 2.0724e-04 / 7.0612e-07 bit-for-bit.

### 2. Refractor — normalized, and the algebra's prediction confirmed

Innermost `RP` and per-layer `RP1` restored to Born & Wolf; the stale
`!RP1=-RP1 ! testing` line dropped.  The packet's r·r argument holds
exactly: the transmitted Ex/Ey/Ez are **BIT-IDENTICAL** before and after,
on a Bench-emitted singlet with a 3-layer absorbing stack on the powered
face, 209 rays, for BOTH 45° linear and circular input.

That A/B could have been vacuous (if the coated branch never ran, or if
`RP` never reached the output), so the inconsistent flip was built too —
`RP1` flipped while `RP` stayed standard, which is what a careless edit
produces.  It moves transmitted x-power by **−3.24%**, max per-ray
|ΔEx|/|Ex| = 1.87e-2.  So the path is live and the invariance is a
result, not an absence.

**New finding, recorded not fixed.**  Coated and uncoated `Refractor`
transmission are on DIFFERENT amplitude normalizations: the uncoated
branch applies the radiometric factor `sqrt(n₂cosθ₂/(n₁cosθ₁))` (the
`S1` at elemsub.F ~:1147), the coated branch omits it entirely.  Measured
with an index-matched single layer (n=1.5 on an n=1.5 substrate —
optically a bare interface): coated/uncoated |Ex| = **0.816442** at
normal incidence, exactly 1/√1.5, falling to 0.804789 off-axis.  A coated
lens under-transmits by ~18% in amplitude (~33% in power) relative to the
same surface uncoated.  Fixing it changes results and deserves its own
decision + gate, so it is logged in `PLAN_POLARIZATION`, the engine
`CLAUDE.md`, and the report's coverage section instead.  Related gap:
the coated Refractor branch has **no analytic gate at all** — the
Fresnel gate covers `Reflector` only.

### 3. Re-runs and the report

| Gate | Result |
|---|---|
| mmacos full suite (gfortran) | see the run recorded in the commit below |
| GMI regression | **6/6, `vs-ref = 0.000e+00`** — bit-identical, as required for an ifPol-off consumer |
| polval regen | all gate thresholds pass, including 5 new ones |

The two-mirror numbers ARE unchanged, which was the point of running
them: every published 2a/2b/polval result stands.

Report: new section `polval/50_sp_sign.md.in` (§4) with the closed-form
comparison, the radial law, the pre-fix A/B table, and the Refractor
scope note; gate-index rows 4.1/4.2; engine-defect entry #4 in the
conventions section, including the note that the §2.2 fold gate could
not have caught this because its analytic was transcribed from the
engine (now textbook).  The pre-fix numbers went to `external.json` as
historical, with the command that reproduces them.

### Corrections to earlier text in this packet

* §5's suite prediction was already corrected in the Fable-lane section
  (the scratch patch was half-applied); the basis-artifact assertions in
  `tJonesPupil` are KEPT and were not touched here.
* The scope list said the coated `Reflector` branch and `Refractor` were
  "unaudited".  Both are now audited: `Reflector` fixed in cb29ea5,
  `Refractor` normalized in 25c4386 with the transmitted-field A/B above.
