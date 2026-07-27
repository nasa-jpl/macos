# Review packet — Phase 3 polarizing elements (polarizer + waveplate)

**Date:** 2026-07-27 · **Branch:** `pol-core` (both repos) · **Lane:** Opus
(worklist item 4) · **Reviewer ask:** one convention decision, plus a normal
line review of the diff.

---

## What landed

`TrPolarizerElt(15)` — an ideal linear polarizer, finished from a
name-table-only stub — and `WavePlateElt(18)`, a new linear retarder
(`mEltTypes` 17 → 18). Both are transmissive, both ride a new `PolElt`
routine in `elemsub.F`, both are gated on `ifPol`.

| Layer | What |
|---|---|
| Engine physics | `PolElt` (`elemsub.F`), modelled on `RefSrf` verbatim for geometry |
| Trace dispatch | `tracesub.F` **and** `propsub.F` (see the finding below) |
| Storage | `PolAxisElt(3,mElt)`, `PolRetElt(mElt)`, and `JmatElt(2,2,mElt)` — previously dead — now filled |
| Rx keywords | `PolAxis=`, `Retardance=` in the `msmacosio.inc` element chain + the `macosio.F` MOD chain |
| Defaults | `ChkDf2` requires both on the element types that use them (strict, on purpose) |
| SAVE | writer inverts the load-time `Wavelen` scaling, like the `Coating=` writer |
| API | `polelt_set` / `polelt_get` / `jmat_elt_get`, all codegen Path A |
| mmacos | `macos.polarizer`, `macos.waveplate`, `macos.elt_jones` + `Session` methods |
| Gates | `tPolElement`, 23 tests, SUITE_FAST; fixtures `Rx_PolElt.in` + `Rx_PolElt_Ref.in` |
| Docs | manual §4 "Polarizers and Waveplates", cmdref NOTES ×3, both `CLAUDE.md`s |

`RfPolarizerElt(14)` is deliberately **still a stub** — see the decision
request.

---

## THE DECISION REQUEST — the off-normal axis convention

This is the one item I stopped short of settling, per the standing rule that
convention arbitration in the polarization core goes to the Fable lane.

**The question.** An ideal polarizer is defined by one axis in the plane
transverse to the ray. The element's axis, however, is declared in global
coordinates and must be *projected* into that plane. Off normal incidence,
two constructions that agree exactly at normal incidence diverge:

* **(A) — what I implemented.** Declare the PASS axis; project it
  orthographically, `â = unit(a − (a·r̂)r̂)`; the blocked direction is
  `b̂ = r̂ × â`.
* **(B) — the alternative.** Declare the BLOCK axis (the wire direction of a
  wire-grid polarizer, say); project *that*; transmit the complement.

They differ because orthographic projection does not preserve orthogonality.

**Measured size.** With the ray tilted by `a` and the declared axis at
azimuth `φ` from the plane of incidence, the angle between the two candidate
pass axes has the closed form, at `φ = 45°`,

```
Δ = acos( 2 cos a / (1 + cos² a) )
```

which is `0.90°` at 10° AOI, **`3.56°` at 20°**, and `7.9°` at 30°; bounded
above by `sin² a` throughout. It is **identically zero** at normal
incidence, and — the part worth flagging — **also identically zero when the
declared axis lies in, or perpendicular to, the plane of incidence.** The
obvious way to test this (axis along x, tilt in the x–z plane) is exactly
that degenerate case; my first version of the gate did precisely that and
reported a clean zero, which is how I noticed. `Δ(aoi, φ)` is now pinned in
`tPolElement/test_offnormal_convention_magnitude`, including the two
vanishing cases, so the number is on the record either way the decision goes.

**What I did NOT do.** I did not pick (A) as settled physics. Everything
gated in this landing is at strict normal incidence, where the two
constructions coincide exactly and no choice is being exercised. The
docstring, cmdref NOTES, manual section and engine `CLAUDE.md` all state
which construction the code uses and that it is provisional.

**Be clear about what the gate measures.** `test_offnormal_convention_magnitude`
computes both constructions in MATLAB from the geometry; it does **not** drive
the engine off normal. That is deliberate — it is a statement about the size of
the ambiguity, not a measurement of the engine — but it means the off-normal
code path is **unvalidated as well as unsettled**. No fixture in the suite
sends an off-normal ray through a polarizing element at all. Whichever way the
convention goes, landing it should come with a fixture that actually tilts the
beam (a converging source, or an off-axis field point — tilting the *element*
does not tilt the ray in a collimated on-axis bundle, which is why the existing
fixture cannot be adapted).

**Why `RfPolarizer` is still a stub.** A reflective polarizer is *inherently*
off-normal — at normal incidence it would send light back along the input
ray. So it cannot be implemented at all without first settling the above.
Implementing it "for completeness" would have meant landing exactly the guess
the rule exists to prevent.

**My read, offered as input, not as a decision.** (A) is the better default:
it makes the keyword mean the same thing for a polarizer and a waveplate (the
transmitted / fast axis), and the physically-motivated case for (B) — a wire
grid — is arguably not what an "ideal polarizer" element should be modelling
anyway. But there is a real argument that a tilted ideal polarizer should not
be offered at all, and that anything at a substantial AOI belongs on a coated
`Reflector`/`Refractor`, where s and p are set by the physical plane of
incidence and the thin-film recursion is already gated against Fresnel. If
that is the call, the right change is to make `PolElt` *refuse* rays beyond a
small AOI rather than silently apply rule (A) — a small edit, and I would
rather do it than have the rule quietly become precedent.

---

## THE FINDING — there are two element dispatch chains, and I only wired one

Worth a paragraph because it is a general trap, not a detail of this slice.

`tracesub.F` dispatches elements for the ray trace. `propsub.F`
(`CPROPAGATE`) **re-traces the rays that seed the diffraction grid through
its own, separate `ELSE IF (EltID(iElt).EQ.n)` chain.** I wired only
`tracesub.F` first. The result passed every ray-level gate — Malus, crossed
extinction, QWP, HWP, composition, unitarity, all of it — and was
**invisible to `intensity` / `complex_field`**.

Measured, before the fix: with crossed polarizers, the ray power at the
detector was `3.6e-33` while the detector *plane* held the full `9.69e-01`.
That is the failure mode to be afraid of: it does not look like a crash or a
wrong number, it looks like "polarization has no effect on the image", which
is a plausible thing for a reviewer to believe.

Caught only because I wrote `test_grid_carries_the_polarizing_train` as a
*scope* claim about Tranche 1 — I expected it to pass and was checking that
the element ordering in the fixture was right. It is now the tripwire; any
future field-touching element needs both chains.

(`srtrace.F` has a third chain, but it is inside `SRTRACE_Test` under
`#if 0` — dead, consistent with the Phase-0 finding. Nothing to add there.)

---

## Conventions — all four extend the pinned set; none is new law

1. **Axis as a 3-vector**, not an angle in an element frame. Sidesteps
   inventing a "which in-plane direction is zero degrees" convention. Follows
   the plan's own wording ("keyword for the transmission-axis vector").
2. **Orthonormalize the partner axis** (`b̂ = r̂ × â`) rather than projecting a
   second declared axis. This one is **forced, not chosen**: a lossless
   retarder must be unitary, and only an orthonormal eigenbasis makes
   `diag(1, e^{-iδ})` unitary. Projecting two axes independently gives a
   non-orthogonal pair and a non-unitary "retarder".
3. **Retardance sign read off the engine.** `C1 = exp(−i·2πLN/λ)`
   (`elemsub.F:395`) ⇒ the slow axis takes the more negative phase ⇒ with the
   declared axis as FAST, `J = diag(1, e^{-iδ})`. Not legislated; derived
   from the same line the conventions table already cites.
4. **Retardance stored physically** as `(n_slow−n_fast)·d`; the Rx value is in
   waves at parse-time `Wavelen` and is scaled at load, the trace divides by
   the current λ. Identical treatment to `Coating=` thickness — a plate is
   fixed glass, and a wavelength sweep is chromatic.

---

## Gates — 23, and how each avoids being vacuous

Every physics prediction is written from the Jones algebra in the test
comments, **not transcribed from the engine**. That is the standing lesson
from the `r_p` sign defect: the Fresnel fold gate could not catch a sign
error because its "analytic" came from the engine's own expression.

| Gate | Claim | Result |
|---|---|---|
| `malus_law` | `I(θ) = I₀cos²θ`, 13 angles | 1e-12 of `I₀` |
| `crossed_polarizer_extinction` | exactly-orthogonal axes → **exactly 0** | `0` (all three components) |
| `qwp_linear_to_circular` | `S1 = S2 = 0`, **`S3/S0 = −1`** | 1e-14 |
| `qwp_circular_to_linear` | QWP→QWP: DoLP 1, DoCP 0, *and* the intermediate really circular | 1e-14 |
| `hwp_rotates_by_2theta` | orientation `= 2θ`, fitted slope | slope `2` to 1e-10 |
| `two_qwp_equal_one_hwp` | cascade identity, field for field | 1e-15 |
| `waveplate_is_unitary` | power conserved, `J′J = I`, linear ×3 + circular | 1e-14 / 1e-15 |
| `zero_retardance_is_identity` | `e^{-i0} = 1` exactly | 1e-17 |
| `retardance_sign_matches_engine` | `arg(Ey/Ex) = −2πR` | 1e-12 |
| `retardance_is_chromatic` | R halves when λ doubles | 1e-12 rel |
| `unpolarized_bit_identical_to_reference_twin` | pol-off ≡ `Reference` surfaces | **bit-identical** (OPD, intensity, ray status) |
| `grid_carries_the_polarizing_train` | detector plane obeys Malus too | 1e-10 |
| `offnormal_convention_magnitude` | the ambiguity above, incl. both zero cases | closed form to 1e-12 |
| + 10 plumbing gates | keywords, round-trip, SAVE, guards, dirty-trace | — |

**Non-vacuity is measured in-suite, not argued**, on the principle that the
pre-Phase-3 behaviour of these EltIDs was "silently do nothing" — an element
that still did nothing would pass a badly-built suite:

* Malus: `std(I)/mean(I) > 0.5` and `min/max < 1e-25` — a pass-everything
  element gives a flat curve and fails both.
* QWP: the same rig at `R = 0` gives `S3/S0 = 0` and `S1/S0 = 1`.
* HWP: at `R = 0`, rotating the plate gives fitted slope `0`, not `2`.
* Composition: a *single* QWP at the same axis differs from the pair by >10%.
* Unitarity: the polarizer's Jones is checked to **fail** `‖J′J − I‖ > 0.5`,
  so the unitarity test cannot pass on an engine that stopped applying
  anything.

**Circular input states are load-bearing**, for the reason `tPolContrast` and
`tVecChain` established: linear-only suites are blind to conjugation and
retardance-sign errors. The signed `S3/S0 = −1` is the sharpest of these —
it flips outright if the retardance convention flips, where a `|S3|` gate
would accept either.

---

## Scope stated honestly

* **Normal incidence** is the gated regime, for the convention reason above.
* **Tranche 1 still applies.** A polarizing element placed *after* the first
  physical-optics leg transforms rays and never reaches the diffraction grid.
  `Rx_PolElt.in` therefore puts all four polarizing elements before its single
  `NFPlane`→`Geometric` leg, and its header says why. `Rx_Coro`-style layouts
  with a polarizer between legs need Tranche 2 (`J_run`), unchanged.
* **Idealizations:** no ray splitting (a PBS is two traces, or better a coated
  `Reflector` at 45°), no walk-off, no face Fresnel loss, no substrate. The
  output is purely transverse — any longitudinal component at the surface is
  discarded, which is what a 2×2 Jones element means (exactly zero for the
  collimated normal-incidence fixture, O(NA) otherwise).
* **Not attempted:** VVC, `RfPolarizer`, the polarization phase-shifting
  Twyman-Green example (that is Phase 2d / `pol-ifo`, and needs the Bench
  `add_polarizer`/`add_waveplate` emitters, which are on the other branch).

---

## Regression

| Gate | Result |
|---|---|
| `tPolElement` | **23 pass, 0 fail** |
| mmacos full suite (gfortran) | **463 pass, 0 fail** — fast 312 / masks 62 / freeform 60 / proper-512 10 / pol-512 6 / proper-1024 13 (was 440; +23 new) |
| ifx build (`makems.sh release`) | **SUCCEEDED** — the standing ifx smoke on new fixed-form code |
| pymacos suite (ifx-linked) | **6652 pass, 0 fail** |
| pymacos PROPER comparison | **all phases passed**; committed residuals reproduced |
| GMI regression | **6/6, `vs-ref = 0.000e+00`** — bit-identical, as required for an ifPol-off consumer |
| polval regen | all gate thresholds pass, including the 15 new ones at model 128 |

**The polval regen is itself a no-op proof.** Regenerating the whole report
across this change moved **no existing number**: 119 → 136 tokens, 17 added,
none removed, **none changed**; and all ten pre-existing figures are
pixel-identical (byte differences are PNG metadata only). Every published
Phase 0–2c result stands.

pymacos bindings for the three new API routines are **not** in this landing
(mmacos-only, as with the `view_rx` engine leg). The routines are in
`macos_api_mod`, so the pymacos side is a shim addition whenever it is wanted.
The pymacos suite above is therefore a regression check that the engine change
is inert for that binding, not coverage of the new elements.

---

## Fable-lane decision (2026-07-27 pm, appended post-review)

Independent verification: tPolElement re-run 23/23 on this box; PolElt
engine block line-read (geometry, the four conventions, both dispatch
chains); packet cross-checked against the diff.  Line review: PASS, no
findings.  The four conventions are correctly derived rather than
legislated — (2) unitarity-forced orthonormalization and (3) the
retardance sign read off elemsub.F:395 are exactly right.

**THE DECISION — neither (A) as landed nor refusal.  The rule is:
PROJECT THE MATERIAL AXIS.**

The invariant physical object in a tilted polarizing element is the
material direction fixed in the element — the absorbing direction (wire
grid: the wires; dichroic sheet: the dipole chains) for a polarizer, the
crystal fast axis for a waveplate.  The transverse-plane behaviour
follows by projecting THAT vector:

* **Waveplate: the landed implementation is CORRECT as is.**  The
  declared axis IS the material (fast) axis, so projecting the declared
  axis — construction (A) — is the material-axis rule for this element.
  No change.
* **Polarizer: flip to the block-axis projection, keyword semantics
  UNCHANGED.**  `PolAxis=` keeps declaring the transmission axis (the
  natural user handle, and Opus's uniformity argument survives intact).
  The material absorbing direction is its in-element-plane complement,
  `ŵ = unit(ψ̂ × â_pass)`; per ray, project ŵ into the transverse plane,
  extinguish it, transmit `r̂ × ŵ_proj`.  At normal incidence this is
  identical to (A), so every gate in this landing stands bit-for-bit.

This is not a taste call — the exact ambiguity has been adjudicated by
experiment: Korger et al., *"The polarization properties of a tilted
polarizer,"* Opt. Express **21**, 27032 (2013), measured a tilted
dichroic polarizer and found it follows the dipole (material-axis)
model, not the geometric pass-axis projection.  Pull the paper during
the landing and anchor the polval section on it (external.json entry,
same pattern as the other external anchors).

**Refusal (the packet's option C) is rejected**: a polarizer in a
converging beam legitimately sees a few degrees of AOI (f/10 → ~3°,
Δ ≲ 0.1°), and refusing would break those layouts to guard against an
ambiguity that is now settled.  The coated-Reflector recommendation for
LARGE deliberate AOI (45° PBS-style) stays as written in the scope
notes — that advice is right regardless.

**Landing spec (next Opus slice — now fully templated):**
1. Flip `PolElt`'s polarizer branch to the block-axis projection above
   (waveplate branch untouched).
2. An off-normal fixture that actually tilts the BEAM (converging
   source or off-axis field — the packet is right that tilting the
   element does nothing in a collimated on-axis bundle).
3. Engine-driven A/B gate: the Δ(aoi,φ) closed form already pinned in
   `test_offnormal_convention_magnitude` becomes an engine-vs-closed-form
   assertion (drive the engine at 20° AOI, φ = 45°, assert the
   transmitted axis matches the material-axis prediction and differs
   from (A) by the predicted 3.56°).  Include the two degenerate-zero
   cases as vacuity guards, as the MATLAB-side gate already does.
4. Normal-incidence gates must stay bit-identical (they should by
   construction — assert it).
5. Docs flip from "provisional (A)" to the settled material-axis rule,
   everywhere the provisional flag is currently visible.

**`RfPolarizer` stays a stub** — the convention no longer blocks it,
but a reflective wire grid carries additional physics (grid reflection
efficiency, s/p of the substrate) beyond the axis rule, and nothing
needs it yet.  The loud load failure stays; do not build it for
completeness.

**The dual-dispatch-chain finding — endorsed, and it is the session's
most valuable output.**  "Passes every ray-level gate, invisible to the
image, reads as 'polarization has no effect'" is the exact failure
shape the r_p defect taught us to fear, caught this time by a scope
claim written to PASS — which is the review methodology working.  The
tripwire test + both-CLAUDE.md recording is the right durable form; any
future field-touching element needs tracesub AND propsub, and now the
record says so.

**Idealization note to add at the waveplate (one line, ride the next
landing):** the ideal retarder ignores the AOI-dependence of retardance
in a real crystal plate (field-of-view effects); anyone needing that
needs a birefringent-plate model, out of scope.

**Ack, no action:** the signed S3/S0 = −1 gate (a |S3| gate would
certify a conjugated retarder — same lesson as the coherency bug, now
applied preemptively); the Reference-twin bit-identity gate; non-vacuity
A/Bs against the "silently do nothing" prior; pymacos shim deferral
correctly scoped as regression-only.

---

## Landing — the material-axis rule, implemented and gated (2026-07-27 pm, Opus)

**macos `216c56c` + MACOS_resources `a3417ce`, both on `pol-core`, pushed.**
All five points of the landing spec are in. The decision itself was not
re-opened.

### 1. The flip

`PolElt`'s basis construction now selects the vector to project by element
type. For a `WavePlate` it is the declared (fast) axis — unchanged, since
that IS the material axis. For a `TrPolarizer` it is `mvec = psi x
PolAxis`, the absorbing direction; that is projected into the transverse
plane, becomes `bhat` (extinguished), and the transmitted `ahat = bhat x
rhat` completes it. Keyword semantics are untouched: `PolAxis=` still
declares the transmission axis.

This is exactly Eq. (5)–(6) of Korger et al. — `â = unit(P_A − (P_A·ẑ)ẑ)`,
`E_out = E_in − (E_in·â)â`, i.e. `T_A = 1 − ââᵀ = t̂t̂ᵀ` with `t̂ = ẑ × â`.
Worth recording from reading the paper: **both** constructions are in the
literature, and the one that shipped first is Fainman and Shamir's
transmitting-axis projection (*Appl. Opt.* **23**, 3188, 1984), which
Korger et al. cite as the model they are testing against. Their §3 states
the two "coincide only for normal incidence" and their Fig. 2(c) is the
measurement: *only the absorbing model explains the drastic change of the
transmitted state of polarization observed when the polarizer is tilted*.
So the rejected option was a defensible published model, not a slip — which
is the reason the ambiguity was real enough to need arbitration.

Two consequences of `psi` entering the construction, both documented at
every surface: a polarizer's declared axis is now taken **modulo its
component along the element normal**, and a pass axis **parallel** to the
normal extinguishes (it leaves no in-element blocked direction). The
pre-existing degenerate case — material axis along the ray — is unchanged
and shares the same guard.

### 2–3. The fixture, and the engine-driven A/B

`tests/Rx/Rx_PolElt_Tilt.in` tilts the **beam**: `ChfRayDir = (sin20, 0,
cos20)`, source frame following it, two `TrPolarizer` elements left at
`psiElt = −z` so they see exactly 20° while every other element stays
normal to the beam. The packet's warning was right — tilting the element
would have measured nothing.

Section F of `tPolElement`, four new gates (23 → 27):

| Gate | Result |
|---|---|
| AOI is 20° ray by ray, from the engine's own direction cosines and normal | 7.1e-15° residual |
| transmitted axis **is** the material rule (all 12453 rays) | 1.2e-16 rad |
| and **misses** the pass-axis projection by the closed form | 3.5616°, matching to 1.3e-14° |
| detector-plane power at the material-rule crossed null | 9.1e-33 relative |
| detector-plane power at the pass-axis crossed null | 1.526e-02, = the predicted cos² to 4.5e-16 |
| degenerate azimuths (0°, 90°): constructions agree, engine on both | 7.6e-18 / 5.6e-17 rad |

(the numbers are the polval driver's, which re-measures all of it
independently of the test class; the test class's own tolerances are looser
by design)

Two remarks on how these were built. First, the transmitted axis is read
**from the field** (a polarizer's output is `t̂` times a complex scalar) with
the ray direction and surface normal taken from `rayfield_get`, so nothing
assumes the fixture's declared numbers — the discipline the odd-mirror PEC
gate established. Second, the angle is computed as `atan2(|cross|,|dot|)`,
not `acos(dot)`: the first version used `acos` and reported 8.5e-7 rad for a
residual that was really 1e-16, which would have forced a meaninglessly loose
tolerance.

**The grid-side gate needed a different idea than the ray-side one.** A
component-ratio comparison on the detector planes would have imported a
question about their frame. Instead it uses the crossed-analyzer null: off
normal, the analyzer azimuth that extinguishes is a *different number* under
the two rules — 131.445° vs 138.555°, 7.11° apart — so pointing the engine at
each in turn discriminates them with **total intensity alone**. The
separation is ~30 orders of magnitude, and the leak at the wrong setting is
asserted against its predicted `cos²`, so it is a number rather than a
"greater than zero".

**Non-vacuity, measured not argued.** The pass-axis engine was rebuilt (both
`EltType.EQ.TrPolarizerElt` basis branches made unreachable, full
`makems.sh` + mex relink) and the suite re-run: the two new engine gates
**fail** — transmitted axis off by the full 3.5616°, and the null/leak
values swap to 1.53e-2 / 1.86e-32. All 23 original gates pass on **both**
engines, as does the degenerate-azimuth gate — which is precisely its job.
Captured in `external.json` as `X_MATAXIS_PREFIX_*`.

### 4. Normal incidence is bit-identical, and now asserted

Values were captured from the pre-flip engine *before* the flip (macos
`52c7669`) and are asserted literally in
`test_normal_incidence_unchanged_by_the_material_flip`: the per-ray complex
field after polarizer→two waveplates→analyzer at generic axes, and a
four-point Malus curve at the **detector** — all matching to all 17
significant digits, ray side and grid side.

That is not automatic. `ahat = bhat x rhat` returns the declared axis to the
last bit only if nothing renormalizes it afterwards, so `PolElt` deliberately
omits the `DUNITIZE` that the waveplate branch applies to its partner axis
(the vector is already unit to a rounding error; the redundant division would
perturb it by an ULP). The reason is now a comment in the routine, because it
looks like an oversight.

### 5. Docs

Every provisional flag is gone: engine `PolElt` header + `elt_mod.F`,
`macos_f90/CLAUDE.md` (conventions list renumbered to five, the OPEN block
replaced by the settled rule), `mmacos/CLAUDE.md`, `PLAN_POLARIZATION.md`
(worklist item 4 and the Phase 3 status block), the two veneer docstrings,
cmdref NOTES ×2 (regenerated), manual §4, `tests/README.md`, and polval §6.7
— rewritten from "the one convention question" to the settled rule plus its
engine gate, with the frontmatter and coverage sections following. The
waveplate idealization note (retardance does not vary with AOI; a real
crystal plate's does — the field-of-view effect) is added at the engine
header, the veneer, cmdref, the manual, and the polval idealization list.
Korger et al. is in `external.json` as `X_KORGER`, with the Fainman–Shamir
citation recorded alongside it as the rejected construction.

`RfPolarizer` stays a stub, and its *reason* changed: no longer the axis
convention, now grid reflection efficiency and the substrate's own s/p
response. Stated that way in all five places that mentioned it.

### Regression

| Gate | Result |
|---|---|
| `tPolElement` | **27 pass, 0 fail** (23 → 27) |
| mmacos full suite (gfortran) | **467 pass, 0 fail** — fast 316 / masks 62 / freeform 60 / proper-512 10 / pol-512 6 / proper-1024 13 |
| ifx build (`makems.sh release ifx`) | **SUCCEEDED** — the standing ifx smoke |
| GMI regression | **6/6, `vs-ref = 0.000e+00`** |
| pymacos suite (ifx-linked) | **6694 passed, 1 skipped** — identical to the pre-flip capture |
| polval regen | all gate thresholds pass, including the 10 new ones at model 128 |

**The polval regen is again a no-op proof, and this time it is the
bit-identity claim at report scale**: 136 → 146 tokens, 10 added, **none
removed, none changed**. Every published Phase 0–3 number is untouched by the
flip, which is what "identical at normal incidence" has to mean if it means
anything.

The engine was rebuilt for both compilers after the last (comment-only) edit
and `tPolElement` re-run on that final binary, so the shipped mex and the
committed source are the same code.

**What the cross-binding suites do and do not cover here.** No fixture in
pymacos or GMI instantiates `TrPolarizer` or `WavePlate` — checked, not
assumed — so `PolElt` is unreachable from either and this change cannot
affect their results. Both were re-run anyway and are unchanged; they confirm
inertness rather than covering the new behaviour, and the pymacos bindings for
the three polarizing-element API routines are still deferred as in the first
landing. The **PROPER-comparison driver** (`run_proper_tests.sh`) was *not*
re-run separately this time, for the same reachability reason; its committed
residuals are the ones captured at the `r_p` landing earlier today.

One process note, since it cost a wasted 9-minute run: `pymacos/tests/context.py`
resolves `src/` as `Path(".").absolute().parent/'src'`, so pytest **must** be
invoked from `tests/`. From the repo root every test errors on
`ModuleNotFoundError: No module named 'pymacos'` — 6584 "failures" that are
entirely the invocation. `external.json`'s command string now says so.
