# Gate-review packet — external anchor for the protected-metal machinery

Opus lane. Closes the small worklist item Fable opened in
`REVIEW_POL_2C_2026-07-27.md`:

> **§3 last paragraph — the 151× stays model-relative; correct instinct.**
> The MgF₂/Al number additionally stacks the thin-film recursion through a
> quarter-wave overcoat — that is the part wanting an external anchor (a
> vendor retardance curve or a published protected-Al polarization-aberration
> result) before it appears in any budget document.

One item, one session, per the discipline rules. **No engine change** — the
comparison did not demand one, and §5 records the moment it looked like it
might.

---

## Verdict in one paragraph

The machinery is **anchored and correct**: driven with a publication's own
indices, film thickness, wavelengths and incidence angles, the engine
reproduces that publication's protected-aluminium model to **2.8e-14** in
diattenuation and **4.9e-14** in retardance, against a stated measurement
accuracy of **±0.01** per normalized Mueller element. The 151× is therefore
not corrupted by a thin-film bug — and, checked a second way, it is
*arithmetically right*. But the anchor exposed something else: **the
configuration it was measured on is not the configuration it is described
as**, and the design rule drawn from it has the **wrong sign**. That is the
finding; §4.

---

## What shipped

| Piece | Where |
|---|---|
| `pol_external_anchor` tool (6 files + README) | `mmacos/tools/pol_external_anchor/` |
| `tPolExternal` — 8 gates, added to `SUITE_FAST` | `mmacos/tests/` |
| 16 `external.json` entries with full provenance | `mmacos/tools/pol_validation_report/` |
| polval **section 8** + 7 gate-index rows | `macos/docs/macos-manual/polval/85_external_anchor.md.in` |
| polval **§5.5 correction notice** | `macos/docs/macos-manual/polval/60_contrast_floor.md.in` |
| `PLAN_POLARIZATION` worklist item closed | `macos/` |

**Suite:** `tPolExternal` 8/8, ~5.8 s at model 128. polval re-rendered, all
tokens resolve.

---

## 1. The anchor, and why this one

G. van Harten, F. Snik & C. U. Keller, *"Polarization properties of real
aluminum mirrors I. Influence of the aluminum oxide layer"*, PASP **121**,
377 (2009), doi:10.1086/599043; arXiv:0903.2740v1.

Chosen because it has **both** halves, which is rare: a full Mueller matrix
measured by ellipsometry at 14 incidence angles over 6–70° and four
wavelengths, *and* every model input stated numerically — both indices, the
film thickness, and the film model.

**The method is the point.** Rather than compare our numbers against a
published *result*, reproduce the publication's *configuration*: their
indices, their 4.12 nm oxide, their 500/550/600/650 nm, their angles.
`macos.coating` takes arbitrary n, κ and physical thickness, so this needed
no engine change and no fixture surgery. This isolates the machinery from
**index-table disagreement**, which is real and is not a machinery error —
the paper says the aluminium k values "vary widely throughout the
literature" and *fits* k rather than adopting a table. A comparison that
conflated the two would have been worthless in exactly the direction that
matters.

Their model (Eqs 1–6) is the Macleod / Born & Wolf characteristic-matrix
form. `vh_thinfilm.m` implements it from Macleod ch. 2 and their stated
equations, **never transcribed from `elemsub.F`** — the r_p-sign lesson.

Two conventions were settled rather than assumed, and one was settled by
evidence:

- **Index sign.** They use `N = n − ik`, k ≥ 0 = loss — the same convention
  the engine stores, so nothing translates. Load-bearing, not cosmetic: they
  report that the opposite sign made their fitted oxide come out ~50 nm
  instead of ~4 nm.
- **p̂ frame — measured, not derived.** The perfect-conductor case has a known
  answer, and the engine returns r_s/r_p = **+1.000000, imaginary part
  exactly 0**, at every angle. Bridge = 0. An earlier version of this harness
  *assumed* a π bridge from the ray-following doctrine and was wrong by
  exactly 180° everywhere; that is why this is now a gate
  (`test_phat_convention_bridge_is_zero`) and not a comment.
- **Their Eqs (5)–(6) print in an order PDF extraction scrambles**, leaving
  the s/p admittance assignment ambiguous. Settled by the *sign* of the [1,2]
  Mueller element: it must be positive for a metal, and their Fig. 1a plots
  it on a 0.00–0.15 axis. The Macleod-standard assignment gives 0.087–0.115
  at 70° across their wavelengths; the swapped one gives the wrong sign.
  Pinned by `test_admittance_assignment_sign`.

---

## 2. (a) The machinery check — PASS

Four wavelengths × six angles inside their measured range, compared **per
ray** at each ray's own incidence cosine (no single-AOI approximation):

| quantity | worst deviation | their accuracy |
|---|---|---|
| diattenuation ([1,2] element) | **2.828e-14** | ±0.01 |
| retardance, same Mueller units | **4.937e-14** | ±0.01 |

Both sides implement the same thin-film theory, so round-off agreement is
the *expected* result. The publication's ±0.01 is what the comparison would
have had to beat to mean anything, and is quoted that way rather than as an
achievement.

**Tolerance provenance, as requested.** ±0.01 is their Sec. 2: "an absolute
accuracy of ~1% of element [1,1]". Their fitted oxide error (±0.08 nm)
propagates to only 1.0e-4 in [1,2] and 2.6e-3 in the retardance block, so
their stated inputs determine the reference curves well inside their own
measurement bar — the comparison is not limited by input uncertainty.

---

## 3. Non-vacuity — the two ways it can fail are their own physics

- **Omitting the 4.12 nm oxide** (engine oxidized, analytic bare) disagrees
  by more than ±0.01. This independently reproduces the paper's central
  claim — that a ~4 nm film must be modelled at all.
- **The historical ~50 nm oxide**, which the paper attributes to a sign error
  on the imaginary index, is separable from 4.12 nm at their accuracy.

Plus two structural guards: the frame-free estimator cross-checks itself
(`M11/M22` vs `−M12/M21`, agreeing at ~1e-15), and the fold rig is asserted
to present the incidence angle it was asked for — see §5.

---

## 4. THE FINDING — the 151× is arithmetically right and verbally wrong

**This is the item that wants a second opinion.**

With the machinery anchored, §5.5's ladder becomes checkable *internally*.
For an on-axis rotationally symmetric train the mirrors' s/p axes share one
azimuth, so `J = c₀I + c₂·[cos2φ, sin2φ; sin2φ, −cos2φ]` and the
cross-polarized amplitude is set entirely by `ε_tot = c₂/c₀`; the cross power
ratio between coatings is |ε_tot|². Three independent routes agree:

| source | MgF₂/Al ÷ bare Al, cross power |
|---|---|
| independent analytic at the fixture's actual λ | 6.35 |
| the engine's own Jones pupil | 6.42 |
| `macos.pol_contrast_floor`, as §5.5 reports it | 5.42 |

So **the number is correct** (the spread is dark-zone annulus weighting vs a
pupil mean). What is wrong is what it was called.

*A note on what is deliberately absent.* I tried twice to quote the
fixture's per-element incidence angles and both routes are unusable through
this binding — `ray_field`'s `.nx/.ny/.nz` is the element **axis** normal
broadcast to the grid, not the local surface normal (tell: primary and
secondary return identical statistics), and the successive-element deviation
route returns exactly 90° everywhere. Rather than publish a third guess, the
claim is dropped: the ratio is flat at **6.34 … 6.36** across the whole
plausible near-normal range, so nothing here depends on it. `vh_cass_probe.m`
records both dead ends so the next person does not re-walk them.

`Rx_Cass_FarField` runs at **`Wavelen = 1.0E-06` m**. The 2c coating
constants are labelled `% Al at 632.8 nm` and `% ~quarter wave at 632.8 nm in
MgF2`. At the fixture's own wavelength, 110 nm of MgF₂ is **0.607**
quarter-waves — not a quarter wave, and not an ordinary protected-aluminium
recipe.

That distinction decides the **sign** of the design rule, because the
overcoat trade reverses across the quarter-wave condition. Same 110 nm film,
same aluminium, same few-degree incidence:

| film | optical thickness | cross power vs bare Al |
|---|---|---|
| 110 nm MgF₂ at 632.8 nm | 0.96 λ/4 | **0.18** |
| 110 nm MgF₂ at 1 µm *(the fixture)* | 0.607 λ/4 | **6.35** |
| true quarter wave, 181.2 nm at 1 µm | 1.00 λ/4 | **0.0157** |

A genuine quarter-wave overcoat **suppresses** the polarization floor by
about 1.8 decades. The report's sentence — that the protective overcoat
"costs most of another" decade — is **withdrawn**: true of the film as
configured, false of the recipe it was described as.

**Handling, for review.** Per the instruction not to reconcile silently:

- §5.5 keeps its measured numbers and gains a **correction notice** pointing
  at §8.3. The numbers are right; only the generalisation was wrong.
- §8.3 states the reversal with all three thicknesses.
- `tPolExternal/test_overcoat_trade_reverses_with_optical_thickness` pins it
  — the same film reducing the driver at 632.8 nm and raising it at 1 µm,
  plus a true quarter-wave suppressing it — so this cannot silently regress.
- **`tPolContrast`'s asserted 27.898 / 151.31 were NOT touched.** They are
  correct measurements of the configuration as built. Changing the fixture
  would be a different item, and it would move gated numbers.

**Recommended follow-on (not done here, deliberately):** decide whether the
2c coating ladder should be re-run with a self-consistent configuration —
aluminium indices at the fixture's own wavelength and a genuine quarter-wave
overcoat — or whether `Rx_Cass_FarField` should move to 632.8 nm. That is a
fixture-and-gated-numbers decision, which is a Fable call, not a mechanical
one.

---

## 5. Two harness bugs worth recording — both mine, both instructive

Neither reached a deliverable, but both are the *class* of error this
programme keeps paying for, and both were caught by physics rather than by a
test.

1. **The fold rig swept the complement of the intended angles.** Mirror
   deviation is `180° − 2·AOI`, not `2·AOI`. Getting it backwards is
   **self-cancelling at exactly 45°** — the one angle the pre-existing
   Fresnel gate runs at — so nothing in the suite would have caught it, and
   the first run "agreed" at 45° and nowhere else. Now gated
   (`test_fold_rig_aoi_is_what_was_requested`).

2. **The single-trace extraction hard-codes the launch frame.** tJonesPupil's
   Fresnel gate divides out the input state's projection on the s and p axes,
   which needs the engine's per-ray launch frame — in practice a hard-coded
   `xGrid = [0;1;0]`. Also safe at 45° and quietly wrong elsewhere. The
   symptom was diattenuation coming out nearly **flat in AOI** (−3.1e-3 at 2°
   vs −4.7e-3 at 45°) where an isotropic surface must give D ∝ θ² — the same
   "flat where physics demands a power law" signature that exposed the r_p
   sign defect. Replaced with a **frame-free two-trace construction**:
   `M = diag(r_s,r_p)·R(φ)` ⇒ `r_s/r_p = M11/M22 = −M12/M21`, φ cancelling
   identically and the two estimates cross-checking each other.

**The near-miss worth flagging.** Between those two bugs, the harness
produced a coherent-looking story in which the engine contradicted the
literature and the 151× was inverted — a plausible "engine defect" narrative
built on two of my own errors. What broke it was not a test but a units
check: noticing `Wavelength= 1.0000000D-06` in the trace log while the
coating constants said 632.8 nm. Worth recording that the *first* coherent
explanation was the wrong one, and that the discipline rule which saved it
was "stop and check the configuration before believing a contradiction".

---

## 6. What this does NOT claim

- It compares against the publication's **model** at their stated inputs.
  Their measured matrices are published as *figures*; digitising them would
  inject a transcription error larger than anything being measured. Their own
  model-to-measurement agreement (±0.01) is the outer envelope and is quoted
  as such.
- One wavelength region, one metal, one dielectric. It anchors the
  characteristic-matrix recursion on dielectric-on-metal — not absorbing
  overcoats, many-layer stacks, or the transmission side (§7 covers that).
- **The aluminium indices used in §5.5 remain a model choice**, not an
  anchored one. This pins the machinery, not the material data.
- Not re-run this session, and stated as such: pymacos/ifx, PROPER-compare,
  GMI. Nothing here touches the engine — the tool and gate are pure
  binding-layer MATLAB — so those cannot have moved.

---

## What to look at first, if tokens are short

1. **§4** — the 151× disposition. The number survives; its description and
   the design rule drawn from it do not. The decision it invites is whether
   to re-configure the 2c fixture (and move gated numbers) or leave it as a
   documented one-off datum.
2. **§5's near-miss** — two of my own errors briefly produced a convincing
   "the engine disagrees with the literature" story. Worth a sanity check
   that §4's conclusion is not the same failure one level up.
3. **§1's p̂ bridge** — measured as 0 where the ray-following doctrine
   predicts π. That is worth an independent look; the machinery agreement is
   exact either way for diattenuation, but retardance depends on it.
