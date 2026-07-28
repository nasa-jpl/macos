# REVIEW — the overcoat quarter-wave reversal, measured on both sides

**Date:** 2026-07-28 · **Branch:** `pol-core` (both repos) · **Lane:** Opus
**Scope:** evidence run only. **No engine change. No fixture change. No
previously-gated number moved.**

Closes the follow-on Fable opened in `REVIEW_POL_EXTERNAL_2026-07-28.md`.

---

## 1. What was open

The external-anchor packet found that the 110 nm MgF₂ overcoat the Phase-2c
coating ladder applies is **0.607 quarter-waves** at `Rx_Cass_FarField`'s own
1 µm — not the 0.96 that its `% ~quarter wave at 632.8 nm` comment describes
— and that the overcoat polarization trade **reverses** across the
quarter-wave condition. The design-rule sentence in §5.5 inverted as a
result.

That correction rested on an **independent analytic** (`vh_thinfilm`). The
ruling was: gated fixtures do not move; instead put **engine** numbers on both
sides of the reversal, via a companion run at 632.8 nm, so the chromatic
statement is a measurement rather than an assertion.

---

## 2. What was done

A companion ladder on the **unmoved** `Rx_Cass_FarField`, model 256,
x-polarized, both mirrors coated. `Wavelen` stays `1.0E-06` on disk; the
companion wavelength is applied at **runtime** with `macos.set_src_wvl` after
`load_rx`. Four ladder points per wavelength — uncoated → 200 nm Al →
110 nm MgF₂/Al → a **true** λ/(4·1.38) MgF₂/Al — plus one control.

**This is a measurement, not a unit conversion.** `macos.coating` takes
*physical* thickness and the engine divides by the *current* `Wavelen` when
it applies the layer phase, so a film is fixed glass under a wavelength
change and its optical thickness moves. That is the entire mechanism under
test, and the control in §4 is what proves it is the mechanism.

Cross-polarized **total power** is the reported quantity, never an annulus
mean: a fixed pixel annulus subtends a different λ/D range at the two
wavelengths, so an annulus statistic would mix the coating effect with the
diffraction scale.

Two ratios are reported because they answer different questions:

* **total** — MgF₂/Al ÷ bare Al on the full cross-polarized power. What a
  designer sees; it carries the irreducible **geometric** cross term.
* **excess** — the same ratio formed from the coating excess over the
  uncoated baseline (identically the ratio of the two `d_cross_rel`). The
  coating-only quantity, and the one directly comparable with a pure-Fresnel
  analytic.

---

## 3. Result — the reversal is engine-measured on both sides

| run | film, quarter-waves | MgF₂/bare, **total** | MgF₂/bare, **excess** | analytic |
|---|---|---|---|---|
| 1 µm — the fixture's own λ | 0.6072 | **5.2707×** *costs* | 5.4238× | 6.35× |
| 632.8 nm — the companion | 0.9595 | **0.2035×** *suppresses* | 0.17496× | 0.18× |
| a true λ/4, either λ | 1.0000 | **0.053207×** | 0.019269× | 0.0157× |

**The reversal, as one number: 25.899×.**

**The 1 µm side reproduces the recorded ladder exactly, in the same
harness:** `d_cross_rel` comes back at **27.8977** and **151.311** against
`tPolContrast`'s asserted 27.898 / 151.31. That is what makes the companion
run comparable rather than a nearby but separate measurement.

**Agreement with the independent analytic** (excess column, which is the
like-for-like one): −15% at 1 µm, **−2.8%** at 632.8 nm, +23% at the true
quarter wave. All three land on the correct side. The spread is the pupil's
AOI range (1.3–5.1°) against the analytic's single near-normal angle, and it
is widest at the quarter-wave point because that is where the coating term is
nearly nulled and a ratio near a null is most sensitive to angle spread.

**Why the engine's true-λ/4 total does not follow the analytic down.** At
that point the coating term is nearly extinguished and what remains is the
**geometric** cross term, which no coating removes: the true-λ/4 total sits
at **1.5376×** the *uncoated* floor. The overcoat has taken the coating
contribution down to the geometry, and the total cannot go below it. The
excess column, which subtracts that floor, is the one that tracks the
analytic — and does.

---

## 4. Non-vacuity — the reversal can be made to disappear, and the gate notices

The counterfactual is specific and physical, not a knob: an engine that
pinned coating thickness in **waves** rather than in metres — i.e. that
evaluated the "632.8 nm" coating constants *achromatically* — would at
632.8 nm be tracing 110 nm × (632.8/1000) = **69.608 nm** of MgF₂, the film
with the same optical thickness in waves that 110 nm has at 1 µm. It is
reachable from the real engine by asking for that film.

Measured, that control gives **5.2707×** — it lands back on the 1 µm answer
to **2.1e-08** and shows **no reversal at all**. Run against it
(`oc_nonvacuity`), **3 of 3** of the gate's 632.8 nm assertions fail:

```
  [FAIL] ratio_mgf2 < 1 (the film must SUPPRESS at 632.8 nm)  got 5.27071, wanted 1.00000
  [FAIL] ratio_mgf2 == 0.20351 within 2%                      got 5.27071, wanted 0.20351
  [FAIL] reversal == 25.899 within 2%                         got 1.00000, wanted 25.89900
```

Three further structural identities localise the effect **in the film and
nowhere else**. Each is exact by construction — the coating coefficients
depend on λ only through `n·d/λ`, and the indices are held fixed — so each is
pinned at round-off rather than at the 2% measurement tolerance:

| identity | measured |
|---|---|
| the metal-only leg (no film) does not move with λ | 1.3e-08 |
| the quarter-wave condition is λ-invariant (181.2 nm at 1 µm vs 114.6 nm at 632.8 nm — different glass, same answer) | 1.9e-08 |
| the uncoated geometric floor is λ-invariant (`cross_over_co`) | 1.6e-15 |

---

## 5. Tolerances, and why they are what they are

* **RelTol 0.02** on the four ladder pins (5.2707, 0.20351, 25.899,
  0.053207). This is `tPolContrast`'s existing convention for coating numbers
  — the 27.898 / 151.31 already asserted there carry the same tolerance —
  chosen to survive a BLAS/compiler reshuffle. It is three orders of
  magnitude tighter than the factor-26 effect being gated, so it cannot mask
  a physics change.
* **RelTol 1e-6** on the three invariances and on the achromatic control's
  agreement with the 1 µm answer. These are algebraic identities, not
  measurements; the observed residuals are ~2e-8, so 1e-6 leaves ~50×
  headroom for round-off without admitting anything real.
* **Bare inequalities** (`< 1` at 632.8 nm, `> 1` at 1 µm, `> 1` for the
  control) carry the qualitative claim independently of every pin, so a
  future retune of the pins cannot quietly delete the finding.

---

## 6. What landed

**macos** (`pol-core`)
* `docs/macos-manual/polval/85_external_anchor.md.in` — new **§8.3.1**, the
  engine measurement on both sides, the three invariances, the non-vacuity
  control, the reconciliation with the analytic, and the design rule restated
  chromatically. §8.3's table gains an engine column beside the analytic one.
* `docs/macos-manual/polval/60_contrast_floor.md.in` — §5.5's correction
  notice now cites **engine** numbers on both sides and states the design
  rule (specify overcoats at λ/4 of the working wavelength).
* `docs/macos-manual/polval/90_coverage.md.in` — three §8.3.1 rows in the
  gate index.
* Rendered `.md` regenerated; guards green.

**MACOS_resources** (`pol-core`)
* `mmacos/tools/pol_overcoat_chromatic/` — `oc_ladder.m` (the measurement),
  `oc_nonvacuity.m` (the counterfactual demonstration), `README.md`.
* `mmacos/tests/tPolContrast.m` —
  `test_overcoat_trade_reverses_across_the_quarter_wave_condition`
  (13 assertions, 3.7 s). Class setup gains the tool on the path, same
  arrangement as `tPolExternal`'s `vh_*` helpers.
* `mmacos/tools/pol_validation_report/run_pol_validation.m` — the C26 token
  block in the model-256 group, plus 8 new entries in that size's gate table.

**What did NOT move, deliberately:** `Rx_Cass_FarField.in`,
`tPolContrast/test_coating_sensitivity`'s 27.898 / 151.31, the engine, and
`tPolExternal`'s analytic reversal gate (which remains the independent check
this run is compared against).

---

## 7. Verification

| check | result |
|---|---|
| `tPolContrast` | **15/15 pass** (was 14) |
| the new gate | pass, 3.7 s |
| `oc_nonvacuity` against the achromatic counterfactual | **3/3 assertions fail** — non-vacuous |
| `make_polval.sh` (models 128 / 256 / 512, merge, render) | green — see §8 |

Not re-run this session, and stated as such: pymacos/ifx, PROPER-compare,
GMI. Nothing here touches the engine — the tool, the gate and the report are
pure binding-layer MATLAB and Markdown — so those cannot have moved.

---

## 8. For review, if tokens are short

1. **§3's excess-vs-total distinction.** The two ratios differ by a factor
   ~2.8 at the true quarter wave, entirely because the total carries the
   geometric floor. If that reading is wrong, the design-rule magnitude
   (\"buys ~1.7 decades\") is the number that moves — the *sign* of the
   reversal does not depend on it.
2. **§4's control.** It is the whole non-vacuity argument, and it is a
   counterfactual reached by asking the real engine for a different film
   rather than by modifying the engine. Worth checking that this is a fair
   stand-in for an achromatic implementation and not a weaker claim.
3. Whether §5.5's residual prose should shrink further now that §8.3.1 exists
   — the correction notice is long, and two sections now carry the same
   physics.
