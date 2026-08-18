# Review packet — IFO polarization slice 1 (Twyman-Green, coated BS)

**Date:** 2026-07-27 · **Branch:** `pol-ifo` (MACOS_resources; off
`bench-builder`, merged `ifo-l2` + `pol-core`) · **Engine:** macos
`pol-core`, gfortran · **Lane:** Opus (IFO design, slice 1 of 3) ·
**Reviewer ask:** line review of the example + confirm the two findings
below; no convention decision.

---

## What landed

A polarization-honest Twyman-Green scoring harness, entirely in one
example directory — **no engine change, no shared design-layer change**
(the whole slice is `examples/design/bench_ifo_pol/`). It puts a real
Al coating on the beam-splitter face, turns `ifPol` on, and scores the
arm-differential diattenuation / retardance at recombination through the
**already-gated Phase-2 layer** (`jones_pupil` + `pol_maps`), ray-level
only.

| Layer | What |
|---|---|
| Rig | `macos.design.twyman_green` (unchanged), emitted for both arms |
| Coating | Al `n=1.45, κ=7.54`, 200 nm physical, applied via `macos.coating` **after** `load_rx` |
| Jones | `jones_pupil` (double-pole basis) at each arm's Recomb; ref arm forced into the test arm's exit basis (`'axis'`,`'xref'`) |
| Differential | 2×2 per-ray transfer `M = J_test · inv(J_ref)` → `pol_maps` → D / retardance |
| PSI error | pupil variation (piston-removed) of the co-pol fringe phase `arg⟨E_test,E_ref⟩` |
| Gates | (1) single-surface 45° Fresnel analytic on the BS reflect; (2) pol-off bit-identity |

### Scoring rule honored (Tranche 1, §3)

The BS sits **past the first propagation leg**, where Tranche 1 caps the
diffraction grid. Every number here comes from `macos.ray_field` via
`jones_pupil` — **no vector-diffraction intensities are used to score.**

---

## Gates (with numbers)

**Gate 1 — coating machinery vs textbook Fresnel.** Trace the test arm
with the BS coated, s+p lit; pull the per-ray reflected s/p amplitudes at
the BS and the **measured incident field** (a separate trace to the
element just before the BS — see finding 2); form the
convention-independent ratio `(Es/Ep)·(qp/qs)` and compare to the
Born&Wolf bare-interface coefficients (`r_p = (N₂cᵢ−N₁cₜ)/(N₂cᵢ+N₁cₜ)`,
`r_s = (N₁cᵢ−N₂cₜ)/(N₁cᵢ+N₂cₜ)`, N₂ = n−iκ). Analytic written from the
textbook, **not** transcribed from the engine, so the r_p sign is pinned
non-circularly (this is exactly the gate-blindness lesson from
`REVIEW_POL_SP_SIGN`). Mirrors `tJonesPupil/test_fold_fresnel_analytic`.

| quantity | value |
|---|---|
| BS mean AOI | 45.000° |
| RS/RP magnitude residual | **1.13e-14** |
| RS/RP phase residual | **2.97e-14** |

PASS (threshold 1e-12).

**Gate 2 — pol-off bit-identity.** With polarization OFF the coating must
be inert, so the coated train's OPD equals the uncoated train's to
round-off. Measured on both arms, uncoated vs coated, pol OFF:

| arm | OPD max-diff |
|---|---|
| test | **0.000e+00 mm** |
| ref | **0.000e+00 mm** |

PASS (threshold 1e-12 mm). Exactly zero, which is the strong form — the
coating touches nothing on the scalar path.

---

## Result (the physics the harness produces)

Arm-differential polarization at recombination, x-pol input, 3001 common
rays:

| quantity | mean | pupil variation (rms) | max |
|---|---|---|---|
| Diattenuation D | 7.21e-02 | 2.0e-08 | 7.21e-02 |
| Retardance (rad) | 8.35e-02 | 3.6e-08 | 8.35e-02 |

**PSI phase-error contribution** (piston-removed co-pol fringe): **2.3e-6
nm RMS** at 632.8 nm.

Reading: the coated BS face is hit as an **external** air→Al reflection in
the test arm but an **internal** glass→Al reflection (plus two glass
transits) in the reference arm — genuinely different Jones — so the
arm-differential D≈0.072 / retardance≈0.084 rad is a real common-mode
asymmetry, not identity. The pupil **variation** is at round-off (2–4e-8)
because both arms recombine near-collimated and near-common-path: the
polarization aberration is essentially pure piston, so the
spatially-varying PSI error is negligible (2.3e-6 nm) in this compensated
null. That is the physically-expected story for a balanced TG; slice 2
(BS AOI vs clearance) and slice 3 (polarizing-PSI variant) are where the
variation is meant to grow.

---

## Two findings worth the reviewer's eye

**Finding 1 — the Bench builder stamps the mirror perfect-conductor idiom
on transmitting Refractors, which is opaque glass under `ifPol`.**
`macos.design.Bench` sets `Extinc=1e22` on **every** element. On a
`Reflector` that is the standard neutral perfect-conductor trick
(RS=−1, RP=+1, polarization-neutral). On a **transmitting `Refractor`**
(compensator faces, BS substrate transits) an extinction of 1e22 means
perfectly **absorbing** glass, so with `ifPol` on the Fresnel
transmission kills the beam — the field collapses by ~1e-22 at the first
compensator face and is zero everywhere downstream. Scalar tracing
ignores extinction entirely, so this is **invisible without
polarization** and every prior (scalar) Bench IFO result is unaffected.

*Fixed in `Bench.m` (Dave's call, 2026-07-27):* the four transmitting
Refractor sites in `add_bs_transmit` and `add_bs_reflect_return` now stamp
`Extinc=0` (transparent glass); the Reflector in `add_bs_reflect_return`
keeps `1e22`. `blank()` already defaulted `extinc=0`, so `add_lens` was
always fine — this was purely the Reflector idiom copied onto the BS
transmit/return faces. **Verified scalar-path bit-identical:** the scalar
`bench_ifo` example, regenerated with vs without this change, produces
identical output (the only field that moved, `R_hat` = the recovered
2e7 mm ≈ flat radius, drifts 4e-2 mm = 2e-9 relative between ANY two runs
— an ill-conditioned-fit artifact present with the change reverted too).
`tBench` 5/5 green. The polarization gates are unchanged from the
example-workaround version (bit-identical numbers).

**Finding 2 — Gate 1 must use the field incident *on the BS*, not the
source launch state.** My first Gate-1 pass built the reference input
from the source launch frame and saw a 1.3e-3 **magnitude** residual with
a **round-off phase** residual (3e-14). That split is diagnostic: the
field arriving at the BS has already passed L1's two glass refractions,
which carry a ~1e-3 s/p **diattenuation** and (real index) **zero
retardance** — exactly a magnitude-only error. Using the measured
incident field (`ray_field(iBS−1)` from a separate trace; RayE is
per-element overwritten) as the reference makes L1 cancel in the ratio
and drops the residual to 1.1e-14. Noted because it is a general trap for
any single-surface Fresnel gate placed downstream of other optics.

---

## Files

- `mmacos/examples/design/bench_ifo_pol/example_bench_ifo_pol.m` — the harness
- `mmacos/examples/design/bench_ifo_pol/CURRENT_SLICE.md` — resume state
- `mmacos/examples/design/bench_ifo_pol/ifo_{test,ref}_pol.in` — emitted rigs
- `mmacos/examples/design/bench_ifo_pol/bench_ifo_pol_results.mat` — gate numbers + Jones pupils

No files outside that directory changed. `Bench.m` and `twyman_green.m`
are at their `bench-builder`/`ifo-l2` baseline (verified `git diff`
clean).

---

## Deferred (post-review)

- **Slice 2** — BS AOI vs mechanical-clearance trade (the AOI is where the
  arm-differential D/retardance and their pupil variation grow).
- **Slice 3** — a polarizing-PSI variant (ideal polarizer/waveplate in the
  collimated normal-incidence legs) + comparison against this baseline.

(The Bench `Extinc` policy from finding 1 is now settled and fixed in
`Bench.m`, not deferred.)

---

## Fable-lane review (2026-07-28) — PASS, with an exact independent check

**The headline numbers are independently CONFIRMED to all printed
digits.**  A from-scratch textbook Fresnel computation of the actual
emitted trains (test: external air→Al reflect × three 45° transit
pairs; ref: two transit pairs × internal glass→Al reflect; n_g = 1.5,
N_Al = 1.45−7.54i) gives D_diff = 0.0721, ret_diff = 0.0835 rad —
matching the harness's 7.21e-02 / 8.35e-02 exactly.  (A naive
one-reflection-per-arm model gives D = 0.103 and the same retardance;
the D difference is the net extra transit pair on the test side —
worth knowing when sanity-checking future variants.)  This is the
engine + jones_pupil + pol_maps + harness composing to the closed-form
prediction on a seven-surface differential train: the strongest
verification this slice could have.

**Finding 1 (Bench Extinc) — the fix is right and complete.**
1. `Extinc=0` IS the correct transparent-glass default.  Optical glass
   at visible wavelengths has κ ~ 1e-8..1e-7 (a negligible intensity
   effect over bench path lengths) and a REAL index contributes exactly
   zero retardance; stamping a fake small κ would inject arbitrary
   unphysical numbers into every Jones result with no measurement
   behind them.  Absorbing glass, if ever needed, is a per-design
   property, not a builder default.  The fix also makes Bench
   self-consistent with `blank()`/`add_lens`, which always used 0.
2. Sweep performed: every remaining `extinc = 1e22` stamp in the design
   layer sits on a Reflector (`add_bs_reflect`, `add_mirror`, the
   internal-reflect face of `add_bs_reflect_return`, `add_oap`,
   `add_relay`) — the conductor idiom where it belongs.  No other
   method carries the leaked idiom.  The scalar-invisibility point is
   worth its place in the record: scalar tracing ignores Extinc
   entirely, so this class of error is detectable ONLY under ifPol —
   one more reason the polarization-honest baseline was worth building
   before the trades.

**Finding 2 — endorsed, and the diagnostic is the keeper.**  The
magnitude-residual-with-round-off-phase split (1.3e-3 vs 3e-14) is
precisely the signature of an upstream REAL-index transit
(diattenuation without retardance), and reading it that way is what
turned a "gate tolerance problem" into "use the measured incident
field."  General rule, now on the record: a single-surface Fresnel
gate downstream of ANY optic must reference the measured incident
field, not the launch state.

**The interpretation ("correct for a collimated common-path null") —
endorsed.**  The means are real (external vs internal reflection +
unequal transit counts are genuinely different Jones) but common-mode:
they cost fringe VISIBILITY (~(1−cos 0.0835)/2 ≈ 0.2% modulation, plus
an amplitude-imbalance term ~D²/2 ≈ 0.3%) and shift the null piston —
neither is a PSI phase error.  The variation IS the PSI error term and
sits at round-off because a collimated planar bench gives every ray
the same 45.000° and the same plane of incidence.  2.3e-6 nm RMS is
the right answer for this configuration.

**One instruction for slice 2 (BS AOI vs clearance):** score BOTH
terms — the pupil-variation PSI error (which will stay small while
collimated) AND the visibility cost of the arm-differential MEANS
(which is what the AOI sweep actually moves at this bench's
f-numbers).  A trade scored on variation alone would conclude "AOI
doesn't matter," which is true only of the term that was already
negligible.

Slices 2 and 3 may proceed.
