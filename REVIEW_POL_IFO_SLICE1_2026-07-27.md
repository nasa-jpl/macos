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
