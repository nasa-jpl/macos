# Review packet — IFO polarization slice 2 (BS-AOI vs mechanical-clearance trade)

**Date:** 2026-07-27 · **Branch:** `pol-ifo` (MACOS_resources; off
`bench-builder`, merged `ifo-l2` + `pol-core`) · **Engine:** macos
`pol-core`, gfortran · **Lane:** Opus (IFO design, slice 2 of 3) ·
**Reviewer ask:** line review of the builder parametrization + the
slice-2 harness; confirm the trade curve and its knee; confirm the
curve gate is non-vacuous. No convention decision.

Slice 1 (Twyman-Green, coated BS, ray-level Jones) passed with an exact
independent check and is the anchor for this slice. Slice 2 is the
BS-AOI vs mechanical-clearance trade: sweep the fold angle down from
45°, score what smaller AOI buys against what it costs, and report the
curve and its knee — not a prior.

---

## What landed

Two files change, both regressions-clean:

| Layer | What |
|---|---|
| Builder | `twyman_green.m` gains one parameter, `BS_AOI` (deg, default 45). The fold turn is `180 − 2·BS_AOI`; the reflect direction is `[cos turn; −sin turn; 0]`. **45° is pinned to the exact `[0;-1;0]` literal** so the default rig emits BIT-IDENTICALLY to the pre-slice-2 baseline (verified: both arms `diff`-clean). Every downstream face — compensator, BS-substrate transits, the internal-reflect return — already tracks the BS normal through the shared `bs` token and the recomputed chief, so only the reflect turn is set here. |
| Harness | `examples/design/bench_ifo_pol/example_bench_ifo_pol_slice2.m` — the sweep, the curve gate, the three scores, the trade curve, the non-vacuity check. Entirely in the example dir; no shared-layer or engine change beyond the one builder parameter. |

The builder change is the only edit outside the example directory. It is
additive and default-preserving; `tBench` 5/5 green.

---

## The sweep is gated, not just run

**Per-AOI, before any score:**
1. **Re-emit both arms** at the BS angle (compensator + substrate faces
   track automatically — verified by the Gate-1 measured AOI below).
2. **Trace-clean check** — both arms must keep ≥90% of the 45° reference
   successful-ray count (they keep 100% across the whole 45→15° sweep;
   the geometry stays feasible — it is *clearance*, not ray failure,
   that flags the crowding).
3. **Gate 1 (per angle)** — the coated single-surface Fresnel gate from
   slice 1, generalized to the measured AOI: engine s/p reflected
   amplitudes vs textbook Born & Wolf bare-interface coefficients at the
   actual AOI. `|RS/RP|` residual ≤ 1.1e-14, phase ≤ 3e-14 at every
   angle. Confirms the coating machinery follows the re-emitted geometry.

**Curve gate (the slice-1 Gate-1 idea extended to the whole sweep).**
At each AOI the engine's arm-differential mean D and retardance must
match a **full-train textbook closed form**, written from Born & Wolf
(not transcribed from the engine, so the r_p sign is pinned
non-circularly). Construction, per the Fable review:

> external-reflect × (net transit pair) vs (transit pair) × internal-reflect

The planar fold makes the whole coated path DIAGONAL in one common s/p
basis, so the differential reduces to two complex ratios
`M_{s,p} = (t_{s,p}^pair)^Δtp · r_{s,p}^ext / r_{s,p}^int`, with:
- `r^ext` = external air→Al reflection at AOI θ (test arm),
- `r^int` = internal glass→Al reflection at the refracted angle θ_g (ref arm),
- one **net glass transit pair** (Δtp = 1): the arms' glass OPL balances
  (compensator), but each transmissive plate crosses two air/glass
  boundaries per thickness while the internal-reflect return crosses only
  one, leaving one net transit pair in the differential.

Result: **all 13 AOI points pass the curve gate at relative-D ≤ 3e-13
and absolute-retardance ≤ 8e-15** — the engine reproduces the full-train
closed form to round-off across the entire sweep.

**Non-vacuity (the trap the Fable review documents).** The single net
transit pair is load-bearing. It is itself a strong s/p diattenuator, so
**dropping it (Δtp = 0, the one-reflection model) misses the engine's
measured D by 358%** at 45° (D0 = 0.0158 vs 0.0721). That mis-modeled
point is asserted to FAIL the same curve gate — so the gate has teeth.
(Retardance is transit-pair-invariant, glass being real-index, so the
D channel is the one that catches the omission — the gate checks both.)

---

## The trade (the physics the harness produces)

x-pol input, 3001 common rays, model 256, HeNe, bare Al BS coating
(n=1.45, κ=7.54, 200 nm). Three scores per AOI:

**(1) Fringe visibility from the arm-differential means.** Reported as
visibility, not raw D/ret: `V = sqrt((1 + sqrt(1−D²)·cos ret)/2)`, the
worst-case (45°-to-axes input) polarization fringe contrast set by the
differential. Falls monotonically as AOI drops.

**(2) PSI pupil-variation phase error.** Co-pol fringe phase, piston
removed. **Stays at 6e-7 … 2e-6 nm across the whole sweep** — round-off,
as expected while both arms recombine collimated and common-path. The
polarization aberration here is essentially pure piston; confirming the
pupil variation stays negligible **is a result, not a failure** (slice 3,
which puts polarizing elements in, is where it is meant to grow).

**(3) Mechanical clearance (MIN_SEP, from `ray_hist`).** Minimum
beam-ENVELOPE separation between the test-arm excursion (compensator →
test optic → return, the nodes that swing back as the fold tightens) and
the incoming source→BS collimated beam. Same physicality standard as the
MET launcher work (a stated 10 mm reference floor, not a prior). Falls
monotonically once the fold exceeds ~40° turn.

| AOI (deg) | D | ret (rad) | vis. cost 1−V | clearance (mm) | curve gate |
|---|---|---|---|---|---|
| 45.0 | 0.07213 | 0.08352 | 1.52e-03 | 43.6 | pass |
| 40.0 | 0.05451 | 0.05965 | 8.16e-04 | 44.8 | pass |
| 35.0 | 0.04013 | 0.04185 | 4.20e-04 | 42.3 | pass |
| 30.0 | 0.02850 | 0.02855 | 2.03e-04 | 36.7 | pass |
| 25.0 | 0.01923 | 0.01863 | 8.96e-05 | 28.1 | pass |
| 22.5 | 0.01539 | 0.01470 | 5.66e-05 | 22.8 | pass |
| 20.0 | 0.01202 | 0.01134 | 3.42e-05 | 16.9 | pass |
| **17.5** | 0.00912 | 0.00851 | **1.94e-05** | **10.4** | pass |
| 15.0 | 0.00664 | 0.00614 | 1.02e-05 | 3.5 | pass |

(Full 13-point table in `bench_ifo_pol_slice2_results.mat` and the run log.)

**Knee: AOI ≈ 17.5°** at the 10 mm clearance floor, where the fringe-
visibility cost is **1.9e-5 — a ~78× improvement over the 45° value
(1.5e-3)**. Reading: the polarization payoff (visibility cost, and its
driver the differential retardance) falls monotonically as AOI drops;
the mechanical clearance falls with it as the larger fold swings the
compensator/test-optic excursion back toward the incoming beam. The knee
is the smallest AOI that still clears the physicality floor — below it
the folded arm crowds the front end (3.5 mm at 15°). The closest-approach
node is the compensator return face (`Comptxbu`) at every angle.

The knee AOI scales with the reference floor and the bench's collimated
beam radius (~22 mm here); it is a *worked* datum on this rig, not a
universal number. A tighter bench (smaller beam, more clearance headroom)
pushes the knee lower and buys more polarization payoff.

---

## Figure

`bench_ifo_pol_slice2_trade.png` — two stacked panels, AOI axis reversed
(sweep direction 45°→15°, left→right):
- top: fringe-visibility cost (log, left) and differential retardance
  (right), both falling — the polarization payoff of smaller AOI;
- bottom: min beam-envelope clearance with the 10 mm MIN_SEP reference
  and the 17.5° knee marked.

---

## Files

- `mmacos/src/+macos/+design/twyman_green.m` — `BS_AOI` parameter (default 45, bit-identical)
- `mmacos/examples/design/bench_ifo_pol/example_bench_ifo_pol_slice2.m` — the sweep harness
- `mmacos/examples/design/bench_ifo_pol/bench_ifo_pol_slice2_trade.png` — trade curve
- `mmacos/examples/design/bench_ifo_pol/bench_ifo_pol_slice2_results.mat` — all scores + gate residuals
- `mmacos/examples/design/bench_ifo_pol/s2_{test,ref}_aoi{45,15}.in` — the sweep endpoints (evidence)
- `mmacos/examples/design/bench_ifo_pol/CURRENT_SLICE.md` — resume state
- `.gitignore` — the numbered sweep intermediates + `fort.6`

No engine change. `Bench.m` / `twyman_green.m` are otherwise at their
`bench-builder`/`ifo-l2` baseline (the `BS_AOI` addition is the only
delta, verified default bit-identical). `tBench` 5/5.

---

## Deferred (slice 3)

Polarizing-PSI variant (ideal polarizer / waveplate in the collimated
normal-incidence legs) + comparison against this baseline. That slice is
where the PSI pupil-variation (score 2) is designed to grow above
round-off.

---

## Fable-lane review (2026-07-28) — PASS; the judgment call is RIGHT with two riders

**The trade curve is independently CONFIRMED at every table point.**  My
slice-1 full-train construction, swept over 45/30/17.5/15°, reproduces
BOTH columns to all printed digits: D = 0.0721/0.0285/0.0091/0.0066 and
1−V = 1.52e-3/2.03e-4/1.94e-5/1.02e-5 (normalized state-overlap,
both-components-lit input — the same construction the harness scores).
The curve gate's non-vacuity (one-reflection model off by 358%) matches
the 40%-in-D error class my slice-1 review documented, now measured as
a gate failure.  The knee's honest labelling ("a worked datum on this
rig, not a universal prior" — it scales with beam radius) is the right
epistemics for a design-rules number.

**The flagged judgment — excluding the scalar throughput-imbalance
term — is CORRECT for a polarization trade**, for the reason CCMac
gave plus one more: the scalar imbalance is trimmed by conventional
means (split-ratio choice, ND trim) and including it would bury the
polarization signal the sweep exists to measure.  Two riders, one
paragraph each, ride slice 3:

1. **Context row, not a swept metric:** report the TOTAL fringe-contrast
   budget once at 45° and once at the knee — scalar-throughput
   imbalance and polarization-differential terms side by side — so no
   reader mistakes the swept metric for the whole budget.
2. **Record the input-state dependence.**  The scored metric assumes
   both s and p lit (the harness's 45°-ish input).  For input ALIGNED
   with s or p the differential Jones is diagonal: the state overlap
   reads ~1 and the entire polarization cost migrates into the
   amplitude-imbalance term the metric excludes — the metric would read
   "no cost" while the cost is real.  One sentence in the report
   prevents that misreading; a worst-case-over-input-states variant is
   the future-proof form if slice 3 wants it.

Also noted with approval: the Python-prototype-before-MATLAB step
(construction verified against slice-1's numbers before any harness
code), the per-angle trace-clean + Fresnel gates rather than
endpoint-only, and the BS_AOI builder default pinned to the exact
legacy literal with bit-identity verified both arms.

Slice 3 (polarizing-PSI variant) may proceed, with the two riders and
the slice-1 review's placement rule (ideal polarizer/waveplate in
collimated normal-incidence legs) as standing constraints.
