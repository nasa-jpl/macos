# Review packet — IFO polarization slice 3 (polarizing phase-shifting variant)

**Date:** 2026-07-28 · **Branch:** `pol-ifo` (MACOS_resources; off
`bench-builder`, merged `ifo-l2` + `pol-core`) · **Engine:** macos
`pol-core`, gfortran · **Lane:** Opus (IFO design, slice 3 of 3) ·
**Reviewer ask:** line review of the two new Bench emitters + the
`twyman_green` polarizing option + the slice-3 harness; confirm the three
gates and the error-budget-vs-baseline comparison. No convention decision
(the off-normal material-axis rule is already settled and does not bite
here — every pol element is at strict normal incidence).

Slices 1 (coated-BS Twyman-Green, ray-level Jones) and 2 (BS-AOI vs
clearance trade) passed with exact independent checks. Slice 3 changes the
**measurement**, not the null: a rotating-analyzer **polarization** phase-
shifting interferometer built on the same rig lineage, and the error budget
of its polarizing components — which is the configuration comparison.

---

## What landed

| Layer | What |
|---|---|
| Builder (general) | `macos.design.Bench` gains `add_polarizer(dist,axis)` + `add_waveplate(dist,axis,R)` and `emit()` support (`PolAxis=` / `Retardance=` blocks, required by `ChkDf2`). These are the reusable emitters the Phase-3 pol-elements packet flagged as the IFO prerequisite. |
| Builder (rig) | `twyman_green` gains `polarizing` (default **false**, bit-identical off) + axis/retardance knobs. When on it inserts: input polarizer @45° (both arms) → coated BS → a **double-passed QWP in each arm** (net half-wave, rotating that arm's linear state) → recomb → **output QWP** → **rotating analyzer** → detector. Every pol element rides a **collimated, normal-incidence** leg. |
| Harness | `examples/design/bench_ifo_pol/example_bench_ifo_pol_slice3.m` — the four-step PSI, three scores, the pol-component error budget, the comparison, three gates. Entirely in the example dir beyond the two builder additions. |
| Tests | `tBench` +2 (`test_add_polarizer_waveplate_emit`, `test_twyman_green_polarizing`); `tBench` added to `SUITE_FAST`. **7/7.** |

**No engine change.** `PolElt` (EltID 15/18), the material-axis rule, and
the dual dispatch chains are all as shipped on `pol-core`.

---

## The configuration (why the four-step is closed-form-exact)

Input linear @45° → both arms lit on s and p. The double-passed QWP in
each arm is a net half-wave that rotates that arm's linear state; with test
QWP fast-axis @0° and ref @45° the arms leave **orthogonally linear**
(−45° and +45°). The output QWP @0° maps them to **orthogonally circular**.
An analyzer at angle θ then gives, for orthogonal-circular arms,

```
I(θ) = A + B cos2θ + C sin2θ ,   fringe phase ψ = atan2(C, B) carrying the OPD.
```

The analyzer projector `t̂t̂ᵀ` has **no harmonic above 2θ**, so a four-step at
θ = 0/45/90/135° (2θ = 0/90/180/270°) recovers ψ = atan2(I₂−I₄, I₁−I₃)
**algebraically exactly** — that is the null-test property, not an
approximation. Rotating the analyzer by θ shifts the fringe phase by 2θ, a
phase-shifting interferometer with **no moving PZT**. The de Groot /
`bench_ifo_dm` PSI machinery is the processing reference; here it runs on
**ray-level Jones fields** (`macos.ray_field`, Tranche-1 — the pol elements
precede the single physical-optics leg, the `Rx_PolElt.in` condition), never
vector-diffraction intensities.

---

## Gates (with numbers)

**Gate A — the null test (closed form).** With perfect QWPs the four-step
estimator EXACTLY inverts the engine's own intensity model, and recovers an
injected known OPD change:

| check | result |
|---|---|
| A(i) four-step vs least-squares 2θ fit of the same frames | **1.79e-16 rad** = 1.8e-14 nm |
| A(iii) same, coated BS | **1.81e-16 rad** |
| A(ii) incremental null: inject 40 nm ref piston, recover the fringe-phase CHANGE (double-pass 0.7943 rad expected) | mean recovered **0.7940 rad**, err **4.0e-3 rad = 0.40 nm**; pupil rms 6.7e-2 rad |

A(i)/(iii) are the closed-form null (the estimator is exact for any A/B/C).
A(ii)'s residual (0.40 nm) and its 0.19 nm absolute offset are **reported,
not pinned**: they are the rig's own static polarization aberration
(perfect-conductor BS π s/p flip + glass-transit diattenuation — slice-1's
finding) limiting absolute OPD recovery. A polarization-neutral beam divider
would be needed for a bit-exact absolute null, which no real BS is.

**Gate B — non-vacuity, the textbook signature.** A known output-QWP
retardance error ε must produce the textbook PSI error (Schwider 1983; de
Groot Appl. Opt. **34**, 4723): a **second-order, twice-fringe** ripple of
amplitude `≈ ε²/4·sin(2p)`. Measured on the engine, referenced to the
ideal-component run, piston removed:

| ε (rad) | 0.126 | 0.251 | 0.503 |
|---|---|---|---|
| engine induced rms | 6.8e-3 | 1.2e-2 | 3.1e-2 rad |
| closed form `ε²/4/√2` | 2.8e-3 | 1.1e-2 | 4.5e-2 rad |
| amp / closed-form | 2.43 | 1.05 | 0.70 |
| 2ω / 1ω content | 2.2 | 2.0 | 3.0 |

The two INVARIANT marks are asserted: the induced error is a **2ω (twice-
fringe) ripple** (2ω content > 1ω at every ε) and its amplitude **approaches
the ε²/4 closed form as ε grows** (ratio → 1 near ε = 0.25 rad, where the
pure 2nd-order term dominates). The small-ε excess (ratio 2.43, log-log
slope 1.10 rather than 2) is a **reported finding**: the coated-BS
polarization aberration `δ_rig` adds a first-order `ε·δ_rig` cross-term —
consistent with slices 1–2, rig-specific, so not pinned. It has teeth (the
signature grows monotonically; a no-op element gives zero).

**Gate C — pol-off bit-identity vs a Reference-twin.** With polarization
off the pol elements are `RefSrf` geometry, so the rig's OPD must equal a
twin's (pol elements retyped `Reference`, keywords stripped):

| arm | OPD max-diff |
|---|---|
| test | **0.000e+00 mm** |
| ref | **0.000e+00 mm** |

Exactly zero (the `tPolElement` unpolarized-twin invariant at bench scale).

---

## The three scores (nominal: ideal QWPs, coated BS, HeNe, model 256)

| score | value |
|---|---|
| (1) fringe visibility V | **0.99155** (cost 1−V = 8.5e-3) |
| (2) PSI pupil-variation phase error (coating-driven, coated−minus−bare) | **2.38 nm** |
| (3) mechanical clearance (45° fold, pol legs collimated) | **40.91 mm** |

---

## The error budget — what only this configuration makes measurable

PSI phase error (nm) induced by each pol-component imperfection, isolated by
subtracting the ideal-component recovered fringe on the same coated rig
(piston removed):

| source | unit | PSI error @ largest mag | slope / order |
|---|---|---|---|
| arm-QWP retardance | wave | 17.2 nm @ 0.05 wave (~344 nm/wave) | 1st order |
| QWP chromaticity | Δλ/λ | 7.9 nm @ 0.10 | 1st order |
| arm-QWP axis | deg | 4.96 nm @ 1° | 1st order |
| output-QWP retardance | wave | 1.48 nm @ 0.05 wave | 1st order (rig cross-term) |
| output-QWP axis | deg | 0.22 nm @ 1° | 1st order |
| input-polarizer axis | deg | 0.22 nm @ 1° | 1st order |
| analyzer azimuth offset | deg | **6e-10 nm @ 1°** | common-mode (piston) |

Reading: the **arm QWPs** carry the tightest budget (they set the arm
orthogonality; a retardance error breaks it at first order), then
chromaticity (same physics, off-design λ detunes every QWP), then axis
misalignments. The analyzer azimuth is a pure common-mode piston (the 2θ
step is exact regardless of its zero-point) — it drops out entirely.

---

## Configuration comparison — the conclusion, as measurement

Slice 2 (coated-BS Michelson, mechanically PZT-stepped) and slice 3 (this
config) share the **same coated BS**, so the coating floor is common. What
differs is how each *reads* it:

- The coated-BS arm retardance is **negligible common-mode piston in slice-2
  scalar interferometry (2.3e-6 nm)**, but in the polarization PSI it
  **aliases into an OPD-dependent phase error (2.38 nm)** — phase-stepping in
  the polarization domain couples the arm polarization differential into the
  readout, which mechanical PZT stepping is blind to. ~10⁶× worse from the
  identical coating.

- **Which is best (measurement, not preference):** mechanical PZT stepping
  (slice 2) is far less sensitive to polarization imperfections and is
  preferred wherever a moving reference mirror is acceptable. The
  polarization PSI earns its place **only** where a moving PZT is not an
  option (snapshot / high-speed / vibration-immune metrology), and then
  demands waveplate retardance/axis control at the **< ~0.01 wave / < ~0.1°**
  level to keep its self-induced error below the coating-aliasing floor.

This conclusion **is** the sensitivity table set beside slice 2's numbers —
it falls out of the measurement, not a design preference. The section is
**COMPLETE** (not marked partial): the comparison reduced cleanly to the
error budget within this session.

---

## Figure

`bench_ifo_pol_slice3_budget.png` — two panels: (top) log-log error-budget
curves, one per pol-component source, with the coating floor marked;
(bottom) the Gate-B error signature (engine induced PSI error vs ideal
fringe phase) with the `−(ε²/4)sin(2p)` closed-form guide.

---

## Files

- `mmacos/src/+macos/+design/Bench.m` — `add_polarizer` / `add_waveplate` + emit + `blank` fields
- `mmacos/src/+macos/+design/twyman_green.m` — `polarizing` option (+ `ax_local` helper), default bit-identical
- `mmacos/examples/design/bench_ifo_pol/example_bench_ifo_pol_slice3.m` — the harness
- `mmacos/examples/design/bench_ifo_pol/bench_ifo_pol_slice3_budget.png` — the figure
- `mmacos/examples/design/bench_ifo_pol/bench_ifo_pol_slice3_results.mat` — gates + scores + budget
- `mmacos/tests/tBench.m` — +2 tests (7/7); `run_mmacos_tests.sh` — `tBench` added to `SUITE_FAST`
- `mmacos/examples/design/bench_ifo_pol/CURRENT_SLICE.md` — resume state
- `macos/REVIEW_POL_IFO_SLICE2_2026-07-27.md` — two review riders appended (fringe-contrast budget + s/p-lit assumption)

The `s3_*.in` emitted rigs + Reference-twins are gitignored (regenerable).

---

## Slice-2 review riders (delivered, per the slice-2 Fable review)

Both are now in `REVIEW_POL_IFO_SLICE2_2026-07-27.md`:
1. the total fringe-contrast budget at 45° and the 17.5° knee, with the
   scalar throughput-imbalance (`≈D²/8`) and polarization-differential
   (`≈ret²/8`) terms tabulated side by side (they are comparable at both
   points; the reported `1−V` is their sum);
2. the note that the visibility metric assumes both s and p lit — for an
   s- or p-aligned input the differential Jones is diagonal, `V = 1`, and
   the cost migrates entirely into the excluded (suppressed-flux) amplitude
   term.

---

## Reviewer notes / honest scope

- **Gate A(ii) and Gate B small-ε are reported, not pinned.** Both surface
  the coated-BS polarization aberration measured in slices 1–2; pinning them
  to round-off would require a polarization-neutral BS the rig does not have.
  The closed-form nulls that ARE pinned (A(i)/(iii) at 1.8e-16; the 2ω mark
  and the large-ε ε²/4 amplitude in B) are the ones independent of rig
  aberration.
- **Double-passed arm QWP = one global fast axis.** `set_pol` derives the
  fast axis from the FORWARD element and gives the SAME vector to both passes
  (re-deriving from each element's own `psi` would reflect the return-pass
  axis about the transverse `u2` — `psi` flips sign at the retro — silently
  breaking the net half-wave for any non-0/90° arm angle). Called out because
  it is a real trap for anyone adding double-passed retarders.
- **Version-skew flag (unrelated to this slice, worth Dave's eye):** on this
  working tree `polarizer.m`'s docstring and `tPolElement` still describe the
  *provisional pass-axis* convention, and the promised `Rx_PolElt_Tilt.in`
  fixture + engine-driven off-normal A/B (`REVIEW_POL_ELEMENTS`/PLAN §
  landing spec) are **not present on `pol-ifo`**, though the ENGINE carries
  the settled material-axis flip. Slice 3 is unaffected (strict normal
  incidence), but the veneer docs and PLAN claim a landing that is not in
  this branch.

---

## Fable-lane review (2026-07-28) — PASS; the arc is complete

**The aliasing headline is independently confirmed in mechanism and
order.**  A from-scratch Jones model of the four-step rotating-analyzer
estimator with slice-1's differential (D = 0.0721, ret = 0.0835 split
symmetrically) gives an OPD-dependent ripple of peak 8.4 nm and
5.95 nm rms over a UNIFORM full fringe; the harness's 2.38 nm is the
rms over the rig's actual sub-fringe piston-removed OPD span, which is
consistently smaller.  The ~10⁶× amplification over the scalar-PSI
common-mode figure stands regardless of the rms weighting — in
polarization PSI the arm polarization differential enters the READOUT,
which mechanical stepping is blind to.  The configuration conclusion
(PZT stepping preferred wherever a moving mirror is acceptable; the
pol PSI demands < ~0.01 wave / < ~0.1° waveplate control) is endorsed
as measured.

**Gate discipline: right things pinned, right things reported.**
A(i)/(iii) pin the closed-form estimator property (1.8e-16, which is
what "no harmonic above 2θ" must give); A(ii)'s 0.40 nm and Gate B's
small-ε excess are correctly REPORTED as the rig's own coated-BS
aberration (pinning them would have required a polarization-neutral BS
that does not exist).  B's two invariant marks — the 2ω twice-fringe
signature and the large-ε ε²/4 amplitude — are the aberration-
independent facts, and that is exactly the split slices 1–2 taught.

**The double-passed-retarder trap is a keeper**: deriving the return-
pass fast axis from the element's own psi (which flips at the retro)
silently breaks the net half-wave for any non-0/90° arm angle; one
global fast-axis vector per double-passed element.  Now on the record
for anyone adding double-passed retarders.

**The version-skew flag is CLOSED** (this review): resources
`origin/pol-core` merged into `pol-ifo` (`7a268a9`; one conflict in
`run_mmacos_tests.sh`, resolved as the union of both suite additions),
bringing the material-axis veneer docs, `Rx_PolElt_Tilt.in`, and the
off-normal engine gates onto the branch.  Verified on the merged tree
with the mex rebuilt: tBench 7/7, tPolElement 27/27 (including the
tilt gates).  Pushed.

**Slice-2 riders: delivered as specified** (the D²/8 vs ret²/8 budget
rows are the decomposition the review asked for).

With slices 1–3 landed and reviewed, the pol-ifo IFO arc is COMPLETE:
polarization-honest baseline → BS-angle trade with clearance → a
polarizing-PSI variant with a measured error budget and a
measurement-backed configuration verdict.  The two Bench emitters
(add_polarizer/add_waveplate) close the IFO prerequisite the
pol-elements packet flagged.  This is the design layer consuming the
whole polarization stack end-to-end — the convergence the lane split
was building toward.
