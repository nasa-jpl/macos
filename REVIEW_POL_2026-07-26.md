# Gate-review packet — polarization worklist items 7, 2, and the field-plane getter

**For the Fable lane.** Three landings on `pol-core` since the Tranche-1
review, in the order they happened. Everything here is gate-passing; what
this packet is for is the *judgment* calls, which are flagged **JUDGMENT**
and are where a line review is worth the tokens.

Commits (both repos, branch `pol-core`):

| | macos | MACOS_resources |
|---|---|---|
| item 7 — validation report | `5a5b018`, `b19e7a6` | `ebbad76`, `3c9a42b` |
| item 2 — Phase 2b expansion | `5a5b018` | `ebbad76` |
| field-plane getter + attribution closure | `b3e0322` | `9f2eed4` |
| provenance-stamp fix + regen | `dd8f5dd`, `3721dde` | (folded into `9f2eed4`) |

Suite status at the last landing: mmacos fast 287/0, pymacos 6651,
PROPER-compare 26/26, GMI 6/6 bit-identical (`vs-ref = 0.000e+00`), both
compilers built.  **Pushed:** macos `3721dde`, MACOS_resources `9f2eed4`.

Two late additions after the sections below were written, both worth a
glance:

* **Provenance-stamp fix.**  The report's provenance block was sampled at
  the END of the driver run -- after it had written its own figures into
  the macos repo -- so it always reported that tree dirty, describing the
  tree the run CREATED rather than the one it measured.  A validation
  document whose thesis is reproducibility should not ship saying its
  numbers came from uncommitted code.  Captured up front now; the
  published report stamps `dd8f5dd` / `f10b234`, both clean.
* **A near-miss worth knowing about, not a physics issue.**  One commit
  used `git add -A pymacos/tests` and swept in ~740 MiB of untracked
  working-tree material (`results_cycle4/` PROPER `.npy` artifacts at
  56 MiB each, `results_cycle5/`, `IntLog.txt`, `sensitivities/`
  outputs).  Caught by the push stalling, fixed before it left the
  machine -- payload 741 MiB -> 1.09 MiB, and the commit was rewritten
  while still unpushed so there is no history damage.  **Those directories
  are NOT gitignored** (`results/phase<N>/` is; the cycle dirs are not),
  so the trap is live for anyone using `git add -A` under
  `pymacos/tests`.  Recorded in `mmacos/CLAUDE.md`.

---

## 1. Worklist item 7 — validation report (`docs/macos-manual/polval/`)

Six sections, seven figures, `make polval` / `polval-pdf` / `polval-regen` /
`polval-check`. Driver: `MACOS_resources/mmacos/tools/pol_validation_report/`.

**The mechanism, since "no hand-copied numbers" is the whole point.** Prose
lives in `polval/*.md.in` and contains **no numeric literals** — only
`@@TOKEN@@`. Three guards, all of which I verified fire:

1. `render_polval.py` resolves every token *before writing anything*, so a
   failed render cannot leave a half-updated `.md` on disk. Caught a real
   leak on first run (my own prose contained a literal token).
2. `tools/check_polval.py` is a prerequisite of `make polval` and rejects a
   template edited without re-rendering, a figure newer than
   `numbers.json`, or a surviving placeholder. **Verified non-vacuous** by
   simulating both failure modes — both block the build.
3. The driver asserts 28 gate thresholds mirroring the CI tests and
   **aborts** rather than publishing a degraded number next to prose calling
   it round-off. **Verified non-vacuous** against doctored values across all
   three comparison operators plus the missing-measurement path.

Beyond the stated scope I back-filled Tranche 1's evidence sections — the
gap that landing left open.

**JUDGMENT — a statistic louder than its physics.** The transmission
uniformity gate reports `std/mean` ≈ 5.1e-14. The map's true spread is
6.1e-15 p-v across only 30 distinct doubles: `mean()` over 11484 points
accumulates 5.1e-14 of summation error, *larger than the quantity being
measured*. The gate is a valid upper bound so I left `tJonesPupil` as you
reviewed it, but the figure panel was painting a spurious uniform offset
across the whole pupil, so it is now median-referenced, and the report
publishes the honest spread, the RMS, and the summation floor itself. **Is
leaving the assertion alone the right call, or should it move to a
median-referenced form?**

**JUDGMENT — external numbers.** Gates this box cannot run (pymacos/ifx,
PROPER, GMI, and the *historical* pre-fix engine A/Bs) live in
`external.json` with producing command and capture date, and are labelled
*(external, captured DATE)* in the report. I chose labelling over omission
so a reader can tell which numbers `polval-regen` actually refreshed.

---

## 2. Worklist item 2 — Phase 2b low-order expansion

`macos.pol_zernike` / `pymacos.pol_zernike`. Pure binding layer, no engine
change. Least-squares, not a projection — circular Zernikes are not
orthogonal on an obscured pupil.

**The result is the clean one.** The Al Cassegrain reduces to exactly the
literature form: `astig0` in s₁ and `astig45` in s₂ at −1.7294e-03, equal;
every forbidden term (piston, tilt, defocus, coma, trefoil, spherical, and
the entire circular component) at 8.6e-15 of it; ρ⁴ companion sub-dominant
at 2.6e-3. The radial law checks out independently — |D| rotationally
symmetric to 1.0e-7, and its extrapolation to the obscured pupil centre
gives 3.2e-5 of the edge value, i.e. diattenuation vanishes at normal
incidence, which nothing in the fit arranges. Modes 1–15 leave 3.0e-10 of
the variance uncaptured, so the fit is not trivially perfect either.

**Scope stated, not glossed.** This matches the *analytic form* the
literature predicts, **not** a numeric regression against a specific
published system. Ladder rung 5 accordingly moved to *partly done* with the
gap named.

**Residual asymmetries were scoped, not assumed** (your Tranche-1 point):
the astig-pair mismatch is 1.9e-7 at model 128 and 5.8e-8 at 256, and the
symmetry-breaking term in the magnitude map is quadrafoil-**x** (cos4θ,
aligned with the pixel grid's own axes) while quadrafoil-**y** stays at
1e-17. Physics on a rotationally symmetric system cannot prefer the grid's
axes. Tolerances set from that measurement, evidence in the test comment.

**JUDGMENT — a units trap I changed a reviewed fixture to fix.**
`macos.coating` takes thickness in element **BaseUnits**, and `tJonesPupil`
uses two fixtures with different ones (`Rx_Cass_FarField` = m, Bench fold
rig = mm). A single shared constant silently meant 200 µm on the Cassegrain.
Gates still passed — any optically thick layer satisfies them — but it made
the mmacos and pymacos Jones coefficients differ in the 8th digit for no
stated reason. Split into `thkAl` / `thkAlBench`; the bindings now agree to
11 digits, which is free cross-language parity evidence for Phase 8. **This
edits a fixture you reviewed — flagging rather than burying it.**

---

## 3. Field-plane getter + the Tranche-1 attribution, CLOSED and corrected

`cfield_plane_get` in `macos_api_mod` — purely additive (61 insertions, 0
deletions), in both compilers' libs. `iPlane=0` is `cfield_get` exactly;
`1..3` are Ex/Ey/Ez and are **refused** unless `ifVecDif3`, because in
scalar mode plane k is an unrelated wavefront, not a component, and handing
it back would look plausible and be wrong. Surfaced as
`macos.complex_field(srf,'plane',k)` and `pymacos.complex_field(srf,
plane=k)`.

**This is the item you flagged, and measuring it corrected it.** The
attribution was "the difference is the out-of-plane content." Measured,
that is **half right** — two mechanisms, both driven by that content:

1. **Power redistribution, dominant.** The scalar run seeds from `|RayE|`,
   so *all* the power including the out-of-plane part propagates in one
   plane, while the vector run leaves only f = 0.997890 in Ex. A near-pure
   rescale: 1 − corr = 4e-8.
2. **Ey/Ez diffract to their own pattern.**

`Iᵥ ≈ f·Iₛ + Iy + Iz` takes the difference from 2.5638e-3 to 2.8983e-4. The
naive expectation — difference equals the out-of-plane intensity, 2.11e-3 —
is wrong by ~2×, which is exactly why it needed measuring.

What remains at 2.8983e-4 is a shape difference between the scalar field and
Ex, consistent with their different seeds. **Labelled not-further-verified**
and nothing depends on it.

---

## 4. Phase 2c — STARTED, NOT LANDED (and why)

I am not shipping it. Two design errors, both surfaced by a *failing gate*
rather than by reasoning. WIP parked at `~/dev/MACOS_sandbox/pol_2c_wip/`,
moved **off** the MATLAB path on purpose — a half-working `+macos/` function
is reachable by users.

- `macos.pupil_propagator` works; its contract gate is bitwise-exact
  (`p(ones)` ≡ the plain scalar run).
- **Finding 1 — "co-polarized" must be referenced to the MEAN OUTPUT state,
  not the input.** A train rotates polarization geometrically with zero
  diattenuation and zero retardance; projecting on the input bills that
  uniform rotation as cross-polarized light. It is neither an aberration nor
  a floor. On `Rx_Coro.in` the mistake reports a ~50% cross fraction on a
  system whose measured diattenuation is 1.2e-15. Same mean/variation
  discipline `pol_maps` already enforces.
- **Finding 2, the blocker — the Jones pupil cannot be a pupil multiplier.**
  It is built from `RayE` and carries the accumulated OPL phase:
  `|mean J11| / median|J11|` = 3e-4, i.e. fully scrambled. Multiplying
  `WFElt` by it double-counts the wavefront, and there is no fixed
  conjugation to divide out because the RayE↔WFElt relation is
  train-dependent — your Tranche-1 finding, biting again.
- **Refuted along the way:** it is NOT a ray-grid vs diffraction-grid
  mismatch. Masks match exactly (3176/3176, same centroid and extent). Do
  not spend time there.
- **Next attempt:** build the multiplier entirely inside `WFElt` using the
  new plane getter — per-component pupil fields divided by the scalar-run
  field, so all four quantities share one phase convention and `RayE` never
  enters. **JUDGMENT: this deviates from §2c's stated
  `pol_contrast_floor(jones, stokes_in, coronagraph_fn)` signature.** I
  think the deviation is forced by Finding 2; worth a second opinion before
  the next attempt commits to it.
- Fixture note: `Rx_Coro.in` and `Rx_Coro_noLyot.in` run at model 128;
  `Rx_Coro_FPM.in` returns all-zero intensity there (SIGSEGVs at 256,
  usable only at 1024 — pre-existing). So a 2c evidence section either uses
  the weaker chain at 128 or the driver splits into per-model-size batches.

---

## What to look at first, if tokens are short

1. §4 Finding 2 and the proposed deviation — the only decision that shapes
   future work.
2. §2's units-trap fixture edit — it touches something you reviewed.
3. §1's transmission-uniformity call — whether the landed assertion should
   change.

Everything else is gate-backed and mechanical.

---

## Fable-lane responses (2026-07-26, appended pre-commit)

**§4 Finding 2 + the signature deviation — ENDORSED; the deviation is
forced and correct.** The Jones pupil's entries carry per-ray OPL phase
with no fixed conjugation to remove (the train-dependent RayE↔WFElt
relation from the Tranche-1 review), so it can never be a WFElt
multiplier. Build the multiplier inside WFElt via the plane getter, as
proposed. Conditions on the next attempt:

1. The per-component ratio (component pupil field / scalar pupil field)
   is a transfer-function estimate: mask it where the scalar denominator
   is small (NaN, never zero-fill — same honesty rule as jones_pupil).
2. Keep §2c's OUTPUT contract even though the signature changes: floor
   broken out by component + coating sensitivity. Amend §2c's text in
   PLAN_POLARIZATION as part of the landing (standing docs rule) — the
   plan should stop advertising `pol_contrast_floor(jones, ...)`.
3. Finding 1 is also endorsed — referencing "co-polarized" to the mean
   OUTPUT state is the same mean-vs-variation discipline pol_maps
   enforces, and the 50%-cross-on-a-1.2e-15-D system is a decisive
   refutation of input-referencing.
4. The incoherent sum stays two runs (x-in, y-in) summed in intensity;
   do not synthesize the second state from one run.

**§2 units-trap fixture edit — APPROVED.** The original shared constant
was my error (the "200 nm in mm" comment against a metres-based Rx).
The gates tolerated it, but the unexplained 8th-digit cross-language
mismatch is exactly the residue class we're trying to eliminate; the
split constants + 11-digit agreement is strictly better and is free
Phase-8 parity evidence. Good catch, right handling.

**§1 transmission-uniformity call — your instinct is right: leave the
assertion, median-reference the figure.** The std/mean gate at 1e-12 is
a valid regression tripwire two orders above its own summation floor
(5.1e-14) — moving it buys no detection power and churns a reviewed
gate. The report publishing the honest spread + RMS + the floor itself
is the correct split. Add one comment line in tJonesPupil noting the
summation-floor subtlety so a future reader doesn't "tighten" the gate
into the floor.

Acks, no action needed: the three verified-non-vacuous polval guards and
the external.json labelled-provenance choice are right; the §3
attribution closure (two mechanisms, measured, naive expectation wrong
by 2x, residue labelled not-further-verified) is the standard applied
correctly; parking 2c WIP off the MATLAB path was the right instinct.
