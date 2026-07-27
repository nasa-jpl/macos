# Gate-review packet — Phase 2c, the polarization contrast floor

Opus lane, worklist item 3. One item, one session, per the discipline rules
in `PLAN_POLARIZATION.md`. Built to the §2c analyzer-at-detector spec as the
Fable lane endorsed it (`REVIEW_POL_SP_SIGN_2026-07-27.md`, "Left for Opus"
item 4); both pupil-multiplier sketches stayed retired.

**Commits:** MACOS_resources `848245d`, macos `b726125` (docs + plan) and
`3bdb531` (regenerated report). Session A (the r_p sign fix tail) was not
touched.

---

## What shipped

| Piece | Where |
|---|---|
| `macos.pol_contrast_floor(pupil, det, ...)` + `Session` method | `mmacos/src/+macos/` |
| `tPolContrast` — 14 tests, model 256 | `mmacos/tests/` |
| `tPolContrastCoro` — 6 tests, model 512 | `mmacos/tests/` |
| `SUITE_POL_512` group + `polfloor` shortcut | `run_mmacos_tests.sh` |
| polval driver split per model size + `merge_numbers.py` | `mmacos/tools/pol_validation_report/` |
| report §5 (`60_contrast_floor.md.in`) + 7 gate-index rows + 2 figures | `macos/docs/macos-manual/polval/` |
| cmdref NOTES, manual §5 subsection, `PLAN_POLARIZATION` §2c amended | `macos/docs/`, `macos/` |

**Suite:** `polfloor` 20/20; full `run_mmacos_tests.sh` results below.

---

## 1. The design, as built

The spec was followed without deviation, so the interesting content is what
implementing it measured. Briefly, for the reviewer's frame:

- The split is at the detector on `complex_field(det,'plane',1..3)`. No pupil
  multiplier exists anywhere in the code.
- The analyzer is the dominant eigenvector of the pupil coherency matrix,
  computed per input state, so "co-polarized" is referenced to the mean
  output state.
- `'unpolarized'` is two traces summed in intensity, each with its own
  analyzer. Nothing is synthesized.
- Ratio maps NaN-mask small denominators (gated with `floor_tol` driven high
  enough that the mask actually bites — otherwise that assertion is vacuous
  on a full FFT grid, where nothing is below 1e-12 of peak).

Two API decisions that were NOT in the spec, both forced:

1. **`.scope`** — the function measures how much of the ray-level
   cross-polarized fraction the diffraction grid actually carries, and warns
   when it is not all of it. See §3; without it the Rx_Coro number would read
   as a floor when it is a lower bound.
2. **Coating sweep sets must all cover the same elements**, and the sweep
   leaves the last set applied. `coat_set` takes ≥1 layer, so a coating can be
   overwritten but never cleared; sets over different elements would silently
   accumulate. Validated with an error rather than documented away.

---

## 2. Gates, and whether they are vacuous

| Spec gate | Result | Fixture |
|---|---|---|
| x-pol reduces to the scalar contrast curve at round-off | cross **exactly 0**; `I_co` vs scalar 1.33e-15 of peak; radial contrast curve 2.36e-15 over 180 bins | `Rx_VecChain` |
| Parseval on the split | 1.34e-24 of peak | `Rx_Cass_FarField` |
| energy bookkeeping | closure 1.80e-16; `I_total == I_co+I_cross+I_long` bitwise; unpolarized == x + y bitwise | `Rx_Cass_FarField` |
| floor reported by component | co / cross / long, `cross/co = 7.0612e-07`, plus dark-zone mean/median/peak per channel | `Rx_Cass_FarField`, `Rx_Coro` |

**Where the reduction gate is run, and why not on Rx_Coro.** "At round-off"
is only defensible on a train where polarization is a no-op by construction.
`Rx_VecChain` is that train (the tVecChain argument, unchanged). On Rx_Coro
the polarized total legitimately differs from the scalar run — the scalar
seed puts all the power in one plane — so no tolerance there would be
defensible, and the gate would have been a tolerance-fitting exercise.

**Non-vacuity, checked rather than argued:**

- The **circular input state** caught a real bug (§4). Without it the
  analyzer gate passes with a conjugated coherency matrix.
- The **full-carry check** on `Rx_Cass_FarField` is measured against a
  ray-level cross fraction of 7.06e-07 — there is real cross-polarized light
  to fail to carry, so `carried == 1` is a result, not a 0/0.
- The **NaN-mask gate** asserts the mask is non-empty before asserting its
  contents, and asserts the unmasked entries equal the raw ratio (nothing
  zero-filled).
- The **polval gate guard** was re-verified on the new tables: degrading
  `C22_PARSEVAL`'s limit and adding an unmeasured token both abort the run
  with the report unwritten (`C22_PARSEVAL = 1.34293e-24 (needs < 1e-30)`;
  `C99_MISSING: not measured`). Reverted afterwards.

---

## 3. The finding — Tranche 1 caps the floor, and inverts the coating trade

This is the item that most wants a second opinion.

`Rx_Coro` is a coaxial chain of seven reflectors with propagation legs
interleaved. Tranche 1 seeds the component planes from `RayE` at the FIRST
physical-optics leg and thereafter multiplies all three by a common scalar
phase, so six of the seven mirrors transform the rays but not the grid. That
is the documented Tranche-1 validity condition; what was not known is the
size of the effect on a real chain:

| | grid (what diffraction carries) | rays (the full train) | carried |
|---|---|---|---|
| uncoated | 4.78973e-09 | 5.69403e-09 | **0.8412** |
| all 7 mirrors Al | 5.1066e-09 | 9.0336e-09 | **0.5653** |

So the reported coating sensitivity on that chain is **−3.2%** where the
ray-level truth is **+59%** — the wrong *sign*, because only the first
mirror's coating precedes the seed leg. A floor computed there and quoted
without this caveat would say "coating the mirrors slightly improves
contrast," which is backwards.

Handling, for review:

- The function computes both fractions itself, per input state, reports
  `.scope.carried` / `.worst` / `.full_chain`, and raises
  `macos:pol_contrast_floor:tranche1`.
- Per-state, not aggregated: the power-weighted mean over the x and y runs on
  Rx_Coro comes to 1.02 and would have declared the chain healthy. The worst
  state is reported instead.
- `tPolContrastCoro/test_coating_sensitivity_is_not_trustworthy_here` pins
  the wrong-sign number and the carried fractions **deliberately**, so
  Tranche 2 has to come back and change them.
- Report §5.6 states it; §5.7 labels the coronagraph floor a lower bound;
  the coverage section carries it as an open item.

The substantive deliverable therefore comes from `Rx_Cass_FarField`, where
both mirrors precede the single far-field leg and `carried == 1.00000`:

| configuration | annulus-mean cross contrast (10–40 px) | cross power vs uncoated |
|---|---|---|
| uncoated (perfect conductor) | 9.59e-12 | — |
| bare Al, 200 nm | 3.34e-10 | 27.9× |
| MgF₂ 110 nm / Al 200 nm | 1.99e-09 | 151.3× |

The design-rules line reads: *on this train the coating, not the geometry,
sets the polarization floor* — bare metal costs about a decade and a half
over the geometric term, and a quarter-wave protective overcoat costs most of
another. **Judgment wanted:** those are the model's numbers with tabulated Al
at 632.8 nm; I have not cross-checked the MgF₂/Al retardance against a vendor
curve or a published protected-Al polarization-aberration result, so I would
not put the 151× in a budget document without that step.

---

## 4. One bug, and what found it

The coherency matrix was first written `C = [Ex'*Ex, Ex'*Ey; Ey'*Ex, Ey'*Ey]`.
In MATLAB `'` conjugates its LEFT operand, so that is `conj(C)`, whose
dominant eigenvector is the CONJUGATE analyzer. For any LINEAR input state
the analyzer is real and the two are identical — x, y and 45° all passed. For
a circular input it is exactly ORTHOGONAL: reported cross/co jumped from
1.4e-06 to 7.13e+05.

It was found by the calibration run that fed the test tolerances, not by
reasoning, and the circular state is now a CI gate with its own non-vacuity
companion (`test_analyzer_would_be_wrong_if_conjugated`). Same lesson as
tVecChain's 45°/circular states: a polarization gate run only on linear
states can be blind to a conjugation.

**Corollary worth recording — Finding 1's original evidence was an
artifact.** The 2026-07-26 packet justified output-referencing with "a ~50%
cross fraction on a system whose measured diattenuation is 1.2e-15." On the
corrected engine that is gone: the r_p sign defect *was* the 50%. On these
coaxial fixtures the mean output state now equals the input state to
~1e-16 (`|1 − |a·v||` at 0 … 1.1e-16). Finding 1's discipline is still right
and still what the code does — it just cannot be demonstrated on any current
fixture, and the report does not claim it can. A fixture with a real
geometric rotation (an off-axis fold train) would demonstrate it; not built
here.

---

## 5. Infrastructure change — polval is now per-model-size

2c cannot run at model 128: `Rx_Cass_FarField` and `Rx_VecChain` declare
`nGridpts=256`, `Rx_Coro` declares 511. The driver's "one model size, one
session" constraint was therefore hit head-on, and it was split the way its
own README said to:

- `run_pol_validation(polvalDir, model)`, model ∈ {128, 256, 512}, each
  writing `generated/parts/numbers_<model>.json`;
- `gate_limits(model)` — per-size threshold tables, so a missing measurement
  in one group cannot be masked by another group's tokens;
- `merge_numbers.py` — refuses to merge parts with colliding tokens or with
  disagreeing provenance (different session), and stamps the model sizes as a
  list;
- `make_polval.sh` loops `MODELS="128 256 512"`.

Existing 128 numbers are unchanged by the split (33 gate thresholds, same
values).

---

## 6. Re-run results

| Gate | Result |
|---|---|
| `run_mmacos_tests.sh polfloor` | 20 pass, 0 fail (256: 14, 512: 6) |
| mmacos full suite | **440 pass, 0 fail** (fast 289 / masks 62 / freeform 60 / proper-512 10 / **pol-contrast-512 6** / proper-1024 13) — was 420 before the 20 new tests |
| polval regen | 48 gate thresholds pass across the three sizes (33 + 10 + 5); 119 tokens merged; `check_polval` clean; `make polval` builds docx + HTML |
| polval split is a no-op | all **81** pre-existing tokens bit-identical to the committed values; the eight pre-existing figures pixel-identical |

Not re-run this session, and stated as such: pymacos/ifx, PROPER-compare and
GMI. **Nothing in this session touches the engine** — `pol_contrast_floor` is
pure binding-layer MATLAB, the only Fortran-adjacent change is the polval
driver's invocation pattern — so those three suites cannot have moved. If the
reviewer wants them anyway, the commands are in `external.json`.

---

## What to look at first, if tokens are short

1. **§3** — the Tranche-1 shortfall and the wrong-sign coating sensitivity.
   That is the finding, and the decision it invites is whether the Rx_Coro
   floor should appear in the report at all or wait for Tranche 2.
2. **§4's corollary** — Finding 1's justification came from the pre-fix
   engine. Nothing depends on it now, but the record should not keep citing
   it as evidence.
3. **§3's last paragraph** — whether the 151× MgF₂/Al number is safe to state
   as a design-rules result without an external cross-check.

Everything else is gate-backed and mechanical.
