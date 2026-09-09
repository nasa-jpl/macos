# The "center-channel speckle" in `dwdsurf_*_5fov` — measured attribution

Terminal Opus for Dave, 2026-09-08.  Answers `BRIEF_to_sens_noise_center.md`
(CCL, same day).  Scope as briefed: dw/dsurf, one deck, centre zoom, centre
field; stop at the first measured attribution.

**Verdict in one line.**  The speckle is not a response.  It is the
finite-difference noise floor, and the floor exists because **the engine's
re-traced OPD is not idempotent on this deck at its nominal field**: two
identical traces alternate between two results differing by up to 6 ulp of the
24 459 mm accumulated optical path.  Element 4 (CenterSegment) is simply the
only channel whose real response never rises above that floor.

---

## 0. Provenance

| | |
|---|---|
| engine | macos `dev-candidate` **990ee5e** (contains `44fc362`, the STOP-handedness fix, and the 4-probe FEX), gfortran release, `build_release_gfortran/lib/libsmacos.a` 2026-09-08 07:28 |
| binding | MACOS_resources `dev-candidate` **3034fff**, `mmacos/src/mmacos.mexa64` 2026-09-08 07:28 (linked against the above; engine sources clean at HEAD) |
| deck | `~/dev/MACOS_sandbox/sens_noise/jwst_ote_designc.in` — Luis's copy.  `diff` against the committed fixture `mmacos/templates/50_sensitivities/zoom_5x5/jwst_ote_designc.in` shows **exactly three added header lines and nothing else**: `PgplotImage= Color`, `UseChfRay4OPD= Y`, `ApStop= -1.499707248121 84.72440986213 36939.79268433` |
| run state | model 512, `nGridPts` 63, stop elt 25 (FSM), OPD read at elt 27 (`wf_elt_auto` → nElt−1, the ExitPupil `Return`), 2301 rays traced / 2184 passing / 117 obscured, `delta` 1e-6, central differences |
| scripts | `~/dev/MACOS_sandbox/sens_noise/work/` (`s0`…`s11`, `cli_toggle_probe.py`) — not committed |

Nothing was regenerated in `zoom_5x5`; no engine or test change was made.

---

## 1. Step 0 — which "center"?  **Both, and they are different centres**

The brief offered two readings.  Measured, they are *both* right and they name
two different things, which is why the symptom reads as one channel.

### (a) The element: it is the CenterSegment's channel

Single-field, centre-zoom `dw_dsurf` on elts 4 / 5 / 23 / 24 (`s0`), raw
columns, OPD in BaseUnits (mm) per mm of Kr / per unit Kc:

| channel | norm | rms | max abs |
|---|---|---|---|
| **Elt 4 Kr** | 3.43880e-04 | **7.35836e-06** | 1.81899e-05 |
| **Elt 4 Kc** | 5.57692e-05 | **1.19335e-06** | 1.09139e-05 |
| Elt 5 Kr | 9.54901e-02 | 2.04330e-03 | 1.52722e-02 |
| Elt 5 Kc | 4.25050e+00 | 9.09524e-02 | 9.40448e-01 |
| Elt 23 Kr | 7.41092e-01 | 1.58579e-02 | 3.26254e-02 |
| Elt 23 Kc | 7.41572e+00 | 1.58682e-01 | 5.01021e-01 |
| Elt 24 Kr | 1.13957e-01 | 2.43846e-03 | 6.35373e-03 |
| Elt 24 Kc | 9.06222e-01 | 1.93914e-02 | 6.68460e-02 |

Reconstructed on the canvas (`work/s0_columns_raw.png`, each panel on its own
colour scale, as Luis's `show_dw_channels` draws them): elts 23 / 24 are smooth
full-pupil responses, elt 5 is a clean single-segment response, and **elt 4 is
pure salt-and-pepper over the whole pupil** — the reported speckle.

A roughness statistic makes it quantitative — `R = rms(nearest-neighbour
difference) / rms(column)`; `R → 0` for a smooth response, `R → √2 = 1.414`
for uncorrelated pixel noise:

```
Elt 4 Kc  R = 1.3958      <- pixel noise
Elt 4 Kr  R = 0.2264      (its energy is 97% piston; the rest is the same noise)
Elt 5 Kr  R = 0.4901      Elt 5 Kc  R = 0.5872
Elt 23 Kr R = 0.0594      Elt 23 Kc R = 0.1043
Elt 24 Kr R = 0.0513      Elt 24 Kc R = 0.0588
w_nom     R = 0.1012      (the nominal wavefront, for reference)
```

### (b) The block: the floor exists only at the CENTRE FIELD

Measured directly — a **zero-poke** difference (two traces, no perturbation,
divided by 2·delta) at every one of the 25 blocks the driver harvests (`s9b`):

| | f0 (centre) | fUL | fUR | fLL | fLR |
|---|---|---|---|---|---|
| z0  | **379 px, rms 1.19335e-06** | 0 | 0 | 0 | 0 |
| zUL | **364 px, rms 1.25401e-06** | 0 | 0 | 0 | 0 |
| zUR | **363 px, rms 1.18315e-06** | 0 | 0 | 0 | 0 |
| zLL | **377 px, rms 1.22715e-06** | 0 | 0 | 0 | 0 |
| zLR | **353 px, rms 1.20472e-06** | 0 | 0 | 0 | 0 |

`0` is exact — at every off-axis field the trace is bit-reproducible and the
floor is identically zero.  So the 5 centre-field blocks (one per zoom state)
are the only ones carrying speckle, in **every** channel; and within them elt 4
is the only channel where the speckle is the whole content.  That is exactly
Luis's "fixed for a given (field, zoom), in the center channel only".

The stop-call count is not a confound: 1, 2 or 3 `stop(25)` calls give the same
379 toggling pixels at the centre field and the same exact 0 at the corner
(`s10`).  Nor is it "on axis" as such — a field scan along +x (`s11`) gives
toggling at 0, 1e-9, 1e-8 rad, **none** at 1e-7…3e-5, toggling again at 1e-4,
none at 2.90888e-4.  Whether a given source state lands on the flip is
effectively arbitrary; it is fixed once the state is fixed.

---

## 2. The attribution

### 2.1 The floor is the same in every column, including channels the poke cannot reach

Segment 5's footprint is 124 of the 2184 valid pixels (5.7% ≈ one of 18
segments).  A poke of segment 5 cannot change the path of a ray that misses
segment 5, so **everything outside that footprint is artifact by construction**
(`s5` §2):

| channel | on-footprint rms | OFF-footprint rms |
|---|---|---|
| Elt 5 Kr | 8.5753e-03 | **1.1916e-06** |
| Elt 5 Kc | 3.8171e-01 | **1.1916e-06** |
| Elt 4 Kc (whole pupil) | — | **1.1933e-06** |

The off-footprint content of a live column equals the whole of the elt-4
column.  One floor, under every channel.  It is invisible in elt 5 / 23 / 24
only because their own colour scale is 3–5 decades larger — in the picture it
shows as the scattered saturated dots outside segment 5.

### 2.2 A zero poke reproduces the elt-4 column

`SurfChannel.apply(0)` twice — set the parameter to nominal, `modify()`,
trace; repeat — then divide by 2·delta exactly as the FD engine does (`s1` B):

```
elt 4 Kr ZERO-poke column: rms 1.1933e-06  max 1.0914e-05  nnz 379
elt 4 Kc ZERO-poke column: rms 1.1933e-06  max 1.0914e-05  nnz 379
   vs the real elt-4 (ptt-removed) column: |corr| 0.999154, rms ratio 1.0009
```

The elt-4 "sensitivity" is the zero-poke difference to within 0.1%.

Corroborating: the Kr and Kc columns, piston removed, agree to
**max abs diff 1.6941e-21 on an rms of 1.1933e-06** — 1.4e-15 relative, i.e.
bit-identical.  Two physically independent pokes cannot give the identical
residual; it is not a residual of the pokes.

### 2.3 The mechanism: the re-traced OPD is not idempotent — a strict 2-cycle

Ten nominal traces in one process, `modify()` between each, nothing perturbed
(`s1` A):

```
trace  1 : max|dW| 0.0000e+00   nnz     0 / 2184
trace  2 : max|dW| 2.1828e-11   nnz   379 / 2184   rms 2.3867e-12
trace  3 : max|dW| 0.0000e+00   nnz     0 / 2184
trace  4 : max|dW| 2.1828e-11   nnz   379 / 2184
   ... alternating to trace 10.   V2==V4, V3==V5: a strict 2-cycle.
```

`2.1828e-11 mm = 6.00 × eps(2.4459e4 mm) = 6 ulp` of the accumulated optical
path; the typical toggling pixel moves 1–2 ulp.  Divided by `2·delta = 2e-6`
this is a column grain of `eps(OPL)/(2·delta) = 1.81899e-06` — and the measured
floor rms is 1.1933e-06, max 1.09139e-05 (= 6 grains).  The elt-4 columns take
integer multiples of that grain in the range −10…+1 (Kr) and −6…+5 (Kc); a live
column (Elt 5 Kc) spans ~5.2e5 of them.  (Integrality itself is trivial — every
OPD value is a multiple of the path's ulp — so the *count* of grains is the
statistic, not the integrality.)

Because the toggle is deterministic, the resulting column is a **fixed**
pattern for a fixed (field, zoom).  That is what made it look physical.

### 2.4 It is engine-side, not a binding artifact

Driven through the interactive `macos` CLI over a pty (model 512,
`mod nGridpts=63`, `stop elt 25 0,0`, then `OPD 27` four times with `MOD`/`quit`
between — the interactive path to the same flag reset `MODIFY` performs;
`cli_toggle_probe.py`):

```
RMS OPD error      Average OPD
6.846100093D-06    2.8293106046702711D-07
6.846100131D-06    2.8293109045036889D-07
6.846100093D-06    2.8293106046702711D-07
6.846100131D-06    2.8293109045036889D-07
strict 2-cycle: True
```

Bit-for-bit the same two states the mex alternates between.

### 2.5 Scope of the non-idempotency — what it is NOT

* **Not Luis's three header lines.**  The committed fixture (no
  `UseChfRay4OPD=`, no `ApStop=`), stop elt 25, toggles identically: 6.00 ulp,
  2-cycle confirmed.  What the chief-ray reference *does* change is the spatial
  extent — with `UseChfRay4OPD= Y` only the 379 rays whose own path toggles
  differ; with the default aperture-mean reference **all 2184** differ, because
  the mean itself moves and shifts every ray by a common amount.  The chief
  reference is the better of the two here.
* **Not the stop.**  Stop elt 25, object-space stop from the deck's `ApStop=`,
  and no stop command at all all toggle (379 / 370 / 379 pixels, 6.00 ulp).
* **Not the FreeForm promotion.**  The pre-promotion deck (`8316b68`, all 27
  surfaces `Conic`) toggles too — 2184 pixels, 5.787e-11 max, 2-cycle.
* **Not universal.**  `e5hex1` (7 hex segments) and `Rx_Cass_NS` are exactly
  idempotent — 0 toggling pixels, max abs difference 0.000e+00.

**What actually alternates inside the trace was not determined.**  That is the
open engine question and is left for Dave / CCL: a deterministic 2-cycle points
at state carried from one trace into the next (a solver warm-start, or a
retained intermediate) rather than at random round-off, and it is deck- and
source-state-dependent.

---

## 3. Recommendation

**Set the finite-difference step from a measured floor; concretely, raise the
`dw_dsurf` step on this rung from `delta = 1e-6` to `1e-4`.**

The floor scales exactly as 1/delta (the numerator is delta-independent) while
the signal does not move.  Delta scan, elts 4 and 5, column RMS (`s3`):

| method | delta | Elt 4 Kr | Elt 4 Kc | Elt 5 Kr | Elt 5 Kc |
|---|---|---|---|---|---|
| central | 1e-07 | 1.19335e-05 | 5.57126e-05 | 2.04384e-03 | 9.09501e-02 |
| central | 1e-06 | 7.38795e-06 | 1.19335e-06 | 2.04336e-03 | 9.09524e-02 |
| central | 1e-05 | 7.38795e-07 | 1.19335e-07 | 2.04334e-03 | 9.09523e-02 |
| central | 1e-04 | 1.12769e-06 | 2.16293e-08 | 2.04335e-03 | 9.09523e-02 |

The elt-4 columns fall by exactly ×0.1 per decade of delta — **including the
Kr column's piston**, so that is noise too, not a chief-ray-reference response
as it first looks.  The live columns are flat to five figures across four
decades.  The zero-poke floor confirms the law with the same mantissa:
`1.19335e-06` at 1e-6 → `1.19335e-08` at 1e-4 (max `1.09139e-05` →
`1.09139e-07`).

Column-wise, not just in RMS (`s8` A) — `max|col(1e-4) − col(1e-6)| / rms`:

| channel | rms @1e-6 | rms @1e-4 | max abs dcol / rms |
|---|---|---|---|
| Elt 5 Kr | 2.04330e-03 | 2.04335e-03 | 5.929e-03 |
| Elt 5 Kc | 9.09524e-02 | 9.09523e-02 | 1.316e-04 |
| Elt 23 Kr | 1.58579e-02 | 1.58542e-02 | 1.019e-03 |
| Elt 23 Kc | 1.58682e-01 | 1.58682e-01 | 7.783e-05 |
| Elt 24 Kr | 2.43846e-03 | 2.43808e-03 | 4.238e-03 |
| Elt 24 Kc | 1.93914e-02 | 1.93915e-02 | 4.887e-04 |

Every one of those differences is the **1e-6 run's own floor**: the floor's max
is 1.09139e-05, so 1.09139e-05/2.04e-03 = 5.3e-03 against the measured 5.9e-03
for Elt 5 Kr, 6.9e-05 vs 7.8e-05 for Elt 23 Kc, and so on for all six.  Nothing
non-linear appears at 1e-4.  Scope: measured on this deck, this rung; the
durable form of the rule is to *measure* the floor (two extra traces per block,
2.4% of an 85-trace block) and pick delta and a reporting threshold from it,
rather than to hard-code 1e-4.

### It also fixes the dead-channel bookkeeping, with no code change

`flag_zero_norm_channels` flags a block at `1e-6 × median live block RMS`.
On the dw/dsurf rung today (`s3`, `s8` C):

| delta | elt-4 block RMS | median live | ratio | flagged |
|---|---|---|---|---|
| 1e-6 | 1.1922e-06 | 4.5562e-01 | 2.617e-06 | `[]` — **not flagged** |
| 1e-4 | 1.1922e-08 | 4.5562e-01 | 2.617e-08 | `4` — flagged |

So at the current step the committed `run_dwdsurf_5zoom_5fov.m`'s
`flag_zero_norm_channels` → `drop_channels` step is a **no-op**: elt 4 sits 5
decades below the median live optic, not 6, so its two noise columns ship in
the Jacobian.  (The README's number-free claim is stated for the dw/dx rung,
where elt 4 really is ~5e-7 of live and the flag does fire.)  Raising delta
puts elt 4 back under the existing threshold on its own merits — response-keyed
still, no element number hard-coded.

---

## 4. Two things noticed en route (neither is what Luis saw)

1. **The committed `dwdsurf` artifacts in `zoom_5x5` are stale.**
   `find_powered_elts` now returns `[4 5 … 24]` — 21 elements, **42 channels** —
   since `Segment` became powered-capable (Dave's 2026-09-05 ruling), and
   `run_sensitivities` calls `dw_dsurf_multi` with no `'elts'`, so the harvest
   auto-discovers all of them.  The committed
   `dwdsurf_5zoom_5fov_..._sens_report.txt` records `54585x4` — four channels,
   and the pre-`stop-enforced-chief` row count.  The driver header and the
   README row ("SM (M2) + TM (M3) … = 4 channels") say the same stale thing.
   Not regenerated, per the brief.

2. **A one-shot stale first trace at elt 27, immediately after a trace to elt
   26** (`s5` §3, `s6`, `s7`).  `trace(26); opd()` then `trace(27); opd()`
   returns rms 3.72982e-05 mm where the correct value is 6.85038e-06 — a factor
   5.4.  The next `trace(27)` recovers the correct map exactly.  Only that pair:
   23→24, 24→25, 25→26 and 27→28 are all clean, as is a fresh load.  **The
   sensitivity harvest never sees it** — `local_wf` always traces to the same
   `wf_elt`, and six `modify()`+`trace(27)` cycles stay inside the 6-ulp floor
   (verified, `s7` B).  Flagged for Dave / CCL with the reproducer; not chased.

---

## 5. Files

Scripts kept (not committed) under `~/dev/MACOS_sandbox/sens_noise/work/`:

| file | what it measures |
|---|---|
| `s0_which_center.m`, `s0b_analyze.m` | Step 0: column norms, roughness, `s0_columns_raw.png` |
| `s1_repeat.m` | 10-trace repeatability, zero-poke FD, Kr/Kc identity |
| `s2_localize.m` | deck / stop / OPD-reference variants, `s2_toggle_map.png` |
| `s3_delta.m`, `s8_reco.m` | delta scan, linearity, `flag_zero_norm_channels` |
| `s5_where.m`, `s6_elt27.m`, `s7_stale.m` | read-element scan; the elt-26→27 anomaly |
| `s9b_blocks.m`, `s10_confound.m`, `s11_fieldscan.m` | floor per block; confound and field scans |
| `cli_toggle_probe.py` | the pty-CLI cross-check (engine-side confirmation) |
| `jwst_ote_designc_prepromo.in` | the pre-FreeForm-promotion deck, from `8316b68` |
| `jwst_ote_designc_meanref.in` | Luis's deck with `UseChfRay4OPD= Y` stripped |
