# Polarization Phase 0 — Audit & Scoping Findings

> Deliverable for `PLAN_POLARIZATION.md` Phase 0 (CCMac lane). No code changes.
> Line references are against `pol-core` (macos) as of this note. Each finding
> gates a downstream decision; the **bottom-line** for each is called out.

---

## 0.1 Vector-propagation leg coverage

Enumerated every physical-optics propagation leg reachable from `CPROPAGATE`
(`propsub.F`) and the DFT legs (`dftsub.F`). Dispatch is by `PropType(iEm1)`.

| Leg | PropType | Routine | Field arg | Vector? |
|---|---|---|---|---|
| Near-field sphere→plane | 2 | `NFPROP` | `WFElt(1,1,iWF)` | **scalar only** |
| **Far-field sphere→plane** | 3 | `PFFPROP` (vec) / `FFPROP` (scalar) | `WFElt(1,1,1)` all 3 planes / `(1,1,iWF)` | **VECTOR** (only leg) |
| Near-field plane→plane | 4 | `PPPROP` | `WFElt(1,1,iWF)` | scalar only |
| Near-field 2×RS sphere→sphere | 5 | `NFPROP` | `WFElt(1,1,iWF)` | scalar only |
| Spherical/Fresnel legs | 6,7,8 | `SFPROP`/`FRPROP` | `WFElt(1,1,iWF)` | scalar only |
| NF DFT sphere→plane | (PropType iEm1) | `NFPropDFT` | `WFElt(1,1,iWF)` | scalar only |
| NF DFT Siegman-Sziklas | 14 | `NFPropDFT` | `WFElt(1,1,iWF)` | scalar only |
| FF DFT | 15 | `FFPropDFT` | `WFElt(1,1,iWF)` | scalar only |

- The **only** `IF (ifVecDif3)` propagation branch is at `propsub.F:1604`, calling
  `PFFPROP` (`propsub.F:2273-2297`), which loops `K=1,3` and runs an independent
  `DSWAP3`/`DFOURN`/`DSWAP3` on each of `WFElt(:,:,1..3)` = Ex/Ey/Ez.
- The field-assembly branch (`propsub.F:1373-1426`) does load all 3 planes when
  `ifVecDif3`, but on **any non-PropType-3 leg** the propagator only touches
  `WFElt(:,:,iWF)` (a single plane, `iWF=iCurWFElt`). The other two component planes
  are carried but **not propagated** — they are stale after the first scalar leg.
- `FFObscure` (Lyot/mask apply on the diffraction grid, `propsub.F:2301-2473`, fires
  only for `EltID==7`) writes `WF(iWF,jWF)=0` on **one plane** — it does not zero the
  other two component planes either.

**BOTTOM LINE.** Vector diffraction closes **only for a single far-field FFT hop**
(pupil → focus). It does **not** close a coronagraph chain
(pupil → FPM → Lyot → focal), which requires vector fidelity across near-field and
plane-to-plane legs plus a vector-aware `FFObscure`. **The VVC acceptance test (Phase 3)
is gated: "vectorize the chain" is a real Phase 3a**, scoped as: (a) per-component
propagation in NFPROP/PPPROP/SFPROP/FRPROP/NFPropDFT/FFPropDFT, (b) 3-plane `FFObscure`,
(c) the field-assembly/`CumLStart` bookkeeping per component. **Track A (IFO PSI) is
unaffected** — its detector-mode fields close with a single far-field hop, so it needs
only Phase 1.

---

## 0.2 `srtrace.F` `ifPol` scoping

- `srtrace.F` contains exactly one routine, `SRTRACE_Test` (`srtrace.F:25-449`), a
  **test driver**. Its `Logical,parameter :: ifPol=.false.` (`srtrace.F:67`) is local to
  it. Its **only caller** is `macos_cmd_loop.inc:6295`, which is inside an `#if 0` block
  (dead — the `SRTRACE` command is disabled).
- The **production** single-ray / chief-ray path is `SRTrace` and `CRTrace` in
  `tracesub.F` (`SRTrace` at `tracesub.F:2805`), both of which take `ifPol` as an
  **argument** and honor it. All live callers (SPOT/LocalCoord/STOP/FDPN via
  `macos_cmd_loop.inc:1680`, `stop_set.F`, `tracesub.F`) pass the real `ifPol`.

**BOTTOM LINE.** Non-issue. The hardcoded `ifPol=.false.` lives only in dead test code;
the Jones-pupil pupil-grid trace uses `CTRACE`, and single-ray probes use
`CRTrace`/`SRTrace`, all of which already thread `ifPol`. **No lifting required.** The
plan's expected answer ("not blocking") is confirmed — and it's not even a
declaration/threading change, because the production path never had the parameter.

---

## 0.3 Coating / wavelength semantics

- **Parse** (`msmacosio.inc:2648-2666`): each `Coating=` layer reads
  `(n, κ, thk_waves)`; thickness is converted to **physical** immediately via
  `EltCoatThk = thk_waves · Wavelen / IndRefArr` (`:2660-2661`), using the `Wavelen`
  **current at parse time**. Boundary media `IndRefArr(0,·)`/`(nCoat+1,·)` are snapshotted
  from the **adjacent elements' `IndRef`/`Extinc`** (`:2663-2666`) — so `Coating=` must
  appear **after** `IndRef=` in the element block.
- **Trace** (`elemsub.F:512-516`): layer round-trip phase uses the **current**
  wavelength: `C2 = i·2·(2π/λ_now)·d·Nb`, `d` physical. So a loaded stack is physically
  fixed and **broadband sweeps are already correct** — λ_now varies, `d` doesn't.
- **SAVE** (`iosub.inc:1754-1765`): inverts to waves via
  `thk = EltCoatThk · IndRefArr / Wavelen` using `Wavelen` **current at save time**.

**BOTTOM LINE.** Semantics are correct and confirmed. The one genuine wrinkle: the
`.in` "thk in waves" number is *waves at whatever `Wavelen` was current when the line was
parsed/written* — so a load→(change `Wavelen`)→save round-trip preserves the **physical**
stack but writes a **different waves number** than the original file. **Phase-1 `coat_set`
takes physical thickness and sidesteps this entirely.** *Test to add in Phase 1:*
load → save → reload at a different `Wavelen`, assert the physical `EltCoatThk` arrays
match to round-off (pins the SAVE inversion + documents the coupling). Also document the
`Coating=` after `IndRef=` ordering requirement.

---

## 0.4 Precision audit

Swept the entire polarization field chain for any single-precision stage:

| Stage | Type | Source |
|---|---|---|
| Per-ray field `RayE` | `Complex*16` | `elt_mod.F:450` |
| Diffraction grid `WFElt`, `JmatElt`, `WFbuff` | `Complex*16` | `elt_mod.F:448,450` |
| Coating recursion locals (`RP,RS,RPnew,...,NaCmplx,ccfb_arr`) | `Complex*16` | `elemsub.F:275-278` (Reflector), `780-783` (Refractor) |
| FFT `DFOURN` data + trig accumulators | `Complex*16` / `Real*8` | `sunsub.F:149-152` |
| `DSWAP3` | `Complex*16` | `sunsub.F:56` |
| Output buffers `RaySpot`,`OPDMat` | `Real*8` | `macos_mod.F:37` |
| `PixArray`, `Cmatrix`, `DrawRayVec` | `SREAL` = **`REAL*8`** | `realtype.h:12` ("SREAL replaced REAL*4, 09-2007") |

**BOTTOM LINE.** No single-precision stage anywhere in the polarization path. `SREAL` is
a `REAL*8` macro throughout (the 2007 conversion). At 1e-10 contrast (field amplitudes
~1e-5) the chain is safe. Nothing to fix.

---

## 0.5 Coating subsystem inventory — TWO models

| | Model A (polarization path) | Model B (NS refractive path) |
|---|---|---|
| Keywords | `Coating`, layer lines `(n,κ,thk)` | `nCoatElt`, `CoatIndxElt`, `CoatThkElt` |
| Storage | `EltCoat(0:mElt)`, `IndRefArr(0:mCoat,·)`, `ExtincArr(0:mCoat,·)`, `EltCoatThk(mCoat,·)` (`mCoat=10`) | `nCoatElt(mElt)`, `CoatIndxElt(20,·)`, `CoatThkElt(20,·)` |
| Complex index? | **Yes** (`ExtincArr` = κ per layer) | **No** (real index only) |
| Parse | `msmacosio.inc:2648-2666` | `msmacosio.inc:2855-2869`, `macosio.F` MOD |
| Applied | `Reflector`/`Refractor` Fresnel + thin-film recursion (`elemsub.F:432-547`, `1080-1155`) | `AirGap` transmittance (`elemsub.F:1675-1735`), via `Refractor`/`NSRefractor` |
| SAVE | `iosub.inc:1754-1765` | `iosub.inc:1740-1744` (+ a pre-existing `nCoatElt` size ToDo at `:1740`) |
| Gated on `ifPol`? | Yes | Yes |

Both are disjoint (separate keywords, storage, code paths). Model A is the one the
Jones-pupil work rides (full multilayer, complex index, R and — dead-coded — T). Model B
is a real-index single-gap/stack transmittance helper for non-sequential refraction only.

**BOTTOM LINE / DAVE DECISION.** Phase 1 `coat_set`/`coat_get` should target **Model A**
(the polarization path with complex-index support). **Open for Dave:** unify B into A, or
leave B as the NS-refractive helper and extend only A (incl. the Phase-4 spatial map)?
Recommendation: leave B alone, extend A — unification is a separable cleanup with its own
regression surface (the `nCoatElt` SAVE ToDo lives there).

---

## 0.6 Phase-convention consistency

Derived the single (ω,k) convention from the code, not legislated:

- **Propagation to a surface** (`elemsub.F:387-390`): `S1 = -2π·L/λ`,
  `C1 = exp(i·S1·Na) = exp(-i·2π·L·N/λ)` with `N = n − iκ`. Over path `L`:
  - real part `exp(-i·2π·n·L/λ)` → phase = **−kL**, decreases as `L` grows;
  - imaginary part `exp(-2π·κ·L/λ)` → **decay** for κ>0. ✓
- **Coating layer** (`elemsub.F:514-516`): `C1 = exp(-i·2·(2π/λ)·d·Nb)` — same sign, the
  factor 2 is the reflection round trip through the layer.

Both consistent with a traveling wave `exp[i(ωt − kz)]` — i.e. the **exp(+iωt)
time-harmonic / exp(−ikz) spatial** convention, with absorbing index **`N = n − iκ`,
κ>0 = loss**. This matches (a) the 2026-07-25 IFO empirical finding that *engine field
phase advances as OPL shortens* (shorter L → larger −kL → phase increases ✓), and (b)
the pymacos↔PROPER `opd_sign_flip=True` reconciliation (PROPER uses the opposite
exp(−iωt) convention).

**BOTTOM LINE.** Convention is self-consistent and pinned:
- Time-harmonic: **exp(+iωt)**; spatial propagator **exp(−ikz)**.
- Absorbing index: **n − iκ**, κ>0 = loss.
- Jones storage basis: **linear (x,y)** (Phase 2 decision), circular by unitary change.
*Test to add:* a lossy-Al bare reflectance check against tabulated (n,κ) to confirm the κ
sign end-to-end. This convention block goes into a new engine-polarization section of
`macos_f90/CLAUDE.md` in Phase 1.

---

## Summary — what Phase 0 changes about the plan

1. **Phase 3a (vectorize the coronagraph chain) is confirmed real and is the plan's
   largest schedule item.** Only the far-field FFT hop is vector-capable today. Scope
   recorded in 0.1. VVC acceptance depends on it; Track A does not.
2. **`srtrace.F` `ifPol` is a non-issue** — dead test code; production paths already
   thread `ifPol`. Remove this from the Phase-1 worry list.
3. **Precision is clean** end-to-end (`SREAL`=REAL*8; all field/coating/FFT stages
   Complex*16/Real*8). No fix needed.
4. **Coating semantics correct**; the only wrinkle is the waves-at-parse-λ coupling,
   which `coat_set`(physical) sidesteps. One round-trip test to add in Phase 1.
5. **Two coating subsystems**; Phase 1 targets Model A. Unify-vs-extend is a Dave
   decision (recommendation: extend A, leave B).
6. **Phase convention pinned** (exp(+iωt)/exp(−ikz), n−iκ). Goes into `CLAUDE.md` + a
   lossy-Al sanity test in Phase 1.

**Gate:** ready for review. No engine code touched. On approval, proceed to Phase 1
(`pol_set`/`pol_get`, `vecdif_set`, `coat_set`/`coat_get`, `rayfield_get`, `Ex0Ey0` fix).
