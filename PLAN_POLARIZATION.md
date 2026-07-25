# Expose & extend MACOS polarization physics — IFO metrology + coronagraph modeling

> Revised from the Opus draft 2026-07-25 (CC review): engine claims re-verified against
> source, two factual corrections folded in (coating-wavelength semantics, Ex0Ey0
> keyword state), the **interferometer/bench application added as a co-equal track**
> (Dave: "add polarization to this model, AND to the coronagraph models"), binding
> mechanics aligned with the actual mmacos codegen path, and the validation ladder
> extended with an IFO self-measurement cross-check.

## Context — two applications, one core

**A. Interferometer / bench metrology** (`examples/design/bench_ifo*`, the VSG-class
testbeds). The rigs run 45° plate beamsplitters and folds. The 2026-07-25 AOI sweep
found scalar metrology AOI-neutral in the compensated design — **the entire case for
smaller AOI lives in polarization** (R_s/R_p diattenuation, coating retardance), and
today none of it is modeled. Polarization sets fringe-visibility loss, s/p phase
splitting (a PSI systematic), and enables the classic instrument architectures we
cannot currently build: PBS + quarter-wave-plate Twyman-Green, polarization
phase-shifting.

**B. Coronagraph contrast floor.** `CORONAGRAPH_DESIGN_RULES.md` §3/§9 flag
polarization aberration as *"a hard floor for 1e-10 contrast … does not respond to
scalar DM control."* The budget needs a **number** — the polarization-limited contrast
floor and its sensitivity to coating choice — not just maps.

MACOS already contains substantial polarization physics in the engine, but **none of it
is reachable from pymacos/mmacos**, so neither application can use it. This plan
exposes the existing physics first (Dave's directive), builds the shared Jones-pupil
machinery, then delivers per-application analyses, then (optionally) extends the engine
with polarizing elements and spatially-variable coatings.

**Deliverable framing.** Phases 1+2 are the committed, shippable core: *given a
prescription, return the Jones pupil, the polarization-aberration decomposition, the
IFO visibility/PSI-systematic numbers (track A) and the contrast floor (track B).*
Phases 3 and 4 add modeling reach, carry most of the schedule risk, and are
independent of each other — sequence them by need (PBS/waveplate IFO and VVC trades →
Phase 3 first; segmented-primary coating non-uniformity → Phase 4 first).

## Branching & execution sequencing (agreed with Dave 2026-07-25)

**Branches — split at the consumer boundary** (single-issue rules; stacked-branch
pattern per bench-builder/expose-beam):

1. **`pol-core`** — single-issue branch in **both repos** (macos: engine API +
   `Ex0Ey0` fix; MACOS_resources: mmacos/pymacos bindings + tests). Scope: Phases
   0–1 plus the application-neutral part of Phase 2 (Jones pupil, `pol_maps`,
   unitarity/2θ tests). **CUT 2026-07-25 and pushed:** macos `pol-core` off
   `expose-beam` (742ce29 — carries the `beam_set` template the plan cites; same
   commit as macos bench-builder), MACOS_resources `pol-core` off **`bench-builder`**
   (cea3d60 — Dave's call: Phase 1's track-A smoke test drives
   `macos.design.twyman_green`, which lives there, and bench-builder ⊇ expose-beam).
   Stacked PRs reduce as parents merge. Later phases (elements, spatial coatings)
   become their own follow-on single-issue branches/PRs; the phase boundaries are
   PR-sized on purpose.
2. **`pol-ifo`** — thin consumer branch **stacked on `bench-builder`**: Phase 2d, the
   Bench coating emission, `add_polarizer`/`add_waveplate` wiring, the
   `bench_ifo_pol` example. Needs both parents; lands after bench-builder and
   pol-core merge (stacked PR auto-reduces as parents land).

⚠ The public history rewrite landed ~2026-07-24: confirm the current base branch
name (`dev` vs `sls-dev`) before cutting branches, and cut only from a fresh
post-rewrite clone (never push from a stale pre-rewrite clone).

**Execution — split at the judgment boundary, not the phase boundary:**

- **CCMac (Opus 4.8): Phases 0 and 1, later Phase 4, and the polarizer/waveplate
  half of Phase 3.** Template-following with objective, machine-checkable gates
  (`beam_set`/`beam_get` pattern, `gen_mex_wrappers.py` codegen, AmplMat model,
  Reflector s/p projection to copy; round-trip identity, geometry invariance,
  crossed-polarizer extinction, GMI 6/6). Delegation is safe where passing the gate
  *is* the definition of done. The Phase-0 audit note doubles as CCMac's onboarding
  vehicle for the engine polarization internals.
- **Claude 5 (this environment): Phase 2's physics core, the VVC, and any Phase 3a.**
  Where the plan itself warns things go quietly wrong: double-pole basis
  construction, retardance branch handling, mean-vs-variation separation, convention
  arbitration, and un-templated engine surgery (vectorizing propagation legs). The
  2026-07-25 bench_ifo_dm diagnosis (a plausible "9 nm algorithm floor" that
  decomposed into a branch-cut artifact + a sign convention + real retrace physics)
  is the reference case for why these stay with the stronger model: an executor that
  stops at "tests pass, plausible number" ships a wrong budget.

**Sequence:** CCMac starts Phase 0+1 on `pol-core` now → review at the phase gate
here → Phase 2 lands on top of their Phase 1 → `pol-ifo` follows here (the
bench_ifo machinery and its 3e-5 pm validation history live in this environment).

---

### What already exists in the engine (verified 2026-07-25, line-checked)

- **Per-ray vector E-field** `RayE(3,mRay)`, `Complex*16` end-to-end
  (`elt_mod.F:450,846` — double precision confirmed; updated per-ray in `CPROPAGATE`
  at `propsub.F:507-509` and via the `Reflector`/`Refractor` `Evec→Eout` argument pair,
  `elemsub.F:206-216`).
- **Fresnel + recursive multilayer thin-film coatings**, complex index stored as
  `DCMPLX(n, −κ)`, in `Reflector` (`elemsub.F:432-547`) and `Refractor`
  (uncoated Fresnel + coated recursion with **live TP/TS transmission**).
  All **gated on `ifPol`** — inert when polarization is off.
  ⚠ `Reflector`'s transmittance sub-blocks are `if (.false.)` dead code — reflective
  elements return R only (fine for mirrors; relevant if anyone expects a coated
  reflector to report transmitted leakage).
- **Vector diffraction** (3 independent FFTs of Ex/Ey/Ez), **far-field FFT leg only**
  (`propsub.F:1380-1426, 1604`, gated `ifVecDif3`).
- **Source polarization** `Ex0/Ey0` (`src_mod`); CLI `POLARIZATION`/`NOPOLARIZATION`/
  `VECTOR`/`SCALAR` (`macos_cmd_loop.inc:1072-1114`). `VECTOR` requires `ifPol` first;
  the `POLARIZATION` path enables `ifVecDif3` only when `mWF≥3` — **all stock models in
  `macos_param.txt` have `mWF=3`**, so this is available at every model size, but the
  API must still assert it.
- **Coating wavelength semantics (verified — corrects the earlier draft).** The parser
  converts the Rx thickness spec (units of waves) to **physical thickness** at load
  time via `Thk·Wavelen/IndRef` (`msmacosio.inc:2648-2675`); the trace-time recursion
  then applies phase `2·(2π/λ_now)·d·N` with the **current** wavelength
  (`elemsub.F:512-518`). So a loaded stack is physically fixed and **broadband sweeps
  are already correct**. The SAVE writer already inverts the normalization
  (`iosub.inc:1747-1766`), so save/load round-trips preserve the physical stack. The
  real remaining wrinkles: (a) the Rx spec means "waves at whatever `Wavelen` was
  current at parse" — a file-order / cross-λ-interpretation coupling to document, and
  (b) `Coating=` must follow `IndRef=` in the file (the parser snapshots the boundary
  indices). Phase 1's `coat_set` API takes **physical thickness** and sidesteps all of
  it.

### Gaps this plan closes

- Zero binding exposure (can't even turn `ifPol` on from Python/MATLAB).
- `RfPolarizerElt(14)`/`TrPolarizerElt(15)` are **name-table-only stubs** (verified: no
  dispatch anywhere in `tracesub.F`/`propsub.F`/`elemsub.F`); `JmatElt(2,2,mElt)` is
  allocated/zeroed/deallocated and never otherwise touched (dead — confirmed).
- **`Ex0Ey0` prescription keyword is silently swallowed** (worse than dead): a live
  `ELSE IF (LCMP(VAR_NAM,'Ex0Ey0',6)) THEN GO TO 50` no-op at `msmacosio.inc:176`,
  next to the commented-out `PolBeam`/`PolSrc` block (150-174). A user setting source
  polarization in a `.in` gets no error and no effect.
- No Jones/Stokes/Mueller math, no polarization metrics, no contrast floor, no IFO
  polarization diagnostics.
- No waveplate/retarder, no polarizer, no vector-vortex mask.
- Coatings are per-element uniform; the one spatial map (`AmplMat`) is real-amplitude
  only.
- **`srtrace.F` hardcodes `ifPol` — as a `Logical,parameter :: ifPol=.false.`**
  (`srtrace.F:67`), so lifting it is a declaration + argument-threading change, not a
  flag flip. Whether it matters is a Phase-0 scoping question (srtrace serves the
  single-ray/chief-ray paths — SPOT beam frames, LocalCoord, sensitivity probes — not
  the pupil-grid trace the Jones pupil uses).
- **Vector propagation is not known to close the coronagraph chain.** Only the
  far-field FFT leg is vectorized. A VVC test needs pupil → FPM → Lyot → focal with the
  vector field preserved at every leg. Audited in Phase 0; if it doesn't close,
  "vectorize the chain" becomes its own Phase 3a.

---

## Conventions to pin down before any code

Cross-code polarization disagreements are almost always convention mismatches. Fix
once, document in `macos_f90/CLAUDE.md` (new engine-polarization section), assert in
tests:

| Convention | Decision |
|---|---|
| Time-harmonic phase | Derive from the engine, don't legislate: the coating recursion uses `N = n − iκ` with propagation factor `exp(−i·2kd·N)` (`elemsub.F:512-518`), and the 2026-07-25 IFO work established empirically that **engine field phase advances as OPL shortens** (frame-ratio diagnostic, exact to 1e-4 pm; same convention the pymacos↔PROPER `opd_sign_flip` reconciles). Phase 0 writes down the single consistent (ω,k) convention these imply and pins it with a test. |
| Absorbing index | n − iκ as implemented; κ > 0 is loss under the engine's phase convention (verify sign once in Phase 0 with a lossy-metal reflectance check against tabulated Al). |
| Jones storage basis | **Linear (x,y)**, with a basis-change utility for circular. Circular is a fixed unitary transform; storing linear and converting is strictly simpler than the reverse. |
| Reference frame | **Double-pole (dipole)** by default; local ŝ/p̂ and global as options. See Phase 2a. |
| Retardance sign | Fast axis leads. State it; VVC handedness depends on it. |
| Diattenuation | Intensity-based, D = (T_max − T_min)/(T_max + T_min) ∈ [0,1]. |

---

## Phase 0 — Audit & scoping (no new features; gates everything downstream)

Short, cheap, can reshape the plan. Deliverable is a findings note, not code. Items 3
and 6 were pre-answered during this review and only need the recorded confirmation.

1. **Vector propagation leg coverage.** Enumerate every propagation leg (`PFFPROP`,
   near-field/Fresnel, angular-spectrum, pupil-plane products, Lyot application) and
   record vector-capable vs scalar-collapse. Output: legs × capability table. **If the
   coronagraph chain does not close, Phase 3's VVC acceptance test is rescoped and
   "vectorize the chain" becomes Phase 3a.** Note for track A: the IFO detector-mode
   PSI needs only trace-to-detector fields — it closes with Phase 1 alone, so track A
   does not wait on this finding.
2. **`srtrace.F` `ifPol` scoping.** Which analyses route through srtrace (chief-ray/
   SPOT/LocalCoord/single-ray probes), and do any polarization use cases need them?
   Expected answer: not blocking for Jones-pupil work; record it.
3. **Coating/wavelength semantics — CONFIRMED during review** (see above). Remaining
   Phase-0 work: one empirical round-trip test (load → save → reload at a different
   `Wavelen`, compare physical stacks) to pin the SAVE-writer inversion, plus document
   the parse-order couplings.
4. **Precision audit.** `RayE` confirmed `Complex*16`; sweep the remaining chain (the
   vector FFT arrays, coating recursion locals — spot-checked `COMPLEX*16` in
   `Refractor`) and record. At 1e-10 contrast we track field amplitudes near 1e-5; a
   single-precision stage is fatal and invisible.
5. **Coating subsystem inventory — two models confirmed:** the polarization-path stack
   `EltCoat/IndRefArr/ExtincArr/EltCoatThk` and the separate `nCoatElt/CoatIndxElt/
   CoatThkElt` (+`AirGap`) model (both appear in the SAVE writer, `iosub.inc:1738-1766`
   — note the pre-existing `nCoatElt` size ToDo there). Document both; decide with
   Dave whether to unify or extend only the polarization-path one.
6. **Phase-convention consistency** (from the conventions table): one note deriving
   (ω,k) sign from the recursion + the IFO empirical finding, plus a lossy-reflectance
   sanity check.

---

## Phase 1 — Expose existing polarization physics (FIRST deliverable)

**Binding mechanics (corrected to match the actual toolchain).** New API routines go
in `macos_api_mod.F90` following the `beam_set`/`beam_get` template
(`macos_api_mod.F90:4482-4541`): `SystemCheck()` guard → marshal → drive `SMACOS` or
read `use`d module state → `OK=PASS`. Then:

- **mmacos:** most routines ride the **codegen Path A** — add the routine, re-run
  `python3 gen_mex_wrappers.py` from `mmacos/src/`, done (no `mmacos_mex.F` edits).
  Only `rayfield_get` needs a **hand-written Path B helper** (complex split + many
  outputs; follow `do_complex_field`, add to `HAND_WRITTEN`). Veneers in `+macos/*.m`
  + `Session.m` per the naming/validation conventions in `mmacos/CLAUDE.md`.
- **pymacos:** bind in `pymacos_f2py.f90` (NOT "pymacosf90.pyf" — the binding is an
  f2py Fortran shim) + `macos.py` methods.
- Every routine: `matlab.unittest` test + pymacos pytest. Work lands on `sls-dev`
  (both repos); full `run_mmacos_tests.sh` + pymacos pytest green pre-commit; engine
  rebuild chain per `reference_build_test_workflow` (rm mex to force relink).

**Engine API (`macos_api_mod.F90`):**

- `pol_set(OK, on, ExRe,ExIm, EyRe,EyIm)` / `pol_get(OK, on, ...)` — drive the
  `POLARIZATION`/`NOPOLARIZATION` logic and set `Ex0/Ey0`; read back `ifPol`,
  `ifVecDif3`, `Ex0/Ey0`. Mirror `macos_cmd_loop.inc:1072-1097`. Assert `mWF≥3`
  before enabling vector diffraction and return a distinct failure code if not (all
  stock models pass, but a custom param file may not).
- `vecdif_set(OK, on)` — `VECTOR`/`SCALAR`. Enforce the CLI's ordering rule
  (polarization must be on) with a clear error rather than the CLI's silent reset.
- `coat_set(OK, iElt, nLayer, IndRef(:), Extinc(:), Thk(:))` / `coat_get(...)` —
  read/write `EltCoat/IndRefArr/ExtincArr/EltCoatThk` **in physical thickness units**
  (element length units); the waves-normalization is an internal storage detail applied
  on set and inverted on get. `coat_get∘coat_set == identity` to round-off is a
  required test. Must set `modified_rx` (grid-setter dirty-the-trace convention). Must
  also populate the boundary-layer entries (`IndRefArr(0,·)`, `(nCoat+1,·)`) the way
  the parser does — the recursion reads them.
- `rayfield_get(OK, ExRe,ExIm, EyRe,EyIm, EzRe,EzIm, kx,ky,kz, nx,ny,nz, status, N,
  iElt)` — harvest `RayE(3,:)` per ray **plus ray direction cosines, surface normal
  at `iElt`, and the ray-status mask** (buffer-getter pattern; `ray_status_get` is the
  precedent for per-ray `elt_mod` harvests). Phase 2 cannot build a correct Jones
  basis without the geometry, nor honest statistics without the mask; adding them later
  changes a published signature.

**Fix `Ex0Ey0`** — replace the silent swallow at `msmacosio.inc:176` with a real parse
(complex pair ×2), and add the SAVE-writer line (cross-ref `SAVE_KEYWORD_AUDIT.md`) so
source polarization round-trips through a `.in`. Decide with Dave whether to also
revive `PolSrc=` (the commented block) as the on-switch keyword or leave state control
to the API/CLI only.

**Bindings:** `macos.polarization(...)`, `macos.vector_diffraction(...)`,
`macos.coating(...)`, `macos.ray_field(...)` + Session methods + pymacos equivalents.

**Tests/examples:**
1. **Geometry invariance (prerequisite for Phase 2).** On a coating-free system,
   polarization-off and polarization-on traces produce **bit-identical `ray_status`
   and OPD**. Phase 2 assumes ray-by-ray correspondence between traces with different
   input states; this test makes that assumption safe.
2. **Round-trip + wavelength stability.** `coat_get∘coat_set` identity; a stack set at
   λ₁ and evaluated at λ₂ reproduces the analytic Fresnel result for that physical
   stack at λ₂ (should pass already — the semantics are correct; the test pins them).
3. Turn `ifPol` on for `Rx_Cass_FarField.in`, Al coating (tabulated n,κ) on an oblique
   fold; vector-diffraction PSF differs from scalar; coating changes throughput;
   reconcile with the PROPER far-field baseline where polarization-neutral.
4. **Track A smoke:** `ifPol` on the bench_ifo Twyman-Green with a bare-metal Al BS
   coating — fringe visibility drops below the scalar 0.868 by the s/p reflectance
   split; number checked against the single-surface Fresnel analytic at 45°.
5. **`ifPol` off is a no-op.** GMI regression stays 6/6 (its Rx are polarization-off).

---

## Phase 2 — Jones pupil, polarization metrics, and the two application deliverables

The engine has no Jones matrix. Build it the standard PRT way **in the binding layer**
(no engine change beyond Phase 1). This is where the physics can go quietly wrong, so
it is specified in more detail.

### 2a. Jones pupil — and the reference-frame problem

`macos.jones_pupil(srf, 'basis', ...)` (mmacos + pymacos): trace twice with orthogonal
input states ((1,0), (0,1)) via `pol_set`, harvest `RayE` + geometry at the pupil via
`rayfield_get`, assemble the 2×2 Jones per pupil point.

**The basis is not a detail.** The 3×3 PRT matrix is basis-independent; the 2×2 Jones
pupil is not. The naive local ŝ/p̂ basis carries a coordinate singularity that shows up
as **spurious retardance with a tilt/astigmatism-like signature** — at 1e-10-relevant
levels the artifact can exceed the physics. Standard remedy: the **double-pole
(dipole) coordinate system** (Chipman/Lam/Young).

- `basis`: `'double-pole'` (**default**), `'local-sp'`, `'global'`; docstring states
  `'local-sp'` is diagnostic-only, never for budget numbers.
- Return N×N×2×2 complex **plus the vignetting mask**; vignetted points are `NaN`, not
  zero (zero-filling seeds RMS statistics with a ring of perfect nulls).

**Gating acceptance test — unitarity.** Reflectivity forced to unity, no coating: the
Jones pupil is **unitary at every unvignetted point** and reduces to identity (within
the system's geometric rotation) in the correct basis. Cheap; catches basis,
normalization, and sign errors in one check. Regression on every commit touching this
path.

### 2b. Polarization-aberration metrics

`macos.pol_maps(jones)` — pure MATLAB/NumPy:

- **Polar/SVD decomposition** J = H·U per point; diattenuation from singular values,
  retardance from the eigenphases of U.
- **Retardance branch handling:** eigenphase differences are mod 2π; unwrap across the
  pupil and flag ambiguous points instead of silently choosing.
- **Separate the pupil-mean from the spatially-varying part and report both.** Uniform
  retardance is a state change, not a contrast term; only the variation drives the
  floor (track B) or the PSI systematic (track A). One conflated RMS is the likeliest
  way to produce a plausible-looking wrong number.
- **Optional expansion into the standard polarization-aberration terms** (piston/tilt/
  defocus/astig in diattenuation and retardance) — makes two-mirror results directly
  comparable to the published literature; worth it for validation alone.

### 2c. Coronagraph contrast floor (track B deliverable)

`macos.pol_contrast_floor(jones, stokes_in, coronagraph_fn)`:

1. Decompose the Jones pupil into co-polarized (J_xx, J_yy) and cross-polarized
   (J_xy, J_yx) components — the cross terms are the ones scalar DM control cannot
   touch.
2. Propagate each component independently through the coronagraph.
3. Sum intensities **incoherently** across orthogonal output states.
4. Return the floor **broken out by component**, plus its sensitivity to the coating
   parameters.

The answer becomes *"polarization sets the floor at N×10⁻¹¹ and here is how it moves
with coating choice"* — the design-rules line item.

### 2d. IFO polarization metrology (track A deliverable — NEW)

`bench_ifo_pol` example (+ small `+macos/+design` helpers), built on the existing
bench_ifo/bench_ifo_dm machinery (differential PSI, `PSI_MODE`, per-ray pupil map):

1. **Visibility budget.** Coated 45° BS + folds: predicted fringe visibility from the
   s/p amplitude split vs the measured visibility in the simulated interferograms.
2. **PSI polarization systematic.** Run the PSI pipeline with the source polarization
   at 0°, 45°, 90°: the recovered surface shifts by the s/p retardance difference of
   the train. Report the apparent surface error vs input state — this is the number a
   metrology error budget wants (and it is measured with machinery already proven to
   3e-5 pm, so anything it shows is real polarization physics, not pipeline noise).
3. **The AOI trade, closed.** Repeat at BS AOI 45°/30°/22.5° (the `ifo_aoi_sweep`
   scaffold exists): scalar metrology was AOI-neutral; this quantifies the actual
   polarization payoff of smaller AOI — the open question from 2026-07-25.
4. **Bench builder wiring:** coating options on `add_mirror`/`add_bs_reflect`/`plate`
   (emit `Coating=` blocks + `IndRef` ordering handled), `twyman_green('coating',...)`
   pass-through, so polarization studies run through the same one-call rig builder.

**Defer:** populating dead `JmatElt(2,2,mElt)`. Per-element Jones is not what either
application needs; the per-pupil-point export is.

**Tests/examples (beyond the two application deliverables):**
- Unitarity (2a) — gating.
- Folded 2-mirror at ~12° AOI, Al: diattenuation scales with AOI, **vanishes at normal
  incidence**; RMS retardance vs analytic single-surface Fresnel.
- **Rotational-symmetry invariant.** On-axis rotationally-symmetric system → the
  characteristic 2θ symmetry in diattenuation/retardance maps. Strong, cheap,
  basis-sensitive; a broken frame breaks the symmetry.
- Basis comparison `'double-pole'` vs `'local-sp'` on the same system, artifact
  magnitude documented in the worked example.

---

## Phase 3 — New polarizing elements (engine)

**Gated on Phase 0 finding 1 for the VVC only** — the polarizer and waveplate serve
track A through the ordinary trace path and do not wait on vector-propagation
closure.

Each element needs: EltID dispatch in `tracesub.F`/`propsub.F`/(as needed) `srtrace.F`,
a surface routine in `elemsub.F` transforming `Eout` from `Evec`, keywords in
`msmacosio.inc` + defaults in `iosub.inc`, SAVE writer, Phase-1-style API/binding
exposure — and, for track A, a Bench `add_*` emitter.

- **Ideal linear polarizer** — finish `RfPolarizerElt(14)`/`TrPolarizerElt(15)`.
  Keyword for the transmission-axis vector; `Eout` = projection onto the axis, reflect
  variant sends the rejected component. Model on the s/p projection in `Reflector`
  (`elemsub.F:385-428`). **Document:** sequential tracing yields one output port per
  run — a **PBS requires two traces** (the Twyman-Green driver already runs per-arm
  traces, so this fits the existing example structure exactly).
- **Waveplate / retarder** — new element type (extend `mEltTypes`): settable retardance
  + fast-axis orientation; applies the retarder Jones in the local transverse basis.
  **Documented as a thin, non-ray-splitting idealization** — no o/e walk-off. Also the
  primitive for bounding **stress birefringence** in transmissive elements (e.g. the
  wedged Fresnel plate upstream of DM1).
- **Vector vortex (VVC)** — charge-l geometric-phase mask coupling circular states
  with exp(±i·l·θ). **Retardance is a first-class parameter, not idealized away:**

  ```
  J(θ) = cos(δ/2)·I − i·sin(δ/2)·(n̂(θ)·σ),  fast axis at l·θ/2
  ```

  The identity term is the non-vortex leakage, amplitude cos(δ/2); with δ = π + ε,
  leakage intensity = sin²(ε/2) ≈ (ε/2)² falls out automatically. Make δ spatially
  variable and wavelength-dependent so vendor retardance-vs-λ curves load directly and
  broadband null depth is computed honestly.

**Tests/examples:** crossed polarizers → ~0 transmission; QWP converts linear →
circular (Stokes check); **polarization phase-shifting Twyman-Green** (PBS + QWPs —
the classic instrument, now buildable; recovered surface must match the PZT-stepping
result on the same rig); VVC null depth vs ε against the analytic (ε/2)², with and
without a circular analyzer built from polarizer + waveplate.

---

## Phase 4 — Spatially-variable coatings

Extend per-element uniform coatings to spatially varying, modeled on the `AmplMat`
ray→pixel projection (`elemsub.F:3348-3375`, `AmplSrfdx`/`nAmplMat`) — **with two
deliberate departures**:

**Departure 1 — representation.** Per-pixel full layer stacks run to gigabytes and
force a recursion per ray per surface. The physically-motivated variation is nearly
always smooth-and-low-order (radial grade) or **piecewise constant** (segment-to-
segment witness-sample variation across a segmented primary — the actual coronagraph
case). So: **primary form = per-pixel zone-index map + small coating table**; escape
hatch = per-pixel stacks, off by default; optionally precompute per-zone r_p/r_s at
the working wavelength to skip the recursion in the piecewise case.

**Departure 2 — interpolation.** Do **not** inherit `AmplMat`'s nearest-neighbor
lookup: quantization steps in a smooth phase map generate speckles exactly where we
can't afford them. Bilinear minimum; grid-resolution sensitivity an explicit test.

Mechanics otherwise parallel AmplMat: `CoatMapdx`, `ifCoatMap`/`iEltToCoatSrf`,
`CoatFile` keyword parallel to `AmplFile`/`AmplInit` (`surfsub.F:3278-3335`), SAVE
writer, API/bindings. In `Reflector`/`Refractor`, when `ifCoatMap(iElt)` is set, build
`nb_arr/kxb_arr` per-ray from the map — **the recursion is unchanged, only its inputs
become spatially varying.**

**Tests/examples:** radial grade → Jones pupil varies as designed; per-segment coating
→ Jones pupil shows the segment pattern; grid-resolution sweep shows no speckle-floor
dependence on map sampling.

---

## Validation ladder (polarization ground truth)

Reconciling with PROPER only proves the scalar path is intact. In increasing effort:

1. **Analytic single-surface Fresnel Jones** — exact, available now; the unit-level
   check (both tracks).
2. **Unitarity** of a lossless Jones pupil (2a) — the single most diagnostic test.
3. **Rotational-symmetry 2θ invariant** — cheap, basis-sensitive.
4. **IFO self-measurement cross-check (NEW).** Put a known retarder/coating in the
   test arm of the simulated Twyman-Green and *measure* its retardance with the PSI
   pipeline (proven to 3e-5 pm); compare against the same coating's Fresnel analytic
   and against the Phase-2 Jones-pupil decomposition. Three independent paths to one
   number, entirely in-house.
5. **Published two-mirror polarization-aberration results** — compare form and
   magnitude of the diattenuation/retardance patterns (the 2b expansion exists to make
   this direct).
6. **One cross-check against a commercial PRT code (Polaris-M) if a JPL license is
   reachable** — the only item here that satisfies an outside reviewer on its own.

**Per-phase mechanics:** `run_mmacos_tests.sh fast` between edits; full suite +
pymacos pytest pre-commit; every `matlab -batch` ends `exit(0)`. Build both compilers
(`makems.sh release` + `makems.sh release gfortran`) — new fixed-form code passes
gfortran's stricter checks (LOGICAL `.eqv.`, ≤72-col, cpp `//` gotcha per
`mmacos/CLAUDE.md`). GMI regression stays 6/6 (its Rx are polarization-off; Phase-1
gating must be a no-op there).

---

## Files touched (representative, per phase)

- **Engine API:** `macos_f90/macos_api_mod.F90` (all phases).
- **Engine physics:** `macos_f90/elemsub.F` (Ph3 surface routines, Ph4 map lookup),
  `macos_f90/tracesub.F` + `propsub.F` (Ph3 dispatch, Ph4 per-ray inputs; possibly
  Ph3a vector legs), `macos_f90/elt_mod.F` (element constants + arrays),
  `macos_f90/surfsub.F` (Ph4 CoatFile loader), `macos_f90/srtrace.F` (Ph0/1 `ifPol`
  parameter, if in scope).
- **Parse/IO:** `macos_f90/msmacosio.inc` (Ex0Ey0 + Ph3/4 keywords), `macos_f90/
  iosub.inc` (defaults + SAVE writers; note the existing coating-SAVE audit comments),
  `macos_f90/macosio.F` (interactive MOD), `macos_f90/SAVE_KEYWORD_AUDIT.md`.
- **mmacos:** `src/gen_mex_wrappers.py` rerun (+ one Path-B helper in `mmacos_mex.F`),
  `src/+macos/*.m` + `Session.m`, `tests/t*.m`; `src/+macos/+design/Bench.m` +
  `twyman_green.m` (coating emission, Ph3 `add_polarizer`/`add_waveplate`).
- **pymacos:** `src/cmake/source/pymacos_f2py.f90` + `src/pymacos/macos.py` + tests.
- **Docs:** `macos_f90/CLAUDE.md` (new engine polarization/coating section — home of
  the conventions table), `mmacos/CLAUDE.md` / `pymacos/CLAUDE.md` cross-refs, a
  worked polarization example under the coronagraph design docs, and the bench_ifo_pol
  example README.

## Reuse (do not reinvent)

- `beam_set`/`beam_get` — API setter/getter template; `gen_mex_wrappers.py` Path A for
  everything mechanically regular.
- `cfield_get`/`cfield_apodize` — split-real/imag buffer + WFElt-multiply patterns.
- `ray_status_get` — per-ray `elt_mod` harvest precedent (for `rayfield_get`).
- `AmplMat`/`AmplInit`/`ChkRayTrans` — Phase 4 model (with the two departures).
- Existing s/p projection in `Reflector` — Phase 3 Jones-apply model.
- `modified_rx` dirty-the-trace convention for all setters.
- bench_ifo/bench_ifo_dm PSI machinery (proven to 3e-5 pm) — track A test vehicle.

---

## Open items to confirm with Dave

- **Reference-frame default** — double-pole proposed; confirm, and whether any JPL
  coronagraph polarization results we'll be compared against used a different frame.
- **Coating subsystem unification** — unify the two coating models or extend only the
  polarization-path one? (Phase 0 finding 5 informs.)
- **`PolSrc=`/`Ex0Ey0=` Rx semantics** — should a `.in` file be able to turn
  polarization ON (revive `PolSrc`), or only carry the source state (`Ex0Ey0`) with
  on/off left to API/CLI?
- **Phase 3 vs Phase 4 priority after 1+2** — PBS/QWP + VVC trades (Ph3) vs segmented
  coating non-uniformity (Ph4). Note track A gives Phase 3 a customer immediately
  (polarization phase-shifting IFO).
- **Scope of Phase 3a** (vector propagation through the full coronagraph chain),
  pending the Phase 0 audit — the plan's largest schedule unknown.
