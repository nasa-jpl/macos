# Design Layer — Plan v3

A high-level design surface on top of MACOS. The user expresses
*design intent* (topology, free parameters, merit function); a builder
generates a concrete MACOS prescription directly (no CodeV step);
MACOS handles ray trace, diffraction, and inner-loop alignment; an
outer optimizer (MATLAB `fmincon`) explores the design space.

Revised 2026-06 (v3) from the original draft + two design reviews.
Companion to `PLAN.md`.  The analysis surface is **front-end-agnostic**
(§1.0): the expected dominant user *imports* a CodeV/Zemax-converted
Rx and runs sensitivities (Sprint 2A-i), and the de-novo *builder*
(telescope first, Sprint 2A-ii; spectrograph third) populates the same
surface.  The evaluation core is proven coronagraph-first on the
existing Rx corpus (Sprint 1).
**MATLAB-first** — the JPL user base is MATLAB-heavy.  The Python
port stays cheap by construction (see §3 state-as-data rule).

> **Related plans.**  This file owns the MATLAB **design layer**
> (builder, importer, `vary`/`evaluate`/`sensitivities`, compensators,
> metrology orchestration).  **Engine + wrapper** work (Fortran,
> `macos_api_mod`, pymacos/mmacos command surface) is owned by
> `PLAN.md`.  Ownership rule: an item has exactly one home and one
> checkbox; this plan **cross-references `PLAN.md` by section** for
> engine dependencies and never duplicates their checkboxes.  The
> §9.1 queue lists engine asks *motivated by* the design layer, each
> pointing at its `PLAN.md` owner where one exists.

---

## 0. Conventions for Claude Code (CC) sessions

Notes addressed to CC appear inline as blockquotes:

> **CC:** like this.  They carry implementation guidance, guardrails,
> and pointers that aren't derivable from the surrounding text.

Standing rules for any CC session working this plan:

> **CC:** Read `~/dev/macos/CLAUDE.md` and
> `~/dev/MACOS_resources/mmacos/CLAUDE.md` before touching code.
> They are working-memory cheatsheets of gotchas (fixed-form Fortran
> conventions, cpp `//` trap, batch-mode `exit(0)` rule, mixed
> tabs/spaces — do not reformat, codegen workflow).  This plan does
> not repeat them.

> **CC:** Standing regression rule (PLAN.md §5.4): every new wrapper,
> helper, or design-layer function lands WITH a matlab.unittest test,
> even if the motivating task didn't need one.  Dev loop:
> `./run_mmacos_tests.sh fast` (~10 s) between edits; full suite
> pre-commit.  End every `matlab -batch` invocation with `exit(0)`.

> **CC:** New Fortran follows macos/CLAUDE.md conventions: IMPLICIT
> NONE, sequential DO only, fixed-form cols 7–72, `!-->` / `!<--`
> markers around added blocks.  New MATLAB follows the `+macos/`
> conventions: `arguments` validation, struct returns, SI metres at
> the user surface, split `get_X`/`set_X`.

> **CC:** All work lands on `sls-dev` (macos) + companion `sls-dev`
> (MACOS_resources).  Tag sprint completions (`design-sprint-N`).

> **CC:** Each sprint tag ships with a **runnable worked-example
> script** (e.g. `example_sensitivities_from_rx.m` at 2A-i) — that is
> what turns an intermediate stage from "tagged" into "adoptable."
> The example runs end-to-end against a committed fixture and ends in
> `exit(0)`; it doubles as the sprint's smoke demo.

---

## 1. Architecture: three tiers, two loops

### 1.0 Two front-ends, one analysis core

The analysis surface — `vary`, `evaluate`, `sensitivities`,
compensators, `optimize`, the ray-loss guard, parallel FD — is
**front-end-agnostic**.  Two front-ends populate it:

- **Import** (expected dominant entry point): `macos.design.System
  .from_rx(path)` loads the Rx via SMACOS and reads element
  parameters back through the existing getter surface (engine
  readback — there is NO fixed-format text parser in MATLAB).  The
  real-world user arrives with a CodeV/Zemax-converted prescription
  and wants sensitivities, not a de-novo layout.  For imported
  systems, design vars map to in-session perturbations and element
  setters — no re-emit, no reload per outer step.
- **Builder** (de-novo): the family / derivation path (§§4–7).
  Emits an `.in`; reloads on geometry change.

Converters (CodeV script, future Zemax) are **importers** that feed
`from_rx`; `describe()` / `check_clipping()` / `ValidatePrescription`
become the diagnostic surface those "crude" conversions never had
(see §10 Made #10).

Everything downstream of the front-end — §§1.1–1.3 mechanics, the
merit/optimizer/sensitivity machinery — is shared.  The split below
(§1.1) gains a third caching mode for the import path.

**Tier 1 — MATLAB design layer** (`+macos/+design/` in
MACOS_resources/mmacos).  Builder API + outer optimizer.  Owns design
intent: topology, free vs. fixed parameters, merit, bandwidth, FoV.

**Tier 2 — Prescription artifact.**  `.in` file emitted by the
builder whenever *geometry* changes.  Language-neutral handoff;
readable and reproducible on its own; the complete record of every
derived quantity.

**Tier 3 — MACOS execution.**  mmacos wrappers over `macos_api_mod`.

Two nested loops:

- **Inner (MACOS, existing strength).**  Fixed topology + dimensions;
  aligns DMs / rigid-body DOFs against an in-instrument target via
  SLSQP.  Already in place (`calib_run`, `calib_set_var_elt`,
  `calib_set_target`, Phase 1a–1c).
- **Outer (MATLAB `fmincon`).**  Topology fixed; *dimensions free*
  (focal lengths, distances, conics, asphere/freeform coefficients,
  mask radius, Lyot undersize, despace…).  Each step regenerates
  geometry via the builder, runs the inner loop if applicable,
  evaluates the band-and-FoV-averaged merit.

### 1.1 Per-evaluation state vs. geometry — the caching split

| Changes per… | Examples | Mechanism |
|---|---|---|
| outer step (builder geometry) | spacings, Kr/Kc, coefficients, mask radii | builder re-derives → re-emit `.in` → `load_rx` (once per step) |
| outer step (imported geometry) | despace, tilt, conic, Zernike on an imported Rx | perturbation / element setters in-session — **no re-emit, no reload** |
| inner trace (state) | wavelength, chief-ray direction (field point) | setters only: `set_src_wvl`, `set_src_fov('src_dir',…)` — never re-emitted |

```
┌─ outer (fmincon over normalized design params) ───────────────┐
│  x_phys = unnormalize(x)                                      │
│  spec   = builder.update(spec, x_phys)                        │
│  rx     = builder.build(spec)        % derive + validate + emit│
│  m.load_rx(rx)                       % ONCE per outer step    │
│  merit = 0                                                    │
│  for λ in band:                                               │
│    m.set_src_wvl(λ)                  % + dispersive-elt update│
│    for f in field_points:                                     │
│      m.set_src_fov('src_dir', f)                              │
│      [inner CALIB if in-instrument DOFs]                      │
│      m.trace(); assert_no_ray_loss()                          │
│      merit += w(λ,f) * user_merit(read_outputs(), λ, f)       │
│  return merit  (penalized if rays lost / clipping)            │
└───────────────────────────────────────────────────────────────┘
```

> **CC:** The private `evaluate_(spec, x)` function is the ONLY path
> to the engine during `optimize()`.  It owns the call sequence —
> trace-dependent commands (`opd`, `intensity`, `complex_field`)
> silently return zeros without a prior `trace_rays`, and
> `modified_rx` wipes trace state.  No user-level hand-sequencing.

### 1.2 Session model

- One Fortran session per MATLAB process.  The design class takes
  `model_size` at construction, calls `init` once, and errors on any
  change attempt within a study.  (The model_size heap-corruption
  bug is fixed engine-side, but one-size-per-study remains the
  policy: geometric studies and 512-class diffraction studies have
  different sampling needs and should not interleave.)
- Interactive `macos.*` calls during a running `optimize()` corrupt
  the evaluation — document, don't guard.
- **Parallelism:** `fmincon` `'UseParallel',true` distributes FD
  gradient probes over `parfor` workers.  Each worker process loads
  its own mex → its own libsmacos state → no shared-state hazard.
  `evaluate_` must be worker-safe: takes the Rx *path* + parameters,
  performs its own load, relies on no desktop-session state.

> **CC:** Design for worker-safety from Sprint 2A-i even if
> UseParallel ships off by default.  Each worker pays one cold
> init+load; wins when traces are slow (diffraction evals), not for
> 10 ms geometric traces.

### 1.3 Evaluation-pipeline hardening (review findings — non-optional)

These are where outer-loop optimizations live or die:

1. **Emission precision.**  The Rx emitter writes ALL numeric
   parameters at full double precision (`ES24.16E3`-class).
   fmincon's FD step (~1e-8 relative) must survive the text
   round-trip or gradients are quantization noise.
2. **Variable normalization.**  Design vars span mm, mrad,
   dimensionless conics, 1e-12 asphere coefficients.
   `vary(...)` bounds are used to normalize internally to
   [0,1]; `builder.update` unnormalizes.  fmincon never sees raw
   units.
3. **Inner-loop determinism.**  Inner CALIB runs with fixed
   iteration budget + tolerance much tighter than the outer FD
   step's merit effect; warm-start the inner solution from the
   previous outer iterate (cost AND smoothness).  If gradient noise
   persists anyway, `patternsearch` is the fallback, not tighter
   SQP knobs.  For multimodal landscapes — coronagraph
   mask-parameter merits can have multiple basins where two-mirror
   conics do not — `MultiStart` (wrapping fmincon over a set of
   seeds) is the complementary fallback.
4. **Ray-loss guard (Sprint 2A-i, not Sprint 4).**  RMS-WFE-over-
   surviving-rays *improves* when rays vignette, and nRays changing
   is a step discontinuity.  Every evaluation asserts nRays
   constant; lost rays → penalized merit (or `nonlcon` violation).
   Backend: the per-ray status machinery (RayStatus / nBadRays with
   corrected accounting) — expose a summary via `macos_api_mod` if
   not already wired.
5. **Emitted-Rx validation.**  `.build()` runs the Phase-1
   prescription validator (`ValidatePrescription`) on every emitted
   `.in`.  Emitter bugs surface as clean messages, not fixed-format
   READ crashes mid-optimization.
6. **Bandwidth default.**  Geometric OPD is λ-independent for
   all-reflective systems: default nλ=1 when (all-reflective AND
   geometric merit); the λ loop turns on for refractive elements or
   diffraction merits.  Free 3× on telescope iteration speed.

> **CC:** Items 2, 4 are Sprint 2A-i acceptance criteria (analysis
> core); items 1, 5 are Sprint 2A-ii (emitter).  Write the FD-
> survivability test explicitly: emit Rx at x and x+1e-8·x, diff the
> files, assert the parameter actually changed.

---

## 2. User interaction flow

The intended session shape, end to end.  Stages 1–3 are iterative
and cheap; stage 4 is the expensive loop; stages 5–6 close out.

```matlab
%% Stage 1 — declare intent (seconds; no engine calls)
% Fixed-topology families (Cass/RC/Gregorian/DK) auto-populate M1/M2/FP
% with standard names — no add_* calls needed.  (N-mirror and
% instruments use add_* explicitly; see §7.2.)
t = macos.design.Telescope('family','RC', ...
        'aperture_diameter_mm', 6000, 'primary_fnum', 2.0, ...
        'system_fnum', 20.0, 'BFD_mm', 1000, 'model_size', 256);
t.set_field_points(macos.design.hexgrid(deg2rad(0.01), 1));
t.set_bandwidth(632.8e-9);            % nλ=1; all-reflective default

%% Stage 2 — build + sanity (first engine touch)
rx = t.build();        % derive §4 table → validate → emit .in
t.describe();          % table: every derived value + its source
t.diagram();           % 3-D side view, element names, beam radii
r  = t.check_clipping();   % worst margin per element

%% Stage 3 — single-point evaluation (trust before optimizing)
s = t.evaluate();      % baseline merit + per-field/λ breakdown
% user iterates 1↔3 until the baseline design makes sense

%% Stage 4 — optimize
% vary(elt, param, ...) shares addressing with override(elt, param, value):
% a design var IS an optimizer-driven override (§5.4).
t.vary('M2','despace', 'bounds',[-2 2], 'unit','mm');
t.vary('M2','conic',   'bounds',[-1.2 -0.8]);
t.vary('M2','refocus', 'role','compensator');  % solved inner-loop, not by fmincon
t.set_outer_merit('rms_wfe');        % per-(λ,field); loop owns the averaging
[opt, hist] = t.optimize('algorithm','sqp', 'MaxIter',60, ...
                         'UseParallel',false);

%% Stage 5 — inspect
t.plot_history(hist);       % merit vs iter; var trajectories
t.describe();               % updated derived values
s = t.evaluate();           % final per-field/λ breakdown

%% Stage 6 — export / handoff
t.save('rc_design_v1.in');          % frozen, reproducible Rx
t.save_spec('rc_design_v1.mat');    % design spec struct (re-loadable)
% future: t.export_proper(...) for FALCO/PROPER contrast runs
```

Design principles encoded in the flow:

- **No silent derivation.**  `describe()` shows every §4-table value
  with provenance (`user | derived(family) | derived(layout) |
  override | optimized`).  The user can audit before optimizing.
- **`evaluate()` before `optimize()`.**  A single trusted evaluation
  is the unit the optimizer multiplies; the flow makes the user
  exercise it standalone first.
- **The Rx is the deliverable.**  `save()` produces a `.in` that
  stands alone — loadable in interactive macos, shareable, diffable.
- **Spec round-trip.**  `save_spec`/`load_spec` re-creates the
  builder state exactly; the optimizer can resume.

> **CC:** `describe()` and the provenance tags are not cosmetic —
> they are the debugging surface for "why did the builder pick
> that".  Implement provenance as a field carried in the spec struct
> alongside each resolved value.

---

## 3. State-as-data rule (Python port insurance)

The `Telescope` / `Coronagraph` / `Spectrograph` classes hold a
**plain struct** (the design spec); all derivation is pure functions:

```
spec ──resolve──▶ geometry ──emit──▶ Rx text
```

Methods only edit the spec.  Consequences:

- The eventual Python port is a struct translation + function
  transliteration, with **byte-identical `.in` output on the shared
  golden specs** as the cross-language parity criterion.
- The §4 resolution table is testable in isolation (no engine).
- `save_spec` / `load_spec` are trivial.

> **CC:** Resist putting state in class properties beyond the spec
> struct + session handle.  If a method needs scratch state, it
> belongs in the spec or in function locals.

---

## 4. What the user provides (input inventory)

Rule: the user states intent; everything the Rx format needs is
derived; `override('<elt>','<param>',value)` disables the
corresponding derivation for that one quantity.  Minimal Cassegrain
≈ ten numbers + three `add_*` calls.

### 4.1 System level (constructor)

| Input | Status |
|---|---|
| family | required (Cass / RC / Gregorian / DK / TMA / FourMirror / Freeform / Coronagraph / Spectrograph) |
| aperture diameter D | required |
| system f/# (or EFL) | required |
| primary f/# | required for 2-mirror families |
| BFD | required (or `'derive'` where family equations allow) |
| optical axis | default `[0 0 1]` |
| object distance | default infinity (zSource ~1e22) |
| model_size | default by merit type; fixed per study |
| BaseUnits / WaveUnits | builder convention, overridable |

### 4.2 Source (mostly derived)

| Input | Status |
|---|---|
| wavelength list + weights | required (≥1); per-eval state, not re-emitted |
| λ_ref (layout wavelength) | = first λ unless dispersive elements present, then required |
| field points (θx,θy) + weights | required (on-axis = `[0 0]`); per-eval state; helpers `hexgrid` / `cross_grid`; no default grid imposed |
| sampling (GridType, nGridpts) | default by policy (modest circular for geometric; 511-class odd for diffraction); overridable |
| chief ray start point / direction | derived (standoff upstream of M1, along axis) |
| source beam aperture / obscuration | derived: fills M1 + margin incl. chief-ray field walk |
| flux, IndRef, polarization | defaults 1 / vacuum / off |

### 4.3 Per element (ordered topology — the user's main job)

`add_mirror`, `add_lens`, `add_grating`, `add_OAP`, `add_DM`,
`add_mask`, `add_lyot`, `add_focal_plane`, `add_detector` — each with
a name, in propagation order.

| Input | Status |
|---|---|
| spacing_after | required for N-mirror (one `'derive'` slot allowed); derived for 2-mirror families |
| surface type | default `conic`; `asphere_monomial` / `freeform_zernike` / `grid` with order or mode list |
| per-mirror f/# | default family equations or equal-power split; overridable |
| aperture shape | default circular auto-sized; override hex / annular (Cass M1 hole) / segmented (SegMirMaker hook) / polygon |
| off-axis distance (OAPs / unobscured) | explicit value in MVP; clearance-derived later |
| lens: material + form | required (glass name from catalog; curvatures or focal length + bending; thickness) |
| grating: ruling density, order, orientation | required |
| mask params (FPM radius, Lyot undersize) | required for those elements (often design vars) |
| DM actuator / mode count | required |
| obscurations (M2 shadow, spider) | shadow derived for on-axis families; spider opt-in (legs, width) |
| Vpt, Rpt, psi | ALWAYS derived (§5 table) |
| Kr, Kc | derived per family; user-facing only as design-var deltas |
| `Element=`, SrfType, PropType, IndRef, ZernType | ALWAYS derived — institutional knowledge lives here |

### 4.4 Evaluation & optimization

| Input | Status |
|---|---|
| design vars + bounds (+ natural units) | required for `optimize()`; declared via `vary(elt, param, …)` (shares addressing with `override`) |
| outer merit | per-(λ,field) built-in (`rms_wfe` / `rms_spot` / `contrast`); callback override; averaging owned by the λ×field loop + §4.2 weights |
| inner target (DarkZone params) | required for coronagraph; absent for telescope |
| constraints (clearance margins, ray-preservation) | defaults ON; thresholds overridable |
| optimizer options (algorithm, MaxIter, UseParallel) | defaults |

---

## 5. Parameter resolution — Aperture / Vpt / Rpt / psi / Kr / Kc

User states intent; the builder derives the rest closed-form where it
can; the outer optimizer refines what closed-form can't pin down
(TMA / 4-mirror / freeform).  Anything derived is overridable.

### 5.1 Common to all families

| Parameter | Source | First-cut value |
|---|---|---|
| **M1 Vpt** | layout origin | `[0 0 0]` (lab frame) |
| **M1 Rpt** | layout | coincident with M1 Vpt |
| **M1 psi** | optical axis | `[0 0 1]` |
| **M1 Aperture** | user D | circular `D` (override → hex, segmented, polygon) |
| **Ck Vpt** (k > 1) | layout: cumulative spacings along chief ray | `Ck_Vpt = C(k−1)_Vpt + spacing_after(C(k−1)) · dir(k−1→k)` where `dir` comes from each component's exit-direction rule (§6.2) |
| **Ck Rpt** | layout | coincident with Ck Vpt |
| **Ck psi** | layout | mirror: `−propagation_direction_at_Ck`; lens surface: `+propagation_direction`; grating: per mount geometry |
| **Ck Aperture** | beam footprint + field walk + margin | `D_k = 2 · (r_marginal(z_k) + h_chief,maxfield(z_k)) · margin` (default margin 1.05) |
| **FP Vpt** | layout | along exit chief ray at `spacing_after(C_last)` = BFD |
| **FP psi** | propagation direction | exit-axis sign at FP |
| **FP Aperture** | user or derived | detector size, or `2·f_sys·tan(max_field)` (+ dispersion spread for spectrographs) |

> **CC:** The aperture rule now includes the chief-ray walk term —
> the original `2·r_marginal·1.05` clips wide-FoV TMAs at the field
> edges by construction.  Source-block derivation (standoff, beam
> sizing, grid policy) is part of this table's scope: it is the
> first thing the emitter writes and must be unit-tested like the
> rest.

### 5.2 Closed-form Kr / Kc by family

| Family | Pinning strategy | Kr first-cut | Kc first-cut |
|---|---|---|---|
| **Classical Cass** | full closed-form from `D`, `f_M1`, `f_sys`, `BFD` | `Kr_M1 = 2·f_M1`; `Kr_M2 = 2·f_M2`, `f_M2 = p·q/(p+q)`, `p = d_M2 − f_M1`, `q = d_M2 + BFD` | `Kc_M1 = −1`; `Kc_M2 = −((m+1)/(m−1))²` |
| **Ritchey-Chrétien** | Cass layout; Kc for zero spherical + coma | as Cass | aplanatic forms **including the back-distance (β = BFD/f_M1) term** — see CC note |
| **Classical Gregorian** | requires `d_M2 > f_M1` (intermediate focus) | as Cass with sign of `f_M2` flipped | `Kc_M1 = −1`; `Kc_M2 = −((m−1)/(m+1))²` (prolate ellipsoid, −1 < Kc < 0) |
| **Dall-Kirkham** | spherical secondary | as Cass | `Kc_M2 = 0`; `Kc_M1` from the aplanatic-spherical condition **with β term** — see CC note |
| **TMA (3-mirror)** | layout from user spacings; equal-power split default `f_i ≈ N·f_sys`; `override('Mk','fnum',…)` to repower | `Kr_k = 2·f_k` | `Kc_k = 0`; outer loop drives Kc + asphere coeffs |
| **4-Mirror** | as TMA | `Kr_k = 2·f_k` | `Kc_k = 0`; outer loop drives Kc + asphere + freeform |
| **Freeform** | as TMA/4-mirror; any mirror may carry monomial / Zernike / grid terms | `Kr_k = 2·f_k` (base sphere) | `Kc_k = 0` + freeform coefficients as design vars |

`m = f_sys / f_M1` (magnitude convention) for the 2-mirror families.

> **CC:** The original draft's Gregorian secondary
> (`+((m+1)/(m−1))²`, labeled prolate) was wrong on sign and shape —
> corrected above.  The original RC (`Kc_M1 = −1 − 2/(m²(m+1))` etc.)
> and DK (`Kc_M1 = −1 + 2/m³`) forms drop the back-distance
> dependence; implement the exact β-dependent forms (Schroeder,
> *Astronomical Optics*, 2-mirror aberration chapter) and DO NOT
> trust any transcription — including this one — without the §5.3
> raytrace test.

### 5.3 Family validation tests (Sprint 2A-ii acceptance)

Per family: build the Rx, trace on-axis and at a small field, Zernike-
decompose the WFE, assert:

- third-order spherical ≈ 0 (all families, on-axis),
- coma ≈ 0 (RC — aplanat check),
- known nonzero residuals match textbook scaling (classical Cass coma,
  DK field coma) at loose tolerance.

This catches every formula error in §5.2 at once AND pins the builder
against MACOS sign conventions (the classic 2-mirror failure mode).

### 5.4 Spacings, sanity checks, overrides

- `'spacing_after', value` is the default; a single `'derive'` slot is
  solved from the f_sys / BFD constraint (multiple `derive` → error;
  zero `derive` + inconsistent f_sys → warning).  For 2-mirror
  families both spacings can be derived.
- `.build()` sanity checks: monotone non-negative spacing cascade;
  per-element clearance (full `check_clipping()` report deferred to
  Sprint 4, but pass/fail margin check runs from Sprint 2A-ii); BFD > 0;
  Cass-class FP reachable through the M1 hole; Gregorian intermediate
  focus exists (`d_M2 > f_M1`); N-mirror power sum consistent with
  f_sys; emitted Rx passes `ValidatePrescription`.
- Override rule unchanged: any table entry is overridable by
  `t.override('<elt>','<param>',value)`; an override disables that
  one derivation, everything else proceeds.
- **`vary('<elt>','<param>', …)` is the same addressing as
  `override`** — a design var IS an optimizer-driven override.
  Declaring `vary` on a derived quantity takes precedence over the
  derivation for the duration of `optimize()`.
- **Compensators.**  `vary('<elt>','<param>', 'role','compensator')`
  marks a DOF solved in the *inner* loop per outer iterate
  (warm-started, §1.3.3) rather than driven by fmincon.  Rationale:
  a despace design var without a refocus compensator lets defocus
  dominate the merit and looks like a broken tool; the same concept
  is what tolerancing needs (compensated-residual budgets).

---

## 6. Component model

### 6.1 Component ≠ Rx element

The builder's topology unit is the **component**; MACOS Rx elements
are the expansion product.  Mirror → 1 element; singlet → 2 elements
(`L1_S1`, `L1_S2`); cemented doublet → 3; grating → 1 with grating
keywords; DM → 1; mask / Lyot → 1 `Element= Obscuring`.

### 6.2 Component contract (what every class supplies)

1. **Paraxial behavior** — power + principal-plane locations, so the
   first-order layout chain works generically.  Mirrors: principal
   planes at the vertex.  Thick lens: lensmaker + principal-plane
   offsets; internal thickness owned by the component
   (`spacing_after` is vertex-to-vertex along the chief ray).
2. **Exit-direction rule** — `exit_dir = component.direction(in_dir, λ_ref)`.
   Mirror: reflect.  Lens: pass.  Grating: grating equation at the
   design order and λ_ref.  This generalizes the original
   "sign alternates after each reflection" Vpt-chaining rule.
3. **Rx expansion** — emit N elements with ALL derived keywords
   (`Element=`, SrfType, IndRef/glass name, ZernType, PropType…).
4. **Perturbation grouping** — design vars and perturbations declared
   on a component map to a rigid group move of its elements
   (CPERTURB_GRP / `prb_elt_grp`, already wired in mmacos), never to
   loose per-surface moves.
5. **Aperture derivation** — per-surface footprint via §5.1 rule.

> **CC:** Implement the contract in Sprint 2A-ii with Mirror, Mask,
> FocalPlane only (all trivial cases) — the point is the interface,
> so Lens/Grating (Sprint 2C) slot in without re-architecting the
> Vpt chain or the design-var mapper.

### 6.3 Institutional-knowledge guards baked into expansion

- **Masks emit `Element= Obscuring`, unconditionally.**  A mask on a
  `Reference` element clips rays but leaves the diffraction grid
  untouched (the Phase-5 silent failure).  `.build()` hard-errors if
  an override makes a masking element any other type.  Model:
  `CoroExample.in` Elt 6.
- **Explicit ZernType, restricted to the dispatch-handled set.**
  ZernTypeL = 0 silently no-ops the Zernike-apply path; Noll(10) /
  ExtFringe(11) parse but silently no-op in the propsub/srtrace
  chains.  The builder only emits handled types.

> **CC:** the ELSE-with-error on the propsub.F / srtrace.F ZernTypeL
> dispatch chains **already landed** (PLAN.md §0, commit b2c2eb8;
> §9.1 Q3) — Noll is handled, unhandled types error loudly.  No
> Fortran work remains even though the design layer is about to
> become the first heavy user of SrfType 14 + FFZernType.

### 6.4 Glass catalog (refractive components)

MACOS resolves glass → n(λ) via a separate text-file catalog
(user-extensible).  Design-layer consequences:

- `add_lens(..., 'material','N-BK7')` emits the glass *name*; no
  dispersion code in MATLAB.
- **Catalog expansion work item:** write an `.agf → MACOS catalog`
  converter (Python, offline) rather than hand-curating — Schott /
  Ohara / CDGM / Hoya / IR vendors all publish .agf.  Shipped
  "usual set" is a generated artifact: Schott main line, fused
  silica, CaF₂ / BaF₂ / MgF₂, and the IR set (Si, Ge, ZnSe, ZnS,
  sapphire).  First: check which dispersion formulas the MACOS
  catalog parser accepts; refit foreign formulas if needed (refit
  residual ~1e-7 is fine).  dn/dT out of scope.
- Validation: per-glass n(λ) vs published values to ~1e-6; one traced
  doublet vs CodeV.  Side benefit: the CodeV converter can start
  emitting glass names instead of baked numeric indices.
- **Pre-flight check (Sprint 0):** when does the engine re-resolve
  glass → n(λ)?  If at trace / wavelength-MOD time, `set_src_wvl`
  suffices; if at Rx-load only, the eval loop needs a per-λ index
  push or reload.  One glass element, two wavelengths in-session,
  compare against fresh loads.

### 6.5 Dispersive systems make layout λ-dependent

Spectrograph geometry closes at `λ_ref` (chief ray follows the
grating equation at the design order); other wavelengths walk across
the detector — that's the instrument.  Hence: `λ_ref` in the spec
distinct from evaluation bandwidth (§4.2); detector sizing adds the
dispersion-direction λ-spread to the field-walk margin; spot-based
merits over λ × field replace pupil-referenced WFE as the natural
objective.  Bandwidth helper provides named line sets (`'FdC'`, …)
for achromat-type merits.

HOEs: parked behind the same component interface as an advanced
pass-through (construction-point keywords supplied raw, no
derivation) until a real design needs more.

### 6.6 Metrology (three tiers)

Optical metrology (edge sensors, laser gauges, fiducial truss) enters
the design layer as a **measurement-space Jacobian** — `dMeasurement /
dDOF` — that plugs into `sensitivities()` (Sprint 2A-i)
alongside the merit-space Jacobian.  Three tiers, increasing trust
required:

1. **Interface (now).**  Measurement space in `sensitivities()` is a
   Jacobian *contract*: any backend producing `dMeasurement/dDOF`
   plugs in.  Backend #1 is the SegMirMaker `Hx.m` edge-sensor model
   — MATLAB-native, trusted, no dependence on engine met functions.

2. **Validation (separable — §9.1 Q8; VALIDATION/CHARACTERIZATION,
   not implementation).**  John Lou's MACOS met functions have useful
   outputs but are unvalidated, incomplete, and lack truss-
   optimization support.  Acceptance tests are closed-form geometric
   truths, not regressions against the code's own output:
   - gauge change = projection of relative fiducial motion onto the
     line of sight (analytic equality);
   - perturbation orthogonal to LOS → zero reading (null test);
   - launcher/target swap → reciprocal response;
   - FD vs analytic consistency where a linearized path exists.
   Deliverables: a pass/fail map per met function + an inventory of
   missing capability (that inventory is the spec input for joint-
   side completion work).

   > **CC (guard):** the design layer must NOT consume the engine met
   > functions for budget work before Q8 passes.  "Has produced
   > useful outputs" is exactly the trust level at which a silent
   > sign / frame error survives longest.

3. **Truss optimization (Sprint 6+).**  Vars = fiducial / launcher
   placement + beam topology; merit = observability of the DOF set
   (min singular value / condition number of the restricted met
   Jacobian; worst-case unobservable mode); constraints = launcher
   count, mounting real estate, beam clearance.  First cut is **pure
   MATLAB kinematics** (fiducial coordinates + LOS projections →
   Jacobian → merit; no engine in the loop), per the logic-starts-
   MATLAB-side principle.  MACOS's role: validate selected
   configurations through the real gauge functions, and trace
   met-beam clearance through structure + optical train.
   Dependency chain: tier-1 interface (above) → Q8 → kinematic
   optimizer → MACOS-backed refinement.

---

## 7. Builder API examples

### 7.1 Classical 2-mirror Cassegrain-class

```matlab
% Fixed-topology family: M1/M2/FP auto-populated by the constructor.
t = macos.design.Telescope( ...
        'aperture_diameter_mm', 6000, 'primary_fnum', 2.0, ...
        'system_fnum', 20.0, 'BFD_mm', 1000, ...
        'optical_axis', [0 0 1], 'family', 'Cassegrain');
t.vary('M2','despace', 'bounds',[-2 2],  'unit','mm');
t.vary('M2','tilt',    'bounds',[-1 1],  'unit','mrad');
t.vary('M2','conic',   'bounds',[-1.2 -0.8]);
```

### 7.2 N-mirror — TMA / 4-mirror / freeform (wide FoV)

```matlab
t = macos.design.Telescope( ...
        'aperture_diameter_mm', 6000, 'system_fnum', 20.0, ...
        'BFD_mm', 500, 'family', 'TMA');
t.add_mirror('M1', 'spacing_after_mm', 6000);
t.add_mirror('M2', 'spacing_after_mm', 3000);
t.add_mirror('M3', 'spacing_after_mm', 'derive');
t.add_focal_plane('FP');
t.set_surface('M1', 'type','conic');
t.set_surface('M2', 'type','asphere_monomial', 'order',6);
t.set_surface('M3', 'type','freeform_zernike', 'modes',4:15);
t.vary('M2','asphere_coeffs', 'bounds',[-1 1]*1e-12);
t.vary('M3','zern_coeffs',    'bounds',[-1 1]*1e-9);
t.set_field_points(macos.design.hexgrid(deg2rad(0.05), 2));
t.set_outer_merit('rms_wfe');        % per-(λ,field); loop owns the averaging
```

### 7.3 Coronagraph

```matlab
c = macos.design.Coronagraph('input_focal_plane', ota_fp_elt, ...
                             'pupil_diameter_mm', 50, 'model_size', 512);
c.add_OAP('OAP1', 'recreates','pupil');     % conjugate solve — the one smart keyword
c.add_DM ('DM1',  'modes', 60);             % modal control (see §8 Sprint 3 scope)
c.add_DM ('DM2',  'spacing_after_mm', 300); % explicit distance in MVP
c.add_OAP('OAP2', 'recreates','focus');
c.add_mask('FPM', 'shape','disk');          % emits Element= Obscuring
c.add_OAP('OAP3', 'mirrors','OAP1');
c.add_lyot('LS1');
c.add_OAP('OAP4', 'recreates','focus');
c.add_detector('SciCam');
c.vary('FPM','radius',     'bounds',[2 8], 'unit','lambdaD');
c.vary('LS1','undersize',  'bounds',[0.7 0.95]);
c.set_bandwidth([550 600 650]*1e-9);
c.set_inner_target('DarkZone', 'inner_lamD',7, 'outer_lamD',10);
c.set_outer_merit('contrast', 'inner_lamD',7, 'outer_lamD',10);  % built-in per-(λ,field)
[opt, hist] = c.optimize();
```

> **CC:** MVP DSL discipline: the conjugate-recreating solve
> (`'recreates','pupil'|'focus'`) is the only smart placement
> keyword — it is essential for pupil/focal-plane fidelity.  Talbot-
> fraction placement etc. are explicit distances plus a helper
> function until they earn DSL status.

### 7.4 Spectrograph (Sprint 2C worked example)

```matlab
s = macos.design.Spectrograph('slit_fnum', 8, 'lambda_ref', 600e-9, ...
                              'band', [450 750]*1e-9, 'order', 1);
s.add_lens   ('COL', 'material','N-BK7', 'focal_length_mm', 200);
s.add_grating('GR',  'lines_per_mm', 600, 'order', 1);
s.add_lens   ('CAM', 'material','N-BK7', 'focal_length_mm', 150);
s.add_detector('DET');
s.set_outer_merit('rms_spot');       % per-(λ,field); loop owns the averaging
```

Exercises in one example: 2-surface expansion, glass-name emission,
λ-dependent exit direction, λ_ref layout, dispersion-aware detector
sizing, chromatic (nλ>1) merit.

---

## 8. Sprint sequence

### Sprint 0 — pre-flight experiments (hours each; do before 2A)

- [x] **Endurance test (Q5).**  **Landed + PASS 2026-06-12.**  Loop
      `load_rx` → `trace` on Rx_Cass_FarField (model_size 128) in one
      session.  Result over 500 iters: **rmsWFE bit-identical** (unique
      = 1, max|diff| = 0 — zero state drift across repeated load); RSS
      warms up through ~iter 175, takes one ~14 MB arena event, then
      **+0 kB / 25-iter — flat steady state** (no linear leak; the
      initial naive (end−start)/N "240 kB/iter" was a warm-up
      artifact).  Load-only and trace-only in isolation are each ~flat
      (1.5 / 0.5 kB/iter, allocator noise).  Test: `tEndurance.m`
      (on-demand, `./run_mmacos_tests.sh tEndurance`, ~31 s; not in
      fast/all).  **Foundation certified for the 2A-i evaluate_ loop.**
- [ ] **Glass re-resolution test** (§6.4) — determines the eval
      loop's per-λ step on refractive systems.

(The CALIB-through-diffraction timing probe from the original draft
is absorbed into Sprint 1's E3.)

> **CC:** Each pre-flight lands as a (skippable-if-slow) test or a
> dated note in this file with the measured numbers — they are
> design inputs, not throwaway probes.

### Sprint 1 — coronagraph-driven wrapper parity + macos ↔ MATLAB split

**Vehicle: the existing coronagraph Rx corpus, no builder.**
`Rx_Coro.in` (HCIT-style), the Phase-5 FPM + Lyot configuration
(known baselines: ~3.2e6 peak suppression, ~3e-10 dark-zone contrast
at 7–10 λ/D), and `CoroExample.in`.  Sprint 1 proves the *evaluation*
surface on hand-written prescriptions; the builder (Sprint 2)
later drives the same surface.

Wrapper items:

- [ ] `intensity()` / `complex_field()` exercised on the coronagraph
      corpus specifically — README marks them MVP-wired, but at
      512-class grids on `Rx_Coro.in`, not just `Rx_Cass`.  Verify +
      regression-test.
- [x] **Already present in mmacos:** `set_src_wvl(λ)` (direct state
      setter, no `modify()` round-trip) and `set_src_fov('src_dir',
      DIR, …)` (absolute chief-ray pointing setter — docstring
      specifically calls out "tip source pointing between field
      points without round-tripping through perturb_src", which is
      exactly the inner-loop need).  Verification only.
- [x] Expose a **ray-loss summary** (per-category counts from
      RayStatus) — backend for the §1.3 guard.  **Landed 2026-06-11**
      (Q2): `ray_status_get` in `macos_api_mod` + `get_ray_status(N)`
      mmacos wrapper + tests.
- [ ] Verify CALIB wrappers handle DM-element variables on `Rx_Coro.in`.
- [ ] Tests: load `SegDemo3.in`, set field point, trace, read WFE;
      expected values from pymacos.

**Division-of-labor experiments** — how much lives in macos (Fortran)
vs MATLAB is an open question, resolved empirically here, on the
coronagraph examples, with measured numbers.  Each experiment is a
short MATLAB script + a dated note in this file:

- [x] **E1 — MATLAB-side outer merit.**  Read `intensity(det)`,
      compute annular dark-zone contrast in MATLAB (port the
      `contrast.py` λ/D machinery from pymacos).  Measure per-eval
      wall time at nλ = 3 on the Phase-5 configuration; reproduce the
      ~3e-10 baseline.
      **DONE 2026-06-13.**  Script
      `mmacos/examples/design/coro/E1_darkzone_contrast.m`; ported λ/D
      machinery `{radial_profile, first_airy_null,
      lambda_over_D_pixels, radial_contrast}.m` in the same folder;
      synthetic-Airy unit test `tCoroContrast` (5/5, in SUITE_FAST).
      Config: `Rx_Coro_noLyot.in` (no-mask Strehl ref) +
      `Rx_Coro_FPM.in` (FPM=400 µm + Lyot=14 mm), detector Elt 21,
      model_size 1024 (nGridpts 511), λ=850 nm.
      **Measured numbers:**
      - λ/D = 8.61 px (first Airy null at 10.50 px); no-mask peak 0.3042.
      - **Peak suppression factor = 3.21e6** — matches the documented
        ~3.2-million Phase-5.2 baseline to 3 sig figs (the hard
        cross-check; confirms macos physics + the λ/D port are both
        correct).
      - Dark-zone radial contrast (Strehl-norm to no-mask peak):
        ring-**mean** over 7–10 λ/D = 1.5e-9, dipping to **7.9e-10 @
        9.9 λ/D**; the e-10 floor noted as "~3e-10" is the band
        *minimum* (finer sampling between bins), not the annular-ring
        mean.  Bright outer-ring artefact at ~15 λ/D present as
        documented (1.0e-9, above the 10 λ/D dip).
      - **Per-eval wall time = 2.4 s/eval; nλ=3 merit = 7.24 s total**
        (model 1024, full CPROPAGATE chain to Elt 21 via
        `set_src_wvl` + `intensity(21)`).
      **Takeaway for E2/E3:** a MATLAB-side outer merit is ~2.4 s per
      (λ,eval) at 1024.  An fmincon outer loop doing O(10²–10³)
      function evals × nλ × n_field then costs minutes–hours purely in
      the diffraction trace — so the *merit* staying MATLAB-side is
      fine, but the *inner* DM solve (E2) cannot be a naïve
      MATLAB-driven fmincon over many modes at this per-eval cost.
      That is the bar the Sprint-3 Fortran `DarkZone` target must beat.
- [x] **E2 — DM control entirely in the outer loop.**  DM Zernike /
      Fourier modes as fmincon design vars, no inner loop, MATLAB
      merit from E1.  Find the mode count where cost becomes
      prohibitive — this sets the bar a Fortran DarkZone target must
      beat.
      **DONE 2026-06-13.**  Script
      `mmacos/examples/design/coro/E2_dm_modes_cost.m` on `Rx_Coro_DM.in`
      (DM = Elt 4, a 1024×1024 grid-data Reflector driven via the
      `elt_srf_grid_data` setter — the SAME carrier influence-function /
      measured DM maps will use; smoke-verified the setter moves the
      focal-plane PSF, not a silent no-op).  Fourier "ripple" mode basis,
      `fmincon` (sqp, forward-diff), objective = selectable dark-zone
      metric (`dark_zone_metrics.m`: mean | peak | floor | median |
      energy; default mean).  Full workspace (incl. 1024² images) saved
      to `coro/results/` for resume; pruned keep-last-2 + weekly cron
      `clean_results.sh`.
      **Measured cost** (model 1024, ~3.4 s/trace):

      | nModes | evals (MaxIter=2) | s/eval | s/iter = (nModes+1)×3.41 |
      |---|---|---|---|
      | 3  | 12 | 3.46 | 20.7 |
      | 8  | 27 | 3.41 | 46.1 |
      | 15 | 48 | 3.34 | 80.1 |

      - **FD-gradient cost scales exactly as (nModes+1) traces/iter.**
        Extrapolated prohibitive count (K=30 iters, 1-hr budget):
        **~34 modes single-λ, ~11 at nλ=3** — i.e. a GLOBAL modal basis
        is hopeless for a real DM (hundreds–thousands of actuators) under
        MATLAB-driven FD.
      - **Multiplexed Jacobian — the influence-function payoff.**
        Actuator influence has compact support (< ~3 act spacings), so a
        separable ~6×6 (36-phase) poke coloring measures the FULL
        Jacobian in **~37 traces ≈ 126 s, INDEPENDENT of actuator count**
        (vs 2001 traces ≈ 6826 s naïve for a 2000-act DM).  This is the
        reason to model DMs as influence-function / grid surfaces, and it
        reframes the Sprint-3 bar: **the cost is the EFC linear algebra +
        iteration, not the trace count.**  (Applies only to LOCAL
        actuator DOFs; global Zernike/Fourier modes stay (nMode+1).)
      - **Dark-zone metrics** (flat-DM, full annulus 7–10 λ/D, 11864 px):
        mean 6.0e-5, peak 2.0e-4, **floor 2.9e-10** (= the documented
        "~3e-10", confirming that figure is the *floor*, not the mean),
        energy 0.71.  FD-contrast control made no progress in 2 iters
        (flat gradient) — confirming naïve FD-contrast is the wrong
        inner loop; **E-field/EFC Jacobian control is what's needed.**
      - **Region caveat:** scored over a full 360° annulus, but
        `Rx_Coro_DM` has ONE DM → the fair region is **one-sided**
        (deeper contrast, smaller area).  Cost numbers are
        region-INDEPENDENT so the headline stands; the achievable
        *contrast* would improve one-sided.  `dark_zone_metrics` now
        takes a `'side'`/`'sector'` arg for this (Sprint-3 DarkZone
        target must tie region geometry to DM count).
      Tests: `tCoroContrast` (8/8, pure-math, in SUITE_FAST) pins the
      ported λ/D machinery + `dark_zone_metrics` incl. the one-sided
      region.
- [x] **E3 — CALIB mechanics through diffraction.**  Exercise the
      existing CALIB targets on DM elements of `Rx_Coro.in`; time the
      objective evaluation when it must run the CPROPAGATE chain at
      nGridpts = 511.  (Absorbs the Sprint-0 CALIB-timing pre-flight.)
      Determines whether a Sprint-3 inner loop is seconds or hours.
      **DONE 2026-06-13.**  Script
      `mmacos/examples/design/coro/E3_calib_timing.m` on `Rx_Coro_DM.in`
      (DM Elt 4), model 1024 / nGridpts 511.  Workspace saved to
      `coro/results/`.
      **Key finding — the shipping CALIB targets are all RAY-TRACE
      objectives** (WFE / WFE_ZMODE / BEAM / SPOT / OPL compute from the
      ray trace at the target element; none runs the diffraction FFT
      chain).  Per-eval costs (forced real re-trace via modify()):
      - ray-trace (modify+trace) = **1.03 s**
      - diffraction-propagate (INT only) = **1.07 s**
      - full diffraction objective (trace+INT) = **3.14 s** (ratio 3.1×
        over ray-trace, *not* the 1000× a naïve no-op probe suggests).
      - **Existing CALIB (WFE, 2-DOF DM tip/tilt recovery): 33 s,
        converged** (WFE 6.3e-5 → 2.4e-12).  So the existing ray-trace
        inner loop is **seconds**.
      **Projected diffraction-scoring DarkZone inner loop** (3.14 s/eval,
      K=30 iters):
      - naïve FD: **0.3 h @ 12 DOF, 2.6 h @ 100, 52.4 h @ 2000 actuators**
      - multiplexed/EFC Jacobian (E2 separable poke): **~116 s (37
        traces), independent of actuator count.**
      **Verdict:** a Sprint-3 DarkZone target driven by naïve FD is HOURS
      at real DM scale and must instead use the multiplexed-poke / EFC
      E-field Jacobian (then it's ~2 min).  The Fortran win is avoiding
      the (nDOF+1)-trace FD penalty, not raw per-trace speed.
- [x] **E4 — λ-loop placement.**  Measure the load / trace / propagate
      cost split per wavelength to quantify what a future MACOS-side
      `spectral_run` amortization would actually buy over a MATLAB
      loop of `set_src_wvl` calls.
      **DONE 2026-06-13.**  Script
      `mmacos/examples/design/coro/E4_lambda_loop_cost.m` on
      `Rx_Coro_FPM.in`, model 1024, nλ=5.  Workspace saved to
      `coro/results/`.  Cost split per wavelength:
      - load_rx (one-time, amortized away) = 5.4 s
      - set_src_wvl + modify = **0.001 s** (negligible)
      - ray-trace = **0.62 s**;  diffraction-propagate = **2.77 s**
        → **trace/propagate = 0.22** (propagate dominates).
      - **MATLAB↔mex per-λ overhead = 0.5 ms** (negligible).
      **What `spectral_run` would buy:** for nλ=5, a MATLAB
      `set_src_wvl` loop (re-traces each λ) = 16.9 s vs a `spectral_run`
      that traces ONCE (reflective → λ-independent geometry) and loops
      only the diffraction = 14.5 s — **saves ~2.5 s (15%)**, and that
      15% is the redundant re-trace, NOT mex overhead (which is ~0).
      **Verdict: don't build `spectral_run` for reflective
      coronagraphs** — the MATLAB loop is within 15% and the mex
      round-trip is free.  It only pays off when trace >> propagate
      (huge ray counts / cheap propagation) or for REFRACTIVE systems
      where per-λ glass re-resolution makes the trace dominate (§6.4 /
      Sprint-0 glass test).

**Sprint 1 division-of-labor experiments E1–E4 COMPLETE (2026-06-13).**
**Net architectural resolution: diffraction-based optimization is
driven in the MATLAB layer; NO diffraction metric calcs or
diffraction-scoring inner loop need to be added to the Fortran code.**
What the measurements support:
- The merit/metric calc (contrast from `|field|²`) is trivial → MATLAB
  (E1; metrics selectable in `dark_zone_metrics`).
- The expensive piece — the diffraction propagation — is *already* in
  Fortran (`CPROPAGATE`); MATLAB reaches it via existing `intensity()`
  / `complex_field()`.  Nothing new in Fortran.
- MATLAB↔mex per-call overhead is **0.5 ms** (E4) → no performance case
  for relocating the loop into Fortran.
- The only real cost driver is the Jacobian trace count; the fix is the
  multiplexed separable-poke / EFC **E-field** Jacobian (~37
  `complex_field` evals, actuator-count-independent — E2/E3), which is
  MATLAB orchestrating existing Fortran calls + linear algebra.
- **This revises the earlier "Fortran `DarkZone` CALIB target"
  assumption (Sprint 3):** a Fortran CALIB DarkZone target would still
  use FD gradients → (nDOF+1) penalty → hours.  The right design is a
  **MATLAB EFC controller over `complex_field`**, not new Fortran.
- **One Fortran-side assembly primitive DOES belong in the loop:
  `COMPOSE`.**  Each λ has a different native focal sampling
  (focal dx ∝ λ); macos's `COMPOSE` + `ADD`/`DADD` (intensity) /
  `CADD` (complex amplitude) assembles the per-λ PSFs onto a FIXED
  pixel grid in Fortran — the natural broadband-PSF builder, so MATLAB
  scores the composed broadband PSF instead of resampling each λ by
  hand.  `COMPOSE` + `ADD` is **now wrapped** as `compose()` in pymacos +
  mmacos (2026-06-13, macos `6cd1372` + MACOS_resources `1a3fc13`):
  the incoherent-intensity path is live; coherent `CADD` is still a
  Fortran stub.  Broadband dark-zone scoring can build on it directly.
- Residual Fortran-pull futures (only if they materialise): REFRACTIVE
  systems (per-λ glass re-resolution makes the ray trace dominate —
  E4 `spectral_run` regime), or letting macos skip the redundant
  upstream re-trace per DM poke — both optimizations of *existing*
  Fortran, not new diffraction-metric code.
- Dark-zone region geometry (one-sided for 1 DM, annulus for 2; E2/Dave)
  is a MATLAB-side scoring parameter (`dark_zone_metrics` 'side'/'sector').

Resolution principle: merit and control logic start MATLAB-side
(flexible, zero Fortran risk); a capability moves into macos only
when the measurements show the inner loop needs it at speed (the
DarkZone target, Sprint 3) or when it must participate in propagation
itself (the `Element= Apodizer` work in PLAN.md §2.1).  E1–E4's
numbers make that call, and they scope Sprint 3 before any Fortran
is written.

> **CC:** The wrapper items are almost certainly codegen Path A (add
> to `macos_api_mod`, re-run `gen_mex_wrappers.py`) — see mmacos
> CLAUDE.md "Adding a new command".  Record E1–E4 results as dated
> notes under this section; they are design inputs for Sprint 3, not
> throwaway probes.

### Sprint 2A-i — import / analysis core (lands first)

The front-end-agnostic analysis surface (§1.0), proven on an
**imported** Rx — the deliverable the existing user base actually
wants.  No family math, no emitter required.

> **CORE COMPLETE 2026-06-12** (tag `design-sprint-2A-i`).  Package
> `+macos/+design/System` shipped over four slices on MACOS_resources
> `sls-dev`: from_rx (fc61db6) → sensitivities (22d133a) → vary
> (6be28b7) → evaluate/optimize (9240e55).  Tests tDesignSystem (5),
> tDesignVary (12), tDesignSensitivities (2, bitwise vs standalone
> drivers), tDesignOptimize (4) — all green in SUITE_FAST.  Examples
> in `mmacos/examples/design/`.  End-to-end proof: import e5hex1 →
> 0.1 mm despace (WFE 7.0e-2) → optimize → 3.1e-5 (~2000×), despace
> recovered to −2.6e-5 mm.  **Deferred to follow-ons** (not blockers):
> nested λ×field loop (single-field nλ=1 today); compensator
> inner-loop *solve* (role is declarable; operationalizing it waits on
> the compensator pass); worker-parallel validation; Zernike/surface
> *optimization* (sensitivities cover Zernike; optimize() is rigid-MVP).

- [x] New package `+macos/+design/` under `MACOS_resources/mmacos/`.
- [x] `macos.design.System.from_rx(path)` — load via SMACOS, read
      element parameters back through the existing getter surface
      (engine readback; no MATLAB text parser).
- [x] `vary(elt, param, …)` mapped to in-session perturbations /
      element setters (no re-emit, no reload per outer step —
      §1.1 imported-geometry row).  Name-based DOFs (no 0-based leak),
      local/EltCoord frame, aliases (despace/tilt/decenter), Zernike,
      `'role','compensator'` *declarable*.  Compensator inner-loop
      *solve* + nominal snapshot/restore (§9.1 Q9 → PLAN.md §11.7) are
      follow-ons; `evaluate` restores via `load_rx` today (Q5-certified
      bit-stable).
- [x] `evaluate_` with canonical call sequence, ray-loss guard
      (nRays drop → penalised merit; the Q2 per-category getter is the
      richer backend, wired when needed), internal [0,1] design-var
      normalization.  Worker-safe by construction (own load); parallel
      not yet exercised.
- [x] `sensitivities()` — **builds on the existing Phase 7 `dw_dx`
      channel / supervisor machinery** (PLAN.md §5.4 Phase 7.a/7.b,
      already bitwise mmacos==pymacos): `RigidBodyChannel`,
      `SourceChannel`, `FocalPlaneChannel`, `GroupedRigidBodyChannel`,
      driven by `dw_dx.m` / `dw_dx_multi.m` over a field set.  Maps
      `vary(...)` DOFs onto those channels and returns the rigid +
      Zernike Jacobians as separate matrices (bitwise vs standalone) —
      NOT re-derived FD.  FD-from-scratch via `evaluate_` is the
      fallback only for a var type no channel covers yet (conic /
      asphere until a `ConicChannel` lands).  Measurement-space when a
      metrology backend is attached (§6.6).
- [x] Outer fmincon loop (single-field, nλ=1 default for
      all-reflective); merit `rms_wfe`.  Nested λ×field + callback
      merits are the follow-on.
- [x] First result: align-and-recover on e5hex1 (`example_align_from_rx`)
      — runs unchanged on a CodeV-converted Rx (the import path is
      Rx-agnostic).

### Sprint 2A-ii — `macos.design.Telescope` 2-mirror builder

The de-novo builder, landing onto an analysis core that already
works (2A-i).

> **CORE LANDED 2026-06-16** (MACOS_resources `sls-dev`: builder
> `e61549b`, first result `771de28`; reference set `2b7ab06`).
> `+macos/+design/Telescope.m` builds all four 2-mirror families →
> emits a MACOS `.in` → validates by SMACOS load.  `tDesignTelescope`
> 10/10 green (fast suite 143/0).  **The convention was de-risked in
> the runnable (pymacos) path against the shared fixtures BEFORE
> transcribing to MATLAB**: `KcElt=K`, `KrElt=-|R|`, `psiElt`→CoC
> (one rule, all surfaces — concave M1 & convex Cass secondary -z,
> concave Gregorian secondary +z), and the SMACOS-load-required
> `nOutCord`/`Tout` block (interactive CLI defaults it; SMACOS does
> not — `load_rx`→nElt=0 without it).  See the agent reference memory
> `reference_macos_rx_emission_convention`.

- [x] Spec struct + pure-function resolve/emit (§3); component model
      with Mirror / FocalPlane (Mask is a coronagraph element — defer
      to Sprint 3).  Fixed-topology families auto-populate M1/M2/FP.
- [x] Closed-form layout for Cass / RC / Gregorian / DK with the §5.2
      β-dependent forms; §5.3 raytrace validation.  Validated vs the
      shared fixtures (`optical_design/fixtures/`): R/K to ~1e-5;
      on-axis spherical-free (Cass & Gregorian machine-zero, the
      latter exercising the concave-secondary psi-flip); RC aplanat
      coma-squared signature (Cass/RC field-WFE ratio grows as
      field→0); DK largest coma.  `make_fixtures.py` is the ported
      closed-form source.
- [x] Full-precision (`%.16E`) `.in` emission; validate-by-load through
      SMACOS on every `build()`.  (Standalone `ValidatePrescription`
      wrapper + the emit-at-x-vs-x+δ FD-survivability diff: deferred —
      load-validation covers structural correctness today.)
- [x] `describe()` with provenance; `save_spec`/`load_spec`;
      deterministic-emit test (same spec → byte-identical `.in`, the
      §3 parity property).  (Committed-golden cross-language anchor:
      defer until the emitter stabilises in 2B/2C, when apertures /
      obscurations / field handling land.)
- [x] First result: a BUILT Cassegrain imported via `System` recovers
      an M2 despace error through `optimize`
      (`tDesignTelescope.test_built_telescope_feeds_alignment` +
      `examples/design/example_telescope_align.m`) — closes the
      builder→analysis loop.

**Deferred to follow-ons** (not blockers): the `ValidatePrescription`
standalone wrapper + FD-survivability diff; the committed-golden
byte-identical-Rx anchor (emitter still grows in 2B/2C); arbitrary
`optical_axis` (MVP is +z); `diagram()` + full `check_clipping()`
(Sprint 4); nested λ×field evaluation (inherited from 2A-i).

### Sprint 2B — N-mirror (TMA / 4-mirror / freeform)

> **DE-RISK + PROOF COMPLETE 2026-06-17** (pymacos scratch in
> `/tmp/{derisk_tma,korsch_search,proof_korsch}.py`; nothing built into
> the design layer yet — that is the remaining BUILD).  Every *mechanism*
> 2B needs is now proven; see [[project_design_layer]] +
> [[reference_macos_rx_emission_convention]] in the agent memory and the
> findings below.  Reason multi-mirror = WIDE FIELD: the value of each
> added mirror only shows under MULTI-FIELD optimization.

**Multi-field optimization — decided: MATLAB-driven, SELECTABLE engine**
`optimize('engine','native'|'fmincon')` (Dave 2026-06-17):
- **`native` (default fast path)** drives MACOS's *existing, mature*
  multi-field optimizer via **`calib_run`** (the CALIB wrappers, already
  in `macos_api_mod` + pymacos `m.calib`/mmacos `+macos/calib*`).
  `nls_optim_dvr` (LM) + SLSQP/NPSOL do **multi-field × multi-λ,
  FoV-weighted** least-squares over per-element DOFs (**DOF 7 = radius,
  DOF 8 = conic**), Zernike modes, and aspheric coeffs, targeting
  WFE/SPOT/WFE_ZMODE (up to 12 FOV × 6 λ).  `m.calib()` returns
  per-(FOV,λ) WFE before/after.
- **`fmincon` (flexible/research path)** = MATLAB field loop
  (`set_src_fov` → trace → WFE) + FoV-weighted merit + `dw_dx` analytic
  gradients; for multimodal freeform (MultiStart/patternsearch).
- **Structure constraint (Dave):** conics/radii/spacings are ONE shared
  physical system; only the chief ray (`OptChfRayDir`) varies per field.
  The native optimizer enforces this (one `aparams` vector; the
  `Do ifov` loop only changes the field).
- **Engine gap:** programmatic field-set setters
  `calib_add_fov`/`calib_set_wavelens` (Phase 1d, deferred).
  Today the FOV set comes from the `.in` AFOV keywords.

**Native-optim `.in` config keywords** (bake into the emitted Rx, then
`m.calib()`): header `OptTarget= WFE`, `OptWFElt= <FP elt>`,
`OptMaxItrs= N`, `OptFEX= No`; per-element `VarElt= 0 0 0 0 0 0 0 1`
(8-DOF mask, here DOF 8 = conic); field set = repeated `OptChfRayDir=` /
`OptChfRayPos=` + one `OptFOVWt=` of weights.  **GOTCHA:** the nominal
`ChfRayDir=` SHARES the `OptChfRayDir` parse block, so it IS field 1 —
emit `OptChfRayDir` for the OFF-axis fields only and size `OptFOVWt` to
`1 + n_off`, else a list-directed-read EOF crash.

**Emission convention (coaxial multi-mirror), validated:** all mirrors
`psiElt=(0,0,-1)` (ONE shared axis, NOT alternating-toward-CoC),
`KrElt=-|R|`, `KcElt=K`, vertices folded along z; matches the working
on-axis TMA `optiixonaxisz1_v4.in`.  Free-form layouts are the exception
(per-element psi — see `e5mono.in`).

**Credible design + proof (the result that matters):** the
`tma_fixture.json` is a 3rd-order-**Seidel regression toy** (R2=0.28,
f/1.5 primary) — huge higher-order, NOT buildable; don't demo on it.  A
**balanced Korsch** (R1=8,t1=3,R2=2,t2=4.5,R3=4,D=1, **f/8**), conics
seidel-solved, traces **0.17λ on-axis** — proving the *emission* is
correct.  `calib_run` over its 3 conics drove FoV-RMS WFE **1.4µm → 1.9nm
(0.003λ)** across 0–2.4′ — diffraction-limited.  Acceptance bar (Dave):
**WFE generally << 0.1λ AT THE EXIT PUPIL**, over a 2–20 arc-min field;
the convincing "each mirror earns its keep" demo is the wider field
where a 2-mirror hits its coma/astig wall and the TMA does not.

**Optimization target = exit-pupil WFE; is `add_pupil` needed in the
loop? (Dave's open question — likely NO):** MACOS's OPD at the focal
plane, evaluated over the ray grid, IS the exit-pupil-referenced
wavefront W(x,y) (each ray's OPD vs the reference sphere = the pupil
aberration).  The proof optimized exactly that (`OptWFElt = FP`, no
`add_pupil`) → 0.003λ.  So `add_pupil` is most likely for the
**deliverable** (an explicit, accessible exit pupil for downstream
instruments / coronagraph / Lyot), NOT a prerequisite of the
optimization flow, which already minimizes exit-pupil WFE via the FP
OPD.  CONFIRM during the build: compare the optimizer's FP-WFE against
the WFE read at an explicit `ExitPupil` Return surface (should match);
if MACOS's FP reference removes a term you need referenced to a SPECIFIC
pupil, that's the case where `add_pupil` would enter the loop.

**`add_pupil` (exit-pupil reference surfaces, wanted by Dave):** a 2-pass
op — (1) emit optics→FP, trace at a field, `m.fex()` to find the exit
pupil; (2) re-emit `Return@image → ExitPupil Return (= elt nElt-1) →
Detector FocalPlane (= elt nElt)`; `m.sxp()`/`fex` sets the EP radius at
trace time (positional convention: EP = nElt-1, FP = nElt).  Matches
`Rx_Cass_FarField.in` (Return1@image → ExitPupil@-4.06 → Detector).

**The BUILD — CORE SHIPPED 2026-06-18** (MACOS_resources `sls-dev`
`dbbb21f`; fast suite 149/0; `tDesignTelescope` +4 = 14):
- [x] N-mirror `Telescope`: `add_mirror`/`add_focal_plane` + spacing
      resolution (last = `'derive'` paraxial focus); all-psi=−z fold
      emission; **seidel-seed conics** — `macos.design.seidel_seed.m`
      (ported from `optical_design/seidel.py`+`make_tma_fixture.py`,
      validated: `seidel_seed([8 2 4],[3 4.5],1)`→`K=[-0.622 0.148
      -3.904]`, EFL=8).  (Deferred: equal-power-split default +
      `override('Mk','fnum',…)`; the builder takes explicit radii today.)
- [x] `add_pupil` (FEX-located EP; flat image-Return + spherical EP-Return
      before the FP, FP preserved; works for 2-mirror AND TMA).
- [x] `optimize('engine','native')`: emits the `Opt*` block + per-mirror
      `VarElt` conic DOF, runs `calib`, reads optimised conics back, emits
      the clean design.  **Generalised to 2-mirror AND N-mirror.**  Drives
      the Korsch seed 0.155λ → **0.003λ across 0–2.4′**.  (`fmincon` path
      still the deferred alternate.)
- [ ] `set_surface` → SrfType 2/4/12/14 dispatch; design vars →
      MonCoef/FFZernCoef emission (sparse FFZernModes).  **DEFERRED** —
      the freeform-surface vars; not needed for the conic-DOF TMA loop.
      (ZernTypeL dispatch ELSE-with-error already landed — §9.1 Q3.)
- [x] **Demo** `mmacos/examples/tma_widefield/example_tma_widefield.m`:
      RC (2-mirror) vs Korsch TMA optimised over 0–20′ — the RC **walls
      at 0.45–1.66λ** (2 conics can't null field astigmatism) while the
      TMA holds **< 0.1λ to ~15′** (11–26× gain).  Saves `.in`+`.mat`+a
      WFE-vs-field PNG, calls `add_pupil`, no `exit(0)`.  RC + align
      examples retrofitted to the same rules.
- [x] First wide-FoV **freeform** result — **LANDED 2026-06-24** as
      `mmacos/design/sz_tma/` (sphere+Zernike unobscured TMA, e5mono-derived,
      f/21). Builder gained `set_base_sphere` (Kc=0) + `add_mirror`'s `convex`
      flag (psi→downstream CoC) + convex-aware `paraxial_focus_` + a ZernModes
      single-line emit fix; `set_freeform`/`optimize_freeform` (CALIB OptZern)
      do the staged center→field correction. On-axis diffraction-limited
      (0.044λ from a 35700λ all-sphere start); 2-D ±2′ is field-limited (0.78λ
      area-wtd — narrow-field e5mono geometry + 15 modes). See
      `mmacos/design/README.md` + agent memory `project_sz_tma`. NOTE: Seidel
      mis-models a convex secondary's focus AND conics — focus fixed
      (paraxial_focus_, gated on any-convex; fixtures untouched), the conic fix
      is the approved follow-on (`project_seidel_convex_bug`). HELD/uncommitted.
      Next: pupil-referenced `optimize_freeform('engine','jacobian')` (linear
      dW/dZern solve, reuses the GMI `dw_dz_zernike` sensitivity) + widen the
      2-D field (M2 closer to input beam / intermediate focus closer to M2).
- [x] **Exit-pupil-surface evaluator** — **LANDED 2026-06-26.**  New
      engine command `XPS` (`tracesub.F`: FEX generalized to the whole
      grid — per-ray field-differential `FindCrossPt`; vertex = FEX to
      0 m) + `macos.pupil_quality()` (mmacos) fit the pupil SURFACE as
      low-order Zernikes in the FLAT entrance-pupil coordinate.  Pure
      geometry, no OPD (the field↔aperture dual of image-surface finding:
      hold the entrance position, vary direction — per Dave).  On
      `sz_tma` the 3-mirror exit pupil is **+1.67 mm defocus + 1.77 mm
      astig** — the curved+astigmatic pupil a **+1 near-focus field
      mirror** must null for a flat Lyot/DM/apodizer conjugate.  This is
      the measurement substrate for **simultaneous focal+pupil
      optimization** (pupil residuals stack into the same weighted vector
      as the focal OPD).  Engine: sls-dev `5bd1804`; tool: mmacos sls-dev
      `80645af`.  Next: in-engine Zernike fit + reference-element
      placement (the "place a pupil" capability for the FSM/EP/apodizer/
      Lyot stations), then feed pupil residuals to the jacobian/lsqnonlin
      engine.

> **Example-building rules (Dave 2026-06-17, apply to ALL real examples
> incl. the 2A-ii 2-mirror ones — retrofit them):** live in a
> `mmacos/examples/<dir>/`; `save()` the `.in` + `save_spec()` the `.mat`
> at the end; call `add_pupil` so the exit pupil is accessible; and do
> **NOT** end in `exit(0)` (that batch-mode rule is only for the
> `test_*.m` smoke scripts + matlab.unittest, never a user demo).

### Sprint 2B+ — the 3+1 coronagraph front end (SHIPPED 2026-07-04)

> **PRODUCT FRAME (Dave 2026-07-04, standing):** the real product of
> the design layer is **utilities and examples users adapt** for
> related design studies, in a telescope PROGRESSION — 2-mirror →
> 3-mirror → **3+1** → N-mirror — then an INSTRUMENT-building
> sequence.  Layout: `mmacos/design/src/` (script-level utilities) +
> `mmacos/design/examples/` (all drivers moved there, `bcccea6`).
> Future: refit the 2-mirror examples onto the utilities structure.

**Standing design constraints (Dave 2026-07-04):** (1) *packaging* —
the telescope must fit a **cylindrical launch shroud**; keep M2 close
to the incoming beam as PM–SM grows (`packaging_report` measures the
shroud envelope + train length); (2) *coronagraph polarization* — the
per-mirror **AOI SPREAD across the rays of the beam** < 15° preferred
(`aoi_report`; the spread at M1 is the primary's own convergence
≈ D/R1, so the knob is a slower primary = longer PM–SM — directly
against packaging); (3) relays serve small per-instrument patches;
the shared wide field lives at the TMA focus (a near-focus relay
cannot carry ±2.5′ — the image walks ±96 mm across M4); (4) the HWO
wide field (10×20′, multi-instrument) needs a 4th POWERED imaging
mirror — field-quadratic astigmatism is what walls the 3-mirror
(measured: conic wall 3′-square/4′-circular DL @1µm, 5′-circular
0.095λ; freeform M2+M3 buys only 1.4× because the astig orientation
varies over the field; per-field astig removal → 0.036λ = the proof).

Shipped (`MACOS_resources` sls-dev `f817f61` + `04295fa`, tests
48/48): `examples/tma_centered` (obscured-vs-unobscured A/B — the
"why does j18 do better" answer: symmetry intact 0.065λ vs section
0.098λ on the 5′ circular field @1µm, + the 2.3µm yardstick);
`examples/tma_3plus1` per Dave's 3-file structure — the j18 DEMO
(0.84D compact, 30″ patch 0.034λ, pupil relayed ~10× flatter, AOI
spread 21/24° documented out-of-preference), the **AOI constraint
finder** (steps PM–SM, verifies the FULL 4-mirror chain: f/2.0 =
1.7×j18 separation MEETS 15°, but decenter grows 0.71→1.50D →
**shroud 1.6→3.2×D: AOI-safe eccentric sections are
shroud-expensive**), and the polsafe optimize (0.061λ patch, 5/5
clear, AOI 14.2°).  Utilities: `tma_conic_recipe`, `wfe_field_diag`,
`aoi_report`, `packaging_report` (design/src); `field_ring`,
`Telescope.trace_at_field`, `add_mirror 'conic'` seeds,
`optimize 'elts'` subset, psi parity ≥4th mirror (package).
**Next in this thread:** the tilted-fold 3+1 (Bauer folds hug the
incoming beam — the shroud-cheap unobscuring; fold path exists);
pupil fine-tune (K4 scan vs `pupil_quality`); then Sprint 2D
segments this parent's M1.

### Sprint 2C — refractive + dispersive components

- [ ] Lens (2-surface expansion, glass names, group perturbation) and
      Grating (λ-dependent exit direction) per §6.
- [ ] Glass-catalog expansion + .agf converter (§6.4).
- [ ] λ_ref layout machinery + dispersion-aware detector sizing.
- [ ] Worked example: §7.4 spectrograph; chromatic merit exercises
      nλ>1 for real.

### Sprint 2D — segmentation (SegMirMaker orchestration)

> Sequenced here for readability; per the v3 review it may land
> after 2A-ii / parallel to 2B at Dave's discretion.

**Orchestrate, do not reimplement.**  The psi-flip detection, FF
replication with ZerntoMon dispatch, and per-segment Mon-frame
conventions are debugged logic that must not be duplicated in MATLAB.

- [ ] **Q7 — SegMirMaker batch mode** (§9.1): the nine interactive
      dialog answers driven from a control file or args; interactive
      mode retained for standalone users.
- [ ] Design-layer `segment('M1','rings',2,'gap',…,'dofs',6)`:
      emit unsegmented Rx → run SegMirMaker batch naming the parent
      element → splice `.presc` segment blocks in place of the parent
      with downstream renumbering (the emitter owns numbering) →
      ingest `Hx.m`.  Segments carry provenance
      `derived(segmirmaker)`.
- [ ] Post-splice validation: `ValidatePrescription` + the
      SEGRAYTRACE center-ray check (segment-center trace endpoint vs
      `RptElt`) run automatically.

### Sprint 3 — `macos.design.Coronagraph` + DarkZone target

- [ ] `Coronagraph` class per §7.3; conjugate-recreating OAP solve.
- [ ] **New CALIB target type `DarkZone`** — annular image-plane
      integral (`inner_lamD`, `outer_lamD` or absolute radii),
      evaluated through the diffraction chain (per E3).
- [ ] `macos_api_mod` wrapper `calib_set_target_darkzone(...)` +
      mmacos surface + tests.
- [ ] **Scope guard:** inner-loop DM control is MODAL (tens of
      Zernike/Fourier modes) — SLSQP over thousands of raw actuators
      via FD is not viable.  High-actuator-count EFC (Jacobian-based)
      belongs to FALCO integration (PLAN.md §2.3); DarkZone is scoped
      accordingly.
- [ ] Warm-started inner loop between outer iterates (§1.3.3).

### Sprint 4 — diagnostics + LAYOUT REALIZABILITY (pulled forward 2026-06-19)

The 2B coaxial TMA emit SELF-OBSCURES (all psiElt=(0,0,-1), vertices on
the z-axis: M1@0, M2@-3, M3@+1.5, FP@-0.5 -> M1 and the FP sit ON the
M2->M3 beam).  It "traces fine" only because MACOS walks elements in LIST
order and never checks body-occults-beam.  So the layout tooling came
early.  Off-axis reference: `~/dev/MACOS_sandbox/old_Rx/dmt6mono.in`
(on-axis conic TMA -> decenter/tilt + flat FOLD mirrors -> beam clears;
nECoord=6/TElt frames).

- [x] **Real-ray viewer via a DATA-only DRAW path** (macos `f3e98e5`;
      MACOS_resources `caf3b86`; both PUSHED).  Engine: `src_mod`
      `DrawDataOnly` + capture arrays; DRAW handler copies its bundle +
      skips render under the flag; `macos_api_mod` `draw_rays_cmd`/
      `draw_rays_get` (real ray bundle as data, no Giza -- DRAW already
      does plane YZ/XZ/XY + first/last-elt slicing).  mmacos:
      `+macos/draw_rays.m`, `Telescope.diagram()` (cheap marginal pass)
      and `Telescope.view_layout(plane,opts)` -- real rays + conic-sag
      surfaces drawn to the actual beam FOOTPRINT, opts plane/istart/
      iend/hide/nrays.  Reveals the coaxial obscuration directly.
- [x] **`check_clipping()`** (DONE, 2026-06-19) -- reconstructs the real
      ray bundle in full 3-D from TWO orthogonal DRAW projections (YZ->Y,Z
      ; XZ->X,Z ; shared Z = integrity check -- orthographic axis-pick
      confirmed in `pgplotsub.F`, so 2 plane-calls suffice, not 3), then
      tests each physical body (disk Vpt/psi/ap_r) for piercing a beam
      segment between two OTHER elements.  Judged in 3-D (a 2-D projection
      paints FALSE conflicts).  Per-element struct (name/kind/ap_r/foot_r/
      margin/obstructs/ok); prints a table; wired into `build('check',true)`
      as a warn-only gate.  Correctly flags the coaxial TMA (M1+FP in the
      converging M2->M3 beam) and reports Cassegrain central obscuration
      (secondary in the incoming beam).  2 tests in `tDesignTelescope`
      (16/16 green).  Pure MATLAB -- rides the existing draw_rays getter,
      no engine rebuild.
- [x] **Off-axis design generation a+b** (DONE 2026-06-20, sls-dev,
      pushed).  Take the on-axis TMA off-axis WITHOUT folds, vertices
      pinned + axes aligned (e5mono/dmt6mono "design on-axis then move
      off-axis" recipe):
      - `set_field_bias(arcmin)` (tilt ChfRayDir; optimize() centers its
        eval fields on the bias) `1ff9790`; `set_aperture_decenter(m)`
        (off-axis used beam) `634a0bd`.  Off-axis mirrors emit `ApType=
        None` (no clip); `aperture_full_field()` = per-elt clear aperture
        (radius,xc,yc) over the FoV, emit-ready as Circular `ApVec=(r,xc,
        yc)` `ba19f93`.
      - `optimize('dofs',mask)` over `[TIP TILT CLOCK DX DY PIST ROC
        CONIC]` (dopt_mod.F:43): CALIB tilt/decenter+conic, bakes into
        psiElt/VptElt; clean re-emit reproduces (conics rot-sym so TElt
        roll is trace-neutral).  +6' TMA 1.6e-6 -> 3e-9 diff-limited.
        Sensible `TElt` emitted = local surface frame (Z=outward normal,
        XY tangent; trace-neutral PERTURB/sensitivity interface) `2867370`.
      - `check_clipping` clearance: signed clearance per elt (`d1664d4`)
        using the off-axis footprint PATCH as the body (`da3af58`); a
        ~1/2-aperture decenter clears M2+FP.  ENGINE BUG noted (not fixed):
        `macos.save_rx` malformed-format crash `macos_IO.f90:200`.
- [x] **Off-axis SECTION (RptElt != VptElt) + unobscured RC + AREA-mode
      optimize + report views** (SHIPPED 2026-06-21, sls-dev `bc21082`).
      Engine-true off-axis primitive: VptElt = parent figure VERTEX,
      psiElt = parent AXIS, RptElt = section POLE (used sub-aperture center;
      carries the TElt / perturbation center).  `ConSrf` takes the conic sag
      from VptElt ONLY, so an off-axis section is just RptElt != VptElt on
      the SAME parent figure (the j18sc segmented-PM pattern).
      - `set_offaxis(clear,...)`: REQUIRED named target ('M3'|'all'|cellstr|
        'none'); `clearance_solve_` bisects the decenter on `check_clipping`;
        `resolve_section_poles_` sets poles + analytic section normals
        (trace-neutral).  RC binding body = M1 (return beam), clears d=0.89D.
      - `optimize('fields',(:,2))`: arbitrary 2-D (thx,thy) field set,
        auto-drops the on-axis row (= implicit field 1); the +y
        `fields_arcmin` path is byte-unchanged (direction cosines).
        `macos.design.field_grid` / `field_cross` build the AREA + CROSS
        field-set modes (CALIB's 12-FoV cap -> area-optimize practical to a
        3x3 grid = 9 FoV; finer grids are evaluation-only).
      - `view_field_map` (WFE contour/surf over the 2-D field) +
        `view_orthoviews` (multi-plane layout report; view_layout core
        refactored into a shared `draw_plane_` helper).  LAYOUT DEPTH FIX:
        the conic sag uses the FULL transverse radius sqrt(h^2 + woff^2),
        the out-of-plane center per element from an ORTHOGONAL DRAW fan
        (`beam_offplane_`), so an off-axis section sits at the right depth in
        the cross-plane view (the M1-in-XZ offset).  FP/EP/Return curves
        drawn to the real beam FOOTPRINT; vertical label collision-avoidance.
      - `realize_apertures` sizes clear apertures to the full-field 2-D
        footprint (both fans); `add_pupil`'s flat image-Return renamed
        `PupImg` -> `FP_return` (it sits at a focal plane, not a pupil image).
      - Worked example `examples/design/rc_offaxis/example_rc_offaxis.m`:
        ONE FoV input -> on-axis 0.0019lam -> set_offaxis('all') d=0.89D
        0.75lam -> axial optimize(roc+conic) 0.0033lam -> 3x3 AREA optimize
        (tip/tilt/dy+roc+conic) worst-field 0.49 -> 0.046lam diff-limited &
        CLEAR -> 7x7 WFE map + x/y cross-sections + orthoviews.  Rigid DOFs
        (not conics) correct the linear off-axis astigmatism.
        `tDesignTelescope` 32/32, fast suite 167/0.
- [ ] **Off-axis follow-ons** (NEXT thrust): OAP relay (FP -> collimating
      OAP -> DM/filter space -> refocusing OAP -> FPM ...) -> JWST-style
      field-bias TMA (reuse `set_offaxis` named-clearance + the field-bias
      knob) -> `set_surface` freeform (SrfType 14 + FFZernCoef) -> fold flats
      (TElt frames).  Everything downstream already runs on whatever `emit_`
      produces.  (The old `examples/design/tma_offaxis/` is superseded by
      rc_offaxis -- left untracked, pending delete.)
- [x] Examples reorg (DONE, `d9978d8`): `mmacos/examples/` ->
      sensitivities / design / coronagraph.
- [ ] `plot_history`; all worked examples in the manual; `diagram()` 3-D;
      tag `design-layer-v1`.

### Sprint 5 — simultaneous focal + pupil optimization (added 2026-07-02)

**Goal (Dave):** n-mirror telescopes with diffraction-limited multi-field
WFE at an ACCESSIBLE focus AND a flat, well-imaged, field-stable,
ACCESSIBLE exit pupil (a flat DM/Lyot/apodizer/FSM conjugate).

> **Deferred here from the engine FEX rework (Dave, 2026-07-03):**
> *should the pupil finders (FEX/XPS/pupil_quality) include OBSCURED
> rays?*  Raised by `eac2_7seg` (coronagraph whose mask obscures all
> rays; its FEX legs disagree −289.7 vs +52597.6).  Engine state when
> parked: FEX's chief-pair vertex/radius and the new beam-footprint
> guard already use obscured rays (geometric trace continues past
> obscuration, `LRayOK` gating); FEX's CENTROID mode honors `iObsOpt`
> (default = unobscured only) and, when every ray is obscured, builds
> `psip` from an uninitialized `CentroidSpot` — latent garbage-normal
> path, same family as the §0 ORS chief-ray bug.  Needs more work —
> take up alongside J2–J5 pupil metrics, where coronagraph pupils
> (annular/occulted beams) are first-class.

**Structural insight (proven by sz_tma):** Zernike departures don't move
chief-ray geometry, so the pupil objectives are functions of the
0th-order layout alone (sphere radii, spacings, tilts, stop) while
focal WFE is recoverable by the Zernike layer for any packaging-feasible
layout.  → nested optimization: OUTER loop over layout knobs scores
pupil + accessibility, INNER CALIB Zernike polish supplies best-achievable
WFE per layout.  **The lever that decouples the two:** a +1 field mirror
at/near the real intermediate focus (3+1 config) — power at an image
plane has ~zero first-order effect on focal WFE but full authority over
pupil curvature; its astigmatic figure controls pupil astig.  sz_tma's
baseline error to null: **+1.67 mm pupil defocus + 1.77 mm astig**
(`macos.pupil_quality`, landed 2026-06-26 — see Sprint 2B tail).

**Objective set** (each → a measurable; weights carry unit conversion,
focal residuals in waves, pupil residuals in mm):

| # | Objective | Measurement | Status |
|---|---|---|---|
| J1 | Area-weighted multi-field RMS WFE (+ worst-field check) | CALIB / `optimize_freeform` staged S0→S1→S2 | exists |
| J2 | Pupil surface flatness: pupil defocus Z4 + astig (or `sag_rms`, piston/tilt removed) | `pupil_quality` | exists |
| J3 | Pupil image sharpness: crossing-cloud thinness | `pq.fit_resid` | exists |
| J4 | Pupil mapping linearity: RMS nonlinear residual of (U,V)→transverse-crossing linear fit | data in `pq.uv`+cloud; metric NOT coded | slice 1 |
| J5 | Pupil stability over field: vertex walk + sag-coef spread across the field set | needs multi-field loop | slice 2 |
| C1 | Accessibility constraints: pupil vertex real/external w/ clearance margin; FP out of beam; unobscured; 100% ray pass | `check_clipping` + ray-loss guard; pupil-station body NOT wired | slice 3 |

**Slices (all MATLAB over existing engine commands — no Fortran):**

- [ ] **Slice 1 — distortion metric.**  Extend `pupil_quality`:
      fit the best linear map (U,V)→transverse crossing position;
      return `pq.distortion_rms` (fraction of pupil radius), `pq.mag`,
      anamorphism.  Pin sz_tma values in a `tPupilQuality` test.
- [ ] **Slice 2 — `pupil_quality_multi`.**  Loop `set_src_fov` →
      `pupil_quality` → restore nominal (mirror the `dw_*_multi`
      supervisor pattern).  Returns per-field struct array + variation
      stats (lateral/axial vertex walk, per-coef spread).  Test on
      sz_tma over the S2 field set.
- [ ] **Slice 3 — pupil-station accessibility.**  Add the pupil
      vertex + a placeholder flat as a body in `check_clipping`
      (reuse the off-axis footprint-patch machinery); signed clearance
      → outer-loop constraint.
- [ ] **Slice 4 — objective wiring.**  Stack pupil residuals into the
      outer optimizer as ONE weighted least-squares vector
      `[w_f·OPD(fields); w_p·(Z4,astig)_pupil; w_d·distortion;
      w_s·field-stability]` (the "pupil residuals stack into the same
      weighted vector as the focal OPD" intent recorded in the 2B
      tail).  Default architecture = nested (outer fmincon/SLSQP over
      layout knobs, inner CALIB Zernike); single-vector jacobian/
      lsqnonlin engine as the follow-on.  Weights from PHYSICAL
      tolerances, not tuning: WFE 0.07 waves; pupil sag from the
      Talbot/Fresnel phase→amplitude budget (coronagraph) or
      registration blur ≪ actuator pitch (metrology); distortion as
      pitch fraction.
- [ ] **Slice 5 — 3+1 builder support.**  `add_mirror` placement at/
      near the intermediate focus (`'recreates','pupil'`-style
      conjugate solve, §7.3 keyword) so the +1 field mirror is a
      first-class design variable (radius → pupil defocus; Zernike
      astig figure → pupil astig).
- [ ] **Slice 6 — worked example `design/tma_3plus1/`.**  Extend
      sz_tma with the +1 field mirror; staged run; before/after
      `pupil_quality` (1.67 mm defocus / 1.77 mm astig → ~0) with WFE
      held diffraction-limited; per example rules (save .in + .mat,
      `add_pupil`, figures, no `exit(0)`).  (An untracked
      `design/tma_3plus1/tma_3plus1.m` stub exists — absorb it.)
- [ ] Tests grow alongside each slice (standing rule).

**Acceptance:** a 3+1 derived from sz_tma with worst-field WFE
< 0.07 waves over ±2′, pupil sag (piston/tilt removed) reduced ≥10×
vs the sz_tma baseline, distortion + field-stability reported, pupil
AND focus accessible per `check_clipping`, committed with tests.

**Open questions (record answers in §10):** (a) pupil-sag tolerance
normalization — Talbot budget vs registration budget, per use case;
(b) do pupil terms ever enter the CALIB inner loop (engine work) or
stay outer-only (default: outer-only until measured to matter);
(c) in-engine pupil-Zernike fit + "place a pupil" reference-element
placement (already queued in the 2B tail "Next") — pull only when
slice 4 wants speed.

### Sprint 6+ (deferred)

- Polychromatic CALIB (broadband inner-loop merit in one solve).
- `spectral_run(λ_list, …)` amortization convenience.
- Clearance-derived OAP off-axis distances; Talbot DSL keyword.
- Surrogate outer optimizers (Bayesian / SMAC) if eval cost demands.
- `export_proper()` → PROPER/FALCO prescription for end-to-end
  contrast (the system-model schema makes this an emitter, §3).
- Python port of the builder — struct translation + transliteration,
  parity = byte-identical Rx on golden specs.
- Segmented-primary integration is **no longer deferred** — see the
  Sprint 2D slice (SegMirMaker orchestration).
- **Truss / metrology optimization** (§6.6 tier 3): fiducial +
  launcher placement and beam topology optimized for DOF
  observability.  First cut pure-MATLAB kinematics; MACOS validates
  selected configs through the real gauge functions + traces met-beam
  clearance.  Gated on §9.1 Q8.

---

## 9. MACOS-side support work — when each piece is needed

| Sprint | macos / `macos_api_mod` change | Why |
|---|---|---|
| 0 | (none — experiments only) | design inputs |
| 1 | `set_src_wvl` / `set_src_fov` already present in mmacos (verification only); ray-loss per-category breakdown (Q2 — small).  E1–E4 experiments (no other Fortran — measurements scope Sprint 3) | FoV loop; §1.3 guard; macos↔MATLAB split |
| 2A-i | (none in Fortran) | import / analysis core entirely MATLAB over existing getters |
| 2A-ii | (none in Fortran) | builder + outer loop entirely MATLAB |
| 2B | (none — ZernTypeL dispatch ELSE-with-error already landed, §9.1 Q3 / PLAN.md §0) | silent no-op guard closed |
| 2C | catalog formula check (possibly parser extension); .agf converter (offline tool) | refractive components |
| 2D | Q7 — SegMirMaker batch mode (control-file / arg-driven) | segmentation orchestration |
| 3 | `DarkZone` CALIB target + wrapper | coronagraph inner loop |
| 4 | (none) | docs + diagnostics |
| (validation) | Q8 — met-function validation harness (closed-form geometric tests) | gate before metrology budget work (§6.6) |

### 9.1 Separable engine queue vs. joint development

Split criterion: **where the item's spec comes from.**  Specs pinned
by existing engine semantics → separable, work now in macos/smacos
with their own acceptance tests (pymacos pytest is the fastest
engine-level harness; mmacos picks the routines up via a codegen
re-run).  Specs downstream of Sprint 1 measurements or design-layer
API shape → joint, after E1–E4.

**Separable queue (spec-ready today):**

| # | Item | Acceptance |
|---|---|---|
| Q1 | ~~`set_src_dir` + `set_wavelength` in `macos_api_mod`~~ — **already landed:** `set_src_fov` (absolute pointing) and `set_src_wvl` are in mmacos via the existing `macos_api_mod` setters.  Sprint 1 verifies only. | pymacos pytest vs journal-driven references |
| Q2 ✅ **landed 2026-06-11** | Per-category ray-status getter in `macos_api_mod` (`ray_status_get`) exposing `RayStatus(:)` + `RayFailElt(:)` verbatim from `elt_mod`.  mmacos surface `get_ray_status(N)` returns the integer-coded category per ray (`RayStat_OK/Obscured/Miss/Bracket/MaxIter/Undef`) + per-category counters — complements the binary `get_ray_info`.  Codegen Path A.  (Also bundled the latent `libslsqplib.a` mmacos-link fix.) | **Met:** `tMacosPkg` pins Rx_Cass_FarField counts (12850 rays = 11484 OK + 1366 Obscured; matches engine "Obscured: 1366") + cross-check vs `get_ray_info` |
| Q3 ✅ **landed (verified 2026-06-11)** | ZernTypeL dispatch ELSE-with-error (propsub.F / srtrace.F / tracesub.F).  **Owned by PLAN.md §0** (commit b2c2eb8) — cross-ref, not a separate track.  Both propsub.F and srtrace.F now handle Noll (ZerntoMon6) + carry an ELSE-with-error for unhandled types. | **Met:** `ZernType= Noll` produces non-zero OPD response; unhandled types error loudly |
| Q4 | Glass catalog: formula-coverage audit → parser extension if needed → .agf converter + generated usual set | n(λ) vs published to 1e-6; one doublet vs CodeV |
| Q5 ✅ **PASS 2026-06-12** | Endurance test (load/trace, one session) — `tEndurance.m`.  No findings to fix: rmsWFE bit-identical across 500 iters, memory plateaus (no linear leak). | **Met:** bit-identical rmsWFE each iter; flat steady-state memory |
| Q6 | `Element= Apodizer` (independently motivated, PLAN.md §2.1 Thrust B) | participates in subsequent trace; PROPER cross-check |
| Q7 | SegMirMaker batch mode — drive the nine interactive dialog answers from a control file / args; interactive mode retained for standalone users | batch run on `test_in/` parents reproduces interactive output byte-identically |
| Q8 | Met-function validation harness (**VALIDATION/CHARACTERIZATION, not implementation**) — John Lou's gauge functions checked against closed-form geometric truths (LOS-projection equality, null test, launcher/target reciprocity, FD-vs-analytic).  **The PLAN.md §4.5 PERTURB coverage gaps are expected Q8 failures** — `SrfMetPos` is updated by only 2 of 5 perturbation paths (`CPRead` / `CPERTURB_2` / `LnkEltCPERTURB` don't), so any Q8 test that perturbs through those paths reads stale metrology.  Q8 should surface them as known failures pointing at PLAN.md §4.5, not silently pass. | pass/fail map per function + missing-capability inventory (spec input for joint completion) |
| Q9 → **PLAN.md** | Nominal snapshot/restore (`+macos.snapshot` / `+macos.restore`).  **Owned by PLAN.md §11.7 (Phase 7) / tracked MISSING in §11.6**; spec = GMI `ObtainNominalSettings` field list (PLAN.md §11.5).  Cross-ref only — no duplicate checkbox here.  **`evaluate_`'s `from_rx` setter path depends on it:** the §1.1 imported-geometry row ("no reload per outer step") needs a lightweight restore-to-nominal between outer iterates; without it, restore falls back to full `load_rx` (as `dw_dx.m`'s `reload_rx=true` does today), defeating the import-path caching win. | (acceptance owned by PLAN.md) |

**Held for joint development (spec depends on measurements / API):**

- DarkZone CALIB target + wrapper — E2 sets the bar it must beat,
  E3 its cost budget; modal-vs-actuator scope shapes the interface.
- `spectral_run` — exists only if E4 shows a real win.
- Phase 1d CALIB setters (`calib_add_fov`, `calib_set_wavelens`) —
  spec-ready from existing AFOV semantics; pull the moment E3/E4
  want them, don't pre-build.

> **CC:** Queue items Q1–Q6 are self-contained work packages: each
> gets its own branch-local test before any design-layer code
> consumes it.  Keep the queue this size — every engine item must
> justify itself independently of the design layer; anything that
> only makes sense IF the design layer wants it is joint-side by
> definition.

---

## 10. Decisions

### Made

1. **MATLAB-first**; Python port deferred, kept cheap via §3.
2. **Rx text emission** (not API-driven element construction) —
   debuggable, version-controllable, language-neutral.
3. **DarkZone in both tiers:** MACOS-side target for the inner loop;
   MATLAB merit for the outer (flexibility > speed there).
4. **Outer optimizer `fmincon`**; `patternsearch` as the
   gradient-noise fallback; surrogate methods deferred.
5. **Coronagraph-driven from Sprint 1** — the evaluation surface and
   the macos↔MATLAB split are proven on the existing coronagraph Rx
   corpus (no builder).  Telescope remains the first *builder*
   example (Sprint 2A-ii: layout math, no DMs); spectrograph third.
   Import (`from_rx`) analysis lands even earlier, in Sprint 2A-i.
6. **macos↔MATLAB division of labor is decided by measurement**
   (Sprint 1 E1–E4), not upfront.  Default: logic starts MATLAB-side;
   it moves into Fortran only for inner-loop speed (DarkZone) or to
   participate in propagation (Apodizer).
7. **nλ=1 default** keyed to *all-reflective*, not to merit type.
8. **Component as topology unit** from Sprint 2A-ii.
9. **Modal inner-loop DM control**; actuator-level EFC → FALCO.
10. **Two front-ends, one analysis core** (§1.0); **import
    (`from_rx`) is the expected dominant entry point.**  Converters
    are importers feeding `from_rx`; `describe()` / `check_clipping()`
    / `ValidatePrescription` are the diagnostic surface crude
    conversions never had.  Zemax completion stays deprioritized —
    the future path (populate the component model directly + share
    the emitter) is smaller than finishing a text-to-text script.
11. **Units policy: bare SI canonical** at the user surface (matches
    the `+macos/` convention).  `_<unit>` suffixed argument names and
    `'unit',<name>` options accepted as documented sugar;
    `'unit','lambdaD'` for mask radii.  Examples keep mm / mrad
    sugar for readability; wavelengths bare SI.  *(Confirmed
    2026-06-11.)*
12. **Metrology measurement space = a Jacobian contract** (§6.6);
    backends pluggable; SegMirMaker edge sensors (`Hx`) first; engine
    met functions gated on Q8 validation.
13. **Segmentation by SegMirMaker orchestration** (batch + splice,
    Sprint 2D), never by reimplementation in MATLAB.

### Open

1. FoV weighting model (equal-weight default; Strehl-weighted and
   worst-case as built-ins).
2. Bandwidth spec sugar (`'500-900 nm @ 5'`) over the λ-list
   primitive; named line sets (`'FdC'`).
3. Family coverage beyond the §5.2 seven (aplanatic Gregorian,
   Schiefspiegler, off-axis Cass) — add when a need arrives.
4. Better N-mirror first-cut than equal-power split (Korsch-class
   TMA equations) — evaluate after Sprint 2B against a real wide-FoV
   prescription.
5. Power-split sign convention default (alternate, concave M1) —
   revisit if a real layout fights it.
6. MACOS catalog dispersion-formula coverage (Sprint 2C check).

---

## 11. Branch routing

- All work on **`sls-dev`** (macos) + companion `sls-dev`
  (MACOS_resources), per PLAN.md §9.
- Design package: `MACOS_resources/mmacos/src/+macos/+design/`.
- DarkZone target (Sprint 3) → `sls-dev`; fast-forward to `opt-dev`
  at release-worthy sprint ends.
- Tags: `design-sprint-0` … `design-sprint-4`, `design-layer-v1`.
