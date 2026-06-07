# Design Layer — Plan v2

A high-level design surface on top of MACOS. The user expresses
*design intent* (topology, free parameters, merit function); a builder
generates a concrete MACOS prescription directly (no CodeV step);
MACOS handles ray trace, diffraction, and inner-loop alignment; an
outer optimizer (MATLAB `fmincon`) explores the design space.

Revised 2026-06 from the original draft + design review.  Companion
to `PLAN.md`.  The evaluation surface is proven coronagraph-first on
the existing Rx corpus (Sprint 1); the telescope is the first
*builder* example (Sprint 2); a simple spectrograph (lens + grating)
the third.
**MATLAB-first** — the JPL user base is MATLAB-heavy.  The Python
port stays cheap by construction (see §3 state-as-data rule).

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

---

## 1. Architecture: three tiers, two loops

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
| outer step (geometry) | spacings, Kr/Kc, coefficients, mask radii | builder re-derives → re-emit `.in` → `load_rx` (once per step) |
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

> **CC:** Design for worker-safety from Sprint 2A even if
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
   `declare_design_var` bounds are used to normalize internally to
   [0,1]; `builder.update` unnormalizes.  fmincon never sees raw
   units.
3. **Inner-loop determinism.**  Inner CALIB runs with fixed
   iteration budget + tolerance much tighter than the outer FD
   step's merit effect; warm-start the inner solution from the
   previous outer iterate (cost AND smoothness).  If gradient noise
   persists anyway, `patternsearch` is the fallback, not tighter
   SQP knobs.
4. **Ray-loss guard (Sprint 2A, not Sprint 4).**  RMS-WFE-over-
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

> **CC:** Items 1, 2, 4, 5 are Sprint 2A acceptance criteria.  Write
> the FD-survivability test explicitly: emit Rx at x and x+1e-8·x,
> diff the files, assert the parameter actually changed.

---

## 2. User interaction flow

The intended session shape, end to end.  Stages 1–3 are iterative
and cheap; stage 4 is the expensive loop; stages 5–6 close out.

```matlab
%% Stage 1 — declare intent (seconds; no engine calls)
t = macos.design.Telescope('family','RC', ...
        'aperture_diameter_mm', 6000, 'primary_fnum', 2.0, ...
        'system_fnum', 20.0, 'BFD_mm', 1000, 'model_size', 256);
t.add_mirror('M1');
t.add_mirror('M2');
t.add_focal_plane('FP');
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
t.declare_design_var('M2_despace_mm', 0, 'bounds',[-2 2]);
t.declare_design_var('M2_conic',     [], 'bounds',[-1.2 -0.8]);
t.set_outer_merit('RMSWFE_FoV_band_averaged');
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
| design vars + bounds (+ natural units) | required for `optimize()` |
| outer merit | default `RMSWFE_FoV_band_averaged`; callback override |
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

### 5.3 Family validation tests (Sprint 2A acceptance)

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
  Sprint 4, but pass/fail margin check runs from Sprint 2A); BFD > 0;
  Cass-class FP reachable through the M1 hole; Gregorian intermediate
  focus exists (`d_M2 > f_M1`); N-mirror power sum consistent with
  f_sys; emitted Rx passes `ValidatePrescription`.
- Override rule unchanged: any table entry is overridable by
  `t.override('<elt>','<param>',value)`; an override disables that
  one derivation, everything else proceeds.  Design vars declared on
  a derived quantity take precedence over the derivation (an
  optimizer-driven override).

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

> **CC:** Implement the contract in Sprint 2A with Mirror, Mask,
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

> **CC (Fortran, Sprint 2B):** add the missing ELSE-with-error to the
> propsub.F / srtrace.F ZernTypeL dispatch chains — the design layer
> is about to become the first heavy user of SrfType 14 + FFZernType.

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

---

## 7. Builder API examples

### 7.1 Classical 2-mirror Cassegrain-class

```matlab
t = macos.design.Telescope( ...
        'aperture_diameter_mm', 6000, 'primary_fnum', 2.0, ...
        'system_fnum', 20.0, 'BFD_mm', 1000, ...
        'optical_axis', [0 0 1], 'family', 'Cassegrain');
t.add_mirror('M1');
t.add_mirror('M2');
t.add_focal_plane('FP');
t.declare_design_var('M2_despace_mm',   0.0,  'bounds',[-2 2]);
t.declare_design_var('M2_tilt_xy_mrad', [0 0],'bounds',[-1 1]);
t.declare_design_var('M2_conic',        [],   'bounds',[-1.2 -0.8]);
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
t.declare_design_var('M2_asphere_coeffs', zeros(1,3), 'bounds',[-1 1]*1e-12);
t.declare_design_var('M3_zern_coeffs',    zeros(1,12),'bounds',[-1 1]*1e-9);
t.set_field_points(macos.design.hexgrid(deg2rad(0.05), 2));
t.set_outer_merit('RMSWFE_FoV_band_averaged');
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
c.declare_design_var('FPMRadius',     4.0,  'bounds',[2 8],     'unit','lambdaD');
c.declare_design_var('LyotUndersize', 0.85, 'bounds',[0.7 0.95]);
c.set_bandwidth([550 600 650]*1e-9);
c.set_inner_target('DarkZone', 'inner_lamD',7, 'outer_lamD',10);
c.set_outer_merit(@(I_stack, lam, f) my_avg_contrast(I_stack, 7, 10));
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
s.set_outer_merit('RMS_spot_lambda_field');
```

Exercises in one example: 2-surface expansion, glass-name emission,
λ-dependent exit direction, λ_ref layout, dispersion-aware detector
sizing, chromatic (nλ>1) merit.

---

## 8. Sprint sequence

### Sprint 0 — pre-flight experiments (hours each; do before 2A)

- [ ] **Endurance test:** loop `load_rx` → `trace` 500× on Rx_Cass in
      one session; assert bit-identical rmsWFE every iteration + flat
      memory.  This is the design loop's exact load profile; no
      current test exercises it.
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
- [ ] Expose a **ray-loss summary** (`nBadRays` + per-category counts
      from RayStatus) — backend for the §1.3 guard.
- [ ] Verify CALIB wrappers handle DM-element variables on `Rx_Coro.in`.
- [ ] Tests: load `SegDemo3.in`, set field point, trace, read WFE;
      expected values from pymacos.

**Division-of-labor experiments** — how much lives in macos (Fortran)
vs MATLAB is an open question, resolved empirically here, on the
coronagraph examples, with measured numbers.  Each experiment is a
short MATLAB script + a dated note in this file:

- [ ] **E1 — MATLAB-side outer merit.**  Read `intensity(det)`,
      compute annular dark-zone contrast in MATLAB (port the
      `contrast.py` λ/D machinery from pymacos).  Measure per-eval
      wall time at nλ = 3 on the Phase-5 configuration; reproduce the
      ~3e-10 baseline.
- [ ] **E2 — DM control entirely in the outer loop.**  DM Zernike /
      Fourier modes as fmincon design vars, no inner loop, MATLAB
      merit from E1.  Find the mode count where cost becomes
      prohibitive — this sets the bar a Fortran DarkZone target must
      beat.
- [ ] **E3 — CALIB mechanics through diffraction.**  Exercise the
      existing CALIB targets on DM elements of `Rx_Coro.in`; time the
      objective evaluation when it must run the CPROPAGATE chain at
      nGridpts = 511.  (Absorbs the Sprint-0 CALIB-timing pre-flight.)
      Determines whether a Sprint-3 inner loop is seconds or hours.
- [ ] **E4 — λ-loop placement.**  Measure the load / trace / propagate
      cost split per wavelength to quantify what a future MACOS-side
      `spectral_run` amortization would actually buy over a MATLAB
      loop of `set_src_wvl` calls.

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

### Sprint 2A — `macos.design.Telescope` 2-mirror MVP

- [ ] New package `+macos/+design/` under `MACOS_resources/mmacos/`.
- [ ] Spec struct + pure-function resolve/emit (§3); component
      interface (§6.2) with Mirror / Mask / FocalPlane.
- [ ] Closed-form layout for Cass / RC / Gregorian / DK with the §5.2
      corrected, β-dependent forms; §5.3 raytrace validation tests.
- [ ] Full-precision `.in` emission + `ValidatePrescription` on every
      build; FD-survivability test (§1.3).
- [ ] `evaluate_` with canonical call sequence, ray-loss guard,
      worker-safety; internal [0,1] design-var normalization.
- [ ] Outer fmincon loop, nested λ × field (nλ=1 default for
      all-reflective); merit built-in + callback path.
- [ ] `describe()` with provenance; `save_spec`/`load_spec`;
      golden-spec → byte-identical-Rx regression.
- [ ] First result: M2 alignment minimizing FoV-averaged WFE vs a
      known-good Cass baseline (ZGD_test_files).

### Sprint 2B — N-mirror (TMA / 4-mirror / freeform)

- [ ] Generalized `add_mirror` + spacing resolution (single `'derive'`).
- [ ] Equal-power split default; `override('Mk','fnum',…)`.
- [ ] `set_surface` → SrfType 2 / 4 / 12 / 14 dispatch; design vars →
      MonCoef / FFZernCoef emission (sparse FFZernModes form).
- [ ] ZernTypeL dispatch-chain ELSE-with-error (Fortran; §6.3).
- [ ] First wide-FoV result: TMA, M2 monomial + M3 freeform-Zernike,
      band-and-FoV-averaged WFE over several arcmin.

### Sprint 2C — refractive + dispersive components

- [ ] Lens (2-surface expansion, glass names, group perturbation) and
      Grating (λ-dependent exit direction) per §6.
- [ ] Glass-catalog expansion + .agf converter (§6.4).
- [ ] λ_ref layout machinery + dispersion-aware detector sizing.
- [ ] Worked example: §7.4 spectrograph; chromatic merit exercises
      nλ>1 for real.

### Sprint 3 — `macos.design.Coronagraph` + DarkZone target

- [ ] `Coronagraph` class per §7.3; conjugate-recreating OAP solve.
- [ ] **New CALIB target type `DarkZone`** — annular image-plane
      integral (`inner_lamD`, `outer_lamD` or absolute radii),
      evaluated through the diffraction chain (per Sprint 0 timing).
- [ ] `macos_api_mod` wrapper `calib_set_target_darkzone(...)` +
      mmacos surface + tests.
- [ ] **Scope guard:** inner-loop DM control is MODAL (tens of
      Zernike/Fourier modes) — SLSQP over thousands of raw actuators
      via FD is not viable.  High-actuator-count EFC (Jacobian-based)
      belongs to FALCO integration (PLAN.md §2.3); DarkZone is scoped
      accordingly.
- [ ] Warm-started inner loop between outer iterates (§1.3.3).

### Sprint 4 — diagnostics + docs

- [ ] Full `check_clipping()` (worst margin per element; the
      pass/fail version ships in 2A) and `diagram()` (3-D side view,
      names, distances, beam radii).
- [ ] `plot_history`.
- [ ] All three worked examples in the manual; tag `design-layer-v1`.

### Sprint 5+ (deferred)

- Polychromatic CALIB (broadband inner-loop merit in one solve).
- `spectral_run(λ_list, …)` amortization convenience.
- Clearance-derived OAP off-axis distances; Talbot DSL keyword.
- Surrogate outer optimizers (Bayesian / SMAC) if eval cost demands.
- `export_proper()` → PROPER/FALCO prescription for end-to-end
  contrast (the system-model schema makes this an emitter, §3).
- Python port of the builder — struct translation + transliteration,
  parity = byte-identical Rx on golden specs.
- Segmented-primary integration (SegMirMaker as an aperture-shape
  provider for M1).

---

## 9. MACOS-side support work — when each piece is needed

| Sprint | macos / `macos_api_mod` change | Why |
|---|---|---|
| 0 | (none — experiments only) | design inputs |
| 1 | `set_src_wvl` / `set_src_fov` already present in mmacos (verification only); ray-loss per-category breakdown (Q2 — small).  E1–E4 experiments (no other Fortran — measurements scope Sprint 3) | FoV loop; §1.3 guard; macos↔MATLAB split |
| 2A | (none in Fortran) | builder + outer loop entirely MATLAB |
| 2B | ZernTypeL dispatch ELSE-with-error | silent no-op guard |
| 2C | catalog formula check (possibly parser extension); .agf converter (offline tool) | refractive components |
| 3 | `DarkZone` CALIB target + wrapper | coronagraph inner loop |
| 4 | (none) | docs + diagnostics |

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
| Q2 | Per-category ray-status getter in `macos_api_mod` exposing `RayStatus(:)` (+ optionally `RayFailElt(:)`, `RayFailMsg`) verbatim from `elt_mod`.  mmacos surface: `get_ray_status(N)` returning the integer-coded category per ray (`RayStat_OK/Obscured/Miss/Bracket/MaxIter/Undef`) — complements the existing binary `get_ray_info` (`ok_trace`, `ok_pass`).  Codegen Path A; ~30 min. | known-vignetting Rx → known per-category counts |
| Q3 | ZernTypeL dispatch ELSE-with-error (propsub.F / srtrace.F) | `ZernType= Noll` errors loudly, no silent no-op |
| Q4 | Glass catalog: formula-coverage audit → parser extension if needed → .agf converter + generated usual set | n(λ) vs published to 1e-6; one doublet vs CodeV |
| Q5 | Endurance test (500× load/trace, one session) → fix findings | bit-identical rmsWFE each iter; flat memory |
| Q6 | `Element= Apodizer` (independently motivated, PLAN.md §2.1 Thrust B) | participates in subsequent trace; PROPER cross-check |

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
   example (Sprint 2A: layout math, no DMs); spectrograph third.
6. **macos↔MATLAB division of labor is decided by measurement**
   (Sprint 1 E1–E4), not upfront.  Default: logic starts MATLAB-side;
   it moves into Fortran only for inner-loop speed (DarkZone) or to
   participate in propagation (Apodizer).
7. **nλ=1 default** keyed to *all-reflective*, not to merit type.
8. **Component as topology unit** from Sprint 2A.
9. **Modal inner-loop DM control**; actuator-level EFC → FALCO.

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
- Design package: `MACOS_resources/mmacos/+macos/+design/`.
- DarkZone target (Sprint 3) → `sls-dev`; fast-forward to `opt-dev`
  at release-worthy sprint ends.
- Tags: `design-sprint-0` … `design-sprint-4`, `design-layer-v1`.
