# Design Layer — Proposed Workflow

A high-level design surface on top of MACOS. The user expresses
*design intent* (topology, free parameters, merit function); a builder
generates a concrete MACOS prescription; MACOS handles ray trace and
inner-loop alignment; an outer optimizer (MATLAB `fmincon` or
SciPy `minimize`) explores the design space.

Drafted 2026-06; companion to `PLAN.md`.  Telescope is the first
worked example; coronagraph the second.  MATLAB-first because the
JPL user base is MATLAB-heavy; Python port is a future possibility
once the MATLAB API stabilizes.

---

## 1. Architecture: three tiers, two loops

**Tier 1 — MATLAB design layer.** Builder API + outer optimizer.
Owns design intent: topology, free vs. fixed parameters, merit
function, bandwidth, field-of-view.

**Tier 2 — Prescription artifact.** `.in` file emitted by the
builder each iteration.  Language-neutral handoff; remains readable
and reproducible on its own without the builder.

**Tier 3 — MACOS execution.** mmacos wrappers over `macos_api_mod`.
Per evaluation: load Rx, set wavelength, set chief-ray direction
(field point), run inner CALIB if there are in-instrument DOFs to
align, read intensity / WFE, return.

Two nested loops:

- **Inner (MACOS, existing strength).**  Fixed topology + dimensions;
  aligns DMs / rigid-body DOFs against an in-instrument target via
  the SLSQP back end.  Already in place via CALIB + the Phase 1
  wrappers (`calib_run`, `calib_set_var_elt`, `calib_set_target`).
- **Outer (MATLAB `fmincon`).**  Topology fixed; *dimensions are
  free* (OAP focal lengths, distances, mask radius, Lyot undersize,
  M2 despace, conic, etc.).  Each step regenerates geometry via the
  builder, calls the inner loop, evaluates the bandwidth-and-FoV-
  averaged merit, returns scalar.

```
┌─ outer (MATLAB fmincon, design params) ───────────────────┐
│   for each fmincon step:                                  │
│     builder.update(design_params)                         │
│     rx = builder.build()                                  │
│     merit = 0                                             │
│     for λ in bandwidth:                                   │
│       for field in field_points:                          │
│         m.load(rx); m.set_wavelength(λ); m.set_src_dir(field)│
│         ┌─ inner (MACOS SLSQP, if applicable) ──┐         │
│         │   m.calib(...)                       │         │
│         └────────────────────────────────────────┘         │
│         eval = m.intensity(detector) or m.rms_wfe(...)    │
│         merit += w(λ,field) * user_merit_fn(eval, λ, field)│
│     return merit                                          │
└───────────────────────────────────────────────────────────┘
```

Per-evaluation cost scales as nλ × nField.  3 wavelengths × 7 hex
field points = 21 traces per outer step.  At ~0.1–1 s per trace, a
60-step outer loop is a few minutes to ~20 minutes.  Cache the Rx
across the inner double-loop so only the chief-ray direction and
wavelength change between traces.

---

## 2. Builder API — Telescope (first example)

```matlab
t = macos.design.Telescope('aperture_diameter_m', 6.0);

% Topology
t.add_primary('M1', 'fnum_solve','from_fov_and_M2');
t.add_secondary('M2');
t.add_focal_plane('FP');

% Outer-loop design vars
t.declare_design_var('M2_despace_mm',   0.0,  'bounds',[-2 2]);
t.declare_design_var('M2_tilt_xy_mrad', [0 0],'bounds',[-1 1]);
t.declare_design_var('M2_conic',       -1.0,  'bounds',[-1.2 -0.8]);

% Bandwidth
t.set_bandwidth([500 700 900]*1e-9, 'weights',[1 1 1]);

% Field of view — two equivalent forms
t.set_field_points(deg2rad([0 0; 0 0.01; 0.01 0; -0.01 0; 0 -0.01]), ...
                   'weights',[1 0.5 0.5 0.5 0.5]);
% or
t.set_field_grid('shape','hexagonal', 'max_half_angle',0.01*pi/180, ...
                 'n_rings',2);

% Merit
t.set_outer_merit('RMSWFE_FoV_band_averaged');     % built-in
% or
t.set_outer_merit(@(wfe_stack, lam, field) ...);   % full callback

% Sanity check, optimize, visualize
report = t.check_clipping();
[opt, history] = t.optimize('algorithm','sqp', 'MaxIter',60);
t.diagram();
```

Telescope is the simpler topology: no DMs, no FPM, no Lyot, no inner
loop (no in-instrument DOFs to align).  This lets Sprint 2 produce a
working design loop with **no new Fortran code** — validating the
architecture before investing in coronagraph-specific MACOS features.

## 3. Builder API — Coronagraph (second example)

```matlab
c = macos.design.Coronagraph('input_focal_plane', ota_fp_elt, ...
                             'pupil_diameter_mm', 50);

c.add_OAP('OAP1', 'recreates','pupil', 'fills','DM1');
c.add_DM ('DM1',  'actuators', 64);
c.add_DM ('DM2',  'at_fraction_of','talbot', 'fraction', 0.25, 'ref','DM1');
c.add_OAP('OAP2', 'recreates','focus');
c.add_mask('FPM', 'shape','disk');
c.add_OAP('OAP3', 'mirrors','OAP1');
c.add_lyot('LS1');
c.add_OAP('OAP4', 'recreates','focus');
c.add_detector('SciCam');

c.declare_design_var('FPMRadius',    4.0,  'bounds',[2  8],  'unit','lambdaD');
c.declare_design_var('LyotUndersize',0.85, 'bounds',[0.7 0.95]);

c.set_bandwidth([550 600 650]*1e-9, 'weights',[1 1 1]);

% Inner-loop target (MACOS CALIB aligns DMs each evaluation)
c.set_inner_target('DarkZone', 'inner_lamD',7, 'outer_lamD',10);

% Outer-loop merit (user callback, sees intensity per λ × field)
c.set_outer_merit(@(I_stack, lam, field) my_avg_contrast(I_stack,7,10));

[opt, history] = c.optimize();
```

Coronagraph stresses every part of the pipeline: builder topology,
design vars in outer loop, DMs in inner loop, multi-λ averaging,
image-plane merit, DarkZone target.  The on-axis case has degenerate
FoV; the FoV machinery from the telescope sprint still applies for
off-axis science targets at known separations.

---

## 4. Sprint sequence

### Sprint 1 — wrapper parity in mmacos (no Fortran)
- Port `intensity()` and `complex_field()` from pymacos to mmacos
  (the buffers already live in `macos_api_mod`; just need mex wrapping).
- Expose `set_wavelength()` (likely already there as `set_src_lambda`
  or `MOD WAVELEN=`; verify and document).
- Expose `set_src_dir(θx, θy)` — chief-ray direction setter for the
  field-point loop.  Probably reachable via `set_src_csys`'s field-
  angle path; verify and wrap.
- Verify CALIB wrappers handle DM-element variables on `Rx_Coro.in`.
- Tests: hand-written mmacos test that loads `SegDemo3.in`, sets a
  field point, traces, reads WFE; expected values from pymacos.

### Sprint 2 — `macos.design.Telescope` end-to-end MVP
- New MATLAB package `+macos/+design/` under `MACOS_resources/mmacos/`.
- `Telescope` class with the API above; closed-form layout math;
  `.in` emission.
- MATLAB outer loop with nested λ × field double loop.
- Single-λ first, then enable bandwidth (one extra `for` loop).
- Merit: `RMSWFE_FoV_band_averaged` built-in plus user-callback path.
- First real result: M2 alignment that minimizes band-and-FoV-averaged
  WFE.  Compare to known-good baseline (Cass-like Rx in ZGD_test_files).

### Sprint 3 — `macos.design.Coronagraph` + DarkZone target type
- `Coronagraph` class with the API above.
- **New MACOS target type:** `DarkZone` — annular integral over an
  image-plane region.  Takes `inner_lamD`, `outer_lamD` (or absolute
  pixel radii) as parameters.  Adds to the existing CALIB target
  enum alongside WFE / Zern / Beam.
- `macos_api_mod` wrapper: `calib_set_target_darkzone(inner, outer, units)`.
- mmacos surface: `m.calib_set_target_darkzone(...)`.
- Inner CALIB now optimizes DMs against a real coronagraph merit
  rather than relying on outer-loop loops over individual DM modes.
- Bandwidth-and-FoV cost is real here; budget accordingly.

### Sprint 4 — diagnostics + docs
- `check_clipping()` — traverse layout, walk the beam through each
  element, report worst margin per element.
- `diagram()` — MATLAB plot/patch 3-D side view, annotated with
  element names, distances, beam radii.
- Both examples written up in the manual.
- Tag the snapshot (e.g., `design-layer-v1`).

### Sprint 5+ (deferred)
- Polychromatic CALIB — inner loop optimizes DMs against a multi-λ
  merit in one solve (broadband EFC), not just outer-loop averaging.
- `m.spectral_run(λ_list, weights, target)` convenience to amortize
  Rx-load cost across wavelengths when geometry doesn't change.
- Surrogate-modeled outer optimizers (Bayesian opt / SMAC) if
  `fmincon` eval cost gets prohibitive.
- Python port of the builder — mechanical translation from MATLAB
  once the API stabilizes.

---

## 5. MACOS-side support work — when each piece is needed

| Sprint | MACOS / `macos_api_mod` change | Why |
|---|---|---|
| 1 | `set_src_dir` exposed in mmacos | FoV loop |
| 1 | `intensity` / `complex_field` ported pymacos → mmacos | Image-plane merits |
| 2 | (none in Fortran) | Builder + outer loop entirely MATLAB |
| 3 | New `DarkZone` CALIB target type | Coronagraph inner loop |
| 3 | `calib_set_target_darkzone(inner, outer, units)` wrapper | mmacos surface |
| 4 | (none) | Docs + diagnostics |

Sprint 1 is mmacos-only — fast, no Fortran risk.  Sprint 2 produces
a working telescope-design loop with no new Fortran, validating the
architecture before Sprint 3 invests in `DarkZone`.

---

## 6. Decisions made

1. **Single-λ vs. bandwidth-aware first build.**  Single-λ in
   Sprint 1; turn on bandwidth in Sprint 2.  Keeps the integration
   test sharp.
2. **DarkZone target — MACOS or MATLAB merit?**  Both.  MACOS-side
   DarkZone for the inner loop (so DM CALIB actually optimizes the
   right thing); MATLAB merit for the outer loop (flexibility >
   speed at that tier).
3. **Outer optimizer.**  `fmincon` for now.  Revisit SLSQP-as-a-
   service only if outer-loop merit eval cost blows up.
4. **Order of examples.**  Telescope first (simpler, validates
   plumbing without DM CALIB), coronagraph second.  Telescope
   performance metric is over a **defined band AND a defined
   field of view.**

## 7. Decisions still open

1. **Field-grid default form.**  Hex grid (`n_rings=1` = 7 points)
   vs. explicit list as natural input.  Probably both; hex as default.
2. **FoV weighting model.**  Equal-weight default; Strehl-weighted
   built-in option.  Worst-case merit vs. averaged — both available.
3. **Bandwidth specification.**  Discrete λ list as primitive;
   spec-band sampling (`'500-900 nm @ 5 samples'`) as friendlier
   wrapper.
4. **`fnum_solve` semantics.**  How aggressive should the closed-
   form layout solve be in the Telescope builder?  Need to balance
   "user wires up M1 explicitly" vs. "builder derives from desired
   FoV + M2 placement."  Decide after Sprint 2 first cut.

---

## 8. Branch routing

- All work lands on **`sls-dev`** (new features) per the revised
  branch model in `PLAN.md` §9.
- MATLAB design package lives in `MACOS_resources/mmacos/+macos/+design/`
  on the companion `sls-dev` branch.
- DarkZone target type (Sprint 3) is a real new feature → `sls-dev`.
- At the end of each sprint, if release-worthy, fast-forward-merge
  `sls-dev` → `opt-dev`.
- Tag sprint completions (`sprint-1-wrappers`, `sprint-2-telescope`,
  `sprint-3-darkzone`, `sprint-4-docs`) for rollback points.
