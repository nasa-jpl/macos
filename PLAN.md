# MACOS Development Plan

Working document. Tasks grouped by thrust; checkboxes track state. Items not
yet started have empty `[ ]`; completed `[x]`. Discrete subtasks under each
task.

> **Related plans.**  This file owns **engine + wrapper** work
> (Fortran, `macos_api_mod`, pymacos/mmacos surface).  The MATLAB
> **design layer** (builder, importer, `vary`/`evaluate`/
> `sensitivities`, metrology orchestration) is planned in
> `PLAN_DESIGN_LAYER.md`.  Ownership rule: an item has exactly one
> home and one checkbox; the other plan **cross-references by
> section**, never duplicates the checkbox.  Engine items the design
> layer depends on stay here; the design plan points at them.

---

## 0. Hygiene / quick fixes

- [x] Zernike trace-dispatch ELSE branches in `propsub.F`, `srtrace.F`, `tracesub.F`; bring all three into parity, modernize to `ZernType_*` constants. Test: any Rx with `ZernType= Noll` and a non-zero coefficient should produce non-zero OPD response.
- [x] `fk.c` unconditionally in `libsmacos.a` (drop `BUILD_GMI` gate). Makes `makems.sh + makegmi.sh` work standalone.
- [x] Build scripts auto-bootstrap bundled readline on first clone. `makems.sh`, `makegfortran.sh`, `makeall.sh` detect missing `libreadline.a` and run `./configure -q && make` in `macos_f90/readline-8.2/`.
- [x] `dopt_init_vars` stop clobbering meaningful defaults of `nitrs_dopt`, `dopt_tol`, `SvdSvCut`.
- [x] Commit `ZGD_test_files/` consolidated corpus (joint-dev pickups + local additions, drop `macos_param.txt` symlink).
- [x] Drop marker-based gate from `run_regression.sh`; trust MATLAB exit code now that `-reentrancy=none` is on the ifx mex link.
- [x] `define_local_csys` follow-up from the `develop_STOP` branch head — already on opt-dev + release-candidate via the IRIS STOP-rewrite merge; only `develop_STOP`-vs-rest delta is author-credit comments.
- [x] Audit `dopt_init_vars` for any other meaningful-default-clobbered-by-zeros collisions beyond the three found today.  Closed 2026-06-04 (commit 0ee4b23): `OptTgtElt` and `OptAlg` now initialize to the sentinel values the LOAD path expects (-1 and `NonLin`), not 0.  The pre-existing zeroing was benign when `dopt_init_vars` only ran during init, but with `macos_realloc=.true.` now firing on every size-change init (the §0 fix), a re-init after a load would clobber these, disabling the `OptTgtElt < 0 → nElt-1` resolution and landing OptAlg outside its enum.
- [x] Quiet the `MBFile6: Unidentified string` warning when a parser hits a target-specific keyword (e.g. `OptBeamPos=`) under a non-matching `OptTarget`. Closed 2026-06-04 (commit 0ee4b23): added a 7-char `OptBeam` catch-all in `msmacosio.inc` immediately after the explicit BEAM_TARGET-gated handlers.  Keys recognized as valid syntax but silently skipped when OptTarget doesn't match.  Same pattern can be replicated for other target-specific keys if more surface.
- [x] **Heap corruption mistakenly attributed to engine — actually a codegen dim-arg ordering bug.**  Closed 2026-06-03 (commit 5890710 mmacos).
  The "trace+perturb on opt_example crashes" pattern was NOT a
  macos engine bug.  Root cause: the mmacos mex codegen's
  `do_calib_set_var_elt` and `do_calib_set_target` dispatchers
  allocated their variable-length array args BEFORE reading the
  size scalar from prhs.  `n_zern` was uninitialized garbage at
  the `allocate(zern_modes(max(n_zern,1)))` site, and the
  subsequent `mxCopyPtrToReal8` read past the source array,
  corrupting the next heap allocation.  Crash surfaced on the
  next mex call (`load + calib_clear_var_elts + calib_set_var_elt`
  was the minimal repro; the corruption tripped after `clear`
  printed ok but the SIGSEGV was actually inside set_var_elt's
  prologue).
  The codegen intended to reorder reads so dim-scalar prhs args
  come first, but its dim-name extractor used regex
  `^[A-Za-z_]\w*$` -- only bare identifiers.  Dim expressions
  like `max(n_zern, 1)` slipped through, so n_zern never got
  added to dim_names and was read in declaration order (last),
  AFTER zern_modes was allocated using it.
  Fix: extract every bare identifier from each dim expression via
  `re.findall(r'[A-Za-z_]\w*', d)` and exclude Fortran intrinsics
  (max/min/size/kind/len).  Audited all 89 emitted dispatchers
  with a verifier script; only the two CALIB setters had the
  bug.  tCalib (9 tests) restored and green.
  Diagnosis path: matlab -Dgdb with `handle SIGSEGV stop print`
  caught the exact crashing function (`do_calib_set_var_elt`)
  before the malloc arena even noticed the corruption.  Future
  similar bugs: this is the fastest diagnostic path -- skip
  valgrind, go straight to gdb on matlab.
  The ORIGINAL §0 bug (model_size transitions triggering malloc
  aborts) was separately fixed -- see next entry.
- [x] **Original §0 bug: `init()` re-init heap corruption on
  model_size transitions.**  Closed 2026-06-03 (opt-dev commit
  e2e8bf6, release-candidate cherry-pick 1d54dd9).  Two
  compounding causes:
  (1) `macos_api_mod.init()` updated `curr_model_size` and called
  `macos_init_all(new_size)` but never set
  `macos_realloc = .true.` -- so SMACOS's module-saved scratch
  buffers (L1, R1, R2, D2, DV1, DV2, CD1, CD2, DrawEltVec,
  DrawRayVec, PertVec, DWF in smacos_vars_mod) stayed at the
  OLD size on the next SMACOS dispatch.
  (2) `sunsub.F`'s NR_FFT-branch `DFOURN` allocated its scratch
  `DATA(:)` ONCE on first call (gated by `first_entry` SAVE'd
  LOGICAL) and never grew it.  After init(new_larger_size),
  `DFOURN` was called with bigger NN(:), `SzData = 2*NN(1)*NN(2)`
  exceeded the existing allocation, and the inner FFT butterfly
  wrote past the end of `DATA` -- crashing in the next process
  malloc.
  Diagnosis path: pymacos repro walked
  `init(128) -> intensity -> init(512) -> intensity ->
  init(1024) -> intensity`; gdb backtrace placed the SIGSEGV
  inside `dswap2_ <- cpropagate_ <- int_cmd_`; tracing
  `macos_realloc` + `first_enter` through SMACOS via temporary
  printfs confirmed the SMACOS scratch buffers reallocated
  (after the macos_api_mod fix), so the remaining write was
  in a routine NOT in smacos_vars_mod's realloc list -- which
  led to `DFOURN`'s `first_entry` gate.
  Fix lands the macos_api_mod half on opt-dev only (the file
  doesn't exist on release-candidate / main yet); the sunsub.F
  half lands on both branches and is universally correct for
  GMI / interactive macos / smacos_dvr too.
  Followup cleanup: the per-size-group split in
  `run_mmacos_tests.sh` and the per-phase pytest split in
  `run_proper_tests.sh` can both be collapsed back to one
  process now -- left for a separate commit so the fix can be
  reviewed independently.  Also bring `main` along when it
  next syncs.
- [~] Port worth-keeping IEEE / accuracy patterns from
  `docs/Archive/dev_optimization_surfsub/` into opt-dev's surfsub.
  **Partial close 2026-06-04 (commit 0ee4b23):** Pattern 4
  (near-tangent fallback) ported into `surfsub.F ConSrf` -- when
  `k2 = b² - 4ac < TOL_TANGENT·b²` the standard two-root selection
  amplifies noise into the spread between near-identical roots;
  collapse to the vertex formula `L = -b/(2a)` instead.
  `TOL_TANGENT` (~1d-14) already in `Constants` from the earlier
  pull.  Patterns 2 (discriminant pre-check) and 3 (cancellation-
  aware quadratic-root selection) already present in current
  `ConSrf`.  **Still open:** pattern 1 (IEEE flag clear+check at
  `ConSrf` and `MonSrf` boundaries) -- the most invasive of the
  four (multiple sites, requires use of `ieee_arithmetic` and
  `ieee_exceptions`); land separately when there's a concrete
  silent-NaN report to test against.
  The five archived patches (Sigrist, 2026-01) never merged; the
  branch has been deleted from origin.  README.md in that directory
  enumerates four bounded patterns worth lifting by hand:
  (1) `ieee_arithmetic`/`ieee_exceptions` flag clear+check at ConSrf
  and MonSrf boundaries — catches silent NaN propagation;
  (2) discriminant pre-check before `SQRT(k2)` — distinct failure
  mode from the bracket/max-iter cases the per-ray status work
  already handles;
  (3) cancellation-aware quadratic-root selection (conjugate form
  when `|b| ≈ √k2`);
  (4) near-tangent fallback `L=-b/(2a)` when `k2 < TOL_TANGENT·b²`.
  Pull `TOL_TANGENT` / `TOL_CANCEL` / `TOL_GEOM` from patch 0003's
  Constants prelude.  Skip the PAUSE-removal patch (already done on
  opt-dev via 51bbf1f), the `SolveConicIntersection` extraction
  (conflicts with FreeForm dispatch), and the GridSrf rewrite
  (overlaps too much with the jGridSrf + FreeForm grid-component
  work).
- [x] Renormalize `psiElt` after the `Q·psi` rotation in `CPERTURB_PROG` (funcsub.F:349-350).  Closed 2026-06-04 (commit 0ee4b23): one-line normalization after the matrix multiply.  `sin²(θ) + cos²(θ) ≠ 1` exactly in IEEE 754 for some θ (1e-6, 3e-5 notably) used to leave psi off by 1 ULP and drift slowly under repeated perturbs, producing a ~3e-14 OPD round-trip residual.  Regression probe at `MACOS_resources/mmacos/tests/tPerturbRoundtrip.m` was written defensively to allow both pre-fix (within 4*eps) and post-fix behaviour; now post-fix.

---

## 1. Thrust A — Test + document, sensitivity workflow as headline

### 1.1 Test infrastructure

- [ ] macos-CLI regression harness, modeled on `MACOS_resources/GMI/regression/run_regression.sh` but driven by `smacos_dvr` or by piping `.jou` to `macos`. Compares `OPD` / `SPOT` / `WFE` to reference snapshots at absolute tolerance.
  - [ ] Choose harness language (bash + Fortran driver) and layout
  - [ ] Define reference-snapshot format and bootstrap mechanism
  - [ ] Add tolerance comparison utility analogous to GMI's `compare_within.m`
  - [ ] Add per-test pass/fail reporting
  - [ ] Hook into CI (eventually)
- [ ] Wire `MACOS_resources/GMI/test_ff/` corpus (21 files) into a regression harness
  - [ ] Choose canonical subset
  - [ ] Generate references for each
  - [ ] Add to GMI regression suite alongside existing 6 tests
- [ ] Wire `macos/ZGD_test_files/` corpus (48 files) into the macos-CLI harness
  - [ ] Audit which prescriptions exercise distinct code paths
  - [ ] Build journal + reference for each canonical one
  - [ ] Treat realistic-system files (`afta_v2`, `hcitv2_lou`, `access121310_v8test`, `iris_dp_v14`, `iris_dp_ZGD`) as end-to-end integration tests
- [ ] Pymacos test expansion in `MACOS_resources/pymacos/tests/`
  - [ ] Optimization paths (LM, multi-field once available)
  - [ ] Polarization (once elements land)
  - [ ] Sensitivity-matrix end-to-end
  - [ ] Group perturbations
- [ ] Decide cross-repo test-corpus consolidation strategy
  - [ ] Identify duplicate prescriptions between `macos/ZGD_test_files/` and `MACOS_resources/GMI/test_ff/`
  - [ ] Decide: keep deliberate duplication, factor into shared location, or split by purpose
  - [ ] Some coronagraph Rxes will remain private — they belong in `MACOS_resources` only
  - [ ] Manual examples and reader exercises must be in `macos` directly (public-only consumers)

### 1.2 Function-by-function exercise pass

Each function → regression test → manual entry → fix what surfaces. Order roughly by user-visibility.

- [ ] Surface types
  - [ ] Flat (SrfType=1)
  - [ ] Conic (SrfType=2)
  - [ ] Aspheric (SrfType=3)
  - [ ] Monomial (SrfType=4)
  - [ ] Interpolated (SrfType=5)
  - [ ] Anamorphic (SrfType=6)
  - [ ] UserDefined (SrfType=7)
  - [ ] Zernike (SrfType=8)
  - [ ] GridData (SrfType=9)
  - [ ] Toric (SrfType=10)
  - [ ] AsGrData (SrfType=11)
  - [ ] MonGrData (SrfType=12)
  - [ ] ZrnGrData (SrfType=13)
  - [ ] FreeForm (SrfType=14)
- [ ] Element types
  - [ ] Reflector
  - [ ] FocalPlane
  - [ ] Reference
  - [ ] HOE
  - [ ] Grating
  - [ ] Refractor
  - [ ] Obscuring
  - [ ] Return
  - [ ] NSRefractor
  - [ ] LensArray
  - [ ] Segment
  - [ ] NSReflector
  - [ ] TrGrating
  - [ ] RfPolarizer (existing, not yet exercised)
  - [ ] TrPolarizer (existing, not yet exercised)
  - [ ] CGHNullPlate
  - [ ] DoeTrGrating
- [ ] Commands
  - [ ] OLD
  - [ ] MODIFY
  - [ ] TRACE
  - [ ] OPD
  - [ ] SPOT
  - [ ] INT
  - [ ] PIX
  - [ ] STOP (LM-based ChiefRayAiming)
  - [ ] FEX
  - [ ] SXP
  - [ ] PERTURB
  - [ ] GPERTURB
  - [ ] CPERTURB
  - [ ] LPERTURB
  - [ ] METcalc (metrology beam lengths; see §4.5)
  - [ ] NOMINAL / SetToNominalSettings
  - [ ] CALIB (covered by Thrust C regression matrix)
  - [ ] VARS / FOVS / WLENS (configuration inspection)
  - [ ] VALIDATE
  - [ ] JOURNAL / JOU (recording vs playback)
  - [ ] EXE / EXECUTE / LQW (journal playback)
  - [ ] LOG (intensity log10 display, not transcript)
  - [ ] IMGMODE (negative/positive polarity)
  - [ ] PGP / NEWPAGE
  - [ ] COMPOSE
  - [ ] SAVE / RESTORE / WFITS
- [ ] Channel-input forms
  - [ ] Dense ZernCoef
  - [ ] Sparse `nZernCoef` + `ZernModes` + `ZernCoef`
  - [ ] FFZernCoef (FreeForm composite)
  - [ ] MonZernCoef (FreeForm composite)
  - [ ] ZernType variants (all 11 named conventions)
  - [ ] AsphCoef (dense + sparse via nAsphCoef)

### 1.3 Sensitivity-matrix workflow (Thrust A headline)

> **Design-layer convergence.**  The design layer's first deliverable
> — a sensitivity table on an imported (CodeV-converted) Rx via
> `from_rx` / `vary` / `sensitivities` (PLAN_DESIGN_LAYER Sprint 2A-i)
> — **is** this Thrust A headline recipe at the package level.  Its
> required runnable worked example (`example_sensitivities_from_rx.m`)
> doubles as the headline example here; the underlying engine is the
> bitwise-verified Phase 7 `dw_dx` channel machinery (§5.4).  Keep the
> CLI/GMI/pymacos recipes below as the engine-level reference; the
> `+macos` design-layer surface is the recommended user entry point
> (see §11.7).

- [ ] Document the workflow as a first-class user recipe
  - [ ] GMI path: `param.zernSrf` / `param.dmSrf` / `param.rbSrf` / `param.gridSrf` channels; multi-call pattern; output `OPD` / `WFE` / `SPOT`
  - [ ] Pymacos path: `dw_dx` centered-difference; group perturbations; SXP / SRS / ORS wrappers; source-vs-group cross-check pattern
  - [ ] Forthcoming mmacos path (placeholder until §5 lands)
- [ ] Sensitivity-matrix regression test
  - [ ] Build a small reference system (Cassegrain or similar)
  - [ ] Compute sensitivity matrix via both paths
  - [ ] Freeze matrix as reference
  - [ ] Test for bit-stable reproduction across macos changes
- [ ] Worked end-to-end example in manual
  - [ ] Cassegrain or 3-mirror system
  - [ ] 6-DOF rigid-body sensitivity matrix
  - [ ] Demonstrate WFE residual under perturbation + correction

---

## 2. Thrust B — Coronagraph beyond PROPER, accessible to PROPER users

### 2.1 Capabilities beyond PROPER's reach

- [ ] Hard-edged amplitude apodizer
  - [ ] New `Element= Apodizer` prescription element
  - [ ] Parametric amplitude profiles (top-hat, Gaussian, prolate-spheroidal)
  - [ ] Participates in subsequent ray trace (closes `pymacos.apodize` post-hoc-only limitation)
- [ ] Complex grid-data element
  - [ ] New SrfType for complex amplitude + phase from FITS file
  - [ ] Supports characterized metasurfaces, AOMs with phase imperfections, measured apodizers
  - [ ] Real + imag from paired FITS, or single complex FITS
- [ ] Mask library expansion
  - [ ] Vortex masks (topological-charge parameter)
  - [ ] Half-plane / knife-edge masks
  - [ ] Complex amplitude masks via FITS input
- [ ] Polarizing elements (per-ray Jones-matrix plumbing already exists; Fourier vector diffraction already exists)
  - [ ] Quarter-wave plate (QWP), arbitrary fast-axis orientation
  - [ ] Half-wave plate (HWP)
  - [ ] Linear polarizer, configurable angle
  - [ ] Vortex retarder, topological-charge parameter
  - [ ] Polarizing beam splitter (PBS) — design sketch first; ray-branching is the hardest piece
- [ ] Polarization metrics
  - [ ] Degree of polarization (DOP) at arbitrary surface
  - [ ] Stokes vector (S0, S1, S2, S3) at arbitrary surface
  - [ ] Polarization-resolved PSF at focal plane
  - [ ] Diattenuation / retardance maps
- [ ] Multi-field-point metrics (also blocks Thrust C completion)
  - **Ownership update (2026-06-11):** the **outer-loop** field-of-view
    sweep (evaluate a merit over a fan of source positions, aggregate
    band-and-FoV) is now owned by the **design layer** — its λ×field
    loop drives `set_src_fov` per field point with no engine change
    (PLAN_DESIGN_LAYER §1.1, §4.2).  What remains engine-side here is
    narrower: **multi-field awareness inside the inner CALIB solve**
    (a single optimization that scores against several fields at once),
    and that is **gated on E3/E4** (PLAN_DESIGN_LAYER §8 Sprint 1) —
    don't build it until those measurements show the inner loop needs
    it.  The subtasks below are scoped to that inner-CALIB case.
  - [ ] Extend `design_optim_mod` to evaluate at a fan of source positions (inner-CALIB multi-field; gated on E3/E4)
  - [ ] Aggregation modes: RMS-of-RMS, worst-field-point, weighted sum
  - [ ] GMI channel for multi-field input
  - [ ] Pymacos wrapper
  - [ ] Regression: 3-mirror system, optimize for 5-field-point average WFE

### 2.2 Accessibility to PROPER users

- [ ] `pymacos.proper_compat` Python module
  - [ ] Inventory FALCO-and-PROPER API surface (`prop_begin`, `prop_circular_aperture`, `prop_lens`, `prop_propagate`, `prop_dm`, `prop_pixellate`, etc.)
  - [ ] Implement facade calling through to pymacos primitives
  - [ ] Cover the subset used by FALCO testbed scripts
- [ ] MATLAB PROPER-compat layer over mmacos (parallel to Python; depends on §5)
- [ ] PROPER → macos migration tutorial
  - [ ] Pick a representative PROPER example (e.g. Krist's tutorials)
  - [ ] Convert to macos prescription side-by-side
  - [ ] Highlight each `prop_*` call and its macos equivalent
- [ ] Contrast-curve scoring API
  - [ ] Generalize `tests/proper_compare/contrast.py` machinery into `pymacos.contrast`
  - [ ] Load result + plot vs λ/D + compare against design target

### 2.3 FALCO integration

- [ ] Confirm FALCO current maintained-codebase status
  - [ ] MATLAB-primary vs Python-secondary state
  - [ ] Last-release vs ongoing-internal cadence
  - [ ] Willingness to accept patches / cooperate on a macos backend
- [ ] Inventory FALCO's PROPER usage
  - [ ] Grep FALCO sources for `prop_*` calls
  - [ ] Identify the subset that needs compat-layer coverage
  - [ ] Tightens §2.2 scope from "all of PROPER" to "what FALCO needs"
- [ ] Integration design doc
  - [ ] Pick path: (i) `pymacos.proper_compat` if FALCO-Python is the target; (ii) `mmacos` + MATLAB PROPER-compat if FALCO-MATLAB is the target
  - [ ] Sizing the API surface in either direction
- [ ] End-to-end FALCO + macos demo
  - [ ] Take one FALCO EFC run with PROPER as the model
  - [ ] Swap in macos via the compat layer
  - [ ] Confirm same dark-zone contrast achieved within numerical tolerance

### 2.4 Sub-wavelength / FDTD coupling (out-of-scope as embedded; document handoff)

- [ ] Identify which coronagraph features genuinely need FDTD
  - [ ] Roughened mask edges
  - [ ] Sub-10λ pinholes
  - [ ] Metasurfaces
  - [ ] Form-birefringent stops
- [ ] Pymacos handoff pattern
  - [ ] `m.export_field(srf)` — write complex field + spatial sampling at any element
  - [ ] `m.import_field(srf, field, dx)` — read returned field back
  - [ ] Worked example with MEEP on a small case (FPM throat)
- [ ] Document the workflow as an appendix; don't promise embedded FDTD

---

## 3. Thrust C — Optimizer robustification

### 3.1 Target-type fixes

- [ ] `WFE_TARGET` default to flat reference
  - [ ] Treat missing `OptTgtWF=` as `objfun_nom=0`, not a fatal file-open
  - [ ] Preserve file-load path when one is specified
- [ ] **Implement or remove `WFE_ZMODE_TARGET` and `OPL_TARGET` — now urgent.**
  As of §3.5 (CALIB wrappers), pymacos and mmacos users can pick these
  targets programmatically (`m.calib_set_target('WFE_ZMODE', modes=...)`)
  and the wrapper accepts the call without complaint.  Under the hood
  the optimizer either no-ops or fails silently.  Before §3.5 these
  paths were buried in `.jou` flows — now they're a documented part
  of the binding API.  Resolve before §3.4 regression matrix lands.
  - [ ] Either wire up `objfun_nom` branches for both
  - [ ] Or remove from the enum + parser so users can't silently pick a non-functional mode
  - [ ] Either way, document the resolution in pymacos's `m.calib_set_target` docstring and mmacos's `+macos/calib_set_target.m` help text
- [ ] Diagnose `SPOT_TARGET` reporting 0 iterations even with `dopt_init_vars` fixed
  - [ ] Find the early-exit gate
  - [ ] Fix and add regression case

### 3.2 Parser / Rx-keyword fixes

- [ ] Add `OptTol` Rx keyword
  - [ ] Mirror `OptMaxItrs` plumbing
  - [ ] Default 1e-12 (current source default after dopt fix)
- [ ] Fix `OptMxItrs` / `OptMaxItrs` keyword-prefix mismatch
  - [ ] Existing parser matches `'OptMaxItrs', 10` chars
  - [ ] Common typo `OptMxItrs` (9 chars) silently drops the user override
  - [ ] Either accept `OptMxItrs` as alias or match shorter prefix
- [ ] Standardize `OptAsph` vs `OptAsphCoef` keyword naming
  - [ ] Parser matches 7-char prefix `OptAsph` so both work today
  - [ ] Document the canonical form
- [ ] Cosmetic: silence `MBFile6: Unidentified string` warning when target-specific keywords appear under non-matching `OptTarget` (e.g. leftover `OptBeamPos=` under `OptTarget=SPOT`)

### 3.3 Algorithm / control-flow fixes

- [ ] Move convergence check to every iteration
  - [ ] Current behavior: fires only at `itr ∈ {0.21·tot_itrs, 0.5·tot_itrs, 0.71·tot_itrs}`
  - [ ] For `tot_itrs ≤ 4` the check never runs
- [ ] Consolidate two `nls_optim_dvr` call sites
  - [ ] `macos_cmd_loop.inc:555` nominal-sensitivity path
  - [ ] `macos_cmd_loop.inc:3630` main path
  - [ ] Same setup logic, two copies
- [ ] `quit` after `CALIB` exits cleanly
  - [ ] Currently save prompt drops back to `MACOS>` after answering
  - [ ] Needs second exit; `EXE` journals therefore can't terminate clean
- [ ] Audit `dopt_init_vars` for any other meaningful-default-clobbered-by-zeros collisions

### 3.4 Regression test matrix

8-12 new tests covering every `OptTarget` × DOF-type combination.
**Re-scoped (2026-06-02):** with §3.5 (CALIB wrappers) landed, these
no longer need to be `.in` / `.jou` file pairs — pymacos's
`tests/test_calib.py` and mmacos's future `tCalib` can express them
as pytest / matlab.unittest methods that drive the optimizer
programmatically.  Faster to author, easier to maintain, CI-ready.
The `opt_example.in` style fixtures stay (as the loadable Rxes
under test) but the test logic moves to Python / MATLAB.

|  | rigid-body | Zernike | Aspheric | mixed |
|---|---|---|---|---|
| BEAM | exists (`opt_example`) | new | new | new |
| SPOT | new | new | new | new |
| WFE (flat reference) | new | new | new | new |
| WFE_ZMODE | (only if implemented) | new | new | new |

- [ ] BEAM × rigid-body — done (`opt_example`)
- [ ] BEAM × Zernike
- [ ] BEAM × Aspheric
- [ ] BEAM × mixed
- [ ] SPOT × rigid-body
- [ ] SPOT × Zernike
- [ ] SPOT × Aspheric
- [ ] SPOT × mixed
- [ ] WFE × rigid-body (flat reference)
- [ ] WFE × Zernike (flat reference)
- [ ] WFE × Aspheric (flat reference)
- [ ] WFE × mixed (flat reference)

### 3.5 Optimizer wrappers (pymacos + mmacos) — Phase 1

Expose the CALIB optimizer through the language-neutral wrapper
layer so Python / MATLAB users can drive it without `.jou` files.
Single-pass thin wrap first; broader features (FOV list, multi-
wavelength, NPSOL constraints, beam-target sub-config) layered
on as needed.

- [x] **Phase 1a — bare `m.calib()` invocation** (commits d9a874c
  macos / 3ec2c28 MACOS_resources)
  - [x] `macos_api_mod` adds `calib_run` (returns rtn_flag +
    per-FOV/wavelength `old_wfe`/`new_wfe`) and `calib_buffer_dims`
    (reports max_fov × max_wl from dopt_mod so callers can size
    output buffers correctly)
  - [x] pymacos's `m.calib() -> dict` wraps both via the f2py
    forwarder.  3 pytest cases against `opt_example.in`: baseline,
    perturbed (1 mrad on Elt 1), schema check.  All pass.
  - [x] Full pymacos suite still 6604/6604; PROPER suite green.
- [x] **Phase 1b — programmatic setters** (commits bce9dab macos /
  d62b242 MACOS_resources)
  - [x] `calib_clear_var_elts`, `calib_set_var_elt(iElt, dofs[8],
    zern_modes[])`, `calib_set_iter`, `calib_set_tol`,
    `calib_set_target(t, wf_zern_modes[])` -- programmatic AVAR /
    MVAR / DVAR / SETC equivalents.
  - [x] Pymacos wrapper accepts both 8-int positional masks AND
    DOF name lists (TIP/TILT/CLOCK/DX/DY/PIST/ROC/CONIC); target
    accepts integer enum AND name strings (WFE/WFE_ZMODE/ZWF/
    BEAM/SPOT/OPL).  6 new pytest cases (total 9).
  - [x] Pymacos suite 6610/6610.
- [x] **Phase 1c — port to mmacos** (one commit on MACOS_resources)
  - [x] All 7 calls in mmacos_gen_cmds.txt; codegen emits clean
    dispatchers (with 3 codegen fixes: paren-aware namelist
    splitter, decl-section line-continuation joiner, dim-symbol
    use_only_map consultation).
  - [x] `+macos/calib.m`, `+macos/calib_clear_var_elts.m`,
    `+macos/calib_set_var_elt.m`, `+macos/calib_set_iter.m`,
    `+macos/calib_set_tol.m`, `+macos/calib_set_target.m`.
  - [x] `macos.Session` delegators added.
  - [x] `tCalib` matlab.unittest class (9 tests) -- restored
    2026-06-03 after the codegen dim-arg ordering fix (commit
    5890710) closed the "heap-corruption" attribution.  All
    9 tests pass: 3 baseline (load+perturb+calib), 6 setters.
- [ ] **Phase 1d — FOV / wavelength / NPSOL setters** (deferred)
  - [ ] `calib_add_fov(dir[3], pos[3], weight)` /
    `calib_clear_fovs()` -- programmatic AFOV.
  - [ ] `calib_set_wavelens([wl1, wl2, ...])`.
  - [ ] `calib_set_constraints(...)` -- NPSOL path (opt-dev only;
    ifdef-gated).
  - [ ] `calib_set_beam_target(...)` -- sub-config for BEAM_TARGET
    (beamPosElt, ifOptBeamPos, beamRefRayElt, etc.).
- [ ] **Phase 1e — introspection** (deferred)
  - [ ] `calib_get_var_elts() -> dict` -- list current variable
    elements with their DOFs + Zernike modes.
  - [ ] `calib_get_config() -> dict` -- snapshot of target, iter,
    tol, FOV list, wavelength list.

---

## 4. New core capabilities

### 4.1 Generalize FreeForm to non-conic bases

- [ ] Refactor `FreeFormSrf` so the conic-base assumption is replaced by a `BaseSrfEval` interface
- [ ] Extend composition to Toric, Anamorphic, Aspheric, UserDefined bases
- [ ] Per-base regression tests confirming current FreeForm behavior bit-stable
- [ ] Manual chapter documenting the composition rules

### 4.2 Polarization elements (per-ray Jones state + Fourier vector diffraction already exist)

See §2.1 above. Elements + metrics belong to Thrust B's coronagraph work but live architecturally here.

### 4.3 Vector diffraction

Already exists in macos. No new work needed beyond documentation in Thrust A.

### 4.4 Multi-field-point optimization metrics

See §2.1 and §3 above.

### 4.5 Metrology completion

Data model in place today: per-element arrays `iEltToMetSrf`, `nMetPos`, `tMetSrf`, `ntMetPos`, `SrfMetPos(3,mMetPos,mMetSrf)`, `SrfMetMea(mMetBeam,mMetSrf)`, `metBeamFlg(mMetPos,mMetPos,mMetSrf)`, flat output `metMeasBuf(mSysMetBeam)` + `nMetMeas`; edge-sensor pair `EdgeSensor(mElt)` + `EdgeSensVec(9,mES,mElt)`; limits `mMetSrf=20`, `mMetBeam=48`, `mMetPos=48`, `mSysMetBeam=128`, `mES=10`. `METcalc` command (3-char min-match) dispatches `SrfMetCalc` ([utilsub.F:1783](macos_f90/utilsub.F#L1783)) — pure 3-D Euclidean distance between point pairs in the global frame. Rx keywords `nMetPos` / `tMetElt` / `EdgeSensors` / `ShowMetData` ([msmacosio.inc:1894](macos_f90/msmacosio.inc#L1894), [3257](macos_f90/msmacosio.inc#L3257)). GMI mex exposes output #9 `MetMeas` via `pflg(27)` triggering `SMACOS('METcalc',...)` ([GMI.F:1448](../MACOS_resources/GMI/GMI.F#L1448)). Example Rx exercising the full path: [optiixonaxisz1_v4_pmsm_met.in](../MACOS_resources/GMI/optiixonaxisz1_v4_pmsm_met.in).

- [ ] PERTURB coverage gap — `SrfMetPos` updated by only 2 of 5 perturbation paths
  - [x] `CPERTURB` ([funcsub.F:184](macos_f90/funcsub.F#L184))
  - [x] `CPERTURB_GRP` ([funcsub.F:1254](macos_f90/funcsub.F#L1254))
  - [ ] `CPRead` (funcsub.F) — add rotate-then-translate SrfMetPos block
  - [ ] `CPERTURB_2` (macos_ops.F) — same pattern
  - [ ] `LnkEltCPERTURB` (lnk_pert.inc) — same; linked-element metrology silently wrong without this
  - [ ] Regression: per-path test pinging SrfMetPos before/after perturbation, verify they all match
- [ ] `EdgeSensor` / `EdgeSensVec` parsed-but-unused dead path
  - **Customer-of-record: the design-layer metrology tier**
    (PLAN_DESIGN_LAYER §6.6).  Its backend #1 is the SegMirMaker
    `Hx.m` edge-sensor model — MATLAB-native, trusted, and
    independent of this engine path.  **Recommendation: remove (or
    defer) the engine `EdgeSensor` keyword** rather than build a
    consumer; the metrology Jacobian-contract interface (§6.6 tier 1)
    has no need for it, and an unconsumed-but-parsed keyword is
    exactly the silent-acceptance trap §6.6 tier-2 (Q8) is meant to
    catch.  Revisit only if a non-design-layer consumer appears.
  - [ ] Decide (per the above): remove the keyword to stop silently accepting data — recommended — OR implement consumer (edge-sensor-distance measurement command, Hx-matrix output)
  - [ ] If keep: define semantics matching SegMirMaker convention (position + two direction vectors → differential displacement)
  - [ ] Add output formatter if kept
- [ ] `SPFcalc` dangling reference
  - [ ] [GMI.F:1454](../MACOS_resources/GMI/GMI.F#L1454) calls `SMACOS('SPFcalc')` for Steward Platform measurements; no matching case in `macos_cmd_loop.inc`
  - [ ] Decide: implement (companion to METcalc — leg-length measurements modeling actual Stewart-platform geometry) OR remove from GMI mex side
  - [ ] If implement: mirror metrology arrays as `nPfPos` / `tPfElt` / `pfMeasBuf` / `nPfMeas`; flag wired to `pflg(28)`
- [ ] Prescription output formatter missing
  - [ ] [iosub.inc](macos_f90/iosub.inc) `PrtSingleEltInfo` emits no `nMetPos` / `tMetElt` / `EdgeSensors` fields
  - [ ] Load → modify → save round-trip silently drops metrology config (and edge-sensor config once consumer lands)
  - [ ] Add output blocks following per-element keyword format already used elsewhere
  - [ ] Regression: round-trip `optiixonaxisz1_v4_pmsm_met.in` and confirm METcalc output identical
- [ ] Line-of-sight obstruction handling
  - [ ] Today `SrfMetCalc` reports geometric Euclidean distance regardless of intervening obscurations
  - [ ] Optional: opt-in `ShowMetObscur=Y` flag that ray-traces each beam through prescription obscuration geometry, reports blocked beams as such instead of silently reporting straight-line length
- [ ] Pymacos exposure
  - [ ] Implement `met_calc()` / `met_measurements()` Python wrappers following the `intensity()` / `complex_field()` `_cmd` + `_get` template
  - [ ] Resolve the `# [ ] perturbElt_METROLOGY_NODES` TODO at [macos.py:2795](../MACOS_resources/pymacos/src/pymacos/macos.py#L2795) — `SrfMetPos` updates flow through automatically once the PERTURB gaps above are closed
  - [ ] Pymacos regression: 6-pt metrology config on a Cassegrain-like Rx, perturb SM, confirm beam lengths track analytically
- [ ] GMI regression coverage
  - [ ] `optiixonaxisz1_v4_pmsm_met.in` already exercises METcalc but isn't in the GMI regression suite — add it (alongside §1.1 GMI test expansion)
- [ ] Sanity-check capacity limits for HWO-class systems
  - [ ] Re-evaluate `mMetSrf=20`, `mMetBeam=48`, `mSysMetBeam=128` against realistic segmented-mirror metrology layouts
  - [ ] Document the limits in the manual either way
- [ ] Manual chapter (see §6.2)

---

## 5. Wrapper layer (Python + MATLAB sharing one Fortran API)

The factoring that keeps maintenance cost bounded if both wrappers exist.

- [x] Factor `pymacos.f90` into language-neutral `macos_api_mod.F` + thin `pymacos_f2py.f90` bridge
  - [x] Codegen script `gen_pymacos_refactor.py` parses pymacos.f90, emits the two split files mechanically (kept as a one-shot tool; not in-tree)
  - [x] `macos_api_mod.F` (4280 lines): pure Fortran 2008 module `macos_api_mod`, no `!f2py` annotations, holds all module state and routine bodies. This is what the future mmacos MATLAB mex calls into.
  - [x] `pymacos_f2py.f90` (1295 lines): defines `module api` (the Python-facing name `lib.api`), `use macos_api_mod, only: NAME_impl => NAME` for every routine, then thin forwarding wrappers carrying the `!f2py` annotations and calling `NAME_impl(...)`. Each wrapper is just a signature copy + one `CALL` line.
  - [x] `CMakeLists.txt` wires both into the f2py build — `macos_api_mod.F` compiled first (its module is consumed by the wrapper), only `pymacos_f2py.f90` parsed by `numpy.f2py` to derive the Python signatures.
  - [x] pymacos regression: 6601/6601 pass in 131 s (vs 207 s baseline; same numerics, faster from warm cache).
  - [x] PROPER-compare suite green — physical-optics agreement preserved at the same numeric precisions as before.
- [x] Build `mmacos_mex.F` MATLAB bridge over `macos_api_mod.F90`
  - [x] Single `mmacos.mexa64` (~2.2 MB) with command-string dispatch, lives at [`MACOS_resources/mmacos/`](../MACOS_resources/mmacos/)
  - [x] `mmacos_mex.F` ~530 lines: parses PRHS(1) as the command name, `SELECT CASE` to per-command helpers; each helper validates args, copies in via `mxCopyPtrToReal8`, calls the corresponding `macos_api_mod` routine, copies out via `mxCopyReal8ToPtr`
  - [x] MVP commands wired: `init`, `load_rx`, `save_rx`, `modified_rx`, `n_elt`, `opd`, `intensity`, `complex_field`, `dx_at`, `base_unit_to_metres`, `apodize`, `perturb_elt`
  - [x] `Makefile` follows GMI's pattern — MATLAB auto-detected (lex-latest under `/usr/local/MATLAB`), ifx + gfortran arms (ifx default for now; macos_api_mod uses ifx-only idioms today), `-reentrancy=none` on the ifx link to avoid libifcoremt parked-thread SIGSEGV at MATLAB process exit (same workaround as GMI)
  - [x] Built clean under ifx against `build_release/` libsmacos.a
  - [x] [`test_mmacos.m`](../MACOS_resources/mmacos/test_mmacos.m) smoke test green via `matlab -batch`: init + load_rx (Rx_Cass_FarField.in, 6 elements) + n_elt + base_unit_to_metres + modified_rx + perturb_elt round-trip → **6/6 pass, exit 0, clean MATLAB teardown**
  - [x] `trace_rays` command added; trace-dependent commands `opd`, `intensity`, `complex_field`, `dx_at`, `apodize` now all exercise green in the smoke test
  - [ ] Additional commands as needs surface: `sxp`, `xp_get/set/fnd`, `spot_cmd/get`, `ray_info_get/set`, `stop_info_*`, `ors_run`, `srs_run`, source-side `set_src_*` / `get_src_*`, all the `elt_*` Rx-inspection / -modification routines
  - [x] **macos_api_mod promoted into libsmacos.a** — file moved from `MACOS_resources/pymacos/src/cmake/source/macos_api_mod.F90` to `macos/macos_f90/macos_api_mod.F90`; added to `SMACOS_ONLY_SOURCES` in top-level `CMakeLists.txt`. Both pymacos and mmacos now just link `libsmacos.a` and pull the `.mod` from `mod_smacos/`. Eliminates the cross-repo `API_MOD_SRC` Makefile path entirely. Inlined the only constant `mPix2` used from the dropped `pymacos.inc`. Verified clean rebuild of macos + pymacos + mmacos, pymacos pytest 6601/6601, PROPER-compare all phases green, mmacos smoke 11/11.
  - [x] **macos_api_mod ported to gfortran-clean idioms** — fixes applied: replaced `(/'m', 'cm', 'mm', ...'/)` with explicit `character(len=4)` arrays (gfortran rejects the implicit-length form for line-truncation reasons); converted 16 `do concurrent` constructs to plain `DO` loops (some inside `BLOCK` for index-decl scoping) per the project no-DO-CONCURRENT convention; bulk-replaced 50 `LOGICAL == PASS/FAIL` sites with `.eqv.`; per-site-patched 16 `INTEGER == PASS/FAIL` sites to integer equality (`var == 0` / `var /= 0`) where the variable was declared INTEGER not LOGICAL; rewrote `prb_elt_grp`'s `ifGlobal` LOGICAL-array tests (`ifGlobal==0`, `ifGlobal==1`, `ifGlobal(i)>0`) to use `.not.` / `.and.` / direct LOGICAL tests; rewrote `set_src_csys`'s `filter/=0` (LOGICAL) to bare `filter`. Verified clean compile + numerical parity under both ifx and gfortran.
  - [x] **mmacos `FC ?= gfortran` default flipped** (matches GMI); ifx still works via `make FC=ifx` and remains the GMI-style fallback. gfortran smoke test green at 11/11 with clean MATLAB teardown.
- [ ] Side-by-side documentation
  - [ ] "Pymacos call → mmacos call → underlying Fortran" reference card
  - [ ] Helps users move between languages
- [ ] Long-term decision: does mmacos replace `call_GMI` or coexist as a finer-grained alternative?
  - [ ] Survey existing GMI users
  - [ ] Plan a backward-compat shim if mmacos absorbs GMI workflows
  - [ ] Document migration path

### 5.4 mmacos → pymacos feature parity (CodeV/PROPER regression + dw/dz, dw/dx)

Goal: bring mmacos to the same user surface as pymacos so the CodeV and
PROPER regression suites and the dw/dz_Zernike and dw/dx sensitivity
drivers can be run from MATLAB with bit-identical numerics. Two layers
exposed throughout — power users can call `mmacos('cmd', ...)` directly,
casual users get the `MacosSession` class veneer. Phases ordered so each
unlocks the next.

**Standing rule (added 2026-05-30):** treat regression-suite growth as a
side activity in every subsequent phase. When a new +macos wrapper,
helper, or mex command lands, add a `tCodeV*` (or equivalent) test
that exercises it — even if the immediate motivating task didn't
require it. Goal is a continuously expanding `make unittest` covering
the realistic mmacos surface, so by the time Phase 8 (cross-language
verification) runs there's already substantial mmacos-side coverage
to compare against pymacos.

- [x] **Phase 1 — Command-surface parity (codegen)**
  - [x] Codegen script `MACOS_resources/mmacos/gen_mex_wrappers.py`
    parses `macos_f90/macos_api_mod.F90`, emits `do_<name>` mex
    helpers and a `gen_dispatch` fallback into `mmacos_gen.F`.
    Picks up arg types (logical/integer/real/character), `intent`,
    array dim symbols (including dim args read first so subsequent
    allocations have sizes), local `integer, parameter` aliases for
    elt_mod dim symbols, multi-line subroutine arg lists.
  - [x] Generated dispatcher: 78 codegen routines + 13 hand-written
    cases in `mmacos_mex.F`. Hand-written cases (init, load_rx,
    save_rx, modified_rx, n_elt, opd, intensity, complex_field,
    dx_at, base_unit_to_metres, apodize, perturb_elt, trace_rays)
    keep their bespoke logic; everything else falls through to
    `gen_dispatch`.
  - [x] Smoke test extended: 15/15 pass under gfortran (R2026a),
    incl. four codegen-routed commands (`get_src_sampling`,
    `elt_grp_max_all`, `src_wvl` getter, `elt_vpt` getter).
  - [x] Low-level surface: `mmacos('cmd', args...)` covering 91
    commands total (was 13). Skipped: `elt_csys_get` (rank-3
    array — hand-write if exercised), single-element `perturb_elt`
    api (collides with hand-written array form; expose as a
    distinct cmd later if needed).

- [x] **Phase 2 — `+macos/` package + `macos.Session` class** (high-level layer)
  - [x] `+macos/` MATLAB package — one `.m` per public function.
    Improvements over pymacos's surface (with user buy-in to "make it
    better"): split getters/setters into `get_*` / `set_*` (MATLAB-
    idiomatic, autocomplete-friendly), single-element vs array
    perturbations split into `perturb` / `perturb_many` (was confusing
    `perturb` vs `prb_elt` in pymacos), `trace` returns a struct
    (`s.nRays`, `s.rmsWFE`), `dx_at(srf, unit)` does unit conversion
    in MATLAB (`'m'` default, `'mm'`, `'cm'`, `'um'`, `'native'`).
    Coverage v1: 23 functions — init / load_rx / save_rx / modify /
    num_elt / has_rx / cbm / sys_units / get_src_sampling /
    set_src_sampling / get_src_wvl / set_src_wvl / get_elt_vpt /
    set_elt_vpt / get_elt_psi / set_elt_psi / get_elt_rpt /
    set_elt_rpt / perturb / perturb_many / perturb_src / trace /
    opd / intensity / complex_field / dx_at / apodize.
  - [x] `macos.Session` classdef in the same package — thin handle
    over the function layer.  Constructor takes `model_size`,
    every method delegates to the same-named package function.  No
    per-instance Fortran state (libsmacos.a owns it).
  - [x] Mex layer cleanup: hand-written cmd `'perturb_elt'` renamed
    to `'prb_elt'` (matches the api routine it actually calls —
    the array form).  Codegen now emits cmd `'perturb_elt'` mapping
    to the single-element api routine; surface is more discoverable.
  - [x] Smoke test `test_macos_pkg.m`: 25/25 pass.  Three sections:
    (A) function-style end-to-end, (B) Session-style end-to-end,
    (C) cross-surface state coherence (function-set → class-get and
    class-set → function-get both observe the same value).
  - [x] **Both layers remain first-class.** Power users call
    `mmacos('cmd', ...)` directly; routine work calls `macos.*`
    functions; OO-flavor code uses `m = macos.Session()`.

- [x] **Phase 3 — Test infrastructure**
  - [x] `matlab.unittest` skeleton in
    `MACOS_resources/mmacos/tests/` with shared utility
    `tests/private/rx_fixture_path.m` (pulls from the pymacos Rx
    corpus so no duplication).  Five test classes:
    `tMmacosCmd` (raw mex layer), `tMacosPkg` (+macos functions),
    `tMacosSession` (class delegation), `tCrossSurface` (mutate via
    one surface, observe via another — proves shared backend),
    `tPerturbRoundtrip` (regression-pins the ULP residual finding
    so a future psi-renormalize fix in CPERTURB_PROG doesn't break
    cleanly — see §0 follow-up).
    50 tests, all pass.  ~6 s cold + 5 s suite.
  - [x] `run_mmacos_tests.sh` analog of `run_proper_tests.sh` — auto-
    rebuilds mex if stale, supports `-k <substring>` and direct class-
    name filter args.
  - [x] Wired into `make unittest` alongside the existing `make test`
    quick-smoke target.  Quick-smoke (`test_mmacos.m`,
    `test_macos_pkg.m`) kept as readable `fprintf`-style diagnostics —
    different intent from unittest (CI / assertion layer).
  - [x] Tolerance / `.mat` reference loader deferred — not needed yet.
    Phase 4 (CodeV regression port) will add them when the suite
    starts comparing numerically against pymacos reference outputs.

- [ ] **Phase 4 — CodeV regression port**
  - [x] **Slice 1**: `test_api_rx_grating.py` (17 tests) → `tests/tCodeVGrating.m`.
    Added 12 +macos grating wrappers (`elt_grating_any`,
    `elt_grating_fnd`, `get/set_elt_grating_params`,
    `get/set_elt_grating_order`, `get/set_elt_grating_rulewidth`,
    `get/set_elt_grating_type`, `get/set_elt_grating_dir`).
    Reference values transcribed from
    `pymacos/tests/rx_data.py::Rx_Grating_001()` into
    `tests/private/rx_grating_001_data.m` (single fixture; tiny).
    17/17 pass; full suite now 67/67.  Surfaced and fixed one
    arg-list-order bug in `set_elt_grating_params` — the api signature
    is `(ok, iElt, Spacing, Diff_Order, h1HOE_, reflective, setter)`
    but the declarations list them in `(ok, iElt, Diff_Order, Spacing, ...)`
    order, which had tricked a swap.  Audit of other +macos wrappers
    found no other ordering mismatches.
  - [x] **Slice 2**: `test_masks.py` (8 classes, ~10K sub-cases) → 8
    `tCodeV{Ape,Obs}Masks{Circ,Ellipse,Rect,Polygon}.m` classes.
    Surprise on review: this file does NOT consume `rx_data.py` fixtures
    (only the grating slice does).  Tests instead **mutate the .in text
    file directly** (overwriting specific lines with new `ApType=` /
    `ApVec=` / `nObs=` values), reload via `macos.load_rx`, trace, and
    check that ray-fall geometry sits inside the analytic mask shape.
    So no `.mat` export pipeline needed; the .m transcription path
    continues to work.
    Pymacos's matrix counts 6584 sub-tests because pytest expands
    each parametrize axis individually; mmacos collapses each
    (Rx, srf, sub-test variant) into a single matlab.unittest method
    that LOOPS over the geometric parameters, with per-iteration
    diagnostics on failure.  62 methods total, ~10K assertions under
    the hood, all green; full suite 129/129.
    Landed:
    - Hand-written `do_elt_csys_get` in `mmacos_mex.F` (rank-3
      `csys(6,6,N)` is outside the codegen ≤2D ceiling; uses
      `mxCreateNumericArray` for the 3D output).  Added to
      `HAND_WRITTEN_CMDS` so the codegen inventory lists it.
    - Two new +macos wrappers: `get_ray_info(N)` and
      `get_elt_csys(srf)`, both returning structs.
    - 5 geometry helpers in `tests/private/`: `rectangular_polygon`,
      `hexagon`, `poly_lines`, `chk_polygon_pts`,
      `ray_pos_at_srf_in_tangent_plane`.  The last writes a mutated
      Rx text via MATLAB's `writelines`, loads it, runs two traces
      (to FocalPlane then to srf), and projects rays into srf's
      local 2D tangent plane via `csys_local^T · (ray_pos − vpt)`.
    - `tests/private/rx_mask_params.m` transcribes
      `test_masks.py:42`'s RX_PARAMS — but with **1-based** line
      indices (pymacos stored 0-based slice indices).  Off-by-one
      caught during port.
    - `tests/private/tolerances.m` mirrors pymacos `_Tol` with named
      fields (`.P`, `.r`, `.L`, `.v`, `.eps`).  Available for future
      tests even though the mask suite uses geometric assertions.
    Bug surfaced + fixed during port: the polygon-aperture global-
    vertex variant (`run_hex_glb`) was using `pts - [dx, dy]` in the
    inside-check, but pymacos uses `pts - [dx*dx_fact, dy]` because
    the PolyApVec write applies `dx*dx_fact` to the shift.  For
    `dx_fact = +1` the two expressions agree; for `dx_fact = -1`
    (parabola_glb srf3) they don't.  Failed 3 cases, fixed to match.
  - [x] **Bit-identical pass — partial:**
    - Slice 1 (grating) — **met directly**.  Reference values are
      transcribed verbatim from `rx_data.Rx_Grating_001()`; tests
      assert `verifyEqual(value, ref_value, AbsTol 1e-15)`.  Round-
      trip setters return the set value at machine precision.
    - Slice 2 (masks) — **met indirectly**.  Tests are geometric
      (every ray lies inside/outside the analytic mask shape), same
      assertion form pymacos uses.  Both languages drive the same
      libsmacos.a, so the trace outputs ARE bit-identical, but the
      tests don't directly compare the ray arrays — they assert the
      geometric consequence.  Direct cross-language ray-position
      comparison belongs in Phase 8 (export `ray_info.pos` from
      pymacos to `.mat` for a canonical (Rx, srf) combo, have mmacos
      load it and `verifyEqual(pos_mm, pos_py)`).

- [ ] **Phase 5 — PROPER regression port**
  - [x] **MATLAB PROPER installed** at `~/dev/proper_matlab/`
    (v3.3.1 from sourceforge.net/projects/proper-library, MATLAB
    translation by Gary Gutt of Krist's IDL original).  Smoke-tested:
    `prop_begin` returns a wavefront struct.
  - [x] **Slice 1 — Phase 1 (Cass-FF)** ported.
    `tests/proper_compare/tProperCompareCassFF.m` — 4 tests:
    `test_proper_cass_ff_runs`, `test_macos_cass_ff_runs`,
    `test_compare_cass_ff_psf` (analytical PROPER aperture; max|a-b|
    aligned < 0.1), `test_compare_cass_ff_psf_with_opd` (macos OPD +
    amplitude mask passed through to PROPER; max|a-b| ≈ 1.1e-11,
    matching pymacos's reported precision).  All 4 pass.
    Infrastructure:
    - `tests/proper_compare/+geometries/cass_farfield.m` — geometry
      params struct.
    - `tests/proper_compare/proper_run_cass_ff.m` — PROPER driver
      (supports macos_opd / macos_amplitude pass-through, opd_sign_flip
      to reconcile macos's OPD convention with PROPER's prop_add_phase).
    - `tests/proper_compare/macos_run_cass_ff.m` — mmacos driver.
    - `tests/proper_compare/private/compare_and_record.m` — minimal
      comparison harness with metrics struct + optional 3-panel PNG
      (macos | PROPER | difference) written to `results/phase<N>/`.
    - `tests/proper_compare/private/{centroid_loc, crop_center,
      embed_macos_array_in_proper_grid}.m` — supporting helpers.
    - `run_mmacos_tests.sh` adds PROPER + proper_compare/ to MATLAB
      path, walks subfolders so `TestSuite.fromFolder` picks up the
      new class.
    - `.gitignore` excludes `tests/proper_compare/results/phase*/`
      so per-run PNGs don't churn the repo (analog of pymacos's
      gitignored results dir).
    - **`run_mmacos_tests.sh` now splits the full-suite run into
      per-model_size matlab -batch invocations** (§0 follow-up logs
      the underlying macos engine bug).  Without the split, Phase 5
      tests (model_size=512) running after Phase 3/4 tests
      (model_size=128) in the same MATLAB session trigger
      `malloc()` heap corruption on the next trace.  Same issue
      pymacos works around via per-phase pytest processes.
    - Named groups for dev-loop speed: `./run_mmacos_tests.sh fast`
      (size=128 non-masks, ~10 s), `masks` (size=128 masks, ~10 min),
      `proper` (size=512 PROPER, ~15 s).  Full `./run_mmacos_tests.sh`
      runs everything (~11 min).  Iterating on a Phase 6+ slice that
      doesn't touch masks should use `fast` between edits and reserve
      the full suite for pre-commit.
  - [x] **Slice 2 — Phase 2 (Coro NF)** ported.
    `tProperCompareCoroNFprop.m` — 3 tests covering NFPlane Elt 2→3
    of Rx_Coro.in (774 mm Fresnel).  Both engines fed the same
    `complex_field` at Elt 2; sum-norm comparison.  Pymacos hits
    2.5e-10 max; mmacos hits **4.836e-13 max** — 500× tighter (see
    FFT-backend curiosity below).  model_size=1024, new third group
    in run_mmacos_tests.sh.
  - [x] **Slice 3 — Phase 3** ported.
    `tProperCompareCoroPhase3.m` — 6 tests covering the rest of the
    Coro chain: 3a NFPlane (5→6), 3b sphere→plane (8→9), 4a pupil
    →pupil through focus (8→10), 4b NFPlane (13→14, no Lyot),
    5.1 ExitPupil→SciFP no mask, 5.2 ExitPupil→SciFP with FPM+Lyot.
    All numerics match pymacos's reported bounds.  Added geometries
    `coro_sphere_to_plane` + `coro_pupil_to_pupil` and the matching
    macos_run / proper_run drivers.
  - [x] **Slice 4 — Phase 6a (apodiser)** ported.
    `tProperCompareCoroApodizer.m` — 1 test, soft Gaussian-edge
    apodiser at Elt 5 of Rx_Coro_noLyot.in, then NFPlane to Elt 6.
    Same mask handed to both engines (macos via `macos.apodize`,
    PROPER via `prop_multiply`).  max=3.971e-08 — matches pymacos's
    4e-8 exactly.  Added `+apodizer/` package
    (`build_apodised_mask`, `circle`, `gaussian_edge_taper`).
  - [x] **Slice 5 — Cass-FF aberrations** ported.
    `tProperCompareCassFFAberrations.m` — 6 parametrized cases
    (nominal + Tx/Ty/Tz ±1-5 µm perturbations to M2).  Each
    perturbs macos, captures OPD at the exit pupil, hands it to
    PROPER via prop_add_phase, compares focal-plane PSFs.  max=
    9.337e-12 — same precision as the unperturbed Phase 1 result.
  - [x] **Slice 6 — band-limited masks + DM phase imprint** ported.
    Two new test classes:
    - `tBandLimitedMask.m` — 3 tests (pure math, no macos/PROPER):
      BL disc integral equals area / dx² to 6+ digits; super-sampled
      disc integral to ~1e-4; BL central peak + radial mean at r=r0
      invariant across N ∈ {128,256,512,1024} to better than 1-2%.
      New +apodizer helpers: `band_limited_circle`,
      `build_band_limited_mask` (analytic FT of unit disc, ifft2 to
      real space, fftshift).
    - `tProperCompareCoroDMPhase.m` — 3 parametrized cases (Z4
      defocus, 5-cycle sinusoid, filtered noise k=12) at λ/20 RMS
      over a 16 mm support disk.  macos via `apodize_complex` (pure
      phase imprint), PROPER via `prop_add_phase`.  All three hit
      sum-norm RMS ~1.2e-9, max ~3.7e-8 — matches pymacos exactly.
      New +macos wrapper: `apodize_complex` (splits MATLAB complex
      mask into real+imag for the cfield_apodize_complex api).
      Codegen now emits `do_cfield_apodize_complex` (removed from
      the HAND_WRITTEN skip list).
  - [ ] **Slices 7+** still pending:
    - `test_coro_dm_grid_self.py` (6 parametrized = 3 shapes × 2
      detectors; macos self-consistency between apodize_complex and
      elt_grid).  Needs a +macos.elt_grid wrapper — non-trivial
      because pymacos transposes the grid array for the Fortran
      column-major API and MATLAB is already column-major (easy to
      get wrong).  Test has no hard assertion, only a diagnostic
      print — value is mostly in the wrapper + PNG artefacts.
    - `test_coro_aberrations.py` (5 parametrized chain runs).
      Needs ~500 LOC of port: `SystemState`/`CORO_LAYOUT`/
      `run_chain_with_state` from pymacos's `aberrations.py`, plus
      `lambda_over_D_pixels`/`plot_contrast_curves` from
      `contrast.py`.  Mirrors the dw/dx state-vector model that
      Phase 7 will build separately — defer until Phase 7 lands the
      driver code, then both share the chain runner.
    - `test_psf.py` (1 real, 1 skipped) — duplicates Cass-FF sanity,
      low priority.
    - `run_broadband_*.py` are not tests; report generators.
  - [x] Validate at the same RMS / max numeric tolerances pymacos hits
    — met for Slice 1 (Cass-FF nominal_with_opd: 1.1e-11 vs pymacos's
    reported ~1e-11) and Slice 2 (Coro NF: 4.836e-13 sum-norm max —
    500× TIGHTER than pymacos's 2.5e-10, see curiosity below).
  - [ ] **Curiosity (logged 2026-05-30, defer):** mmacos's Phase-2
    NF comparison hits 4.836e-13 sum-norm max vs pymacos's reported
    2.5e-10 — same macos backend, same geometry.  Both PROPER
    implementations are the same algorithm but different FFT
    backends: MATLAB PROPER (Gutt's port) calls MATLAB's `fft2`,
    which on R2026a dispatches to Intel MKL; PyPROPER3 (Krist's
    Python port) uses NumPy's pocketfft by default.  Likely
    explanation: MKL's accumulation strategy is more precise than
    pocketfft on 1024² complex transforms.  Confirm by saving
    `complex_field` at Elt 2 to .mat and feeding it through both
    PROPER backends directly (strips the macos hand-off).  Practical
    consequence: mmacos-side PROPER tolerances can be tighter than
    pymacos's — better regression detection at the cost of
    cross-language reproducibility.
    **Remediation path (closes the 500× gap):** install `mkl_fft`
    (`pip install mkl-fft` or `conda install -c conda-forge mkl_fft`)
    and monkey-patch numpy.fft before importing PyPROPER3:
    ```
    import mkl_fft.interfaces.numpy_fft as _mkl
    import numpy.fft as _np_fft
    _np_fft.fft2  = _mkl.fft2
    _np_fft.ifft2 = _mkl.ifft2
    import proper  # picks up the MKL versions transparently
    ```
    PyPROPER3 imports `numpy.fft.fft2` / `ifft2` at module-load, so
    the patch has to happen BEFORE the `import proper`.  Suggested
    integration: an env flag (e.g. `PYMACOS_USE_MKL_FFT=1`) that
    triggers the patch from pymacos's `tests/proper_compare/conftest.py`
    before any test imports PROPER.  Default off so existing pymacos
    numerical references don't shift; opt-in for tighter agreement.
    Caveats: pin `MKL_NUM_THREADS=1` for the deterministic path
    (multi-thread MKL reductions vary at the ~1-ULP level); literally
    bit-identical to MATLAB requires matching MKL version + alignment
    + thread count, but "comparable" (within a few ULP) is realistic.
    Other options ranked: scipy.fft + MKL backend (same as mkl_fft,
    newer API) > pyFFTW (FFTW3 — better than pocketfft, not quite
    MKL's accumulation order).  Skipping for now: mmacos suite is
    already at the FFT-precision limit; pymacos would need this for
    its own tighter regression bars.

- [ ] **Phase 6 — `dw_dz_zernike` driver**
  - **Build pymacos first, then port to mmacos.** pymacos already has
    the single-field worker (`tests/sensitivities/dw_dz_zernike.py`,
    346 LOC); we're only adding the multi-field supervisor + the
    sensitivity-bundled `.mat` format.  mmacos backport is then a
    mechanical translation of a validated design.
  - **Risk to mitigate:** pymacos-first can leak Python idioms into
    the mmacos port (dicts, list comprehensions, numpy fancy
    indexing, Python-specific .mat field shapes).  Mitigation:
    state the supervisor API contract in language-neutral form
    BEFORE either implementation — name the field-set datatype
    (Nx2 array, not list-of-tuples), the stack-output names
    (`dwdxall`, `w0_stacked`, `indxall`, `field_table`), the
    .mat layout, the worker→supervisor calling convention — then
    implement the same shape in both.  The m2v idiom is already
    shared; the worker/supervisor surface stays small.
  - **Canonical output: `dwdxall`.**  The supervisor's bundled
    output is the full state-vector control model
    `wall = dwdxall * x + w0_stacked`, where `wall` stacks per-field
    wavefronts, `x` is the channel coefficient vector, `w0_stacked`
    is the per-field nominal-OPD baseline.  `dwdxall` is the
    LOAD-BEARING output (not just "a convenience stack" — it's what
    feeds future control-loop calcs).  The .mat must name it
    `dwdxall` for downstream consumer scripts.
  - **Dev-cycle speed: reduce rays.**  During development run
    `macos.set_src_sampling(64)` (or smaller) so each per-field
    trace finishes in ~0.5 s; full prescription sampling only for
    the validation run.  Same trick saves ~10x wall time on the
    multi-field iteration loops.
  - **Multi-field design (worker + supervisor split):**
    The same prescription is loaded once; the source is tilted to
    each field point in turn via `set_src_fov` (or `perturb_src` of
    Elt 0 by the y-field-angle).  Default field set is 5 points —
    center + 4 corners (upper-left, upper-right, lower-left,
    lower-right of the chief-ray FoV).  Variants: N×N grid (every
    1/n step in field x and y) or 5-points-with-FSM-angle-tweak to
    expose a wider field of regard.  One Rx for all field points.
    Per-field workflow:
    ```
    load -> STOP+FEX -> single-field dw_d* -> perturb_src by next field -> repeat
    ```
    The dwdz/dwdx matrices stack vertically (one block per field):
    `w_stacked = dwdx_stacked * x + w0_stacked`.  Use the same m2v
    bookkeeping as the single-field case — build one big nominal-OPD
    matrix `OPDall` with per-field tiles separated by `zeros(nsize)`
    spacers, then `[wall, indxall] = m2v(OPDall)` gives the linear
    index map that handles both the math (stacking) AND the display
    (`plot(OPDall)` shows the spatial layout — center middle, corners
    at corners — for free, no GridSpec / subplot positioning needed).
  - **Phase 6 mmacos sub-tasks** (after pymacos lands):
  - [ ] `+macos/dw_dz_zernike.m` (single-field worker) plus the
    `+channels/` submodule (ZernikeCoefChannel only).
  - [ ] `+macos/dw_dz_zernike_multi.m` (supervisor; takes field set,
    loops worker + perturb_src, stacks dwdz + builds OPDall).
  - [ ] Driver takes a `MacosSession` instance (mirrors pymacos's
    `MacosModel` arg); a free-function form taking no session also
    exposed for raw-mex users.
  - [ ] Writes the same `.mat` format pymacos's `dw_dz_zernike` does
    so downstream consumers don't fork.

- [x] **Phase 7.a — `dw_dx` single + multi-field, per-element + source + FP**
  - Per-element `RigidBodyChannel`, `SourceChannel`,
    `FocalPlaneChannel` (fp_mode ∈ {track, srs, sxp, none}; fex
    deferred — no Fortran wrapper).
  - `+macos/dw_dx.m` single-field driver + `+macos/dw_dx_multi.m`
    supervisor (same 5-field default / NxM grid / fields-file
    pattern as Phase 6).
  - Numerical agreement on e5hex1: BITWISE mmacos == pymacos for the
    full 6 DOFs (60 non-FP channels were bitwise from the start; FP
    rotation channels carried a 7e-7 m spatial residual root-caused
    to MATLAB `norm()` vs NumPy `np.linalg.norm` differing by 1 ULP
    on the EP-psi normalize in `FocalPlaneChannel.rotate_ep_about_fp_rpt`
    — fixed by switching to `sqrt(sum(v.^2))`; commit 345dbe8).
  - tDwDx 5/5 pass; `examples/sensitivities/e5hex1/run_dwdx.m` end-to-end.
- [x] **Phase 7.b — Group channels (GPERTURB) + dw/dx integration**
  - `GroupedRigidBodyChannel` (handle classdef), `prb_grp` /
    `get_elt_grp` / `set_elt_grp` / `del_elt_grp` foundation,
    `grouped_rigid_body_channels` builder, `parse_rx_groups`
    EltGrp= parser, dw/dx driver `'groups'` + `'groups_auto'` opts.
  - Numerical agreement on e5hex1 'Cam' group ([9,10,12,13]):
    BITWISE mmacos == pymacos standalone (6 DOFs) and integrated
    (33 per-elt Tx,Ty,Tz + 3 group Tx,Ty,Tz = 36 channels,
    max|diff| = 0).
  - tDwDxGroups 4/4 pass (commit 1a80e0e).
  - Quirk preserved (mirrors pymacos): GroupedRigidBodyChannel
    translation increment is passed straight to `prb_grp` in the
    Rx's BaseUnits, NOT SI metres — per-element `RigidBodyChannel`
    does the CBM conversion via `macos.perturb`, group does not.
    Locked in for bit-identical output; if pymacos ever changes
    this convention, mmacos follows in lockstep.
- [ ] **Phase 7.b.2 — Group consistency tools (deferred)**
  - Port `predict_global_rigid_response` (analytical prediction of
    a global rigid-body group response from per-element local-frame
    columns + the rigid-body kinematics) and `group_synthesis_matrix`
    (the column-combination weights that test whether per-element
    sensitivities correctly reconstruct the group column).  Cross-
    check tool, not on the dw/dx hot path — defer until a user
    workflow needs it.

- [ ] **Phase 8 — Cross-language verification**
  - [ ] Run identical Rx + perturb sequences through pymacos and
    mmacos; assert numerical equality at machine precision.
  - [ ] **Anchor cases to add first (called out by Phase 4 slice 2
    deferral):** for a small canonical set of (Rx, srf, mask) combos
    pulled from `tCodeVApe/ObsMasks*`, export the post-trace
    `ray_info.pos` (and `.dir`, `.opl`) from pymacos to `.mat` and
    have an mmacos test class `tCrossLangRayInfo` load + assert
    bitwise equality.  Closes the gap between "both engines pass
    geometric assertions on the same Rx" and "both engines produce
    identical ray-trace arrays."
  - [ ] This is the actual proof of "same backend." Add to CI once a
    CI substrate exists.

Total estimate: 12–16 days, ordered so no back-pressure between phases.
Phases 4, 6, 7 are the user-facing milestones; 1, 2, 3 are scaffolding.

---

## 6. Manual rewrite

### 6.1 Infrastructure

- [ ] Source-of-truth decision
  - [ ] Keep `macosMan4_01.docx` as canonical
  - [ ] Or switch canonical to a less-fragile format (LaTeX, Markdown + pandoc)
  - [ ] Likely stay with Word given user familiarity
- [ ] Formatting cleanup pass on existing chapters
  - [ ] Consistent section / heading / cross-ref / equation styles
  - [ ] Fix the broken cross-references and ToC
  - [ ] Probably requires real Word, not pandoc round-trip
- [ ] Image / figure handling
  - [ ] Decide source-of-truth for figures (SVG? PNG? vector PDF?)
  - [ ] Set up figure regeneration workflow

### 6.2 New chapters and sections

- [ ] Group perturbations (GPERTURB, CPERTURB_GRP)
- [ ] Optimization: `nls_optim_dvr` (always available)
- [ ] Optimization: `np_optim_dvr` (USE_NPSOL only, opt-dev branch)
- [ ] Optimization: target types, DOF taxonomy, convergence controls (see §3)
- [ ] Sensitivity-matrix workflow (Thrust A headline) — engine-level recipes + the design-layer `from_rx`/`vary`/`sensitivities` surface (PLAN_DESIGN_LAYER Sprint 2A-i) as the recommended entry point
- [ ] Polarization ray trace + metrics
- [ ] Vector diffraction
- [ ] Multi-field-point metrics
- [ ] FreeForm surface authoring (`SrfType=14`, FF/Mon/grid composition, `FFZernType` / `MonZernType`)
- [ ] Generalized FreeForm on Toric / Anamorphic bases
- [ ] Metrology + edge sensors (data model, `METcalc`, Rx keywords; see §4.5)
- [ ] Segmented-mirror prescription generator (`SegMirMaker`; see §7)
- [ ] `ParseZernType` modern names + legacy aliases
- [ ] Per-ray status tracking (`RayStatus`, `RayFailMsg`, WARN breakdown)
- [ ] SXP vs FEX (chief-ray-to-FP geometry)
- [ ] IMGMODE / annotation system
- [ ] LOG command (intensity log10 display — name collision with JOURNAL)
- [ ] VALIDATE command
- [ ] Build choice (gfortran default, `-reentrancy=none` for ifx, USE_NPSOL flag, opt-dev branch)
- [ ] Apodizer + extended mask library (after Thrust B)
- [ ] PROPER → macos migration appendix (after Thrust B)
- [ ] FDTD coupling appendix (see §2.4)

### 6.3 Worked examples

- [ ] Cassegrain sensitivity matrix
- [ ] Optiix single-field optimization
- [ ] Optiix multi-field optimization
- [ ] HCIT-style coronagraph (FPM + Lyot)
- [ ] FreeForm segment design + assembly
- [ ] DM phase-imprint validation
- [ ] IRIS observatory model walkthrough
- [ ] FSM beam-steering (from `opt_example`)
- [ ] Asphere coefficient optimization (in progress from `opt_example_asph`)
- [ ] Zernike-coefficient optimization
- [ ] Constrained optimization (USE_NPSOL only)
- [ ] Cached sensitivities (`OptSavePinv` / `OptUseSavedPinv`)
- [ ] Source-side optimization (`VarSrc`, `OptSrcRpt`)

---

## 7. Infrastructure / existing-work absorption

These are independent of the thrusts above; land when bandwidth permits, each lands with a regression test as it touches its surface.

- [ ] Surfsub robustness pass — re-implement engineering intent on current FreeForm code shape
  - [ ] IEEE NaN/INF guards in `ConSrf`, `ConicDeltaS`, `FreeFormSrf` conic-solve sites
  - [ ] Secant-Brent hybrid root finder (Phase 1 secant, Phase 2 Brent fallback)
  - [ ] `SolveConicIntersection` factoring to unify standalone-ConSrf and FreeForm-composite conic solves
  - [ ] Conic-radicand-check documentation (numerical-stability notes)
- [ ] `define_local_csys` robustness / speed improvement
- [ ] Glass dispersion model updates from `dev_optimization_surfsub` branch — pick out landable hunks
- [x] **SegMirMaker** — modernize legacy `SMPGe.for` (VAX-ism `TYPE`/`ACCEPT`/CR-only line endings) and extend it to support a FreeForm parent surface
  - [x] Standalone tool at [`MACOS_resources/segmirmaker/`](../MACOS_resources/segmirmaker/) — `SegMirMaker.f` + `CMakeLists.txt` + `makesegmirmaker.sh`, linking against pre-built `libsmacos.a`
  - [x] Portability pass: VAX `ACCEPT`/`TYPE`/`REPEAT`-`EXIT` → standard prompt+`READ` helpers (`DQUERY`/`IQUERY`/`CQUERY`); LF line endings
  - [x] `SurfCoordFF` wrapper around MACOS `FreeFormSrf` — conic-only path degenerates correctly when parent has `lMon=lFF=nGridMat=0`
  - [x] Parent prescription loader via `SMACOS('OLD',...)`; reads `KcElt` / `KrElt` / `psiElt` / `VptElt` and FF/Mon/grid blocks from chosen `iParent`; sentinel filename falls back to legacy interactive conic dialog
  - [x] Per-segment output writer extended to copy parent FF + grid blocks into every segment block (Mon stays empty per design); conic-parent path preserves legacy SMPGe block format
  - [x] Edge-sensor Hx matrix path now uses composite (conic+FF+grid) tangent plane at midpoint
  - [x] VAX original retired to `Archive/SMPGe.for`
  - [x] Regression corpus in `test_in/` (e5hex2, ff_hex2, ffparent, monparent, psiparent, e5pie, e5mono, ...)

---

## 8. Cross-thrust dependencies + sequencing

### Dependencies

| Thrust | Depends on | Enables / unblocks |
|---|---|---|
| §1 (test + doc) | — | All other thrust work documented + tested |
| §2 (coronagraph) | §5 mmacos MVP for FALCO MATLAB path | Apodizers / masks / pol elements used by §3 regression |
| §3 (optimizer) | — | §2.1 multi-field-point metrics; §1.3 sensitivity-matrix |
| §4 (core caps) | §5 wrappers for exposure to users | §2 element work via shared backend |
| §4.5 (metrology) | PERTURB-coverage gap (3 of 5 paths) | Sensing-chain validation; HWO laser-truss workflows |
| §5 (wrappers) | — | mmacos consumed by §2 FALCO integration |
| §6 (manual) | each preceding thrust's content | User-facing documentation |
| §7 (infrastructure) | — | Underlies §1.2 regression robustness |

### Suggested sequencing

Coronagraph is the priority for elapsed time. Thrusts run in parallel where possible.

**Weeks 1-2**
- §0 remaining items (`define_local_csys`, `dopt_init_vars` audit, parser-noise quieting)
- §5.1 factor `pymacos.f90` → `macos_api_mod.F` + `pymacos_f2py.f90`
- §5.2 mmacos MVP

**Weeks 3-6**
- §2.1 coronagraph element pass (apodizer, polarization, masks, complex grid)
- §3 optimizer robustification (parallel)
- §1.2 function exercise on surface types

**Weeks 4-8**
- §2.2 PROPER-compat module (Python first, MATLAB after §5)
- §2.3 FALCO conversation + inventory
- §1.2 function exercise on element types + commands

**Weeks 6-10**
- §2.1 multi-field-point metrics (rendezvous of §3 and §2.1)
- §1.3 sensitivity-matrix recipe + regression
- §6 manual chapter writing (concurrent)

**Weeks 8-12**
- §4 generalize FreeForm to non-conic bases
- §2.3 FALCO end-to-end demo
- §6 manual cleanup pass

### Test + doc sequencing pattern

For every new capability:

1. **Feature landing (Day 0)**: Fortran implementation + CLI regression test + pymacos wrapper + pymacos regression test in one PR. "Shipped" means: works in CLI, has tests in pymacos.
2. **mmacos wrapper (lag ≤ 2 features)**: each new feature gets a parallel mmacos wrapper batched with the next mmacos update.
3. **Documentation (lag ≤ 3 features)**: manual chapter with CLI prescription + journal + pymacos + mmacos examples side by side.

---

## 9. Branch / repo policy (revised 2026-06)

- `sls-dev` — integration branch for new work; SLSQP default + NPSOL opt-in coexist here
- `opt-dev` — release target; bug fixes only; NPSOL source tree gets removed at promotion to `main`
- `release-candidate` — **frozen** at `19bfbf8` (mZern cherry-pick); no new commits
- `main` — public release surface; promoted from `opt-dev` on each release
- Promotion gate: `sls-dev` fast-forward-merges into `opt-dev` when ready; `opt-dev` strips NPSOL and promotes to `main` for public release
- Test files duplicated where needed: not all macos users access `MACOS_resources`, but all `MACOS_resources` users see `macos`; manual examples and reader exercises must be in `macos` directly; some coronagraph Rxes remain private and stay in `MACOS_resources` only

---

## 10. Strategic decisions already made

- **C-rewrite**: no. Modernize in-place via Fortran 2008/2018 idioms (BLOCK, RESHAPE, derived types with FINAL, explicit interfaces, IMPLICIT NONE). Legacy stays. Python+Fortran split (pymacos) is already serving the "modern frontend" role.
- **Wrappers strategy**: language-neutral Fortran API as single source of truth; Python and MATLAB bridges over it. Pays back the factoring cost after ~3 backports.
- **FDTD**: out of scope for embedding in macos. Document handoff pattern with external tools.
- **GMI build default**: `gfortran` (clean MATLAB exit). `ifx` available via `makegmi.sh ifx` with `-reentrancy=none` linking the single-threaded `libifcore.so.5`.

---

## 11. GMI replacement surface

Inventory of what mmacos must cover to retire GMI as the de-facto
sensitivity / sensing entry point.  Source: survey of
`MACOS_resources/GMI/{GMI,GMIG}.F`, `GMI.inc`, `call_GMI.m`,
`optiixInit_jzlou.m`, and the `test_ff/` driver scripts on 2026-05-31.

GMI is monolithic: one MEX (`GMI.mexa64`), one work routine
(`GMI_DVR`), 14 positional PRHS inputs, 9 named PLHS outputs, and
~30 named behaviours toggled via a packed `pflg` flag array.  The
mmacos equivalent doesn't have to mirror the call shape — it can
expose the same capabilities via the existing `+macos.*` function-
package idiom — but the **set of capabilities** below is what
"retire GMI" means.

### 11.1 Top-level invocation shape (today)

```matlab
[PIX, CE, OPD, OPDMask, SPOT, WFE, c, metMeas, USER] = ...
    call_GMI(prb, pzern, pgrid, pdm, pfa, prad, pimg, ...
             InfFcnZern, InfFcnGrid, param[, winfil])
```

`param` is a struct (~30 fields, table below).  Each call: snapshot
nominal state → apply perturbations → trace → assemble outputs →
restore.  Per-call cost is dominated by the trace; the snapshot/
restore book-keeping is what makes a sensitivity loop coherent.

### 11.2 Perturbation channels

| Channel | Surface list | Coefficient vector | What it writes | Reset behavior |
|---|---|---|---|---|
| Zernike | `param.zernSrf` (Nx2: ids + pass-flag) | `pzern` (N·mzern) | `ZernCoef(4..nzern+3, iElt)`; **forces `SrfType=8`** | Idle ELSE zeros coefs + sets `SrfType=2` |
| MonZern (FreeForm) | `param.monzernSrf` | `param.pmonzern` (N·nmonzern) | `MonZernCoef(4..nmonzern+3, iElt)`; caller must already have set SrfType=14 | Zeros coefs |
| Rigid body (per-element) | `param.rbSrf` (Nx2: ids + global(0)/local(1) frame) | `prb` (mprb-vec) | `PERTURB`/`GPERTURB` via SMACOS — updates `TElt`, `psiElt`, `vptElt`, `pMon/xMon/yMon/zMon`, on FreeForm also `pFF/xFF/yFF/zFF`, `pData/xData/yData/zData` | `SetToNominalSettings` restores from snapshot |
| Rigid body (group) | `rbSrf` last col == 1 → group via `EltGrp` | `prb` | `GPERTURB` propagates to every element in the group | as above |
| Grid (surface deformation) | `param.gridSrf` | `pgrid` (N·mgrid²) | `GridMat(:,:, iEltToGridSrf(iElt))` per actuator pattern | Not snapshotted (overwritten each call) |
| DM (poked actuators) | `param.dmSrf` | `pdm` (mpdm-vec) | Influence-function-weighted grid deformation per actuator | as above |
| InfFcnZern | scoped by zernSrf | 15-vec Zernike influence function (input) | Scales the per-segment Zernike apply | — |
| InfFcnGrid | scoped by gridSrf | `mgrid x mgrid` grid influence (input) | Scales the per-segment grid apply | — |

Surface lists are 2-D matrices.  The **last column** is a flag (for
`rbSrf` it's the global/local-frame switch); a 1-column list silently
disables the channel.  This is the source of the common "wrote
`[1:N]'` and got no perturbation" footgun.

Compile-time sizes (`GMI.inc`): `numseg=7`, `mzern=45`, `mgrid=256`,
`mprb=(numseg+55)·6`, `mpzern=numseg·mzern`,
`mpgrid=mgrid²·((numseg+1)·2+numacf)`, `mpdm=(numseg+2)·2`, `mpfa=7`,
`mprad=numseg`, `mpimg=100`, `mpflg=2000`.  mmacos can size
dynamically — no recompile to grow.

### 11.3 Outputs (PLHS array, 10 entries returned by the mex; `call_GMI.m` re-packs to 9)

| Idx | Name | Shape | Meaning |
|---|---|---|---|
| 1 | `PIX` | nPix×nPix | Detector pixel array (after cPix rebin + jitter + shot/read noise + QE + bias + crosstalk) |
| 2 | `RealEF` | nPix×nPix | Real part of pixel-plane complex EF |
| 3 | `ImagEF` | nPix×nPix | Imag part — `CE = complex(RealEF, ImagEF)` MATLAB-side |
| 4 | `OPD` | nGridPts×nGridPts | OPD map at `ifOPD` reference element |
| 5 | `OPDMask` | nGridPts×nGridPts | Per-pixel mask (1 if traced, 0 if blocked) |
| 6 | `SPOT` | 2×iSpot | Spot diagram (x,y per ray hit) |
| 7 | `WFE` | 1×1 | RMS WFE in WaveUnits |
| 8 | `cent` | 2×1 | Centroid (xcent, ycent) |
| 9 | `metMeas` | nMetMeas×1 | Metrology sensor measurements (when `ifMetCalc`) |
| 10 | `USER` | mElt×1 | User-defined per-element scratch (currently RptElt diagnostic) |

### 11.4 Modes / flags (the `pflg` payload — parsed by `ExtractFlagParameters`)

Fixed-position scalars (`pflg(1..29)`):

| pos | param field | Meaning |
|---|---|---|
| 1 | `ifFEX(1)` | Find exit pupil; `==2` also writes XP geometry to pflg(31-37) |
| 2 | `ifPupilImg` | Compute pupil image |
| 3 | `cGrid` | Output grid sampling (≤ nGridpts) |
| 4 | `cPix` | Output pixel-plane sampling |
| 5 | `DMlim` | DM stroke clip |
| 6 | `ifOPD` | Element id where OPD is reported (also drives mask) |
| 7 | `ifShotNoise` | Add Poisson shot noise to PIX |
| 8 | `sigReadNoise` | Gaussian read noise σ |
| 9 | `sigJitterX` | Pixel jitter σ along x |
| 10 | `sigJitterY` | Pixel jitter σ along y |
| 11 | `sigCrosstalk` | Pixel-to-pixel crosstalk σ |
| 12 | `StartSeed` | RNG seed |
| 13 | `transMaskThreshold` | Translation-mask threshold |
| 14 | `rotMaskThreshold` | Rotation-mask threshold |
| 15 | `pixelSize` | Physical detector pixel pitch |
| 16 | `mzern` | # Zernike modes per surface (must match `GMI.inc:mzern`) |
| 17 | `QE` | Quantum efficiency |
| 18 | `DBias` | Detector bias |
| 19 | `wlens` | Wavelength selector |
| 20 | `ifPIX` | Run PIX command (else PIX returned empty) |
| 21 | `ifRetRefSrf` | Return reference-surface info |
| 22 | `ifSPOT` | Element id for SPOT diagram (0=off) |
| 23 | `ifPIXflip` | Flip PIX y-axis |
| 24 | `ifPIXSpotDetCheck` | Spot-vs-pixel-plane consistency check |
| 25 | `ifSysCalib` | Run `SystemCalib` pre-pass |
| 26 | `ifPIXElt` | Element id where PIX is sampled |
| 27 | `ifMetCalc` | Run `MetCalc` (metrology) |
| 28 | `ifSpfCalc` | Run `SpfCalc` |
| 29 | `ifRetUserSrf` | Populate `USER` array |

Variable-length blocks (starting at `pflg(30)`, sentinel `9999` =
absent): STOP(4) → iFSM+TFSM → iFDP → gridSrf → zernSrf → dmSrf →
rbSrf → RptSrf+RptElt → ifFEX(2..8) tail → monzernSrf → RefSurfs →
INTsrf.  Order is positional; adding fields without updating both
sides shifts every later index.

### 11.5 Internal work routines (the "what GMI does") — coverage status

| Subroutine | Purpose | mmacos coverage today |
|---|---|---|
| `ReadDebugEnv` | Read `GMI_DEBUG` env var | trivial — env-based debug flag |
| `GMI_DVR` | Top-level driver | what dw_dx (Phase 7) starts to replace |
| `ObtainNominalSettings` | First-call snapshot of: Wavelen, Flux, zElt, ChfRayPos/Dir, xGrid/yGrid/zGrid, Tout; per-elt KrElt/KcElt/nObs/ObsType/ObsVec/psiElt/vptElt/rptElt/TElt/pMon-xMon-yMon-zMon/pFF-xFF-yFF-zFF/pData-xData-yData-zData/SrfMetPos | **MISSING** — need `+macos.snapshot()` / `+macos.restore()` |
| `SetToNominalSettings` | Restore from snapshot at start of every later call | same |
| `ExtractFlagParameters` | Unpack `pflg` into named locals | Absorbed by named-arg `+macos` API — the flag-packing layer goes away |
| `ApplyPerturbationToOpticalSystem` | Loop the three channels (Zern, MonZern, rigid-body) | Channel-wise: `+macos.perturb`, `+macos.perturb_many`, `+macos.perturb_grp`, and `mmacos('elt_srf_zrn_coef', ...)` cover the apply; the loop+reset structure is what dw_dx builds |
| `ApplyPerturbationToOpticalSystem_Prb` | Per-element rigid-body loop | `+macos.perturb_many` ✓ |
| `SystemCalib` | Pre-pass calibration (when `ifSysCalib`) | **MISSING** — investigate what it actually does |
| `MetCalc` | Metrology sensor calc | **MISSING** |
| `SpfCalc` | "SPF" calc (what?) | **MISSING** — figure out scope |
| `ExecFDPCmd` | Run FDP command | **MISSING** |
| `ComposePIX` | Assemble PIX (cPix rebin + noise + jitter + QE + bias + crosstalk) | **MISSING** — substantial port; raw mex has the buffer but no noise pipeline |

### 11.6 Capability coverage matrix

> **Parity audit 2026-06-12.**  mmacos has **full engine-level parity**
> with pymacos: every implemented pymacos function maps to a shared
> `macos_api_mod` routine that mmacos also exposes (the raw command
> layer is 104 commands vs pymacos's ~90 implemented functions; the two
> pymacos `*ActivePointSrc` functions are `#ToDo` stubs).  Remaining
> gaps are **convenience `+macos` veneers over already-wired raw
> commands**, not capability.  First batch closed 2026-06-12 (spot,
> fex/xp, Kc/Kr — below).  Still veneer-only: grid-surface family
> (`elt_srf_grid_*`), `src_size`/`src_csys`/`src_info`, surface-csys,
> `elt_zrn_type`/`norm_rad`, `set_ray_info`, `elt_grp` query helpers.

| Capability | Status in mmacos |
|---|---|
| Perturbation channels: rigid-body | `+macos.perturb`, `+macos.perturb_many`, `+macos.perturb_src` ✓ |
| Perturbation channels: Zernike | Raw `mmacos('elt_srf_zrn_coef', ...)` ✓; no `+macos.elt_zrn_*` veneer yet |
| Perturbation channels: MonZern (FreeForm) | Raw `mmacos('elt_srf_mon_zrn_coef', ...)` ✓; no veneer |
| Perturbation channels: grid surface | Raw `mmacos('elt_srf_grid_data', ...)` ✓; no `+macos.elt_grid` veneer (column-major transpose convention needs care) |
| Perturbation channels: DM (influence-function-weighted grid) | **MISSING** — needs InfFcnGrid + per-actuator grid build helper |
| Group perturbations (EltGrp) | `+macos.perturb_grp` planned — uses `prb_elt_grp` cmd that's already wired |
| Nominal save/restore (sensitivity-loop coherence) | **MISSING** — top of Phase 7 list |
| OPD output | `+macos.opd()` ✓ |
| Complex EF output | `+macos.complex_field(srf)` ✓ |
| OPD mask | Returnable via `macos.get_ray_info().ok_pass` reshape; no `macos.opd_mask()` shortcut yet |
| Spot diagram | `+macos.spot(srf)` ✓ (2026-06-12) — struct: pts/centroid/shift/csys/nSpot; over `spot_cmd`+`spot_get` |
| PIX output (raw) | Hand-written `mmacos('apodize', ...)` writes to WFElt; no PIX pull-back yet |
| PIX with noise model (ComposePIX) | **MISSING** — port from GMI's ComposePIX (cPix rebin, shot/read noise, jitter, crosstalk, QE, bias) |
| WFE | Returned by `+macos.trace()` as `.rmsWFE` ✓ |
| Centroid | Computable from `get_ray_info`; no `+macos.centroid()` shortcut |
| FEX / exit pupil | `+macos.fex(mode)` + `+macos.get_xp` / `+macos.set_xp` ✓ (2026-06-12) over `xp_fnd`/`xp_get`/`xp_set` |
| Conic / radius (Kc, Kr) | `+macos.get_elt_kc` / `set_elt_kc` / `get_elt_kr` / `set_elt_kr` ✓ (2026-06-12) |
| Metrology calc | **MISSING** |
| FDP / FSM | **MISSING** |
| System calibration | **MISSING** |
| Reference-surface info return | partial — `ifRetRefSrf` triggers a structured output not yet exposed |
| Per-element USER scratch | low priority |
| **Design optimization (CALIB)** | **§3.5 Phase 1c ✓** — `+macos.calib` + `calib_set_*` family wraps `nls_optim_dvr` (the unconstrained LM path).  Constrained NPSOL path is opt-dev-only via macos's existing dispatch; surfacing it through wrappers is Phase 1d. |

### 11.7 Phased roadmap to retire GMI

| Phase | Delivers | GMI coverage |
|---|---|---|
| **Phase 7 (dw/dx) — must include:** | Nominal save/restore (`+macos.snapshot`/`+macos.restore`); Zern/MonZern/rb perturb channels (already exist as commands; build the loop+reset structure); OPD/CE pulls (exist); SystemState / chain-runner that doubles as the `test_coro_aberrations.py` infrastructure deferred from Phase 5 slice 7 | ~60% of GMI's everyday use (sensitivity loops, EFC inputs) |
| **Phase 7+ slice A** | Spot veneer; PIX-with-noise (port ComposePIX); FEX wrapper; OPD mask veneer; centroid veneer | ~25% more |
| **Phase 7+ slice B** | DM influence-function wrapper; grid-surface `elt_grid` veneer (with the column-major transpose convention pinned in tests) | ~5% more |
| **Phase 7+ slice C** | SystemCalib, MetCalc, SpfCalc, FDP, FSM | Final ~10% — specialized (closed-loop calibration demos, metrology workflows) |

At end of Phase 7 + slices A–C, mmacos covers GMI's full surface; the
GMI mex can be deprecated.  Existing GMI callers either migrate to
the `+macos` package directly or use a thin compatibility shim
(callable `mmacos_gmi_compat(prb, pzern, ...)` that maps the 14-arg
positional form to the modern `+macos` calls — useful for grandfathering
the corpus of test_gmi*.m scripts).

> **Recommended GMI-replacement user surface:** the design layer's
> `System.from_rx` → `vary` → `sensitivities` flow
> (PLAN_DESIGN_LAYER §1.0 / Sprint 2A-i) is the package-level
> successor to the GMI sensitivity-loop workflow — import a
> prescription, declare DOFs, get a Jacobian, all over the same
> bitwise-verified `dw_dx` channel machinery the `+macos` package
> already exposes.  The raw `+macos` channel calls and the
> `mmacos_gmi_compat` shim remain for power users and grandfathered
> scripts; new users get pointed at the design-layer surface.

### 11.8 Open questions

- `SpfCalc` purpose — comment in GMI.F is sparse.  Need to read the body or ask jzlou.
- `SystemCalib` operating mode — is it a closed-form pre-pass or an iterative refinement?  Affects whether mmacos can do it inline or needs its own driver.
- DM influence-function format — `InfFcnGrid` is `mgrid × mgrid`.  Single-actuator pattern that gets convolved per actuator?  Or per-actuator stack?  Affects the DM apply API shape.
