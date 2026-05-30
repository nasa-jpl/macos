# MACOS Development Plan

Working document. Tasks grouped by thrust; checkboxes track state. Items not
yet started have empty `[ ]`; completed `[x]`. Discrete subtasks under each
task.

---

## 0. Hygiene / quick fixes

- [x] Zernike trace-dispatch ELSE branches in `propsub.F`, `srtrace.F`, `tracesub.F`; bring all three into parity, modernize to `ZernType_*` constants. Test: any Rx with `ZernType= Noll` and a non-zero coefficient should produce non-zero OPD response.
- [x] `fk.c` unconditionally in `libsmacos.a` (drop `BUILD_GMI` gate). Makes `makems.sh + makegmi.sh` work standalone.
- [x] Build scripts auto-bootstrap bundled readline on first clone. `makems.sh`, `makegfortran.sh`, `makeall.sh` detect missing `libreadline.a` and run `./configure -q && make` in `macos_f90/readline-8.2/`.
- [x] `dopt_init_vars` stop clobbering meaningful defaults of `nitrs_dopt`, `dopt_tol`, `SvdSvCut`.
- [x] Commit `ZGD_test_files/` consolidated corpus (joint-dev pickups + local additions, drop `macos_param.txt` symlink).
- [x] Drop marker-based gate from `run_regression.sh`; trust MATLAB exit code now that `-reentrancy=none` is on the ifx mex link.
- [ ] `define_local_csys` follow-up from the `develop_STOP` branch head — small robustness / speed improvement not yet on release-candidate.
- [ ] Audit `dopt_init_vars` for any other meaningful-default-clobbered-by-zeros collisions beyond the three found today.
- [ ] Quiet the `MBFile6: Unidentified string` warning when a parser hits a target-specific keyword (e.g. `OptBeamPos=`) under a non-matching `OptTarget`. Treat as "not relevant in this mode" rather than "unparseable junk."
- [ ] **Fix `init()` re-init heap corruption on model_size transitions.**
  When `macos_init_all(new_size)` is called in a process that already
  ran a trace at a different size, the next FFT-bearing trace aborts
  in `malloc()` / `free()` (`invalid size (unsorted)`, `unaligned tcache
  chunk`, `munmap_chunk: invalid pointer`).  Surfaced in mmacos's
  matlab.unittest run when Phase 5 PROPER tests (model_size=512)
  followed Phase 3/4 tests (model_size=128) in the same process.
  Pymacos has the same bug; its `run_proper_tests.sh` works around
  it by invoking a separate pytest process per phase.  mmacos's
  `run_mmacos_tests.sh` does the analogous split into per-size matlab
  -batch invocations (see §5.4 Phase 5 notes).  Workaround works but
  is unsatisfying — the underlying engine reallocation in
  `macos_init_all` (and the routines it cascades into:
  `src_mod_init_vars`, `elt_mod_init_vars`, the FFT-buffer allocators
  in `fftsub`, etc.) is leaving dangling pointers somewhere that get
  freed at a wrong size on the next trace.  Likely candidates:
  (a) module-level allocatables in some `*_mod_init_vars` that
  `deallocate` based on the OLD size; (b) `iEltToiWF` / FFT plan
  caches that key on size but free with `mFFT` from a stale module
  param.  When fixed, drop the split in `run_mmacos_tests.sh` and
  put the full suite back in one matlab -batch.
- [ ] Renormalize `psiElt` after the `Q·psi` rotation in `CPERTURB_PROG` (funcsub.F:349-350). Currently a round-trip `perturb(+θ) + perturb(-θ)` along a single axis leaves `psi` off by 1 ULP because `sin²(θ) + cos²(θ)` ≠ 1 exactly in IEEE 754 for some specific θ values (e.g. 1e-6, 3e-5). The artifact is at the eps × |coord| floor — invisible against any practical signal — but causes psi to drift slowly under many repeated perturbs and produces a ~3e-14 OPD round-trip residual that briefly confused the §5.4 Phase 2 +macos smoke-test author. One-line fix: `psi = psi / norm2(psi)` after the rotation. See `MACOS_resources/mmacos/test_state_after_roundtrip.m` for a regression probe.

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
  - [ ] Extend `design_optim_mod` to evaluate at a fan of source positions
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
- [ ] Implement or remove `WFE_ZMODE_TARGET` and `OPL_TARGET`
  - [ ] Either wire up `objfun_nom` branches for both
  - [ ] Or remove from the enum + parser so users can't silently pick a non-functional mode
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

8-12 new tests, each becoming an `.in` / `.jou` pair in `ZGD_test_files/` and a regression entry. Covers every `OptTarget` × DOF-type combination.

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
  - [ ] Decide: implement consumer (edge-sensor-distance measurement command, Hx-matrix output) OR remove the keyword to stop silently accepting data
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
  - [ ] **Slices 6+** for the remaining proper_compare files:
    - `test_coro_aberrations.py` (2 parametrized tests; Coro chain
      under SM tilts/translations).
    - `test_coro_dm_phase.py` (1 parametrized; DM phase imprint).
    - `test_coro_dm_grid_self.py` (1 parametrized; DM-grid
      self-consistency).
    - `test_band_limited_mask.py` (5 tests; band-limited FPM
      construction).
    - `test_psf.py` (1 real test, 1 skipped; PROPER-only baseline
      that duplicates the Cass-FF sanity we already have — low
      priority).
    - `run_broadband_*.py` are not tests; they generate reports.
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
  - [ ] `+macos/dw_dz_zernike.m` plus the `+channels/` submodule
    (ZernikeCoefChannel only).
  - [ ] Driver takes a `MacosSession` instance (mirrors pymacos's
    `MacosModel` arg); a free-function form taking no session also
    exposed for raw-mex users.
  - [ ] Writes the same `.mat` format pymacos's `dw_dz_zernike` does
    so downstream consumers don't fork.

- [ ] **Phase 7 — `dw_dx` driver**
  - [ ] `+macos/dw_dx.m` plus the full channel system as MATLAB
    classdef: `RigidBodyChannel`, `FocalPlaneChannel`,
    `GroupedRigidBodyChannel`, `SourceChannel`,
    `predict_global_rigid_response`, `group_W` synthesis.
  - [ ] Biggest port — channels carry real semantic content beyond a
    thin mex wrapper.
  - [ ] Same dual surface: `MacosSession`-driven default, raw-mex
    constructor for power users.

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
- [ ] Sensitivity-matrix workflow (Thrust A headline)
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

## 9. Branch / repo policy

- `release-candidate` is the slim public-facing branch (no NPSOL, no in-flight work)
- `opt-dev` adds the USE_NPSOL constrained-optim path
- New feature work lands on feature branches → merges into `release-candidate` after testing
- `main` lags `release-candidate`; sync periodically
- Test files duplicated where needed: not all macos users access `MACOS_resources`, but all `MACOS_resources` users see `macos`; manual examples and reader exercises must be in `macos` directly; some coronagraph Rxes remain private and stay in `MACOS_resources` only

---

## 10. Strategic decisions already made

- **C-rewrite**: no. Modernize in-place via Fortran 2008/2018 idioms (BLOCK, RESHAPE, derived types with FINAL, explicit interfaces, IMPLICIT NONE). Legacy stays. Python+Fortran split (pymacos) is already serving the "modern frontend" role.
- **Wrappers strategy**: language-neutral Fortran API as single source of truth; Python and MATLAB bridges over it. Pays back the factoring cost after ~3 backports.
- **FDTD**: out of scope for embedding in macos. Document handoff pattern with external tools.
- **GMI build default**: `gfortran` (clean MATLAB exit). `ifx` available via `makegmi.sh ifx` with `-reentrancy=none` linking the single-threaded `libifcore.so.5`.
