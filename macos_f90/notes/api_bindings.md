# MACOS engine notes — macos_api_mod and the bindings' tests: the wrapper layer, pymacos/PROPER regression, the Obscuring-mask rule

> **Read this file before touching:** `macos_api_mod.F90`, pymacos/mmacos veneers, coronagraph fixtures.  These notes are the engine record
> (what moved, the measurements, the gates); they are NOT auto-injected — the index is
> `macos_f90/CLAUDE.md`, the root rules `../../CLAUDE.md`.  Text is moved verbatim from the
> former monolithic cheatsheet (2026-10-15 split); add new entries here, not to the index.

## Regression tests via pymacos (~/dev/MACOS_resources/pymacos/)
- CodeV-comparison tests in `tests/test_api_rx_grating.py` and `tests/test_masks.py`
  cover the geometric / ray-trace paths (6601 tests, all passing).
- PROPER-comparison tests in `tests/proper_compare/` cover the physical-optics paths
  (INT/PIX/DFT-propagation) that CodeV can't validate. Organised by phase:
  - **Phase 1** (`test_cass_ff*.py`, results in `results_phase1/`): far-field
    image-plane PSF comparison on `Rx_Cass_FarField.in`. Nominal + SM Tx/Ty/Tz
    perturbations.  Peak-normalised agreement at 1e-11 with OPD pass-through.
  - **Phase 2** (`test_coro_nfprop.py`, results in `results_phase2/`):
    near-field plane-to-plane propagation between Elt 2 and Elt 3 of
    `Rx_Coro.in` (HCIT-style coronagraph, simplified-conic version).
    Sum-normalised agreement at 2.4e-14 RMS, 4.8e-13 max -- at
    double-precision FFT round-off (needs the runtime `dx_at()` query
    described below + odd `nGridpts=511`).
  - **Phase 3 + 4 + 5** (`test_coro_nfprop_phase3.py`, results in
    `results_phase3/`): six additional Coro steps:
      - 3a NFPlane 5→6 (sampling-limited 3.7e-8)
      - 3b sphere→plane 8→9 (4.2e-11)
      - 4a pupil→pupil 8→10 through focus (5.9e-13)
      - 4b NFPlane 13→14 (7.4e-5; localized post-focus diffraction ring)
      - 5.1 ExitPupil 20→FocalPlane 21, no mask (2.3e-9)
      - 5.2 same with FPM=400 um + Lyot=14 mm coronagraph
        (2.6e-6 Strehl-norm against a peak suppressed by 3.2 million
        relative to the un-coronagraphed baseline).
  - **Phase 6a** (`test_coro_apodizer.py` + `apodizer.py` + new
    `pymacos.apodize` wrapper): pupil apodisation via a user-supplied
    NxN amplitude mask, multiplied into macos's WFElt and PROPER's
    wavefront via the same numpy array.  Apodiser construction uses
    KxK super-sampling for sub-pixel area-weighted aperture edges so
    low-N results converge to high-N as the same physical shape.
    First test (Gaussian-edge pupil at Elt 5, NFPlane to Elt 6):
    4.0e-8 max, 1.1e-9 RMS (same sampling-limited regime as 3a,
    confirms the wrapper integrates correctly with downstream
    propagation).
  - **Contrast scoring** (`test_coro_contrast_curve.py` +
    `contrast.py`): radial dark-zone contrast vs lambda/D at the
    science focal plane.  lambda/D derived empirically from the
    no-mask PSF's first Airy null (1.22 lambda/D with a fractional-
    depth guard).  Phase 5.2 dark-zone contrast hits ~3e-10 in the
    7-10 lambda/D range; bright outer-ring artefact at ~15 lambda/D
    (Lyot edge diffraction).  Baseline for Phase 6b/c apodised
    designs to score against.
- Run via `./run_proper_tests.sh` at the pymacos root.  It rebuilds pymacosf90
  (if needed) and invokes each phase in its own pytest process to dodge a
  pymacos state-leak across model_size transitions (512 ↔ 1024).  Artefacts
  go into per-phase `results_phaseN/` directories.
- The pymacos<->PROPER coupling needs three reconciliations to reach the
  observed agreement:
  - **Aperture match.** PROPER takes amplitude DIRECTLY from macos's mask
    via `prop_multiply`, NOT from an analytical circular_aperture +
    obscurations model.  Mismatched apertures put light where macos has
    zeroed it out, breaking the wavefront tilt and halving the apparent
    PSF shift under perturbation.
  - **Sign flip.** macos OPD sign convention is opposite to PROPER's
    `prop_add_phase` input.  Default `opd_sign_flip=True` reconciles them.
  - **Normalisation choice.** `compare_and_record` takes `norm_kind=`
    `'peak'` (Strehl-norm, default; right for image-plane PSFs) or `'sum'`
    (flux-norm; right for pupil-plane / near-field intensities).  Peak-norm
    inflates the residual on flat-top NF PSFs (where peak position inside
    a uniform region is noise-dominated); sum-norm gives the physically-
    meaningful flux-conservation precision.
- Centroid-based (not peak-based) alignment in `compare_and_record`:
  intensity-weighted center of mass is robust for both sharp Airy peaks
  (Phase 1) and flat-top NF pillars (Phase 2).

## macos_api_mod: language-neutral wrapper layer for pymacos + mmacos
*(The routine catalog is cmdref Part II (generated from the api source); the bindings program is PLAN_CONSOLIDATION §1.4a/b, BRIEF_to_bindings.md.  Below: the layer's shape and its porting record.)*

`macos_f90/macos_api_mod.F90` (compiled into `libsmacos.a`) is the
shared SMACOS-call wrapper that both `pymacos` (Python) and `mmacos`
(MATLAB) bind into.  Lives in `MODULE macos_api_mod`; holds the SMACOS
call-line buffer (`command`, `CARG/DARG/IARG/LARG/RARG`), the shared
output arrays (`PixArray`, `OPDMat`, `RaySpot`, `OPD`, `SPOT`, `PIX`,
`USER`), the package-level state flags (`firstEntry`, `rxLoaded`), and
~96 wrapper subroutines covering init, load_rx, save_rx, modified_rx,
n_elt, trace_rays, opd_val, spot_cmd/get, int_cmd/get, cfield_cmd/get,
cfield_apodize, elt_dx_get, base_unit_to_metres, perturb_elt,
prb_elt(_grp), the elt_* Rx-inspection/modification family, the
source-side set_src_*/get_src_*, the elt_grp_* family, stop_info_*,
xp_get/set/fnd, sxp_fnd, ors_run, srs_run, plus internal helpers
(`SystemCheck`, `EltRangeChk`, `StatusChk1`, `translateSurfaceID`,
`translateEltID`, `checkSurfaceID`, `checkEltID`).

`ray_status_get(OK, RayStatus_, RayFailElt_, nRays)` (added 2026-06-12,
Q2) exposes `elt_mod`'s per-ray `RayStatus(:)` / `RayFailElt(:)`
(`RayStat_OK/Obscured/Miss/Bracket/MaxIter/Undef`) verbatim from the
last trace — the per-category complement to `ray_info_get`'s binary
LRayOK/LRayPass.  mmacos veneer: `get_ray_status(N)`.  Backs the design
layer's ray-loss guard (`PLAN_DESIGN_LAYER.md` §1.3.4).

The MATLAB **design layer** (`macos.design.System`,
`MACOS_resources/mmacos/src/+macos/+design/`) is the first heavy consumer of
this wrapper surface — import (`from_rx`) → `sensitivities` (harvests the
Phase 7 `dw_dx` channels) → `vary`/`evaluate`/`optimize`.  It adds no
Fortran; it's a pure MATLAB layer over these routines.  Note: SLSQP
objects live in `libslsqplib.a`, a separate archive from `libsmacos.a`
— any binding that links libsmacos must also link slsqplib (the SLSQP
back end references `slsqp_`).

Promotion history: until §5.2 the file lived in
`MACOS_resources/pymacos/src/cmake/source/macos_api_mod.F90` and was
compiled separately by each binding's own CMake/Make build.  Now it
lives in `macos_f90/` and the binding builds just link `libsmacos.a`
+ pull the `.mod` from `mod_smacos/`.  pymacos's `CMakeLists.txt`
takes `MACOS_BUILD_DIR` as a cache var (default
`~/dev/macos/build_release`); mmacos's Makefile takes the same.

gfortran cleanliness: the file is portable across ifx and gfortran
as of §5.2.  Notable porting fixes — replaced `(/'m','cm',...'/)`
implicit-length char arrays with explicit `character(len=4)` forms;
converted all `do concurrent` to plain `DO` loops (some inside
`BLOCK` for index-decl scoping); 50× `LOGICAL == PASS/FAIL` sites
swapped to `.eqv.`; 16× INTEGER-comparison sites switched to integer
equality (`var == 0` / `var /= 0`) where the variable was INTEGER;
LOGICAL-array `ifGlobal` tests rewritten in `prb_elt_grp`; bare
LOGICAL `filter` test in `set_src_csys`.  pymacos pytest 6601/6601
+ PROPER-compare green under both compilers; mmacos 11/11 under
both.


## pymacos intensity() / complex_field() / dx_at() / apodize() wrappers
> These span `macos_api_mod.F90` + `pymacos_f2py.f90` + `macos.py`.
> The Python-facing detail already lives in
> `../MACOS_resources/pymacos/CLAUDE.md` — keep it there (single source
> of truth) and reference it; do not duplicate. The buffer-wrapper
> *template* (spot_cmd/spot_get pattern) belongs with macos_api_mod above.

## Coronagraph mask element type (silent failure mode)
*(Documented: Manual §6 "Setting up Diffraction Propagations" (masks) and §4 Reference Surfaces.  Below: how it was found.)*
Circular obscurations declared in a prescription (`nObs`, `ObsType=Circle`,
`ObsVec=...`) are **only applied to the diffraction-grid wavefront when the
element is `Element= Obscuring`**.  If the element is `Element= Reference`
(or any non-Obscuring type), macos still parses the obscuration metadata
and applies it to geometric rays during ray tracing, but the
diffraction-grid `WFElt` array sails through untouched -- a hard-edge
"mask" that's invisible to the diffraction propagation.

Caught while building the Phase 5 coronagraph test
(`tests/proper_compare/test_coro_nfprop_phase3.py`).  Original
`Rx_Coro_FPM.in` had Elt 9 (the FPM) as `Element= Reference` with a
132 um circular obscuration -- it looked like a working coronagraph
in the prescription, but `complex_field(9)` and `intensity(9)` were
byte-identical with and without the mask.  Suppression at the science
focal plane was ~17%, all of which came from flux trimming at the
Lyot stop (Elt 14) -- not the FPM.  After changing to
`Element= Obscuring`, the FPM finally bites: on-axis suppression at
Elt 21 became factor ~1e5 from FPM alone, factor ~3e6 with the Lyot.

The right model to follow is `MACOS_resources/docs/macos-manual/examples/CoroExample.in`
Elt 6 (the working CoroMask).  Any new coronagraph test prescription
should use `Element= Obscuring` for the mask element.

