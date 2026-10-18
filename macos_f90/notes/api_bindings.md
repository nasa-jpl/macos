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

`srt_trace(OK, xoff, yoff, iEndElt, Pos_, Dir_, IndPrev_, IndCur_, RayOK_, Status_,
FailElt_, nElt_)` (2026-10-17, BRIEF_to_bindings B.2, the first single-ray routine) is the CLI's
SRT made SMACOS-dispatched (LoadStack already had its entry; only `#ifdef MACOS_CMD` kept it out):
one ray from a normalized source-aperture point, its crossing and outgoing direction at each
element from traceutil_mod's `thisRay`.  SRTrace's CSMACOS dump (a `status="replace"` write to
`filnam` = the deck's path without `.in`) was removed with it.  The surface normal SRTrace saves
is NOT returned: CTRACE updates `Nsrfm` only where it reflects or refracts, so other columns
hold the last ray's.  An OBSCURED ray keeps `LRayOK` (it is traced on; `Status_` = Obscured).
The veneers (`macos.srt`, pymacos `srt`) mark the trace modified: SRTrace re-seeds the source
ray set.  Gates tSrt / test_srt.py (the CLI's printed numbers + the paraboloid sag + the axial
focus + a full trace unchanged after it).

`ray_trace(OK, iRay, ...)` (the CLI's RAY) and `segray_trace(OK, iSeg, ...)` (SEGRAYTRACE; the
ray through Segment element iSeg's `RptElt`, forced onto that segment) return the same buffers
from `thisRay`: under SMACOS each handler calls `ThisRayAlloc(nElt)` and sets `ifSRT` around its
one CTRACE; the CLI path prints as before.  **A SMACOS handler must not loop back to its prompt**:
the prompt re-reads an EMPTY stack and keeps the old value, so RAY's out-of-range re-prompt and
every `GO TO 11` in SEGRAYTRACE hung the host -- both now re-prompt under `MACOS_CMD` only.
segray_trace refuses a non-Segment element before the call and clears `thisRay`, so a refusal
inside the handler cannot return the previous call's ray.  For the Segment elements the ray does
not meet, a column repeats the previous state (the source point before the hit segment, the hit
segment's after).  Found on the way: CTRACE clears `ifRayPrint` on return, and SEGRAYTRACE set it
once before its loop -- a second segment in one command traced silently (journal leg added).
Gates tSegRay / test_segray.py (the CLI's numbers for segments 3 and 7; every segment's crossing
on the chief-parallel line through its `RptElt`, from the deck's lines).

`ref_rays(OK, iKind, RefRay_, Off_, nRef_)` (RRAY iKind=1, PRAY iKind=2; nRef_ = mPupRay = 9)
returns `refRay` and the new `refRayOff(2,mPupRay)` (elt_mod; each found ray's XY offset from the
chief ray in the source grid frame -- FSR returns it in `srp`, which the handlers overwrote).
`-1` = no ray found, `0` = a slot the kind does not fill.  Engine defects fixed with it: PRAY's
first ring point read `srp(1)=9.0` (ray 2 was RRAY's +x edge ray); RRAY left slots 6..9 holding a
previous PRAY's rays and TPR traced them; FSR left `ifRayPrint=.TRUE.`, so the next OPD printed
every ray's segments (CLI and bindings alike) -- FSR now restores the caller's value.  **The
single-ray api routines restore `nGridPts`**: ColSource sets the global `nGridPts = npts` (one
short) and only `SUBROUTINE OPD` puts it back, so after RAY / SEGRAYTRACE / RRAY / PRAY the
`get_src_sampling` / `opd_val` / `opd_mask_get` callers saw nGridPts-1 and pymacos `opd()` /
mmacos `macos.opd` returned a crop of the last map (SRT does not leak it).  ray_trace,
segray_trace and ref_rays save it before the SMACOS call and restore it after.  Gates tRefRays /
test_ref_rays.py (the CLI's ids and offsets on Rx_Cass_FarField; each ray within two grid
spacings of its edge / 0.9-ring point from the deck's Aperture=; the ring's mirror symmetry;
each offset tied to its index through `ray()`; RRAY's cleared slots; the map read WITHOUT a
re-trace; no ray printing after) and the `*_leaves_the_session_grid` legs in tSrt / tSegRay /
test_srt.py / test_segray.py.

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

## One default per option: a veneer never carries a default an api getter can supply (2026-10-16)

Dave's rule (2026-10-10, PLAN_CONSOLIDATION 1.4a; landed by BRIEF_to_bindings A.6):
every session/user option the bindings expose has ONE default, held by the
ENGINE (set at allocation where it is engine state); the CLI, mmacos and pymacos
inherit it.  A veneer that needs the value reads the api getter (`modes_get`,
`opd_ref_get`, `ffcut_get`, ...); an omitted argument means "the engine's", never
a literal restated in the veneer.  Where the api has NO getter (spot centre
`spcOption`, the OBS option `iObsOpt`, model size) the veneer's literal is a
restatement that can drift: the table marks it **agrees** and the fix is an api
getter (a new routine -- Dave's call), not a second literal.

The table is cmdref `25_defaults.md`, GENERATED by
`MACOS_resources/mmacos/src/gen_mex_wrappers.py --targets defaults` from
`binding_overrides.OPTIONS`: every quoted default is READ from its source (engine
allocation site, mmacos `.m` arguments block, pymacos signature) at generation
time, and a pattern that stops matching stops the generator -- the table cannot
keep describing a default the code no longer has.  Gates: mmacos `tDefaults` +
pymacos `test_defaults.py`, the same decks and the SAME numbers in both.

Drift found and fixed 2026-10-16 (each by REMOVING the restatement):
- the eight mmacos dw_* harvests defaulted `opd_ref='chief'` and re-applied it
  after their reload: on `tst_opdref_mean.in` (`UseChfRay4OPD= N`) the harvest's
  `w_nom` was the CHIEF map (|w_nom - opd()| up to 2.518e-12 against the deck's
  mean) and the session was left on chief.  Default now '' = the deck's own.
- `fex` / `sxp` (both bindings) defaulted `mode=1`; the CLI's FEX reads the
  session's CENTROID/CHIEFRAY mode.  Default now = `modes_get().centroid`.
  (mmacos `sxp` also refused `mode=0`: `mustBePositive`.)
Dave's rulings (2026-10-16) on the three routed items: (1) model size defaults to
256 everywhere -- the CLI's start-up prompt already took 256 on Return
(`default_problem_size`, macos_mod.F; the first table said "no default", a misread
of the commented-out lines), mmacos `init` gains `= 256`, pymacos `init` was 512;
(2) api getters `spc_option_get` and `obs_opt_get` added, the spot centre and the
OBS option now inherit, and `spot_cmd` restores `spcOption` after its SPOT (the
`xp_fnd` / `ifCentroid` pattern: a per-call centre never changes the session); (3)
pymacos `spot()` returns the api's 4-vector shift [chief xy, vertex xy], the
mmacos convention.

