# Current Slice — in-flight working state

> **CC:** This file is the *ephemeral* working memory for the ONE sprint
> slice in progress.  It is the deliberate complement to the permanent
> record: `PLAN.md` / `PLAN_DESIGN_LAYER.md` hold *landed* state (sprint
> checkboxes, `CORE COMPLETE` blockquotes, the §10 Decisions ledger);
> the agent `MEMORY.md` holds durable learnings; `CLAUDE.md` (+ nested)
> holds mechanical gotchas.  THIS file holds only the half-done middle
> that compaction throws away — and it is **promoted then cleared** when
> the slice lands.  It is never a second source of truth.

> **CC (post-compaction / session resume):** read this file FIRST, then
> the plan section it anchors, then the CLAUDE.md set per the root
> directive.  If this file is at the empty template below, no slice is
> in flight — pick the next unchecked item from the active sprint.

> **CC (standing conventions, do not relax here):** work lands on
> `sls-dev` (macos) + companion `sls-dev` (MACOS_resources); every new
> function ships with a matlab.unittest test; `./run_mmacos_tests.sh
> fast` between edits, full suite pre-commit; every `matlab -batch`
> ends in `exit(0)`; each sprint tag ships a runnable worked-example.

---

## Active slice

**2026-07-05: centered-TMA FP EXTRACTION day (MACOS_resources only, no
engine change) — see [[project-fold-extraction]] memory.**  Dave's
recipe, all real-ray verified: the coaxial Korsch FP (on axis at
z=−4.07, swallowing the M2→M3 science beam) is extracted by (1) FIELD
BIAS = source tilt (10′ clears the honest 0.25 m FP body off the
science beam; 5′ doesn't), (2) a FLAT FOLD **in the M2→M3 FEED** behind
the PM, 90° into +x with psi in the X-Z plane → M3 + image + FP on a
**flat X-Y bench** at z≈0.87 behind the primary, (3) M3 back 2 m (j18's
own spacings put the exit pupil 1 m in FRONT of the PM — no fold fits).
- Builder (Telescope.m): `add_fold` (reflect-downstream isometry,
  WFE-neutral 6e-16), `set_hole` (perforated-primary clearance),
  `center_focal_plane`, `add_focal_plane('ap_r')`, ApType=None policy
  for folds+FP (ap_r = check_clipping BODY not a stop — honest FP as a
  stop made CALIB rigid trials lose all rays), optimize writes conics/R
  back to the mirror list (re-resolve keeps them), powered-only DOF
  enrollment, `view_orthoviews` 'zoom' detail panel + XY cross-plane
  ray reconstruction (both fans via fan_pt_ — beam extents measurable),
  tighter emit_ source standoff.  design/src: `fold_station_report`.
- Examples (Dave's finder/optimize pattern): `tma_centered_fold_search`
  (bias×dM3 ladders; placement rules M1_KEEPOUT 0.6 m / ARM_MIN 0.4 m /
  shortest backbone; CHOSE dM3=2, 15′, fold z=0.871, gap 0.221, shroud
  0.99×D, all-clear but M2's own shadow) + `tma_centered_foldfp`
  (ROC+conic ring balance — rigid DOFs re-point INTO the bias and eat
  the clearance, check_clipping catches it; −tilt 1.87 waves @2.3 µm =
  the honest extraction price; buy-downs noted).  6 new tests.
- **Astig round (pm, `382a3d6`):** margin→0.08 (Dave OK) → **10′ bias**
  geometry; the astig Dave saw after FEX = axis-solved conics evaluated
  at the bias; **fix = conics+ROC solved AT the bias field only
  (annular anastigmat): 30.6 → 0.061 waves @1 µm**; then **M1 stop-
  surface freeform (Dave OK'd M1 Zernike), modes 5:15 → 0.046 @1 µm /
  0.020 @2.3 µm, FOV unchanged** (`91d47af`).  Freeform traps printed
  in the example: ring-balance trades the center; 2-mirror single-field
  over-fit; **lMon must be the BEAM footprint, not body ap_r —
  `'lmon'` option**.  FOV: blur (−tilt) DL @2.3 µm over ~0.5′-dia core.
- **Plot rounds (`91d47af` + `4f422a9`):** Return retro verified exact
  (dir·−dir = 1.0; rhat=−ihat per Dave) — the XY "separation" was the
  beam-center fan reconstruction → XY panels now render TRUE positions
  from **`Telescope.ray_bundle`** (Dave's DRAW rethink, NO engine
  change).  Slice physics: the fold preserves y → **pupil-Y slice
  shows bench beam width** (pupil-X spread maps into z, invisible);
  slices = grid column within one pitch, evenly subsampled (nearest-N
  → gap); XY panels of folded designs start AT the fold.  **Figures
  saved in the example dir = usual practice.**  Suite 192/192.
  **ALL PUSHED (both repos sls-dev @ macos 9809af2 / MACOS_res
  4f422a9).**
- **2026-07-06 follow-up day (MACOS_res sls-dev `be03e49`, LOCAL —
  suite 194/194): next-time items 1+2+5 CLOSED, 3 half-answered.**
  1. ~~Broken saved .in~~ **FIXED**: `add_pupil` seeded BOTH inserted
     Returns with literal z-facing psi `[0 0 1]` and probed FEX about
     the unbiased +z axis; on the fold-in-feed bench the FP_return
     flat was ray-parallel → chief died.  Seeds now derive from the
     chief line prev→FP (identical to legacy on axial trains); probe
     offsets from the biased chief.  CLI `ray 1` traces the full
     7-elt train, retro legs equal to 7e-9.  **Vintage puzzle**: the
     earlier verify that traced OK was the RETURN-LEG-fold topology
     (axial beam at FP) — z-assumption only fatal for fold-in-feed.
     **FEX EP placement was RIGHT all along**: independent two-field
     chief-crossing = 0.616 m from image (0.5 mm agreement); the
     "paraxial 1.5–1.7 m" estimate was wrong.
  2. ~~Defocused spot~~ **REAL and FIXED**: FP station was 0.227 mm
     off best focus; removed by [2d] below.
  5. ~~FP tilt~~ **`Telescope.align_focal_plane`** (Dave: grid of
     foci, 2×2 prelim → 5×5/7×7 final): maps best-focus points (3-D
     closed-form LSQ point per field, no scan) over an N×N field
     grid + center anchor, fits the detector plane, sets FP Vpt+psi,
     returns tilt/defocus/sag map (`.map` ready-to-plot).  Guarded
     to run BEFORE add_pupil.  **Folded TMA: FP tilt 7.32° wrt
     chief, field-curvature sag ±1.1 µm over ±0.25′ (≪ f/20 DoF) —
     and the −tilt FOV ladder extends DL @2.3 µm from the ~0.5′ core
     to the FULL 2′-dia mapped field (0.25′ 0.016 / 0.5′ 0.018 /
     1.0′ 0.027 waves).**  Demo [2d] stage + fpmap figure; 2 tests.
  - Item 3 CLOSED by item 1 (Dave: "ORS is resolved by Item 1").
    Item 4 (metric) deferred; PSF metrics now available (below).
- **2026-07-06 pm (MACOS_res sls-dev `11ff118` + `7f54151`, LOCAL,
  suite 196/196):**
  - **`design_report`** (design/src) — Dave's one-page report: mirror
    list, first-order (EFL/f/#-at-M1/f/#-at-FP/plate scale/λD — **EFL
    measured LIVE from chief displacement per field angle**;
    spec.derived.EFL is the Seidel seed and was 4× off on the folded
    j18: the M3 pushback changed the relay — the folded design is
    actually EFL 34.4 m f/5.2), WFE ladder **with Strehl**, FP tilt,
    EP handle, shroud/clearance/AOI.  **Strehl = exact coherent sum
    over the de-tilted EP-referenced OPD** — NOT the INT pixel peak:
    an off-axis PSF walks across the FarField window and the pixel
    ratio measures sampling (1′ ring 0.98 exact vs 0.28 INT).
    **add_pupil now emits `PropType=FarField` on the EP ONLY** (Dave:
    everything else Geometric) → INT at the FP = the PSF (metric hook).
    FP-vs-own-retrace-legs excluded from the clearance verdict.
    tma_centered README (step-by-step adaptation guide) + demo [6]
    report stage.  Folded TMA: Strehl 0.954/0.860@1′/0.106@2.5′.
  - **`tma_unobscured` example** (Dave: visible 500 nm coronagraph +
    imager + spectrometer front end, slower M1, M2 close to the
    source–M1 beam): finder walks the slower-M1 ladder at **CONSTANT
    feed f/# (f1·m2 = 10)** — holding m2 fixed drives the feed toward
    the system f/# as M1 slows, the M3 relay degenerates to 1:1 and
    the decenter blows past 2.5·D (first failure mode).  **Design
    point f/2.5 (Dave)**: decenter 1.66·D, AOI spread 13° < 15,
    all-clear, shroud 3.64×D.  Demo: 7.66 → 0.011 waves ([2b] off-axis
    refigure!) → 0.030 field-balanced → freeform → align (FP tilt
    4.77°, defocus 0.605 mm) → **EFL 132.1 / f/20.01, Strehl 0.954
    center / 0.852 @0.5′ @500 nm, UNOBSCURED**.  Robustness: all-lost
    sentinel (9.9999e36) now NaNs wfe_field_diag/Strehl rows instead
    of crashing (a ring outside the realize_apertures envelope).
  - **VERDICT (Dave 2026-07-06): this design approach — the conic
    eccentric-pupil section — will NOT meet our requirements.  Closed
    out; the example stays as the recorded trade study (constant-feed
    ladder, AOI-vs-shroud price).  DIRECTION: return to the
    SPHERE+ZERNIKE approach for 3+n mirrors** (the sz_tma lineage:
    all-sphere base + Zernike departures, real intermediate focus;
    plus n relay/imaging mirrors — the 3+1 progression), targeting
    the visible coronagraph + imager + spectrometer front end.
- **2026-07-06 eve — the pivot slice (MACOS_res `afce4a0`, LOCAL,
  suite 196/196): `design/examples/freeform_unobscured/`** (name =
  Dave's) — sphere+Zernike at 500 nm on sz_tma's e5mono tilted-fold
  geometry, full toolkit (staged S0/S1/S2 → align_focal_plane →
  add_pupil FarField → **standalone reload verification as a stage,
  1305/1305 VERIFIED** → design_report).
  - **S0 centers DL at the visible bar (0.030 waves)** — the strategy
    holds at 500 nm; the ±1′ S2 field solve trades the center away
    (center 0.16 −tilt / worst 0.54 waves; Strehl 0.38 center).
    **NEXT design conversation: mode depth / staging weights / field
    size / the +1 mirror.**  Packaging: UNOBSCURED, **shroud 1.86×D
    (half the eccentric section's 3.6×D)**, AOI spreads
    8.9/10.2/0.9° all under 15°, EP 1.02 m from the image, FP tilt
    1.95°, 15 mm defocus removed by align.
  - **Dave's staging rule: NO apertures in the first design steps** —
    add them when the design approaches objectives.
  - **REAL BUG exposed (pre-existing, queued): `realize_apertures`
    measures footprint centers in GLOBAL XY (draw_rays) but emits
    them as LOCAL ApVec offsets** — correct only while the element
    origin sits at the global origin (coaxial/eccentric-section
    parents); saved TILTED-FOLD designs lose EVERY ray on reload.
    **`sz_tma.in` carries this latent** (verified: 0/1305 pass on
    standalone reload), tma_offaxis.in likely too.  Stopgap
    `Telescope.clear_realized_apertures` documents it; proper fix =
    ray_bundle 3-D footprints projected into the engine's aperture
    frame, when apertures re-enter the flow.
- Weak-POWER fold option = requested follow-on (seam noted in add_fold).

**Previously (2026-07-04): TWO thrusts, ALL COMMITS LOCAL (Dave has
NOT said push; four commits pending his word):**
- **macos sls-dev `662e86e`** — SAVE audit element-data bucket
  (PLAN §0 item 2 follow-through).  All gates green (mmacos 30/0,
  GMI 6/6, e5hex1 legacy SAVE byte-identical, fixture
  `ZGD_test_files/tst_save_keys.in` round-trips byte-identical under
  ifx AND gfortran AND ifx==gfortran).
- **MACOS_resources sls-dev `bcccea6` + `f817f61` + `04295fa`** —
  the design-layer day (details below).

### SAVE element-data bucket (macos `662e86e`) — as landed
All emission gated on value-set-in-memory:
- All 18 element-data keys round-trip: Coating (thickness un-scaled
  ×IndRef/Wavelen), Grad*, Doe*+OrderHOE (new DoeTrGrating CASE),
  Ampl*, LensArrayIndRef (covers XYIndRefFile), ArrIndRef +
  ArrWaveLen (λ2..λn), SegAp*, ZernCenter/XDir/YDir/Rad,
  ZernAnnularRatio (**writer misspelled `ZernAnnualRatio` since
  birth — ratio silently lost every reload**; + FF/Mon emission),
  ZCOZernType, GridSrfOrder, lData, nMetPos/tMetElt, EdgeSensors.
- Structural: pData..zData widened to any nGridMat>0 (SrfType 9/11);
  nGridMat/GridFile/GridSrfdx emitted for grids on NON-grid SrfTypes
  (iris_dp_ZGD Conic NSReflector segments lost their whole grid).
- **lensarr_indexes.inc heap stomp fixed**: 107×107 rec table vs
  mLenslet=250 → LensArrayIndRef/XYIndRefFile parse corrupted memory;
  table now sized from mLenslet + clamp + parse-site count guard.
- FmtD: '-0.0E+00' → '0.0E+00' (cross-compiler byte-identity).
- REMAINING for Dave: Opt*/CALIB-family + trace-state singles SAVE
  policy; Lou-UpdateNotes cross-check.

### Design-layer day (MACOS_resources, 3 commits) — as landed
Dave's product frame: **utilities + adaptable examples**; progression
2-mirror → 3-mirror → 3+1 → N-mirror, then instruments.  Layout:
`design/src/` utilities, `design/examples/` examples.
- `bcccea6` restructure: all 8 driver dirs → `design/examples/`
  (git mv), progression README.
- `f817f61` both TMA families: `tma_conic_recipe` + `wfe_field_diag`
  (design/src), `field_ring` + `Telescope.trace_at_field` (package),
  `examples/tma_centered` — the **j18 A/B answer**: same parent, 5′
  circular field @1µm → CENTERED 0.065λ (DL, symmetry intact) vs
  SECTION 0.098λ (excess = induced binodal field astig; freeform
  M2+M3 can't fix field-varying orientation).  Field ladder: conic
  wall 3′ square / 4′ circular DL; 5′ circular 0.095λ.
- `04295fa` **the 3+1 coronagraph front end**, three files per Dave:
  `tma_3plus1.m` (j18 DEMO: 0.84D compact section, 30″ patch 0.034λ,
  pupil relayed ~10× flatter −0.13/0.20mm, 5/5 clear; AOI spread
  21/24° on M1/M2 documented OUT of the 15° preference),
  `tma_3plus1_aoi_search.m` (constraint finder: steps PM–SM via
  tma_layout, verifies the FULL 4-mirror chain — AOI spread +
  clearance + shroud; **f/2.0, t1=11.8m=1.7×j18 MEETS 15°, but
  decenter grows 0.71→1.50D → shroud 1.6→3.2×D — AOI-safe eccentric
  sections are SHROUD-EXPENSIVE**), `tma_3plus1_optimize.m` (polsafe
  deliverable: 0.061λ patch, 5/5 clear, AOI 14.2°).
- Builder capabilities (each tested; tDesignTelescope 48/48):
  **psi parity ≥4th mirror** (mirrors 1-3 keep legacy all-(0,0,−1);
  a 4th mirror concave-to-a-−z-beam needs +1 — it traced CONVEX
  before, diverging the relay), `add_mirror(...,'conic',K)` explicit
  seeds (seidel can't seed a relay-past-focus chain),
  `optimize(...,'elts',SET)` subset (image vs field-mirror DOF split).
- Standing design rules recorded (Dave, all 2026-07-04): product =
  utilities+examples; packaging = fit a CYLINDRICAL LAUNCH SHROUD
  (keep M2 near the incoming beam as PM–SM grows); coronagraph
  polarization = per-mirror **AOI SPREAD across the beam** < 15°
  (not chief-ray absolute); the wide field (HWO 10×20′) needs a 4th
  POWERED imaging mirror, not more freeform on three; relays serve
  small per-instrument patches, the shared field lives at the TMA
  focus.

### Follow-on candidates (next session picks with Dave)
- **Tilted-fold 3+1**: Bauer folds hug the incoming beam → the
  shroud-cheap alternative to the eccentric section (the fold path
  `resolve_nmirror_fold_` already exists); the natural resolution of
  the packaging↔AOI trade the finder quantified.
- Pupil fine-tune: K4 / M4-conjugate scan against `pupil_quality`
  (current polsafe pupil −0.29mm defocus / 0.37mm astig).
- **Refit the 2-mirror examples onto the utilities structure**
  (Dave's note, recorded in design/README roadmap).
- Sprint 2D: segment the 3+1's M1 (SegMirMaker orchestration) — the
  reference parent now exists.
- N-mirror stage / HWO 10×20′ shared field (4th powered imaging
  mirror); instrument-building sequence.
- SAVE follow-ons: Opt*/trace-state SAVE policy (Dave); `Glass=` BK7
  fixture (zero coverage); mmacos opd-veneer clean error on the
  9.9999e36 sentinel; FEX footprint-autoswitch arm unexercised /
  journals unregressed; eac2 + obscured-rays parked at
  PLAN_DESIGN_LAYER Sprint 5 head.  Dave: NO opt-dev cherry-picks
  for now (incl. the older 60f886d).

### Previous state (2026-07-03, all pushed)
PLAN §0 items 1–5 ALL LANDED + PUSHED: #1 FEX `1139f78`, #2
SAVE/ApStop `96696fa`, #3 comments `d2da3cc`, #5 DXCALC `65a9c7f`,
#4 glass `65b2a67` (all macos sls-dev).

### Items 5 + 4 as landed (2026-07-03 pm)
- **Item 4 (glass tables)**: `tools/gen_glass_builtin.py` (new) →
  generated `macos_f90/glass_builtin.f90` (`LoadBuiltinGlass`, 200
  glasses from the canonical macos_glass_list.txt; verified
  programmatically identical, all 200×6).  All THREE loader
  fallbacks (macos_glass / smacos_glass / rl_macos_glass .inc) call
  it instead of hardcoded Air/BK7/LAK9; a found file still
  overrides.  CMakeLists: glass_builtin.f90 in COMMON_SOURCES.
  Tested both branches (file present → unchanged; file hidden →
  "using built-in catalog: 200 glasses").  NOTE: no Rx anywhere uses
  `Glass=` — zero test coverage for glass-by-name refraction; a
  BK7-refractor fixture is a worthwhile follow-on.
- Pending gate: final combined full mmacos suite RUNNING (GMI
  already all-pass; ifx+gfortran builds + mex rebuilt with both
  fixes).  On green: ready to commit as two commits (#5 utilsub
  DXCALC guard; #4 glass bake-in + codegen + CMake) — awaiting Dave.
- **Root cause (gdb, ifx build): `DXCALC` (utilsub.F) — when EVERY
  pupil slice has zero valid rays, the slice-search loop never
  assigns `iimax` (only set when `l>lmax`, and all l=0), and the
  second pass indexes `LRayOK` from `j=(iimax-1)*npts+1` =
  uninitialized garbage.**  ifx stack garbage → wild index → SIGSEGV
  (`forrtl severe(174)`, Dave's repro); gfortran's garbage happened
  benign (CLI survived) — why it looked ifx-only.  The earlier
  `dxRef=1d0` l==0 guard was necessary but not sufficient.
- Fix: `lmax=0; iimax=1` init + skip the second slice pass entirely
  when `lmax==0` (l stays 0 → existing warning + degenerate-beam
  guards take over).  Repro Rx: `scratchpad/optiix_badstop.in`
  (optiix_dp_conic with ApStop moved to 1e6 → all rays obscured at
  elt 1); driver `scratchpad/drive_macos_opd.py` (pty+gdb, answers
  the giza device prompt with /null).
- Verified: ifx CLI `OPD 55` → "All rays are obscured or lost" +
  9.9999e36 sentinel + prompt returns (was SIGSEGV); mmacos opd on
  the same Rx leaves MATLAB alive (wrapper returns non-struct on the
  sentinel — optional mmacos-veneer polish, not engine).  GMI
  all-pass.  Full mmacos suite RUNNING.  Both builds + mex rebuilt.

### Item 3 as landed (d2da3cc):
  - Capture: `GET_EQ` whole-line-'%' branch → `RxCommentCapture`
    (iosub.inc; banner-furniture filter for idempotency); storage in
    `elt_mod` (`RxCommentLine/Anchor`, cap 500, `RxCommentCtx`,
    `LRxCommentCapture` gate — MOD dialog also calls GET_EQ and must
    not capture; disarmed at both hosts' MBFile6 exits).
  - Anchor: `RxCommentCtx` = 0 in header, iElt per element block
    (set in msmacosio.inc's element DO loop).
  - Emit: `RxCommentEmit(0)` at `PrtSourceInfo` end, `(iElt)` at
    `PrtSingleEltInfo` end — comment before `iElt= k+1` is read in
    block k → end-of-block-k emission reproduces position (verified:
    "%Verify fElt, eElt for the following…" lands right before its
    ellipsoid element).
  - Rode along: `FmtD` 16→**17 sig digits** (bit-exact real reload);
    skip `BaseUnits/WaveUnits= (default-unit)` placeholder (GET_EQ
    eats parens on reload → came back empty; idempotency breaker).
- Verified: eac2 (14 comments) all survive at correct positions;
  pass2==pass3 byte-identical EXCEPT pre-existing 2-ulp `ChfRayPos`
  re-aim oscillation (ApStop-driven recompute each load — NOT a SAVE
  defect, out of scope).  ifx/gfortran CLI SAVE byte-identical.  GMI
  all-pass.  Full mmacos suite RUNNING (commit+push on green — Dave
  gave advance permission).
- **NEXT: #5 (opd all-rays-lost SIGSEGV), #4 (glass tables)**; work
  SAVE_KEYWORD_AUDIT element-data bucket (Dave review); eac2/
  obscured-rays parked at PLAN_DESIGN_LAYER Sprint 5.

### Last landed (2026-07-03) — FEX EP-radius rework + guards (§0 item 1)
- **FEX EP radius = chief-ray distance EP→(iElt+1) plane, ALWAYS** —
  whatever iElt+1 is (FP, coronagraph mask, ZWFS…): the far-field
  propagation distance for physical optics (Dave's spec).  Sign per EP
  EltID (Return −cr1dir / Reference +cr1dir).  Legacy `zp_iEm1` =
  fallback + autoswitch alternative.  Guards (all noisy): telecentric
  (parallel probe chief rays → FindCrossPt 0/0; keep station, radius =
  station→iElt+1, FLAT 1d22 fallback); beam-footprint sanity
  autoswitch; Rx-order flag (Return before EP return should usually be
  a Reference — fires on most legacy Rx, deliberate nudge).
  `tracesub.F` FEX block + decls; SXP/XPS untouched (SXP redundant now).
- **Compat pass green:** conforming double-pass Rx → legs equal by
  construction (pre-EP Return AT the focus) → round-off no-op.
  `Cassegrain.in` −1.2295 → +6.7907 == CassWithExitPupil exactly
  (right answer where legacy was degenerate).  **`eac2_7seg` −289.7 →
  +52597.6 — material shift, Dave to review.**  GMI all-pass; mmacos
  30 classes / 0 verification failures (tProperCompareCassFF crash
  pre-existing); pymacos has no FEX tests.
- **fex_sweep tool** saved at `MACOS_resources/mmacos/tools/fex_sweep/`
  (per-Rx process isolation; corpus = old_Rx sandbox + manual
  examples).  Corpus notes: ape.in trips a legacy loader parse bug
  ("Bad real number", segfaults under ifx build — pre-existing);
  several old Rx fail load (SegDemo/FFSegDemoAll/Luneberg/eac1/…)
  or the sweep's heuristic stop (j18sa/j18sc/6MST_segV3/dmt6seg).
- **Not yet exercised:** the footprint-autoswitch arm (no corpus Rx
  has a degenerate DEFAULT leg).  Journals not regressed.
- **Parked → design layer (Dave 2026-07-03):** the eac2_7seg leg
  discrepancy (−289.7 vs +52597.6; coronagraph, mask obscures all
  rays) and the "should FEX include obscured rays?" question — see
  the deferral note at PLAN_DESIGN_LAYER §8 Sprint 5 head (FEX
  centroid mode's all-obscured path builds psip from uninitialized
  CentroidSpot — latent; vertex/radius/footprint already use
  obscured rays).  Do NOT chase eac2 further in the engine for now.
- Negative-L answer (Dave's Q): ordinary train rejects L<0 (GO TO 98
  → miss); near Reference/Obscuring/Return/LensArray ifLNsrf=TRUE
  allows it and ConSrf picks the root by |L²−mpr| PROXIMITY to the
  ray's vertex distance (surfsub.F:~148) — no flow-of-light sense.
  That heuristic is WHY a too-small reference sphere fails silently.

opt-dev cherry-pick of 60f886d still pending.

### Last landed (2026-07-02 pm) — NS grid-frame fix + ns_griddata decomposition
- **Grid figures on NSReflector segments collapsed to a center-pixel
  per-segment piston**: the NS-block Reflector call sites passed the grid
  frame `pData/xData/yData/zData` indexed by `iElt` (NS-group entry elt,
  frameless → zeros) instead of `imin` (elt actually hit).  Two-line fix,
  tracesub.F:3714 + propsub.F:983.  **macos sls-dev `60f886d`** (pushed;
  **opt-dev cherry-pick NOT yet done** — it has the same slip via 1b535a5).
- **ns_griddata example redone as a 5-variant decomposition** (Luis Q&A):
  Conic/Zernike/GridData/ZrnGrData+flat/ZrnGrData variants; superposition
  3e-8 rms (real 2nd-order slope-resampling cross-term, mostly double-pass);
  flat grid bit-identical to none.  GridSrfdx=1.1 (=280/255 — grid span is
  (nGridMat−1)·dx centered on pData; Luis's 0.2 covered only ±25 of the
  140-radius segment).  **MACOS_resources sls-dev `e6f4f17`** (pushed).
- mmacos regression: 30 classes, 0 failures; tProperCompareCassFF heap
  crash pre-existing/known.  Both release builds + mmacos mex rebuilt.
- Remaining NS gap: the `NSRefractor` ROUTINE (refractive NS) still
  null-frames GridData — iris (all-reflective) unaffected.

### Landed earlier 2026-07-02 — conforming Reference (PASSIVE)
- **`Element=Reference` now accepts `Surface=Zernike`/`Aspheric`** so a
  conforming Reference CARRIES a Zernike basis definition (segment shapes) for
  GS-basis dev, but has **no effect on the light** (RefSrf unchanged; coeffs
  stored, never injected — I first built it ACTIVE = WRONG, Dave caught it,
  reverted; verified with-ref==no-ref OPD 9e-12).  Engine = 3 files
  (`EltSurfCompat` gate + 2 shared-parser fixes: `ZernModes` single-vs-wrapped
  IOSTAT read; `SrfTypeName(EltID→SrfType)` warning mislabel).  **macos sls-dev
  `c9fa767`** (pushed).  Example `e5hex2_refzern` (passivity + `make_gs_basis`
  + `run_dwdgrid_multi` split + `verifyall`) + `segment_grid_basis`/
  `grid_channels` exclude refs/non-segments: **MACOS_resources sls-dev
  `52d688d`** (pushed).  Other-session fixes: `b0d044c`.  **Rx GOTCHA: segment
  grid frame `pData/xData/yData/zData` must = clocked `pMon/xMon/yMon/zMon` or
  pokes don't localize.**  [[project-conforming-reference]]
- **NOT run:** pymacos regression (needs ifx build_release + pymacos rebuild;
  mmacos suite was green for the change).  Dave's `macos` alias (ifx
  build_release) rebuilt to current/passive.

### Landed 2026-07-01 — pointers for the next session
- **dwdgrid segmented examples** (MACOS_resources sls-dev `c07ed65` + `8beb774`,
  pushed): `run_dwdgrid_multi_multisegbasis` (per-segment `segment_grid_basis`
  struct fed as `dw_dgrid_multi` `influence`; verified 72 chan / 0.0% inter-seg
  overlap / 92px spread) + `run_dwdgrid_multi_singlesegbasis` (single shared
  `gs_zernike_segment_basis`) + generic `run_dwdgrid_multi` + `run_dwd{x,z,surf}
  _multi`; both segmented library drivers (generic RX='' flavor) +
  `plot_dw_per_element` + `gs_zernike_segment_basis` committed.  Experimental
  SegDemo scratch dirs deleted; unique FreeForm/Zernike fixtures preserved in
  `~/dev/MACOS_sandbox/segdemo_fixtures/`.  [[project-gridmat-generator]]
- **GMI build fix** (macos opt-dev `796dc51` + sls-dev `2f4948e`, pushed): opt-dev
  `makeall` `-DBUILD_GMI=ON`→`OFF` (Makefile is sole GMI builder) + MATLAB_ROOT
  autodetect (both branches, was hardcoded R2025b) — fixes Scott's GMIG.F
  `fintrf.h` fail (opt-dev double-build).  [[project-build-fixes-issue56]]
- **iElt fix** cherry-picked to MACOS_resources opt-dev (`0f071f2`→`15603a8`).
- All four branch-refs (both repos × sls-dev/opt-dev) in sync with origin.

### Open loose ends (next candidates)
- **Design layer NEXT (planned 2026-07-02): Sprint 5 — simultaneous
  focal + pupil optimization** (`PLAN_DESIGN_LAYER.md` §8 Sprint 5).
  Six slices, all MATLAB over existing XPS/`pupil_quality`/
  `check_clipping`: distortion metric → `pupil_quality_multi` →
  pupil-station clearance → stacked-residual objective wiring →
  3+1 field-mirror builder support → `design/tma_3plus1/` worked
  example (null sz_tma's +1.67 mm pupil defocus / 1.77 mm astig while
  holding WFE diffraction-limited).
- ~~Deferred engine: `Surface=Zernike` for `Element=Reference`~~ DONE
  2026-07-02 (passive), see above.  Remaining:
  Rx-collapse (modes×coefs in the engine); more Zernike types.
- (Pre-existing) `tma_onaxis` designer + the convex-secondary REVERT question
  ([[project-design-drivers]]); `test7.in` / `zmode_end` rename
  ([[project-engine-fixes-lega-shipped]]); layout realizability
  ([[project-layout-realizability]], local unpushed commits).

## In-session state NOT yet committed
**2026-07-07: freeform_unobscured "+1" — add M4 (sphere+Zernike
relay/field mirror) for a wider DL FOV + an accessible flat pupil
(Dave: "add mirrors … then work layout details. Then add
instruments").**
- Geometry PROVEN (probe A): 4-mirror folded chain builds/traces via
  `add_mirror` ×4 + `'derive'` (seidel any-convex signed-curvature
  path is N-general); M4 1.0 m past the TMA focus, R4=2·m·d4/(1+m)
  (m=1.6 → R4=1.2308, final ~f/20-34); `add_pupil` relays the EP to
  ~0.80 m past M4 on the 1.6 m M4→FP leg (mid-leg, accessible,
  ~36 mm — DM-scale).  tilt4=−12° keeps AOI ~12° and the back end
  compact behind M1.
- Builder edit (LOCAL): `optimize_freeform` `'lmon'` accepts a
  VECTOR (one per ELTS entry; scalar broadcasts; stable-unique
  pairing) — a near-focus M4 needs lMon ~0.09 m vs the 4 m body
  default (100× — fatal degeneracy).  +1 test
  `test_optimize_freeform_lmon_vector_validates` (tDesignTelescope).
- **Staging laws learned (probes A-D; engine-crash class confirmed):**
  (1) joint multi-mirror Zernike solves at a SINGLE field are
  DEGENERATE (surfaces near-redundant at one field) — CALIB lmlsq
  "DOFs correlated" then SIGSEGV on the diverged coefficients.
  (2) **center-only solves cannot pin the physical focus**: the OPD
  reference sphere absorbs pure defocus, so with modes 3/4 (BornWolf
  tilt/defocus) in the set CALIB parks the REAL focus anywhere (331 m
  measured) while reporting 0.03-wave "DL"; align_focal_plane then
  shows tilt≈90° + metres of "defocus" (foci strung along the axis).
  (3) dropping modes 3/4 does NOT fix it: the all-sphere start
  focuses at 20.5 m (not the paraxial 3.0 m!) — **e5mono's mm-scale
  Z3/Z4 terms are LOAD-BEARING power that pulls the focus to the FP;
  only MULTI-FIELD stages pin it honestly** (shipped 3M S2 → −15 mm).
  (4) M1-M3 modes cannot null the center THROUGH an uncorrected M4
  (stalls 1.1-1.3 waves — M4's term maps to the pupil distorted,
  beyond the mode set).
  → Correct flow: **full 3M staged solve (S0/S1/S2) → align (honest
  focus) → save_spec → fresh 4M build with M3 spacing = measured
  focus + d4 → carry M1-M3 freeform (capture spec.elt(k).freeform,
  re-apply set_freeform) → M4-alone FIELD solve (multi-field pins
  M4's power) → joint 4-mirror FIELD refine.**  Fix M4's lMon ONCE
  (field zone: 1.05×(walk+foot)) before any M4 solve — changing lMon
  rescales coefficient meaning.
- **ROOT CAUSE FOUND (probes E2–E7, 2026-07-07): the SHIPPED 3M
  solve is pathological — metre-scale canceling Zernike coefficients
  on M3 (mode 3 = +1.171 m! modes 10/11/21/22 = ±3.2–5.9 m) from the
  lMon ill-conditioning trap (M3's 0.12 m footprint on the 4-m body
  normalization = 3% of the disk), PLUS BornWolf mode 3 (tilt) is a
  pure GAUGE for CALIB's chief-referenced OPD metric (the chief
  follows the wedge — tilt is never pinned at ANY field count).**
  The 3M chain is self-consistent (OPD/align follow the wedged
  chief) but the real beam axis leaves the geometric chief: the M4
  append (placed on the resolve chief) saw the beam 0.6 m off-vertex
  → wrong shell conjugates → "relay image at 0.946 m vs paraxial
  1.6" → every M4 solve diverged (SFFZP bracketing, all rays lost).
  Chain of evidence: e6 (feed identical 3M vs 4M, .in diff = only
  nElt+zElt), python independent sphere trace (image at 1.5998 ✓
  textbook), e7 dir-free two-point waist + |hit−V4| = 0.55–0.64 m,
  coefficient dump (the smoking gun).  **Fix (probe E8, iterated): per-mirror lMon = the FIELD-ZONE
  radius (footprint + field walk), fixed ONCE before S0 — center-
  footprint lMon ([4.2 0.70 0.10]) fixes S0 (0.34 waves, sane
  coefficients) but stalls S1 at 23 waves: fields walk ~58 mm/′ on
  M3 (chief pivots at the M1 stop), outside a center-sized disk →
  [4.2, 0.70, 0.20] for ±1′.  Tilt mode 3 STAYS — dropping it
  stalls S0 at ~10 waves (it is the merit's tilt-removal channel);
  with a well-conditioned basis LM grows no wedges, and the flow
  VERIFIES that (solved max|coef| print + the focus' lateral offset
  from the geometric chief).  Mode 4 = real power (the all-sphere
  start does NOT focus at the paraxial station), pinned by the
  multi-field stages.  lMon changes between solves are FORBIDDEN
  (coefficients are tied to their normalization radius — the
  read-back now stores lmon and re-solves inherit it).**  Also fixed en route: optimize_freeform read-back now
  stores lmon with the solved figure + inherits stored lmon on
  re-solves (was silently re-emitting solved coefficients on the
  body radius = a different surface).
- 3M reference to beat: center 0.29 raw/0.164 −tilt, 1′ ring
  0.56/0.457 waves, Strehl 0.376 (shipped report) — likely partly
  ARTIFACTS of the pathology; e8 re-solves the 3M under the
  corrected doctrine before appending M4.
- **e8/e9 verdicts:** honest 3M closes S0 0.46 / S1(±0.5′) 3.95
  waves with mm-scale coefficients (M3 max 1.3 cm ✓) and small
  chief wedge (lat 41 mm); a 3M S2 at ±1′ dives into a bad basin
  (21/58 waves — the outer ring is beyond the honest 3-mirror
  basis; the pathological 0.54 was bought with the insane terms).
  M4 append on the honest state is geometrically CLEAN at last
  (beam 7.7 mm off-vertex; exit waist 1.563 vs paraxial 1.600) —
  BUT the M4-ALONE ±1′ field solve stalls (130/346): in the honest
  state the exit pupil hugs the image (chief pivots ~0.1 m out →
  M4 field walk 258 mm/′ vs 9 mm footprint = DISJOINT per-field
  patches; 11 fields × local DOFs ≫ 15 modes on the union zone =
  underdetermined in M4's favor... against).  **Architecture fix
  (e10): solve the 3+1 as ONE system — S0 center M1-M3 (through
  the M4 sphere) → S1/S2 JOINT all-four FIELD solves.**
- **e10 verdict (the third architecture, same wall): S0 1.41 /
  S1 joint 7.2 / S2 joint ±1′ = 23.5 center / 36.9 worst — all
  coefficients mm-scale (sane).  Across append/M4-alone,
  append/joint, and single-system joint: the 15-mode e5mono basis
  does NOT reach the ±1′ ring at 500 nm on this geometry.  Levers:
  mode depth (BornWolf extends only to ~23 via indices 26-33 —
  ZerntoMon3's permutation table ends there; ANSI reaches 45),
  pupil-aware re-layout (Sprint 5), or field/λ spec.**
- **Optimizer thread (Dave's question 2026-07-07): CALIB's FD-LM
  isn't the root problem but is the wrong engine for figure solves
  — no bounds (trial steps SIGSEGV the engine), no degeneracy
  visibility ("DOFs correlated" aborts), scalar merit with gauge
  modes, 200 iterations × FD-Jacobian cost.  Built
  `design/src/zern_jacobian_solve.m` (Dave's 2026-06-24 idea):
  poke-per-mode OPD Jacobian (engine setters, no re-emit per poke)
  → per-field piston/tilt PROJECTED out of residual+columns (kills
  the tilt gauge by construction; merit = tilt-removed WFE = the
  blur metric) → truncated-SVD minimum-norm solve (spectrum
  printed — degeneracy visible; canceling-coefficient null-space
  solutions excluded by construction) → apply, re-emit, iterate
  2×.  ~100× cheaper than a CALIB run.  probe_jac1 (running) =
  first live test on e10's S2 state.**
- Diagnostic toolkit built (scratchpad probes): dir-free two-point
  waist (LSQ of ray lines), feed/exit crossing checks, coefficient
  dumps, independent python sphere trace.  `get_ray_info` after
  `trace(k)`: `.dir` = OUTGOING at elt k (doc says incoming — fooled
  E5(a)).
- NEXT after probe: worked example (same dir, staged script per the
  finder/optimize pattern) + README + fast suite; then layout
  details (clearance/AOI/shroud/fold stations); then instruments.

## Just tried / ruled out (with why)
- lMon policy alone as the S0 stall cause — RULED OUT (probe B:
  default lMon also stalls at 1.10 waves, but footprints collapse
  clean, so the relay seed geometry is right).
- Joint M1-M4 center solve — RULED OUT (degenerate → lmlsq fail →
  engine SIGSEGV; probe B).

## Next concrete step
Commit the shippable set (Telescope.m lmon fixes + 4 tests +
design/src/zern_jacobian_solve.m) on fast-suite green; then the
design fork is DAVE'S CALL: (a) mode-depth test — one ANSI-45
jacobian assembly + column-space projection answers "does depth
reach ±1′?" in ~15 min; (b) pupil-aware re-layout (Sprint 5
simultaneous focal+pupil); (c) re-scope field/λ.  The 3M example
also needs re-solving under the doctrine (its shipped numbers rest
on pathological surfaces) — Dave to confirm before his example's
recorded results change.  The drafted 4M example is parked in the
scratchpad (freeform_unobscured_4m_DRAFT.m) until the fork closes.

## Open micro-questions (slice-local)
—

## Promote-on-land  →  then CLEAR this file
> Same commit as the `design-sprint-N` tag: move each item to its
> permanent home, then reset this file to the empty template below.
- [ ] PLAN checkbox(es) ticked: <which>
- [ ] `CORE COMPLETE <date>` blockquote added to the sprint: <one-line summary>
- [ ] §10 Decisions: <new "Made" entry / resolved "Open" item>, dated
- [ ] CLAUDE.md / nested gotcha captured (if a new trap was found): <where>
- [ ] agent `MEMORY.md` learning (if workflow-level): <one line>
- [ ] worked-example committed + named: <example_*.m>
- [ ] **reset CURRENT_SLICE.md to empty template**

---

## Empty template (reset state — copy over the above on land)

```
## Active slice
- Sprint / item: —
- Plan anchor: —
- Branch / worktree: sls-dev + sls-dev @ —
- Definition of done (honest): —

## In-session state NOT yet committed
—

## Just tried / ruled out (with why)
—

## Next concrete step
—

## Open micro-questions (slice-local)
—

## Promote-on-land → then CLEAR this file
- [ ] PLAN checkbox(es)
- [ ] CORE COMPLETE blockquote
- [ ] §10 Decisions entry
- [ ] CLAUDE.md / nested gotcha (if any)
- [ ] agent MEMORY.md learning (if any)
- [ ] worked-example committed
- [ ] reset this file
```
