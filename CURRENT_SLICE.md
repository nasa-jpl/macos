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

> **CURRENT STATE (2026-07-21) — READ THIS FIRST.**  The active thread
> is the **end-to-end worked-example series** `MACOS_resources/mmacos/
> design/examples/e2e/` (s1 design → s2 relay → s3 segmentation →
> s4 sensitivities → s5 MET → s6 compare → **s7 closed-loop sim**).
> Everything through **s7 is SHIPPED + COMMITTED** on `sls-dev`
> (MACOS_resources `224891b`, macos `e3905ec`; both PUSHED 2026-07-21).
> The s7 design + physics + numbers are recorded in the "2026-07-20 /
> -21 (s7 SIMULATOR SESSION)" block further down this file.
>
> **NEXT = s7b:** upgrade the RBCS pose estimator from the static
> weighted-LS/BLUE form to the **steady-state Kalman filter** (Tesch
> *RBCS Algorithms* §2.3.3 eq 12-14, predict/update with the Riccati
> gain) and add **figure states** to the measurement model via the
> `dmdz`/`dmdgrid` blocks (macos.design.dmet_dfig) so the loop can
> SENSE and CORRECT the figure floor itself — not just the periodic
> image-based WFC.  The OSE single-step static estimator is this with
> converged gains.  Background PDFs in `MACOS_sandbox/Documents/`:
> `Tesch_RBCS_algorithms.pdf` (read), `OSE_Eqns_2019.pdf` (read),
> `2025_JATIS_HWO_Special_Issue-2.pdf` (PENDING — read before s7b).
>
> **Resume protocol:** read root+nested `CLAUDE.md`, then memories
> `[[project_recast_runners]]` (runners + RBCS loop + s7),
> `[[project_e2e_example]]` (s1–s6 heritage),
> `[[feedback_demo_plot_conventions]]` (movie plot rules), then the s7
> block below.  Runner = `design/runners/run_simulator.m`; driver =
> `examples/e2e/s7_simulate.m`; test =
> `tRunCompare/test_run_simulator_time_history` (SUITE_FAST).
> Everything OLDER than the s7 block (below, from 2026-07-16 on) is
> LANDED HISTORY — context, not in-flight work.

**2026-07-16: VISUALIZATION + MET-layout physicality/v2.** Session
record (chronological):
1. **Luis GridData transpose — DONE+PUSHED** (MACOS_res `17f1239`):
   `macos.read_grid_file`; engine GridInit reads file line=COLUMN;
   mmacos `elt_grid_add`=[x,y] vs pymacos=[y,x] OPPOSITE (memory
   `reference_grid_orientation_convention`).
2. **General visualizer — DONE, COMMITTED NOT PUSHED** (macos
   `032ddb8` + MACOS_res `d8011c6`): Dave: "work with any prescription
   — beam, optics, MET paths if present" (Lou's unfinished 3-D
   visualizer, modernized).  Engine: `Draw3DVec` 3-D DRAW capture
   (traceutil_mod; CTRACE fills at the 3 DrawRay sites via
   iDrawRay_global) + `draw_rays3d_get` + `met_geom_get` (endpoints ==
   met_get order; ride perturbations) + **OPD batch-hang fix: GO TO 6
   reprompt → abort — `trace(k)` on a Segment/NS elt spun FOREVER in
   SMACOS (same class as SPOT ed48d4f; sweep of remaining IACCEPT_S
   reprompt sites still OPEN)** + draw_rays_cmd vestigial stack args
   removed ("Unknown command" noise; load_rx still emits some —
   pre-existing).  mmacos: `view_rx` (DRAW-fan harvest, NOT per-elt
   trace(k)) / `draw_rays3d` / `met_geom`; design layer: `met_view`
   (3-D + face-on panels) / `hex_tile` / `seg_boundary` (hex + PIE
   wedges, arc-length `sample`); `segment_rx` exposes width/gap/grid.
   Examples: `examples/view_rx_demo` (Cass/Coro/e5mono+met) + e5_seg
   figures.  Tests: tMet 6/6, tMetView 4/4 (**findall can't reach
   sgtitle Text — title mirrored to fig.Name**), tSegmentRx 4/4,
   tReadGridFile 5/5.
3. **Dave's design-review constraints (figures drove these — the
   PURPOSE of the viz)**: boundary-true tiles (width/2 apothem, ONE
   global clocking per manual §segments — NOT per-seg face frames);
   launchers AT segment edges + edge_off 5 mm (add_met DEFAULT; edge
   placement alone: e5_seg edge+MET 5.98→4.34 nm); **MIN_SEP 50 mm
   between ANY two launchers (corner junctions!)**; **fiducials must
   MOUNT ON M2 ≤~25 mm inside its ApVec rim (615 mm here; r_fid≈590)**;
   **nf ≥3, likely 6**; **extra/M3 launcher ring hugs that element's
   physical radius + edge_off (was floating at r_fid)**; met_view hub
   disc now drawn at REAL hub radius (was fiducial-fit — masked the
   rim violation).  metopt v2: spread + CLUSTER-PAIR families,
   per-beam fiducial assignment enumeration (add_met `pair_map`),
   rim-zone RFID, MIN_SEP gate, hierarchical passes.  **DONE
   (3bdcda1+1103de0+369de7a)**: MIN_SEP killed the v1 corner ring
   (3 segments' launchers within 10 mm at every corner junction =
   INFEASIBLE); physical as-built prior 9577 → edge 3641 → MET 450 →
   edge+MET **232 nm** (MC 1.3%) — **aft/M3 ring at its physical
   100 mm = THE bottleneck**; optimizer (42-DOF SEGMENT sub-merit)
   → **3.358 nm**, winner spread [30 80 160]° + opposite-jump 6-fid
   map [1 4 2 5 3 6] at the rim (rfid=615), engine-FD 0.00%; the
   cluster family (64,800 layouts; kept FIRST-CLASS — deformation/
   rigid decoupling, Dave) lost to spread on rigid DOFs.
   PATTERN_FRAME 'radial'|'segment' knob (ring-uniform builder
   parts).  met_view readable views: face-on projected beams +
   color=segment association + **M2-M3 face-on inset** under the
   legend.
4. **e5_pie POLY-APERTURE ASSESSMENT DONE** (design/examples/e5_pie,
   6391b00; supersedes the sandbox copy): verdict YES — PolyApVec +
   explicit xObs emission loads clean; segments clip ZERO rays at
   nominal (fail_elt histogram); round-trip 6e-8 mm → 'rxpoly' reader
   trivial; pie seg_boundary/edge placement validated on the REAL
   fixture.  **Ray-loss question RESOLVED (Dave OPD review, MACOS_res
   f52ba4c): the RayFailElt=14 rays fall into the inter-segment GAPS —
   correct physical clipping, NOT an engine bug (elt-14 = cosmetic
   return-leg attribution).  Same commit: center segment (Elt 1) now
   emits a generated circumscribed 24-gon PolyApVec (Polygonal, no
   Circular special case); loss 456→442; all 7 polygons round-trip
   6e-8 mm.  rxpoly productization UNGATED.**
   Gotchas: macos.trace().nRays = SOURCE count (use ok_pass/fail_elt
   for parity); .presc segment blocks carry no ApType/nObs (append).
5. **Census for Dave's presentation**:
   `~/dev/MACOS_sandbox/improvements_census_draft.md` (~500 commits
   FreeForm→today, categorized, presentation arcs).  MEMORY.md index
   compacted 22→10 KB (no content loss).

**2026-07-16 (second session): the whole NEXT queue LANDED (all
local/unpushed on sls-dev — macos d09adc0+19d2c08, MACOS_resources
f52ba4c+10a955c+1dde42b+1dd4a64):**
1. **e5_pie manual example + aperture productization** (`10a955c` +
   `1dd4a64`): Dave's OPD review closed elt-14 (gap rays, correct
   clipping).  Center pie segment = HEXAGON from the traced footprint
   (poke-diff, NOT a disc/24-gon; corners at k·60°); ring-1 wedges
   abut it along straight CHORDS (flat (w−g)/2 + gap → chord (w+g)/2;
   obscuration = apex TRIANGLE, no inner arc — Dave's "circle around
   Seg 1").  `design.seg_apertures` (hex corners exact / pie hexagon +
   chorded sectors; xObs from psiElt NOT zMon) + `segment_rx
   emit_apertures/ap_pad/ap_obs` + `seg_boundary source=rxpoly`
   (auto when every segment declares PolyApVec; boundary = polygon
   minus obscuration via polyshape, LARGEST region — %.10E rounding
   leaves slivers, boundary #1 blind-take lost a tile).  Manual-grade
   runner `design/examples/e5_pie/e5_pie.m`, figure per step; pupil
   overlay needs affine centroid calibration (OPD grid transposed +
   mirrored; xGrid=(−1,0,0)).  met_view M2-M3 inset got axis labels.
   tSegmentRx 8 checks green.
2. **hub/aft DOFs in the analytic MET merit** (`1dde42b`):
   `design.met_bodies` = engine-truth frames (RptElt pivot + TElt
   rotation block via elt_rpt/elt_csys_get; segment triads reproduced
   EXACTLY).  metopt v3 merit = full 54 DOFs; fiducials RIDE the hub,
   aft ring rides the aft body.  Winner FLIPS to CLUSTER pairs at the
   corners (pmap [6 3 1 4 2 5], rim fiducials): rms 3.429 nm,
   engine-FD 0.00%; worst-mode 184 nm exposes the weakly-observed aft
   direction the 42-DOF merit hid.  tMet pins the 54-col identity.
3. **IACCEPT_S reprompt-loop sweep** (macos `19d2c08`): 21 sites in
   macos_cmd_loop.inc (BUILD/ORS/SRS/FEXIT/SXP/XPS/SPOT/MVAR/DVAR/
   GPERTURB/LPERTURB + AVAR's host-killing `stop`) now abort to the
   main loop; kept the 0=quit loops/MRESET/zernRange (self-
   terminating).  Both engines rebuilt, mmacos relinked, fast suite
   25 classes + tSegmentRx + tMet green; BUILD/SPOT at a Segment now
   abort in 0.01 s.

**2026-07-16/17 (third session): view_rx v2 — SOLID layout viewer
(Dave-planned + 6 review rounds; macos 30f4125+0c05406, MACOS_res
ac57bf4+b7c4c4b, all LOCAL):**
- Engine api: `ray_hist_set`/`ray_pos_hist_get` (traceutil RayPosHist
  — Lou's Vis3D substrate, now API-visible; slot 1 = source),
  `elt_info_get` (EltID/ApType/ApVec/xObs/lMon/PolyApVtx — IN-PLANE
  about VptElt, 5e-8 mm round-trip), `src_seg_get` (GridType/nSeg/
  width/gap).
- view_rx: rings-and-spokes FILLED bundle from the full traced grid
  (polylines connect the slots each ray REACHED — segmented ok is
  sparse), THIN flat two-tone sag-following shells (LightTools ref
  deck; no lighting), meridian profile curves both faces, sag sign
  calibrated vs crossings, consecutive-Refractor lens JOIN, EXACT hex
  Segment tiles (src_seg_get truth — no overlap, gaps read),
  source-plane ring (collimated-source cue), 'show'
  beam|beam+met|met, per-channel ray colors.
- macos.view_std (NEW): standard beam-aligned 4-panel figure (front
  from behind the SOURCE at the optic's face / back / iso / side),
  SOURCE AT LEFT, light → right, per-panel [az el] fine-tune;
  manual campos + camva('auto') + axis-off panels (hand-computed
  CameraViewAngle misframes; labels collide otherwise).
- Gotchas: ray_hist('on') must macos.modify() (grid-setter-retrace
  class); DRAW fan resamples its own rays (validate to grid pitch).
- Demo: 4 cases incl. e5hex1 via view_std; tViewRx 5/5 (SUITE_FREEFORM
  256), tMetView 4/4, tMet 6/6, quick 67/0.

**NEXT SESSION (Dave's objective, 2026-07-17): COMPLETE END-TO-END
worked example for users to hack, built from the parameterized
runners/utilities.  Dave 2026-07-17: EVERY stage runner must produce a
THOROUGH design report (saved text alongside the artifacts — identity,
first-order, field performance, and the stage's own metrics: Jacobian
condition/rank for stage 4, MET observability + post-control residual
for stage 5, time-history stats for stage 6) AS WELL AS graphics
(view_std/view_rx + the stage's metric figures).  Dave's SEQUENCE
(verbatim):
1. telescope design (TMA+FF feeding back end), with views and
   performance report;
2. add imaging instrument (3-4 mirrors to widen the field), new views
   and performance report;
3. segmentation, new views;
4. generate dwdx, dwdz, dwdgrid;
5. MET, MET-optimized performance report, dxdl, dxde, dwdl, dwde; new
   views with MET;
6. a SIMULATOR that generates PSFs using mmacos OR the linear model
   (user switch), driven by an x, z, grid TIME HISTORY.
Building blocks already shipped: Telescope/sz_tma builders +
design_report (1-2), segment_rx + emit_apertures (3), dw_d*_multi
harvests (4), add_met + metopt v3 + dmet_dx/dldx_analytic/met_bodies
+ edge_sensors (5; dxdl/dxde = the estimator gains from H, dwdl/dwde
= dwdx·gains), COMPOSE/psf + linear w = dwdx·x + dwdz·z + dwdgrid·g
(6).  view_std/view_rx figures at every stage.

Smaller queue: GMI mex not relinked against the new engine (next
makeall); pymacos export of the 4 new api routines deferred; rxpoly
for IMPORTED segmented Rx untested on an external fixture;
PLAN_DESIGN_LAYER promote-on-land.  Full record + gotchas: memory
`project_met_visualizer.md` + `project_e5pie_apertures.md`.**

**2026-07-17/18 (e2e sessions 2-3): STAGES 1+2 LANDED + VIS REDESIGN
(MACOS_res `e20d6b1` -> `5e94320` (s2 v1) -> `cbbb7c9` (VIS), all
PUSHED except cbbb7c9 LOCAL).  ARC: s1 v1 (f/1.25, m2=16, fold, bias
5') was NOT a VIS imager (Dave ran the .in) -- bias was the killer
(~bias^2).  VIS design point (Dave): f/1.75, m2=8, M3 EXTRACTION TILT
1.2 deg (return off the feed axis -> fold clears at bias 2').  s1 =
DL telescope (+-1' 0.0157 -tilt waves @500nm); s2 = 3-mirror bench
relay (M4 corrector/M5 collimator/M6 camera, radii DERIVED from the
collimator condition), +-2' at 0.26 -tilt, distortion 0.28" (M4
blur-guarded distortion stage stands down).  PROCEDURE (Dave: record
everything) = README 'design procedure' section, 11 rules each from a
failed run: WFE-based bias pick (small-bias conic basin K3->-3..-4 is
repeatable, continuation does NOT rescue), guarded M1 common-mode
null, SVD engine for degenerate bases (CALIB SIGSEGVs), joint-then-
null order, detector-plane re-fit last, distortion!=blur (M4 near
focus = the reflective distortion corrector; affine-projected chief
metric; blur-guarded).  DEAD ENDS on record: per-field Zernike patch
corrector (rank collapse), 4th near-pupil mirror (common-mode), full
BornWolf 3:25 (+13% only), relay tilts 6 deg (in/out beams overlap).
m2=12 A/B in ~/dev/MACOS_sandbox/e2e_m2_12 (bias 4', 0.25 waves).
OFFNER RELAY SHIPPED (MACOS_res 4ab5229): DL over the FULL +-2'
(Strehl 0.99->0.93, -tilt 0.015-0.043, pure spheres) via NEW
design/src/offner_layout.m (concentric chief solve -> Bauer chain;
spheres held, convex stop mirror out of solves, M4 = flat routing
fold).  Zigzag variant kept selectable (DL to 1.5', a7ed2f6).

**2026-07-18 (e2e session 4): FIELD SHIFT + HOLE CHAIN + s3 LANDED.**
(1) Field center -0.7' off the s1 bias (Dave read the s2 WFE map;
P.inst.field_dy_arcmin applied to set_field_bias -> artifact chief =
science center, s3-s6 inherit): worst +-2' 0.043 -> 0.0231 -tilt,
Strehl floor 0.965; [4f] center scan in runner+report; full -1.05'
re-solve = flatter interior / worse 2' edge, BOTH kept in
s2_variants/dy-0.70|dy-1.05 (Dave: keep both).  (2) M1-hole chain
(Dave: show the hole in all layout views): set_hole EMITS a real
ObsType=Circle obscuration (trace clips the 5 central rays); NEW
engine api elt_obs_get (macos_api_mod, codegen'd; pymacos export
deferred); view_rx renders circular obscurations (view_std inherits);
segment_rx carry_obs -> center segment.  (3) s3_segmentation.m: pie
(7, e2e_pie.in) + hex2 (19, e2e_hex2.in), physical apertures, 128-pt
grid (Dave: 41 too coarse), P.seg.variant="pie" feeds s4-s6; pie
12090/12520 pass (392 gap/rim), hex2 9838/9876 (hex tiling has no gap
rays); both VERIFIED.  TWO segment_rx product bugs fixed (tests
green): Surface=Zernike parent FIGURE dropped by SegMirMaker (15 um!)
-> carried into each segment's FF channel; SegMirMaker<->engine
tiling contract needs xGrid=(-1,0,0) AND SegXgrid=(-1,0,0) in the
merged Rx (design-layer (+1,0,0) left ray->segment 180 deg off the
frames; PSEG ignores SegXgrid = dead code, HSEG anchors to it).
**LATENT e5 FINDING for Dave: the e5 hex corpus has the SAME 180 deg
frames<->tiling offset (ray_hist-verified) — invisible bare, but
e5_seg's forward model pairs dwdx (tiling identity) with dedx/dldx
(frame identity) point-reflected; fix = SegXgrid -1 / regenerate,
touches committed Sprint-2D references — Dave's call.**
center_focal_plane now engine-truth (trace+get_ray_info; was
draw_rays plot-projection, U-sign follows grid handedness — why the
emitter KEEPS +1 and only segment_rx's output flips).  Fast suite
229/0; tSegmentRx 8/0 (+figure/obs-carry tests), tDesignTelescope
67/0 (+hole test).
**2026-07-18 (cont.): GENERAL FIX LANDED IN THE GENERATOR (Dave:
"fix it in e5 -- preferably in a general way").**  SegMirMaker.f, two
fixes (MACOS_res d54405f): (1) header SegXgrid now emits the in-plane
basis ACTUALLY USED — the back-facing-mirror 180-deg basis flip
(segment-numbering convention) negated xs/ys for frames+SegCoord but
the header kept the pre-negation vector → engine HSEG tiled rays
point-reflected from the frames (EVERY back-facing fixture, e5
included).  (2) Zernike-aware LoadParent: Surface=Zernike parent
figure merged into the FF channel (segments present the as-designed
surface; RptElt placement sees the figure sag).  Byte-identity refs
regenerated (delta = SegXgrid line only; Hx identical).  NEW
tSegmentRx gate test_state_consistency_rays_on_frames (ray_hist truth
vs frames, Pie+Hex, every off-center segment) = THE invariant: seg
k's DOFs move seg k's w/e/l only.  e5 corpus REGENERATED consistent:
e5_seg edge+MET 229.5 nm (was 232), metopt 3.421 nm (was 3.429),
engine-FD 0.00%; e5_pie unchanged behavior.  segment_rx text-carry
stands down when SMM already carried the figure.  Suites: tSegMirMaker
3/0, tSegmentRx 9/0, tMet 6/0, tEdgeSensors 3/0, tMetView 4/0, fast
230/0.

**s3 GRAPHICS + PIE-GEOMETRY FIXES (Dave review, 2026-07-18).**
Dave flagged the s3 views (EP dashed circles E15/E16, unclear pie
center segment) — pulling the thread found TWO real geometry bugs in
the pie aperture emission, not just rendering:
- **view_rx polish:** Return elements HIDDEN by default (`'returns',
  true` restores; tViewRx test); Segment plates draw with alternating
  face tints, center cell (Seg1, carries the hole) distinctly darker,
  and NO per-tile profile curves (7–19 overlapping meridian-curve
  sets were the spoke mush hiding the center segment).
- **Pie wedge apertures were apex SECTORS to r=0** (+ triangle
  obscuration faking the chord): every declared-polygon consumer
  (view_rx plates, step-3 figures) drew wedges COVERING the center
  hexagon with edges converging at the center.  A ring-1 wedge is
  CONVEX (disc ∩ 3 half-planes) → now emitted as the TRUE chorded
  polygon, no obscuration.  **Gaps are uniform-width slots** (side
  edges PARALLEL to the sector rays at ±g/2, not angular offsets that
  converge at the center — Dave).  New shared helper
  `macos.design.pie_wedge_geom` used by seg_apertures + seg_boundary
  + view_rx (ONE geometry source, no third desync).
- **pie_rings classification robustness:** e2e's Zernike-figured
  parent scatters ring-1 wedge radii ~5 µm (frame-tilt leak in the
  tiling-plane projection); the old 1e-6·max(rc) tolerance split the
  6 wedges into degenerate 1–2 member "rings" → 2π/nnz sector spans +
  go/sin(π) 1e14-vertex blowups in the emitted .in (the broken
  s3_footprints_pie Dave saw).  New `macos.design.pie_rings`
  (width-scaled tolerance, shared by all three sites) + jitter test.
  e5 never tripped it (exact symmetry) — the e2e example did its job.
- Regenerated: e5_pie (wedges 15 verts, 502 clips, no obscurations)
  + e2e s3 (pie 12106/12520, 376 gap/rim clips ≈ prior 392; hex2
  unchanged 9838/9876; both VERIFIED).  tSegmentRx 11/0 (+2 geometry
  tests), tViewRx 6/0 (+returns test).
- MET is NOT in s3 by design — metrology config + optimization is s5.

**s3 GRAPHICS REVIEW → PIE GEOMETRY FIXES + HELPER HOIST (Dave,
2026-07-18).**  Dave's s3-view flags exposed real emission bugs:
- **view_rx:** Return elements HIDDEN by default (E15/E16 EP dashed
  circles; `'returns',true` restores, tViewRx test); Segment plates
  alternate face tints (center cell darkest, carries the hole), NO
  per-tile profile curves (the spoke mush).
- **Pie wedges were apex sectors to r=0 + triangle obscuration** →
  declared-polygon consumers drew wedges covering the center hexagon.
  Ring-1 wedge is CONVEX → emitted as the TRUE chorded polygon, no
  obscuration; **gaps = uniform-width slots** (side edges PARALLEL to
  sector rays at ±g/2 — angular offsets converge at the center,
  Dave's rule).  Shared `macos.design.pie_wedge_geom` (seg_apertures
  + seg_boundary + view_rx = one geometry source).
- **`macos.design.pie_rings`:** ring classification needs a
  WIDTH-scaled tolerance — the e2e Zernike-figured parent scatters
  ring-1 wedge radii ~5 µm (frame-tilt leak in the tiling-plane
  projection); the old 1e-6·max(rc) split them into 1–2 member
  "rings" → 2π/nnz spans + go/sin(π) 1e14 vertices in the .in.  e5's
  exact symmetry never trips it — the second consumer found it.
- **Hoist (Dave: general-purpose runners):** `macos.design.
  seg_footprints` (poke-diff engine-truth footprint measurement) +
  `seg_footprint_view` (calibrated overlay figure) — e5_pie.m and
  s3_segmentation.m now thin narratives over them; duplicated
  pupil_axes_ copies deleted.  Regens BEHAVIOR-IDENTICAL pre/post.
- Regenerated: e5_pie (15-vert wedges, no obs, 502 clips) + e2e s3
  (pie 12106/12520, 376 gap/rim clips; hex2 9838/9876; both
  VERIFIED).  Suites: tSegmentRx 12/0 (3 new tests), tViewRx 6/0,
  tMet 6/0, tMetView 4/0, **fast 232/0**.
- **PLAN §0 model-transition heap crash REOPENED** (Dave asked; the
  full no-arg mmacos suite SIGSEGVs at tViewRx setup after ~26
  classes in one process — evidence + repro + suspects recorded in
  PLAN.md; `fast` subset is the green signal until refixed).
- MET is NOT in s3 by design — metrology config + optimization = s5.
- **s2 AUTO FIELD CENTER LANDED (Dave's ask):** sections [2]–[4e] now
  run in a 2-pass loop; [4g] maps ±3′ (13×13), finds the centroid of
  the raw-WFE<0.02 region, and if the chief is >0.15′ off it, pass 2
  re-solves there.  Result: chief moved another −0.71′ (bias 1.3′ →
  **0.59′ total**), worst ±2′ −tilt **0.0231 → 0.0215**, <0.02-region
  31→52/169 pts, centroid residual +0.03′; [4f] scan confirms the new
  center is the optimum (both ±0.35′ neighbors worse).  Clearance
  UNCHANGED by engine-truth ray-to-body margins (M2 −230 mm input-on-
  M2-back obscuration + FM mm-graze, both documented; M7 physical
  +126 mm — the report's new M7 flag is pupil-retrace bookkeeping,
  not light).  s2_wfe_field.png now ±3′ with the 0.02 contour +
  science patch + chief/centroid markers.  s3 rerun on the new
  parent: pie 12098/12520 (376 gap clips), hex2 9830/9876, both
  VERIFIED.  Old hand-picked variants preserved in s2_variants/.
- **s4 LANDED (MACOS_res 89f06f8): `s4_jacobians.m`** — dwdx (60412×90,
  rank 90; spectrum shows EXACTLY 21 strong modes = 7 segs ×
  piston/tip/tilt then a 2-decade cliff), dwdz (×56, FFZern 4..11 per
  seg, cond 4.4e2), dwdgrid (×42, 6 aperture-confined pokes/seg on
  grid-augmented `e2e_pie_grid.in` — flat 256 grids in each segment's
  CLOCKED Mon frame; cond 1.7), all via the production dw_d*_multi
  supervisors over C+4 corners at ±2′.  Per-segment column-norm table
  (piston dominant, clocking ~null) in s4_report.txt; s4_jacobians.mat
  carries the three outputs for s5/s6.
  POST-s4 FIXES (Dave review, e701f87): SMM pie frames = axis of symmetry (wedge bisector / center flat-normal; hex heritage unchanged; refs regen, tSegMirMaker 3/0); s4 dwdz = monzern (segment-LOCAL, cond 438→5.4; ffzern = parent-aperture basis, wrong channel); s4 dwdgrid = segment_grid_basis G-S (cond 1.26); segment_grid_basis ray-hist fallback for Segment pm_ref_elt; LESSON = copy the tested sensitivities runners CONFIGURATION (run_dwdz, run_dwdgrid_multi_multisegbasis), not just the API.
  NEXT: s5_met.m — MET config with the SHAPE-CLASS launcher constraint
  (below) + dedx/dldx joining s4's dwdx; then s6 simulator.

**2026-07-19 (s6 SESSION, post-compact resume point): run_compare
SHIPPED — the compare stage runner + s6 driver + tRunCompare.**

- `design/runners/run_compare.m` (Dave's spec): pokes each rigid DOF
  of the sensed bodies (100 nrad / 100 nm) on the optimized met Rx;
  per poke TWO graphics — mmacos ENGINE | LINEAR MODEL — each the
  center-field OPD change (parula, shared clim) above stacked bars of
  l / e_piston / e_gap / e_shear (shared ylims); dwell 1.6 s (Dave:
  0.25 read too fast — runner default + driver + GIF DelayTime all
  1.6); frames in <name>_frames/ + <name>_compare.gif + agreement
  report + <name>_compare.mat (exports dwdu = segments+SM columns for
  s7; ~7 MB, committable).
- Engine truth per channel: w = RigidBodyChannel apply → trace(nElt−1)
  → OPD delta (same convention as the s4 harvest; plain bodies need no
  FP tracking); l = macos.met() (met points ride the element under
  programmatic perturb); e = FINITE-ROTATION kinematics at the Hx
  SensorPos points — axis recovered from the row's own translation
  block (a_w = T_s·blkᵀ), arm = SensorPos − rpt_s, full Rodrigues;
  first order == the dedx row, difference = true linearization error.
- **STALE-dedx catch (the debugging story):** first run showed e
  disagreeing O(1) on in-plane DOFs while translations are
  algebraically forced to agree → the s5 e2e_pie_met.mat dedx was ONE
  Hx GENERATION BEHIND (s5 ran before the gap/shear axis-relabel
  regen; piston rows agreed, in-plane pair rotated).  run_compare now
  builds dedx FROM THE Hx SIDECAR (pad + rot cols ×cbm, same as
  run_met) and only cross-checks the met .mat, warning
  'run_compare:stale_dedx'.  s5 rerun queued/in-flight to refresh the
  .mat estimator products (dxde/dwde) for s7 — WEM numbers are
  rowspace-invariant, layout unchanged.
- **s6 RESULT (e2e pie, 54 DOFs): engine == linear to ≤9.1e-5 (w),
  ≤1.2e-6 (l), ≤6e-8 (e rot; translations to machine 1e-16).**
  Null-response DOFs get a floor (opts.w_floor 1e-12 m): segment Rz
  (clocking through the near-flat parent leaves ~0.01 nm — REAL, above
  floor; true nulls are the flat fold's in-plane DOFs + FPA-class
  bodies at ~5e-15 m FD noise) → reported 'null', not a noise ratio.
- Bodies without dwdx channels (aft ring riding the FPA in the e5
  fixture) get ZERO linear columns — matching the engine's null plain-
  trace response; control bodies are asserted actuatable.
- tRunCompare (SUITE_FAST): e5 pie fixture, REAL single-field coarse
  harvest ('grid','1x1', ngridpts 15, fp_mode none — no ApStop in the
  SMM corpus) + run_met-with-no-jac met struct; gates = the physics
  (w<5%, l<1%, e<1e-3, FPA null, dwdu 48 cols).  36 s.  PASSES.
- **s7 NOTES (Dave, this session): single-step STATIC estimator** —
  the OSE form x̂ = x̄ + K(m−m̄), converged steady-state gains (= the
  dxde/dxdl MMSE partitions run_met exports), m = [l; e]; OSC
  controller u = −[c_wu·I + (DwDu)ᵀDwDu]⁻¹(DwDu)ᵀ·dwdx·x̂.
  Background papers: MACOS_sandbox/Documents/OSE_Eqns_2019.pdf +
  2025_JATIS_HWO_Special_Issue-2.pdf (read OSE; JATIS pending).
- s5 rerun DONE: e2e_pie_met.mat refreshed on the current Hx — WEMs
  reproduce EXACTLY (as-built 15.57/3.918, optimized 3.766/1.777, FD
  0.00%, MC 1.9%) confirming rowspace invariance of the relabel.

**z + grid EXTENSION + FIGURE SENSING (Dave, same session):** "add
dwdz and dwdgrid visualization to s6 — keep going after dwdx, z then
grid" and "evaluate the effect of grid and z DOFs at each sensor
location (l and episton) → generate dmdz and dmdgrid":
- run_compare now phases x → z (MonZernChannel pokes on the met Rx) →
  grid (GridChannel pokes on og.rx_path, basis REBUILT via the same
  segment_grid_basis call — mismatch fails the agreement gate).  One
  GIF, frames p%03d, per-channel worst summary.
- **macos.design.dmet_dfig** (+ macos.design.zern_seg_eval): dmdz /
  dmdgrid = mode shape at each Hx SensorPos + met_geom launcher point
  (src_elt gives the launcher→segment map engine-truth), projected on
  the measurement axis: dl = −û·n̂·f(p_src), de = (â·n̂)·f(p_q).
  Piston + l rows carry it; gap/shear are SLOPE-ORDER small (~2-3% on
  the e5 f/1.75 parent), NOT zero — in-plane axes are ⊥ the LOCAL
  normal, Mon sag displaces along the FACE normal zMon (real model
  response; flagged for Dave's model review).  Exported in
  s6_compare.mat as [l;e]-ordered dmdz/dmdgrid for the s7 H.
- **zern_seg_eval convention PINNED BY ENGINE GATE** (tRunCompare
  test_zern_grid_engine_equivalence): the same MonZern mode sampled
  onto a grid channel + poked via elt_grid_add reproduces the
  MonZernCoef poke's OPD (scale within 2%, corr 0.998 = grid
  discretization) → lMon-normalized, UN-normalized ANSI (MonZernType=
  ANSI; NORM_RMS gates to Norm* types only — the E1 fix), Mon frame =
  clocked face frame.
- **UNIT GOTCHA (found via 999×≈1/cbm mismatch on the mm fixture):**
  the dwdz/dwdgrid supervisors (dwdz_for_current_source) return
  OPD-BaseUnits per coef-BaseUnits — NO cbm scale, UNLIKE dwdx's
  OPD-metres convention.  run_compare ×cbm on the figure-poke linear
  maps; invisible on e2e (m), 1000× on e5 (mm).
- The engine METcalc/Hx hold met + sensor points RIGID → z/grid
  engine bars read zero vs the dmdz/dmdgrid linear bars — the movie
  shows the model-vs-engine gap explicitly (labeled 'l mdl'/'e_pist
  mdl' in the report; natural follow-up = engine met points riding
  figure, Dave's call).
- tRunCompare EXTENDED (z leg with a real monzern mini-harvest — the
  e5 FreeForm REFRACTOR elt 9 legitimately sweeps in, counts dynamic
  — + the engine equivalence gate): 2/2 green, ~60 s.
- **.doc REPORTS (Dave):** MACOS_sandbox/e2e_reports/make_reports.py
  (pandoc) → s1–s6_report.docx: front matter + stage summary + stage
  figures + verbatim report; s5 carries the MET-OPTIMIZER DEEP DIVE
  (merit/WEM algorithms, shape-class search flow, code flow, results
  tables — s5_met_optimizer_deepdive.md).  s1–s5 built; s6 after its
  rerun.
- **GRID-BASIS PERSISTENCE (the deep find):** run_compare's grid
  agreement exposed that the s4 og dwdgrid block could NOT be
  reproduced by a rebuilt basis — modes 1–5 bit-exact, but each
  segment's LAST G-S mode came out ORTHOGONAL across sessions (|Δ|=√2;
  in-session double-builds bit-identical; raw detrended modes well-
  conditioned sv≥0.19 — a discrete tie-break below the G-S waterline,
  root cause unchased).  Fix = ARCHITECTURAL: the influence basis is
  PART OF THE JACOBIAN'S DEFINITION → run_sensitivities persists it
  (og.sgb) and run_compare/dmet_dfig consume it VERBATIM (rebuild
  path kept with a loud warning).  Also: run_compare must RELOAD the
  grid Rx after any segment_grid_basis call (it traces/ray-hist-pokes;
  leftover state biased grid FDs ~20% — the harvest supervisor
  reloads, manual loops must too).
- **FINAL s6 (persisted basis): x ≤2.8e-3, z ≤1.8e-3, grid ≤0.235.**
  The grid residual is an **OPEN ENGINE QUESTION for Dave**: the
  engine's grid-surface response is NOT proportional across sub-µm
  amplitudes on METRE-based Rx — central FD at the harvest delta
  (1e-6 BU) reproduces the Jacobian to 1.7e-5, but forward pokes give
  dW(h)/h drifting 0.15→0.24→0.42→0.98 over h=1e-8..1e-6 and BLOWING
  UP 84× at h=1e-5 (10 µm figure → ~570 µm OPD response?!); fex
  referencing ruled out; the mm fixture at the same physical poke
  agrees <5% → units-dependent (absolute-EPS class?) behavior in the
  FreeForm grid intersection path (GSZPSolve/SFFZPSolve suspects).
  x/z channels are clean at the same amplitudes (z fwd 1e-7 vs 1e-6
  agree to 1.5e-3).  run_compare report carries the note.
- Fast suite FINAL: **242/0** (tRunCompare = 2 tests, x/z/grid legs +
  engine convention gate).
- .doc reports (Dave): s4 doc gained the SENSITIVITIES GENERATION
  STORY deep dive (channels/referencing/grid substrate/basis
  persistence/unit conventions/validation); standalone CONCISE
  met_optimizer_concise.docx (with the optimizers-used note: discrete
  block coordinate descent + exhaustive per-class enumeration +
  top-K refinement, closed-form MMSE inside, no NLP solver).
  All in ~/dev/MACOS_sandbox/e2e_reports/ (make_reports.py).

**2026-07-20 (s7 SIMULATOR SESSION): run_simulator SHIPPED (Dave's
two-part spec) + the grid non-proportionality question ANSWERED
(numerical, not small-angle).**
- **Q2 answered (Dave asked if the grid rel 0.235 is just too-large-
  poke nonlinearity): NO — three discriminators.** (1) 100 nm on the
  8 m aperture is ~12 ppb sag and the z path is proportional to
  1.5e-3 at the same amplitudes; (2) the FD ratio does NOT converge
  as h→0 (0.15 at 1e-8 → 0.98 at 1e-6) — real nonlinearity improves
  as h shrinks, a granularity/EPS floor doesn't; (3) the mm fixture
  at the SAME PHYSICAL poke is clean <5% — physics can't see the Rx
  units, an absolute solver tolerance can.  Engine investigation
  still open (GSZPSolve/SFFZPSolve).
- **run_simulator (design/runners/, SHIPPED): time-series x/z/grid →
  movie, UNCORRECTED + CORRECTED legs** (Dave 2026-07-20: history
  opens with µm-to-mm misalignments; initial image-based WFC
  u = −pinv(dwdu)·w(frame 1) solved from the ENGINE wavefront and
  HELD; then nm-to-µm random-walk steps; 1000 s at 1–10 s steps).
  Frames: OPD unc | OPD corr | PIX psf | COMPOSE broadband psf +
  m=[l;e] bars (piston/gap/shear colors) + ACCUMULATING rms-WFE
  (log, both legs) and Strehl curves (unc λ0/corr λ0/corr bb, peak
  ratio to nominal).  m bars = the validated linear model
  [dldx;dedx]·(x+u) + dmdz·z + dmdgrid·g (dedx rebuilt from Hx;
  dmet_dfig blocks; grid Rx has no metrology — engine-l cross-check
  runs when the history has no grid states).  Per-frame corrected-leg
  engine-vs-linear w_rel = the s6 gate extended to MIXED states,
  computed on the DRIFT INCREMENT (frame t − frame 1): the absolute
  frame-1 state is deliberately nonlinear at µm-mm amplitudes and the
  control compensates it, so an absolute comparison would only
  re-measure that compensation.
- **WFC ITERATES (wfc_iters, default 3):** one linear solve at 202 µm
  left an ENGINE residual of 1.3 µm vs 40 nm predicted — the ~0.5%
  per-column nonlinearity of µm-mm states, not a solve error.  Real
  image-based WFC iterates; each Gauss–Newton refinement re-measures
  the engine wavefront at the corrected state and ridge-solves the
  update (monotone state path — never toggles back, respecting the
  two-pass non-closure rule).
- **TWO-PASS ENGINE SCHEDULE (correctness find):** toggling ±u per
  frame through the incremental perturb path does NOT close — fixed-
  order single-axis rotation increments leave a SYSTEMATIC ~|u_rot|²
  non-closure per cycle that accumulates LINEARLY (~µrad phantom over
  100 frames at 50 µrad u).  run_simulator plays the whole
  uncorrected history (storing per-frame OPD + psf peak), reloads the
  Rx, then plays the corrected history — within a pass the large
  state applies once and increments are nm-scale.
- **WFC SOLVE = TIKHONOV RIDGE (both failure modes hit, fixture +
  e2e):** plain pinv (tol 1e-6) noise-amplified near-degenerate dwdu
  combos (piston-like Tz directions) into huge canceling commands —
  engine honors them only to first order → WORSE psf (0.72→0.21)
  despite lower fitted rms; a hard SVD cutoff at 1e-2 over-truncated
  (624 nm of the 202 µm e2e initial error uncorrected → corr Strehl
  0.25 at 500 nm).  Ridge u = −V·(sv/(sv²+λ²))·Uᵀw, λ = wfc_tol·s1
  (default 1e-3; e2e driver tunes 3e-4 → predicted residual ~30 nm)
  corrects every direction sv ≫ λ and bounds each command by
  |w|/(2λ).  This IS the OSC controller form — carries straight into
  the estimator/controller loop.
- tRunCompare/test_run_simulator_time_history (SUITE_FAST, fixture
  reuse): µm-scale opening state + drift + z + g; gates = control
  collapse (corr < 0.2×unc), Strehl improvement, mixed-state w_rel
  <0.15, artifacts/m_hist/dmdz/dmdgrid dims.  3/3 green.
- **METROLOGY LOOP = RBCS estimator/controller (Tesch "RBCS
  Algorithms" ch 2.3 + 3.3; several live-review corrections):**
  `met_loop` (default on): the post-WFC state IS the control TARGET;
  each frame the sensed drift δm = m − m_ref drives a POSE ESTIMATOR
  then a CONTROLLER (Dave: "estimate the state without weighting the
  WF impact").
  - **Estimator = weighted-LS / BLUE (Tesch §2.3.2 eq 11):**
    δx̂ = R_meas·δm, R_meas = (Hᵀ N⁻¹ H + R_x⁻¹)⁻¹ Hᵀ N⁻¹, H = dmdx,
    N = sensor-noise cov (default [1 pm laser truss, 1 nm edge]),
    R_x = state/disturbance prior (default "auto" from the drift
    std, PTT ≫ lateral ≫ pinned-aft).
  - **Controller = min pose error (Tesch §3.3.1 eq 16-17):**
    u_t = u_{t−1} − k_p·δx̂(control DOFs), k_p 0.5 (<1 for margin;
    the TCE integrator makes the loop robust to gain error).
  - **THE BUG (Dave: "MET control is not correctly implemented" →
    pointed at Tesch):** my first loop used a RAW pinv(dmdx) = the
    basic-LS estimator (Tesch §2.3.1 eq 10).  It amplifies un-modelled
    δm content (figure drift, linearization residual) by 1/σ_min in
    the weakly-observed rigid directions; the integrating loop RAN
    AWAY — the engine e2e went 0.02 nm at t=10 → 2e5 nm at t=20 →
    5.3e6 nm.  Standalone linear reproduction: basic-LS max 5.6e6 nm
    WF residual / 34 mm commands vs BLUE max 6.98 nm / 74 nm commands
    on the identical drift (open-loop drift was 81.6 nm).  The R_x
    prior is the fix: weak/pinned DOFs fall back to prior-0 instead of
    being inverted.  This is STATE weighting (noise + disturbance
    stats), NOT wavefront-impact weighting (Tesch §3.3.2 eq 19,
    deliberately unused per Dave).
  - Bars show the SENSED DRIFT δm, not absolute m (µm-scale post-WFC
    offsets swamped the nm drift — Dave's "MET results not changing").
    Figure drift ALIASES into x̂ via the l/e_piston rows (no figure
    states in the simple estimator) — negligible at realistic figure
    drift; s7b's H adds dmdz/dmdgrid.
  - Loop-STABILITY gate added to tRunCompare: corrected rms ≤
    uncorrected rms every frame (the old pinv went 100× above).
- **DOF-class statistics (Dave):** support structure puts 10× more
  error in PTT (local Tz/Rx/Ry — the class MET+edge find and correct)
  than lateral (Tx/Ty/Rz); equal allocation had loaded the ridge-
  dropped weak directions and left a pessimistic ~40 nm one-shot
  residual.  **Figure drift = few nm TOTAL** (1 nm/step had the 98
  coefs each walking to ~10 nm → >100 nm rigid-uncontrollable figure
  dominating the corrected leg's drift-away).
- e2e driver s7_simulate.m (final): ~1 µm PTT initial (0.25 µrad
  tip/tilt, 1 µm piston; lateral /10), drift 4 nrad-nm/step (5×
  reduced; lateral /10), figure 0.004 nm/step, T=100 @ 10 s; M3→FPA
  + aft ring pinned (truly stable structure, Dave).  Running WFE
  printed in the OPD panels' xlabels (Dave).  Stale-frame cleanup
  built into the runner (s6 lesson).
- **WF-MAINTENANCE RECONTROL + TWO SCENARIOS (Dave 2026-07-21):**
  `wfc_reset_times` re-runs the image-based WFC mid-history (Tesch's
  WF Maintenance Activity) + `wfc_on_frame` delays the initial WFC so
  the movie opens uninitialised ("no system starts perfect" — first
  two data points at the as-deployed ~100 µm, then control turns on).
  Two 500s runs, reset@400s, one GIF each (s7A/s7B):
  - **s7A metrology-bias:** the metrology zero-point drifts; the loop
    holds the BIASED reading so the true wavefront walks off UNSEEN
    (`meas_bias = −dmdx·p(t)`); the 400s image-based reset
    re-references and knocks it back.  ENGINE: corr holds ~1 nm →
    **44 nm by 390s → reset → 1.0 nm** → 12 nm by 500s.
  - **s7B focus/astig figure:** per-segment focus(5)+astig(4,6) trend
    (60 nm, 2× for visibility per Dave); the truss reads RIGID POSE
    only (`loop_senses_figure=false`) so figure accumulates unseen; the
    tight-ridge reset (`wfc_reset_tol=1e-5`) engages the LATERAL DOFs.
    ENGINE: 20 nm → **45.6 nm by 390s → reset → 24.4 nm** → 31 nm.
  - **PLOT conventions (Dave, [[feedback_demo_plot_conventions]]):**
    delayed init (first 2 pts at as-deployed 100 µm); Strehl EXACT from
    OPD `|<exp(i2πW/λ)>|²` (≤1, not psf-peak >1); autoscale the
    corrected OPD panel (shared 100 µm scale hid the nm structure);
    broadband Strehl trace dropped; legends non-obscuring (hide the
    reset-marker auto-legend entry); running WFE in the OPD xlabels.
  - **KEY PHYSICS (Dave, corrects my earlier error):** RB control CAN
    counter SEGMENT focus/astig on a parabolic parent — a segment x
    move changes its local best-fit radius (focus + a bit astig),
    y/twist add astig.  Verified on the s4 dwdx: **focus RB-residual
    0.017, astig 0.31–0.49 (via lateral DOFs), higher order 0.6–0.98
    (uncorrectable)**; the reset needs tol 1e-5 to reach the weak
    lateral DOFs (3e-4 leaves astig 0.78).  My earlier per-segment
    reset failed 24→24 because the loop sensed+pre-consumed the figure
    AND the loose reset tol missed the lateral DOFs.
  - There is NO wavefront error BOTH MET-invisible AND rigid-
    correctable in a well-instrumented segmented system (rigid-
    correctable ⟹ segment pose ⟹ MET-visible); the WFC reset earns its
    keep against (A) MET bias drift and (B) figure the MET is DEFINED
    not to sense.
- NEXT (s7b): upgrade the estimator to the STEADY-STATE KALMAN form
  (Tesch §2.3.3 eq 12-14, predict/update with the Riccati gain) +
  figure states via dmdz/dmdgrid in H + sensor noise; the OSE
  single-step static estimator is this with converged gains.

**2026-07-19 (RECAST SESSION, earlier): the
sensitivity diagnostic + runner recast + SMM EDGE-SENSOR REWORK all
LANDED (pushed: MACOS_res 20d7f03+b41503d, macos 3a89432).**

DIAGNOSTIC (Dave's four questions, all closed):
- s4 == the runners BIT-IDENTICAL (same dw_d*_multi supervisors).
- Segment coordinates CORRECT (engine-truth poke footprints match
  frames on e2e AND e5pie; apparent reflection = canvas convention).
- Nominal subtraction CORRECT (reference frozen per field; median
  col-vs-w0 corr 1-3%; top col Elt 17 Rx 0.63 = physical FP signature).
- "e5 doesn't compare": the E5 CORPUS is the trap, not s4 --
  SegMirMaker replicates the PARENT's grid channel (pData=0, full-
  aperture span) into every segment block; segment-frame bases poked
  against those paint a CENTRAL DOT and rank-collapse (e5pie dwdgrid
  rank 15/42 cond 1e7 vs healthy 42/1.26).  ALSO: SMM corpus ships NO
  ApStop (exit-pupil machinery fails).  Last-key-wins parsing means
  appending correct grid lines cannot fix stale ones.

RECAST LANDED (MACOS_resources, all local):
- macos.design.grid_augment_rx: per-segment grid channels in the
  CLOCKED Mon frames, REPLACING stale lines.  Span default = the
  dxGrid convention (Dave): GridSrfdx = Aperture/(ng-1) (span = the
  beam; e5mono heritage 31.25 = 8000/256).  NEVER size from lMon:
  pie-wedge lMon is the hex 'length', not a circumscribing radius --
  a 2.2*lMon span clipped wedge corners and made s4 dwdgrid
  non-physical (Dave caught it; segment_grid_basis now WARNS when
  footprint rays fall outside the grid span).
- design/runners/run_sensitivities.m: general stage runner (dwdx/
  dwdz/dwdgrid + dwdsurf opt-in; grid_basis multi|single; influence
  passthrough for DM maps; ApStop injection preflight 'stop' opt;
  zkinds; segment-only conditioning; per-segment column norms).
- design/runners/run_segmentation.m: parent .in -> verified segmented
  .in (parity, engine-truth footprints, apertures, reload gate).
- s3/s4 = thin drivers; regen VALIDATED (e2e_pie.in/e2e_hex2.in
  byte-identical; s4 dwdx/dwdz numerics unchanged, dwdgrid healthy).
- PLOT CONVENTIONS (Dave): piston removed in the shared
  plot_dw_per_element/plot_dw_channels (parula stays); per-element
  pages center+multi in <name>_pages/ subfolder (page flood).
- sensitivities/run_dwd*_multi.m x6: PRESERVED as same-name,
  same-CONFIG thin wrappers over run_sensitivities (Dave: people use
  them -- do NOT delete); examples/run_* dirs = thin drivers too.
- Tests: tRunSensitivities (full-rank + poke-localization regression
  gates on the SMM trap corpus) + tRunSegmentation; both in
  SUITE_FAST.

SMM EDGE-SENSOR REWORK (Dave's spec, landed + validated):
- Per SHARED EDGE: 2 sensor locations at +/-SensorOff (new prompt,
  default 0.25*width; new stdin answer BEFORE the final Y --
  fixtures + segmirmaker_run 'sensor_off' updated) x 3 axes:
  1 piston = surface normal, 2 gap = in-plane PERP to edge, 3 shear =
  in-plane ALONG edge.  NO absolute-piston anchor row (not a
  measurement -- Dave).  Pie = 72 rows = 24 locations x 3 axes; hex2
  config-2 = 252 (42 edges x 6).  Prescriptions BYTE-IDENTICAL.
- QUIRKS FIXED (flagged for Dave): legacy rows reused rhoi for BOTH
  segments + projected triad-before-cross; new rows use each
  segment's OWN arm: dm/d(del_s)=+/-a'T_s, dm/d(th_s)=+/-(rho_s x
  a)'T_s.
- Hx.m now carries MeasAxis/MeasLoc/SensorPos; edge_sensors ingests
  (axis/loc/sensor_pos/has_anchor; legacy files still parse).
- tEdgeSensors REWRITTEN (axis recovery both sides, per-segment
  sensor-point coincidence, on-edge at +/-SensorOff; rows determine
  the point only PERP to the axis -- all checks project the axis
  out; in-plane axes are point-local).  tSegMirMaker refs
  regenerated.  3/0 + 3/0.
- s5 RERUN with 72-row dedx: as-built edge+MET WEM 10.81 -> 3.918,
  optimized 4.59 -> 1.777 (in-plane DOFs now edge-observable);
  MET-only optimized 3.766 unchanged (truss stands alone).  dldx FD
  1.1e-7 PASS, merit FD 0.00%.  (Hx radhat/tanhat->gap/shear relabel
  afterward leaves WEM invariant -- same rowspace, isotropic noise;
  e2e sidecar regenerated with final semantics.)

DAVE'S s6/s7 SPEC (VERBATIM INTENT, next to build):
- Linear model: w = dwdx*x + dwdz*z + dwdgrid*grid + dwdu*u + w0;
  u = control = SEGMENT + SM rigid-body DOFs (dwdu = those dwdx
  columns).  Measurement m = dmdx*x + m0, m = [l; e], dmdx = [dldx;
  dedx].
- s6 = run_compare: poke each DOF in turn (default 100 nm / 100
  nrad), display 2 graphics -- mmacos vs linear -- EACH = OPD map
  above stacked bar charts of l, e_piston, e_gap, e_shear; settable
  dwell (default 0.25 s); + saved frames + per-poke agreement report.
  Engine side: reuse dw_dx/met_calc paths for w/l; e engine-truth
  from finite rotations at es.sensor_pos.  s7 = run_simulator
  (estimator/controller, uncontrolled + controlled).

COMMIT GATE IN FLIGHT: s3 final Hx regen + fast suite running; on
green commit BOTH repos (macos: SegMirMaker untouched -- SMM lives in
MACOS_resources).  NO PUSH until Dave asks.  DEFERRED: hole->SM
chain rerun; run_design (s1/s2 recast); deletions list (inventory in
hand: freeform_unobscured retired-candidate, examples/design/oatma
orphan, e5hex1 old runners -- Dave must sign off; NEVER
coro_planet_demo); GMI mex relink; pymacos api exports.

**2026-07-19: RUNNERS DOCTRINE (Dave) + s5 BUILT AS THE GENERAL
run_met RUNNER.**  Dave: the PRODUCT is a small set of reusable stage
runners handing the .in from stage to stage — design → segmentation →
sensitivities → met → **compare** (NEW stage: step x/z/grid states,
emit w/e/l from mmacos AND the linear model) → simulate (needs an
ESTIMATOR/CONTROLLER; uncontrolled + controlled outputs).  Decision:
build s5/s6 as general runners from day one; recast s1–s4 into
run_design/run_segmentation/run_sensitivities in ONE pass after s6
(interfaces settle when the last consumer exists), then rethink the
mmacos file structure for new users + DELETE obsolete code/examples.
Also queued: run_dwdgrid_multi should adopt the ray-history footprint
approach.  Landed this session (mmacos, LOCAL until suites green):
- `design/runners/` (on the mmacos_setup path) + README (pipeline
  table, handoff contract = .in + declared sidecars) + **run_met.m**:
  as-built add_met → trace parity + aft-beam HOLE-margin table → dedx
  (Hx→SI: rot cols ×cbm ONLY) → dldx engine-FD vs analytic GATE →
  merit table → gains dxde/dxdl/dwde/dwdl → met_layout_opt →
  realized winner + engine-FD validation → MC acceptance → report +
  met_view/view_rx/metric figures + <name>_met.mat.
- `macos.design.met_layout_opt` = e5_seg_metopt v3 hoisted, SHAPE-
  CLASS aware: classes by boundary CONGRUENCE (vert count + edge
  lengths of seg_boundary polygons — pie → hexagon+wedge, hex → one,
  rxpoly imports classify free), one pattern per class about each
  member's own frame-x (pattern_frame 'segment'; 'radial' = e5
  heritage), coordinate-descent sweeps + top-K cross refinement.
  e5_seg_metopt.m is now a thin consumer — regen reproduces the v3
  winner EXACTLY (3.421 nm, FD 0.00%, same cluster angles; pmap/fclock
  differ by the fiducial-ring rotation symmetry; v3 baseline was
  ALREADY Inf — infeasible 10 mm corner sep, not a regression).
- `macos.design.seg_from_rx`: rehydrates the seg struct from the
  segmented .in ALONE (elt-type scan + met_bodies triads + lMon +
  src_seg_get tiling) so runners take files, not in-memory structs.
- `dldx_analytic` grew optional `unit_to_m` (default 1e-3 mm
  heritage): **e2e BaseUnits = METRES — the hard-coded mm arm scale
  would have squashed rotation rows 1000×.**  e2e Hx needs NO scaling
  (readings already metres); general rule in run_met: rot cols ×cbm.
- e2e s5 driver `s5_met.m` (thin): hub=M2 (elt 8), **aft ring on the
  FOLD FM (elt 9), NOT M3 — M3 is 2.1 m off-axis, its beams cross the
  M1 plane at ~2 m radius; FM is the on-axis bench body whose 0.10 m
  ring reaches the M2 rim fiducials THROUGH the 0.1225 m hole (run
  verified: crossing radii 0.079–0.10, margin 22 mm CLEAR)**.  54
  DOFs = 7 segs + M2 + FM.  Verified before the OOM (below): trace
  parity exact, FD gate 1.1e-7 PASS, 54 cols found.
- **OOM lesson (the two "VS Code crashed" reports = THIS)**:
  run_met's first merit form trace(Dk·P·Dk') with Dk 60412×54 built a
  29 GB dense → killed the 30 GB box twice.  Fixed to trace(P·G),
  G=D'D (memory feedback_trace_gram_not_outer); long MATLAB runs now
  caged via systemd-run MemoryMax.
- Tests: tMet +test_dldx_analytic_unit_to_m; NEW tRunMet (seg_from_rx
  rehydration identity; shape-class discovery [1 6] + WEDGE-PATTERN
  CONGRUENCE in segment frames; run_met end-to-end on the e5 pie
  fixture with synthetic jac + trimmed grids).  tRunMet added to
  SUITE_FAST.
- IN FLIGHT at session edge: caged s5 e2e run + tMet/tRunMet batch;
  then fast suite → commit both repos.

**DAVE'S s5 DESIGN REVIEW (2026-07-19, mid-session — supersedes the
first s5 numbers):** the first optimized layout (tight ±5° clusters,
WEM 4.9) "scores far worse than others we have evaluated" — it
maximizes neither beam angles nor outer-edge coverage.  Directives,
ALL IMPLEMENTED in run_met/met_layout_opt (rerun pending):
1. HEADLINE METRICS DIMENSIONLESS (WEM = wavefront per unit gauge
   noise, suite ratios held; floor as fraction of prior) — absolute
   nm at arbitrary priors "invites panic" (HWO context).  MET tracks
   the post-WFSC DRIFT state (≪1 nm class), gauge roadmap ~1 pm →
   e2e scenario sigmas now 1e-10 prior / 1e-11 edge / 1e-12 met.
2. Report WEM full / TILT / NON-TILT for MET-only AND edge+MET (tilt
   control absorbs tilt; orthogonal split via per-field piston+tilt
   projection of the s4 per-field blocks — Gram-only, no row order).
3. OPTIMIZE THE TRUSS WITHOUT EDGE SENSORS (must stand alone) on the
   NON-TILT Gram; edge sensors stay in the reported cases + estimator
   products.  Edge-sensor PLACEMENT is never optimized — SMM Hx as-is
   (the figures' open gray circles were the as-built LAUNCHER overlay;
   real sensor midpoints now drawn as gray dots via met_view
   'sensor_pts').
4. MANUAL BENCHMARK (auto 'corner_pairs' extra in met_layout_opt):
   launcher pairs at each class's two outermost boundary corners + a
   pair on the inside edge; evaluated GATE-BYPASSED with its true
   min-sep reported (adjacent-wedge corner pairs sit ~gap apart vs
   the 50 mm rule — a design datum, not hidden).
5. AFT STRUCTURE TO THE SM RADIUS (0.232 m ring on the fold/bench;
   'hole_r' override = SM shadow radius for the clearance check;
   "the hole will be that large anyway").  QUEUED (Dave's call):
   enlarge the REAL Rx hole via P.tel → full s1→s5 chain rerun.
6. AFT LEG SOLVED FIRST (Dave: "solve the simpler M3-SM truss first,
   fix it"): the aft ring's clocking + ITS OWN fiducial map = the
   first coordinate block in the descent, frozen for the class
   sweeps, revisited once.  add_met grew extra_clock/extra_pair_map.
7. PRESET LAYER (Dave: "save the optimized configuration as a new
   as-built for future PIE builds"): met_layout_opt exports every
   winner SCALE-FREE (class pattern angles, fiducial RIM INSET, aft
   block) + 'apply' mode realizes a preset on any same-tiling build;
   run_met 'preset' option makes it the as-built; s5 driver promotes
   to design/runners/presets/pie_met.mat; tRunMet round-trip test.
8. SYMMETRIC (ROTATIONAL) FIDUCIAL ASSIGNMENT + nf=6 (Dave: "with 6
   [fiducials] the segments can be ~interchangeable"): sym_assign
   (default ON) shifts each member's fiducial map by its clocking →
   congruent beam geometry per member; add_met accepts nseg×6
   pair_map; e2e forces nf_grid=6 (the nf=3 raw-merit winner broke
   the 60° symmetry).
**SESSION CLOSE-OUT (Dave 2026-07-19): pushed (MACOS_res 51c2f77 /
macos 5669ffe).  Decisions: (a) hole→SM-radius CHAIN APPROVED but
DEFERRED — P.hole_min_r_m floor is COMMITTED in s1/s2 [5], artifacts
NOT yet regenerated; run the chain (one MATLAB per stage — the §0
model-transition bug punishes mixed-size single processes; script
pattern in the transcript) AFTER the s3/s4 cleanup.  (b) NO corner
hardware sharing (corner_pairs stays a benchmark).  NEXT SESSION =
RECAST: s1–s4 → run_design/run_segmentation/run_sensitivities +
sensitivities/examples/* onto robust runners + file-structure rethink
+ obsolete deletion + run_dwdgrid_multi ray-history; AND investigate
Dave's flag: "improve the pie sensitivity results — they don't
compare well to the e5 examples" (pie dwdx/dwdz/dwdgrid quality vs
the e5 fixtures — start by diffing s4's spectra/conditioning and maps
against run_dwdx/dwdz/dwdgrid_multi outputs on e5).  Then s6
(run_compare first, then run_simulator w/ estimator/controller).**

FINAL RESULTS (s5 v3, nf=6 + sym assignment + solved aft block):
as-built WEM 15.6 MET-only / 10.8 edge+MET → **optimized 3.77 /
4.59, floor 0.28% of prior** — the SYMMETRIC nf=6 config BEATS the
asymmetric nf=3 raw-merit winner (4.86) outright: interchangeability
AND merit.  Dave's corner_pairs: 3.97 (ties the machine) but min-sep
15 mm violates the 50 mm rule at wedge-corner junctions.  WEM_tilt ≡
0 by CONVENTION: s4 dw_dx channels are fp_mode='track' harvested →
global tilt out at the source (prior tilt split 8e-13 nm);
per-segment tilt fully in the merit.  MC 0.7%, FD 0.00%.  Preset
promoted: design/runners/presets/pie_met.mat.  e5 regen under the
new optimizer: **3.421 → 3.255 nm, worst-mode 184 → 171 nm, FD
0.00%** (aft block + sym assignment — improvement, documented).

**s5 DESIGN CONSTRAINT (Dave 2026-07-18): MET-configuration
optimization solves ONE launcher pattern per segment SHAPE CLASS,
expressed in the segment frame, replicated to all same-shape segments
(pie: hexagon class + wedge class; hex2: one hexagon class) — the
tier-3 symmetry-first collapse, now shape-class-aware.  add_met
'launch_pts' + seg_boundary/rxpoly polygons are the realization
hooks.**

NEXT: s4_jacobians.m — dwdx/dwdz/dwdgrid on e2e_pie.in (dw_d*_multi
harvests; per-segment 6-DOF x ordering per Sprint-2D), then s5 MET
(shape-class patterns per above), s6 simulator.  Views + THOROUGH
report every stage.
RESUME RECIPE (post-compaction): (1) re-read root+nested CLAUDE.md,
MEMORY.md (project_e2e_example), this file; (2) artifacts live in
mmacos/design/examples/e2e/ (s1_telescope.in/.mat = DL telescope,
s2_instrument.in/.mat = Offner system, reports+PNGs beside them;
e2e_params.m = every knob; README = 12-rule design procedure);
(3) s3_segmentation.m = NEW runner: segment M1 of the e2e system
(segment_rx parent = s2_instrument.in, elt 1, P.seg block: Hex,
rings 1, gap 25 mm, emit_apertures, model 512 -- the e5_seg pattern,
watch the FF-parent + ApStop notes in project_sprint2d_segmentation);
then s4 dwdx/dwdz/dwdgrid (dw_d*_multi), s5 MET (add_met + metopt v3
+ met_bodies), s6 simulator (COMPOSE/psf + linear model switch).
Views + THOROUGH report at every stage (Dave's standing rule).**

**(superseded record of session 2 below)**
**2026-07-17 (e2e session): STAGE 1 LANDED (MACOS_res `e20d6b1`,
LOCAL).  `design/examples/e2e/` = the 6-stage worked example; all
knobs in `e2e_params.m`.  Dave's decisions this session: (a) user
specifies BASIC telescope parameters, f/# a FREE input (tma_layout,
NOT D-scaling); (b) the case = D=4 m, primary f/1.25, system f/18,
ON-AXIS Korsch taken SLIGHTLY OFF-AXIS; (c) 90-deg FOLD after M2
moves M3 + image + FP BEHIND M1; (d) every stage runner emits a
THOROUGH design report + graphics; (e) stage 2 = JOINT refinement:
as instrument optics are added, keep refining M1-M3 (M2/M3 apertures
may GROW, Zernike terms deepen) to improve field performance.
Stage-1 as-built: m2=16 -> f/20 int focus in front of M1 (met
injection), near-unit M3 relay, ~6% M2 obscuration, fold at z=0.3,
bench M3 [2.1 0 0.3] / FP [-0.5 -0.12 0.3], M1 hole r=0.182, bias
sweep -> 5' (least that fully clears; only M2 obstructs, by design);
solve ladder (each step was a debugged lesson): conics AT the bias
point with FP align FIRST (align between two solves; else -2 mm
field-curvature defocus = 1.6 waves poisons the conic solve) -> joint
FF field solve (1.07 worst) -> M1 stop common-mode null (bias point
0.0024 waves, Strehl 0.97; corners pay: worst +-1' = 1.58 raw/0.97
-tilt = the pure field DIFFERENTIAL stage 2 corrects).  ORDER
MATTERS: static-null-then-joint lands WORSE (1.48) than
joint-then-null -- LM basin path dependence.  NEW
`design/src/field_zone_lmon.m` (doctrine field-zone lMon; tested in
tDesignTelescope).  Reload 1305/1305; fast suite 225/0.
NEXT: s2_instrument (3-4 mirror relay widening toward +-2', joint
M1-M3 re-solve per (e)), then s3..s6 per the sequence below.**

**2026-07-12: PLAN_DESIGN_LAYER Sprint 2D — SEGMENTATION + SENSING
(Dave's objective): segment the design flow incl. per-segment local
coordinate systems, edge sensors, laser metrology; MET truss from M2
to points around M3; future = optimize the MET configuration.**

### The frame (Dave, 2026-07-12 — supersedes §6.6 tier-3 draft merit)
Deliverable = the **linear forward model of a segmented system on one
shared per-segment DOF vector x**, for closed-loop control analysis:
- optics `w = dwdx·x + dwdz·z + dwddm·dm + w0` (channels exist)
- edge sensors `e = dedx·x + e0` (dedx = SegMirMaker `Hx`)
- laser MET `l = dldx·x + l0` (poke → METcalc → FD; NEW channel)
- simple control `dx = −pinv(dwdx)·(w − wtarg)`, w from the ESTIMATE.
**Configuration merit = post-control wavefront residual.**  With
full actuation the correctable part cancels exactly and
`w_post = dwdx·(x − x̂)`, so
`merit = E‖w_post‖² = trace(dwdx · P_δx · dwdxᵀ)`,
`P_δx = X − X·Hᵀ(H·X·Hᵀ + R)⁻¹·H·X`, `H = [dedx; dldx]`, X = prior
cov of x (deploy tolerances), R = sensor noise.  Prior wavefront cov
`W = dwdx·X·dwdxᵀ` is the baseline to beat (report the ratio).  When
actuation ⊄ states, report the uncorrectable projection
`(I − dwdx·pinv(dwdx))·w` separately.  Unobservable directions
saturate at their prior wavefront cost (no singular inverse).  MET
placement optimization (tier 3, future) minimizes this merit over
launcher/fiducial placement + beam topology — pure MATLAB on the S4
outputs.

### Decisions (Dave, 2026-07-12)
1. Dev/test parent = **e5mono import** (SegMirMaker's canonical parent,
   committed references); showcase on the design-layer 3M after.
2. **Hx normal-height edge-sensor model accepted** (relative
   surface-normal displacement at edge midpoints — piston+dihedral, no
   in-plane gap/shear); richer sensor models later = pluggable dedx
   backends behind the same Jacobian contract.
3. **MET default geometry = Stewart-platform trusses**: 6 launchers on
   EACH segment illuminating 3–6 fiducials on M2; 6 launchers around
   M3 illuminating 3–6 fiducials on M2; one measurement per
   launcher/fiducial pair = change in straight-line distance;
   ≥ as many measurements as DOFs.

### Slices
- [x] **S0 — SegMirMaker refresh + batch driver (Q7).  DONE 2026-07-12,
  commit pending fast-suite gate.**  CMakeLists repointed
  build_release_giza→build_release, npsol/lapack/blas→slsqplib;
  **`-fp-model strict` added → a rebuilt binary reproduces the
  committed test_in references BYTE-IDENTICALLY** (first non-strict
  build differed at 1 ulp; strict closed it; verified 2× fresh-dir
  sha256).  Batch mode = scripted stdin (NO control-file mode exists;
  README's "nine questions" is really ~15–17 prompts; DOF default is
  3 not the documented 6 — README fixed).  Committed answer fixtures
  `test_in/e5pie.stdin`/`e5hex2.stdin`.  MATLAB driver
  `macos.design.segmirmaker_run` (fresh scratch dir; copies parent +
  macos_param.txt + GridFile= refs; positional answers; 'Done.' gate)
  + `tSegMirMaker` (3/3: pie + hex2 byte-identity, size-args error)
  wired into SUITE_FAST.  GOTCHA: rerunning in a used dir trips the
  overwrite prompt and shifts the answer stream — always fresh dir.
- [x] **S1 — `segment_rx` splice.  DONE 2026-07-12** (`macos.design.
  segment_rx` + `tSegmentRx` 4/4).  As-built: engine numbers elements
  by READ ORDER (`iElt=`/stale `nElt=` are cosmetic — e5hex2.in loads
  25 blocks despite nElt=24); splice replaces the parent's source
  GridType with the segment tiling + inserts nSeg/width/gap/SegXgrid/
  SegCoord after yGrid; nElt from BLOCK COUNT; downstream iElt
  renumbered cosmetically.  **ApStop never renumbers** (header form =
  3-vector StopPos POSITION msmacosio.inc:90; element form = 2-vector
  offset inside the stop element's own block :2887, last-wins) — but
  if the SEGMENTED element carried the element-form stop it drops
  with the block → warning + out.dropped_apstop.  OptTgtElt/RefElt/
  tMetElt refused pre-run.  Validation = load + trace parity (seg
  0 lost rays, rms 3.564e-5 vs parent 3.920e-5 mm = sampling) +
  frames orthonormal.  **FF-parent gotcha: ring centers are NOT
  3-D-equidistant (mm-scale astig figure) — tiling-plane invariants
  only.**  Original plan line: [ ] Telescope.segment() builder method
  — deferred to the 3M showcase.
  splice .presc blocks in place of the parent elt + downstream
  renumbering (emitter owns numbering; watch ApStop/element-index
  refs) → reload → validate (ValidatePrescription + segment-center
  ray check ≈ SEGRAYTRACE's: center-ray endpoint vs RptElt, tol
  ~1e-11).  Spec carries per-segment frames; `segment_frames()`
  readback.  Segment local csys = RptElt + TElt/pMon,xMon,yMon,zMon
  triad (SEGMENT = full parent surface: same conic/FF/psi/Vpt every
  segment; Mon slot zeroed = reserved per-segment figure channel;
  DOF columns = [rot·x̂ rot·ŷ rot·ẑ | trans·x̂ ŷ ẑ] in the triad).
- [x] **S2 — edge sensors (dedx).  DONE 2026-07-12** (`macos.design.
  edge_sensors` + `tEdgeSensors` 3/3).  Ingests Hx.m via isolated-
  workspace run(); pads row-sparse assignments to nMeas×nState; per-
  segment columns [rot_xyz | trans_xyz] IN THAT SEGMENT'S TRIAD;
  row 1 = master piston.  Validation WITHOUT surface re-eval, from
  the generator's algebra (SegMirMaker.f:588-645): translation
  triplet = normhatᵀ·T → unit norm + T_i·del_i' == −T_j·del_j'
  (same world normal both ways); rotation rows = shared-ρ cross
  terms (the generator uses rhoi for BOTH segments) → [del]ₓρ=th'
  solves ρ = pSeg_i − pr, pr lands on the shared-edge midpoint
  (lateral < 0.05 width; residual limited by Hx TEXT precision
  ~5e-10 — tolerances 1e-8/1e-6, not 1e-12).  Frame-mixing quirk
  noted: Fortran multiplies triad-component rows by world [ρ]ₓ —
  mirror the code, don't "fix" it silently (model review with Dave
  if it matters downstream).
- [x] **S3 — DONE 2026-07-12 (engine leg + add_met Stewart emitter, tMet 4/4; 48-beam truss == geometry to 1e-12)** (`met_calc`/`met_get` in
  macos_api_mod after ray_status_get; capacity mMetSrf 20→64 /
  mSysMetBeam 128→512 in elt_mod; codegen Path A → `macos.met()`
  veneer with SI/native units; `tMet` 3/3 = Q8 closed-form gate:
  exact baseline lengths on a hand-inserted m2→fpa 2-beam fixture,
  gauge Δ == −û·d LOS projection under global perturb (5e-6), ⊥-null
  < 1e-7).  Both engine trees + mex rebuilt.  **Build gotcha: makems
  must run from ~/dev/macos (repo root) — a bad cwd 'succeeds' via
  stale logs; mmacos make needs explicit
  MACOS_BUILD_DIR=~/dev/macos/build_release_gfortran with
  FC=gfortran (its default is still build_release_giza — fix
  sometime).**  REMAINING: builder `add_met(...)`:
  Stewart trusses per decision 3 (fiducial/launcher placement in the
  segment triad → global; emits nMetPos/SrfMetPos rows + tMetElt +
  metBeamFlg).  Engine (small): `met_calc` + `met_get` wrappers in
  macos_api_mod (codegen Path A → mmacos veneer `macos.met()`), and
  CAPACITY BUMP `elt_mod.F:332-336` mMetSrf 20→64, mSysMetBeam
  128→512 (19 segs + M2 + M3 = 21 surfaces; 19·6+6 = 120 beams —
  2-ring barely fits today, 3-ring doesn't).  Ship with Q8-style
  closed-form tests (gauge Δ == LOS projection of relative fiducial
  motion; orthogonal-motion null; reciprocity) = the §6.6 tier-2 gate
  for SrfMetCalc.  Document limitation: straight-line lengths, NO
  LOS/obscuration check (PLAN §4.5 ShowMetObscur deferred).
- [ ] **S4 — SHAPE PER DAVE 2026-07-12: `design/examples/e5_seg/` — a MODIFIABLE runner: e5mono.in -> segmented .in INCLUDING the MET points (segment_rx + add_met) -> dedx (edge_sensors) + dldx (dmet_dx FD channel) -> MET metric performance = trace(dwdx*P_dx*dwdx') vs prior W=dwdx*X*dwdx' ratio.  Example rules: save .in+.mat, figures in dir, README, NO exit(0).  Underlying pieces: `dmet_dx` channel + `forward_model()` + demo.**  FD
  Jacobian via perturb→METcalc→read (dw_dx supervisor pattern; the
  API perturb path CPERTURB_PROG DOES move SrfMetPos).
  `forward_model()` = {dwdx, dedx, dldx, w0,e0,l0} on ONE x ordering
  (per-segment 6-DOF local + M2/M3 rigid).  Worked example: segmented
  e5-class primary → forward model → draw x~X, measure with noise,
  estimate, control dx=−pinv(dwdx)·ŵ → simulated RMS w_post ==
  analytic trace(dwdx·P_δx·dwdxᵀ).  That demo = sprint acceptance.

### Engine facts (2026-07-12 surveys; file:line current that day)
- MET: `nMetPos`/`tMetElt`/`metBeamFlg` parsed msmacosio.inc:1937-1975
  (nMetPos MUST precede tMetElt or parser STOPs); global-frame points
  `SrfMetPos(3,48,20)`; ONE command `METcalc` (3-char min-match,
  macos_cmd_loop.inc:2511) → `SrfMetCalc` utilsub.F:1810 = **pure
  straight-line Euclidean distance** source-point→target-point;
  flat output `metMeasBuf(1:nMetMeas)` (GMI output #9, pflg(27)).
  No LoadStack entry, no api_mod wrapper, no mmacos surface yet.
- SrfMetPos moves under CPERTURB + CPERTURB_GRP + **CPERTURB_PROG
  (funcsub.F:372-380 — the API/mmacos/pymacos path, so FD works)**;
  CPRead / CPERTURB_2 / LnkEltCPERTURB do NOT (PLAN §4.5 gaps — avoid
  those paths for met work).
- Engine `EdgeSensors` keyword = parsed+stored+SAVEd but **dead** (no
  consumer; PLAN §4.5 recommends not building on it) — dedx comes
  from SegMirMaker Hx per §6.6 tier 1.  SAVE round-trips all met
  keys since 662e86e (nMetPos precedes tMetElt in the writer).
- SegMirMaker: standalone fixed-form Fortran, own CMake vs
  build_release; loads the parent through the REAL engine loader
  (needs macos_param.txt + GridFile data in cwd);
  `SEGRAYTRACE` = engine-side QA command (macos_cmd_loop.inc:1391).

### Tier-3 MET-layout optimization plan (Dave's Q 2026-07-12, answered)
Merit evaluation is ANALYTIC once frames are known: dldx rows are
closed-form LOS projections + moment arms (validated FD==analytic in
tMet), dedx fixed, dwdx fixed per optical design -> thousands of
layouts/sec of pure linear algebra, NO engine in the loop (engine FD
only validates the winner).  Approach = HIERARCHICAL, not raw
combinatorics: (1) SYMMETRY first -- one launcher pattern per segment,
mirror-symmetric about the segment centerline, replicated by the hex
rotation group -> collapses the discrete space to one segment's
pattern x global ring params; (2) COMBINATORIC enumeration of that
small symmetric set (affordable analytically); (3) STEP-AND-EVALUATE /
patternsearch polish of the continuous knobs (radii, angles, standoff,
nf) on the shortlist; (4) worst-mode check (max eig of P_w, not just
trace) + a symmetry-BREAKING perturbation stage if a symmetric blind
mode saturates; (5) engine-FD validation of the final layout.  BUILT 2026-07-12: dldx_analytic + e5_seg_metopt.m -- 20160 layouts/12.8 s, as-built 5.156 -> 3.732 nm, engine validation 0.00% off analytic.  Queued: hub/extra bodies in the analytic rows need engine TElt/RptElt frame parse; fid ring on M2 BACK side; beam-clearance check.

## In-session state NOT yet committed
- MACOS_resources: segmirmaker CMakeLists/makesegmirmaker.sh/README/
  CLAUDE.md edits + test_in/*.stdin + mmacos segmirmaker_run.m +
  tSegMirMaker.m + run_mmacos_tests.sh SUITE_FAST entry.  Fast suite
  running as the pre-commit gate → commit on green.
- macos: this CURRENT_SLICE rewrite only (no engine change yet;
  engine work starts in S3).

## Just tried / ruled out (with why)
- Byte-parity chase vs the retired giza-era binary — RESOLVED, not
  ruled out: `-fp-model strict` reproduces the committed references
  exactly; no reference refresh needed.
- Reusing a scratch dir for repeat runs — ruled out (overwrite prompt
  shifts the positional answer stream; driver always makes fresh).

## Parked threads (carry, do not lose)
- **eac5 thread (Dave 2026-07-07):** reflective e5 back end =
  `~/dev/tst_dir/eac5mono.in` m3 (R=0.328 m, lMon 20 mm) + m4
  (R=0.407 m, lMon 30 mm) compact mini-relay ~0.2 m past the f/21
  intermediate focus — the reference for the freeform_unobscured
  "+2" (supersedes single-M4).  File breakage: CRLF + tabs-in-TElt
  (l.146) + `/* */` comment blocks (l.151+); Phase-1 validator
  rejects ("blank line inside TElt block").  NEXT: does the parser
  know `/*`? clean a copy, load/trace, extract m3/m4 conjugates.
- **±1′ design fork (Dave's call pending):** (a) ANSI-45 mode-depth
  test via one Jacobian assembly + column-space projection (~15 min);
  (b) pupil-aware re-layout (Sprint 5); (c) re-scope field/λ.  The
  shipped freeform_unobscured 3M numbers rest on pathological
  surfaces (see [[project_zern_solve_doctrine]]); re-solve under the
  doctrine needs Dave's OK.  Probe scripts:
  `~/dev/MACOS_sandbox/freeform_4m_probes/`.
- Others' actions on the AppScan response: Luis closes #61 + deletes
  ft-dev (both repos); IT rescans + dispositions the 2 conftest FPs.

### Tier-3 AS LANDED (e06c254) + NEXT SESSION (Dave 2026-07-12)
Constrained optimizer SHIPPED: launchers ON the segment hex boundary
offset OUTWARD by EDGE_OFF=5 mm (clearance off the reflecting surface,
by construction; sign-flip for inside-the-rim), 3 MIRROR PAIRS about
each segment's radial centerline, free pair angles (center seg = xhat
line).  add_met 'launch_pts' override realizes any winner.  Hex run:
54400 layouts/34.7 s; edge ring 3.857 -> 2.973 nm rms (worst 215.8 ->
155.9); angular spread [30 90 150] CONFIRMED optimal by search; the
gain came from the HUB FIDUCIAL geometry (rfid 300->150, fclock 105);
engine-FD validation 0.00%.  e5_seg now GRID='Hex'.
**NEXT SESSION (Dave): (1) annotation + graphics for the e5_seg
example (layout views of the truss/launchers/fiducials, metric
figures); (2) more examples.  Queued: hub/extra bodies in analytic
rows (TElt/RptElt parse); fiducials on M2 back side; beam clearance;
full fast suite ran only per-class today (Dave deferred); PLAN
promote-on-land pass.**

### RESUME RECIPE (written for a fresh session/model — Opus next)
1. Read THIS file + memory `project_sprint2d_segmentation.md` +
   `mmacos/CLAUDE.md` (+ root CLAUDE.md per directive).  All Sprint-2D
   code is PUSHED (macos f740cf6 / MACOS_res e06c254 tips).
2. Regenerate working state (~8 min):
   `cd ~/dev/MACOS_resources/mmacos && matlab -batch "run('mmacos_setup.m'); run('design/examples/e5_seg/e5_seg.m'); run('design/examples/e5_seg/e5_seg_metopt.m'); exit(0)"`
   Expect: edge+MET 5.978 nm, MC ~2.3%; optimizer 3.857→2.973 nm,
   engine validation 0.00%.  Tests: `./run_mmacos_tests.sh tMet` (and
   tSegMirMaker/tSegmentRx/tEdgeSensors), 19 green.
3. NEXT TASK (Dave): annotation + graphics for e5_seg, then more
   examples.  Graphics hooks that already exist: `macos.draw_rays` /
   `Telescope.view_layout` (real ray bundle), `am.src_pts/tgt_pts`
   (3xN global truss endpoints — plot3 beams launcher→fiducial),
   `seg.frames` (segment centers/triads for hex outlines + labels),
   `e5_seg.mat`/`e5_seg_metopt.mat` (all Jacobians + metric tables).
   Figures land in the example dir (house rule; no exit(0) in
   examples — the -batch wrapper above supplies it).
4. TRAPS a fresh model WILL hit: macos.modify() after every poke
   (cached OPD); macos.opd() = no args, N×N, WaveUnits; run makems
   from ~/dev/macos ROOT; mmacos mex: make FC=gfortran
   MACOS_BUILD_DIR=~/dev/macos/build_release_gfortran; SegMirMaker
   scratch dirs must be FRESH (overwrite prompt shifts stdin answers);
   engine numbers elements by read order; ad-hoc triads for hub/fpa
   DON'T match engine TElt/RptElt (why analytic rows are seg-only).

## Next concrete step
Fast suite green → commit S0 (MACOS_resources + this file on macos)
→ S1 `segment()` splice: decide splice host (Telescope method vs
System.from_rx path — likely both share one private splice helper),
renumbering rules for downstream element-index keywords, and the
segment-center ray check via mmacos trace.

## Open micro-questions (slice-local)
- S1: which element-index-bearing keywords must renumber on splice
  besides nElt (ApStop index, tMetElt targets, OptTgtElt, RefElt,
  LnkElt/group refs?) — enumerate from msmacosio.inc when building.
- S3: met-point placement API — accept explicit point lists AND the
  Stewart preset; where do M2-fiducial defaults sit (radius on M2
  face)?
- S4: x ordering — segments first (6/seg, triad order rot|trans) then
  M2/M3 rigid?  Must match dw_dx channel DOF naming.

## Promote-on-land  →  then CLEAR this file
> Same commit as the `design-sprint-N` tag: move each item to its
> permanent home, then reset this file to the empty template below.
- [ ] PLAN_DESIGN_LAYER Sprint 2D checkboxes ticked (Q7, segment(),
      post-splice validation) + NEW: S2/S3/S4 sensing items recorded
- [ ] `CORE COMPLETE <date>` blockquote added to Sprint 2D
- [ ] §10 Decisions: Stewart MET default; post-control-residual merit
      (supersedes tier-3 min-σ(H)); Hx sensor scope — all dated
- [ ] PLAN.md §4.5: met wrapper/capacity items ticked when S3 lands
- [ ] CLAUDE.md / nested gotcha captured (segmirmaker batch gotchas →
      segmirmaker/CLAUDE.md — DONE in S0)
- [ ] agent MEMORY.md: sprint-2d project memory
- [ ] worked-example committed + named (S4 demo)
- [ ] **reset CURRENT_SLICE.md to empty template**

---

## Empty template (reset state — copy over the above on land)

```
## Active slice
- Sprint / item: —
- Plan anchor: —
- Branch / worktree: sls-dev + sls-dev @ —
- Definition of done (honest): —

### Tier-3 MET-layout optimization plan (Dave's Q 2026-07-12, answered)
Merit evaluation is ANALYTIC once frames are known: dldx rows are
closed-form LOS projections + moment arms (validated FD==analytic in
tMet), dedx fixed, dwdx fixed per optical design -> thousands of
layouts/sec of pure linear algebra, NO engine in the loop (engine FD
only validates the winner).  Approach = HIERARCHICAL, not raw
combinatorics: (1) SYMMETRY first -- one launcher pattern per segment,
mirror-symmetric about the segment centerline, replicated by the hex
rotation group -> collapses the discrete space to one segment's
pattern x global ring params; (2) COMBINATORIC enumeration of that
small symmetric set (affordable analytically); (3) STEP-AND-EVALUATE /
patternsearch polish of the continuous knobs (radii, angles, standoff,
nf) on the shortlist; (4) worst-mode check (max eig of P_w, not just
trace) + a symmetry-BREAKING perturbation stage if a symmetric blind
mode saturates; (5) engine-FD validation of the final layout.  BUILT 2026-07-12: dldx_analytic + e5_seg_metopt.m -- 20160 layouts/12.8 s, as-built 5.156 -> 3.732 nm, engine validation 0.00% off analytic.  Queued: hub/extra bodies in the analytic rows need engine TElt/RptElt frame parse; fid ring on M2 BACK side; beam-clearance check.

## In-session state NOT yet committed
—

## Just tried / ruled out (with why)
—

### Tier-3 AS LANDED (e06c254) + NEXT SESSION (Dave 2026-07-12)
Constrained optimizer SHIPPED: launchers ON the segment hex boundary
offset OUTWARD by EDGE_OFF=5 mm (clearance off the reflecting surface,
by construction; sign-flip for inside-the-rim), 3 MIRROR PAIRS about
each segment's radial centerline, free pair angles (center seg = xhat
line).  add_met 'launch_pts' override realizes any winner.  Hex run:
54400 layouts/34.7 s; edge ring 3.857 -> 2.973 nm rms (worst 215.8 ->
155.9); angular spread [30 90 150] CONFIRMED optimal by search; the
gain came from the HUB FIDUCIAL geometry (rfid 300->150, fclock 105);
engine-FD validation 0.00%.  e5_seg now GRID='Hex'.
**NEXT SESSION (Dave): (1) annotation + graphics for the e5_seg
example (layout views of the truss/launchers/fiducials, metric
figures); (2) more examples.  Queued: hub/extra bodies in analytic
rows (TElt/RptElt parse); fiducials on M2 back side; beam clearance;
full fast suite ran only per-class today (Dave deferred); PLAN
promote-on-land pass.**

### RESUME RECIPE (written for a fresh session/model — Opus next)
1. Read THIS file + memory `project_sprint2d_segmentation.md` +
   `mmacos/CLAUDE.md` (+ root CLAUDE.md per directive).  All Sprint-2D
   code is PUSHED (macos f740cf6 / MACOS_res e06c254 tips).
2. Regenerate working state (~8 min):
   `cd ~/dev/MACOS_resources/mmacos && matlab -batch "run('mmacos_setup.m'); run('design/examples/e5_seg/e5_seg.m'); run('design/examples/e5_seg/e5_seg_metopt.m'); exit(0)"`
   Expect: edge+MET 5.978 nm, MC ~2.3%; optimizer 3.857→2.973 nm,
   engine validation 0.00%.  Tests: `./run_mmacos_tests.sh tMet` (and
   tSegMirMaker/tSegmentRx/tEdgeSensors), 19 green.
3. NEXT TASK (Dave): annotation + graphics for e5_seg, then more
   examples.  Graphics hooks that already exist: `macos.draw_rays` /
   `Telescope.view_layout` (real ray bundle), `am.src_pts/tgt_pts`
   (3xN global truss endpoints — plot3 beams launcher→fiducial),
   `seg.frames` (segment centers/triads for hex outlines + labels),
   `e5_seg.mat`/`e5_seg_metopt.mat` (all Jacobians + metric tables).
   Figures land in the example dir (house rule; no exit(0) in
   examples — the -batch wrapper above supplies it).
4. TRAPS a fresh model WILL hit: macos.modify() after every poke
   (cached OPD); macos.opd() = no args, N×N, WaveUnits; run makems
   from ~/dev/macos ROOT; mmacos mex: make FC=gfortran
   MACOS_BUILD_DIR=~/dev/macos/build_release_gfortran; SegMirMaker
   scratch dirs must be FRESH (overwrite prompt shifts stdin answers);
   engine numbers elements by read order; ad-hoc triads for hub/fpa
   DON'T match engine TElt/RptElt (why analytic rows are seg-only).

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
