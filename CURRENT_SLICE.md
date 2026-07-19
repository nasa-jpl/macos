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
