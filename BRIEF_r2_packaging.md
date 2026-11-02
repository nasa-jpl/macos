# BRIEF for TO: rodgers2/afocal4 — package the 343 mm design, or price the gap

_Tasking: Terminal Opus.  Supervisor: CC.  Timing: fits AROUND demo
rehearsal support (~2026-09-01) — NOT demo-blocking: deck slide 12
already carries the honest qualifier ("packaging was not driven to
completion; the trade is the result, not the layout"), so Path C below
is pre-satisfied and success only UPGRADES the slide.  STOP at the
timebox.  Origin: Dave 2026-08-30, reviewing slide 12 — "the layout
claims buildability, but the back focal length is longer than the
M1–M2 spacing.  Make it compact by repeated folds, or redesign, or
note it was not driven to completion."
Cold-start reads: memory `project_afocal4_rodgers2` (S4 trade NOT
BUILDABLE — this brief is that finding's fix attempt) +
`project_fold_extraction` (add_fold + fold rules) +
`feedback_frame_before_angle`; `challenges/afocal4/` (FORM_STUDY /
RESULTS / STATUS_S4B / STATUS_S4C records + the b2long/b2_trade
decks); `MACOS_sandbox/slides/deck_rodgers_status.md` slide 5 (its own
footnote already flags "serious repackaging")._

## Non-destruction rule (absolute)

All work in a NEW subdir — `challenges/afocal4/packaging/` — with new
deck names, new figures, its own README.  NOTHING committed under
`challenges/afocal4/` or `challenges/rodgers2/` is modified, moved, or
regenerated: the existing decks, figures and reports are the published
evidence for the trade and for the status/keysight decks, and they
stand regardless of this brief's outcome.

## The problem, precisely

The family-2 winner at the 343 mm pupil distance
(`afocal4_b2long_343mm.in`) holds the LAST MIRROR behind M1 — the
constraint that was enforced — but the back focal path from the fold
to the cold stop/instrument is LONGER than the M1–M2 spacing, so the
"instrument volume behind the primary" is not actually enclosed by
the telescope's own envelope.  The record knew: the trade's own
footnote calls the constraint "a first screen, not a layout".

## Steps, in order

1. **Measure the gap from ENGINE truth** (never .in text-parsing — the
   corpus-indexing lesson): load the committed 343 mm deck, take every
   leg length from the traced chief (`ray_hist`), and state: M1–M2
   spacing, each post-fold leg, the total back focal path, and the
   overhang.  These numbers go in the README and (whatever else
   happens) onto the slide footnote.
2. **Path A — fold it down (preferred):** insert fold flats
   (`add_fold` machinery; folds are FLATS — verify the null: WFE
   before vs after each insertion identical to round-off, the
   check-a-null-operation-is-null lesson) to wrap the back focal path
   into a stated envelope behind M1.  Standing rules: fold AOI < 15°
   unless geometry forces more (then say so); clearances from the
   SIGNED clearance model (the oi_clear pattern in design/src), never
   by eye; beam–mirror floors ≥ a spec you state up front (suggest
   25 mm, the rodgers3 convention, unless the afocal4 record implies
   its own).  Deliverables: `packaging/afocal4_b2pack_343mm.in`,
   view_std layouts (the new content-framed view_std is live),
   a legs/floors/envelope-box table.
3. **Path B — bounded redesign (only if A cannot clear):** allow M3/M4
   spacings and fold positions to move under the SAME constraint set.
   If any powered surface or spacing moves, RE-MEASURE and RE-REPORT
   the trade quantities for the new point (WFE, pupil blur, wander,
   magnification variation) — the slide-12 table numbers are the
   343 mm committed point and must not silently drift; a new point is
   a NEW row, reported beside it.
4. **Path C — price the gap (always delivered if A/B do not close):**
   the step-1 numbers plus one annotated figure showing the long leg
   against the M1–M2 spacing.  The deck footnote then quotes the
   measured overhang instead of a qualitative sentence.

## Gates

- Fold-insertion null: |ΔWFE| at round-off for every added flat,
  asserted, not assumed.
- Signed clearance floors for every leg/obstacle pair in the packaged
  layout, printed in the verdict table (a floor AT the spec
  undershoots — ask the hinge above the gate, the endgame lesson).
- The packaged deck re-traces with zero lost rays and reproduces the
  committed 343 mm design's WFE/pupil numbers (Path A: bit-comparable
  to round-off; Path B: honestly re-reported).
- Any reported fold angle carries its reference frame
  (frame-before-angle).

## Deliverables

1. `challenges/afocal4/packaging/` — deck(s), figures, README
   (numbers-first: the step-1 gap table + the outcome table).
2. Delivery log at the foot of THIS brief (the standing pattern).
3. Memory update extending `project_afocal4_rodgers2`.
4. Slide-12 refresh = CC afterward, under the §5 gate (upgrade the
   footnote on success; quote the measured overhang on Path C).

## Traps (paid for once already)

- matlab -batch: script files + exit(0) in runners; one model size per
  process; MACOS_HOME=/home/dcr/dev/macos/macos_f90.
- Read-tool PNGs cache-stale on overwritten paths — unique filenames,
  verify by printed numbers.
- Work in MACOS_res_dev (the dev worktree); ~/dev/MACOS_resources is
  the demo checkout, detached — leave it alone.
- Commits LOCAL; Dave orders pushes.  No engine work.

---

## Delivery log — 2026-08-30

**Delivered.  Path A closes the depth question; step 1 found a larger,
independent buildability item that no fold can touch.**  Everything in
`MACOS_res_dev/mmacos/challenges/afocal4/packaging/` (new dir, new names,
own README).  Nothing under `challenges/afocal4/` or `challenges/rodgers2/`
was modified, moved or regenerated.  Local commits only; no pushes.

### 1. The gap, from engine truth (`afocal4_b2long_343mm.in`)

| | |
|---|---|
| M1–M2 spacing (the yardstick) | **1.0420 m** |
| deepest optic behind M1 (field mirror) | **1.8866 m** |
| **overhang** | **+0.8446 m = 1.81× the front span** |
| back focal path, M1 plane → interface | 2.8075 m (2.69×) |
| + the stated 1.000 m instrument | 3.8075 m (3.65×) |
| observatory envelope | Ø1.120 m × 3.825 m |
| primary's through-hole | radius ≥ 174.0 mm |

Chief legs: 1.0421 / 2.9249 / 0.5541 / 0.3724 m.  Rays 10 665 over the
0.5°×0.5° box, 0 lost.

### 2. THE BIGGER FINDING — the collimator stands in its own feed beam

Measured over the **field box** (a mirror must carry every field), the
M2 → field-mirror beam runs **through** the collimator's glass:
**−79.9 mm** with a 1.15×hull + 15 mm body, and **−55.4 mm against bare
lit glass with no allowance at all**.  At the box centre alone it is
−6.0 mm, i.e. it was always marginal.

Mechanism: the collimator sits inside the converging feed cone, which
passes it in a 27.8–55.6 mm annulus — fine for ONE field (10.8 mm all
round), but a monolithic collimator must cover its union footprint
(17.0 → **87.0 mm**), and that glass is where the other fields' feed
beams pass.  **A fold is an isometry and carries this across unchanged**
(asserted at 1e-12 mm).  It is a Path-B item, not a packaging one.
`afocal4_pack` could not see it: that gate checks the LAST mirror's exit
leg for beam-to-beam daylight, at the box centre, and never asks whether
a body is standing in a beam.

### 3. What the committed S4b single fold does — rebuilt on this deck

* **It does not touch the depth**: deepest optic 1.8866 m, overhang
  +0.8446 m — identical to unfolded, to the last digit.
* **It costs the shroud**: the instrument leaves radially and reaches
  1.390 m off axis → **Ø2.779 m × 2.626 m** against Ø1.120 × 3.825.
* **The flat does not fit where the gate put it**: `afocal4_pack`
  measured 17.5 mm of daylight against its own 15.0 mm margin and passed,
  but over the field box that flat's footprint is **103.7 mm** in radius
  (body 134.3) and it clips the FM→M3 feed beam by **−73.6 mm**.  The
  gate's margin is a number, not a body.

### 4. Path A — `packaging/afocal4_b2pack_343mm.in`

Four folds in the M2→field-mirror leg, `+z→+x→−z→−x→+z`, normals in x–z
so the y bias maps to y.  Stated: x_step 375, x_out 190, z_front 230,
m3_gap 100 mm.

| | unfolded | 1 fold (S4b recipe) | **4 folds** |
|---|---|---|---|
| deepest optic behind M1 | 1.887 m | 1.887 m | **0.893 m** |
| overhang vs M1–M2 | +0.845 (1.81×) | +0.845 (1.81×) | **−0.149 (0.86×)** |
| optics slab | +1.323…+1.887 | +1.323…+1.887 | **+0.230…+0.893** |
| optics radius | 0.296 m | 0.296 m | 0.529 m |
| instrument reach | 0.465 m, aft | **1.390 m, radial** | 0.518 m, aft |
| observatory envelope | Ø1.120 × 3.825 m | Ø2.779 × 2.626 m | **Ø1.120 × 2.832 m** |
| body floor, pre-existing | −79.9 mm | −79.9 mm | −79.9 mm |
| body floor, new flats | — | **−73.6 mm** | **+5.8 mm** |
| rays lost, field box | 0 | 0 | 0 |

**A metre of length removed for no diameter at all** — the packaged optics
sit inside the primary's own keep-out cylinder (0.529 < 0.560 m) and
inside the M1–M2 slab.  Stated, not smoothed: the new flats clear the
beams by only **+5.8 mm**; each is a ~300 mm flat (footprint radius
127–132 mm) decentred 100–146 mm from the axis it is mounted on, because
the ±76 mm field walk dominates the beam's own 35 mm radius; and four
45° flats in series is a polarization budget this study did not open.

### 5. Gates

* **Fold null, asserted**: every deck re-scored on `afocal4_score`.
  |ΔWFE| 3.07e-8 nm, |Δblur| 4.35e-12 µm, |Δbreathe| 4.14e-14 %,
  |Δwander| 3.58e-12 µm, |ΔM| 1.42e-14; pre-existing clearance floor
  equal to **1.44e-12 mm**, which is the sharper test (a merit can agree
  while the geometry moved).
* **Zero lost rays** on every deck over the whole field box.
* **Signed clearance floors** printed per pair, split pre-existing /
  fold-induced and body / leg-leg.
* **Frames named**: all angles and stations are global, sky at −z.

### 6. Two things the null test caught — both worth keeping

1. **A real hardware bound, found the hard way.** The first route scored
   674 nm off its parent on a deck an isometry cannot change.  Cause: a
   +90°/−90° fold PAIR has perpendicular planes that always intersect (at
   x = x_step/2); a ray landing beyond that line reflects toward a point
   behind it — the engine loses it, and physically the two flats are cut
   into each other.  **Each pair's step must exceed twice the beam's
   half-extent on its first flat, over the whole field box**: measured
   80.4 mm → 180.8 mm.  Empirically exact at 175/200 mm, 1.2 % off at
   150, 6.5 % off at 125.  Now measured and asserted in `pack_route`.
2. **A silent MATLAB trap in this study's own code**: `OK = cat(1, OK,
   hi.ok)` seeded from `[]` returns a DOUBLE array of ones/zeros, so
   `P(:,mask,j)` indexes BY VALUE — N copies of ray 1, right size,
   plausible centroid.  One whole run's clearances were the chief ray's.
   Seed logical accumulators explicitly.

### 7. Open / for Dave

* The collimator-in-the-feed-beam interference is the real buildability
  item and is **bigger than the depth the brief asked about**.  Path B,
  not packaging.  Cheapest candidates: a longer M2→FM leg (buys room,
  costs the depth just removed) or a smaller intermediate image.  Not
  attempted — outside this brief's scope.
* **Record discrepancy, flagged not fixed**: `RESULTS.md` §S4b.4 states
  the folded demonstration's interface pupil at [+0.304, −0.004, +0.614] m
  and an instrument slab +0.464…+0.764 m.  Measured from the committed
  `afocal4_b_final_folded.in`, the fold and the interface plane are both
  at **z = +1.3782 m**, pupil [+0.2483, −0.0051, +1.3782], slab
  +1.228…+1.528 m.  `check_record` reproduces it.  Which is authoritative
  is Dave's call.
* **Slide 12**: the measured replacements for its footnote are in the
  packaging README §1 and §3.  Outward-facing, so it waits on sign-off
  (`doc/STYLE_REPORTS.md` §5).
