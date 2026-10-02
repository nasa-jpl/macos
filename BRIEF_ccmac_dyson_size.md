# BRIEF -- CCMac: how small can the Dyson block be at the same performance?

For **CCMac** (Claude Code on Dave's Mac), written from cleared state
2026-10-02 by CC (Linux).  The dyson5 challenge (a VSWIR Dyson imaging
spectrometer designed from Joe's EMIT-class spec) has a design of record,
R4, whose silica block is 220 mm in radius and therefore 221 mm THICK.
Dave: "200 mm thick glass is super heavy and awkward.  Not to mention
non-uniform.  Can't we make it much smaller with the same performance?"
Nobody is on that question.  Terminal Opus (TO) is mid-telescope (beat 5b)
in the same challenge directory.  Your job: answer the size question with
engine scores, in NEW files only.

## 0. Get current (before reading code)

- Both repos on **`dev-candidate`**, `git pull` after Dave pushes this brief.
  You are current when `git log origin/dev-candidate` shows, in
  `nasa-jpl/macos`, "CALIB: asphere differential step scales with the
  coefficient" (f1d2617) and this brief; in `nasa-jpl/MACOS_resources`,
  "dyson5 beat 5: the telescope that feeds the slit" (e114c99).  If not, STOP
  and say so.
- **Rebuild the engine (gfortran) and relink the mmacos mex.**  Seven engine
  fixes landed 2026-09-30 / 10-01 and this work is wrong without them: the
  glass catalog (a `GlassElt=` element traced as AIR before), medium-aware
  propagation, chord-ruled curved gratings and their path-length jump, short
  `AsphCoef=` lines, the CALIB stride and asphere step.  `rm` the mex to force
  the relink and confirm the binary carries the fix
  (`strings src/mmacos.mex* | grep -c "optics restored to the"` must print 1).
- **Baseline gates, before touching anything:** `./run_mmacos_tests.sh
  tSpectrometerRx` (7), `tGratingImmersed`, `tGlassDispersion`, `tGratingOpl`,
  `tAsphCalib`.  All green on Linux (fast suite 510 / 0).  Report your counts.
- Read, in this order: `mmacos/challenges/dyson5/README.md`;
  `BRIEF_dyson5_beat3b.md` (R4) and `BRIEF_dyson5_beat4b.md` (the meniscus
  search and what the merit weights do); `BRIEF_dyson5_beat4e.md` (the
  closure envelope -- the closest prior art to your task); `dyson_ladder.m`'s
  header; `dyson5_envelope.m` (`point_` is the pattern you will re-use);
  `macos/demo_session/deck_dyson.md` slides "The concentric Dyson and its
  scaling law" and "Where the design closes"; `macos/macos_f90/CLAUDE.md`
  sections dated 2026-09-30 / 10-01.

## 1. What is known

- **The scaling law** (stage s0): a concentric Dyson's corner blur goes as
  h^4.17 / r^3.19 (h = field height on the flat face, r = block radius).  In
  the Dyson form the flat face sits at the centre of curvature, so the
  block's THICKNESS equals its radius.  The 54 mm slit is 2.3x EMIT's.
- **R4 of record** (silica, r 220 mm, F/1.8, 54 mm slit, 18 um pixels):
  keystone 0.0026 px, smile 0.0051 px, CRF 1.327 px, SRF 2.032 px, ensquared
  energy 0.759.  Deck `dyson5_s3_r4.in`.
- **The closure envelope** (beat 4e, `dyson5_s4env.txt`), each point an R4
  re-solve from the record at 30 iterations:

  | point | CRF (px) | EE | note |
  |---|---|---|---|
  | silica r 180 mm | 1.607 | 0.606 | fails CRF |
  | silica r 150 mm | 2.371 | 0.408 | fails CRF |
  | CaF2 r 220 mm | 1.032 | 1.000 | large margin |
  | slit 30 mm, silica r 220 | 1.047 | 0.955 | large margin |
  | F/2.2, silica r 220 | 1.148 | 0.954 | margin; the spec is F/1.8 |

- **Two things the record does NOT know.**  (a) Whether the 180 / 150 mm
  rows are the design's limit or the SOLVER's: every envelope point was
  warm-started from the 220 mm record, and TO has just watched a cold solve
  stall outside its basin on the telescope (a stalled solve reads as a bad
  design).  (b) How far the CaF2 and short-slit margins convert into a
  smaller radius: the radius axis was scanned in silica at the full slit only.
- **A correction to the deck's glass number.**  The trade table's 22.7 L is
  the full hemispherical cap (50 kg in silica).  The beams use a core 147 mm
  in diameter (`blockD` in `dyson5_s3_trade.txt`), so the EDGED part is a rod
  147 mm x 221 mm: 3.7 L, 8 kg.  Report edged volume and mass.

## 2. The work

**New files only**, all in `mmacos/challenges/dyson5/`:
`dyson5_size_trade.m` (the driver), `dyson5_size_fig.m`, records
`dyson5_size.txt` / `.mat` / `.png`, decks `dyson5_size_<family>_r<mm>.in`,
report `BRIEF_dyson5_size.md`.  **Do NOT edit `dyson5_run.m`,
`dyson5_params.m`, `dyson5_envelope.m` or `dyson_ladder.m`**: TO has
uncommitted work in the first two and the others are the record's tools.
Take parameters from `dyson5_params()` and override fields in your own
struct; load the R4 seed the way `stage_s4env_` in `dyson5_run.m` does (read
it, do not call it: it rewrites the envelope's record files).  If you find
you need a change in one of those four files, STOP and report what and why.

**Method: CONTINUATION.**  Each point is one `dyson_ladder(Q, tagp, 'rungs',
5, 'seed', <previous point's r.P>, ...)` call with `point_`'s options and
bounds, warm-started from the PREVIOUS point's solved design, stepping the
radius down ~10 % at a time.  Score in the ENGINE (the ladder's `r.engine`);
chain numbers are not the record.  A point CLOSES when it matches R4 of
record: smile and keystone < 0.1 px, CRF <= 1.33 px, SRF <= 2.05 px, EE >=
0.76, and no variable on a bound.  A solve that ends on a bound is reported
as such and is not closed (name the variables).  Stop a walk two points
after it stops closing.

**Step 1 -- identity, then the solver check (report these first).**
1. Reproduce the record on the Mac: silica, r 220 mm, seeded from itself.
   The engine row must match R4's to solver tolerance; report the
   differences (this is also the Linux-vs-Mac engine identity check).
2. Silica, 54 mm slit, by continuation: 220 -> 200 -> 180 -> 165 -> 150 mm.
   Compare with beat 4e's 1.607 px (180) and 2.371 px (150).  If
   continuation beats them, the envelope's radius axis was solver-limited:
   say so plainly, it changes the deck.

**Step 2 -- the three families (walks downward by continuation):**
| family | glass | slit | radius walk |
|---|---|---|---|
| A | silica | 54 mm | step 1's walk |
| B | CaF2 | 54 mm | 220 -> 100 mm |
| C | silica and CaF2 | 27 mm (1500 pixels: two modules share the swath) | 220 -> 60 mm |
If time allows: F/2.0 at family B's smallest closing radius, as a margin row.
Start family B from the envelope's CaF2 220 mm deck's design
(`dyson5_s4env_glass_CaF2.in` is its emitted deck; re-solve the point to get
the parameter struct) and family C from the 30 mm-slit point walked to 27 mm.

**Per point, report:** the five engine scores; closes / fails / on-bounds;
block radius and thickness; clear diameter; EDGED volume and mass (silica
2.20, CaF2 3.18 g/cm^3); grating radius and diameter; slit plane to grating
vertex; the clearance gate's minimum and worst pair (`spectrometer_clearance`
as stage s3 runs it: the slit-to-detector separation does not shrink with
the block, so the gate may bind before the image does); and two
order-of-magnitude uniformity columns, labelled as estimates: the path error
from an index inhomogeneity of 1e-6 and from a 0.1 K gradient, each over the
double-pass glass path (silica dn/dT ~ 1e-5 /K, CaF2 ~ -1e-5 /K), in waves
at 1 um.

**Pitfalls you will meet.**
- Bounds were written for r = 220 mm.  The meniscus vertex's lower bound
  follows the block radius (20 mm beyond it); its thickness bound is 4 mm
  (the F/1.6 envelope point ended on it).  Do not widen bounds silently: a
  point on a bound is reported, and a SECOND run with that bound scaled by
  the radius is reported beside it.
- Keep the record's merit weights (`ladder_w_dist`, `ladder_w_blur`) so the
  rows compare; beat 4b explains what the distortion weight does.
- The groove period and the focus are re-solved inside the chain at every
  iterate; the band must still span the 9 mm of the detector (the ladder
  checks it -- do not bypass the check at small radii).
- One MATLAB at a time unless you have checked memory; verify a background
  run with a fresh probe before reporting its status.

## 3. Deliverables

1. `dyson5_size.txt`: the table, one row per point, families A-C.
2. `dyson5_size.png`: CRF and ensquared energy against block radius, one
   curve per family, the record's values as lines; a second panel of edged
   mass against radius.
3. `BRIEF_dyson5_size.md`, written for Dave: the answer first (the smallest
   block per family that matches R4, with its mass and thickness), the
   solver-check verdict from step 1, the table, what binds at the small end
   (image, bounds, or clearance), and your recommendation.
4. A short section in `challenges/dyson5/README.md`.
The deck is CC's (Linux): hand over the figure and the numbers, do not edit
`deck_dyson.md`.

## 4. Rules

- `dev-candidate`, commit locally in your clone, **push only on Dave's
  word**; no `--amend`; state SHA + branch when you report.
- No engine change is needed.  If a gate or a point exposes an engine
  defect, report it with a reproducer deck: engine fixes are CC-Linux's lane
  (they need both compilers).
- First report back: section 0's gate counts, step 1's identity and solver
  check.  Then run step 2.
