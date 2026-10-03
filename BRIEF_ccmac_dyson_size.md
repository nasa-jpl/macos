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

---

# Round 2 (2026-10-02) -- the Dyson WITHOUT the meniscus (Jim's point)

Your round 1 is accepted (the slit is the lever; the envelope's radius axis
is the design's limit; deck slides 25-26 carry it, with 130 mm as the
headline two-module block because it keeps the record's distortion).  One
correction for your report: CaF2's index is LOWER than silica's (1.43 vs
1.45), so "higher index narrows the cone" is not the mechanism; say the
mechanism is not established, or test it (family F below is one test).

**Jim (who builds these), to Dave, 2026-10-02:** "We've never built a Dyson
with a meniscus, because of the losses.  Without an AR coating that works
over 400-2500 nm, the losses are ~0.92^2.  You could compensate by making
the system faster, but that brings along a host of other challenges (mass,
volume, image quality, etc).  The meniscus looks very thin, which is a
challenge for fabrication, mounting, and vibe."

So R4's corrector is a liability on three counts: throughput (four extra
uncoated air-glass crossings: 0.966^4 = 0.87 in silica by Fresnel, Jim's
rule 0.92^2 = 0.85), the part itself (4 mm thick, 155 mm across, and the
solve ends ON the 4 mm bound -- it wants it thinner), and, not yet examined,
ghosts (four near-concentric uncoated surfaces between block and grating).

## The question

Which configurations meet the spec with NO meniscus -- ladder rung R3 (the
de-concentred block with its face conic and h^4/h^6 terms; `'rungs', 3`)?

What is known: R3 at silica 220 mm with the 54 mm slit has CRF 2.10 px and
EE 0.48 (fails).  Size alone needs a 341 mm block (R1, radius free).  The
scaling law (blur ~ h^4.17 / r^3.19, quarter-pixel at r = 213 mm for the
54 mm slit's corner, h = 28.2 mm) predicts ~100 mm for a 27 mm slit's
corner (h = 15.7 mm with the 8 mm dispersion offset).  That prediction is
the thing to test: your round 1 shows the short slit closing to 100 mm WITH
the meniscus; if it also closes without, the meniscus goes.

## The work (extend your own driver; still new files only)

Add a `rung` option to `dyson5_size_trade.m` (3 = R3, 5 = R4) and these
families, by continuation exactly as before.  Seed the first point of each
from the R3 of record (the rung-3 entry of the saved ladder, 220 mm / 54 mm
/ silica), walking the SLIT 54 -> 40 -> 27 mm at 220 mm first, then the
radius.
| family | rung | glass | slit | radius walk |
|---|---|---|---|---|
| D | R3 | silica | 27 mm | 220 -> 60 mm |
| E | R3 | CaF2 | 27 mm | 220 -> 60 mm |
| F | R3 | CaF2 | 54 mm | 300 -> 180 mm (does a one-module, no-meniscus Dyson exist at any size up to 300 mm?) |
| G | R4 | silica | 54 mm, 220 mm | the THICK-meniscus basins of the global search (`dyson5_s3_r4global.txt`, starts 12 and 7, ~30 mm thick), engine-scored: a buildable part, if the meniscus is kept at all |

**Two verdict columns per point:** `matches R4` (round 1's rule) and `meets
the SPEC` (smile and keystone < 0.1 px, CRF < 1.5 px, SRF < 2.1 px, no
variable on a bound).  Without the meniscus R4's ensquared energy may be out
of reach while the specification is met; report EE either way, and name the
smallest radius per family under each rule.

**A throughput column:** the number of air-glass crossings on the slit ->
detector path and the uncoated Fresnel product at 1 um for that glass
(labelled as uncoated, normal incidence).  It is the first entry of the
radiometric chain and the direct answer to Jim.

Keep everything else from round 1: engine scores only, bounds reported not
widened silently, edged volume and mass, the clearance gate, the two
uniformity estimates.

## Deliverables

`dyson5_size.txt` / `.png` extended with families D-G (same figure, rung
shown by line style); a "Round 2: no meniscus" section at the TOP of
`BRIEF_dyson5_size.md` with the answer first (the smallest no-meniscus block
per family under each rule, its mass, thickness and throughput against R4
and against round 1's 130 mm block); the corrected CaF2 sentence; README
section updated.  First report back: family D.  Commit locally; push on
Dave's word.

---

# Round 3 (2026-10-02) -- Jim's comparison, spectrometer side; and "do better than Fresnel"

Jim's reply (2026-10-02) reframes the trade.  The VSWIR instrument already
has TWO imaging spectrometers of 3k pixels each (54 mm slits) behind two
telescopes (30 m GSD at ~550 km).  His question: **is four telescopes +
spectrometers with fused-silica lenses and 1.5k pixels better than two with
bigger CaF2 lenses and 3k pixels?**  On CaF2: large pieces exist (50 cm
boules), but a boule is several crystals, your piece is carved from a
single-crystal volume, and the price grows faster than the volume.  On
losses: "the VSWIR spectrometers have always used uncoated Dyson lenses and
eaten the Fresnel losses.  Your estimate is spot-on.  Maybe you and AI can do
better."

## 3a. The comparison table (your families F and D/E, re-cut as Jim asks)

One table, two columns, every row engine-scored from `dyson5_size.txt`:
| | 2 x (3k, 54 mm slit, CaF2, no meniscus) | 4 x (1.5k, 27 mm slit, silica, no meniscus) |
|---|---|---|
| block (radius = thickness), per module | 240 mm (F) | 130 mm headline / 100 mm smallest (D) |
| edged glass per module, total glass | ... | ... |
| single-crystal CaF2 volume needed per block (the carve) | the circumscribing rod + a margin you state | n/a |
| CRF / SRF / EE / smile / keystone | ... | ... |
| air-glass crossings, uncoated throughput | 4 / 0.88 | 4 / 0.87 |
| gratings, detectors, slits | 2 / 2 x 3k | 4 / 4 x 1.5k |
| grating radius and diameter, length | ... | ... |
Add the silica 3k row (220 mm WITH the 4 mm meniscus, 8 crossings) as the
reference Jim knows, and the CaF2 1.5k row (80 mm) as the small end.  Do not
price CaF2 -- state the volume and let Jim price it.

## 3b. "Maybe you and AI can do better" -- the throughput question

Three routes, each a measured number, none a promise:
1. **Fewer crossings.**  Count what each design actually needs: the slit
   mask DEPOSITED on the block's flat face and the detector's window
   cemented to it would leave only the convex face's two crossings
   (0.93 uncoated).  Score that geometry with the chain + engine: slit at
   z = 0 on the face (no air gap) and the detector plane at the face; report
   CRF / EE / clearance against the record.  R5's sweep says air costs
   ~0.5 px per mm at 54 mm; at 27 mm the margin may buy a window.
2. **A broadband AR coating** on fused silica and on CaF2 over 400-2500 nm:
   `+macos/thinfilm_rt.m` (the engine's Abeles stack, gated by
   tPolRadiometric) can score any stack.  Optimise a 2-4 layer design over
   the band for average reflectance (candidate materials: MgF2, SiO2,
   Al2O3, Ta2O5 -- state the indices you use and their source) and report
   the band-averaged transmission for 4 crossings against uncoated 0.87.  A
   single-layer MgF2 on silica is known to do little (index mismatch); say
   so with the number.
3. **The working distance.**  Jim: "If you and AI give up, I will tell the
   tricks to greater working distances with better response functions."
   Before we ask: scan the slit AND detector standoff at the 27 mm slit on
   the 130 mm silica block (0.5 -> 3 mm, both sides equal, R3 re-solved at
   each) and report CRF / EE / clearance.  If the small module tolerates a
   few mm, that is the working distance; if not, the number tells Jim where
   we are stuck.

Deliverables: `BRIEF_dyson5_jim.md` with 3a's table first, then 3b's three
numbers; records `dyson5_jim_*.txt`; new files only, same rules as before.
First report back: the table (3a), which is a re-cut of what you have.

## Round 3 note (2026-10-03, CC): route 1's crossing count -- amend before Jim sees it
Round 3 is accepted (3a is the reply's core; 3b's AR and standoff numbers
stand).  One correction to route 1.  Depositing the slit on the block's face
removes NO crossing: the telescope's beam arrives in air and enters the glass
at the slit plane whatever carries the mask (it removes the mask's standoff,
which is mechanical, not radiometric).  Cementing the detector WINDOW to the
block removes the window's two crossings, which the 4-crossing baseline
never counted (the dewar window is real: block exit, window in, window out
= 3 crossings on that side, versus 1 at the window's vacuum face when
cemented).  So the honest statement is: with the window counted, uncoated
throughput goes 6 crossings / 0.81 -> 4 / 0.87 by cementing the window;
the convex pair's AR then takes it to ~0.90; "2 crossings / 0.934" is not
reachable with a cold detector behind a window.  Please amend
`BRIEF_dyson5_jim.md` 3b route 1 and the bottom line accordingly (one
commit, by path); CC's draft reply to Jim already uses the corrected count.

---

# Round 4 (2026-10-03) -- the three-mirror anastigmat at Jim's numbers, by the design layer

The telescope side of Jim's comparison is open.  TO's three-mirror family
(a coaxial parent with the stop at M2, solved on axis and then moved 9 deg
off axis through the offset_imager ladder) does not image at F/1.8: CC
measured 6-65 px spots at the slit on its best designs.  TO's two-mirror
modified Schwarzschild images (0.5 px on axis) but, at f = 330 mm, scales
to 2-5 px across a 9.4 deg strip (the review's TMS is a short-focal-length
form).  The review's OWN first telescope example for our regime is
**Mouroulis & Green 2018 Fig. 6: a 420 mm F/1.8 three-mirror anastigmat
with a 16 deg linear field** ("dimensions roughly equal to its focal
length"), and sec. 5.2 names the telecentric variant with the stop on the
secondary.  That is Jim's "reasonably long unobscured telescope".  Your
job: reach that class with DIFFERENT tooling from TO's -- the design
layer's `macos.design.Telescope` (TMA layout from the f-numbers, Seidel
seed, `add_fold`, `realize_apertures`, `view_layout`) and the engine's
native multi-field optimizer (CALIB through `optimize`, which closed the
e2e / e2e2 TMAs at f/1.75 primaries), solving the BIASED configuration from
the start rather than an on-axis parent moved afterwards.

**Spec (per telescope; Jim: 30 m GSD at ~550 km, 18 um pixels):** f 330 mm,
D 183 mm at F/1.8, telecentric and flat at the slit, strip 9.4 x 0.3 deg
(3k module) and 4.7 x 0.3 deg (1.5k module), unobscured, every clearance
positive with the Dyson's bodies at the slit (`spectrometer_clearance`,
the pattern in `dyson5_run` stage t2), the exit pupil matched to the
grating (the Dyson accepts chiefs within 0.09 deg of telecentric).  Read
first: the review's sec. 5.1-5.2 and Fig. 6 (`challenges/dyson5/
Mouroulis&Green2018.pdf`, git-ignored, on disk); `templates/10_telescopes/
tma_widefield/example_tma_widefield.m` and `tma_offaxis/`; memory
`project_fold_extraction` (fold rules, "conics solved AT the bias field" --
exactly the lesson TO's route missed); `BRIEF_dyson5_beat5.md` sec. 5 for
the `seidel_seed` PNP first-order defect (EFL 381 vs 126 mm: check the
seed's EFL by exact trace before trusting it).

**Steps, each engine-scored on the slit's 7 fields (rms spot radius and
energy in an 18 um pixel), each a committed record:**
1. The coaxial TMA parent from `tma_layout` at f 330 / F/1.8 with the
   stop on the secondary, conics solved by the native optimizer AT a field
   bias of 3-6 deg along-track (scan the bias), the 9.4 deg strip as the
   field set; report spot per field at each bias.
2. The unobscured section: `add_fold` / off-axis apertures, every clearance
   positive (your own gate, numbers not bodies), the bias that first
   clears.
3. Telecentricity and flatness at the slit (chief-ray angle per field,
   best-focus z per field), the exit pupil's distance.
4. Freeform refinement if conics + aspheres stall above 1 px
   (`optimize_freeform` / `optimize_aspheres` in the Telescope class), one
   rung, with the asphere-step and derivative fixes of 2026-10-01 on the
   engine (pull first; rebuild).
5. The end-to-end row with each module's own Dyson, engine-only join
   (`dyson5_size_F_r240.in` for 3k, `dyson5_size_D_r130.in` for 1.5k),
   scored by the spectrometer's scorer: smile / keystone / SRF / CRF /
   admitted fraction / clearance -- the TMA entry in Jim's table.
New files only (`dyson5_tma_*`, a `BRIEF_dyson5_tma.md`); TO's `tms_*`
and `t4` files are theirs.  First report: step 1's spot-per-field at the
best bias, with the layout rendered.  Length, M2 / M3 diameters and mirror
mass (state the areal density) in every row.
