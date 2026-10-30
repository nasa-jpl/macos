# BRIEF for TO: rodgers3 adjacent-problem run — demo beat 22b build

_Tasking: Terminal Opus.  Supervisor: CC.  Time-critical: the
Keysight/CodeV demo is ~2026-09-01 — wrapper + timed rehearsal +
fallback pre-gen, then Dave's dry run.  Dave's rulings 2026-08-28 are
RECORDED below (§Rulings) — do not re-open them.  Cold-start reads:
this file; `templates/10_telescopes/offset_imager/README.md` +
`oi_walk.m` + `oi_story.m` headers; `t5_walk/t5_walk_REPORT.md` (the
frontier table + record knobs/timings); memory `project_rodgers3`;
`BRIEF_demo_deck.md` slide 22 (the lead-in framing + beat contract);
`challenges/rodgers3/PACKET.md` §B if frontier context is needed._

## The beat (what the audience sees)

CC — live, on screen — takes a field-box spec from the audience,
states the PREDICTED outcome from the committed walk frontier BEFORE
running, then drives one warm-started continuation solve and shows
map + layout + gates confirming (or honestly refusing).  The claim:
the walk table is COMPILED DESIGN KNOWLEDGE — the driver extends it
to an adjacent problem on demand.  Beat (a) = the product in an
engineer's hands (Dave, solo); this beat = AI drives the tool.

## Rulings (Dave 2026-08-28 — settled)

1. **Default/fallback spec = 12×12°** (between committed steps, so
   visibly not canned).
2. **TO builds** (this brief); CC supervises and drives on demo day.
3. **One audience knob** — box full-width in [5, 15]°, offset fixed
   +22.5°, envelope fixed ×1.65, EPD/F# fixed — and the
   **background-solve timing structure** (solve kicked off at the
   top of the demo, TG beat runs while it works, reveal after).
4. **Refusal path held in RESERVE** — script the wording, gate the
   behavior, do NOT stage a plant or spend live time demonstrating
   it.

## Why the knob and warm starts are safe (grounding — don't re-derive)

- `oi_walk` documents `box_deg` as the ONLY supported continuation
  axis (walking offset re-enters the t4 field-walk infeasibility,
  PACKET 2026-08-21).
- Per-step warm starts are committed: `t5_walk/t5_walk_run.mat`
  (re-saved after each step) at 5/8/11/13/15°.  Any asked width
  warm-starts from the largest committed step BELOW it — never a
  cold start (cold = the documented 595 µm failure; that contrast is
  narration, one sentence).
- Out-of-envelope asks: the screen + F8 sentinel refusal answers
  with one accurate sentence (reserve, per ruling 4).  Audience asks
  for offset/EPD variation get the scripted answer: documented
  infeasible axis / full re-instance = overnight run, offer to send
  results.

## Build

**`oi_demo_step.m`** in `templates/10_telescopes/offset_imager/`:
input = box full-width (deg).  Load `t5_walk_run.mat`, pick warm
start = largest committed step ≤ ask, apply the pinned demo knobs,
run screen → ONE full-freedom S5 solve → close/gates → map + layout
figures + a compact printed verdict block: map max, clear floor,
exit err, gates, and one line vs the frontier prediction (the
bracketing committed rows).  Composition of existing calls ONLY (the
`run_s5_budget` leg_ pattern `oi_walk` reuses) — no new solver
machinery, no engine work.

- Go through the bounded stop-pose path — `oi_close` called DIRECTLY
  at large offset is UNGUARDED (known gotcha, memory
  `project_rodgers3`).
- Unique per-run output tag (timestamp) — the stale-PNG trap.
- Metric quoted with every number, the packet contract: strict RMS
  WFE, sphere on the spot centroid on the frozen FPA, anchored at
  the exit pupil, piston-only removal; headline = dense-map MAX.

**Demo knobs** (starting point — rehearsal timing rules):
gn_iters 30→8, nsolve 5×5→3×3; KEEP model 256 / nGridpts 41
(warm-start fidelity) and the 11×11 dense map for the reveal
(trace-only).  A warm-started step is a POLISH (committed neighbor
start qmeans were 22–48 nm), so expect ~2–4 min — but MEASURE.
Record-knob single step is ~20–25 min (121 min / 5 steps): if best
demo knobs still exceed ~8 min, the background structure carries it;
if quality at demo knobs lands far off the frontier prediction,
raise gn_iters and lean on background time instead — never ship a
bad-looking number to save wall clock.  If NEITHER time nor quality
can be made to work, STOP and write it up for Dave — do not degrade
silently.

## Deliverables

1. `oi_demo_step.m` (+ help block per the template's style).
2. **Timed runs at 7, 12, 14°** at the pinned knobs — wall times +
   map max + floor recorded in a short numbers-first section
   APPENDED TO THIS BRIEF (the tg_psi_dm delivery-log pattern) and
   in the template README.
3. **Fallback pre-gen**: the 12×12° default run saved to
   `offset_imager/demo_adjacent/` — PNGs + verdict text.  Fallback
   ladder on demo day: these outputs → committed walk artifacts
   (any step) → frontier table + endgame slide.
4. **Rehearsal script snippet** (a short .md beside the wrapper):
   the beat's spoken structure — spec relay → frontier prediction →
   one visible call → reveal — plus the scripted refusal wording
   and the offset/EPD deflection lines.  Narration seeds:
   - "You asked for 12°.  The frontier brackets that: 11×11 at
     27.3 nm / 25.1 mm floor, 13×13 at 40.0 nm / 24.6 mm — expect
     low-30s and a floor near the spec knee.  Let's find out."
   - ≥13° asks: the clearance knee IS the story — honest deficit,
     priced (endgame: 47.1 nm @ 30.89 mm in the same envelope).
5. Gates (a small test or an assert block exercised by a batch
   runner — SUITE registration optional given demo timeline, Dave's
   call at review):
   - wrapper at 12° reproduces its own pre-gen verdict numbers;
   - an AT-a-committed-step ask (11°) warm-starts from the step
     BELOW and lands within the committed row's neighborhood;
   - out-of-range ask (20°) refuses via the screen with the one
     accurate sentence (the reserve path, gated even though not
     demonstrated live);
   - verdict block's frontier-bracket line quotes the correct
     committed rows.
6. `BRIEF_demo_deck.md` prep-queue 22b row updated with measured
   times; memory update (extend `project_rodgers3` or a small new
   entry linking [[project_keysight_demo]]).

## Traps (paid for once already)

- matlab -batch: script files + exit(0) in the batch RUNNER only —
  never in the wrapper itself; one model size per MATLAB process;
  MACOS_HOME=/home/dcr/dev/macos/macos_f90.
- Read-tool PNGs are cache-stale on overwritten paths — verify
  figures by printed numbers, unique filenames.
- `pgrep/pkill -f` self-match; bound every sentinel loop.
- Sandbox does not block the MATLAB license; run sandboxed +
  background (memory `reference_matlab_test_runs`).

## Standing rules

Commits LOCAL on res dev-candidate; Dave orders pushes.  NO engine
work.  The demo-day driving is CC's job — this task ships the
wrapper, timings, fallback and script, not the performance.

---

# DELIVERY LOG — 2026-08-28 (TO)

## Numbers first

Three widths, run end to end at the pinned knobs.  Every number is the
packet metric: strict RMS WFE, sphere on the spot centroid on the frozen
FPA, anchored at the exit pupil, piston-only removal; headline = 11×11
dense-map MAXIMUM; solve set 5×5 ≠ scoring set.

| ask | warm start | predicted | measured | floor pred → meas | exit err | verdict | wall |
|---|---|---|---|---|---|---|---|
| 7×7° | step 1 (5°) | 18.0 nm | **20.0 nm** (1.11×) | 77.6 → 93.8 mm | 0.002° | PASS | 15.2 min |
| 12×12° | step 3 (11°) | 33.7 nm | **33.6 nm** (1.00×) | 24.8 → 24.9 mm | 0.012° | PASS | 15.2 min |
| 14×14° | step 4 (13°) | 54.9 nm | **51.2 nm** (0.93×) | 21.2 → 24.9 mm | 0.001° | PASS | 15.2 min |

All three PASS both gates.  A solo 12° rerun through the wrapper's own
defaults reproduces the bundle to every printed digit (33.6240 nm /
24.9496 mm / 0.012103°) in 14.87 min — the run is deterministic, and
concurrency is nearly free (three-abreast costs the same wall clock as
one, which is how the bundle is regenerated).

**Pinned knobs: `gn_iters = 1`, `nsolve = 5`, `solve_sampling = 11`**,
model 256, nGridpts 41, 11×11 map.  Wrapper and runner defaults match.

## The knob study (why those knobs, and not the brief's starting point)

12° ask, same warm start, all rows scored identically — they differ only
in solve knobs:

| gn_iters | nsolve | solve_sampling | map max | floor | exit err | wall |
|---|---|---|---|---|---|---|
| 2 | 3 | 21 | 40.96 nm | 24.65 mm | 0.0141° | 13.1 min |
| 8 | 3 | 21 | 38.53 nm | 24.88 mm | 0.0013° | 56.6 min |
| 2 | 5 | 21 | 35.23 nm | 23.89 mm | 0.0427° | 29.5 min |
| **2** | **5** | **11** | **33.55 nm** | 24.94 mm | 0.0005° | 29.5 min |
| 3 | 4 | 11 | 103.85 nm | 18.87 mm | **0.6981°** | 31.5 min |

- **The brief's starting knobs (gn_iters 30→8, nsolve 5×5→3×3) do not
  meet the quality bar, and I did not ship them.**  At `nsolve = 3` the
  12° ask lands at 41.0 nm — *worse than the committed 13° row (40.0)* —
  which on stage reads as the driver doing worse at a smaller box.  The
  brief's suggested remedy (raise gn_iters) is the wrong lever: 2→8
  iterations bought 6% for +43 min, while 9→25 solve fields bought 14%
  for +16 min.  **Spend on fields, not iterations.**  This is the
  recorded S5 solve-field lesson reproducing, with a clean tell in the
  map: the box CORNERS (which a 3×3 grid samples) came out *best* at
  15–22 nm, while the unsampled ridge at YAN 25–27° ran 38–41 nm.
- **`nsolve` must be ODD** — new, and it cost a run to find.  `oi_solve`
  imposes the exit-direction equality on the solve field nearest the box
  centre; an even grid has no centre field, so the exit chief is pinned
  off centre.  Measured at `nsolve = 4`: exit error **0.6981°** (1400×
  the odd-grid rows) with the dense map at 103.9 nm.  Now guarded by a
  warning in the wrapper.
- **`solve_sampling` is NOT a cost lever** — the solve is deck
  write/parse bound, not ray bound: 11 vs 21 nGridpts measured 29.51 vs
  29.49 min.  It is pinned at 11 on *quality* (better on all three
  reported axes at 12°), which is the opposite of why I reached for it.
- `nsolve = 5` at `gn_iters = 8` was abandoned at 512 s, still inside its
  first Jacobian — time-infeasible; recorded rather than left implicit.

## Two brief assumptions that measurement overturned

1. **The ≥13° "honest deficit" narration is wrong — a wide ask PASSES.**
   The brief scripts a wide ask as a priced clearance deficit.  At 14°
   the frontier LINE predicts a 21.2 mm floor (below spec) and the
   re-solve delivers **24.9 mm with both gates PASS**, at 51.2 nm against
   a predicted 54.9.  Same mechanism the endgame found at 15°: the walk's
   widest rows stopped early on the WFE-only plateau break while the
   clearance hinge was still pulling, so the committed line *understates*
   what the envelope holds — now shown to hold in the frontier's interior
   as well.  `REHEARSAL.md` scripts the corrected line.
2. **The timing structure does not fit, and that is Dave's call.**  A
   step at the quality bar is **~15 min**, and no knob buys it back
   without either dropping below `nsolve = 5` (the bad look above) or
   cutting design freedom (which the brief forbids).  This is a
   *scheduling* change, not a degradation: take the spec at the top of
   the Rodgers section (slide ~16) and reveal at 22 — about 15 slides of
   cover — rather than launching "at the top of the demo".  If the spec
   can only be taken at slide 22, run the beat off the pre-generated
   bundle and say so.  Flagging rather than degrading silently, per the
   brief's own instruction.

## A third measured item for the demo (narration, not a defect)

**An ask landing exactly on a committed width is the worst case, and CC
should get out in front of it.**  There the "prediction" is that
committed row itself — solved with 30 GN iterations on a 5×5 grid —
while the live run takes ONE polish step off the width below.  Measured
at 11°: **39.1 nm against the committed 27.3 (1.43×)**, both gates still
PASS, floor 57.9 mm.  Between committed widths there is no such
asymmetry (12° lands 1.00×), which is precisely why the ruled default
spec is 12.  `REHEARSAL.md` carries the pre-emptive line.

## Delivered (all LOCAL on res `dev-candidate`, uncommitted)

- `templates/10_telescopes/offset_imager/oi_demo_step.m` — the wrapper.
  Warm start = largest committed step STRICTLY below the ask (so an ask
  AT a committed width is a real continuation step, not a re-score of an
  already-solved design); frontier prediction printed BEFORE the solve;
  one full-freedom S5 step through the bounded stop-pose path; verdict
  block printed and written beside the figures; unique timestamped tag by
  default (the stale-PNG trap).
- `.../run_oi_demo.m` — batch runner for the background solve (`exit(0)`
  lives here, never in the wrapper).
- `.../demo_adjacent/` — the fallback bundle at **7, 12 and 14°**
  (map/layout/fields PNGs, verdict text, deck, run.mat) plus
  `REHEARSAL.md`: spoken structure, the two scripted refusals, the
  offset/EPD deflections, and the fallback ladder.
- `tests/tOiDemoStep.m` — 11 gates, NOT suite-registered (two full solves
  ≈ 30 min); registration is Dave's call at review.
- README: run line, an "adjacent problem" section, the knob study, and an
  architecture entry.

## Gates

The brief's four, plus one it did not ask for:

- 12° fresh run reproduces the pre-gen verdict numbers (RelTol 1e-6 —
  equality with numerical noise, not a band, justified by the measured
  determinism above);
- an AT-a-committed-step ask (11°) warm-starts from the step BELOW (8°)
  and lands in the committed row's neighbourhood;
- an out-of-range ask (20°) refuses via the screen, in one sentence,
  before any solving, producing no map;
- the verdict block quotes the correct bracketing committed rows, and an
  at-a-step ask brackets onto itself;
- **added: the F8 traceability refusal is actually TRIGGERED** — and
  getting there corrected a wrong assumption of mine, caught by the gate
  itself (it failed on the first run, 10/11).  I had assumed a very wide
  box would trip the no-rays sentinel.  **It does not.  Width alone can
  never trip it on this instance**: the carried design keeps tracing all
  the way out, degrading smoothly rather than losing rays —
  30° → 1765 nm, 60° → 1.25e5, 90° → 2.59e5, 120° → 4.55e6,
  150° → 5.71e6, never the 1e9 sentinel.  So the original gate exercised
  nothing and the wrapper simply fell through to a (slow) real solve.
  Replaced with a genuine trigger via the documented `run_mat` knob — a
  warm start carrying an absurd 50 mm M1 radius, where every ray misses:
  screen returns 1e9, refusal fires, no map produced, and it is fast
  because the screen answers before any solving.  A companion gate pins
  the corollary: **a wide ask is stopped by the RANGE screen, not by
  traceability**, so nobody later "simplifies" the range check on the
  assumption that F8 would catch it.

  Consequence worth knowing for the demo: **the range screen is the only
  operative guard.**  The F8 path is defensive code inherited from
  `oi_walk` (where it earns its keep — it is the t5 *unguided* failure
  mode, cold starts losing 104/121 fields, not a continuation one; the
  committed walk itself took zero halvings).

## Wording bugs caught before they reached a slide

- The block printed "spec >= 25 … PASS" beside a 24.65 mm floor.
  `oi_gates` passes at `min(clear_m)` less a 1.5 mm hinge knee; the block
  now quotes the threshold actually applied and flags WARN.
- The block claimed "the 25 mm clearance hinge is live at this width"
  beside the 7° run's **93.8 mm** floor.  The hinge is dormant there and
  the extra floor is not "bought"; the claim is now conditioned on the
  floor being inside the WARN band.

## Not done / open

- The `.pptx` slide 22b assets — per the standing rule this ships the
  wrapper, timings, fallback and script, not the deck.
- Suite registration of `tOiDemoStep` — Dave's call; runtime is why it is
  out.
- The timing-structure ruling above.
