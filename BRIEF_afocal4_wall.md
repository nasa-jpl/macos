# BRIEF for TO: afocal4 — make the clearance a wall, converge the cleared curve, and write the canonical record

_Tasking: Terminal Opus.  Supervisor: CC.  Timing: START NOW (Dave
2026-08-30: "still a few days to go, and a clean answer here will be
heard by people who will understand it").  NOT demo-blocking — deck
slide 13 ships on the delivered −10° result regardless — but see the
**demo window** note in Task 2: a walled operating point landed by
Sunday evening 08-31 can still improve what the room hears, with fresh
sign-off.  Base: the clearing delivery, `MACOS_res_dev 9902c21` +
`macos 1518cbb`, both LOCAL and unpushed; build on them, do not rebase,
do not push (push is Dave's call)._

_Cold-start reads: `macos/BRIEF_afocal4_clear.md` **including the
DELIVERY LOG at its foot** (the "Open, for Dave" list is this brief's
seed); `challenges/afocal4/clearing/README.md` (the full numbers-first
account — §11 for the leverage-4 architectures); `challenges/afocal4/`
{RESULTS.md, STATUS_S4B.md, STATUS_S4C.md}; the packaging delivery log
at the foot of `macos/BRIEF_r2_packaging.md`; memory
`project_afocal4_rodgers2`.  Precedents that bind here: the S4b
`m3_behind_min` wall + `afocal4_pack_seed` ("**a wall needs a compliant
seed or it is a cage**" — warm-starting lost 4 of 5 trade points), and
the S4c gradient lesson (**a 3e-3 forward difference reads this merit's
gradient 17% low**; the stalls it caused were misread as convergence)._

## Why this slice exists

The clearing delivery proved the collimator interference is structural
(the field-walk ratio law), retired the fold and the station with
measurements, and delivered a −10° extraction tilt that clears the union
gate at +37.82 mm with the wavefront 13.6% better.  Two things keep that
from being the *clean* answer yet:

1. **The clearance is not a wall, and the re-solve spends it.**  At
   −8/−9° the solver walks +23.3/+42.3 mm of margin down to +2.3/+0.7 mm
   because the merit cannot see clearance at all.  So the delivered
   −10° is a point that happens to hold margin, not a point chosen on a
   frontier.
2. **The delivered numbers are budget-limited, not converged** (exitflag
   0 at 427 evaluations, on the forward-difference setting S4c measured
   as bad).  The quoted price (blur 553 µm, breathing 0.816%) may be
   partly solver artifact — in either direction.

The audience for the result includes people who will ask exactly these
questions.  The deliverable is the walled, converged tilt-vs-price
frontier and the record that carries it.

## Task 1 — the union-floor wall in `afocal4_build`, with a compliant seeder

The S4b pattern, exactly: a wall, **never a merit term** (the log-domain
merit doctrine stands untouched).

- Parameter `P.pack.union_min` (mm, declared-body floor; default 0).
  The wall enforces `afocal4_union` floor ≥ `union_min` on the
  **declared body model** (1.15× + 15 mm); print bare-lit-glass beside
  it always, per the gate's own convention.
- **A compliant seeder** so solves start feasible (the
  `afocal4_pack_seed` lesson).  The delivered −10° deck
  (`clearing/afocal4_clear_343mm.in`) is a known-feasible seed; a
  seeder that tilts first and re-poses is presumably the shape.
- **Non-vacuity, asserted in code**: with the wall OFF, the −8/−9°
  re-solves must reproduce the margin-spending behavior (+23.3/+42.3 →
  +2.3/+0.7 mm); with it ON, the same solves hold ≥ `union_min`.
- Cost note: the union gate is sampling-free but not free — if calling
  it inside the solver loop is too slow, the S4b wall precedent (reject
  at evaluation boundary) is acceptable; say which you did and what it
  costs per evaluation.

## Task 2 — the cleared frontier: tilt vs price, walled and converged

The clean answer is a curve, not a point.

- Sweep the extraction tilt (suggest −6…−12° in 1° steps; extend if the
  frontier is still moving at an end) with the wall ON at a real margin
  — run `union_min` = 0 and one positive value (suggest +15 mm, the
  declared allowance's own pad; state your choice) — conics + FM
  standoff + front end re-solved at each point, **central differences**
  (the S4c lesson), converged rather than budget-capped; report
  exitflags and evaluation counts.
- Per point: union floor (declared + bare glass), WFE rung 2 max, blur,
  wander, breathing, max chief AOI, M / exit beam / collimation, rays
  lost.  The same columns the delivery table used, so the two read
  side by side.
- **The operating point** = the minimal pupil price that clears with
  declared margin.  Key question the room will ask: does a walled −8°
  hold real margin at materially less pupil damage than −10°?  If yes,
  that is the number the deck should carry; if no, say so and the −10°
  stands with a converged pedigree.
- **Central-difference polish of the delivered −10° deck** regardless,
  so every quoted number is converged.  If the polish moves a deck-quoted
  number by more than ~1%, flag it to CC immediately (slide 13 carries
  8993 / 553 / 0.82% / +37.8 / 1.24×).
- **Demo window**: the talk is 09-01.  If the walled frontier lands by
  Sunday evening 08-31 with a materially better operating point, CC
  re-cuts slide 13 with Dave's sign-off; after that, the deck ships as
  is and this becomes record work.  Do not rush the convergence to hit
  the window — a budget-capped "better" point is exactly what this slice
  exists to retire.

## Task 3 — the unclaimed-pupil addendum (measure, don't relegislate)

The delivery found a *positive* 4° tilt takes the committed deck's blur
157.0 → 102.6 µm with **no re-solve** — 35% of the pupil metric lying
unclaimed, because a WFE term 130× off target owns the merit's sum of
squares.  Two bounded measurements:

- The signed-tilt blur/wander/breathing curve on the committed deck
  (no re-solve) over ±6° or so — the one-figure statement that the
  committed trade curve was not at its own pupil optimum.
- The same question on the **cleared** operating point: at fixed tilt
  and wall, does a pupil-weighted polish (S4c machinery: `afocal4_score`
  `'zone'`/rim options exist) recover any of the 553 µm without giving
  back WFE or margin?  One or two solves, not a campaign.

Document both as an addendum in the record (Task 5).  The merit doctrine
(log-domain, walls-not-terms) is settled — this task *measures* the
slack, it does not reopen the doctrine.

## Task 4 — leverage 4 stays priced, not built

Only if the walled frontier cannot deliver an acceptable pupil price
does the fifth mirror graduate to a build task — and that would be its
own brief.  Here: keep it priced against the walled frontier's best
point (the delivery's bar: beat blur 553 µm / breathing 0.816% / wander
560 µm at 8993 nm while supplying ≥ 111.6 mm of field-independent
separation).  Update the bar if the frontier moves it.

## Task 5 — the canonical record (Dave's explicit ask)

`clearing/README.md` is the numbers-first account and stays so.  What is
missing is the resolution *in the canonical afocal4 record*, findable by
someone who enters through `challenges/afocal4/` and readable by an
optical designer outside the project:

- A **CLEARING section in `challenges/afocal4/RESULTS.md`** (or a
  `RESULTS_CLEARING.md` beside it — your call, say why): the defect and
  the number that hid it (per-field daylight vs union footprint), the
  field-walk ratio law with the measured terms, the three leverage
  retirements/delivery *as measurements*, the walled frontier from Task
  2 with its operating point, the price table, the packaging
  consequences (deepest 1.81×→1.24×, zero flats, Path A does not close
  — **two ratios, named apart**: deepest/spacing vs overhang/spacing),
  and the standing gate with its non-vacuity statement.
- **Mark the retired remedies retired in place**: the S4b single-fold
  recipe and Path A get a superseded note where they are recorded (do
  not delete them — the gate lesson "a margin is a number, not a body"
  is part of the story), pointing at the clearing record.
- **The §S4b.4 record discrepancy** (folded-demo pupil placement:
  RESULTS.md says z +0.614 m, the committed `.in` measures +1.3782 m;
  `check_record` reproduces).  With the fold recipe retired this is
  historical: suggested handling is a dated correction note in §S4b.4
  stating both numbers and that the recipe is superseded — but it edits
  the existing record, so **flag the exact wording to Dave via CC before
  committing it**.
- Resolved oddities from this slice go into the record **at resolution
  time** (the standing rule), not just the delivery log.

## Ground rules (all established, none new)

- Additive only under `challenges/afocal4/` — nothing overwritten;
  engine untouched; every number engine truth over the field box; no
  text-parsing `.in` files for geometry; `iface` from `zElt`, never
  vertex geometry; an all-failed sweep errors, never returns clean.
- Gates non-vacuous and asserted in code; `tAfocal4` and `tAfocal4Clear`
  stay green; new gates registered in the suite.
- Local commits only, on top of `9902c21`/`1518cbb`; delivery log at the
  foot of THIS brief; nothing pushed; outward-facing text (slides,
  anything leaving the repo) waits on Dave's sign-off via CC.
