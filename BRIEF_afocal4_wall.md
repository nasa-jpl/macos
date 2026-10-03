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

---

# DELIVERY LOG (2026-08-31)

Everything new lives in `MACOS_res_dev/mmacos/challenges/afocal4/wall/` (its
README is numbers-first).  **Nothing under `challenges/afocal4/` was
overwritten**; the additive clauses beside it are listed in §5 below, all
default-OFF and all asserted to leave the committed decks rebuilding byte for
byte.  Local, unpushed.  No engine work.  Base: `MACOS_res_dev 9902c21` /
`macos 1518cbb`, built on, not rebased.

## 1 — the wall, and it is a wall

`afocal4_union`'s floor on the **declared body model** is now a wall in
`afocal4_build` beside the S4b `m3_behind_min` one, via a new
`../afocal4_union_wall.m`.  Four `P.pack` fields: `union_enforce` (default
**false**), `union_min` (0 m), `union_body_k` (1.15), `union_body_pad`
(0.015 m).

* **Default OFF is load-bearing.**  With the wall on, `afocal4_build` cannot
  re-emit the committed 343 mm deck (−79.89 mm) and every S4 / S4b / S4c /
  clearing artifact stops reproducing.  Same reason `P.pack.enforce = false`
  keeps the unbuildable S4 reference alive.
* **`clear_build` DEFERS it past the tilt.**  `afocal4_build` emits the
  *untilted* train; a wall applied there judges the design the tilt exists to
  get away from and rejects every iterate — a cage.  Gated on one `P`, both
  halves.
* **Cost, measured, inside the loop** (the brief asked which and what):
  evaluated INSIDE the build, so every iterate the solver sees is compliant.
  Cold, at solve sampling: `clear_build` 1.61 s + `afocal4_score` 6.56 s =
  8.17 s, wall +4.18 s → **+51 %**.  Warm in a solver loop the whole
  evaluation is ~3.2 s.  Nearly all the wall's cost is the nine-field
  re-trace `afocal4_score` already paid for; the probe count is nearly free
  (314 probes 4.18 s vs 65 probes 3.52 s) so it is **not** tuned down — a wall
  must hold the quantity the gate reports.
* **Sampling, stated not hidden:** the wall is judged at SOLVE sampling and
  every table quotes REPORTING sampling.  A bigger ray grid makes a bigger
  union hull, so the wall's number is the optimistic one — measured **+1.0 to
  +3.2 mm** across the frontier and carried per point as
  `R.sampling_bias_mm`.

## 2 — the compliant seeder, and a finding about the law

`wall_seed`, `afocal4_pack_seed`'s sibling.  Preferences: the tilt alone
(nothing moved — which is what keeps a point comparable with the delivered
row); the smallest field-mirror standoff change that clears, on the parent's
own front end; a different M2 radius; last, the delivered −10° design's DOFs
(flagged as a different basin).  Margin 10 mm over `union_min`.

**THE FIELD-WALK LAW IS 5–6× OPTIMISTIC AS A PREDICTOR OF A STANDOFF CHANGE.**
The obvious cheap ranking is `f̂ = f + 2|α|(d − d₀)`, pure closure arithmetic.
Measured at −6°, moving the standoff −38.6 → +250 mm takes `d` 0.563 → 0.680 m
and the law predicts the floor going −13.00 → **+11.45 mm**; it measures
**−8.25**.  Realised sensitivity **33 mm per metre of `d`** against the law's
**209**.  The law is not wrong — the tilt really does supply a
field-independent `2αd` — but a **standoff** change moves the
field-PROPORTIONAL part at the same time and the two nearly cancel.  *Leverage
2 showing up inside leverage 3.*  So the seeder **probes and bisects on
measurements** instead: 4–6 gate evaluations, no prediction relied on, and
`INFO.slope_mm_per_m` carries the realised slope beside the law's.

Ranking by `afocal4_pack_seed`'s own "weakest field mirror first" is worse
still — at −6° it spent ten gate evaluations walking a −13.0 mm floor down to
**−83.2 mm**.

**And `P.parent` is not the parent design.**  It carries Mike's raw secondary
(R_M2 468.8 mm, t_M1M2 1.0492 m) while the committed 343 mm deck has a
re-solved front end (448.4 mm, 1.0420 m).  Filtering seed candidates through
`P.parent` admitted 21 standoffs of 57 and **not one was the deck's own**;
carrying `D.R2`/`D.t1` into the closure — what `afocal4_build` does — admits
54, spanning `d` = 0.255…0.821 m against the parent's 0.563.

## 3 — the bug that cost an hour, and the rule it earned

**A wall belongs on ITERATES, never on the REPORT that follows them.**
`clear_solve` built its final, quotable deck through the same walled builder
the objective used.  Inside the objective a violation becomes a large finite
residual and the solver backs out of it; in the report path there is nobody to
back out — and because the wall is judged at SOLVE sampling while the report
builds at REPORTING sampling (~2.5 mm lower floor), **a converged design
sitting on its wall throws out of its own report and takes the whole solve
with it.**  Measured twice: `t-80_u15` and `t-90_u15`, both re-run.  Fixed in
`clear_solve` (the report build measures, it does not judge) and guarded in
`wall_point` (a round that throws costs a round, not a night).

Operational note worth keeping: **a running MATLAB caches the function it has
already called**, so a mid-flight fix does not reach processes already in
their solve.  Only points the fleet started AFTER the edit picked it up; the
two that did not were identified from the fleet log's start times and re-run.

## 4 — the addendum: the unclaimed pupil, and a correction to the clearing record

The signed-tilt curve on the committed deck, 1° steps, nothing re-solved
(`wall/afocal4_wall_unclaimed.png`):

* the blur optimum is **101.3 µm at +5°**, not 102.6 at +4° — **35.5 %** below
  the committed design's, for a tilt and no re-solve.  The clearing stage's 2°
  grid straddled it.
* **it is not free.**  What pays is magnification **breathing**, 0.124 →
  0.795 %, and breathing is the one pupil target the committed design *meets*
  (0.4 %).  Wander tracks blur; the wavefront is flat to 0.5 % across ±8°.
* **the free move is a small NEGATIVE tilt:** at −1° breathing reaches
  **0.0381 %**, three times better than committed and ten times inside target,
  for 178.2 vs 157.0 µm of blur and no loss of clearance.

**AND THE REASON THE CLEARING RECORD GIVES IS WRONG.**  It attributes the
unclaimed pupil to a merit owned by a wavefront term 130× off target.  Read
off the residual vector, `afocal4_score` divides the per-field wavefront
residuals by `sqrt(K)`, so the wavefront block is **78 %** of the merit, not
~97 %, and **the merit PREFERS the +5° point** (29.40 against 30.22).  The
committed design is not at its pupil optimum because **an extraction tilt was
never in the DOF set that produced it** (`{conic, standoff, front}`).
Corrected in place in `clearing/README.md` §12.3 and written up in
`RESULTS.md` § C.7a.

**The rung-4 rigid-body tilt is a DIFFERENT operation, with the opposite
sign** — measured, because "the freedom was already there" had to be tested
rather than assumed.  `P.bounds.tilt` allows ±2.86° of FM x-tilt about the
element's **vertex** with the train not re-posed:

| operation on the field mirror | WFE (nm) | blur (µm) | breathing (%) | floor (mm) |
|---|---|---|---|---|
| committed, nothing moved | 10407.0 | 157.0 | 0.1240 | −79.89 |
| rung-4 rigid body +2.86° (vertex, the bound) | **9818.8** | **613.0** | 0.1381 | −86.19 |
| extraction tilt +3.00° (chief, train re-posed) | 10436.9 | **116.7** | 0.5271 | −84.27 |

The rigid-body tilt is a *wavefront* knob that spends pupil; the extraction
tilt is a *pupil* knob that spends a little wavefront.  So the extraction tilt
is a genuinely new degree of freedom — and a second unclaimed quantity sits
beside it: **a rung-4 pass on the committed 343 mm deck is worth ~5.6 % of
wavefront**, on a design whose rung-4 DOFs were never solved.  Neither chased;
both recorded.

## 5 — what changed outside `wall/`, all additive and default-off

| file | change |
|---|---|
| `afocal4_union_wall.m` | **new** |
| `afocal4_build.m` | one deferred-wall clause (§3c) + `'defer_union'` option |
| `clearing/clear_build.m` | one post-tilt wall call; `'defer_union',true` on its inner build |
| `clearing/clear_solve.m` | the final report build measures, does not judge (§3) |
| `afocal4_params.m` | four `P.pack.union_*` fields, wall OFF |
| `tests/tAfocal4Wall.m` | **new**, 8 tests, registered in `SUITE_FREEFORM` |
| `RESULTS.md` | **new `# CLEARING RESULTS` section** (§6) + superseded note on §S4b.4 |
| `README.md` | the packaging → clearing → wall arc made findable from the entry point |
| `packaging/README.md` | Path A marked superseded in place |
| `clearing/README.md` | §12.3 corrected in place |

**Gates: `tAfocal4Wall` 8/8 (new), `tAfocal4` 8/8, `tAfocal4Clear` 8/8.**  And
the delivered cleared deck re-scores to its recorded numbers exactly —
8992.68 / 553.34 / 0.8160 / 559.87 / +37.82, zero rays lost.

## 6 — the canonical record (Task 5)

**Placed as a `# CLEARING RESULTS` section appended to
`challenges/afocal4/RESULTS.md`**, not as a separate `RESULTS_CLEARING.md`.
Why: RESULTS.md is what someone entering through `challenges/afocal4/` reads,
S4b and S4c were both appended as sections rather than split out, and the
superseded notes elsewhere need one findable place to point at.  The two
stage READMEs stay the numbers-first accounts.

It carries: the defect and the number that hid it (per-field daylight vs the
union footprint); the field-walk ratio law with its measured terms and the
`M × iface` pin; the three leverage retirements/delivery as measurements; the
price table; the wall and the seeder; the walled frontier and its operating
point; the packaging consequences with **the two ratios named apart**
(deepest/spacing 1.81× → 1.24× vs overhang/spacing 0.81× → 0.24×, and Path A
does not close on the cleared deck — 96 routes, 15 admissible, every one lost
rays); the standing gate with a **non-vacuity table** naming where each half
is asserted; the addendum; leverage 4's bar; and nine new earned rules
(23–31).

**Retired in place, not deleted:** `packaging/README.md` §3 (Path A) and
`RESULTS.md` §S4b.4 (the S4b single fold) both carry dated superseded notes
that keep the measurements and the "a margin is a number, not a body" lesson.

## 7 — FOR DAVE: the §S4b.4 correction wording, held for sign-off

`packaging/check_record` reproduces the discrepancy.  Measured from the
committed `afocal4_b_final_folded.in`, against what §S4b.4 and
`STATUS_S4B.md` both state:

| quantity | as recorded | as measured |
|---|---|---|
| interface pupil (m) | [+0.304, −0.004, +0.614] | **[+0.2483, −0.0051, +1.3782]** |
| 1000 mm envelope ends (m) | [+1.304, −0.017, +0.614] | **[+1.1990, −0.3151, +1.3782]** |
| instrument z-slab (m) | +0.464 … +0.764 | **+1.2282 … +1.5282** |

On the committed deck the fold and the interface plane are both at
z = +1.3782 m.  **The conclusion those numbers support is unchanged** — every
z is positive, so the instrument sits entirely behind the primary — and the
fold-is-null row, the WFE/blur/breathing/wander table and the 137 mm
closest-approach are measurements of different quantities and are unaffected.

Proposed note, to be inserted at the foot of §S4b.4 and mirrored in
`STATUS_S4B.md`, **not committed until you approve the wording**:

> **DATED CORRECTION (2026-08-31).**  The interface-pupil coordinates and the
> instrument z-slab quoted above for the FOLDED deck do not reproduce from the
> committed file.  Measured by `packaging/check_record` on
> `afocal4_b_final_folded.in`: interface pupil **[+0.2483, −0.0051,
> +1.3782] m**, envelope ending **[+1.1990, −0.3151, +1.3782] m**, z-slab
> **+1.2282 … +1.5282 m**; the fold and the interface plane are both at
> z = +1.3782 m.  The statement the numbers support is unaffected — every z is
> positive, so the instrument is entirely behind the primary — and the
> fold-is-null result is a measurement of a different quantity.  The recipe
> itself is superseded (see the head of this section), so this corrects the
> historical record, not a live design.


**SIGNED OFF (Dave 2026-08-31, via CC): the §S4b.4 correction is
approved as written.**  Commit the note into `RESULTS.md` §S4b.4 and
mirror it in `STATUS_S4B.md` with the handback.
## 8 — THE FRONTIER: both halves of the brief's premise reversed

### 8a — the margin is not spent, and the wall is insurance not the answer

The brief asks for the clearing stage's margin-spending to be reproduced with
the wall off and abolished with it on.  **It does not reproduce.**

| tilt, 0 mm floor | raw tilt | clearing stage (427 ev, forward) | wall OFF, converged | wall ON, converged |
|---|---|---|---|---|
| −8° | +23.34 | **+2.32** | **+28.05** | +28.05 *(identical)* |
| −9° | +42.25 | **+0.69** | **+38.67** | +43.18 |
| −10° | +57.44 | +37.82 | **+45.07** | +45.07 *(identical)* |

Same tilt, same DOFs, same seed, wall off; the only difference is 427
budget-capped forward-difference evaluations against 1209 central-difference
ones.  At −8° and −10° the wall-on and wall-off runs are **identical to the
last digit of round-1 merit** (46.181908, 50.097881) — the wall never rejected
an iterate.  **The margin-spending was a stalled solve on the gradient S4c
had already measured as 17 % low, not a merit blind to clearance.**

The wall is still right to have — nothing else holds the clearance, it refuses
the committed deck at −79.89 mm, and `union_min` is a real threshold — but on
this design it is **insurance**, and **convergence** is what changed the
answer.  Where it DOES bind (+15 mm floor; −9° at 0 mm) it changed the path
and landed a **better** design both times (merit 34.34 vs 36.89; 31.47 vs
32.89), which is not what a cage does.

Determinism checked: the −10° wall-off point ran twice in independent
processes and reproduced round merits 50.097881 / 47.094653 / 35.998867 both
times.

### 8b — the free-standoff sweep is not a tilt-vs-price curve

Sorted by the standoff each solve reached rather than by tilt:

| point | converged s_FM (mm) | K_FM | WFE (nm) | blur (µm) |
|---|---|---|---|---|
| −6° | −229.2 | −3.00 | 12076.2 | 535.8 |
| −8° | +229.7 | −8.25 | 7813.1 | 347.5 |
| −10° (wall off) | +275.9 | −15.78 | 6744.1 | 352.5 |
| −9° | +438.8 | −18.53 | 5288.8 | 288.1 |
| −7° | +535.6 | −25.34 | 3212.5 | 227.8 |

Monotone in the standoff, from a −38.6 mm parent toward a +600 mm bound, with
every solve still descending 24–43 % per round.  **The differences between
tilts are mostly how far each solve walked the standoff.**  Two consequences:
a real tilt curve must hold the standoff fixed (§8c), and **the committed
design is nowhere near its own optimum in the standoff DOF either** — worth
3–4× of wavefront error, where the tilt is worth a few percent.  That is a
second unclaimed quantity beside §4's, and a larger one.

### 8c — the frontier, tilt isolated: THE DELIVERED POINT IS PAST THE KNEE

Standoff pinned at +276 mm, DOFs `{conic, front}`, wall ON at 0 mm, 1628
evaluations over 4 rounds each (round-4 gains 5.4e-2 / 8.3e-4 / 2.1e-2 /
7.7e-3):

| tilt | floor (mm) | bare (mm) | WFE (nm) | blur (µm) | breathing (%) | wander (µm) | AOI | M |
|---|---|---|---|---|---|---|---|---|
| **−8°** | +15.18 | +39.38 | **6513.9** | **279.9** | **0.7210** | **284.2** | 10.68 | 30.0150 |
| **−9°** | **+45.44** | +63.33 | 7794.7 | 352.0 | 0.9190 | 356.4 | 11.01 | 30.0150 |
| −10° | +48.54 | +64.85 | 7682.3 | 456.8 | 1.1912 | 460.9 | 11.29 | 29.9846 |
| −11° | +47.17 | +63.52 | 7464.9 | 482.9 | 1.2044 | 487.0 | 11.34 | 29.9848 |

**The clearance saturates by −9°** (+45.4 / +48.5 / +47.2, flat to ±3 mm) while
the pupil price keeps climbing — so **−10° and −11° are DOMINATED** and the
delivered design sits past the knee.

**The operating point is −9°**, beating the delivered −10° row on four of five
columns: WFE **−13.3 %**, blur **−36.4 %**, wander **−36.3 %**, floor
**+7.62 mm**, breathing +12.6 %.

**And the brief's own question, answered:** *does a walled −8° hold real margin
at materially less pupil damage than −10°?*  **Yes** — **+15.18 mm**, i.e.
exactly the declared allowance's own 15 mm pad, at **49.4 % less blur**, with
breathing and wavefront also better.

Isolated, the tilt's price is small and one-sided: the wavefront moves in a
15 % band with no trend (the clearing stage's "the wavefront is not the price"
survives) while blur, wander, breathing and AOI all grow monotonically with
|tilt|.  **Buy exactly as much tilt as the clearance needs and not one degree
more.**

*Caveat: this curve holds the standoff at one station.  §8b shows the standoff
is worth far more than the tilt, so this is the tilt's price AT a good
station, not the design's optimum.*

## 9 — a second wall defect, found by the addendum

**A WALL IS ONLY A WALL WHILE IT DOMINATES THE MERIT'S OWN SCALE.**
`clear_solve` rejects a wall-violating iterate with a constant residual of 20
per component — merit 5600.  At the study's weights a sound design scores ~30,
so 5600 is an impassable barrier.  Multiply the pupil weights by 16 to measure
§4's slack and a *sound* design scores ~4e4: the same constant now looks
**attractive**, and the solver walks through the wall on purpose.  Measured:
the ×16 run returned a converged `x` whose closure put M3 **1051 mm in FRONT
of the primary** and died in the report build on the S4b packaging wall.

Fixed: the residual scales with the largest merit weight in play,
`20 * max(1, max(P.weights))`.  With the study's own weights that is exactly
20 — **bit-identical to every committed clearing-stage solve**.

**This matters directly for the descent brief**, whose Task 1 introduces a
power-economy regularizer: changing the merit's scale is precisely the
condition under which a fixed-magnitude wall residual stops being a wall.

## 10 — status at handback

**Complete and in the record:** Tasks 1 (the wall + seeder + non-vacuity + the
measured cost), 2 (the frontier, both series, the operating point, the polish
flag), 3's first half (the signed-tilt curve and the corrected explanation), 4
(leverage 4's bar restated against the frontier's own best point), and 5 (the
canonical record, the retirements in place, the §S4b.4 wording held for Dave).

**Task 3's second half, landed after the first commit (§ C.7b):** the
pupil-weighted polish does NOT recover the blur, and the incumbent beats both
re-weighted solves ON THEIR OWN MERITS -- 229.9 vs 301.3 at x4, 3386.5 vs
5886.4 at x16.  The x16 run PLATEAUED there (round-2 gain 1.17e-5); the x4 run
was still descending (round-2 gain 89 %) and is reported as a PROBE, not a
converged point, because quoting a still-descending number as a result is the
exact failure this slice exists to retire.  Both re-weighted solves also broke
the customer interface -- M 0.515 % and 0.98 % off 30 against a 0.1 % target --
because nothing in a re-weighting protects a requirement that was previously
satisfied incidentally.  Earned rule 34.

**Gates at handback:** `tAfocal4Wall` 8/8 (new), `tAfocal4` 8/8,
`tAfocal4Clear` 8/8.  The delivered cleared deck re-scores to its recorded
numbers exactly.  `afocal4_wall()` runs end to end and writes
`afocal4_wall_frontier.png` + `afocal4_wall.mat`.

**Run inventory:** 17 points launched, 13 carry checkpoints and are reported;
the dropped four (−11/−12 past saturation; −9/−10 at the +15 mm wall where it
never binds) and the two re-runs are accounted for in `wall/README.md` §3.
Decks with no checkpoint behind them were deleted — a mid-solve snapshot on
disk is a number waiting to be mis-quoted.

**Not done, deliberately, and named:** no fold-route or envelope work (the
packaging round owns it); no fifth-mirror architecture beyond the bar; the
standoff's own operating point — which § C.4c measures as worth 3–4× of
wavefront, far more than the tilt — is left open and is the largest loose
thread this slice found.  It is the natural first question for the descent.
