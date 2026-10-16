# BRIEF for TO: afocal4 descent — start at seven mirrors, walk back toward buildability

_Tasking: Terminal Opus.  Supervisor: CC.  Timing: START after the wall
slice hands back (TO estimate ~1.5 h from 2026-08-31 am) — its machinery
(union wall, `wall_seed`, walls-on-iterates-never-reports, per-point
sampling bias) is this stage's foundation.  **Time box: 48 hours from
start (Dave 2026-08-31).**  If the box runs out mid-ladder, hand back
the completed rungs — every rung is a finished, checked design (the walk
doctrine's virtue); do NOT rush convergence to finish the ladder (the
budget-cap lesson, learned twice).  NOT demo-blocking: the talk is
09-01, nothing from this slice reaches the deck without Dave's sign-off
via CC; if the top solve or early rungs land before the talk they are
Q&A pocket only ("we're walking down from seven to find out").  Base:
the wall slice's handback commits on `MACOS_res_dev` + `macos`; build on
them, do not rebase, do not push (push is Dave's call)._

_Cold-start reads: `macos/BRIEF_afocal4_wall.md` **including its
DELIVERY LOG** (the wall + seeder, the premise reversal, the
iterates-not-reports rule, the field-walk law's 5–6× optimism as a
STANDOFF predictor, the DOF-set finding, the `P.parent` trap);
`macos/BRIEF_afocal4_clear.md` + delivery log;
`challenges/afocal4/RESULTS.md` end to end — S4b/S4c sections, the
`# CLEARING RESULTS` section, and every earned rule;
`challenges/afocal4/{clearing,wall}/README.md`; `PLAN_AFOCAL4.md`;
memory `project_afocal4_rodgers2`.  Precedent for the method: the
rodgers3 continuation walk (`challenges/rodgers3` + the offset_imager
walk records) — solve the easy end, step toward the hard one
warm-started, every rung finished and gated._

## Dave's rulings (2026-08-31) — the frame, not open questions

1. **Start at 7 or 8 mirrors.**  CC's resolution: start at **7**;
   escalate to 8 only if the 7-mirror top solve cannot meet the full
   requirement set with margin — and that escalation, if it happens, is
   itself a ladder datum, not a setback.
2. **The 71 nm wavefront target JOINS the requirement set** at the top.
   The 4-mirror family never met it (best ~7.5 µm at the operating
   points); the descent starts where it IS met and finds where it
   breaks.
3. **Removal policy is CC's call** — stated in Task 2: predictions rank,
   measurements decide.
4. **A retained flat is NOT a mirror.**  The ladder indexes POWERED
   mirrors only.  Flats kept as folds are listed per rung (count +
   role) but do not count; the full fold-route/envelope packaging
   treatment is a FOLLOW-UP round, not this slice.  The buildability
   WALLS still run at every rung — deferring the packaging round does
   not reopen the S4b hole.

## Why this stage exists

The S4 arc asked "what does the fourth mirror buy" and found a
requirement pair with only half met: pupil control was purchased with
wavefront, and 71 nm was never on the table.  The wall slice then found
the committed design missed even its own pupil optimum because **the
extraction tilt was never in the DOF set** — the signature pathology of
building up from too little freedom.  The descent inverts the whole
question: start where EVERYTHING is met with margin, remove one powered
mirror at a time, and measure what each removal costs.  The ladder
answers, from above, the question the arc has circled for a month: **how
many mirrors does the full requirement set actually need?**  It also
supersedes the leverage-4 "priced, not built" bar — the ladder measures
the fifth mirror's worth directly at N=5, on the same footing as every
other rung.

## Task 0 — the N-mirror closure (`afocal4_build` generalized)

The real new machinery.  Pin the invariants, choose the parameterization
yourself:

- The interface identities stay EXACT closures, not targets: M = 30×,
  collimation, pupil at the interface station (the 4-mirror build holds
  them at 1e-9 — keep that).
- **Extraction tilts are in the DOF set from the start** (per mirror,
  the wall slice's lesson: a pupil knob that costs a little wavefront —
  do not repeat the omission that cost the 4-mirror family 35% of its
  blur).  DOF list explicit in the record: conics, spacings/standoffs,
  front end, tilts.
- Known traps, all pre-paid: the `resolve_nmirror_` psi parity rule
  (starts at k=3); verify every layout against its own paraxial
  prediction ON AXIS before any metric; degenerate-closure walls
  (spacing ≥ 20 mm); expect Mersenne/confocal flat valleys at high N;
  seed from the committed deck's own re-solved front end, never
  `P.parent`'s raw one; `iface` from `zElt`.
- Byte-for-byte: with N=4 and tilts zeroed the generalized build must
  re-emit the committed decks exactly — asserted in code.

## Task 1 — the top of the ladder (N=7, everything met, with margin)

One requirement set, one footing, all rungs (the S4 study's own
targets):

| requirement | target | held as |
|---|---|---|
| WFE rung-2 max over the field box | ≤ 71 nm | target |
| pupil blur | ≤ 47 µm | target |
| footprint wander | ≤ 56 µm | target |
| magnification breathing | ≤ 0.4 % | target |
| interface surface (rim-anchored — the S4c spec rule) | ≤ 0.2 mm | target |
| M / exit beam / collimation | 30× / ~33.6 mm / collimated | closure identities |
| union body-in-beam floor (declared body) | ≥ +15 mm | wall |
| last powered mirror behind M1 | ≥ 500 mm | wall |
| spacing, AOI | ≥ 20 mm, ≤ 15° | walls |
| rays lost over the field box | 0 | gate |

- **The top rung must have SLACK everywhere** — report the margin on
  every row.  A top rung sitting on a target is a cliff, not a start.
- Multi-seed (the S4b basin lesson): ≥ 3 independent seeds; report
  basins and spreads.
- **Power-economy tie-breaker**: solve with a small regularizer
  preferring least Σφᵢ², so power concentrates in the mirrors that
  matter and the near-flat ones identify themselves.  Doctrine note: it
  is a TIE-BREAKER, never a target — assert non-interference by
  re-solving the top without it and showing no requirement margin moves
  materially.  If it does interfere, drop it and say so; the removal
  ranking falls back to probes alone.
- Convergence discipline throughout the slice: central differences,
  converged not budget-capped, exitflags + evaluation counts per solve,
  walls on iterates never on reports, sampling bias carried per point.

## Task 2 — the descent (removal policy: predictions rank, measurements decide)

Per rung, from the current design:

1. **Rank** removal candidates cheaply: |φᵢ| and the power-economy
   solve's own concentration.  Rankings are advisory only (the
   field-walk law's 5–6× standoff optimism is the standing warning
   against trusting closure arithmetic as a predictor).
2. **Probe** the top 2–3 candidates with a short power-to-zero
   continuation (φᵢ → 0 in steps, warm-started re-solves, walls ON,
   budget-capped as probes and labeled so).
3. **Commit** the best candidate's removal with a CONVERGED walk to
   φᵢ = 0; the flat is then deleted if the layout closes without it, or
   retained as a fold (listed, not counted — ruling 4).  Record the
   runners-up's failure modes: which requirement broke first is ladder
   data.
4. **Control** (layered A/B): at each rung, an independent multi-seed
   direct solve at that N.  If the control materially beats the walk
   rung, the walk carried vestigial geometry — adopt the control as the
   rung and SAY SO; that is a finding about the method, and the record
   wants it either way.
5. **Bottom of the ladder**: descend until a rung cannot meet the set —
   converged, walled, multi-seeded.  Autopsy the failure: which
   requirement, at what margin, with what the last feasible rung looked
   like.  Expect the bottom near N=4–5 (the 4-mirror family demonstrably
   cannot meet the wavefront half), but expectation is not a result.

## Task 3 — the ladder (the deliverable)

One table, N descending, the delivery-table columns so everything reads
side by side: WFE rung-2 max, blur / wander / breathing, interface
surface (rim), union floor (declared + bare glass + sampling bias), max
chief AOI, retained flats (count, role), rays lost, exitflag /
evaluations, provenance (walk vs control), and per rung **what the
removed mirror bought** (the delta to the rung above).  Figures: the
ladder chart (each requirement vs N with its target line), a layout
render per rung (the standing review gate — a design is not quotable
until its layout has been looked at), and the failed-rung autopsy.

## Task 4 — the record and the gates

- Stage dir `challenges/afocal4/descent/` (numbers-first README);
  **`# DESCENT RESULTS` appended to `challenges/afocal4/RESULTS.md`**
  (the wall slice's precedent: appended sections, one findable entry
  point); earned rules continue the numbering; resolutions written at
  resolution time; delivery log at the foot of THIS brief.
- `tAfocal4Descent` registered in the suite.  Non-vacuity asserted in
  code: the top rung meets every row of the requirement table; the
  recorded failed rung fails it; the N=4 byte-identity of Task 0 holds.
  `tAfocal4` / `tAfocal4Clear` / `tAfocal4Wall` stay green.
- Everything additive under `challenges/afocal4/`; no engine work; no
  text-parsing `.in` files; committed decks rebuild byte for byte.

## Scope guards

- NO fold-route/envelope optimization (the follow-up packaging round
  owns it).  NO fifth-mirror architecture work beyond what the ladder
  itself measures at N=5.  The merit doctrine (log-domain targets,
  walls-not-terms) stands; the only new doctrine item is the
  power-economy tie-breaker, gated by its own non-interference check.
- Outward-facing text waits on Dave's sign-off via CC — including any
  "how many mirrors does the spec need" headline, which people will
  quote.

---

# DELIVERY LOG (2026-09-01)

Everything new lives in `MACOS_res_dev/mmacos/challenges/afocal4/descent/`
(numbers-first README).  Canonical record: **`# DESCENT RESULTS` appended to
`challenges/afocal4/RESULTS.md`** (§ D.1-D.7), the wall slice's precedent.
Local, unpushed.  No engine work.  Base: the wall slice's handback commits.

## The headline, and it is not the shape the brief expected

**No mirror count in this family reaches the requirement set.**  The stage was
to start at seven mirrors where everything is met with margin and walk down;
**the top rung was never reached**, and what the stage delivers is the reason,
measured three independent ways.

With the pupil requirement **abandoned entirely** and every DOF free -- the
most optimistic wavefront this family can produce -- seven mirrors floor at
**3424 nm against a 71 nm target (48x)**, and three extra mirrors buy **11 %**
over four.  So the ladder's question has an answer of a different KIND than a
rung count.

## 1 - Task 0: the N-mirror machinery, done and gated

The closure generalizes in three lines and they are `afocal4_close`'s own:
`phi_N` recollimates (analytic), `t_{N-1}` sets M (analytic), and only
`phi_{N-1}` needs a root because the chief ray is the residue.  Verified two
ways at N=4 -- against `afocal4_close` (max|dR| 4.4e-16, residuals 0 /
2.2e-16 / 0) and through the builder against the committed deck, **byte for
byte**, under `afocal4_phi4`'s own scan window.

Two things had to be derived, not inherited: the exit marginal height's SIGN
is a property of the layout at general N (both are closed; the one giving a
positive last spacing wins), and the pole problem is worse, so every candidate
is closed and CHECKED.

Machinery: `descent_close`, `descent_build` (tilts in the DOF set from the
start, applied upstream-first so they compose), `descent_seed`,
`descent_require` (TARGETS with margin / WALLS with room left / GATES as
facts; interface surface RIM-anchored per the S4c spec rule), `descent_solve`,
`descent_remove`, `descent_add`, and per-point runners.  Gate
**`tAfocal4Descent` 6/6**, registered in `SUITE_FREEFORM`.

## 2 - Tasks 1-3: what every route said

| route | N | merit | WFE nm | blur um | M err % | verdict |
|---|---|---|---|---|---|---|
| committed | 4 | 30.2 | 10407 | 157 | 0.0221 | the S4b/S4c delivery |
| cold seed | 7 | 70.78 | 12422 | 506 | 3.86 | missed, stalled |
| cold, radii freed | 7 | 70.40 | 11718 | 534 | 3.68 | missed, stalled |
| cold seed 2 | 7 | 707 | 3.7e9 | 1.7e6 | 95 | scrambled, gates caught it |
| cold seed 3 | 7 | 53.57 | - | - | - | missed, worst 177x |
| cold seed | 8 | 146.7 | - | - | - | missed, worst 1768x |
| ascent | 5 | 37.39 | 10775 | 332 | 0.1379 | missed |
| ascent | 6 | 44.20 | 9137 | 721 | 0.0651 | missed |
| ascent | 7 | 42.66 | 7894 | 705 | 0.0632 | missed |

The wavefront-only floors (pupil abandoned, every DOF free):
**N=4 3841.8 / N=5 8077.4 / N=6 5689.0 / N=7 3424.2 nm** -- 54x / 114x / 80x /
48x the target.  Upper bounds (rounds still gaining 18-25 %) and labelled so;
closing 48x would need two orders of magnitude where S4c's 17x-budget solves
moved the same designs 0.02-10.8 %.

Figure: `descent/afocal4_descent_ladder.png`.

## 3 - The findings that are reusable beyond this stage

* **THE PACKAGING STATION OBEYS A PARITY LAW.**  `z_N = sum (-1)^k t_k`, so
  the closure's own last spacing enters with sign `(-1)^(N-1)`.  Compliance
  rate: **88.4 % / 0.03 % / 89.7 % / 0.00014 %** at N = 5/6/7/8 -- a factor of
  ~3000 between adjacent N from one sign.  S4b's "one extra mirror flips the
  parity of the back end", with a rate attached.
* **AND IT DOES NOT TRANSFER TO A SINGLE REMOVAL** (predicted too strongly by
  me, then measured).  Retain clears 3 of 3, delete 1 of 3 -- but M5's delete
  clears, because deleting element k MERGES `t_{k-1}+t_k` and re-signs
  everything after it.  **"N = 6 cannot be built" is false**: at least four
  routes exist from the N=7 base.
* **WHAT THE PUPIL REQUIREMENT COSTS, priced.**  At 343 mm, dropping it takes
  the wavefront 10407 -> 3842 nm, a factor of **2.7**.  S4 ran the same A/B at
  140 mm, got 4 %, and generalized it as "the DOFs do not touch it" -- that
  generalization is operating-point specific and should not be carried.
* **A COLD CLOSURE IS A SPECIFICATION, NOT A DESIGN.**  Four cold seeds landed
  worse than the four-mirror family with twice the freedom; the same machinery
  warm-started produced sound rungs.  A closure holding its conditions at
  1e-16 says nothing about quality -- one such probe traced M = 40.45 against
  a paraxial 30.0000.

## 4 - Two of my own hypotheses, refuted, and one recorded wrong prediction

Recorded because the corrections are the useful part:

* **"the DOF set is the reason"** -- the wall slice's own finding, reached for
  first.  Freeing the radii bought **5.7 %** on a row needing 165x.
* **"the tilts are the missing freedom"** -- the wavefront-only control
  WITHOUT tilts floors at **3841.8 nm, BETTER** than 4497.7 with them.  Tilts
  are a pupil knob, exactly as the clearing stage measured.
* **the parity prediction**, written as a rule and true only as a rate (above).

*A lesson confirmed once does not become the explanation for the next stall.*
Both times, the SEED was the reason.  Earned rules 35-38 in § D.6.

## 5 - Method defects caught, each before it cost a result

* **A grid is a grid**: the N=4 compliance row first read "no compliant
  closure" for the design sitting in the repository complying at +1323 mm.
  The parent's own spacings are now always injected.
* **A wall with only one side is not a constraint**: the first N=7 seed put
  the last powered mirror **10.96 m** behind the primary and was compliant by
  every check.  Bounds now two-sided, in multiples of the M1-M2 spacing.
* **The insertion needs a compliant seed too**: a naive midpoint split put the
  child 479 mm IN FRONT of the primary -- the parity law charging for the
  extra reflection.  `descent_add(...,'search',true)` scans for it.
* **`afocal4_score`'s failure path returns a MINIMAL struct** and
  `descent_solve` reached past it, turning a scored-as-bad iterate into a dead
  run.  Guarded.  `clear_solve` carries the same unguarded access and has
  never been handed that struct -- noted, not chased.

## 6 - Scope not done, and named

No fold-route or envelope work (the packaging round owns it).  No fifth-mirror
architecture beyond what the ladder measured.  **The descent proper -- walking
DOWN from a top rung that meets the set -- was not run, because no such rung
exists**; the removal machinery is built, gated and demonstrated, and the
first removal step is measured in § D.5.

## 7 - FOR DAVE

**The requirement set is not reachable by mirror count in this family, so the
next move is a judgement about what the work is for, not another solve:**

* **the spec** -- 71 nm was set at >= 10x Rodgers' best three-mirror (S3 gate
  review).  Nothing in this arc has come within 48x of it, with the pupil
  requirement abandoned and up to seven mirrors; or
* **the family** -- a coaxial all-reflective afocal whose interface-pupil
  condition consumes the last two powers.  § D.4 now prices that condition at
  a factor of 2.7 of wavefront.

Also still open from the wall slice: the **§S4b.4 correction wording**, held
for sign-off and uncommitted.
