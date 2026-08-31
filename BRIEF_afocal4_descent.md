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
