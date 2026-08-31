# BRIEF for TO: afocal4 stage — get the collimator out of its own feed beam

_Tasking: Terminal Opus.  Supervisor: CC.  Timing: START NOW (Dave
2026-08-30: "we have time to get this done") — but NOT demo-blocking:
deck slide 13 already declares the finding and "redesign queued", so
the demo stands whatever happens here; do not let this collide with
demo-day support.  Origin: the r2-packaging delivery's bigger finding
(delivery log at the foot of `BRIEF_r2_packaging.md` §2).
Cold-start reads: that delivery log + `challenges/afocal4/packaging/`
(README + `pack_route`/`pack_view`/`check_record` machinery); memory
`project_afocal4_rodgers2`; `challenges/afocal4/`
{FORM_STUDY,RESULTS,STATUS_S4B,STATUS_S4C}; the e2e extraction-tilt
precedent (memory `project_e2e_example`: M3 extraction tilt 1.2° took
the return leg off the feed axis); `project_fold_extraction` (fold
near a focus = small flat, small offset — the relay-is-a-small-patch
doctrine)._

## The defect, restated precisely (from the measurement)

In the committed 343 mm family-2 deck — UNFOLDED, before any
packaging — the collimator sits inside the converging M2→field-mirror
feed cone.  Per single field the feed passes it in a 27.8–55.6 mm
annulus with 10.8 mm of daylight; but a monolithic collimator must
cover its UNION footprint over the 0.5°×0.5° box (17.0→87.0 mm), and
that glass is where the other fields' feed beams pass: −55.4 mm
against bare lit glass, −79.9 mm with the 1.15×hull + 15 mm body
allowance, and −6.0 mm at the box centre alone.  Real hardware
vignettes a fraction of the field box; the trace does not show it
because bodies are not obscuring elements.  A fold is an isometry and
carries the conflict unchanged — this needs a DESIGN move, not a
packaging one.

## The leverage map (work these in order; each is priced, none is free)

1. **Extraction fold near the internal focus (first — cheapest
   physics).**  Family 2 was chosen partly BECAUSE it moves the
   internal focus well behind M1; near that focus the beam is
   smallest, so a SMALL flat can take the post-focus leg off the feed
   axis entirely — after which the collimator, field mirror and cold
   stop live beside the feed cone, not inside it.  A flat is
   aberration-neutral (null asserted, the packaging pattern), and
   unlike folding the existing train, inserting the fold BETWEEN the
   conflict partners re-routes one of them — the isometry argument
   does not apply.  Cost: one more reflection, fold-flat clearance
   (use the measured pair-step bound from `pack_route`), and the
   packaging Path-A folds must be re-run on the new leg geometry.
2. **Axial placement + sizing of the collimator within the cone.**
   The union demand (87 mm) and the local annulus both change with
   station; there may be a station where union glass fits inside the
   per-field daylight.  Cheap to SCAN before committing to anything:
   sweep the collimator station, plot union-footprint demand vs
   available annulus — if the curves never cross, say so with the
   figure and move on.
3. **Extraction TILT on the powered leg (the e2e precedent).**  A
   small tilt of the field-mirror/collimator leg separates it from
   the feed geometrically; it puts power on tilted surfaces, so it is
   priced in WFE and pupil metrics by a joint re-solve — report the
   price, do not eyeball it.
4. **A FIFTH mirror (Dave's addition, 2026-08-30).**  Add a powered
   element rather than contorting four: the natural form is a small
   powered mirror near the internal focus that does the extraction
   AND shares the collimation/pupil work — taking the body-in-beam
   conflict away geometrically while ADDING the freedom the pupil
   metrics currently pay for (the 4th mirror bought the pupil at a
   factor-39 image cost; a 5th may re-balance that exchange).  Costs
   to state: one more reflection and alignment, a re-opened element
   count (the "what a 4th mirror buys" story becomes "what a 5th
   buys" — a NEW trade row with the full pupil set), and the
   packaging re-run.  Try it when 1–3 clear only at an ugly price, or
   pursue it in PARALLEL as the counter-design if the station scan
   (2) shows the 4-mirror topology is fundamentally pinched.
5. **The family knob (M2 speed / internal-focus position)** — the
   most expensive lever: it re-opens the family-2 trade at its base.
   Touch it only if 1–4 cannot clear, and then re-report the FULL
   trade row (WFE, pupil blur, wander, magnification variation,
   pupil distance ↔ instrument diameter).

## Where the problem is constrained (hold these fixed)

- The INTERFACE: 30× afocal output (28.7× at the frozen-offset
  variant), ~33 mm collimated beam delivered to the tilted cold stop
  at the stated pose — this is the customer boundary.
- The FIELD: 0.5°×0.5° box, offset 0.6° — spec, not a knob.
- Last powered mirror BEHIND the primary (the constraint the whole
  343 mm trade paid for), and the envelope discipline the packaging
  stage established.
- The committed 343 mm trade row: never overwritten.  Any design
  move that shifts performance produces a NEW row beside it.

## The gate this stage must leave behind

The blind spot that let this ship: `afocal4_pack` checked beam-to-beam
daylight on the last leg at the box centre, and never asked whether a
BODY stands in a BEAM, union-sized, over the field.  Promote the
packaging stage's union body-in-beam measurement into the afocal4
gate set as a standing check: every optic's union footprint (with a
stated allowance — declare it up front; the packaging stage used
1.15×hull + 15 mm) against every beam segment, over the field box,
signed floors printed.  Non-vacuity: the gate must FAIL on the
committed 343 mm deck (it measures −55 mm there) and pass on the
cleared design.

## Deliverables, in priority order

1. The station-scan figure (leverage 2) — even if it only proves the
   scan cannot clear, that is the figure that justifies the fold/tilt.
2. The cleared design: `challenges/afocal4/clearing/` (new dir, new
   names — nothing committed is touched), with the union-clearance
   gate PASSING, zero rays lost, and the interface/field/packaging
   constraints held.
3. Re-run packaging Path A on the cleared design (the four-fold depth
   answer must survive, or its new form reported).
4. The trade quantities for the cleared design as a NEW row beside
   the committed 343 mm row — plus the pupil metrics (blur, wander,
   magnification variation): the 4th mirror exists FOR the pupil, so
   the cleared design must show what it kept.
5. The standing union body-in-beam gate, with its non-vacuity check.
6. README (numbers-first), delivery log at the foot of THIS brief,
   memory update extending `project_afocal4_rodgers2`.
7. Slide refresh afterward = CC, under the §5 gate.

## Traps (paid for once already)

- matlab -batch: script files + exit(0) in runners; one model size
  per process; MACOS_HOME=/home/dcr/dev/macos/macos_f90.
- Fold-insertion null asserted per flat; +90°/−90° pair step > 2× the
  beam half-extent on the first flat (the measured `pack_route`
  bound).
- A gate's margin is a number, not a body — floors from the signed
  model, allowances declared.
- Engine truth for every geometric claim (`ray_hist`), never .in
  text-parsing.
- Frame-before-angle on every reported tilt/fold angle.
- Read-tool PNGs cache-stale — unique filenames, verify by printed
  numbers.
- Work in MACOS_res_dev; ~/dev/MACOS_resources is the demo checkout —
  leave it alone.  Commits LOCAL; Dave orders pushes.  No engine work.

---

# DELIVERY LOG (2026-08-30)

Everything new lives in `MACOS_res_dev/mmacos/challenges/afocal4/clearing/`
(the README there is numbers-first and is the full account).  **Nothing
under `challenges/afocal4/` was overwritten.**  Two additions beside it,
both purely additive: `../afocal4_union.m` (the gate the brief asks to be
promoted) and one clause in `../afocal4_pack.m` that calls it.  Local,
unpushed.  No engine work.

## The answer, in four lines

1. **The interference is structural and obeys a law.**  The collimator's
   footprint and the feed beam's are two scaled copies of the SAME
   off-axis field box, so they can separate only if their scales differ
   by more than `(bias+half)/(bias-half) = 0.85/0.35 = 2.4286`.  Measured
   ratio **1.3025**.  And one of the two scales is **pinned by the
   customer interface**: `c_body(collimator) = M x iface` — 10.290 m/rad,
   measured 11.118 (+8.1 %, the design's own pupil aberration).
2. **Leverage 1 (a flat extraction fold) cannot do it**, exactly.  The
   collimator lands ON the flat at `DIST = LEG - B` (0.808 of the leg
   here) and the two failure modes are complementary about that station.
   Best over every station and both turn directions: **-79.89 mm — the
   parent's own number, to the last digit.**
3. **Leverage 2 (the station) cannot either.**  Over **1.40 m** of
   collimator travel the ratio reaches **1.48**; demand 83.6…100.3 mm
   against 23.0…34.1 mm available.  And **every point of the committed
   trade curve fails the gate** — 50/90/140/220/343 mm read **-29.0 /
   -33.4 / -54.2 / -67.5 / -79.9 mm** with the declared body.  At 50 mm
   the COLLIMATOR pair itself does clear (+36.8 mm, ratio 11.7) — but its
   collimator sits **5.21 m** behind the primary and the deck then fails
   on a different pair.  (Report both columns; the one-column version of
   this table reads "50 mm clears".)
4. **Leverage 3 (an extraction tilt on the field mirror) does it**,
   because a tilt is the one remedy the law does not bind: it separates
   the bundles by a FIELD-INDEPENDENT `2*alpha*d`, measured as a fitted
   offset of **+201 mm**.  Gate **-79.89 -> +37.82 mm**, zero rays lost,
   and the price falls entirely on the pupil.

## The delivered design

`clearing/afocal4_clear_343mm.in` — the committed deck with a -10 deg
extraction tilt on the field mirror and conics + field-mirror standoff +
front end re-solved around it (427 evaluations).

| | committed 343 | cleared -10 deg | ratio |
|---|---|---|---|
| WFE rung 2 max (nm) | 10406.98 | **8992.68** | **0.864** |
| pupil blur rms (um) | 157.02 | 553.34 | 3.524 |
| breathing chief-normal (%) | 0.1240 | 0.8160 | 6.579 |
| wander at refit plane (um) | 161.23 | 559.87 | 3.472 |
| M at box centre | 30.0066 | 30.0148 | 1.000 |
| exit beam (mm) / collimation (urad) | 33.633 / 1477 | 33.571 / 1266 | 0.998 / 0.857 |
| max chief AOI, any mirror (deg) | 12.84 | **10.86** | 0.846 |
| **union body-in-beam floor (mm)** | **-79.89** | **+37.82** | — |
| rays lost over the field box | 0 | 0 | — |

**The wavefront and the interface do not pay** — WFE 13.6 % better, the
exit beam and M held, collimation 14 % better, and the design's *worst*
worked mirror improves (12.84 -> 10.86 deg, inside the 15 deg standing
rule) even though the field mirror itself goes 2.1 -> 10.1 deg.  **What
it costs is the fourth mirror's pupil control**: blur and wander 3.5x,
breathing 6.6x.

## And it packages itself — Path A is not needed and does not close

| | committed | committed + Path A (4 flats) | **cleared, no flats** |
|---|---|---|---|
| deepest optic behind M1 (m) | 1.8866 | 0.8932 | **1.2874** |
| … x the M1-M2 spacing | **1.81x** | **0.86x** | **1.24x** |
| overhang (deepest − spacing, m) | +0.845 | −0.149 | **+0.246** |
| … x the M1-M2 spacing | 0.81x | −0.14x | **0.24x** |
| radius of the optics BEHIND M1 (m) | 0.186 | 0.435 | **0.150** |
| back focal path (m) | 2.808 | 2.808 | **2.192** |
| extra flats | 0 | **4** | **0** |
| union floor (mm) | -79.9 | -79.9 | **+37.8** |

*(Two ratios, named apart deliberately: the packaging stage's headline
"1.81x" is the DEEPEST OPTIC over the M1-M2 spacing, while the overhang
over the same spacing is 0.81x on the same deck.  Quoting one under the
other's name makes the improvement look 3x bigger than it is.)*

The overhang the packaging stage spent four ~300 mm 45-deg flats to
remove, the swing takes from 0.81x to **0.24x** with none.  It does not
reach Path A's negative overhang — but the back-end **girth shrinks**
(0.186 → 0.150 m) where the four folds nearly **tripled** it to 0.435 m,
the back focal path shortens by **0.62 m** (a fold is an isometry and
leaves it), and it opens **no** polarization budget.  Path A was
re-searched over all four of its stated quantities on the cleared deck's
own geometry — 96 routes tried, 15 satisfied both the route algebra and
the plane-intersection bound, and **every one of those lost rays** (best:
-72.74 mm, 592 rays) — **the trombone runs out of leg, because the depth
it was invented to remove is mostly gone.**

And the rest of `afocal4_pack` PASSES on the cleared deck, so this does
not trade one buildability clause for another: last powered mirror
**+774 mm** behind M1 (min 500), fold daylight **27.1 mm** (margin 15),
instrument 233 mm off axis with the **largest that fits Ø421 mm**
against the stated Ø300, union **+37.8 mm**.

## The standing gate (deliverable 5)

`afocal4_union.m`, beside `afocal4_pack`: union footprint over the field
box, convex hull (never a centred disk), **declared** allowance
(1.15 x footprint + 15 mm, printed with every answer, and every table
also carries bare lit glass), sampling-free pierce and distance.
`afocal4_pack` runs it as part 4 and folds the verdict into `K.ok`;
`'union',false` reproduces the old verdict exactly and the three
sub-flags `tAfocal4` asserts are untouched.

**Non-vacuity, asserted in code**: at the same allowance it FAILS on the
committed 343 mm deck (-79.89 mm) and PASSES on the cleared one
(+37.82 mm).  `tAfocal4` stays **8/8**; new `tAfocal4Clear` **8/8**
(`SUITE_FREEFORM`, registered in `run_mmacos_tests.sh`).

## Nulls, all at machine precision

* the design struct recovered from the committed deck rebuilds it
  **byte-for-byte**;
* `clear_build(tilt 0)` == `afocal4_build`, **byte-for-byte**;
* the tilt moves the chief-ray path up to the swung mirror by
  **4.45e-16 m**, the chief still lands on the pivot to **4.45e-16 m**,
  and the beam turns **20.0000 deg** for a 10.0000 deg tilt;
* the measured offset runs **19.35 mm/deg** against the `2*alpha*d` the
  geometry predicts, **19.66** — 1.6 %.

## Things that cost a cycle, recorded

* **`iface` must be recovered from `zElt`, not from vertex geometry.**
  The builder poses the interface plane on the traced chief, so on this
  deck the last mirror's vertex is 359 mm from the interface vertex
  while the standoff is 343.  The vertex reading rebuilds a deck that is
  not the committed one, and it silently shifted a whole scan.
* **An all-failed sweep must be an error, not an empty result.**
  `pack_fold` lives in `../packaging`; a caller without it on the path
  got a clean "a flat cannot clear it" out of nine undefined-function
  exceptions — the answer the study was looking for.  `clear_fold` now
  refuses an empty sweep and prints every station's failure.
* **Do not gate a tilt on the unsigned AOI moving by alpha.**  The field
  mirror is worked at 2 deg, so any larger tilt carries the SIGNED angle
  through zero and |AOI| folds back.  Gate the beam TURN (2 alpha) and
  the pivot instead.
* **Fit the field walk WITH an intercept.**  Forcing it through the
  origin makes a tilted design report a meaningless "walk" that is
  silently absorbing the offset — and the offset IS the remedy.

## Open, for Dave

1. **The clearance is not yet a wall and the re-solve spends it.**  At
   -8/-9 deg the solver walks +23.3/+42.3 mm of margin down to
   +2.3/+0.7 mm, because the merit cannot see it.  The fix is the S4b
   pattern: a union-floor wall in `afocal4_build` beside
   `m3_behind_min` — with a compliant seeder, or it is a cage.  Not done
   here; it is its own slice.
2. **35 % of the pupil blur was lying around unclaimed**: a *positive*
   4 deg tilt takes blur 157.0 -> 102.6 um with no re-solve at all.  The
   committed design is not at its own pupil optimum, and the reason is
   visible in the merit (a wavefront term 130x off target owns the sum
   of squares).
3. **The re-solves are budget-limited, not converged** (exitflag 0 at
   427 evaluations) and ran on the study's default 3e-3 FORWARD
   difference — the setting S4c measured as reading the gradient 17 %
   low.  A central-difference polish is the known next step.
4. **Leverage 4 priced, not built.**  A fifth mirror must supply at
   least **111.6 mm** of field-independent separation, which the -10 deg
   tilt already supplies (201.4 mm).  Its case therefore rests on doing
   it WITHOUT spending the pupil control — i.e. on beating blur 553 um /
   breathing 0.816 % / wander 560 um at 8993 nm.  Three architectures
   and which lever each pulls: clearing README section 11.
5. **Slide 13** can now say what the redesign is and what it costs.
   Outward-facing, so it waits on sign-off (`doc/STYLE_REPORTS.md` §5).
