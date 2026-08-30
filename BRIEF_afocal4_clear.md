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
