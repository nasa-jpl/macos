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
