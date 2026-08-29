# Keysight/CodeV demo — session outline (Dave's working copy)

_Edit freely — this is YOUR session-day copy.  The durable record
stays `BRIEF_demo_deck.md` (linked here); if you change the structure,
tell CC and the brief gets synced.  Demo ~2026-09-01._

## This folder

| entry | what it is |
|---|---|
| `decks/` | the five harvest decks (rodgers3 V3 + walk, CTB, e2e6m ×2) |
| `tg_psi_dm/` | live beat (a): `demo_tg_psi.m`, backup PNG per beat, README |
| `offset_imager/` | live beat (b) machinery: `oi_demo_step.m`, `run_oi_demo.m` |
| `demo_adjacent/` | beat (b) fallback bundle (7/12/14° PNGs + verdicts) + **REHEARSAL.md** (the spoken script) |
| `rodgers3/` | the challenge dir (PACKET, records, deck sources) |
| `BRIEF_*.md` | deck plan + the two live-beat briefs with delivery logs |

## Session-day checklist (before the room fills)

1. Desktop MATLAB session #1 — for the TG beat.  Fresh process
   (model 256; one model size per process).  `cd` into `tg_psi_dm/`.
2. CC terminal session — visible on the projector for beat (b).
   MATLAB license confirmed working (both sessions opened BEFORE).
3. `demo_adjacent/` PNGs + verdict texts known-good (the fallback
   ladder: this bundle → committed walk artifacts in
   `offset_imager/t5_walk/` → frontier table slide).
4. Backup PNGs for the TG beats confirmed present in `tg_psi_dm/`.

## Live-beat structure (REVISED, Dave 2026-08-29 — Challenge-3-first)

**The Rodgers arc now opens with Challenge 3, and the ask+launch
happen at the frontier slide (7)** — the audience chooses a width
WITH the frontier on screen, CC visibly predicts + launches, and the
reveal lands at slide 13, right after the TG beat.  A solve at the
quality bar is ~15.2 min (measured):

- **Slide 7 (design family / frontier):** the frontier table on
  screen — the context for the ask.
- **Slide 8 (demo intro, NEW 2026-08-29):** explains the beat (ask,
  prediction-first, no-AI-in-this-loop, what the reveal shows) with
  the rehearsal verdict verbatim on the right.  **DAVE switches to
  MATLAB here and kicks off `oi_demo_step(<width>)` himself** —
  the without-AI half of slide 4's objective, running live.
- **Slides 9–20 are the cover** (C1, C2 ×2, the IFO block incl. the
  TG live beat, the CTB trio, the e2e6m observatory trio) —
  ~26–32 min against the 15.2-min solve.  COMFORTABLE.
- **Slide 21: the reveal** — map + layout + gates beside the stated
  prediction, right before the closing discussion.
- Refusal path in RESERVE (scripted in `demo_adjacent/REHEARSAL.md`,
  not demonstrated).  NARRATION NOTE: a wide ask (13–15°) now
  PASSES both gates — the committed walk line understates the
  envelope (TO's measured correction; REHEARSAL.md has the wording).
  Do NOT use the old "honest deficit" line for interior widths.

**Beat (a) = Dave SOLO** in desktop MATLAB: `demo_tg_psi` — 7 beats,
pauses between beats, each prints its numbers and has a backup PNG.
OPEN: keep/cut beat 3 (the 11.7% misalignment catch — TO+CC recommend
KEEP; your dry-run call).

## The slides (RESTRUCTURED per Dave 2026-08-28 — DRAFT deck built:
## `deck_keysight.md` → `deck_keysight.pptx`, both in this folder)

**Intro**
1. Title / agenda.  *(the beat-(b) spec ask + VISIBLE launch happen here)*
2. History — one slide (Dave's imagery + sign-off PENDING).
3. One engine, four surfaces.

**AI-driven design (Rodgers arc — themes per the rodgers3 deck)**
4. Three challenges, one method — tasks quoted as posed + ladders.
5. Challenge 3: the problem + all five steps reproduced (the license).
6. Challenge 3: the numbers — 113.6 @ 34.1 mm (0.97×); 45.4 vs 53
   (0.86× at a floor 1.6 mm shy), the solve-field diagnosis.
7. Design family: walk + frontier (the ask's context).
8. **Demo intro — Dave switches to MATLAB, KICKS OFF the solve.**
9. Challenge 1: design-vs-design table (0.97×) + findings.
10. Challenge 2 reproduced: variant table 0.95–1.02×.
11. Challenge 2 extended: pupil metrics (+ convergence-cloud strip)
    + the 4th-mirror trade.

**The IFO demo (expanded 2026-08-29 pm)**
12. Measuring a surface with polarized light — the PRINCIPLE (why
    TG, analyzer angle = phase step, factor 2; sweep-frames figure).
13. The rig, built + calibrated — bilinear sweep + the 11.7% catch.
14. **LIVE: Dave measures the DM** (7 steps with timings, solo,
    desktop MATLAB; poke + recovery figures on the slide).

**The coronagraph testbed (REDONE 2026-08-29 pm, before the e2e
story — highlights from deck_ctb v5, measured results)**
15. CTB end to end: mask families head-to-head table + PROPER
    external validation (correlation 1.000000).
16. The vortex against the Lyot stop — THE CHARGE-4 CONCLUSION
    (8.8e-11 at the band-limited mask's own 36%).
17. The mirrors close the loop: 1.7e-8 → 6.8e-15 (numerical floor);
    physics vs arithmetic floors; pol/bandwidth/vector-vortex
    verdicts.

**The observatory (e2e6m — the coronagraph-bearing design; REDONE
2026-08-29 pm per Dave: testbed proves the models, then the same
modeling builds the flight-like system)**
18. From the testbed to an observatory: light-order sketch + shroud
    packing (7.451 m vs the 8.0 m gate).
19. Segmentation, metrology, controls at 6 m scale (114 MET beams).
20. A randomly disturbed observatory, simulated: 2.5e-7 held vs
    1.7e-6 open — ready to ingest dynamical/thermal time histories.
    WORK IN PROGRESS, declared on the slide.
21. **LIVE REVEAL** — prediction vs the just-solved design.

**Close**
22. AI in design + analysis — the working questions.
23. Summary + public code.

**Backup** (after a plain divider, slide 24)
25. What MACOS does today (moved from main, Dave 2026-08-29 pm).
26. The MATLAB toolbox + design layer (moved from main).
27. Validation anchors (honest PROPER spans).
28. Adjacent-design rehearsal bundle (7/12/14° table).
29. IFO demo backups (one PNG per live step).
30–32. The e2e worked example trio (design / segmentation+MET /
    compare+simulate) — moved from main when e2e6m became the
    featured example.

**TG v2 (delivered 2026-08-29, commit 5cd2112 LOCAL):** slide 13 now
carries the v1-vs-v2 table + the error-budget-inversion figure;
either rig runs the same 7 demo steps (v1 49.8 s / v2 47.0 s) — RIG
CHOICE = Dave's dry-run call; v1 stays the rehearsed default.

**Timing note:** kickoff at slide 8 ≈ 14–16 min in; cover = slides
9–19 including the TG beat ≈ 25–30 min vs the 15.2-min solve —
COMFORTABLE (reveal moved before the discussion).  The
rehearsal-bundle backup slide covers a hard failure.  MATLAB
session for the kickoff must be a SEPARATE process from the TG
session — the solve BLOCKS its session for ~15 min, and the TG
beat runs inside that window.  Have both pre-opened.

## Open items

- Dave's TG dry run + beat-3 keep/cut ruling.
- `tOiDemoStep` suite registration (out on runtime — Dave's call).
- Deck assembly (STYLE_REPORTS §5) — after the dry run; history
  slide + slide 24 stay a Dave+CC pass.
- 22b work is LOCAL/uncommitted in the res tree — commit ordering
  with TO, push on Dave's word.
