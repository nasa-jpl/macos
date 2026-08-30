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

0. **Both demo MATLABs run from `~/dev/MACOS_resources`** — the
   checkout named like the public repo, held DETACHED at
   origin/dev-candidate (it shares one git with MACOS_res_dev, which
   keeps the dev-candidate branch; 2026-08-30).  Refresh procedure
   after any push:
   `git -C ~/dev/MACOS_resources fetch origin`
   `git -C ~/dev/MACOS_resources checkout --detach origin/dev-candidate`
   `make -C ~/dev/MACOS_resources/mmacos`   (mex; no-op when current)
   Verified 2026-08-30: mex built (gfortran engine), tTgPol 9/9 +
   tTgPol2 9/9 + tBench 7/7 from this tree.

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

- **Slide 8 (design family / frontier):** the frontier table on
  screen — the context for the ask.
- **Slide 9 (demo intro):** explains the beat (ask, prediction-first,
  no-AI-in-this-loop, what the reveal shows) with the rehearsal
  verdict verbatim on the right.  **DAVE switches to MATLAB here and
  kicks off `oi_demo_step(<width>)` himself** — the without-AI half
  of slide 5's objective, running live.
- **Slides 10–21 are the cover** (C1, C2 ×2, the IFO block incl. the
  TG live beat, the CTB trio, the e2e6m observatory trio) —
  ~26–32 min against the 15.2-min solve.  COMFORTABLE.
- **Slide 22: the reveal** — map + layout + gates beside the stated
  prediction, right before the closing discussion.  The reveal
  windows AUTO-POP in the solve MATLAB at completion (2026-08-30):
  live-traced layout + WFE map + fields; `oi_demo_show(OUT)`
  re-renders, `oi_demo_show()` no-arg = newest run on disk.
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
4. Inside the toolbox — the mmacos directory map (veneer,
   sensitivities, design layer, runners, templates/tests; NEW
   2026-08-30, fig_mmacos / make_mmacos_fig.py).

**AI-driven design (Rodgers arc — themes per the rodgers3 deck)**
5. Three challenges, one method — tasks quoted as posed + ladders.
6. Challenge 3: the problem + all five steps reproduced (the license).
7. Challenge 3: the numbers — 113.6 @ 34.1 mm (0.97×); 45.4 vs 53
   (0.86× at a floor 1.6 mm shy), the solve-field diagnosis.
8. Design family: walk + frontier (the ask's context).
9. **Demo intro — Dave switches to MATLAB, KICKS OFF the solve.**
10. Challenge 1: design-vs-design table (0.97×) + findings.
11. Challenge 2 reproduced: variant table 0.95–1.02×.
12. Challenge 2 extended: pupil metrics (+ convergence-cloud strip)
    + the 4th-mirror trade.

**The IFO demo (rehearsed with TO 2026-08-29 — good to go)**
13. Measuring a surface with polarized light — the PRINCIPLE (why
    TG, analyzer angle = phase step, factor 2; sweep-frames figure).
14. The rig, built + calibrated — bilinear sweep + the 11.7% catch.
15. **LIVE: Dave measures the DM** (7 steps with timings, solo,
    desktop MATLAB; poke + recovery figures on the slide).

**The coronagraph testbed (highlights from deck_ctb v5)**
16. CTB end to end: layout sketch + mask families head-to-head +
    PROPER external validation (correlation 1.000000).
17. The vortex against the Lyot stop — THE CHARGE-4 CONCLUSION
    (8.8e-11 at the band-limited mask's own 36%).
18. The mirrors close the loop: 1.7e-8 → 6.8e-15 (numerical floor);
    physics vs arithmetic floors; pol/bandwidth/vector-vortex
    verdicts.

**The observatory (e2e6m — testbed proves the models, the same
modeling builds the flight-like system)**
19. From the testbed to an observatory: light-order sketch + shroud
    packing (7.451 m vs the 8.0 m gate).
20. Segmentation, metrology, controls at 6 m scale (114 MET beams).
21. A randomly disturbed observatory, simulated: 2.5e-7 held vs
    1.7e-6 open — ready to ingest dynamical/thermal time histories.
    WORK IN PROGRESS, declared on the slide.
22. **LIVE REVEAL** — prediction vs the just-solved design.

**Close**
23. AI in design + analysis — the working questions.
24. Summary + public code.

**Backup** (after a plain divider, slide 25)
26. What MACOS does today (moved from main, Dave 2026-08-29 pm).
27. The MATLAB toolbox + design layer (moved from main).
28. Validation anchors (honest PROPER spans).
29. Adjacent-design rehearsal bundle (7/12/14° table).
30. IFO demo backups (one PNG per live step).
31–33. The e2e worked example trio (design / segmentation+MET /
    compare+simulate) — moved from main when e2e6m became the
    featured example.

**TG v2 (delivered 2026-08-29, commit 5cd2112 LOCAL):** slide 13 now
carries the v1-vs-v2 table + the error-budget-inversion figure;
either rig runs the same 7 demo steps (v1 49.8 s / v2 47.0 s) — RIG
CHOICE = Dave's dry-run call; v1 stays the rehearsed default.

**Timing note:** kickoff at slide 9 ≈ 15–17 min in; cover = slides
10–21 including the TG beat ≈ 26–32 min vs the 15.2-min solve —
COMFORTABLE (reveal at 22, before the discussion).  The
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
