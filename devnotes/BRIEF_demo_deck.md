# BRIEF: Keysight/CodeV demo deck + live beats — the working plan

_Dave + CC, 2026-08-28.  Demo ~2026-09-01.  This file is the durable
outline (Dave approved the structure 2026-08-28); the deck build runs
under doc/STYLE_REPORTS.md §5 (pre-write gate + Dave's source sign-off
on outward-facing claims).  Cold-start: memory `project_keysight_demo`,
`BRIEF_tg_demo.md` (live beat 1, tasked to TO), this file._

## The 25 slides (approved structure)

> **RESTRUCTURED (Dave 2026-08-28) — the live order is now
> `demo_session/deck_keysight.md` (DRAFT .pptx built beside it):**
> intro (1–3) → Rodgers arc + live reveal (4–10) → TG live gauge
> (11–12) → capability compressed to 2 (13–14) → illustrations
> (15–19) → close (20–21) → backup (anchors, rehearsal bundle, TG
> step backups).  Old slide 9 (Trust) moved to backup; old 4–8
> compressed; section V moved ahead of capability; r1/r2 slides
> harvest `MACOS_sandbox/slides/deck_rodgers_status.md`.  The list
> below is retained as the CONTENT INVENTORY.

I. Intro
 1. Title / agenda.
 2. The family: ONE engine, four surfaces (macos CLI / smacos lib /
    mmacos / pymacos) + the shared macos_api_mod layer, one Rx language.
 3. HISTORY — ONE slide (Dave's ruling): timeline late-1980s NASA tech
    dev -> 1991 Hubble recovery -> JWST modeling + WFSC development ->
    TMT -> validated vs testbed + flight data, many others.
    *(imagery + claims: DAVE'S sources, his sign-off — do not fabricate)*

II. Capability
 4. Current applications: system-level modeling — error budgeting,
    integrated simulation, WFSC for segmented telescopes, coronagraph
    testbeds + instruments.
 5. mmacos capabilities (two columns).  DRIVERS: design layer
    (System/Telescope/Bench), dw_dx/dw_dz/dw_dsurf/dw_dgrid
    supervisors, run_sensitivities / run_met / run_compare /
    run_simulator / run_segmentation, study orchestrators (ctb_study).
    UTILITIES: trace/opd/perturb/fex/pupil_find, view_rx/view_std/
    draw_rays, jones_pupil/polarizer/waveplate, segment_grid_basis +
    grid IO, compose, spot/intensity/complex_field.
 6. The linear-model workflow: w = J·x + w0, channels, multi-field x
    multi-config stacking, groups.
 7. Streamlining: design -> segmentation -> sensitivities -> MET ->
    compare -> simulate; one-call reconstructible studies,
    checkpoint/resume.
 8. The design layer + why: fully public-domain — no proprietary Rx
    needed.

III. Trust
 9. TEST SUITES + ANCHORS (Dave's add — the handshake slide for a
    CodeV audience): geometric <-> CodeV comparison suite (6,601
    tests); physical optics <-> PROPER (1e-11..1e-13); polarization
    <-> textbook closed forms (Born & Wolf, Abeles/Macleod) + published
    lab data (protected-Al Mueller anchor ~1e-14 vs the paper);
    cross-compiler bit-identity (ifx/gfortran, GMI 6/6); the
    NON-VACUITY rule — every gate shown to fail the pre-fix engine.

IV. Illustrations
 10. e2e: telescope + Offner relay, closed-form seed -> DL ±2'.
 11. e2e: segmentation, physical apertures, per-segment channels.
 12. e2e: MET + closed-loop simulator — design-to-control with numbers.
 13. Bench builder -> interferometer surface gauges: the Twyman-Green
     rig (`macos.design.twyman_green`, shipped example =
     `templates/90_polarization/tg_psi_dm/`, NOT 40_benches — that
     name was taken) + the pol-PSI variant = live beat 1.  Include the
     TG-vs-MZ trade note (Dave ruled TG for DM figure; MZ's case =
     isolation/two ports/transmission/dynamics).  MECHANISM NOTE for
     any slide that explains HOW the sweep is free: it is the EXACT
     bilinear analyzer basis (3 traces/arm span every analyzer angle,
     3.5e-10 vs direct) — NOT post-trace Jones at the detector, which
     is wrong at ~3% (field 2.8e-2 non-transverse there).  Delivery
     log at the foot of BRIEF_tg_demo.md.
 14. CTB: vortex EFC to 6.8e-15; band/polarization/vector-vortex layers.
 15. e2e6m: telescope + coronagraph families; closed loop holds the
     dark zone under drift.

V. Rodgers — AI-driven design
 16. The challenges: Mike's problems, the CC-driven method, ground rules.
 17. Rodgers 1 snapshot: his three designs reproduced at 0.99-1.00x his
     metric (the trust anchor).
 18. Rodgers 2 snapshot: the buildability finding — honest negatives.
 19. Rodgers 3: the family + solve discipline.
 20. Rodgers 3: the numbers — 45.4 nm vs his 53, honest rows.
 21. The product: family + driver -> adjacent problems, with or
     without CC.
 22. LIVE DEMO.  Lead-in framing (Dave's ruling 2026-08-28): the two
     beats carry DIFFERENT claims and different drivers — beat (a) is
     the PRODUCT in an engineer's hands (Dave solo, desktop MATLAB;
     the script's pauses/animation gate on the desktop and CC driving
     would collapse the product beat into a second AI beat), beat (b)
     is the designated CC-LIVE moment (AI drives the design tool).
     Say this contrast out loud when introducing the demo — it is the
     "with or without CC" line of slide 21 made concrete.
     (a) pol-PSI Twyman-Green DM gauge — build in visible
     Bench calls, fringes, actuator poke, analyzer sweep (no re-trace),
     4-step PSI map beside truth.  OPEN DECISION (Dave): keep demo
     beat 3 = the 11.7% finding live (the BS diattenuation rotates the
     test arm 7.5° off orthogonal; gauge reads 11.7% high while fringe
     contrast moves 0.17% — a 69× blind spot; fix solved in 5 traces)
     — TO and CC both recommend KEEP: it is the model-catches-what-
     the-instrument-can't beat; (b) rodgers3 adjacent-problem driver
     run on an audience spec (bounded to the driver's validated
     envelope).  Backup snapshots behind BOTH beats; rehearse twice;
     time-box 8 min.
 23. Live results recap (backups if needed).

VI. Close
 24. AI in design + analysis — the questions: verification discipline
     (graphics as gates, adversarial checks, honest rows), where the
     human rules, reproducibility.
 25. Summary + where the public code lives.

## Assets on disk (verified present)

- rodgers3: `mmacos/challenges/rodgers3/deck_rodgers3_V3.pptx` (+
  walk_for_review variant, deck_rodgers3.md/py sources).  DRAFT status
  — Dave's sign-off pass pending.
- CTB: `templates/30_instruments/bench_ctb/deck_ctb.pptx` (v5 content).
- e2e6m: `templates/80_end_to_end/e2e6m_r2/deck_e2e6m_r2.pptx` (27
  slides DRAFT) + `e2e6m/deck_e2e6m.pptx`.
- e2e worked example figures: templates/80_end_to_end + the s1-s4
  dirs; e2e2 deck (challenges-adjacent) STANDS.
- TG rig: `src/+macos/+design/twyman_green.m` +
  `templates/90_polarization/tg_psi_dm/` (DELIVERED 2026-08-28,
  commit 5d5f375: example, 7-beat demo script, backup PNG per beat,
  README with numbers; tTgPol 9/9 in SUITE_FAST).
- ~15 of 25 slides harvest existing material; the NEW builds are the
  history slide (Dave's imagery), slides 5 + 9 (content lists above),
  and the framing slides.

## Prep queue

1. TG pol-PSI demo build — DELIVERED 2026-08-28 (TO, res dev-candidate
   `5d5f375` LOCAL; delivery log at the foot of `BRIEF_tg_demo.md`).
   Remaining: Dave's dry run + the beat-3 keep/cut ruling (slide 22).
2. History-slide materials — DAVE (imagery, claims, sign-off).
3. Slide 22b: **DELIVERED 2026-08-28 (TO, res dev-candidate LOCAL,
   uncommitted; delivery log at the foot of
   `BRIEF_r3_adjacent_demo.md`).**  `oi_demo_step` + `run_oi_demo`
   batch runner + `demo_adjacent/` fallback bundle at **7, 12 AND 14°**
   + `REHEARSAL.md` + `tests/tOiDemoStep.m` (11 gates, unregistered).
   Measured, all PASS both gates: 7° 20.0 nm (pred 18.0), **12° 33.6 nm
   (pred 33.7 — 1.00×)**, 14° 51.2 nm (pred 54.9) — floors 93.8/24.9/
   24.9 mm, exit err ≤0.012°.  Deterministic to every printed digit.
   **Both flagged items RESOLVED:**
   (a) **RULED (Dave 2026-08-28): launch at the TOP of the session**
   (spec asked during the intro, slides 1–3; CC composes + launches
   VISIBLY on the projector — the driving is part of the show — then
   a status glance at ~slide 16, reveal at 22; ~40 min cover for the
   ~15-min solve).
   (b) resolved by TO — wide asks PASS (the walk line understates the
   envelope); `REHEARSAL.md` carries the corrected wording.  Do not
   use the old "honest deficit" line for interior widths.
   Still open: suite registration of `tOiDemoStep` (out on runtime).
   CC drives on demo day.
   **SESSION FOLDER: `macos/demo_session/`** — links to decks, both
   live-beat dirs, fallback bundle, briefs; `OUTLINE.md` there is
   DAVE'S EDITABLE session copy (sync structure changes back here).
4. Deck assembly under STYLE_REPORTS §5 — after 1-3; taskable, but
   history + slide 24 stay a Dave+CC pass.
5. `BRIEF_tg_ifo_v2.md` — **KICKED OFF 2026-08-29 (TO, time-boxed
   sprint):** cube PBS with a published MacNeille stack at 45°, no
   engine work, v2 lives BESIDE v1 with its own seven-beat demo
   script; v1 stays the rehearsed demo default, v2 joins demo day
   only if gated + rehearsed.  QUEUED post-demo:
   `BRIEF_dwd_plot_pagination.md` (TO), `BRIEF_ctb_segmentation.md`
   (campaign), e2e6m parked items.
