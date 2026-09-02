# CCMac — welcome back (2026-09-01 evening, from Dave via CC)

**The one hard constraint first: the Keysight CODE V demo is TOMORROW
MORNING (Sep 2).**  Tonight the box is carrying live compute and the
demo deck is frozen-ish; read the box rules before running anything.

## Where the world moved while you were out (verify in the records, not here)

- **Keysight demo**: `macos/demo_session/` — deck 42 slides
  (`deck_keysight.md` + geo sidecar = source of truth;
  `STATUS_demo_eve.md` = the resume-here file).  Two live elements
  tomorrow: an audience-spec'd telescope solve (`oi_demo_step`) and the
  pol-PSI Twyman–Green DM gauge (RUNBOOK_ifo_demo.md).  **Do not touch
  `deck_keysight_edit.pptx`** — sync is gated (`sync_edit_deck.sh`).
- **e2e6m_r2 CF campaign** (`templates/80_end_to_end/e2e6m_r2/`,
  `e2e6m_r2_LOG.md` foot): the EFC **restart ladder** dug the d=1.10
  apodized train 1.2e-6 → 1.133e-9 (CF3d); **CF5b** re-posed the drift
  series at that dug state (JWST-class 10 nm/hr, 24 h — loop holds
  2e-9); **CF3e closed tonight**: unapodized floor 4.48e-6 at thru
  0.745 → the apodizer buys ~4000× contrast for 8.2× light, and the
  released-la work shows the substrate is apodizer-independent — the
  prolate buys the STROKE PRICE of collecting it.  New controller
  machinery in `cf_efc_lib`: `push` (stroke-released walk) and `bump`
  (beta-bump).  CF3f (Gaussian-feather apodizer) runs overnight.
- **afocal4** (`challenges/afocal4/`, RESULTS.md): the DESCENT ruled
  71 nm unreachable coaxially at any N; the OFFAXIS slice (running
  under TO) seeded off-axis Mersennes — field curvature owns 60–99.7%
  of wavefront variance, an off-axis FOUR-mirror at 2181 nm beats the
  coaxial seven-mirror floor (M-defense pending), N=6 crossover test
  in flight.
- **Branch state**: everything above is LOCAL on dev-candidate
  (both repos), **nothing pushed** — do not push, pull, or rebase.

## Tonight, if you have cycles (both bounded, both de-risk the demo)

1. **Demo smoke test (the useful one).**  Dry-run both live elements
   end-to-end tonight: (a) `oi_demo_step` at 12° — must reproduce the
   rehearsal bundle "to every printed digit"
   (`templates/10_telescopes/offset_imager/demo_bundle/`, backup
   slide 37's claim); (b) walk `RUNBOOK_ifo_demo.md` (fresh CC session,
   `./ifo_dialog.sh`) far enough to confirm the rig runs.  Report
   PASS/FAIL + exact digits to Dave — fix nothing without his word.
2. **Slide-5 source hunt (read-only).**  Slide 5 quotes challenge-1
   reproduction "0.83–1.11×" with no committed source located, and the
   C3 band wording (slide 5 vs 6) needs a records check
   (`challenges/rodgers1/PACKET.md` + addenda, rodgers3 records).
   Deliver the citations (or "not in the record") — the deck edit, if
   any, is CC's through the gated pipeline.

## Box rules tonight (contention killed jobs at 09:11 today)

- 4+ MATLABs are live (CF3e/CF3f ladders + TO's off-axis fleet).  No
  fleets, no model-512 Jacobian sweeps; one MATLAB at a time, and
  prefer `-batch` with a timeout.  If something of yours dies
  silently, suspect seats/memory, not your code.
- No engine rebuilds tonight (the demo runs on the built tree).

## Your September queue (start after the demo, order per Dave)

1. **EP-reset arc** (yours, queued): frozen-EP vs family `reset_xp`
   asymmetry; reproducer = `demo_session/harvest_tel_sens.m` (segment
   tilt on s3_imager_full.in → one-signed domes, E1 Rx neg/pos 0.00 /
   +9.8; Tz clean).  Commits 3130e9e / 7292e93.  Default + regen
   policy = Dave's call — bring him options, not a fait accompli.
2. **Model-transition crash, REOPENED**: 128→512 `macos_init_all`
   heap crash is back and deterministic; full pymacos suite = 308
   tests; fix path = PLAN §0 reopen.
3. VSG2 Stage-3 pupil-reimage (still open from August).

Standing rules unchanged: engine truth over .in text-parsing; write
resolutions at resolution time; meaningful tests only; only push when
Dave asks.
