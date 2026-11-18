# BRIEF to TO (macos-fd lane) — mmacos manual, round 2: Dave's rulings of 2026-10-26 (CC)

Dave is back; the standing push word has lapsed (pushes on his word again).  The defaults are AGREED (e2e; the
helpers kept).  His rulings on Part I, in the order the dependencies dictate -- the library fixes first, because
the chapters regenerate from them; every library change in your lane gated by the mmacos suites (and the pymacos
twin where one exists), committed by explicit path, with the number before and after in the message.

## A. Library fixes (resources), in this order
1. **(kkk) run with CHIEF.**  The sensitivity harvest (`run_sensitivities`, its `opd_ref` default) and
   `run_compare`'s engine truth both use the engine's chief-ray reference -- the one-reference rule of
   notes/opd_reference.md.  Gate: s6's agreement summary returns to the committed August class (x/z/grid w of
   2.6e-3 / 1.8e-3 / 0.29 or better), and tDwDx / the sens-core tests stay green (several pin the 10-10 'chief'
   re-pin already; say which moved and why).
2. **(fff) diagnose the zero gap-ray clip** on the live pie before touching stage 3: is the gap wide enough for the
   source grid (12128 rays live vs 12520 committed -- what changed the grid?), and is OBS ALL in force when the pie
   is built?  Report the cause with the number; fix only if the cause is a library defect, else record the setting
   the template must state.
3. **(hhh) + (iii)** in `run_met` / `met_view`: the metric chart prints every row (no log floor hiding the
   'edge+MET opt' bar), the hub inset's unit label reads metres.
4. **Ch1b: the design_report label** -- "plate-scale EFL at the field point" (72.32 m) beside the system EFL
   (72.14 m), one-word fix, your lane.

## B. Regenerate the e2e template's committed outputs (bbb) (ggg)
After A: rerun s1..s7 of `templates/80_end_to_end/e2e` and commit the committed products the stages write
(report .txt, .in, .mat, PNGs) so a user's rerun reproduces them -- stage 5 with the LIVE search (Dave: that is the
one to show), not the preset.  Say in the commit which numbers moved and by how much (the conics' 4th-5th figure,
the field-edge WFE, the M1 hole, the MET angles, the WEMs).  Then rerun the manual's ch01-ch07 drivers: every pin
either holds or is re-pinned with the mechanism named.

## C. Chapter 2 rewrite: the CTB as the back end (Dave)
The coronagraph behind the e2e telescope is the CTB (`templates/30_instruments/bench_ctb`, our best coronagraph
results), built as a Bench instrument at the telescope's relayed pupil.  Measure first: does bench_ctb's Bench take
the e2e telescope's exit beam (f-number, pupil size, the exit-pupil Return) without a new relay?  If it needs one,
that relay is the chapter's design step.  The Offner imager stays as the imager example; coro_walkthrough and
coro_planet_demo leave the chapter (they stay findings (ccc) (ddd)).  Same three voices, same loop.

## D. A new chapter after 7: coronagraph analyses on the CTB
What the CTB machinery already does, as generated runs: the contrast and dark zone; the polarization impact on
contrast (the CTB's pol legs, PLAN_POLARIZATION's validation ladder); DM sensing and control (the DM gauges of
40_benches as the sensing side, EFC as the control side).  Scope it from what runs today; anything that would need
new library code is a question for Dave, not a chapter.

## Then
Rerun Part II pages whose templates A-B changed (e2e, e5_pie, e5_seg, the sensitivities tools) -- mm_page tells you.
Report in the Part I form: REPORT_PART_I_round2.md.  Findings on the way -> PLAN 1.2 as before.
