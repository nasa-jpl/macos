# STATUS — demo eve (2026-08-31 evening), the resume-here file

**Talk: 2026-09-01 morning.  Deck: 38 slides, current and synced.**
Post-compaction: read this + the md header comment + memory
`project_keysight_demo` before touching anything.

## The deck workflow (all tools in this dir)
- Source of truth: `deck_keysight.md` + `deck_keysight.geo.json`
  (layout/fonts sidecar: x/y/w/h, cx/cy/cw, fs; block keys img:/txt:/code:N).
  Build: `python3 make_brief_slides.py deck_keysight.md` (LOCAL copy only —
  it applies the sidecar).  Render altered slides after every rebuild
  (soffice→pdftoppm, bump the tag; last used wF).
- Dave edits `deck_keysight_edit.pptx`.  **Sync ONLY via
  `./sync_edit_deck.sh`** — refuses on Impress lock or unrecovered diffs
  (`pptx_text_diff.py` + `pptx_geo_diff.py` vs `baseline_pass3.pptx`,
  FULL output, never head/tail-truncated); fold into md/sidecar, rebuild,
  then `--folded`.  `--force` only on Dave's explicit say-so.
- Figure tools: make_ifo_schematic_fig.py (slide 14), make_collision_fig.m
  (slide 13, engine truth), make_sens_figs.m (slides 22-24, r3 harvest
  ONLY), harvest_tel_sens.m (imaging-leg harvest — PULLED from the deck,
  kept as the EP-reset reproducer), make_ecosystem/make_mmacos figs.

## Today's major deck decisions (all Dave-ruled, committed)
slide 11 challenge/pursuit/achieved · slide 13 collision-in-the-part's-
plane figure (three-way layouts → backup w/ marker legend) · slide 14
TG-only schematic + terminal prompt panel + pixelated-camera answer ·
slides 15/16 post-run recaps (v2 numbers) · slide 20 "bench still needs
packaging work" stated forward · slides 22-24 control-basis families
(Rx/MonZern4/Grid4, ONE harvest, corners labeled coronagraph-clipped,
extraction idiom on slide) · fonts/layout per Dave's passes 2-3+.

## Open threads (in priority order)
1. **TO wall-slice handback** — six solves were finishing this morning;
   on handback: polished −10° numbers move slide-13's quoted values
   19–36% → re-cut ONLY with Dave's sign-off (slide currently quotes the
   committed deck truthfully).  §S4b.4 correction SIGNED OFF (in
   BRIEF_afocal4_wall.md foot) — TO commits it.
2. **BRIEF_afocal4_descent.md** — ready; hand to TO after wall handback.
   48 h box.
3. **e2e6m bench repose slice** (post-demo): union gate found bench-in-
   trunk-beam (−5.1 mm bare, Backend) + DM-pocket crossings (−13.1 mm
   both ways) — memory `project_e2e6m` has exact numbers; slide 20
   carries the plan.
4. **EP-reset arc (CCMac, queued)**: new reproducer — segment tilt on
   s3_imager_full.in gives one-signed domes (E1 Rx neg/pos 0.00, +9.8);
   Tz clean.  harvest_tel_sens.m + commits 3130e9e/7292e93.
5. Deck content checks still open from the fresh look: slide-5 C1
   "0.83–1.11×" unsourced; C3 band-vs-measured wording (slide 5 vs 6);
   the walk clearance spec (≥25 mm w/ 1.5 mm hinge) unstated on 8/23.

## Demo day
- IFO demo: RUNBOOK_ifo_demo.md (fresh CC session, ./ifo_dialog.sh, v2).
- Telescope solve: oi_demo_step(width) per slide 9 / rehearsal bundle.
- All numbers in the deck re-derivable from committed records; renders in
  renders/ (wF latest for 22-24).

## IN-FLIGHT at compaction (2026-08-31 late)
- **Slide 10 table change (Dave)**: replace the "(chief)" column with
  MACOS scored at PISTON+TILT REMOVED (the rung-2-style convention).
  Numbers: look in MACOS_res_dev/mmacos/challenges/rodgers1/ (PACKET.md
  + addenda) for per-design piston+tip/tilt-removed RMS; if absent,
  re-score the committed .in decks (S2 verbatim, S3/S4 solves +
  re-traced) — box max/avg over the field map, piston+tilt removed.
  Then: md table edit, rebuild, render p10, gated sync, commit.
- **cf3d deep-dig ladder RUNNING detached** (d=1.10 apl restart ladder,
  monitor bq69ackrd, report cf3d_report.txt, checkpoint cf3d_run.mat,
  wall 6 h; EXTENDED by watcher on completion: resume w/ max_rounds 40, wall 9 h, target 5e-11 -- Dave 2026-08-31 late).  Round 5 = 1.53e-8.  Morning:
  slide-25 re-cut with Dave's sign-off if it digs.
- Dave flag pending: "M2-M4" wording on slides 22/23 (no M4 element).
