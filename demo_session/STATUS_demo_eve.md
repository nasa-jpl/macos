# STATUS — demo eve (2026-08-31 evening), the resume-here file

**Talk: 2026-09-02 morning.  Deck: 40 slides, current and synced
(ladder slide 28 added late 2026-08-31; all decks match the md).**
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
- AUDIENCE Q&A RELAY (Dave 2026-09-01): slide 32 = "The agent's own
  questions for this room" (4 questions, each anchored in a committed
  number: ladder/step machinery, apodizer stroke economy, afocal pupil
  requirement + compressed-beam corrector, conventions as data).
  Protocol: Dave pastes audience answers/questions into THIS supervisor
  session; the agent answers from the records (concise, numbers cited,
  no claim beyond the committed record) or asks follow-ups.  Slide
  numbering after 31 shifted +1 (discussion 33, closing 34, backup 35).
- All numbers in the deck re-derivable from committed records; renders in
  renders/ (wF latest for 22-24).

## IN-FLIGHT (2026-08-31 ~22:00)
- Slide 10 table change: DONE (piston+tilt-removed MACOS column in md).
- **cf3d ladder DONE + FOLDED** (Dave's "redo the e2e6m_r2 reports and
  the deck_keysight* decks" executed): 1.215e-6 → 1.133e-9, 10 digging
  rounds / 4.6 h, plateau PROVEN by 3 no-step rounds (alphas to 1e-10,
  niter 30); la(G) 3.52e-11 at the dug state — the FALCO case.  Record:
  CF3d LOG section + cf3d_* committed (MACOS_res_dev).  Deck: NEW slide
  after the drift slide ("Digging the dark hole: the overnight restart
  ladder", figs/fig_cf3d_ladder.png), drift slide's honest-gap bullet
  re-scoped, closing-slide claim quantified.  Rebuilt 40 slides,
  renders verified (wR).  Edit-deck sync COMPLETE (~00:00 09-01): Dave
  closed with zero diffs; gate passed clean; edit deck + baseline carry
  the 40-slide ladder build.  ALL deck_keysight* decks current.
- "M2-M4" wording: DONE -> M2-M3 (both state-equations slides, 09-01).
