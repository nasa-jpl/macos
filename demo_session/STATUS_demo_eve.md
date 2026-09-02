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

## DEMO-DAY MORNING (written 2026-09-02 ~00:45 — READ THIS FIRST)

**Deck: FINAL. 42 slides, DRAFT retired, export-review marking on
Dave's 17 slides (builder-encoded: EXPORT_MARK_TITLES in
make_brief_slides.py). All deck_keysight* synced via the gate and
PUSHED (nasa-jpl dev-candidate on BOTH repos — pushing is allowed
now, rule Bash(git push*)).** Slide count 42: ladder 27, stations 28,
drift 29, reveal 30, working-questions 31, agent-questions 32,
discussion 33, closing 34, backup 35+.

- **Choreography RULING: live solve on the LINUX box** — slide 30 says
  "in 15 minutes" and the ThinkPad timing run confirmed ~15 min class
  (scratchpad/oi_timing.log, grep THINKPAD TIMING for the final
  number). The Mac (CCMac) = hot spare: 18× faster, byte-identical
  digits — spoken bonus at the reveal, not the driver.
- **Q&A relay live**: slide 32 = five agent questions (afocal pupil +
  compressed-beam corrector; basins/seeding; merit scaling;
  conventions-as-data via macro/API; EFC ladder + stroke economy).
  Dave pastes audience answers/questions into the SUPERVISOR session;
  answer ONLY from committed records, numbers cited.
- **Key quotables**: TO's slice COMPLETE — off-axis beats coaxial at
  ALL FOUR mirror counts; best = SIX mirrors at 33× (2344.5 nm, 3.54×
  gain); mirror economy = four at 42× vs best coaxial 109×; the
  decenter effect is measured/controlled but UNEXPLAINED (4 hypotheses
  dead — slide 15 asks the room); rule 39 = freedom flows to the
  unscored quantity (RESULTS § OFFAXIS wrap, 26 commits local
  unpushed).  DECK now 45 slides: agent questions DISTRIBUTED (9
  basins/merit, 15 afocal+decenter, 32 EFC 'For coronagraph folks',
  35 conventions closer); reveal 33, discussion 36, closing 37,
  backup 38+.  CCMac σ-sweep: σ1 converged 8.8e-5 (too soft), σ3
  wall-limited 1.09e-5 still digging (RESUMED in parallel per Dave),
  σ4 running; depth axis = la@50/laU, c_end wall-limited rows =
  floor-so-far. Apodizer trade: prolate 0.091 thru
  / 1.13e-9 floor; feather σ2 0.607 / **9.527e-7 FINAL** (CONVERGED
  00:17 Sep 2, 22 rounds, laU 1.155e-11 @ 6.3 µm); bare 0.745 /
  4.48e-6 (CLOSED, 18 rounds incl. bump-armed resume 16–18 all
  no-gain, laU 4.4e-11 @ 5.4 µm — linear validity is the
  wall). CF5b: JWST-drift hold at 2e-9 from the dug state.

REVEAL PROTOCOL (dry-run PASS 2026-09-02, 5 s end-to-end): when the
live `oi_demo_step(W)` wraps, run `python3 make_reveal.py <W>` in
demo_session/ — it finds the NEWEST `oi_demo_<W>deg*_verdict.txt`
bundle (live runs are timestamped), crops the figures, splices the
live numbers into slide 33 of a COPY (`deck_keysight_live.pptx`),
and renders `renders/live33-33.png` for pre-flight.  Verified: text
+ geo diffs touch ONLY slide 33.  Dave's ruling: he opens the live
copy on my prompt at the reveal; the presented deck is never edited.
Fallback = the 12° bundle already on the slide.

MORNING HARVEST (in order):
1. cf3f_report.txt tail → final feather floor; write the three-way
   trade table + LOG section; commit + push (MACOS_res_dev).
2. CCMac's cf3h_s* σ-sweep (GATE: their chain must reproduce tag
   f20 / THRU 0.607 / static 9.976e-5 — gate failed ⇒ sweep is noise).
3. TO's § OFFAXIS wrap + N=7 pupil floor + N=6 crossover.
4. oi_timing.log final line → confirm/record the cover time.
5. Monitors: cf3f = b48fbgvkd (keep until DONE); cf3e = b9hse8g8e
   (arm CLOSED at 4.480e-6 — stop the monitor).
6. Deck FINAL v2 (2026-09-02 morning): Dave 08:50 edit pass folded (coordinate-free line, mmacos toolbox retitle, barrier wording, slide-34 sub-bulleted answers) + Marx-informed control≡truth sub-bullet on slide 32; live copy regenerated; gate re-baselined.
