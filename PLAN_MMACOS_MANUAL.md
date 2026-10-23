# PLAN: the mmacos manual -- MATLAB + MACOS with an AI assistant (drafted 2026-10-23, CC for Dave)

**What it is.** A user's guide to mmacos for a reader working with an assistant.  Generated the way the MACOS
manual is: every listing, number and figure comes from a runner the reader re-runs (`make` regenerates them);
nothing is typed.  Each chapter has three voices: what the reader ASKS (the prompt), what the assistant READS
before it acts (the CLAUDE.md set, `optical_design/` rules and guides, `macos_f90/notes/`), and what it RUNS (a
named template or runner, with its transcript and figures).

**Decisions for Dave:** (1) the thread's telescope -- the e2e flow (3 m f/20 TMA, all seven stages, the e5hex1
budget exist) or e2e6m (the 6 m unobscured coronagraph front end); CC recommends e2e.  (2) the forty helpers no
template calls (`design/src` 40 of 75, `sensitivities` 12 of 24 unused): template pages, or retired before
Part II is written -- a catalog of dead functions is worse than none.

## Part I -- the worked thread: one telescope from sketch to error budget

| Ch. | The reader asks for | What runs (exists today) |
|---|---|---|
| 1 Telescope | a telescope from a spec (aperture, field, f-number, obscured or not); the .in file emitted; first-order properties; the pupil found | design layer `Telescope`; `e2e/s1_telescope.m`; the `10_telescopes` family as the answers; `view_rx` |
| 2 Instruments and benches | an OBSERVATORY-CLASS telescope (the TMA of ch. 1) feeding SEVERAL instruments -- an imager, a spectrograph, a coronagraph -- each added with `Bench` at the telescope's focal surface or relayed pupil (the station/standoff rules, folds, relays, the exit-pupil Return); a BENCH on its own (the DM gauges) -- one chapter, the same machinery | `e2e/s2_instrument.m`; `Bench`, `add_fold`, `add_pupil`; `tma_longslit` + `dyson5` (spectrograph), `tma_3plus1` / `coro_walkthrough` / `e2e6m` (coronagraph), `offset_imager` (imager); `40_benches` |
| 3 Segmentation | hex or pie, gaps, edge sensors; SegMirMaker | `s3_segmentation.m`; `e5_pie`, `e5_seg`; `run_segmentation.m` |
| 4 Sensitivities | dw/dx, dw/dz, dw/dgrid, dw/dsurf, multi-field; the OPD conventions | `s4_jacobians.m`; `50_sensitivities`; `doc/SENSITIVITY_TOOLS.md`, `doc/opd_conventions.md` |
| 5 Metrology | laser MET and the truss; the dmet/dx layout optimizer | `s5_met.m`; `run_met.m` |
| 6 Open and closed loop | the time series, the controller, the two-pass rule; compare | `s6_compare.m`, `s7_simulate.m`; `run_simulator.m`, `run_compare.m` |
| 7 Error budget | the coefficients, the WFC and metrology stages, the workbook handed to systems engineering | `run_budget.m`; `50_sensitivities/e5hex1_budget` (the .xlsx) |

**NOTE (Dave 2026-10-23): study `ImageDemo.in` WITHOUT its aberrations as the basis of a design family and
template -- the SMALL DEDICATED IMAGER.**  The Mars Observer Camera is a compact, undersampled, single-purpose
camera: a Cassegrain with a refractive field-corrector relay, the exit-pupil Return and one image plane.  It is NOT
an observatory-class front end feeding several instruments (that is chapter 2's thread, on the TMA), so it is its
own family: the camera template of PLAN_CONSOLIDATION section 4 item 6, a Part II page beside the telescopes, and
chapter 2's counter-example -- what a dedicated imager looks like when the telescope and the instrument are one
design.  Before it becomes a template: strip the Zernike terms (the 2026-10-19 rejoined line showed how sensitive
the deck is to them), check the relay's glasses and apertures, fit its first-order properties (the manual's Section
6 numbers) and its sampling (undersampled by design -- say by how much) to the design.

## Part II -- the template catalog (what we have done)
One page per DESIGN template (35, in the index's order): what it is, the knobs, the figure, the number it
reproduces, and the prompt that makes an assistant adapt it.  The five CHALLENGES get a chapter each, the scoring
convention stated before the number (challenges/README.md's rule).  The 17 tool templates: an appendix table.
Pages are generated from each template's own run (the `make transcripts` pattern); a template with no README gets
one first (23 of 53 have none).

## Part III -- reference
The api catalog = cmdref Part II (generated; cross-reference, do not copy).  The design rules and agent guides
(`optical_design/`) ship on main after the exclusions change; cross-reference.  One chapter the MACOS manual lacks:
how the assistant is expected to work with this library -- Appendix D's content where the MATLAB reader looks, plus
the house rules that bit us (measure before guarding; one parameterized runner per system; no `exit` in examples;
deck figures unmodified; the pin policy).

## Challenges from the literature (candidates; CC 2026-10-23)
- The IODC (ILDC until 1994) lens design problems, one per conference since 1990, entrants' solutions published:
  scored, public, comparable -- the challenge form the README defines.  Most are REFRACTIVE camera lenses: they
  would exercise the glass tables, aspheres, multi-wavelength and the CodeV converter rather than the mirror paths
  the library is built around.  Choose by reading the problem statements (SPIE proceedings), not from memory:
  the candidates, the criteria and the extraction list are `NOTE_to_ccmac_iodc.md` (2026-10-23): 2014 Cinderella
  (ten N-BK7 blanks, fixed glass, monochromatic, spherical, 45 entries) to enter first; 2023 Yin-Yang (two wavelengths,
  free glass) as the stress test; 2010 green lens and 2021 down-under as the alternates; 2017 is a puzzle.
- Korsch's TMA family and Cook's unobscured forms (the telescope templates' ancestry); Bauer / Schiesser / Rolland
  2018 freeform starting geometry (already in `mmacos/doc`): a freeform challenge with a published answer.
- Public prescriptions at real scale: HabEx and LUVOIR final reports (off-axis TMAs, segmented); Roman CGI's public
  PROPER model (already our comparison); JWST, Keck, GMT, TMT first-order designs and segment geometries
  (segmentation, edge sensing, MET at scale).

## Order of work
Part I from the existing stages (TO can start the prose once the telescope is chosen); the ImageDemo study and
camera template; Part II pages generated per template; Part III last.  The manual's own `make`, `src/`, markers and
`reference.docx` are the MACOS manual's, reused (`docs/mmacos-manual/`).
