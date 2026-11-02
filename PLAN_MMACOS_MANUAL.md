# PLAN: the mmacos manual -- MATLAB + MACOS with an AI assistant (drafted 2026-10-23, CC for Dave)

**What it is.** A user's guide to mmacos for a reader working with an assistant.  Generated the way the MACOS
manual is: every listing, number and figure comes from a runner the reader re-runs (`make` regenerates them);
nothing is typed.  Each chapter has three voices: what the reader ASKS (the prompt), what the assistant READS
before it acts (the CLAUDE.md set, `optical_design/` rules and guides, `macos_f90/notes/`), and what it RUNS (a
named template or runner, with its transcript and figures).

**Dave 2026-10-26: defaults AGREED (e2e; the helpers kept).  Round 2 = `BRIEF_to_mmacos_manual_round2.md`: chapter 2 on the CTB as the telescope's back end, a later chapter for the coronagraph analyses (polarization impact on contrast, DM sensing and control), the regenerations and fixes ruled in PLAN 1.2 (bbb) (fff) (ggg) (hhh) (iii) (kkk), ch1's EFL label.**

**Decisions for Dave (RESOLVED 2026-10-26):** (1) the thread's telescope -- the e2e flow (the 4 m f/18 folded Korsch; all seven stages
and the e5hex1 budget exist; e2e2 is the 3 m f/20 design-only flow) or e2e6m (the 6 m unobscured coronagraph
front end); CC recommends e2e, and TO started on it 2026-10-23 as the default.  (2) the forty helpers no
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

**The design / tool split (CC's census of 2026-10-21, written down 2026-10-23 at TO's request; the index carries
no label).**  53 directories = 35 DESIGN + 17 TOOL + 1 planning.  README = has its own README.md today (30 of 53).

| Kind | Directories (README'd in **bold**) |
|---|---|
| design, 10_telescopes (15) | **freeform_unobscured**, **offset_imager**, rc_offaxis, rc_onaxis, rc_unobscured, sz_tma, tma_3plus1, **tma_centered**, tma_freeform, **tma_longslit**, tma_offaxis, **tma_onaxis**, **tma_unobscured**, tma_widefield, wf2_freeform |
| design, 20_segmentation (2) | **e5_pie**, **e5_seg** |
| design, 30_instruments (3) | **bench_ctb**, coro_planet_demo, **coro_walkthrough** |
| design, 40_benches (8) | bench_ifo, bench_ifo_dm, bench_layout, **cophi_dm96**, **pdi_dm96**, **tg_psi_dm96**, **tg_psi_dm96_oap**, **zwfs_dm96** |
| design, 80_end_to_end (4) | **e2e**, **e2e2**, **e2e6m**, e2e6m_r2 |
| design, 90_polarization (3) | bench_ifo_pol, **tg_psi_dm**, **tg_psi_dm_v2** |
| tool (17) | **design_layer_api**, **pupil_id** (10_telescopes: API examples and the pupil-as-surface method); grid_surface (20); **coro_experiments** (30); **dm_gauge_lib** (40); all eleven of 50_sensitivities (**e5hex1_budget**, **e5hex2_refzern**, **gen_segment_gridmat**, **zoom_5x5** README'd; e5hex1, run_dwdgrid_multi x3, run_dwdsurf_multi, run_dwdx_multi, run_dwdz_multi not); **view_rx_demo** (60) |
| planning (1) | **vsg_wip** (not runnable) |

So 20 design templates have a README and 15 do not; 10 tools have one and 7 do not.  Part II pages: the 35 designs
(README'd first, batch by family); the 17 tools as an appendix table; vsg_wip a line.  design_layer_api and pupil_id
are the T1 thread's entry points and get a page each in the tools appendix, not in the design catalog.

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
  TWO glasses N-BK7 + N-SF6, not free glass -- CCMac's correction from Yabe's verbatim statements) as the stress test;
  2010 green lens and 2021 down-under (winner Damien Gawron, EPD 198.4 mm, 90 deg) as the alternates; 2017 is a
  puzzle.  CCMac's record `cli_tests/records/iodc_problems_2026-10.md` (2026-10-23) confirms the prior and adds: the
  entries are CODE V / Zemax files, so `rx_converter` is on the comparison path for any scored year; the SPIE full
  text was unreachable from the Mac too (the JPL library pull is listed line by line in the record); any printed
  prescription is reproduced independently, never copied.  AWAITING Dave: name the problem.
  **Paper IDs for the JPL-library pull (all R. C. Juergens, SPIE proceedings; 2026-10-23):**
  | Year | Proc. SPIE | Paper | DOI |
  |---|---|---|---|
  | 2006 "the lens shuffler" (a 100 mm N-BK7 cylinder divided into elements) | 6342 | 63421U | 10.1117/12.695230 |
  | 2010 "the green lens" (one shared radius, N-BK7, 532 nm; winner EPD 81.3 mm, 43.5 deg, product 3537; 37 entries) | 7652 | 76522I | 10.1117/12.871174 |
  | 2014 "the Cinderella lens" (ten N-BK7 blanks 100 mm dia x 30 mm, 587.56 nm; winner Jon Ehrmann, 25 elements, EPD 33.9 mm, 62.5 deg, product 2119; 45 entries) | 9293 | 929303 | 10.1117/12.2074232 |
  | 2017 "the centennial lens" | 10590 | 1059002 | 10.1117/12.2299564 |
  | 2021 "the down under lens" (winner Damien Gawron, 198.4 mm, 90 deg, 17857) | 12078 | 1207820 | 10.1117/12.2605550 |
  | 2023 "the Yin-Yang lens" (N-BK7 + N-SF6; winner Arnaud Davenel) | 12798 | 1279837 | 10.1117/12.3012062 |
  Correction to the earlier abstract reading: the "81.3 mm / 43.5 deg" winner belongs to 2010, not 2014.
  **Papers read 2026-10-23 (`~/dev/MACOS_sandbox/IODC/`): the verified statements and the build are
  `BRIEF_to_iodc_challenge.md`.  The full texts change the order: 2010 green lens FIRST (its wavefront convention is
  MACOS's native chief-ray reference; one radius + a sign search; distortion <= 5 %), 2014 second (tilt-removed RMS,
  the blank-packing scorer, an intermediate-image winner), 2023 the stress test.  No prescriptions are printed in any
  paper: comparison is against the merit tables.**  **RULED 2026-10-26 (Dave): go with the revised order -- enter 2010 (the
  green lens) first as `mmacos/challenges/iodc2010/`, 2014 second, 2023 the stress test; TO after round 2 of the manual
  (`BRIEF_to_iodc_challenge.md`, the build section).**  **Step 1 DONE (TO, 2026-10-30, LOCAL): `mmacos/challenges/iodc2010/` -- the evaluator with its metric stated first (README), gated against an independent textbook ray trace to 2e-10 wave, the N-BK7 table row, the CODE V `.seq` for Luis, `tIodc2010` (PLAN_CONSOLIDATION 1.2 (hhhh)); the design search is next.**
- Korsch's TMA family and Cook's unobscured forms (the telescope templates' ancestry); Bauer / Schiesser / Rolland
  2018 freeform starting geometry (already in `mmacos/doc`): a freeform challenge with a published answer.
- Public prescriptions at real scale: HabEx and LUVOIR final reports (off-axis TMAs, segmented); Roman CGI's public
  PROPER model (already our comparison); JWST, Keck, GMT, TMT first-order designs and segment geometries
  (segmentation, edge sensing, MET at scale).

## Order of work
Part I from the existing stages (TO can start the prose once the telescope is chosen); the ImageDemo study and
camera template; Part II pages generated per template; Part III last.  The manual's own `make`, `src/`, markers and
`reference.docx` are the MACOS manual's, reused (`docs/mmacos-manual/`).
