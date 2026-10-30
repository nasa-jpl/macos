# BRIEF to TO — the manual's listings and transcripts become GENERATED text (CC for Dave, 2026-10-19)

**Why.** The manual (`docs/macos-manual/src/*.md`, Markdown-first since 2026-10) still carries the 3.2-era
listings and CLI dialogs as hand-pasted text, many of them miscast as pipe tables by the PDF extraction:
Appendix A alone has 1061 table rows that are really the example decks and journals; sections 3-8 hold
dialogs quoted from an engine twenty years old.  Today's finding: the ORS example ran for a day on a deck
whose extracted listing had a Zernike coefficient split at a page margin (PLAN_CONSOLIDATION 1.2 (cc)).
Dave 2026-10-19: regenerate ALL of them with the current engine, and make it repeatable.

**The deliverable.** `make transcripts` in `docs/macos-manual/` regenerates every listing and dialog from
sources, and the chapter sources include the generated files.  After it lands, no listing or dialog in
`src/` is typed by hand.

**Start from a cleared session.  Read in this order:** root `CLAUDE.md`, `docs/macos-manual/CLAUDE.md`,
`cli_tests/cli_cmd_gate.py` (the pty driver and the journal grammar: `# COPY:`, `$TMP`, `$DECK`, `# EXPECT:`),
`docs/macos-manual/tools/regen_engine_figs.py` and `tools/figs/FIG31.jou` (the figure pipeline you are
twinning), `docs/macos-manual/Makefile`, then the ORS/SRS paragraphs of `src/06_diffraction_analysis.md`
(a dialog re-emitted by hand this morning -- the target form).

**Mechanism (build it in this order; commit each step with its gate green).**
1. `tools/transcripts/<name>.jou` -- one journal per quoted dialog, the cmd-gate grammar, run by a new
   `tools/regen_transcripts.py <macos binary>` that writes `src/transcripts/<name>.md`: the dialog as the
   manual prints it (`MACOS\>**cmd**` line, then the engine's lines in a 4-space code block, the prompts
   with the answers as the pty saw them).  Strip the `Graphics device` prompt and `/null` as the FIG
   journals do; keep the engine's spacing.
2. The example decks are included VERBATIM from `docs/macos-manual/examples/<name>.in` (as written, not the
   SAVE-emitted canonical form: the appendix is a tutorial and the emitted form reorders keys), and the
   example journals from `examples/<name>.jou`.  The generator writes `src/transcripts/listing_<name>.md`.
3. A pandoc include step (the `number_sections.lua` filter already handles RawBlocks; add an include filter
   or a pre-pass in the Makefile) replaces a marker line `<!-- transcript: NAME -->` in a chapter source
   with the generated file.  The marker is the only thing a chapter keeps.
4. Appendix A first (`src/91_appendix_a_examples.md`): every listing and journal there becomes a marker.
   Diff the regenerated text against the 3.2 PDF (`pdftotext -layout` of `docs/macosMan3.2.pdf` is in CC's
   scratch; regenerate it) and log every NUMBER that moved in `docs/macos-manual/TRANSCRIPT_DIFF.md` with
   one line of cause each -- that log is the editorial input for Dave.
5. Then the chapter dialogs: 3 (user interface), 5 (ray-trace analysis; RRAY/PRAY/TPR, OPD, SPOT, FEX),
   6 (ORS/SRS/BEAM/INT), 7 (imaging), 8 (differential ray tracing), 11 (polarization).  One journal each.
6. `make all pdf` and the edit copies (`edits/macosMan_{edit,baseline}.docx` and `docs/macosMan.pdf` are
   copies of the build outputs) at the END of the series, not per step.

**Rules.**
- Every journal is also a gate: it carries `# EXPECT:` lines pinned to 5 figures (rule 2) or bounds for
  near-zero quantities (rule 1; PLAN section 4 item 4), and `run_cli_tests.sh cmd` runs it.  A dialog
  that cannot be pinned is a finding, not a transcript.
- An example that no longer illustrates its point is NOT silently replaced: name it in TRANSCRIPT_DIFF.md
  with what it shows now; Dave chooses the replacement (ImageDemo and the seg.in ORS example were two).
- A deck extracted from the PDF is checked for numbers split at the margin before it is trusted:
  `grep -nE '[DE][-+]$' examples/*.in` must print nothing.
- Binaries from the clean worktree `~/dev/macos_cli_base` at a committed SHA, both compilers for the
  gate; the transcripts themselves from ifx (the manual prints one engine).
- Fixed-form rules, explicit-path commits, no `--amend`, no `commit -a`; push only on Dave's word.

**Done when:** `grep -c '^|' src/91_appendix_a_examples.md` is the count of REAL tables (expect 0),
`make transcripts && make all pdf` runs clean from a fresh clone, every journal is green on both
compilers, TRANSCRIPT_DIFF.md is complete, and Dave has the replacement list.
