# MACOS Manual Project

> **Post-compaction / post-upgrade — re-read the docs first.** After a
> context compaction, before resuming manual work, re-read this file,
> the root `CLAUDE.md`, and `README.md` here.

**As of 2026-07-06 the manual is Markdown-first.**  Canonical source:
`src/*.md` + `src/media/`.  The 2026-04/05 docx patch chain
(`patches/`, `scripts/`) was REMOVED from HEAD on 2026-07-07 (it was
frozen provenance, unreferenced by the build, and its files carried
AppScan findings) — recover it from git history if ever needed.  The
`macosMan*.docx` snapshots remain as provenance.  New manual work
edits `src/` only.

## Cycle

```bash
make            # docx + html into build/
make pdf        # via LibreOffice (tools/docx2pdf.py refreshes TOC)
make appendix-c # regen src/93_appendix_c_commands.md from macos_help.inc
make transcripts # regen src/transcripts/ (example listings + CLI dialogs; MACOS_BIN= picks the binary)
make cmdref-regen # regen cmdref Parts I-III + the CLI-suite coverage appendix
                  # (every Part I entry's "*CLI suite:*" line comes from cli_tests/cmd/*.jou)
```

## Rules / gotchas

- **Appendix C is generated** from `macos_f90/macos_help.inc` by
  `tools/gen_appendix_c.py`. Never hand-edit
  `src/93_appendix_c_commands.md`; when a command is added to
  `macos_cmd_loop.inc` + `macos_help.inc`, run `make appendix-c`.
- **Listings and CLI dialogs are generated** (2026-10-19, Dave: "regenerate ALL of
  them ... make it repeatable").  `tools/regen_transcripts.py` writes
  `src/transcripts/listing_<file>.md` (each `examples/` file verbatim) and
  `src/transcripts/<name>.md` (dialogs from `tools/transcripts/*.jou`, the cmd-gate
  grammar plus `# TRANSCRIPT: <name> [+device]`, `# TRANSCRIPT-OFF`,
  `# TRANSCRIPT-BANNER:`); a chapter keeps only the marker line
  `<!-- transcript: NAME -->` (`tools/include_transcripts.lua`; a missing file stops
  the build).  Never hand-edit `src/transcripts/` or type a dialog into a chapter:
  change the journal and rerun.  Every journal pins with `# EXPECT:` (5 figures, bounds
  near zero; a journal with none is refused, a miss writes nothing) and
  `cli_tests/run_cli_tests.sh transcripts` gates both compilers.  Transcripts come from
  the ifx build of the clean worktree.  What moved against 3.2, and the open decisions,
  are in `TRANSCRIPT_DIFF.md`.  Section 4's NEW tutorials are generated too (2026-10-19): a
  `# TRANSCRIPT:` switch between two answers of one NEW cuts the dialog there (`from=` / `upto=`
  regexes); author a long dialog with `tools/prompt_walk.py`, which prints each prompt beside the
  answer it got.  Nothing in `src/` is a typed dialog any more.
- **No TOC in source** — pandoc `--toc` builds it. A plain
  `soffice --convert-to pdf` leaves the Word TOC field EMPTY; the
  Makefile pdf target uses `tools/docx2pdf.py` (UNO: update indexes
  twice, then export) for populated page numbers.
- **Word styling** comes from `--reference-doc=reference.docx` (pandoc
  default restyled). GOTCHA: macosMan4_0_styled.docx is DEFECTIVE as a
  reference-doc (tables collapse to one column) — never switch back.
  Pandoc maps to its own style names (BodyText, SourceCode,
  FirstParagraph, Compact, Heading1-4); restyle by editing those
  styles in the reference docx, not by touching src/.
- **Section order** is the lexical sort of `src/[0-9]*.md`; new
  appendix files slot in at 92_/94_ etc.
- Terminal transcripts = generated (above); other code = indented or fenced blocks; figures = inline
  `![](media/imageNN.png)` with the `**FIGURE n**` caption paragraph
  adjacent. Image `{width=...}` attrs are in fractional inches from
  the docx — keep unless intentionally resizing.
- The migration scripts (`tools/style_map.lua`,
  `tools/split_sections.py`, `tools/convert_from_docx.sh`) are
  one-time; re-running convert_from_docx.sh OVERWRITES src/ edits.
- Legacy back-of-book index was dropped at migration
  (`src/_dropped_legacy_index.md` keeps a copy; stale v3.2 page
  numbers, superseded by --toc + PDF bookmarks + HTML search).
- Known debt: ~73 equations are images/Symbol-font fragments to be
  re-entered as LaTeX math ($...$ renders in HTML; docx gets OMML).
- `RECONCILE_4_01.md` documents folding in the Writer fork
  (`macosMan4_01.docx`); its "(historical)" version-notes sections now
  live at the end of `src/01_introduction.md`.
- New-feature documentation for engine work (FreeForm, SLSQP, XPS,
  COMPOSE, ...) belongs in the matching `src/` section; add a bullet
  to "New Features" in `src/01_introduction.md` per release.

## Product names in prose (Dave 2026-10-20)
`MACOS`, `SMACOS`, *mmacos* and *pymacos* are set so a scanning reader catches them: MACOS and SMACOS in
capitals as written, *mmacos* and *pymacos* in italics in body prose (not inside code spans or paths).  Applied
in the Introduction; apply it to any prose you write or regenerate in the manual and the command reference.
