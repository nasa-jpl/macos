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
```

## Rules / gotchas

- **Appendix C is generated** from `macos_f90/macos_help.inc` by
  `tools/gen_appendix_c.py`. Never hand-edit
  `src/93_appendix_c_commands.md`; when a command is added to
  `macos_cmd_loop.inc` + `macos_help.inc`, run `make appendix-c`.
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
- Terminal transcripts = fenced code blocks; figures = inline
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
