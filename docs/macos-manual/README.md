# MACOS Manual

The manual is maintained as **Markdown source** in `src/` — one file
per section, figures in `src/media/`.  PDF, HTML, and DOCX are build
products; edit the Markdown, never the outputs.

```
macos-manual/
├── src/                        ← THE manual (edit these)
│   ├── 00_frontmatter.md
│   ├── 01_introduction.md
│   ├── ...
│   ├── 91_appendix_a_examples.md
│   ├── 93_appendix_c_commands.md   ← GENERATED from macos_help.inc
│   └── media/                  ← figures (PNG)
├── Makefile                    ← build driver
├── build/                      ← outputs (gitignored)
├── tools/                      ← build + migration scripts
├── examples/                   ← .in/.jou files shown in Appendix A
└── (legacy .docx/.pdf, patches/, scripts/ — see below)
```

## Building

```bash
make            # build/macosMan.docx + build/macosMan.html
make pdf        # build/macosMan.pdf (docx -> LibreOffice headless,
                #  TOC page numbers refreshed via tools/docx2pdf.py)
make appendix-c # regenerate src/93_appendix_c_commands.md from
                #  macos_f90/macos_help.inc  (never hand-edit it)
```

Requirements: pandoc, GNU make, LibreOffice (PDF only), python3.
Section order = lexical order of `src/[0-9]*.md`.  The Table of
Contents is generated at build time (`--toc`) — there is no TOC in the
source.  Word styling comes from `--reference-doc=reference.docx` (pandoc default restyled: Arial headings, Times body). NOTE: do NOT use macosMan4_0_styled.docx as reference-doc — it collapses all tables to single-column.

## Editing guidelines

- Prose: plain Markdown paragraphs.  Terminal transcripts and .in file
  excerpts: fenced code blocks.  Figures: `![](media/imageNN.png)`
  placed inline where they belong; keep the `**FIGURE n** caption`
  paragraph next to the image.
- `src/93_appendix_c_commands.md` is generated — run `make appendix-c`
  after adding a command to `macos_help.inc`.
- Known debt inherited from the PDF extraction: ~73 equations are
  low-resolution images or Symbol-font fragments awaiting re-entry as
  LaTeX math; some diagram label fragments remain near figures.

## History / legacy files

The manual originated in FrameMaker (source lost), survived as
`docs/macosMan3.2.pdf`, was Acrobat-converted to docx, restyled, and
then maintained 2026-04/05 as a chain of Python XML patches
(`patches/`, `scripts/`) producing `macosMan4_0_styled.docx` →
`macosMan4.2beta.docx`.  On 2026-07-06 that endpoint was converted to
the Markdown source in `src/` (see `tools/convert_from_docx.sh`) and
Markdown became canonical.  The patch chain and docx snapshots are
retained for provenance only.  `RECONCILE_4_01.md` records how the
parallel Writer fork `macosMan4_01.docx` was folded in.
