# MACOS Manual

The manual is maintained as **Markdown source** in `src/` — one file per
section, figures in `src/media/`.  PDF, HTML, and DOCX are build
products.  **Edit the Markdown; never edit the outputs.**

```
macos-manual/
├── src/                            ← THE manual (edit these)
│   ├── 00_frontmatter.md           ← title page, contact, authorship
│   ├── 01_introduction.md          ← incl. "New Features" per release
│   ├── 02_technical_overview.md
│   ├── 03_user_interface.md
│   ├── 04_describing_optical_systems.md
│   ├── 05_ray_trace_analysis.md
│   ├── 06_diffraction_analysis.md
│   ├── 07_beam_propagation_imaging.md
│   ├── 08_differential_ray_tracing.md
│   ├── 09_subroutine_macos.md
│   ├── 90_references.md
│   ├── 91_appendix_a_examples.md
│   ├── 93_appendix_c_commands.md   ← GENERATED — never hand-edit
│   └── media/                      ← figures (PNG/JPEG)
├── Makefile                        ← build driver
├── build/                          ← outputs (gitignored)
├── reference.docx                  ← Word/PDF styling (see Gotchas)
├── tools/                          ← build + one-time migration scripts
├── examples/                       ← .in/.jou files shown in Appendix A
├── FIGURE_RESCUE_LOG.md            ← figure-migration status + audit
├── RECONCILE_4_01.md               ← how the Writer fork was folded in
└── (legacy: macosMan*.docx snapshots — frozen provenance)
```

This directory also produces a second document, the **MACOS Command
Reference** (`cmdref/`) — a standalone catalog of every engine command
and every mmacos/pymacos function, written to guide new users of the
public repo.  See "Command Reference" below.

## Building

```bash
cd docs/macos-manual
make            # build/macosMan.docx + build/macosMan.html
make pdf        # build/macosMan.pdf (via headless LibreOffice)
make appendix-c # regenerate Appendix C from macos_f90/macos_help.inc
make cmdref     # build/macosCmdRef.docx + .html
make cmdref-pdf # build/macosCmdRef.pdf
make clean
```

Requirements: `pandoc`, GNU `make`, `python3`; LibreOffice only for the
PDF leg.  Section order is simply the lexical order of `src/[0-9]*.md`.
The Table of Contents is generated at build time (`--toc`) — there is
no TOC in the source, and no index (PDF bookmarks + HTML search
replace it).

For a quick preview while writing, build the HTML — it's the fastest
target and self-contained (open `build/macosMan.html` in a browser).

## How to add or edit content

### Everyday edits

Open the `src/` file for the section, edit, `make`, look at the HTML
(or `make pdf` for final pagination).  The building blocks:

- **Prose** — plain Markdown paragraphs.  Blank line between
  paragraphs.  `*italics*`, `**bold**`, `` `inline code` `` for
  keywords like `VptElt`.
- **Headings** — `##` is a SECTION/APPENDIX, `###` a subsection,
  `####` a command or sub-topic.  Follow the numbering style already
  in the file.
- **Terminal transcripts / .in file excerpts** — fenced code blocks:

  ````markdown
  ```
  MACOS>spot
  Enter element number: 11
  ```
  ````

- **Figures** — put a PNG in `src/media/`, reference it inline where
  the figure belongs, with the caption paragraph directly below:

  ```markdown
  ![](media/my_new_figure.png)

  **FIGURE 67** One-line caption text.
  ```

  Renumber ONLY if inserting between existing figures (grep for the
  neighbors first).  Optional sizing: `![](media/x.png){width="4in"}`.
- **Tables** — pipe tables are easiest to type:

  ```markdown
  | Location  | Frame A | Frame B |
  |-----------|---------|---------|
  | VptElt(1) | 2.5     | 0.0     |
  ```

  (The migrated grid tables `+---+` work too; no need to convert.)
- **Math** — LaTeX between dollar signs: `$W(x,y) = \sum a_j Z_j$` or
  `$$...$$` for display equations.  Renders natively in HTML and
  converts to Word equations in the docx.  This is the preferred
  replacement for the ~73 legacy equation images (see Open items).

### Adding a new command or feature

1. Write the description in the matching section file (e.g. a new
   `####` under Section 5 for a ray-trace command), following the
   pattern of the neighboring commands: short prose, then a transcript
   code block showing the dialog.
2. Add a bullet to **"New Features of Version X"** in
   `src/01_introduction.md`.
3. If the command was added to `macos_help.inc`, run `make appendix-c`
   so the command reference picks it up.  Never edit
   `src/93_appendix_c_commands.md` by hand — it is overwritten.
4. `make` and check the HTML.

### Adding a new section or appendix

Create a new numbered file — the number prefix sets its position
(e.g. `92_appendix_b_papers.md`).  Start it with a `##` heading.  No
registration needed; the Makefile globs `src/[0-9]*.md`.

### Updating the Appendix A examples

The `.in`/`.jou` files live in `examples/`.  When an example changes,
re-run it through the engine and paste the fresh transcript into
`91_appendix_a_examples.md`.  (Automating this — regenerating all
transcripts from `examples/*.jou` at build time, like Appendix C — is
a planned improvement.)

## Command Reference (cmdref/)

`cmdref/` holds the MACOS Command Reference: an orientation chapter
(what MACOS / SMACOS / mmacos / pymacos are and how they relate), a
syntax-conventions chapter per surface with a quickstart, and three
catalogs — Part I engine commands, Part II mmacos/pymacos functions
side by side, Part III the higher-level MATLAB layers.

Parts I–III are generated from the engine and binding sources
(`make cmdref-regen`), so the catalog tracks the code; hand-written
enrichment survives regeneration only inside the entry NOTES markers.
**See [`cmdref/README.md`](cmdref/README.md)** for the file map,
editing rules, and status.

## Gotchas

- **`reference.docx` controls Word/PDF styling** (pandoc-default
  styles, restyled: Arial headings, Times body).  To restyle the
  manual, edit the styles *in that file*; don't touch `src/`.  Do
  **not** substitute `macosMan4_0_styled.docx` — it is defective as a
  reference-doc and collapses every table to a single column.
- **PDF TOC**: a plain `soffice --convert-to pdf` leaves the TOC field
  empty.  `make pdf` uses `tools/docx2pdf.py` (LibreOffice UNO:
  refresh indexes, then export) — always build the PDF through make.
- **Do not re-run `tools/convert_from_docx.sh`** — it was the one-time
  docx→Markdown migration and would overwrite all `src/` edits.
- Legacy image sizes appear as `{width="0.1055…in"}` — fractional
  inches from the docx; harmless, leave unless intentionally resizing.

## History

The manual originated in FrameMaker (source lost), survived as
`docs/macosMan3.2.pdf`, was Acrobat-converted to docx, restyled, and
maintained 2026-04/05 as a chain of Python XML patches (`patches/`,
`scripts/`) producing `macosMan4_0_styled.docx` → `macosMan4.2beta.docx`.
On 2026-07-06 that endpoint was converted to the Markdown source in
`src/` (`tools/convert_from_docx.sh` = pandoc `docx+styles` +
`tools/style_map.lua` + `tools/split_sections.py`) and **Markdown
became canonical**.  The patch chain (`patches/`, `scripts/`) was
removed from HEAD on 2026-07-07 — it was unreferenced by the build
and carried static-analysis findings; recover from git history if
needed.  The docx snapshots are retained for provenance.

Two clean-ups happened at migration:

- **Writer fork folded in** — the parallel LibreOffice edit
  `macosMan4_01.docx` was diffed against the baseline;
  `RECONCILE_4_01.md` records what each side contributed (its
  "(historical)" version-notes sections now end
  `src/01_introduction.md`).
- **Composite figures rescued** — 34 FrameMaker line-art figures
  existed in the docx only as exploded fragments; each was rasterized
  from the vector-clean `macosMan3.2.pdf` into a single
  `src/media/fig_rescued_NN.png` (`tools/rescue_figures.py`; status
  and debris audit in `FIGURE_RESCUE_LOG.md`).

## Open items

1. **Equations**: ~73 remain as low-resolution images or Symbol-font
   fragments inherited from the PDF extraction; re-enter as LaTeX
   math (they were broken in every docx generation too — this is the
   one manual-labor item left).
2. **Figures 31, 56, 59** need hand-tightened crops (see
   `FIGURE_RESCUE_LOG.md`).
3. **New-feature content**: FreeForm details, SLSQP optimization, XPS,
   COMPOSE/broadband, giza graphics, design layer, etc. still need
   sections written.
4. **Appendix A automation**: regenerate example transcripts from
   `examples/*.jou` via the live engine at build time.
