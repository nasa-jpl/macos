# Reconciliation: macosMan4_01.docx (Writer fork) vs macosMan4.2beta.docx

Date: 2026-07-06.  The Markdown source in `src/` was migrated from
**macosMan4.2beta.docx** (the docx patch-chain output).  Dave's
LibreOffice Writer fork **macosMan4_01.docx** was converted with the
same pipeline and diffed against the baseline (~2,000 changed lines).
Findings:

## Baseline (4.2beta) content missing from 4_01 — kept

- `Authorship` block on the title page.
- `VALIDATE` command section (§3).
- The rewritten §9 "Sharing Data via Fortran Modules" cluster
  (SMACOS modules, initializing the runtime, minimal driver,
  save/restore, legacy COMMON support) — 4_01 still had the old
  "Common Blocks and Include Files" text.
- General text cleanup: 4_01 carries roughly 2x the Acrobat
  line-break hyphenation debris ("sur-face", "coni-coid", ...).

## 4_01 content missing from baseline — ported into src/

- "New Features of Version 3.2" and "Changes since MACOS 2.8",
  which Dave had retitled "(historical)" in the Writer fork.  Both
  now live at the end of `src/01_introduction.md`, dehyphenated.

## Not carried

- 4_01's remaining wording deltas are hyphenation/whitespace noise
  or older phrasings of text the patch chain later improved; no
  other substantive content was found on the 4_01 side.

The full machine diff used for this audit can be regenerated with
`tools/convert_from_docx.sh` applied to each docx (into scratch
directories) and `diff` on the outputs.
