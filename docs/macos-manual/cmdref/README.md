# MACOS Command Reference — source

This directory is the source of the **MACOS Command Reference**, a
standalone document cataloging every way to drive the MACOS engine:
the interactive CLI, the SMACOS subroutine interface, and the mmacos
(MATLAB) and pymacos (Python) bindings.  It is written to orient new
users of the public repo; the MACOS Manual (`../src/`) covers
concepts and worked examples, this document covers commands.

## Building

From `docs/macos-manual/`:

```bash
make cmdref       # build/macosCmdRef.docx + build/macosCmdRef.html
make cmdref-pdf   # build/macosCmdRef.pdf
make cmdref-regen # regenerate the Part I-III catalogs from the code
```

The chapters are the numbered `.md` files, assembled in lexical
order:

| File | Content | Maintained |
|---|---|---|
| `00_orientation.md` | what MACOS/SMACOS/mmacos/pymacos are, how they relate, where they live | by hand |
| `01_cli_syntax.md` | CLI conventions (min-match, dialogs, journals) + quickstart | by hand |
| `02_smacos_syntax.md` | SMACOS call convention + init | by hand |
| `03_mmacos_syntax.md` | MATLAB conventions + quickstart | by hand |
| `04_pymacos_syntax.md` | Python conventions + quickstart | by hand |
| `10_engine_commands.md` | Part I: ~120 engine commands by HELP category | **generated** |
| `20_bindings.md` | Part II: mmacos/pymacos functions side by side | **generated** |
| `30_higher_level.md` | Part III: `macos.design`, `macos.channels` | **generated** |

## Editing rules

- The `00_`–`04_` chapters are ordinary Markdown — edit freely.
- The Part I–III catalogs are **generated** by `tools/gen_cmdref.py`
  from `macos_f90/macos_help.inc`, the mmacos `.m` help headers, the
  pymacos docstrings, and `macos_api_mod.F90`.  Running
  `make cmdref-regen` overwrites them — with one exception:
- **Hand-written prose survives regeneration only inside the NOTES
  markers** that every entry carries:

  ```markdown
  <!-- BEGIN NOTES cmd-PERturb -->
  Your expanded description, dialog transcript, related commands...
  <!-- END NOTES cmd-PERturb -->
  ```

  Write enrichment there (the `TODO` marker disappears once a NOTES
  block has content).  Anything outside the markers in a generated
  file will be lost on the next regen.
- New engine commands and binding functions appear in the catalogs
  automatically on regen; a new cross-language naming mismatch shows
  up as a false "*not available*" gap — fix it by extending the
  `SYNONYMS` map in `tools/gen_cmdref.py`.
- Improving a *description* is often better done at the source — the
  `.m` help header or Python docstring — so `help`/`help()` users
  benefit too; the catalog inherits it on regen.

## Status

- Phase A (2026-07-07): scaffolding, syntax chapters, generated
  catalogs — complete; document builds and is usable.
- Phase B (pending): expand Part I entries category by category with
  dialog transcripts and behavior notes.
- Phase C (pending): expand Part II entries cluster by cluster;
  the "*not available*" flags double as a binding-parity worklist.

Part III documents packages under active development; expect it to
change faster than Parts I–II.
