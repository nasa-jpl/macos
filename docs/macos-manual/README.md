# MACOS Manual — Development Project

## Structure

```
macos-manual/
├── macosMan3_2_styled.docx     ← working document (open in LibreOffice)
├── scripts/
│   ├── unpack.py               ← extract docx → working/
│   ├── pack.py                 ← repack working/ → docx
│   ├── docx_helpers.py         ← shared XML utilities (import this)
│   ├── inspect.py              ← explore document without changing it
│   └── transform_template.py  ← copy this to start a new transform
└── working/                    ← unpacked XML (git-ignore this)
```

## Workflow

Every edit cycle follows the same three steps:

```bash
python scripts/unpack.py          # 1. unpack docx into working/
# ... edit working/word/document.xml, or run a transform script ...
python scripts/pack.py            # 2. repack into macosMan3_2_styled.docx
# open in LibreOffice to verify   # 3. check result visually
```

## Exploring the document

```bash
python scripts/inspect.py                    # paragraph style counts
python scripts/inspect.py headings           # full heading tree
python scripts/inspect.py sample CodeBlock   # sample paragraphs of a style
python scripts/inspect.py search "MACOS>"   # find paragraphs by content
```

## Writing a new transform

```bash
cp scripts/transform_template.py scripts/my_change.py
# edit my_change.py — fill in the transform() function
python scripts/unpack.py
python scripts/my_change.py
python scripts/pack.py
```

## Current document style inventory

| Style          | Count | Meaning                                    |
|----------------|-------|--------------------------------------------|
| TableParagraph | 4604  | Content inside tables                      |
| BodyText       | 2954  | Main prose paragraphs                      |
| CodeBlock      | 2917  | MACOS terminal output & command examples   |
| Normal         | ~1354 | Mixed: diagram fragments, math, misc       |
| ListParagraph  | 403   | Bulleted/numbered lists                    |
| Heading4       | 124   | Command-level headings (e.g. INTENSITY)    |
| FigureCaption  | 83    | FIGURE x and TABLE x labels                |
| Heading3       | 59    | Sub-section headings                       |
| Heading2       | 12    | Section headings (SECTION 1–9, APPENDIX)   |
| Heading1       | 1     | Document title                             |

## Known issues requiring manual attention

- **73 Symbol-font paragraphs** — math equations rendered in Symbol charset.
  These need to be rewritten using LibreOffice's equation editor (Insert → Object → Formula).
  Find them with:  `python scripts/inspect.py sample Normal`
  (look for entries with `font='Symbol'`)

- **7 large-font diagram fragments** — label text orphaned from FrameMaker figures.
  They sit near images; decide whether to delete or reformat as captions.

- **~1274 remaining Normal paragraphs** — includes diagram coordinate labels,
  section-TOC remnants, and inline notes. Worth a visual scroll in LibreOffice
  using the Navigator (F5) to catch anything structural.
