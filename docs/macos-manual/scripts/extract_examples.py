"""
extract_examples.py — Extract Appendix A example files (.in / .jou /
.test) from the manual into ~/dev/macos/docs/macos-manual/examples/.

Each example in the manual is introduced by a Heading4 paragraph
whose text matches '*.in', '*.jou', or '*.test', followed by a
sequence of CodeBlock / TOC* / Normal paragraphs containing the file
contents.  We harvest everything from the heading up to the next
Heading2/3/4, drop page-footer artifacts ('Modeling and Analysis...
NN'), and write the concatenated text to disk under examples/.

Usage:
    python scripts/unpack.py macosMan4.1.docx
    python scripts/extract_examples.py
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from docx_helpers import (load_doc, iter_paragraphs, get_style, get_text)

PROJECT_ROOT = Path(__file__).parent.parent
WORKING_DIR  = PROJECT_ROOT / "working"
OUTPUT_DIR   = PROJECT_ROOT / "examples"

EXAMPLE_RE = re.compile(r"\.(?:jou|in|test)$", re.IGNORECASE)
HEADING_STYLES = {"Heading2", "Heading3", "Heading4"}

# Page-footer / running-header patterns to drop.
FOOTER_RES = [
    re.compile(r"^Modeling and Analysis for Controlled Optical Systems\s*\d+\s*$"),
    re.compile(r"^\s*\d{1,4}\s*$"),                 # bare page number
    re.compile(r"^\s*[A-Z][\w\.]+\s+\d{1,4}\s*$"),  # 'Cassegrain.in 167'
]


def is_footer_line(text: str) -> bool:
    text = text.strip()
    if not text:
        return False
    return any(rgx.match(text) for rgx in FOOTER_RES)


def extract_one(paragraphs, start_idx):
    """
    Walk paragraphs[start_idx+1:] until the next Heading2/3/4.
    Return list of text lines for the example body.
    """
    lines = []
    j = start_idx + 1
    while j < len(paragraphs):
        p = paragraphs[j]
        s = get_style(p)
        if s in HEADING_STYLES and get_text(p).strip():
            break
        text = get_text(p)
        if text.strip() and not is_footer_line(text):
            # Each paragraph becomes one logical "line" in the source
            # file.  Multi-line paragraphs (very rare in CodeBlock --
            # codeblocks usually emit one paragraph per source line)
            # already contain literal '\n' from <w:br/> within the run;
            # emit them verbatim.
            lines.append(text.rstrip())
        j += 1
    return lines, j


def main():
    if not WORKING_DIR.exists():
        sys.exit("working/ not found -- run: python scripts/unpack.py macosMan4.1.docx")
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    _, root = load_doc(WORKING_DIR)
    paras = list(iter_paragraphs(root))

    examples = []
    for i, p in enumerate(paras):
        s = get_style(p)
        t = get_text(p).strip()
        if s.startswith("Heading") and EXAMPLE_RE.search(t):
            examples.append((i, t))

    print(f"Found {len(examples)} example headings.\n")

    written = 0
    skipped = 0
    for idx, name in examples:
        lines, end = extract_one(paras, idx)
        if not lines:
            print(f"  skip  {name:30s}  (no content captured)")
            skipped += 1
            continue
        out_path = OUTPUT_DIR / name
        out_path.write_text("\n".join(lines) + "\n")
        size = out_path.stat().st_size
        print(f"  wrote {name:30s}  {len(lines):4d} lines, {size:6d} bytes "
              f"(paras {idx+1}..{end-1})")
        written += 1

    print(f"\nWrote {written} files, skipped {skipped}, into {OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
