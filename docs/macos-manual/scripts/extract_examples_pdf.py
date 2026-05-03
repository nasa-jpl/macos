"""
extract_examples_pdf.py — Re-extract Appendix A example files
directly from the original PDF (using pdftotext -layout).

The .docx round-trip dropped values from many Key= entries in
AOExample.in, SegDemo.in, Luneberg.in, Cassegrain.in.  The PDF text
layer preserves them, so re-extracting from there gives clean files.

Usage:
    pdftotext -layout ~/dev/macos/docs/macosMan3.2.pdf /tmp/man.txt
    python scripts/extract_examples_pdf.py /tmp/man.txt

Writes one file per Appendix-A example into examples/ (overwriting).
"""
import re
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
EXAMPLES_DIR = PROJECT_ROOT / "examples"

# Appendix-A section headers, in document order.  The list is the
# authoritative file roster -- if a heading isn't here the PDF text
# pass skips it.
EXAMPLES = [
    "A.1.1 Cassegrain.jou",
    "A.1.2 Cassegrain.in",
    "A.1.3 CassWithExitPupil.jou",
    "A.1.4 CassWithExitPupil.in",
    "A.2.1 GratingExample.jou",
    "A.2.2 GratingExample.in",
    "A.3.1 CornerCube.jou",
    "A.3.2 CornerCube.in",
    "A.3.3 Luneberg.jou",
    "A.3.4 Luneberg.in",
    "A.4.1 AOExample.jou",
    "A.4.2 AOExample.in",
    "A.4.3 AOmirror.test",
    "A.5.1 SegDemo.jou",
    "A.5.2 SegDemo.in",
    "A.6.1 coroExample.jou",
    "A.6.2 coroExample.in",
    "A.7.1 ImageDemo.jou",
]

# Page-running-header / footer patterns to drop.
DROP_RES = [
    # 'Modeling and Analysis... 191' (page-number trailing)
    re.compile(r"^\s*Modeling and Analysis for Controlled Optical Systems\s*\d*\s*$"),
    # '190    Modeling and Analysis...' (page-number leading)
    re.compile(r"^\s*\d+\s+Modeling and Analysis for Controlled Optical Systems\s*$"),
    re.compile(r"^\s*\d+\s*$"),                          # bare page number
    re.compile(r"^\s*[A-Z][\w\s\.\-]+\s+\d{2,4}\s*$"),   # 'Cassegrain.in 167'
]


def is_dropline(s: str) -> bool:
    return any(rx.match(s) for rx in DROP_RES)


def filename_from_heading(h: str) -> str:
    """'A.4.2 AOExample.in' -> 'AOExample.in'"""
    parts = h.split(None, 1)
    return parts[1] if len(parts) > 1 else h


def find_section_starts(lines, headings):
    """
    Walk lines once; for each heading record the line number where it
    first appears as a section heading (i.e., not inside the TOC).
    The TOC entries have trailing page numbers; section headings do not.
    """
    pos = {}
    for i, line in enumerate(lines):
        for h in headings:
            if h in pos:
                continue
            stripped = line.strip()
            # Match the heading exactly OR with leading whitespace --
            # but NOT with a trailing page number (TOC entries).
            if stripped == h:
                pos[h] = i
    return pos


def extract_block(lines, start, end):
    """Lines from after the heading up to the next section, dropping
    page-runner lines and trimming common leading whitespace."""
    body = lines[start + 1:end]
    # Drop page-runner / footer lines.
    body = [ln for ln in body if not is_dropline(ln)]
    # Drop leading blank lines.
    while body and not body[0].strip():
        body.pop(0)
    while body and not body[-1].strip():
        body.pop()
    if not body:
        return []
    # Determine common leading whitespace among non-blank lines and
    # strip it.
    nonblank = [ln for ln in body if ln.strip()]
    if nonblank:
        common = min(len(ln) - len(ln.lstrip()) for ln in nonblank)
        if common > 0:
            body = [ln[common:] if len(ln) >= common else ln
                    for ln in body]
    # Collapse runs of consecutive blank lines to a single blank.
    collapsed = []
    prev_blank = False
    for ln in body:
        ln = ln.rstrip()
        if not ln.strip():
            if not prev_blank:
                collapsed.append("")
            prev_blank = True
        else:
            collapsed.append(ln)
            prev_blank = False
    return collapsed


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: extract_examples_pdf.py <pdftotext -layout output>")
    text = Path(sys.argv[1]).read_text()
    lines = text.splitlines()

    pos = find_section_starts(lines, EXAMPLES)
    missing = [h for h in EXAMPLES if h not in pos]
    if missing:
        print(f"WARN: not found in PDF text:")
        for h in missing:
            print(f"   {h}")

    # Compute end-of-section for each heading: next heading's line, or
    # next 'A.x' chapter delimiter (e.g. 'APPENDIX B'), or EOF.
    sorted_pos = sorted(pos.items(), key=lambda kv: kv[1])
    enders = [p for h, p in sorted_pos]
    bounds = {}
    for k, (h, start) in enumerate(sorted_pos):
        end = enders[k + 1] if k + 1 < len(enders) else len(lines)
        bounds[h] = (start, end)

    EXAMPLES_DIR.mkdir(parents=True, exist_ok=True)
    written = 0
    for h in EXAMPLES:
        if h not in bounds:
            continue
        start, end = bounds[h]
        body = extract_block(lines, start, end)
        if not body:
            print(f"  skip   {h}  (empty)")
            continue
        fname = filename_from_heading(h)
        out_path = EXAMPLES_DIR / fname
        out_path.write_text("\n".join(body) + "\n")
        size = out_path.stat().st_size
        print(f"  wrote  {fname:30s}  {len(body):4d} lines, {size:6d} bytes")
        written += 1

    print(f"\nWrote {written}/{len(EXAMPLES)} into {EXAMPLES_DIR}/")


if __name__ == "__main__":
    main()
