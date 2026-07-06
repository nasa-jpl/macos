#!/usr/bin/env python3
"""Split the pandoc-converted manual (one big .md) into per-section
source files under src/.  One-time migration helper, kept for
provenance; re-runnable.

Usage: python tools/split_sections.py full_mapped.md src/

Rules:
  - Split on '## ' headings (the docx sections; '#' is the title page).
  - The converted 'Table of Contents' section is dropped (the build
    regenerates a live TOC via pandoc --toc).
  - The legacy back-of-book index (stale v3.2 page numbers, embedded
    after SECTION 7 by the original Acrobat extraction) is dropped;
    a copy is kept in _dropped_legacy_index.md for reference.
  - Image paths './media/...' are rewritten to 'media/...'.
"""
import re
import sys
from pathlib import Path

# heading-prefix -> output file name
FILEMAP = [
    ("SECTION 1", "01_introduction.md"),
    ("SECTION 2", "02_technical_overview.md"),
    ("SECTION 3", "03_user_interface.md"),
    ("SECTION 4", "04_describing_optical_systems.md"),
    ("SECTION 5", "05_ray_trace_analysis.md"),
    ("SECTION 6", "06_diffraction_analysis.md"),
    ("SECTION 7", "07_beam_propagation_imaging.md"),
    ("SECTION 8", "08_differential_ray_tracing.md"),
    ("SECTION 9", "09_subroutine_macos.md"),
    ("REFERENCES", "90_references.md"),
    ("APPENDIX A", "91_appendix_a_examples.md"),
    ("APPENDIX B", "92_appendix_b.md"),
    ("APPENDIX C", "93_appendix_c_commands.md"),
]


def target_for(heading):
    for prefix, fname in FILEMAP:
        if heading.startswith(prefix):
            return fname
    if heading.startswith("Table of Contents"):
        return None  # dropped: regenerated at build time
    return "00_frontmatter.md"  # title/authorship/copyright chunks


def strip_legacy_index(lines):
    """Remove the embedded legacy index: from the '**Symbols**' marker
    (first letter-index block) through the end of the '**Z**' entries.
    Returns (kept_lines, dropped_lines)."""
    letter = re.compile(r"^\*\*(Symbols|[A-Z])\*\*$")
    start = None
    for i, ln in enumerate(lines):
        if ln.strip() == "**Symbols**":
            start = i
            break
    if start is None:
        return lines, []
    # walk forward: index ends at the last letter-block's entries,
    # i.e. just before the next '## ' heading
    end = len(lines)
    for j in range(start, len(lines)):
        if lines[j].startswith("## "):
            end = j
            break
    return lines[:start] + lines[end:], lines[start:end]


def main(src_md, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    text = Path(src_md).read_text()
    text = text.replace("](./media/", "](media/")
    lines = text.splitlines()

    lines, dropped = strip_legacy_index(lines)
    if dropped:
        (out_dir / "_dropped_legacy_index.md").write_text(
            "\n".join(dropped) + "\n")

    chunks = {}          # fname -> list of lines
    order = []
    current = "00_frontmatter.md"
    for ln in lines:
        if ln.startswith("## "):
            current = target_for(ln[3:].strip())
        if current is None:
            continue
        chunks.setdefault(current, [])
        if current not in order:
            order.append(current)
        chunks[current].append(ln)

    for fname, body in chunks.items():
        # trim leading/trailing blank lines
        while body and not body[0].strip():
            body.pop(0)
        while body and not body[-1].strip():
            body.pop()
        (out_dir / fname).write_text("\n".join(body) + "\n")
        print(f"  {fname:42s} {len(body):6d} lines")
    if dropped:
        print(f"  _dropped_legacy_index.md (excised)         "
              f"{len(dropped):6d} lines")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
