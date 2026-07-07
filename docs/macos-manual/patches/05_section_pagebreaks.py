"""
Phase 3: re-insert page breaks before each major section.
Sets the <w:pageBreakBefore/> paragraph property on every Heading2 and
on the document Title -- so SECTION 1..9 and APPENDIX A..C each start
on a new page when Word/LibreOffice paginates.
Word's preferred mechanism is a paragraph property rather than an
explicit <w:br type="page"/> run, because the property survives cut/
paste and copy-of-style operations more cleanly.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from patch_base import apply
from docx_helpers import (iter_paragraphs, get_style, wtag)
import defusedxml.ElementTree as ET

PAGE_BREAK_STYLES = {"Heading2"}    # one per major SECTION / APPENDIX

def ensure_page_break_before(p):
    """Insert <w:pageBreakBefore/> in p's <w:pPr>, creating pPr if needed."""
    pPr = p.find(wtag("pPr"))
    if pPr is None:
        pPr = ET.Element(wtag("pPr"))
        p.insert(0, pPr)
    # Don't double-insert.
    if pPr.find(wtag("pageBreakBefore")) is not None:
        return False
    pbb = ET.Element(wtag("pageBreakBefore"))
    # Spec says pageBreakBefore comes after pStyle, before other
    # properties.  Insert at index 1 if pStyle is at 0; else at 0.
    pStyle = pPr.find(wtag("pStyle"))
    if pStyle is not None:
        children = list(pPr)
        idx = children.index(pStyle) + 1
        pPr.insert(idx, pbb)
    else:
        pPr.insert(0, pbb)
    return True

def transform(root):
    added = 0
    skipped_first = False
    for p in iter_paragraphs(root):
        if get_style(p) not in PAGE_BREAK_STYLES:
            continue
        if not skipped_first:
            skipped_first = True
            continue
        if ensure_page_break_before(p):
            added += 1
    return {"page-break-before set on Heading2": added}

if __name__ == "__main__":
    apply(__file__, transform)