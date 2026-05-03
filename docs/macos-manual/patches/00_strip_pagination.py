"""
Phase 2: Strip pagination artifacts inherited from FrameMaker -> PDF
-> DOCX conversion.

Removes:
  * 326 intermediate <w:sectPr> elements (one per original PDF page).
    Keeps the body-final <w:sectPr> that defines page geometry / margins
    / headers / footers for the document as a whole.
  * All 101 <w:br w:type="column"/> column-break runs.
  * Any <w:lastRenderedPageBreak/> markers (Word will regenerate).
  * Any <w:pageBreakBefore/> paragraph properties.
  * Empty <w:r> shells left behind by removed breaks.

After this patch, the document flows continuously and Word/LibreOffice
re-paginates from the single body-final sectPr.  Phase 3
(05_section_pagebreaks.py) then re-inserts page breaks before each
Heading2 banner.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import wtag


def transform(root):
    body = next(c for c in root if c.tag == wtag("body"))
    body_children = list(body)
    body_final_sectpr = next((c for c in reversed(body_children)
                              if c.tag == wtag("sectPr")), None)

    stripped_sectpr = 0
    stripped_col    = 0
    stripped_lrpb   = 0
    stripped_pbb    = 0
    stripped_runs   = 0

    # 1) Remove intermediate sectPr at body level.
    for c in body_children:
        if c.tag == wtag("sectPr") and c is not body_final_sectpr:
            body.remove(c)
            stripped_sectpr += 1

    # 2) Remove sectPr nested inside <w:pPr> (paragraph-level section
    #    breaks; rare but possible).
    for ppr in list(root.iter(wtag("pPr"))):
        for sect in list(ppr.findall(wtag("sectPr"))):
            ppr.remove(sect)
            stripped_sectpr += 1

    # 3) Remove all column breaks.  A <w:br w:type="column"/> sits in a
    #    run; remove the break, and if the run becomes content-empty,
    #    remove the run too.
    for r in list(root.iter(wtag("r"))):
        for br in list(r.findall(wtag("br"))):
            if br.get(wtag("type")) == "column":
                r.remove(br)
                stripped_col += 1
        # Was the run reduced to <w:rPr/> only?  Drop it.
        non_rpr_children = [c for c in r if c.tag != wtag("rPr")]
        if not non_rpr_children:
            # Find r's parent and remove it.
            for par in root.iter():
                if r in list(par):
                    par.remove(r)
                    stripped_runs += 1
                    break

    # 4) Strip <w:lastRenderedPageBreak/> (Word regenerates).
    for lrpb in list(root.iter(wtag("lastRenderedPageBreak"))):
        for par in root.iter():
            if lrpb in list(par):
                par.remove(lrpb)
                stripped_lrpb += 1
                break

    # 5) Strip <w:pageBreakBefore/> from paragraph properties.
    for ppr in root.iter(wtag("pPr")):
        for pbb in list(ppr.findall(wtag("pageBreakBefore"))):
            ppr.remove(pbb)
            stripped_pbb += 1

    return {
        "intermediate sectPr stripped": stripped_sectpr,
        "column breaks stripped":       stripped_col,
        "empty runs removed":           stripped_runs,
        "lastRenderedPageBreak removed": stripped_lrpb,
        "pageBreakBefore removed":      stripped_pbb,
    }


if __name__ == "__main__":
    apply(__file__, transform)
