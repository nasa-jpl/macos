"""
Phase: remove horizontal-line artifacts left by the FrameMaker -> PDF
-> DOCX conversion.

The original FrameMaker template put a thin horizontal rule across the
top and bottom of every page (header/footer ornaments).  Acrobat
preserved each rule as a DrawingML <w:drawing> element with a
<wp:extent cx="6172200" cy="..."/> at the page's full text width.
There are 1,368 of them in the source document, all with the same
cx and varying cy (1440 .. 25560 EMU = roughly 0.002 .. 0.028 inch).

Strategy
--------
For every <w:drawing> in the document:

  1. Parse its <wp:extent cx="..." cy="..."/>.
  2. If the shape is full-text-width (cx == 6172200 EMU) and thin
     (cy < CY_MAX_EMU), it's a header/footer rule -- remove the whole
     <w:drawing>.
  3. If removing the drawing leaves its parent <w:r> empty of content,
     drop the run too.
  4. If that leaves the paragraph with no content runs and no other
     meaningful children, drop the paragraph as well.

Figure-internal rules (axis lines, plot dividers, leader lines) use
different cx values and are *not* removed.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import wtag

# DrawingML wp namespace.
WP = "http://schemas.openxmlformats.org/drawingml/2006/wordprocessingDrawing"

# Header/footer rule signature.  All 1,368 rules in the source share
# this canonical width; figure-internal lines use other widths.
RULE_CX_EMU = 6172200
CY_MAX_EMU  = 100000        # ~0.11 in -- accept any thickness up to this
                            # so thicker rules (if any) are still caught


def is_header_footer_rule(drawing) -> bool:
    """True if this <w:drawing> wraps a full-text-width thin rule."""
    extent = drawing.find(f".//{{{WP}}}extent")
    if extent is None:
        return False
    try:
        cx = int(extent.get("cx", "0"))
        cy = int(extent.get("cy", "0"))
    except ValueError:
        return False
    return cx == RULE_CX_EMU and 0 <= cy < CY_MAX_EMU


def find_parent_map(root):
    pm = {}
    for parent in root.iter():
        for child in parent:
            pm[id(child)] = parent
    return pm


def run_is_empty(r) -> bool:
    """True if the run has no <w:t>, <w:drawing>, <w:br>, <w:tab> etc."""
    for child in r:
        tag = child.tag
        if tag in (wtag("rPr"),):
            continue
        # Anything else counts as content.
        return False
    return True


def paragraph_is_empty(p) -> bool:
    """True if the paragraph has no content runs and no images / tables."""
    for child in p:
        tag = child.tag
        if tag == wtag("pPr"):
            continue
        if tag == wtag("r"):
            if not run_is_empty(child):
                return False
            continue
        # bookmarks, comments, etc. -- ignore as content for our purposes
        if tag in (wtag("bookmarkStart"), wtag("bookmarkEnd"),
                   wtag("commentRangeStart"), wtag("commentRangeEnd"),
                   wtag("proofErr")):
            continue
        # Real content (sectPr, smartTag, ...): keep
        return False
    return True


def transform(root):
    pm = find_parent_map(root)

    # Collect drawings to remove.
    to_remove_drawings = []
    for d in root.iter(wtag("drawing")):
        if is_header_footer_rule(d):
            to_remove_drawings.append(d)

    # Remove each drawing, then conditionally its run, then conditionally
    # its paragraph.
    drawings_removed = 0
    runs_removed     = 0
    paras_removed    = 0
    for d in to_remove_drawings:
        run    = pm.get(id(d))             # parent <w:r>
        para   = pm.get(id(run)) if run is not None else None
        para_parent = pm.get(id(para)) if para is not None else None

        # 1. Remove the drawing from its run.
        if run is not None:
            try:
                run.remove(d)
                drawings_removed += 1
            except ValueError:
                continue
        else:
            continue

        # 2. If the run is now empty, remove it from the paragraph.
        if run_is_empty(run) and para is not None:
            try:
                para.remove(run)
                runs_removed += 1
            except ValueError:
                pass

            # 3. If the paragraph is now empty, remove it.
            if paragraph_is_empty(para) and para_parent is not None:
                try:
                    para_parent.remove(para)
                    paras_removed += 1
                except ValueError:
                    pass

    return {
        "horizontal-rule drawings removed": drawings_removed,
        "empty runs removed":                runs_removed,
        "empty rule-only paragraphs removed": paras_removed,
        "rule signature":                    f"cx={RULE_CX_EMU} cy<{CY_MAX_EMU}",
    }


if __name__ == "__main__":
    apply(__file__, transform)
