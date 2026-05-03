"""
Patch I: append the SegMirMaker tool note at the end of §4.4.11
'Hex and Pie Segmented Mirrors'.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_LIST, P_H4
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, make_para)


HEAD = "Generating segmented prescriptions: SegMirMaker"

INTRO = (
    "Hand-editing a .in file to define every segment of a 7- or "
    "19-segment mirror is tedious.  MACOS v4.00 ships with a separate "
    "tool, SegMirMaker, in MACOS_resources/segmirmaker/, that "
    "generates a .presc (suitable for inclusion in a .in file) and an "
    "Hx.m edge-sensor measurement matrix from a parent prescription.")

ACCEPT_INTRO = "SegMirMaker accepts:"
ACCEPT_BULLETS = [
    "A MACOS .in parent prescription.  The parent element can be any "
    "conic or FreeForm surface.",
    "The number of rings (nRing = 1 for 7-segment, 2 for 19-segment, "
    "etc.).",
    "The axis orientation psi, segment size or aperture diameter, "
    "inter-segment gap, and segment standoff.",
    "3-DOF (piston/tip/tilt) or 6-DOF segments.",
]

PROPAGATE = (
    "For FreeForm parents, SegMirMaker copies the parent's lFF, "
    "FFCoef, FF coordinate frame, grid data, and grid frame verbatim "
    "into every segment, so all segments share the parent's underlying "
    "FreeForm shape.  Each segment's Mon slot is left empty "
    "(MonZernCoef = 0) but with lMon = L2 (half the segment width); "
    "the user can later add per-segment figure error via the Mon "
    "Zernike coefficients without re-running SegMirMaker.")

LINEAGE = (
    "SegMirMaker is the modernized successor to the SMPGe ('Segmented "
    "Mirror Prescription Generator') VAX Fortran tool from 1992; the "
    "original is preserved in "
    "MACOS_resources/segmirmaker/Archive/SMPGe.for.  See "
    "MACOS_resources/segmirmaker/README.md for full usage.")


def transform(root):
    h = None
    for p in iter_paragraphs(root):
        if get_style(p) == "Heading4" and \
           "Hex and Pie Segmented Mirrors" in get_text(p):
            h = p
            break
    if h is None:
        raise RuntimeError("Couldn't find §4.4.11 Hex and Pie heading")

    parent = find_para_parent(root, h)
    children = list(parent)
    idx = children.index(h) + 1
    # Walk to next H4 (next sub-section) or H3.
    while idx < len(children):
        c = children[idx]
        if c.tag == wtag("p") and get_style(c) in ("Heading3", "Heading4"):
            break
        idx += 1

    new_paras = [
        make_para(P_H4, HEAD),
        make_para(P_BODY, INTRO),
        make_para(P_BODY, ACCEPT_INTRO),
    ]
    new_paras.extend(make_para(P_LIST, b) for b in ACCEPT_BULLETS)
    new_paras.append(make_para(P_BODY, PROPAGATE))
    new_paras.append(make_para(P_BODY, LINEAGE))

    for offset, np in enumerate(new_paras):
        parent.insert(idx + offset, np)
    return {"new paragraphs": len(new_paras)}


if __name__ == "__main__":
    apply(__file__, transform)
