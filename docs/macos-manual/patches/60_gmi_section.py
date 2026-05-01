"""
Patch J1: append new sub-section 'GMI: MATLAB MEX interface' inside §9
SMACOS, just before §10 (Appendix) heading.  If no §10 exists, append
at the very end of the document body.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_LIST, P_H3
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, make_para)


HEAD = "GMI: MATLAB MEX interface (FreeForm updates)"

INTRO = (
    "GMI (MACOS_resources/GMI/) is a thin MATLAB MEX wrapper around "
    "SMACOS, packaging perturbations and pulling back OPD / PIX / "
    "SPOT / complex-EF arrays into MATLAB.  V4.00 changes:")

BULLETS = [
    "NRHS = 14 (was 13).  Old call_GMI.m callers passing 13 "
    "arguments will fail with 'GMI requires 14 input arguments'.  "
    "The new 14th argument is pmonzern, a MonZernCoef perturbation "
    "vector parallel to pzern for FreeForm surfaces.",
    "param.monzernSrf -- column vector of FreeForm element indices "
    "to perturb (analogous to param.zernSrf).  Encoded in pflg with "
    "a 9999 sentinel when absent.",
    "param.pmonzern -- per-mode coefficient vector, length "
    "length(monzernSrf) * nmonzern, where nmonzern defaults to "
    "param.mzern.  Read inside call_GMI.m and forwarded as the 14th "
    "MEX argument; user code does not call GMI() directly.",
    "GMI.inc compile-time sizes: numseg = 7, mzern = 15 (up from "
    "6 / 12) to support 7-segment FreeForm test cases with 15 modes "
    "per segment.  When a runtime param.mzern exceeds the "
    "compile-time mzern, the MEX validation rejects with 'pzern must "
    "be scalar, empty or a mpzern x 1 matrix'.  Bump mzern and force "
    "a rebuild (rm GMI.o GMIG.o; source makegmi.sh -- make does not "
    "detect .inc changes).",
    "Nominal save/restore expanded to FreeForm frames.  "
    "ObtainNominalSettings / SetToNominalSettings now snapshot and "
    "restore pFF/xFF/yFF/zFF and (when nGridMat > 0) "
    "pData/xData/yData/zData alongside the existing pMon / psiElt / "
    "vptElt / TElt entries.  This fixes a v3.x bug where per-call "
    "prb perturbations on FreeForm elements leaked the FreeForm "
    "coordinate frames forward, so the 'nominal-equivalent' OPD "
    "drifted across iterations of a sensitivity loop.",
]

CLOSE = ("See MACOS_resources/GMI/README.md and "
         "MACOS_resources/GMI/CLAUDE.md for full architecture details.")


def transform(root):
    # Anchor: find heading containing 'SECTION 9' (Heading2).  Then
    # insert just before the next Heading2 or at the end of its parent.
    sec9 = None
    for p in iter_paragraphs(root):
        if get_style(p) == "Heading2" and "SECTION 9" in get_text(p):
            sec9 = p
            break
    if sec9 is None:
        raise RuntimeError("Couldn't find SECTION 9 heading")

    parent = find_para_parent(root, sec9)
    children = list(parent)
    idx = children.index(sec9) + 1
    insert_at = len(children)
    while idx < len(children):
        c = children[idx]
        if c.tag == wtag("p") and get_style(c) == "Heading2":
            insert_at = idx
            break
        idx += 1

    new_paras = [
        make_para(P_H3, HEAD),
        make_para(P_BODY, INTRO),
    ]
    new_paras.extend(make_para(P_LIST, b) for b in BULLETS)
    new_paras.append(make_para(P_BODY, CLOSE))

    for offset, np in enumerate(new_paras):
        parent.insert(insert_at + offset, np)
    return {"new paragraphs": len(new_paras)}


if __name__ == "__main__":
    apply(__file__, transform)
