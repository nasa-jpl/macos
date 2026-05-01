"""
Patch B1: replace the §3.2 'Help' body with the v4.00 description of
the rebuilt 15-category HELP output.  The §3.2 heading itself is left
intact; everything between '### Help' and the next ### heading is
swapped for the new content.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_LIST, P_H3
from docx_helpers import (find_para_by_text, find_para_parent, get_style,
                          insert_after_para, make_para, wtag)

HEAD = (
    "Type HELP (or ?) at the MACOS> prompt to see the list of all "
    "commands.  The output is paginated -- hit return to advance to "
    "the next screen.  Commands are grouped into 15 categories:")

CATEGORIES = [
    "Session & files: QUit, END, HELP, REset, STatus, SUMmarize, ELTS, "
    "EXEcute, JOUrnal, MREset, !, PWD, CD, RX/LS/LL, VI/EMAcs",
    "Prescription I/O: NEW, OLD, FID, SAVe, EXPort, SHOw, MODify",
    "Source & wavelength: CHIefray, WLENS, SWL, MULtispec, "
    "NFIlt/RFIlt/SFIlt, ATMosphere, SETC, SAOpt",
    "Ray tracing: RAY, SEGraytrace, PRAy/RRAy/TPR, MAP, COOrd, ACOor, "
    "EFL, BWK, FEXit, FSR, CENter, STOp, CENTRoid, FFP, PFP, FDP, "
    "SPCenter",
    "Behavior toggles: LNEg/NOLNeg, POLarized/NOPol, OBS, "
    "REGrid/NOREGrid, ORS, SRS, NONe",
    "Surface data: SINt, UDSinit",
    "Perturbation: PERturb, GPERturb, PREad",
    "Linear model: BUild, DMBuild, PARtials, LPErturb, LPRead, LREset, "
    "LSPot, LOPd, LPIxilate",
    "Diffraction: PROpagate, BEAm, VECtor, SCAlar, COMpose, ADD, DAdd, "
    "CADD, NOIse, SEEd, STRetch",
    "Outputs: OPD, SPOt, ZABerr, ZCOef, LOS, METcalc, INTensity, PIXel, "
    "AMPlitude, PHAse/PH, REAl, LOG, MTF/CMTF, IMG",
    "Plot style: GRAy, WIRe, SLIce, COLumn, CONtour, IMGmode, CIR, GIR, "
    "PGP",
    "File output: TEXt, BINary, FITs/WFIts, MAT, GETMatlabmatrix",
    "Window & post-processing: PLOcate, NOPLOC, PWIn, SZCo, GBS, "
    "BLUr/GBLur, GAIn, ODRaw, PGD, ROW",
    "System optimization: AVAR, MVAR, DVAR, VARS, AFOV, DFOV, FOVS, "
    "CALib, SFOV",
    "Misc / debug: SRAy, SRT, SRTrace, DOPD, DOPL, ZRM, "
    "JWST_v3d / Vis3d",
]

CASING = (
    "Casing convention.  Commands are shown as <UPPER prefix><lower "
    "tail>, e.g. RAYtrace, SUMmarize, PERturb.  The uppercase part is "
    "the minimum-match abbreviation tested by the dispatcher -- typing "
    "just RAY, SUM, or PER (case-insensitive) is sufficient.")

TAGS = [
    "Prerequisite tags. Lines in the HELP output may end in:",
    "[Rx] -- needs a loaded prescription (run OLD or NEW first)",
    "[BLD] -- needs a built linear model (run BUild or DMBuild)",
    "[DIFF] -- needs a propagated wavefront (run PROpagate first)",
]

CLOSE = ("Commands without a tag work in any state.  A complete "
         "reference of all commands appears in Appendix C.")


def transform(root):
    # Anchor: the '### Help' heading is paragraph styled Heading3 with
    # text exactly 'Help'.
    h = None
    from docx_helpers import iter_paragraphs, get_text
    for p in iter_paragraphs(root):
        if get_style(p) == P_H3 and get_text(p).strip() == "Help":
            h = p
            break
    if h is None:
        raise RuntimeError("Couldn't find §3.2 Help heading")

    # Sweep forward, deleting paragraphs until we hit the next NON-EMPTY
    # H3 (which is 'MACOS Model Size Specification' or 'Entering
    # Commands').  An empty Heading3 immediately after 'Help' is a
    # stylistic placeholder in the source and would short-circuit the
    # delete; skip it.
    parent = find_para_parent(root, h)
    children = list(parent)
    idx = children.index(h)
    deleted = 0
    j = idx + 1
    while j < len(children):
        c = children[j]
        if c.tag == wtag("p") and get_style(c) == P_H3 \
                and get_text(c).strip():
            break
        parent.remove(c)
        deleted += 1
        children = list(parent)

    new_paras = [make_para(P_BODY, HEAD)]
    new_paras.extend(make_para(P_LIST, c) for c in CATEGORIES)
    new_paras.append(make_para(P_BODY, CASING))
    new_paras.extend(make_para(P_LIST, t) if i > 0 else make_para(P_BODY, t)
                     for i, t in enumerate(TAGS))
    new_paras.append(make_para(P_BODY, CLOSE))
    insert_after_para(root, h, new_paras)
    return {"old paragraphs deleted": deleted,
            "new paragraphs": len(new_paras)}


if __name__ == "__main__":
    apply(__file__, transform)
