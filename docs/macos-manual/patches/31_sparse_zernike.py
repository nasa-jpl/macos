"""
Patch F1: append the sparse Zernike-input note to the end of §4.5.4
'Zernike Surfaces'.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_CODE, P_H4
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, make_para)

HEAD = "Sparse Zernike input"

INTRO = (
    "In addition to the dense MonCoef form (120 coefficients in "
    "Cartesian monomial order), MACOS reads Zernike coefficients in a "
    "sparse 'modes + coefs' form:")

EXAMPLE = (
    "    MonZernType= BornWolf\n"
    "    nMonZernCoef= 4\n"
    "    MonZernModes= 4  9  16  25\n"
    "    MonZernCoef= -4.819E+00  4.001E+00  1.197E+00  2.269E-03")

TYPES = (
    "MonZernType selects the Zernike ordering and normalization.  "
    "Supported types: ANSI, BornWolf, Fringe, NormANSI, NormBornWolf, "
    "NormFringe, NormNoll, Noll, NormAnnularNoll, NormHex, ExtFringe.  "
    "The Norm* variants use RMS normalization; the unprefixed variants "
    "are unnormalized.")

CONVERSION = (
    "Conversion to the internal monomial form (MonCoef) happens "
    "automatically at trace time.  The same scheme applies to FF "
    "Zernikes: FFZernType, nFFZernCoef, FFZernModes, FFZernCoef.  "
    "Modes are 1-indexed.")


def transform(root):
    h = None
    for p in iter_paragraphs(root):
        if get_style(p) == "Heading4" and \
           "Zernike Surfaces" in get_text(p):
            h = p
            break
    if h is None:
        raise RuntimeError("Couldn't find §4.5.4 Zernike Surfaces heading")

    parent = find_para_parent(root, h)
    children = list(parent)
    idx = children.index(h) + 1
    # Insert before next H4 (Monomial Surfaces) or H3.
    while idx < len(children):
        c = children[idx]
        if c.tag == wtag("p") and get_style(c) in ("Heading3", "Heading4"):
            break
        idx += 1

    new_paras = [
        make_para(P_H4, HEAD),
        make_para(P_BODY, INTRO),
        make_para(P_CODE, EXAMPLE),
        make_para(P_BODY, TYPES),
        make_para(P_BODY, CONVERSION),
    ]
    for offset, np in enumerate(new_paras):
        parent.insert(idx + offset, np)
    return {"new paragraphs": len(new_paras)}


if __name__ == "__main__":
    apply(__file__, transform)
