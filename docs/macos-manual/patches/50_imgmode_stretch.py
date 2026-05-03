"""
Patch C1: insert new sub-section 'Image polarity and stretch-aware
labels' immediately before §7.3 'Plot Type Commands'.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_LIST, P_CODE, P_H3
from docx_helpers import (iter_paragraphs, get_style, get_text,
                          find_para_parent, make_para)


HEAD = "Image polarity and stretch-aware labels"

POLARITY = (
    "Image polarity (IMGMODE).  Image-type plots (INTensity, "
    "PIXillated, OPD-as-image) can be rendered in either of two "
    "polarities:")

POL_BULLETS = [
    "Astronomical (default, 'ASTRO' or 'NEG'): large pixel values "
    "render dark, small values render bright.  Matches the PGPlot "
    "convention many MACOS users are accustomed to.",
    "Conventional ('CONV' or 'POS'): large values bright, small "
    "values dark.  Matches general scientific imaging.",
]

TOGGLE = "Toggle with the IMGMODE command:"

EXAMPLE = ("    MACOS> IMGMODE\n"
           "    Enter image polarity (NEG/POS, ASTRO/CONV): [NEG]: CONV")

PERSIST = (
    "The setting persists for the session.  It does not affect the "
    "non-image plot types (CONTOUR, SLICE, WIRE, SPOTDIAG, PLOTCOL), "
    "which always render with black ink on white.")

STRETCH = (
    "Stretch-aware colorbar labels.  When STRETCH is active (LOG10 or "
    "SQRT), the colorbar wedge label reflects the active stretch:")

STRETCH_BULLETS = [
    "INT, STRETCH=LIN: label 'Intensity'",
    "INT, STRETCH=LOG10: label 'LOG10 Intensity'",
    "INT, STRETCH=SQRT: label 'SQRT Intensity'",
    "PIX, STRETCH=LOG10: label 'LOG10 Pixel value'",
    "OPD with BaseUnits: label 'OPD (mm)' etc.",
]

LABEL_NOTE = (
    "The label is set by the command handler before the draw routine "
    "emits it, so the polarity / stretch and the label always match.")

ANNOT_INTRO = (
    "Bottom-of-plot annotation.  Most image plot types now print one "
    "or two annotation lines just below the colorbar:")

ANNOT_BULLETS = [
    "OPD: 'OPD=<rms> <unit> RMS, <pv> <unit> P-V'",
    "SPOT: 'RMS spot radius=<r> <unit>'",
    "INT: 'Peak pixel=<MaxInt>'",
    "PIX: 'Peak pixel=<maxval>'",
]


def transform(root):
    h = None
    for p in iter_paragraphs(root):
        if get_style(p) == "Heading3" and \
           "Plot Type Commands" in get_text(p):
            h = p
            break
    if h is None:
        raise RuntimeError("Couldn't find §7.3 Plot Type Commands")
    parent = find_para_parent(root, h)
    idx = list(parent).index(h)

    new_paras = [
        make_para(P_H3, HEAD),
        make_para(P_BODY, POLARITY),
    ]
    new_paras.extend(make_para(P_LIST, b) for b in POL_BULLETS)
    new_paras.append(make_para(P_BODY, TOGGLE))
    new_paras.append(make_para(P_CODE, EXAMPLE))
    new_paras.append(make_para(P_BODY, PERSIST))
    new_paras.append(make_para(P_BODY, STRETCH))
    new_paras.extend(make_para(P_LIST, b) for b in STRETCH_BULLETS)
    new_paras.append(make_para(P_BODY, LABEL_NOTE))
    new_paras.append(make_para(P_BODY, ANNOT_INTRO))
    new_paras.extend(make_para(P_LIST, b) for b in ANNOT_BULLETS)

    for offset, np in enumerate(new_paras):
        parent.insert(idx + offset, np)
    return {"new paragraphs": len(new_paras)}


if __name__ == "__main__":
    apply(__file__, transform)
