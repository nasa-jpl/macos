"""
Patches H + G: insert SEGRAYTRACE command sub-section and per-ray
status reporting sub-section under §5.2 'Ray-Trace Commands'.
Both go just before §5.3 'Beam Set-Up Commands'.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_LIST, P_CODE, P_H4
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, make_para)


SEG_HEAD = "SEGRAYTRACE"

SEG_INTRO = (
    "Trace a single ray that passes through the geometric center of a "
    "chosen segment element (its RptElt), printing the ray state at "
    "every downstream surface.  Useful for verifying that segmented-"
    "mirror prescriptions place segments where you expect.")

SEG_EXAMPLE = (
    "    MACOS> SEGRAYTRACE\n"
    "    Enter segment number (0=quit): [0]: 5\n"
    "    Ray  1 segment from Element 0 (InputRay) to Element 5 (Seg5):\n"
    "       Starting point=  ...\n"
    "       End point     =  ...   (should be RptElt(:,5))\n"
    "       Direction     =  ...\n"
    "    ...\n"
    "    Enter segment number (0=quit): [0]: 0\n"
    "    MACOS>")

SEG_NOTE = (
    "The source ray's position is automatically offset from ChfRayPos "
    "so the resulting ray (parallel to ChfRayDir for a collimated "
    "source, or aimed from ChfRayPos for a point source) lands "
    "exactly on the segment's recorded RptElt.  The segment-dispatch "
    "path is forced to fire only on the chosen segment via "
    "RayToSegMap; other segments are skipped.")

SEG_TOLERANCE = (
    "The trace endpoint at the segment's element should match RptElt "
    "to within SFFZPSolve's tolerance (~1e-11 mm).  Larger mismatches "
    "diagnose generator bugs:")

SEG_BUGS = [
    "~180-deg rotational mismatch (e.g. Seg2's endpoint matches "
    "Seg5's RptElt): the segment-numbering psi-flip in the generator "
    "didn't fire when it should have.",
    "~mm radial mismatch in the surface-normal direction: the "
    "generator computed RptElt on the conic only, missing the FF "
    "Zernike sag.",
]

SEG_PRECONDITION = (
    "The command requires an Rx loaded with at least one Segment "
    "element; otherwise it prints 'Element N is not a Segment' and "
    "re-prompts.")


STATUS_HEAD = "Per-ray status reporting"

STATUS_INTRO = "Every ray in a trace records its terminal status:"

STATUS_VALUES = [
    "OK (0) -- successfully terminated at last requested element",
    "Obscured (1) -- blocked by an obscuration or aperture",
    "Miss (2) -- missed a surface (no intersection found)",
    "Bracket (3) -- root-finding bracket failure on a composite "
    "(FreeForm / aspheric / grid) surface",
    "MaxIter (4) -- root-finder hit its iteration limit",
    "Undef (5) -- initial state, never updated",
]

STATUS_NOTE = (
    "The end-of-trace WARN summary breaks down failures by category "
    "and records the offending element index for the first failure of "
    "each ray.  Bracket / iter messages are throttled after the first "
    "20 to keep the terminal usable when many rays fail.")


def transform(root):
    # Anchor: the §5.3 'Beam Set-Up Commands' Heading3.  We insert
    # immediately before it.
    h = None
    for p in iter_paragraphs(root):
        if get_style(p) == "Heading3" and \
           "Beam Set-Up Commands" in get_text(p):
            h = p
            break
    if h is None:
        raise RuntimeError("Couldn't find §5.3 Beam Set-Up Commands")
    parent = find_para_parent(root, h)
    idx = list(parent).index(h)

    new_paras = [
        make_para(P_H4, SEG_HEAD),
        make_para(P_BODY, SEG_INTRO),
        make_para(P_CODE, SEG_EXAMPLE),
        make_para(P_BODY, SEG_NOTE),
        make_para(P_BODY, SEG_TOLERANCE),
    ]
    new_paras.extend(make_para(P_LIST, b) for b in SEG_BUGS)
    new_paras.append(make_para(P_BODY, SEG_PRECONDITION))

    new_paras.append(make_para(P_H4, STATUS_HEAD))
    new_paras.append(make_para(P_BODY, STATUS_INTRO))
    new_paras.extend(make_para(P_LIST, v) for v in STATUS_VALUES)
    new_paras.append(make_para(P_BODY, STATUS_NOTE))

    for offset, np in enumerate(new_paras):
        parent.insert(idx + offset, np)
    return {"new paragraphs": len(new_paras)}


if __name__ == "__main__":
    apply(__file__, transform)
