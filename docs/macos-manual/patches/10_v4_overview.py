"""
Patch A2: Replace heading '1.3 New Features of Version 3.2' with the
new v4.00 overview.  The original 1.3 list is renamed
'1.3 (historical) Changes in v3.2' below the new content so prior
release-note context isn't lost.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_LIST, P_H3
from docx_helpers import (find_para_by_text, insert_after_para,
                          replace_para_text, make_para)

OVERVIEW_BODY = (
    "MACOS 4.00 (April 2026) is a feature release built on the "
    "Fortran 90 / SMACOS infrastructure introduced in v3.2.  It adds a "
    "new general-purpose composite surface (FreeForm), modernizes the "
    "graphics layer, expands the SMACOS / GMI MATLAB interface, and "
    "adds new utility commands plus a separate segmented-mirror "
    "prescription generator.")

INTRO_LINE = "A summary of changes since MACOS 3.2:"

V4_BULLETS = [
    "FreeForm composite surface (SrfType = 14). Single surface type "
    "combining a conic base with a Cartesian monomial sag (Mon), a "
    "free-form Zernike sag (FF), and an interpolated grid sag, all on "
    "independent coordinate frames.  See §4.5.11.",

    "Per-element Zernike inputs.  Sparse \"modes + coefs\" form "
    "MonZernCoef / MonZernModes and FFZernCoef / FFZernModes is "
    "accepted directly in .in files; conversion to monomial form "
    "happens internally at trace time.  See §4.5.4.",

    "Per-ray status tracking.  Every ray records why it failed "
    "(Miss / Obscured / Bracket / MaxIter) and at which element.  The "
    "end-of-trace WARN summary breaks down failures by category and "
    "throttles bracket / iter messages after the first 20.",

    "Graphics modernization.  Built-in giza graphics backend "
    "replacing legacy PGPLOT (still selectable at build time).  New "
    "IMGMODE command toggles between astronomical (negative, "
    "large->dark) and conventional (positive, large->light) display "
    "polarity.  Colorbar wedge labels reflect the active STRETCH "
    "(LIN / LOG10 / SQRT).  See §7.x.",

    "SEGRAYTRACE command.  Trace a single ray through the geometric "
    "center of a chosen segment element; print ray state at every "
    "downstream surface.  Useful for verifying segmented-mirror "
    "prescriptions.  See §5.2.",

    "SegMirMaker external tool.  Standalone utility, in "
    "MACOS_resources/segmirmaker/, that generates segmented-mirror "
    ".presc and Hx.m files.  Replaces SMPGe (1992) and adds support "
    "for FreeForm parents.  See §4.4.11.",

    "GMI MATLAB interface upgrades.  New 14th MEX argument pmonzern "
    "for FreeForm MonZernCoef perturbations (parallel to pzern).  New "
    "param.monzernSrf element list.  Nominal-state save/restore "
    "extended to cover FreeForm coordinate frames so per-element "
    "rigid-body perturbations no longer leak across GMI calls.  "
    "See §9.",

    "Graceful Ctrl-C.  SIGINT installs a clean-exit handler at "
    "start-up; pressing Ctrl-C at the MACOS> prompt prints "
    "'MACOS: interrupted' and exits cleanly, replacing the "
    "Intel-runtime 'forrtl: severe (69)' traceback.",

    "HELP command rebuilt.  All ~140 dispatched commands are now "
    "listed (vs. ~75 in v3.2), grouped into 15 categories with "
    "one-line effect descriptions and prerequisite tags "
    "([Rx] / [BLD] / [DIFF]).  See §3.2 and Appendix C.",
]


def transform(root):
    h = find_para_by_text(root, "New Features of Version 3.2", style=P_H3)
    if h is None:
        raise RuntimeError("Couldn't find heading 'New Features of Version 3.2'")

    # Rename the existing heading to v4.00.
    replace_para_text(h, "New Features of Version 4.00", style=P_H3)

    new_paras = [make_para(P_BODY, OVERVIEW_BODY),
                 make_para(P_BODY, INTRO_LINE)]
    new_paras.extend(make_para(P_LIST, b) for b in V4_BULLETS)
    # Add a separator heading for the historical list that follows.
    new_paras.append(make_para(P_H3,
        "1.3 (historical) Changes since MACOS 2.8"))

    insert_after_para(root, h, new_paras)
    return {"new heading": 1, "body paragraphs": 2,
            "bullets": len(V4_BULLETS), "historical heading": 1}


if __name__ == "__main__":
    apply(__file__, transform)
