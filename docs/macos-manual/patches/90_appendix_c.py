"""
Patch K1: append 'APPENDIX C — Complete Command Reference' at the end
of the document, regenerated from the live macos_help.inc (single
source of truth).

Each WRITE(*,*) literal in macos_help.inc becomes a line in the
appendix, preserving the same grouping and casing the user sees from
the live HELP command.  Section banners (lines that start with ALL-CAPS
or 'TAGS:' etc.) become Heading4 paragraphs; everything else is a
CodeBlock-style fixed-width paragraph so the alignment stays intact.
"""
import sys
import re
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_CODE, P_H2, P_H3, P_H4
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, make_para)

HELP_INC = (Path("/home/dcr/dev/macos") / "macos_f90" / "macos_help.inc")

INTRO = (
    "The following is the complete command list as printed by the "
    "HELP command, regenerated from macos_f90/macos_help.inc (single "
    "source of truth).  Casing convention: <UPPER prefix><lower tail> "
    "-- the uppercase part is the minimum-match abbreviation tested "
    "by the dispatcher.  Tags: [Rx] needs a loaded prescription; "
    "[BLD] needs a built linear model; [DIFF] needs a propagated "
    "wavefront.")

# Lines that should render as Heading4 (category banners).
CATEGORY_BANNERS = {
    "SESSION & FILES",
    "PRESCRIPTION I/O",
    "SOURCE & WAVELENGTH",
    "RAY TRACING",
    "BEHAVIOR TOGGLES",
    "SURFACE DATA",
    "PERTURBATION",
    "LINEAR MODEL",
    "DIFFRACTION",
    "OUTPUTS (RAY-TRACE & DIFFRACTION)",
    "PLOT STYLE",
    "FILE OUTPUT",
    "WINDOW & POST-PROCESSING",
    "SYSTEM OPTIMIZATION",
    "MISC / DEBUG",
}


def parse_help_inc():
    """Yield (kind, text) pairs from macos_help.inc.

    `kind` is one of 'banner' (Heading4), 'note' (BodyText),
    or 'code' (CodeBlock).  Skips empty WRITE(*,*) ' ' lines.
    """
    pattern = re.compile(r"WRITE\(\*,\*\)'(.*)'\s*$")
    for line in HELP_INC.read_text().splitlines():
        m = pattern.search(line)
        if not m:
            continue
        s = m.group(1).rstrip()
        stripped = s.strip()
        if not stripped:
            yield ("blank", "")
            continue
        if stripped in CATEGORY_BANNERS:
            yield ("banner", stripped)
        elif stripped.startswith("MACOS command summary"):
            yield ("note", stripped)
        elif stripped.startswith("Casing:") or stripped.startswith("Tags:"):
            yield ("note", stripped)
        else:
            yield ("code", s)


def transform(root):
    # Append at very end of body. Find the body element.
    body = None
    for child in root:
        if child.tag == wtag("body"):
            body = child
            break
    if body is None:
        raise RuntimeError("Couldn't find <w:body>")

    # The last child of body is sectPr; insert before it.
    children = list(body)
    sectpr_idx = len(children)
    for i, c in enumerate(children):
        if c.tag == wtag("sectPr"):
            sectpr_idx = i
            break

    new_paras = [
        make_para(P_H2, "APPENDIX C  Complete Command Reference"),
        make_para(P_BODY, INTRO),
    ]
    grouped_codes = []
    for kind, text in parse_help_inc():
        if kind == "banner":
            # Flush any accumulated code lines first.
            if grouped_codes:
                new_paras.append(
                    make_para(P_CODE, "\n".join(grouped_codes)))
                grouped_codes = []
            new_paras.append(make_para(P_H4, text))
        elif kind == "note":
            if grouped_codes:
                new_paras.append(
                    make_para(P_CODE, "\n".join(grouped_codes)))
                grouped_codes = []
            new_paras.append(make_para(P_BODY, text))
        elif kind == "blank":
            if grouped_codes:
                new_paras.append(
                    make_para(P_CODE, "\n".join(grouped_codes)))
                grouped_codes = []
        else:  # code
            grouped_codes.append(text)
    if grouped_codes:
        new_paras.append(make_para(P_CODE, "\n".join(grouped_codes)))

    for offset, np in enumerate(new_paras):
        body.insert(sectpr_idx + offset, np)
    return {"new paragraphs": len(new_paras),
            "source": str(HELP_INC)}


if __name__ == "__main__":
    apply(__file__, transform)
