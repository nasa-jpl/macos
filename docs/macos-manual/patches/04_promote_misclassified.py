"""
Phase: promote misclassified Normal paragraphs to their semantic style.

Audit findings: a few Normal-styled paragraphs were really headings or
misc body text whose style got lost in the PDF round-trip.  The
distinguishing feature is the run formatting (Helvetica 18 bold).

Rules:
  * Normal + Helvetica + 18-pt bold -> Heading3   (e.g. 'Element Types')
  * Normal + Helvetica + 24-pt+ bold -> Heading2  (e.g. 'SECTION 7')
  * Normal + Helvetica + 14-pt bold -> Heading4
  * Other Normal: leave alone (math equations, diagram labels --
    require human review)
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import (iter_paragraphs, get_style, get_text,
                          get_run_info, set_style)


def transform(root):
    promoted = {"Heading2": 0, "Heading3": 0, "Heading4": 0}
    for p in iter_paragraphs(root):
        if get_style(p) != "Normal":
            continue
        text = get_text(p).strip()
        if not text:
            continue
        info = get_run_info(p)
        font = info["font"]
        try:
            sz = int(info["size"]) if info["size"] else 0
        except ValueError:
            sz = 0
        if not info["bold"]:
            continue
        if "Helvetica" not in (font or ""):
            continue
        if sz >= 24:
            set_style(p, "Heading2"); promoted["Heading2"] += 1
        elif sz >= 18:
            set_style(p, "Heading3"); promoted["Heading3"] += 1
        elif sz >= 14:
            set_style(p, "Heading4"); promoted["Heading4"] += 1

    return {f"-> {k}": v for k, v in promoted.items() if v}


if __name__ == "__main__":
    apply(__file__, transform)
