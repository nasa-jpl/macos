"""
Phase: drop empty paragraphs left over from PDF page-break spacers.

Targets:
  * 'Normal'-styled paragraphs whose concatenated text is empty AND
    that contain no inline image / drawing / shape.
  * 'FrameContents'-styled paragraphs that are likewise content-empty.

Keeps a small number of empty BodyText / CodeBlock paragraphs that
provide intentional vertical whitespace between content blocks.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import (iter_paragraphs, get_style, get_text, has_image,
                          wtag)

DROP_STYLES = {"Normal", "FrameContents"}


def transform(root):
    to_remove = []
    for p in iter_paragraphs(root):
        if get_style(p) not in DROP_STYLES:
            continue
        if get_text(p).strip():
            continue
        if has_image(p):
            continue
        # Skip if the paragraph contains any object (math, OLE, etc.)
        if list(p.iter(wtag("object"))):
            continue
        to_remove.append(p)

    counts = {"Normal": 0, "FrameContents": 0}
    target_ids = {id(p) for p in to_remove}
    style_by_id = {id(p): get_style(p) for p in to_remove}
    for parent in root.iter():
        for child in list(parent):
            if child.tag == wtag("p") and id(child) in target_ids:
                counts[style_by_id[id(child)]] += 1
                parent.remove(child)

    return {
        "empty Normal removed":        counts["Normal"],
        "empty FrameContents removed": counts["FrameContents"],
    }


if __name__ == "__main__":
    apply(__file__, transform)
