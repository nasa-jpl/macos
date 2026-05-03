"""
Phase: drop the back-of-book index from the PDF.

The original FrameMaker PDF appended an alphabetical index (e.g.
'macos.glass 32, 46, 47', 'iZern 87', 'NOREGrid 27').  Acrobat carried
this through as Normal-styled paragraphs around para#~8052+.  Word can
regenerate an index from the new headings, so this stale text is just
clutter.

Heuristic: a paragraph is part of the index if it's styled Normal or
FrameContents, fits an "<term> <page-numbers>" shape, and lies in the
contiguous run of such paragraphs (we pick the largest contiguous
cluster to be conservative).
"""
import re
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import iter_paragraphs, get_style, get_text, wtag

INDEX_STYLES = {"Normal", "FrameContents"}

# An index entry: a string ending in one or more "page N" or "N, N, N"
# tokens.  Examples:
#   'macos.glass 32, 46, 47'
#   'NOREGrid 27'
#   'layer thickness 88 number of layers 88'
#   'iZern 87'
INDEX_RE = re.compile(
    r"^[\w\s,\.\(\)\-/'\"]+?\s+\d{1,4}(?:[,\s]+\d{1,4})*\s*$"
)


def looks_like_index_entry(text: str) -> bool:
    text = text.strip()
    if not text or len(text) < 4:
        return False
    if not INDEX_RE.match(text):
        return False
    # Must end with at least one numeric page reference.
    return bool(re.search(r"\d{1,4}\s*$", text))


def transform(root):
    paragraphs = [p for p in iter_paragraphs(root)]
    matches = []
    for i, p in enumerate(paragraphs):
        if get_style(p) not in INDEX_STYLES:
            continue
        text = get_text(p).strip()
        if looks_like_index_entry(text):
            matches.append(i)

    if not matches:
        return {"index entries dropped": 0}

    # Find the largest contiguous run of indices (allowing 1-paragraph
    # gaps for empty separators).
    best_start = best_end = matches[0]
    cur_start = cur_end = matches[0]
    for idx in matches[1:]:
        if idx - cur_end <= 3:   # tolerate small gaps
            cur_end = idx
        else:
            if cur_end - cur_start > best_end - best_start:
                best_start, best_end = cur_start, cur_end
            cur_start = cur_end = idx
    if cur_end - cur_start > best_end - best_start:
        best_start, best_end = cur_start, cur_end

    # Only drop if the run is sizable (at least 20 entries) -- protects
    # against accidentally matching a few legitimate "see page N"
    # cross-references in body text.
    run_len = best_end - best_start + 1
    if run_len < 20:
        return {"index entries dropped": 0,
                "longest candidate run": run_len}

    # Mark every paragraph in [best_start, best_end] for deletion.
    targets = set()
    for i in range(best_start, best_end + 1):
        targets.add(id(paragraphs[i]))

    dropped = 0
    for parent in root.iter():
        for child in list(parent):
            if child.tag == wtag("p") and id(child) in targets:
                parent.remove(child)
                dropped += 1

    return {
        "index entries dropped": dropped,
        "index span": f"para#{best_start}..{best_end} ({run_len} paragraphs)",
    }


if __name__ == "__main__":
    apply(__file__, transform)
