"""
Phase 4 final pass: rename the surviving `TableParagraph` style to
`TableCell` (semantic name) for the few paragraphs that genuinely live
inside <w:tbl>.

After 06_unwrap_layout_tables.py, the only TableParagraph occurrences
should be in real data tables.  Renaming makes the document self-
documenting.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import (iter_paragraphs, get_style, set_style, wtag)


def in_table(p, parent_map):
    cur = parent_map.get(id(p))
    while cur is not None:
        if cur.tag == wtag("tbl"):
            return True
        cur = parent_map.get(id(cur))
    return False


def transform(root):
    # Build parent map.
    pm = {}
    for parent in root.iter():
        for child in parent:
            pm[id(child)] = parent

    renamed_in_tbl     = 0
    renamed_out_of_tbl = 0
    for p in iter_paragraphs(root):
        if get_style(p) != "TableParagraph":
            continue
        if in_table(p, pm):
            set_style(p, "TableCell")
            renamed_in_tbl += 1
        else:
            # Stray TableParagraph outside any table -- promote to BodyText.
            set_style(p, "BodyText")
            renamed_out_of_tbl += 1

    return {
        "TableParagraph -> TableCell (in real tables)":   renamed_in_tbl,
        "TableParagraph -> BodyText (orphaned)":          renamed_out_of_tbl,
    }


if __name__ == "__main__":
    apply(__file__, transform)
