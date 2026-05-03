"""Patch A1: bump cover-page copyright '1992-1999' -> '1992-2026'."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply
from docx_helpers import iter_paragraphs, replace_text_in_para


def transform(root):
    hits = 0
    for p in iter_paragraphs(root):
        hits += replace_text_in_para(p, "1992-1999", "1992-2026")
    return {"runs touched": hits}


if __name__ == "__main__":
    apply(__file__, transform)
