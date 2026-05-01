"""
transform_template.py — Copy this to start a new transformation script.

Usage:
    python scripts/unpack.py
    python scripts/transform_template.py
    python scripts/pack.py

Then open macosMan3_2_styled.docx in LibreOffice to verify.
"""

import sys
from pathlib import Path
from collections import Counter

sys.path.insert(0, str(Path(__file__).parent))
from docx_helpers import (load_doc, save_doc, get_style, set_style,
                          get_text, get_run_info, has_image,
                          iter_paragraphs, print_style_summary, wtag)

PROJECT_ROOT = Path(__file__).parent.parent
WORKING_DIR  = PROJECT_ROOT / "working"


def transform(root):
    """
    Make your changes here. Return a stats dict for reporting.
    
    Common patterns:
    
    # Restyle a paragraph
    set_style(para, "Heading3")
    
    # Remove a paragraph from its parent
    # (collect ids first, then remove in a second pass)
    to_remove = set()
    to_remove.add(id(para))
    for parent in root.iter():
        for child in list(parent):
            if child.tag == wtag("p") and id(child) in to_remove:
                parent.remove(child)
    
    # Replace text in a run
    for t_el in para.iter(wtag("t")):
        if t_el.text:
            t_el.text = t_el.text.replace("old", "new")
    """
    stats = Counter()

    for para in iter_paragraphs(root):
        style = get_style(para)
        text  = get_text(para).strip()

        # ── your logic here ──────────────────────────────────────────────

        pass

    return stats


def main():
    if not WORKING_DIR.exists():
        sys.exit("working/ not found — run:  python scripts/unpack.py")

    print("Loading document...")
    tree, root = load_doc(WORKING_DIR)

    print("Before:")
    print_style_summary(root)

    print("\nTransforming...")
    stats = transform(root)

    print("\nChanges made:")
    for action, count in stats.most_common():
        print(f"  {count:5d}  {action}")

    print("\nAfter:")
    print_style_summary(root)

    print("\nSaving...")
    save_doc(tree, WORKING_DIR)
    print("Done. Run:  python scripts/pack.py")


if __name__ == "__main__":
    main()
