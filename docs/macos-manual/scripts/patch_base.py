"""
patch_base.py — Shared imports and main() runner for v4.0 patches.

Each numbered patch script in this directory imports `apply` from this
module and defines its own `transform(root)` that mutates the parsed
document tree.  The module's `apply()` handles load → transform → save
and prints a one-line status.

Usage from a patch script:

    from patch_base import apply, P_BODY, P_LIST, P_CODE, P_H3, P_H4
    from docx_helpers import make_para, find_para_by_text, ...

    def transform(root):
        ...
        return {"paragraphs added": n}

    if __name__ == "__main__":
        apply(__file__, transform)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from docx_helpers import load_doc, save_doc

PROJECT_ROOT = Path(__file__).parent.parent
WORKING_DIR  = PROJECT_ROOT / "working"

# Style-name shortcuts (from the document inventory)
P_BODY = "BodyText"
P_LIST = "ListParagraph"
P_CODE = "CodeBlock"
P_H1   = "Heading1"
P_H2   = "Heading2"
P_H3   = "Heading3"
P_H4   = "Heading4"
P_FIG  = "FigureCaption"


def apply(script_path: str, transform_fn) -> None:
    """Load, transform, save. Prints a banner with the patch name + stats."""
    name = Path(script_path).stem
    if not WORKING_DIR.exists():
        sys.exit("working/ not found — run: python scripts/unpack.py")
    tree, root = load_doc(WORKING_DIR)
    stats = transform_fn(root) or {}
    save_doc(tree, WORKING_DIR)
    summary = ", ".join(f"{v} {k}" for k, v in stats.items()) or "no changes"
    print(f"[{name}] {summary}")
