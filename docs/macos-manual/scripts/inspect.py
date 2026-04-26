"""
inspect.py — Explore the document structure without changing anything.

Usage:
    python scripts/inspect.py                   # style summary
    python scripts/inspect.py headings          # full heading tree
    python scripts/inspect.py sample <style>    # sample paragraphs of a style
    python scripts/inspect.py search <text>     # find paragraphs containing text

Examples:
    python scripts/inspect.py headings
    python scripts/inspect.py sample Normal
    python scripts/inspect.py search "MACOS>"
"""

import sys
from pathlib import Path

# Allow running from project root or scripts/
sys.path.insert(0, str(Path(__file__).parent))
from docx_helpers import (load_doc, get_style, get_text,
                          get_run_info, print_style_summary, W, wtag)

PROJECT_ROOT = Path(__file__).parent.parent
WORKING_DIR  = PROJECT_ROOT / "working"

HEADING_LEVELS = {
    "Title":    0,
    "Heading1": 1, "Heading2": 2, "Heading3": 3,
    "Heading4": 4, "Heading5": 5,
}


def cmd_summary(root):
    print_style_summary(root)


def cmd_headings(root):
    print("Document heading tree:\n")
    for p in root.iter(wtag("p")):
        style = get_style(p)
        if style not in HEADING_LEVELS:
            continue
        text = get_text(p).strip()
        if not text:
            continue
        level  = HEADING_LEVELS[style]
        indent = "  " * level
        marker = "#" * max(level, 1)
        print(f"{indent}{marker} [{style}] {text}")


def cmd_sample(root, style_id: str, limit: int = 30):
    print(f"Sampling up to {limit} non-empty '{style_id}' paragraphs:\n")
    count = 0
    for p in root.iter(wtag("p")):
        if get_style(p) != style_id:
            continue
        text = get_text(p).strip()
        if not text:
            continue
        info = get_run_info(p)
        print(f"  font={info['font']!r:20s} sz={info['size']:>3s} "
              f"bold={info['bold']}  |  {text[:100]}")
        count += 1
        if count >= limit:
            break
    print(f"\n({count} shown)")


def cmd_search(root, query: str):
    query_lower = query.lower()
    print(f"Paragraphs containing {query!r}:\n")
    count = 0
    for p in root.iter(wtag("p")):
        text = get_text(p)
        if query_lower in text.lower():
            style = get_style(p)
            print(f"  [{style:20s}] {text.strip()[:110]}")
            count += 1
    print(f"\n{count} match(es)")


def main():
    if not WORKING_DIR.exists():
        sys.exit("working/ not found — run:  python scripts/unpack.py")

    _, root = load_doc(WORKING_DIR)

    args = sys.argv[1:]
    if not args or args[0] == "summary":
        cmd_summary(root)
    elif args[0] == "headings":
        cmd_headings(root)
    elif args[0] == "sample" and len(args) >= 2:
        limit = int(args[2]) if len(args) >= 3 else 30
        cmd_sample(root, args[1], limit)
    elif args[0] == "search" and len(args) >= 2:
        cmd_search(root, " ".join(args[1:]))
    else:
        print(__doc__)


if __name__ == "__main__":
    main()
