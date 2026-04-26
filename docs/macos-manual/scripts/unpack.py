"""
unpack.py — Extract a .docx into working/ for XML editing.

Usage:
    python scripts/unpack.py                        # unpacks the project docx
    python scripts/unpack.py path/to/other.docx     # unpacks a specific file

Output lands in working/ (sibling of scripts/).
Always run this before editing; always run pack.py after.
"""

import zipfile
import shutil
import sys
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
DEFAULT_DOCX = PROJECT_ROOT / "macosMan3_2_styled.docx"
WORKING_DIR  = PROJECT_ROOT / "working"


def unpack(docx_path: Path, working_dir: Path) -> None:
    docx_path  = Path(docx_path)
    working_dir = Path(working_dir)

    if not docx_path.exists():
        sys.exit(f"Error: {docx_path} not found.")

    # Clean slate
    if working_dir.exists():
        shutil.rmtree(working_dir)
    working_dir.mkdir(parents=True)

    with zipfile.ZipFile(docx_path, "r") as zf:
        zf.extractall(working_dir)

    parts = list(working_dir.rglob("*"))
    xml_count = sum(1 for p in parts if p.suffix == ".xml")
    print(f"Unpacked {docx_path.name}  →  {working_dir}")
    print(f"  {xml_count} XML files, "
          f"{len(list((working_dir / 'word' / 'media').glob('*')) if (working_dir / 'word' / 'media').exists() else [])} media files")


if __name__ == "__main__":
    docx = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_DOCX
    unpack(docx, WORKING_DIR)
