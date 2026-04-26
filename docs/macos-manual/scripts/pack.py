"""
pack.py — Repack working/ back into a .docx file.

Usage:
    python scripts/pack.py                          # writes macosMan3_2_styled.docx
    python scripts/pack.py output/my_output.docx    # writes to a specific path

Always run unpack.py first, make your edits, then run this.
"""

import zipfile
import sys
import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
DEFAULT_OUTPUT = PROJECT_ROOT / "macosMan3_2_styled.docx"
WORKING_DIR    = PROJECT_ROOT / "working"


def pack(working_dir: Path, output_path: Path) -> None:
    working_dir = Path(working_dir)
    output_path = Path(output_path)

    if not working_dir.exists():
        sys.exit(f"Error: {working_dir} not found. Run unpack.py first.")

    output_path.parent.mkdir(parents=True, exist_ok=True)

    file_count = 0
    with zipfile.ZipFile(output_path, "w", zipfile.ZIP_DEFLATED) as zf:
        for file_path in sorted(working_dir.rglob("*")):
            if file_path.is_file():
                arcname = file_path.relative_to(working_dir)
                zf.write(file_path, arcname)
                file_count += 1

    size_kb = output_path.stat().st_size // 1024
    print(f"Packed {file_count} files  →  {output_path}  ({size_kb:,} KB)")


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_OUTPUT
    pack(WORKING_DIR, out)
