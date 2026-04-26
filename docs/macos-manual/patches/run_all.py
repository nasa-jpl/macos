"""
run_all.py -- run the full v3.2 -> v4.00 patch chain in order.

Cycle:  unpack -> 01..90 patches -> pack -> macosMan4_0_styled.docx

Output: macosMan4_0_styled.docx (sibling of the styled v3.2 source).

Re-runnable: re-running unpacks fresh from the v3.2 source each time,
so partial results from a previous run never leak in.
"""
import sys
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
SCRIPTS      = PROJECT_ROOT / "scripts"
PATCHES_DIR  = Path(__file__).parent
SOURCE_DOCX  = PROJECT_ROOT / "macosMan3_2_styled.docx"
OUTPUT_DOCX  = PROJECT_ROOT / "macosMan4_0_styled.docx"

PATCHES = [
    "01_cover_copyright.py",
    "10_v4_overview.py",
    "20_help_rebuild.py",
    "21_ctrlc_blank.py",
    "30_segmirmaker_note.py",
    "31_sparse_zernike.py",
    "32_freeform_surfaces.py",
    "40_segraytrace_status.py",
    "50_imgmode_stretch.py",
    "60_gmi_section.py",
    "90_appendix_c.py",
]


def run(cmd: list, label: str) -> None:
    print(f"\n--- {label} ---")
    r = subprocess.run(cmd, cwd=str(PROJECT_ROOT))
    if r.returncode != 0:
        sys.exit(f"FAILED: {label} (rc={r.returncode})")


def main():
    if not SOURCE_DOCX.exists():
        sys.exit(f"Source not found: {SOURCE_DOCX}")

    run([sys.executable, str(SCRIPTS / "unpack.py"), str(SOURCE_DOCX)],
        "unpack v3.2 source")

    for patch_name in PATCHES:
        run([sys.executable, str(PATCHES_DIR / patch_name)], patch_name)

    run([sys.executable, str(SCRIPTS / "pack.py"), str(OUTPUT_DOCX)],
        "pack -> macosMan4_0_styled.docx")

    print(f"\nDone.  Output: {OUTPUT_DOCX}")
    print("Open in LibreOffice or Word to verify.")


if __name__ == "__main__":
    main()
