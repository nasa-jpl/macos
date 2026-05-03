"""
clean_4_1.py -- run the cleanup chain against macosMan4.1beta.docx.

Cycle:  unpack -> 00..08 cleanup patches -> pack -> macosMan4.1.docx

Distinct from run_all.py, which is the v3.2 -> v4.00 content-injection
pipeline.  These cleanup patches operate on macosMan4.1beta (the
hand-cleaned-up successor to macosMan4_0_styled) and produce the final
macosMan4.1.docx.

Re-runnable: each pass starts from a fresh unpack of 4.1beta.
"""
import sys
import subprocess
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
SCRIPTS      = PROJECT_ROOT / "scripts"
PATCHES_DIR  = Path(__file__).parent
SOURCE_DOCX  = PROJECT_ROOT / "macosMan4.1beta.docx"
OUTPUT_DOCX  = PROJECT_ROOT / "macosMan4.1.docx"

CLEANUP_PATCHES = [
    # Conservative cleanup pass.  Each patch is surgical and structural
    # only -- no style consolidation, no paragraph reclassification.
    "00_strip_pagination.py",         # remove 326 intermediate sectPr +
                                      # 101 column breaks
    "01_remove_horizontal_lines.py",  # remove 1,368 FrameMaker header/
                                      # footer horizontal-rule artifacts
    "09_smacos_modules.py",           # rewrite obsolete COMMON-blocks
                                      # sub-section to use Fortran modules

    # Disabled (style-changing patches that did too much; preserved in
    # the patches/ directory for reference but not in the pipeline):
    # "02_drop_empty_paragraphs.py"
    # "03_drop_index_pages.py"
    # "04_promote_misclassified.py"
    # "05_section_pagebreaks.py"
    # "08_consolidate_styles.py"
    # 06_unwrap_layout_tables.py was deleted entirely (destroyed real
    # data tables).
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
        "unpack 4.1beta")

    for patch_name in CLEANUP_PATCHES:
        run([sys.executable, str(PATCHES_DIR / patch_name)], patch_name)

    run([sys.executable, str(SCRIPTS / "pack.py"), str(OUTPUT_DOCX)],
        f"pack -> {OUTPUT_DOCX.name}")

    print(f"\nDone.  Output: {OUTPUT_DOCX}")
    print("Re-run scripts/inspect.py to see the cleaned style inventory.")


if __name__ == "__main__":
    main()
