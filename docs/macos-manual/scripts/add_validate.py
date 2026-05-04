"""
add_validate.py — Insert VALIDATE-command documentation into the manual.

Three edits:
  1. Para in the Help-section summary list (line beginning
     "Prescription I/O: NEW, OLD, FID, ..."):  add "VALidate".
  2. After the "OLD and LOAD" Heading4 block (just before the MODIFY
     Heading4):  insert a Heading4 "VALIDATE" + a few BodyText paragraphs.
  3. The "PRESCRIPTION I/O" CodeBlock that mirrors the HELP output:
     extend the listing to include the VALidate line.

Run after `python scripts/unpack.py macosMan4.2beta.docx` and before
`python scripts/pack.py <output>.docx`.
"""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from scripts.docx_helpers import (
    load_doc, save_doc, get_text, iter_paragraphs,
    make_para, insert_after_para, replace_text_in_para,
    replace_para_text,
)

WORKING = PROJECT_ROOT / "working"


def main() -> None:
    tree, root = load_doc(WORKING)
    paras = list(iter_paragraphs(root))

    # 1) Help-section one-liner
    summary_hits = [
        p for p in paras
        if get_text(p).startswith("Prescription I/O: NEW, OLD, FID, SAVe")
    ]
    if not summary_hits:
        sys.exit("Could not find Prescription-I/O summary line")
    n = replace_text_in_para(
        summary_hits[0],
        "NEW, OLD, FID, SAVe",
        "NEW, OLD, FID, VALidate, SAVe",
    )
    print(f"  [1/3] summary line: {n} replacement(s)")

    # 2) New section right after "OLD and LOAD".  Anchor: the bare
    #    "MACOS> loa Cassegrain" line that closes the OLD/LOAD discussion.
    anchor = None
    for p in paras:
        if get_text(p).strip() == "MACOS> loa Cassegrain":
            anchor = p
            break
    if anchor is None:
        sys.exit("Could not find 'MACOS> loa Cassegrain' anchor")
    new_block = [
        make_para("BodyText", ""),
        make_para("Heading4", "VALIDATE"),
        make_para(
            "BodyText",
            "The VALidate command syntax-checks a .in prescription file "
            "without loading it into MACOS.  This is useful when sharing "
            "or hand-editing prescription files: rather than discovering "
            "a typo halfway through a load, the user can check the file "
            "first.  Like OLD, the command takes the filename (with or "
            "without the .in suffix):",
        ),
        make_para("BodyText", "MACOS> validate Cassegrain"),
        make_para(
            "BodyText",
            "If the file is well-formed, MACOS prints \"Cassegrain.in: "
            "OK\".  If a problem is detected, the message identifies the "
            "offending line and key, for example:",
        ),
        make_para(
            "BodyText",
            "Cassegrain.in: line 76: key \"ZernType\" has no value",
        ),
        make_para(
            "BodyText",
            "The same check runs automatically as part of the OLD and "
            "LOAD commands; if it fails, the load is aborted and the "
            "user is re-prompted for a different filename, leaving the "
            "previously loaded prescription untouched.",
        ),
        make_para("BodyText", ""),
    ]
    insert_after_para(root, anchor, new_block)
    print(f"  [2/3] inserted {len(new_block)} paragraphs after OLD/LOAD")

    # 3) HELP-mirror code block.  Find the CodeBlock that begins with
    #    "   NEW         - start a new optical system" and rewrite it
    #    with the VALidate line interleaved.
    help_old = (
        "   NEW         - start a new optical system"
        "   OLD <file>  - load a .in prescription"
        "   FID <id>    - load by file id (from RX listing)"
        "   SAVe        - write current system to .in"
        "   EXPort      - write select results to file"
        "   SHOw <elt>  - print all data for one element [Rx]"
        "   MODify <elt>- interactively edit one element  [Rx]"
    )
    help_new = (
        "   NEW         - start a new optical system"
        "   OLD <file>  - load a .in prescription"
        "   FID <id>    - load by file id (from RX listing)"
        "   VALidate <file> - syntax-check a .in prescription"
        "                     without loading it"
        "   SAVe        - write current system to .in"
        "   EXPort      - write select results to file"
        "   SHOw <elt>  - print all data for one element [Rx]"
        "   MODify <elt>- interactively edit one element  [Rx]"
    )
    fixed = False
    for p in paras:
        if get_text(p) == help_old:
            replace_para_text(p, help_new)
            fixed = True
            break
    if not fixed:
        sys.exit("Could not find HELP-mirror code block to update")
    print("  [3/3] HELP-mirror code block updated")

    save_doc(tree, WORKING)
    print("Saved working/word/document.xml")


if __name__ == "__main__":
    main()
