"""
Patch D: append Ctrl-C, blank-command, and unknown-command notes to
the end of §3.4 'Entering Commands' (just before the next H3 heading).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from patch_base import apply, P_BODY, P_H3, P_H4
from docx_helpers import (iter_paragraphs, get_style, get_text, wtag,
                          find_para_parent, make_para)

NOTES = [
    ("Interrupting MACOS (Ctrl-C)",
     "Pressing Ctrl-C at the MACOS> prompt -- or during a long trace, "
     "build, or propagation -- prints 'MACOS: interrupted' and exits "
     "with status 0.  Internally MACOS installs a SIGINT handler at "
     "startup that uses _exit(0) so it is async-signal-safe and fires "
     "regardless of what the runtime is doing.  This replaces the "
     "Intel-runtime 'forrtl: severe (69)' traceback that earlier MACOS "
     "versions produced on Ctrl-C."),
    ("Empty input",
     "Pressing return at the MACOS> prompt with no command is a no-op "
     "-- useful for clearing the input buffer or scrolling history "
     "without re-running the last command."),
    ("Unknown commands",
     "Unknown commands print a single-line warning rather than the "
     "verbose multi-line traceback that v3.2 produced."),
]


def transform(root):
    # Find '### Entering Commands' heading.
    h = None
    for p in iter_paragraphs(root):
        if get_style(p) == P_H3 and "Entering Commands" in get_text(p):
            h = p
            break
    if h is None:
        raise RuntimeError("Couldn't find §3.4 Entering Commands")

    # Walk forward until next H3 heading; insert before it.
    parent = find_para_parent(root, h)
    children = list(parent)
    idx = children.index(h) + 1
    while idx < len(children):
        c = children[idx]
        if c.tag == wtag("p") and get_style(c) == P_H3:
            break
        idx += 1

    new_paras = []
    for label, body in NOTES:
        new_paras.append(make_para(P_H4, label))
        new_paras.append(make_para(P_BODY, body))

    for offset, np in enumerate(new_paras):
        parent.insert(idx + offset, np)
    return {"H4 sub-headings": len(NOTES), "body paragraphs": len(NOTES)}


if __name__ == "__main__":
    apply(__file__, transform)
