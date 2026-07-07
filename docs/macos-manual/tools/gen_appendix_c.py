#!/usr/bin/env python3
"""Regenerate src/93_appendix_c_commands.md from the live
macos_f90/macos_help.inc (single source of truth — same rule as the
old docx patch 90_appendix_c.py: never hand-edit Appendix C).

Each WRITE(*,*) literal becomes a line of the appendix.  Category
banners become #### headings; explanatory lines become prose;
everything else stays in fenced code blocks so alignment is preserved.

Run from anywhere:  python tools/gen_appendix_c.py
"""
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent
HELP_INC = HERE.parent.parent.parent / "macos_f90" / "macos_help.inc"
OUT = HERE.parent / "src" / "93_appendix_c_commands.md"

INTRO = (
    "The following is the complete command list as printed by the "
    "HELP command, regenerated from macos_f90/macos_help.inc (single "
    "source of truth).  Casing convention: \\<UPPER prefix\\>\\<lower "
    "tail\\> — the uppercase part is the minimum-match abbreviation "
    "tested by the dispatcher.  Tags: \\[Rx\\] needs a loaded "
    "prescription; \\[BLD\\] needs a built linear model; \\[DIFF\\] "
    "needs a propagated wavefront.  For a full per-command catalog — "
    "dialogs, behavior notes, and the mmacos/pymacos programming "
    "interfaces — see the companion **MACOS Command Reference** "
    "(macosCmdRef.pdf; source in docs/macos-manual/cmdref/).")

CATEGORY_BANNERS = {
    "SESSION & FILES",
    "PRESCRIPTION I/O",
    "SOURCE & WAVELENGTH",
    "RAY TRACING",
    "BEHAVIOR TOGGLES",
    "SURFACE DATA",
    "PERTURBATION",
    "LINEAR MODEL",
    "DIFFRACTION",
    "OUTPUTS (RAY-TRACE & DIFFRACTION)",
    "PLOT STYLE",
    "FILE OUTPUT",
    "WINDOW & POST-PROCESSING",
    "SYSTEM OPTIMIZATION",
    "MISC / DEBUG",
}


def parse_help_inc():
    """Yield (kind, text): 'banner', 'note', 'blank', or 'code'."""
    pattern = re.compile(r"WRITE\(\*,\*\)'(.*)'\s*$")
    for line in HELP_INC.read_text().splitlines():
        m = pattern.search(line)
        if not m:
            continue
        s = m.group(1).rstrip()
        stripped = s.strip()
        if not stripped:
            yield ("blank", "")
        elif stripped in CATEGORY_BANNERS:
            yield ("banner", stripped)
        elif (stripped.startswith("MACOS command summary")
              or stripped.startswith("Casing:")
              or stripped.startswith("Tags:")):
            yield ("note", stripped)
        else:
            yield ("code", s)


def main():
    out = ["<!-- GENERATED FILE — do not edit by hand.",
           "     Regenerate with: python tools/gen_appendix_c.py",
           "     Source: macos_f90/macos_help.inc  -->",
           "",
           "## APPENDIX C: Complete Command Reference",
           "",
           INTRO,
           ""]
    code = []

    def flush():
        if code:
            out.append("```")
            out.extend(code)
            out.append("```")
            out.append("")
            code.clear()

    for kind, text in parse_help_inc():
        if kind == "banner":
            flush()
            out.append(f"#### {text}")
            out.append("")
        elif kind == "note":
            flush()
            out.append(text)
            out.append("")
        elif kind == "blank":
            flush()
        else:
            code.append(text)
    flush()

    OUT.write_text("\n".join(out).rstrip() + "\n")
    print(f"wrote {OUT} from {HELP_INC}")


if __name__ == "__main__":
    main()
