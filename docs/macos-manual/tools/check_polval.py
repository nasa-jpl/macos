#!/usr/bin/env python3
"""Refuse to build a polarization validation report that is out of date.

Four failure modes this catches, all of which would otherwise produce a
document that LOOKS authoritative:

  1. an unresolved @@TOKEN@@ in a rendered .md -- the renderer failed or was
     never run, and the report would ship with a placeholder in place of a
     measured number;
  2. a .md.in newer than its rendered .md -- someone edited the prose and
     did not re-render, so the text and the numbers disagree;
  3. a figure newer than generated/numbers.json -- the figures were
     regenerated but the numbers were not (or vice versa), so the panels and
     the tables describe different runs;
  4. the conventions table drifting between its two homes.  The executive
     summary reproduces it verbatim because it is review priority 1 and the
     summary is handed out alone; two copies of one table is exactly the
     arrangement that drifts, so the copies are compared byte for byte.

Usage:  check_polval.py <polvalDir>
Exit 0 = consistent.  Nonzero = do not build; run `make polval-regen`.
"""
from __future__ import annotations

import pathlib
import re
import sys

TOKEN = re.compile(r"@@([A-Za-z0-9_]+)@@")

# A block bracketed by these markers must be byte-identical everywhere it
# appears.  Used for the conventions table, which the executive summary
# reproduces verbatim.
MIRROR = re.compile(
    r"<!-- ([A-Z0-9-]+):BEGIN -->\n(.*?)<!-- \1:END -->", re.S
)


def check_mirrored_blocks(polval: pathlib.Path) -> list[str]:
    """Every marker-bracketed block with the same tag must match exactly."""
    seen: dict[str, tuple[str, str]] = {}   # tag -> (body, first file seen in)
    problems: list[str] = []
    for tpl in sorted(polval.glob("*.md.in")):
        for tag, body in MIRROR.findall(tpl.read_text()):
            if tag not in seen:
                seen[tag] = (body, tpl.name)
            elif seen[tag][0] != body:
                problems.append(
                    f"{tpl.name}: block {tag} has drifted from the copy in "
                    f"{seen[tag][1]} -- the two must be byte-identical"
                )
    return problems


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    polval = pathlib.Path(sys.argv[1])
    problems: list[str] = []

    templates = sorted(polval.glob("*.md.in"))
    if not templates:
        print(f"check_polval: no templates under {polval}", file=sys.stderr)
        return 2

    for tpl in templates:
        md = tpl.with_suffix("")
        if not md.exists():
            problems.append(f"{md.name} has never been rendered from {tpl.name}")
            continue
        if tpl.stat().st_mtime > md.stat().st_mtime + 1.0:
            problems.append(
                f"{tpl.name} is newer than {md.name} -- prose edited without re-rendering"
            )
        left = sorted(set(TOKEN.findall(md.read_text())))
        if left:
            problems.append(
                f"{md.name} carries unresolved tokens: " + ", ".join(left)
            )

    problems.extend(check_mirrored_blocks(polval))

    numbers = polval / "generated" / "numbers.json"
    if not numbers.exists():
        problems.append("generated/numbers.json is missing")
    else:
        nt = numbers.stat().st_mtime
        newer = [
            m.name
            for m in sorted((polval / "media").glob("*.png"))
            if m.stat().st_mtime > nt + 1.0
        ]
        if newer:
            problems.append(
                "figures newer than numbers.json: " + ", ".join(newer)
            )

    if problems:
        print("check_polval: report is NOT consistent --", file=sys.stderr)
        for p in problems:
            print(f"  * {p}", file=sys.stderr)
        print(
            "\nRun `make polval-regen` (needs MATLAB + a built mmacos mex).",
            file=sys.stderr,
        )
        return 1

    print(f"check_polval: {len(templates)} sections consistent with numbers.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
